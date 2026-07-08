"""SDF/MOL2 helpers for the bulk ligand prep pipeline."""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path
from typing import List, Tuple

logger = logging.getLogger(__name__)


OBABEL_TIMEOUT_S = 900
LIGPREP_OBABEL_TIMEOUT_SEC = OBABEL_TIMEOUT_S


def rename_mol2_with_prefix(mol2_path: Path, target_stem: str) -> Path:
    """Rename a MOL2 (and optional .name sidecar) to a deterministic stem."""
    target = mol2_path.with_name(f"{target_stem}{mol2_path.suffix}")
    if target == mol2_path:
        return mol2_path
    if target.exists():
        target.unlink()
    mol2_path.rename(target)
    name_path = mol2_path.with_suffix(".name")
    if name_path.exists():
        name_path.rename(target.with_suffix(".name"))
    return target


def _run_obabel(cmd: List[str], timeout_sec: int) -> bool:
    print("Running Open Babel:", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True, timeout=timeout_sec)
        return True
    except subprocess.TimeoutExpired:
        print("Open Babel timed out.")
        return False
    except subprocess.CalledProcessError as e:
        print("Open Babel failed with code:", e.returncode)
        return False


def _attempt_obabel_series(
    base_cmd: List[str],
    timeout_sec: int,
    threads_list: List[int],
    use_fast_first: bool = True,
) -> bool:
    attempts: List[List[str]] = []
    if use_fast_first:
        for j in threads_list:
            attempts.append(base_cmd + ["--fast", "-j", str(j)])
    for j in threads_list:
        cmd = [c for c in base_cmd if c != "--fast"]
        cmd += ["-j", str(j)]
        attempts.append(cmd)

    total = len(attempts)
    for i, cmd in enumerate(attempts, 1):
        print(f"[OBabel attempt {i}/{total}]")
        if _run_obabel(cmd, timeout_sec):
            return True
        print("Retrying with a more conservative setting...")
    return False


def _sdf_to_mol2(
    sdf_path: Path,
    mol2_path: Path,
    obabel_exe: str,
    timeout_sec: int = LIGPREP_OBABEL_TIMEOUT_SEC,
) -> Tuple[bool, str]:
    """Convert a single SDF to MOL2 using Open Babel."""

    mol2_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        obabel_exe,
        "-isdf",
        str(sdf_path),
        "-omol2",
        "-O",
        str(mol2_path),
        "--gen3d",
    ]

    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout_sec,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        logger.warning(
            "[ligprep] per-ligand sdf->mol2 timed out for %s (timeout=%s s)",
            sdf_path,
            timeout_sec,
        )
        return False, f"timeout: {exc}"

    ok = proc.returncode == 0
    stderr_text = (proc.stderr or "").strip()
    if not ok:
        logger.warning(
            "[ligprep] per-ligand sdf->mol2 non-zero exit for %s: rc=%s stderr=%s",
            sdf_path,
            proc.returncode,
            stderr_text[:200],
        )
    return ok, stderr_text


def count_sdf_records(sdf_path: Path) -> int:
    n = 0
    with open(sdf_path, "r", errors="ignore") as fh:
        for line in fh:
            if line.startswith("$$$$"):
                n += 1
    return n


def split_sdf_into_chunks(
    sdf_path: Path, out_dir: Path, chunk_size: int = 1000
) -> List[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    chunks: List[Path] = []
    idx = 0
    mol_buf: List[str] = []
    mols_in_chunk = 0

    def flush_chunk() -> None:
        nonlocal idx, mol_buf, mols_in_chunk
        if not mol_buf:
            return
        idx += 1
        out_path = out_dir / f"{sdf_path.stem}_chunk{idx:04d}.sdf"
        with open(out_path, "w", encoding="utf-8") as out:
            out.write("".join(mol_buf))
        chunks.append(out_path)
        mol_buf = []
        mols_in_chunk = 0

    with open(sdf_path, "r", errors="ignore") as fh:
        cur: List[str] = []
        for line in fh:
            cur.append(line)
            if line.startswith("$$$$"):
                mol_buf.extend(cur)
                cur = []
                mols_in_chunk += 1
                if mols_in_chunk >= chunk_size:
                    flush_chunk()
        if cur:
            mol_buf.extend(cur)
        flush_chunk()
    return chunks


def convert_sdf_to_mol2_split_parallel(
    sdf_path: Path,
    mol2_output_dir: Path,
    obabel_exe: str,
    threads: int = 32,
    timeout_sec: int = 36000,
    chunk_size: int = 1000,
) -> List[Path]:
    mol2_output_dir.mkdir(parents=True, exist_ok=True)

    chunk_dir = mol2_output_dir / "_sdf_chunks"
    if chunk_dir.exists():
        shutil.rmtree(chunk_dir)
    chunk_dir.mkdir(parents=True, exist_ok=True)

    n_mols = count_sdf_records(sdf_path)
    print(f"SDF has ~{n_mols} molecules")

    chunk_paths = split_sdf_into_chunks(sdf_path, chunk_dir, chunk_size=chunk_size)
    print(f"Split into {len(chunk_paths)} chunk(s) of up to {chunk_size} molecules")

    all_out: List[Path] = []

    t1 = threads if threads >= 16 else 16
    t2 = max(16, t1 // 2)
    threads_list: List[int] = []
    for t in [t1, t2, 16]:
        if t not in threads_list:
            threads_list.append(t)

    for ci, chunk in enumerate(chunk_paths, 1):
        prefix = mol2_output_dir / f"mol2_chunk{ci:04d}_"
        base_cmd = [
            obabel_exe,
            "-isdf",
            str(chunk),
            "--gen3d",
            "-omol2",
            "-m",
            "-O",
            str(prefix) + ".mol2",
        ]
        print(
            f"[ligprep] pre-MOL2-write: chunk={ci} sdf={chunk.name} -> prefix={prefix.name}"
        )
        ok = _attempt_obabel_series(
            base_cmd,
            timeout_sec=timeout_sec,
            threads_list=threads_list,
            use_fast_first=True,
        )
        if not ok:
            print(f"Chunk {ci} failed entirely; moving on.")
            continue

        out_files = sorted(mol2_output_dir.glob(f"mol2_chunk{ci:04d}_*.mol2"))
        print(f"Chunk {ci}: wrote {len(out_files)} mol2 files")
        for _mol2 in out_files[:4]:
            try:
                n_mol2_atoms = 0
                with open(_mol2, "r", errors="ignore") as fh:
                    in_atoms = False
                    for ln in fh:
                        s = ln.strip()
                        if s.startswith("@<TRIPOS>ATOM"):
                            in_atoms = True
                            continue
                        if s.startswith("@<TRIPOS>") and in_atoms:
                            break
                        if in_atoms and s and s[0].isdigit():
                            n_mol2_atoms += 1
                print(f"[ligprep] MOL2_written≈{n_mol2_atoms} file={_mol2.name}")
            except Exception as e:
                print(f"[ligprep] MOL2_count_error {_mol2.name}: {e}")
        all_out.extend(out_files)

    print(f"Total MOL2 files: {len(all_out)}")
    return all_out


def _obabel_convert_chunk(chunk_sdf: Path, out_prefix: Path, obabel_exe: str) -> int:
    cmd = [
        obabel_exe,
        "-isdf",
        str(chunk_sdf),
        "-omol2",
        "-m",
        "-O",
        str(out_prefix) + ".mol2",
        "-j",
        "50",
    ]
    print("Running Open Babel:", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True, timeout=OBABEL_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        print(f"Open Babel timed out on {chunk_sdf.name}")
        return 0
    except subprocess.CalledProcessError as e:
        print(f"Open Babel failed ({e.returncode}) on {chunk_sdf.name}")
        return 0

    written = len(list(out_prefix.parent.glob(out_prefix.name + "_*.mol2")))
    return written


__all__ = [
    "OBABEL_TIMEOUT_S",
    "LIGPREP_OBABEL_TIMEOUT_SEC",
    "_run_obabel",
    "_attempt_obabel_series",
    "_sdf_to_mol2",
    "count_sdf_records",
    "split_sdf_into_chunks",
    "convert_sdf_to_mol2_split_parallel",
    "_obabel_convert_chunk",
]
