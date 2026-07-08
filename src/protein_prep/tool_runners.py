from __future__ import annotations

import os
import shutil
import subprocess
import sys
import sysconfig
from pathlib import Path
from typing import Any, Callable, Sequence, Tuple


def run_subprocess_capture(
    cmd: Sequence[str], **kwargs: Any
) -> subprocess.CompletedProcess:
    options = dict(kwargs)
    if "capture_output" not in options and "stdout" not in options:
        options["capture_output"] = True
    if "text" not in options:
        options["text"] = True
    return subprocess.run(list(cmd), **options)


def run_with_fallback(
    primary_cmd: Sequence[str],
    fallback_cmd: Sequence[str],
    run_cmd: Callable[[Sequence[str]], subprocess.CompletedProcess] = run_subprocess_capture,
) -> Tuple[subprocess.CompletedProcess, bool]:
    first = run_cmd(primary_cmd)
    if first.returncode == 0:
        return first, False
    second = run_cmd(fallback_cmd)
    return second, True


def build_openbabel_add_h_cmd(
    obabel_exe: str, input_pdb: str | Path, output_pdb: str | Path
) -> list[str]:
    return [obabel_exe, str(input_pdb), "-O", str(output_pdb), "-h"]


def run_openbabel_add_h(
    input_pdb: str | Path,
    output_pdb: str | Path,
    openbabel_path: str = "",
    *,
    timeout_sec: int | None = None,
) -> None:
    obabel = openbabel_path or shutil.which("obabel")
    if not obabel or not shutil.which(Path(obabel).name):
        raise RuntimeError(
            "Open Babel check failed: executable unresolved "
            f"(OPENBABEL_PATH='{openbabel_path or ''}', PATH probe='{Path(obabel or 'obabel').name}'). "
            "Install Open Babel or set OPENBABEL_PATH, then run `atlas --verify-tools` or `atlas --doctor`."
        )
    cmd = build_openbabel_add_h_cmd(obabel, input_pdb, output_pdb)
    if timeout_sec is None:
        raw_timeout = str(os.environ.get("OPENBABEL_ADDH_TIMEOUT_S", "600") or "600")
        try:
            timeout_val = int(float(raw_timeout))
        except Exception:
            timeout_val = 600
        timeout_sec = max(5, timeout_val)
    else:
        timeout_sec = max(1, int(timeout_sec))
    try:
        r = run_subprocess_capture(cmd, timeout=int(timeout_sec))
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"Open Babel timed out after {int(timeout_sec)}s for '{input_pdb}'"
        ) from exc
    if r.returncode != 0:
        raise RuntimeError(f"Open Babel failed:\n{r.stderr}")


def build_meeko_base_cmd() -> list[str]:
    """
    Build a Meeko CLI command that is portable across machines/environments.
    Prefer invoking by module with the current interpreter to avoid PATH/script shims.
    """
    try:
        import meeko  # type: ignore[import-untyped]  # noqa: F401

        return [sys.executable, "-m", "meeko.cli.mk_prepare_receptor"]
    except Exception:
        pass

    for name in ("mk_prepare_receptor", "mk_prepare_receptor.py"):
        exe = shutil.which(name)
        if exe:
            return [exe]
    try:
        scripts = sysconfig.get_path("scripts")
        cand = os.path.join(scripts, "mk_prepare_receptor.py")
        if os.path.exists(cand):
            return [sys.executable, cand]
    except Exception:
        pass

    raise FileNotFoundError("Meeko CLI not found via module or script.")


def build_meeko_modern_cmd(
    meeko_base_cmd: Sequence[str], input_pdb: str | Path, output_pdbqt: str | Path
) -> list[str]:
    return list(meeko_base_cmd) + ["--read_pdb", str(input_pdb), "-p", str(output_pdbqt)]


def build_meeko_legacy_cmd(
    meeko_base_cmd: Sequence[str], input_pdb: str | Path, output_pdbqt: str | Path
) -> list[str]:
    # Current Meeko releases reject the historical -r/-o form and request
    # --read_pdb plus -p. Keep this helper as a second retry hook, but make
    # the command compatible with the installed open-source CLI.
    return build_meeko_modern_cmd(meeko_base_cmd, input_pdb, output_pdbqt)


def build_adt_prepare_receptor_cmd(
    mgltools_python: str,
    prepare_receptor_script: str,
    input_pdb: str | Path,
    output_pdbqt: str | Path,
) -> list[str]:
    return [
        mgltools_python,
        prepare_receptor_script,
        "-r",
        str(input_pdb),
        "-o",
        str(output_pdbqt),
        "-A",
        "none",
        "-U",
        "nphs_lps",
    ]
