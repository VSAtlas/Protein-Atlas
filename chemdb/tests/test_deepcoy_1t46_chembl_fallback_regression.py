from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest
import requests

REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = REPO_ROOT / "config.txt"


def _tail(text: str, limit: int = 2000) -> str:
    return text[-limit:] if text else ""


def _set_config_value(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"^({re.escape(key)}\s*=\s*)([^#\n]*)(.*)$", re.MULTILINE)
    if pattern.search(text):
        return pattern.sub(
            lambda match: f"{match.group(1)}{value}{match.group(3)}",
            text,
            count=1,
        )
    return text.rstrip() + f"\n{key}={value}\n"


def _count_unique_smiles(path: Path) -> int:
    unique = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            smiles = stripped.split()[0]
            if smiles:
                unique.add(smiles)
    return len(unique)


def _extract_path(stdout: str, marker: str) -> Path | None:
    for line in stdout.splitlines():
        if marker in line:
            tail = line.split(marker, 1)[1].strip()
            if not tail:
                continue
            path = Path(tail)
            return path if path.is_absolute() else REPO_ROOT / path
    return None


def _extract_audit_path(stdout: str) -> Path | None:
    for line in stdout.splitlines():
        if "[deepcoy.audit.saved]" not in line:
            continue
        match = re.search(r"path=([^\s]+)", line)
        if not match:
            continue
        path = Path(match.group(1))
        return path if path.is_absolute() else REPO_ROOT / path
    return None


@pytest.mark.network
def test_deepcoy_1t46_chembl_fallback_regression(tmp_path: Path) -> None:  # noqa: ARG001
    if os.environ.get("ATLAS_ALLOW_NETWORK") != "1":
        pytest.skip("set ATLAS_ALLOW_NETWORK=1 to run")

    try:
        resp = requests.get(
            "https://data.rcsb.org/rest/v1/core/entry/1T46", timeout=10
        )
        resp.raise_for_status()
    except Exception as exc:  # pragma: no cover - network guard path
        pytest.skip(f"network preflight failed: {exc}")

    original_config = CONFIG_PATH.read_text(encoding="utf-8")
    updated_config = _set_config_value(
        original_config, "DEEPCOY_DECOYS_PER_ACTIVE", "1"
    )
    CONFIG_PATH.write_text(updated_config, encoding="utf-8")
    try:
        env = os.environ.copy()
        env.setdefault("PYTHONUNBUFFERED", "1")
        cmd = [
            sys.executable,
            "DeepCoy_duds/generate_dud_library.py",
            "--pdb",
            "1T46",
            "--sources",
            "chembl",
            "--no-run-deepcoy",
            "--skip-sdf",
        ]
        try:
            result = subprocess.run(
                cmd,
                cwd=str(REPO_ROOT),
                env=env,
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as exc:
            pytest.fail(
                "generate_dud_library.py failed\n"
                f"returncode: {exc.returncode}\n"
                f"stdout tail:\n{_tail(exc.stdout)}\n"
                f"stderr tail:\n{_tail(exc.stderr)}"
            )

        stdout = result.stdout or ""
        stderr = result.stderr or ""

        if "No actives found; using fallback SMILES" in stdout:
            pytest.fail(
                "Unexpected fallback SMILES path triggered.\n"
                f"stdout tail:\n{_tail(stdout)}\n"
                f"stderr tail:\n{_tail(stderr)}"
            )

        actives_path = _extract_path(stdout, "-> Active SMILES written to: ")
        if not actives_path or not actives_path.exists():
            pytest.fail(
                "Actives SMILES path not found in output.\n"
                f"stdout tail:\n{_tail(stdout)}\n"
                f"stderr tail:\n{_tail(stderr)}"
            )

        unique_smiles = _count_unique_smiles(actives_path)
        if unique_smiles < 3:
            pytest.fail(
                f"Too few unique actives: {unique_smiles} in {actives_path}\n"
                f"stdout tail:\n{_tail(stdout)}\n"
                f"stderr tail:\n{_tail(stderr)}"
            )

        audit_path = _extract_audit_path(stdout)
        if not audit_path or not audit_path.exists():
            audit_path = (
                REPO_ROOT
                / "DeepCoy_duds"
                / "deepcoy_work"
                / "1T46"
                / "cache"
                / "source_audit.json"
            )
        if not audit_path.exists():
            pytest.fail(
                "source_audit.json not found.\n"
                f"stdout tail:\n{_tail(stdout)}\n"
                f"stderr tail:\n{_tail(stderr)}"
            )

        audit_data = json.loads(audit_path.read_text(encoding="utf-8"))
        urls = [
            rec.get("url", "")
            for rec in audit_data.get("requests", [])
            if isinstance(rec, dict)
        ]
        fallback_hit = any("/activity?target_chembl_id=" in url for url in urls) or any(
            "/assay?target_chembl_id=" in url
            and "assay_type=B" in url
            and "relationship_type=D" not in url
            for url in urls
        )
        if not fallback_hit:
            pytest.fail(
                "Fallback ChEMBL URLs missing from audit.\n"
                f"stdout tail:\n{_tail(stdout)}\n"
                f"stderr tail:\n{_tail(stderr)}"
            )
    finally:
        CONFIG_PATH.write_text(original_config, encoding="utf-8")
