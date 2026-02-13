from __future__ import annotations

import builtins
import importlib
import sys
from collections.abc import Iterable

import pytest


def _purge_modules(prefixes: Iterable[str]) -> None:
    for name in list(sys.modules):
        if any(name == p or name.startswith(f"{p}.") for p in prefixes):
            sys.modules.pop(name, None)


def _blocked_import(blocked_roots: set[str]):
    real_import = builtins.__import__

    def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
        root = name.split(".", 1)[0]
        if root in blocked_roots:
            raise ModuleNotFoundError(f"blocked optional dependency: {name}")
        return real_import(name, globals, locals, fromlist, level)

    return guarded_import


def test_core_imports_do_not_require_optional_deps(monkeypatch):
    blocked = {"tensorflow", "pymol", "gtk", "gobject", "DeepCoy", "GGNN_DeepCoy"}
    _purge_modules(blocked)
    _purge_modules(
        [
            "main",
            "docking",
            "prep_ligands",
            "post_docking",
            "cap_driver",
        ]
    )

    monkeypatch.setattr(builtins, "__import__", _blocked_import(blocked))

    modules = [
        "main",
        "docking.docking",
        "docking.docking_subruns",
        "docking.docking_vina_multistage",
        "prep_ligands.prep_ligands",
        "post_docking.rescoring.rescoring_scorch",
    ]
    for mod_name in modules:
        importlib.import_module(mod_name)


def test_optional_deps_loaded_only_on_explicit_paths(monkeypatch):
    _purge_modules({"tensorflow", "pymol", "gtk", "gobject", "DeepCoy", "GGNN_DeepCoy"})
    importlib.import_module("main")
    importlib.import_module("docking.docking")
    importlib.import_module("DeepCoy_duds.utils")

    assert "tensorflow" not in sys.modules
    assert "pymol" not in sys.modules
    assert "gtk" not in sys.modules
    assert "gobject" not in sys.modules

    real_import_module = importlib.import_module

    def _raise_missing(name, *args, **kwargs):
        if name == "pymol" or name.startswith("tensorflow") or name == "gtk":
            raise ModuleNotFoundError(name)
        return real_import_module(name, *args, **kwargs)

    monkeypatch.setattr(importlib, "import_module", _raise_missing)

    import cap_driver
    import DeepCoy_duds.utils as deepcoy_utils
    import molprobity_coot

    with pytest.raises(RuntimeError, match="pymol-open-source"):
        cap_driver._require_pymol_cmd()
    with pytest.raises(RuntimeError, match="conda-forge tensorflow"):
        deepcoy_utils._load_tensorflow()
    with pytest.raises(RuntimeError, match="pygobject"):
        molprobity_coot._require_gtk()
