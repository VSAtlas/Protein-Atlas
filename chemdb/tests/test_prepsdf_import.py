import importlib
import sys


def test_prepsdf_import_no_config_load(monkeypatch):
    module_name = "prep_ligands.prepsdf"
    sys.modules.pop(module_name, None)

    def _raise_if_called(*args, **kwargs):
        raise AssertionError("load_config should not run during import")

    monkeypatch.setattr("input_and_export_functions.load_config", _raise_if_called)

    module = importlib.import_module(module_name)
    assert hasattr(module, "decompress_sdf_gz_files")

    sys.modules.pop(module_name, None)
