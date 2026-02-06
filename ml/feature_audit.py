from __future__ import annotations

import csv
from pathlib import Path


class FeatureAuditWriter:
    def __init__(
        self,
        audit_dir: Path,
        tag: str,
        continuous_feature_names: list[str],
        fp_bits: int,
    ) -> None:
        self.audit_dir = Path(audit_dir)
        self.tag = str(tag).strip() or "data"
        self.continuous_feature_names = list(continuous_feature_names)
        self.fp_bits = int(fp_bits)

        self.audit_dir.mkdir(parents=True, exist_ok=True)

        self.features_path = self.audit_dir / f"features_{self.tag}.csv"
        self.onbits_path = self.audit_dir / f"morgan_onbits_{self.tag}.csv"

        self._features_handle = self.features_path.open(
            "w",
            encoding="utf-8",
            newline="",
        )
        self._onbits_handle = self.onbits_path.open(
            "w",
            encoding="utf-8",
            newline="",
        )

        self._features_writer = csv.writer(self._features_handle)
        self._onbits_writer = csv.writer(self._onbits_handle)

        self._features_writer.writerow(
            [
                "row_number",
                "source_index",
                "active",
                "lig_smiles",
                "ex_rec_pdb",
                "pocket",
                "murcko_scaffold",
                *self.continuous_feature_names,
            ]
        )
        self._onbits_writer.writerow(["row_number", "bit"])

    def write_row(self, meta: dict[str, object], continuous_values: list[float]) -> None:
        row = [
            meta.get("row_number"),
            meta.get("source_index"),
            meta.get("active"),
            meta.get("lig_smiles"),
            meta.get("ex_rec_pdb"),
            meta.get("pocket"),
            meta.get("murcko_scaffold"),
            *continuous_values,
        ]
        self._features_writer.writerow(row)

    def write_onbits(self, row_number: int, onbits: list[int]) -> None:
        for bit in onbits:
            self._onbits_writer.writerow([row_number, int(bit)])

    def close(self) -> None:
        self._features_handle.close()
        self._onbits_handle.close()
