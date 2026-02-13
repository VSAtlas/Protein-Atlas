from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy import sparse

from druggability_orchestrator import load_fpocket_metrics_for_ml
from ml.config import FeaturesConfig
from ml.feature_audit import FeatureAuditWriter
from ml.pocket_features import POCKET_FEATURE_COLUMNS


_FPOCKET_FEATURE_NAMES = (
    "fpocket_druggability",
    "fpocket_volume",
    "fpocket_openness",
    "fpocket_polar_fraction",
    "fpocket_has_metal",
    "fpocket_tier",
)


@dataclass(frozen=True)
class FeaturizedRows:
    X: sparse.csr_matrix
    y: np.ndarray
    source_index: list[int]
    murcko_scaffolds: list[str]
    dropped_invalid_smiles: int


def _safe_float(value: object, *, default: float = 0.0) -> float:
    try:
        parsed = float(value)
        if np.isfinite(parsed):
            return parsed
        return default
    except Exception:
        return default


def _safe_optional_float(value: object) -> float | None:
    try:
        parsed = float(value)
    except Exception:
        return None
    if not np.isfinite(parsed):
        return None
    return float(parsed)


def _optional_text(value: object) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.lower() in {"none", "null", "nan"}:
        return None
    return text


def _descriptor_values(mol: Chem.Mol) -> list[float]:
    return [
        float(Descriptors.MolWt(mol)),
        float(Crippen.MolLogP(mol)),
        float(Lipinski.NumHDonors(mol)),
        float(Lipinski.NumHAcceptors(mol)),
        float(rdMolDescriptors.CalcTPSA(mol)),
        float(Lipinski.NumRotatableBonds(mol)),
        float(rdMolDescriptors.CalcNumRings(mol)),
        float(sum(atom.GetFormalCharge() for atom in mol.GetAtoms())),
    ]


def _extra_descriptor_values(mol: Chem.Mol) -> list[float]:
    hbd = float(Lipinski.NumHDonors(mol))
    hba = float(Lipinski.NumHAcceptors(mol))
    try:
        qed = float(QED.qed(mol))
    except Exception:
        qed = 0.0
    return [
        hbd + hba,
        float(rdMolDescriptors.CalcFractionCSP3(mol)),
        qed,
        float(mol.GetNumHeavyAtoms()),
    ]


def _murcko_scaffold_smiles(mol: Chem.Mol) -> str:
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol) or ""
    except Exception:
        return ""


class BigBindFeaturizer:
    def __init__(
        self,
        *,
        bigbind_root: Path,
        features: FeaturesConfig,
        atlas_cfg: dict[str, Any] | None = None,
        fpocket_center_columns: tuple[str, str, str] = (
            "pocket_center_x",
            "pocket_center_y",
            "pocket_center_z",
        ),
        fpocket_variant_column: str | None = "variant",
        fpocket_ph_column: str | None = "pH",
        fpocket_centers_by_pdb: dict[str, tuple[float, float, float]] | None = None,
        logger: logging.Logger | None = None,
    ):
        self.bigbind_root = bigbind_root
        self.features = features
        self.atlas_cfg = atlas_cfg or {}
        self.fpocket_center_columns = fpocket_center_columns
        self.fpocket_variant_column = fpocket_variant_column
        self.fpocket_ph_column = fpocket_ph_column
        self.fpocket_centers_by_pdb = {
            str(key).strip().lower(): tuple(value)
            for key, value in (fpocket_centers_by_pdb or {}).items()
            if key and value is not None
        }
        self.logger = logger or logging.getLogger("ml.featurize")

        self.fingerprint_bits = max(0, int(features.ligand_morgan_fp_bits))
        self._morgan_generator = (
            rdFingerprintGenerator.GetMorganGenerator(
                radius=2, fpSize=self.fingerprint_bits
            )
            if self.fingerprint_bits > 0
            else None
        )

        self._fpocket_metrics_cache: dict[
            tuple[str, str, str, tuple[float, float, float]], dict[str, float] | None
        ] = {}
        self.loaded_fpocket_count = 0
        self.missing_center_count = 0
        self.missing_info_count = 0

        self._continuous_feature_names: list[str] = []
        if features.ligand_descriptors:
            self._continuous_feature_names.extend(
                [
                    "lig_mw",
                    "lig_logp",
                    "lig_hbd",
                    "lig_hba",
                    "lig_tpsa",
                    "lig_rotatable_bonds",
                    "lig_ring_count",
                    "lig_formal_charge",
                ]
            )
        if features.ligand_extra_descriptors:
            self._continuous_feature_names.extend(
                [
                    "lig_hbd_plus_hba",
                    "lig_fraction_csp3",
                    "lig_qed",
                    "lig_heavy_atom_count",
                ]
            )
        if features.pocket_fpocket:
            self._continuous_feature_names.extend(list(_FPOCKET_FEATURE_NAMES))
        if features.pocket_features:
            self._continuous_feature_names.extend(list(POCKET_FEATURE_COLUMNS))
        if features.vina_score:
            self._continuous_feature_names.append("vina_score_placeholder")

    @property
    def feature_names(self) -> list[str]:
        names = list(self._continuous_feature_names)
        if self.fingerprint_bits > 0:
            names.extend(
                f"morgan_bit_{bit_idx}" for bit_idx in range(self.fingerprint_bits)
            )
        return names

    def metadata(self) -> dict[str, object]:
        return {
            "feature_names": self.feature_names,
            "continuous_feature_names": list(self._continuous_feature_names),
            "fingerprint_bits": self.fingerprint_bits,
            "fingerprint_radius": 2,
            "feature_flags": {
                "ligand_descriptors": self.features.ligand_descriptors,
                "ligand_extra_descriptors": self.features.ligand_extra_descriptors,
                "pocket_fpocket": self.features.pocket_fpocket,
                "vina_score": self.features.vina_score,
            },
            "fpocket_stats": {
                "loaded_fpocket_count": self.loaded_fpocket_count,
                "missing_center_count": self.missing_center_count,
                "missing_info_count": self.missing_info_count,
            },
        }

    def _center_from_row(self, row: pd.Series) -> tuple[float, float, float] | None:
        cx_col, cy_col, cz_col = self.fpocket_center_columns
        cx = _safe_optional_float(row.get(cx_col))
        cy = _safe_optional_float(row.get(cy_col))
        cz = _safe_optional_float(row.get(cz_col))
        if cx is None or cy is None or cz is None:
            return None
        return (cx, cy, cz)

    def _center_for_row(
        self, row: pd.Series, pdb_id: str
    ) -> tuple[float, float, float] | None:
        from_columns = self._center_from_row(row)
        if from_columns is not None:
            return from_columns
        return self.fpocket_centers_by_pdb.get(pdb_id.strip().lower())

    def _row_text(self, row: pd.Series, column: str | None) -> str | None:
        if not column:
            return None
        return _optional_text(row.get(column))

    @staticmethod
    def _fpocket_cache_key(
        pdb_id: str,
        variant: str | None,
        ph_label: str | None,
        center: tuple[float, float, float],
    ) -> tuple[str, str, str, tuple[float, float, float]]:
        return (
            pdb_id.strip().upper(),
            (variant or "").strip(),
            (ph_label or "").strip(),
            tuple(round(coord, 3) for coord in center),
        )

    def _load_fpocket_metrics(
        self,
        *,
        pdb_id: str,
        variant: str | None,
        ph_label: str | None,
        center: tuple[float, float, float],
    ) -> dict[str, float] | None:
        cache_key = self._fpocket_cache_key(pdb_id, variant, ph_label, center)
        if cache_key not in self._fpocket_metrics_cache:
            self._fpocket_metrics_cache[cache_key] = load_fpocket_metrics_for_ml(
                cfg=self.atlas_cfg,
                pdb_id=pdb_id,
                variant=variant,
                ph_label=ph_label,
                center=center,
                logger=self.logger,
            )
        return self._fpocket_metrics_cache[cache_key]

    def _row_fpocket_features(self, row: pd.Series) -> list[float]:
        zeros = [0.0] * len(_FPOCKET_FEATURE_NAMES)
        if all(name in row.index for name in _FPOCKET_FEATURE_NAMES):
            injected = [_safe_optional_float(row.get(name)) for name in _FPOCKET_FEATURE_NAMES]
            non_missing = sum(1 for value in injected if value is not None)
            if non_missing > 0:
                self.loaded_fpocket_count += 1
            else:
                self.missing_info_count += 1
            return [_safe_float(value, default=0.0) for value in injected]

        pdb_id = _optional_text(row.get("ex_rec_pdb")) or _optional_text(row.get("pdb_id"))
        if not pdb_id:
            self.missing_info_count += 1
            return zeros

        variant = self._row_text(row, self.fpocket_variant_column)
        ph_label = self._row_text(row, self.fpocket_ph_column)
        center = self._center_for_row(row, pdb_id)
        if center is None:
            self.missing_center_count += 1
            return zeros

        metrics = self._load_fpocket_metrics(
            pdb_id=pdb_id,
            variant=variant,
            ph_label=ph_label,
            center=center,
        )
        if metrics is None:
            self.missing_info_count += 1
            return zeros

        self.loaded_fpocket_count += 1
        return [_safe_float(metrics.get(name), default=0.0) for name in _FPOCKET_FEATURE_NAMES]

    def _row_continuous_features(self, row: pd.Series, mol: Chem.Mol) -> list[float]:
        values: list[float] = []

        if self.features.ligand_descriptors:
            values.extend(_descriptor_values(mol))

        if self.features.ligand_extra_descriptors:
            values.extend(_extra_descriptor_values(mol))

        if self.features.pocket_fpocket:
            values.extend(self._row_fpocket_features(row))

        if self.features.pocket_features:
            values.extend([_safe_float(row.get(name), default=0.0) for name in POCKET_FEATURE_COLUMNS])

        if self.features.vina_score:
            # Placeholder slot only. Docking/import is intentionally out of scope here.
            values.append(_safe_float(row.get("vina_score"), default=0.0))

        return values

    def transform(
        self,
        df: pd.DataFrame,
        *,
        audit_dir: Path | None = None,
        audit_tag: str = "data",
    ) -> FeaturizedRows:
        if df.empty:
            raise ValueError("Cannot featurize an empty dataframe.")

        cont_rows: list[list[float]] = []
        labels: list[int] = []
        source_index: list[int] = []
        murcko_scaffolds: list[str] = []
        fp_rows: list[int] = []
        fp_cols: list[int] = []
        fp_data: list[float] = []
        dropped_invalid_smiles = 0

        audit_writer: FeatureAuditWriter | None = None
        if audit_dir is not None:
            audit_writer = FeatureAuditWriter(
                audit_dir=Path(audit_dir),
                tag=audit_tag,
                continuous_feature_names=list(self._continuous_feature_names),
                fp_bits=self.fingerprint_bits,
            )

        try:
            for idx, row in df.iterrows():
                smiles = str(row.get("lig_smiles") or "").strip()
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    dropped_invalid_smiles += 1
                    continue

                row_number = len(cont_rows)
                continuous_values = self._row_continuous_features(row, mol)
                scaffold = _murcko_scaffold_smiles(mol)
                cont_rows.append(continuous_values)
                labels.append(int(row.get("active", 0)))
                source_index.append(int(idx))
                murcko_scaffolds.append(scaffold)

                onbits: list[int] = []
                if self.fingerprint_bits > 0:
                    if self._morgan_generator is not None:
                        fp = self._morgan_generator.GetFingerprint(mol)
                    else:
                        fp = AllChem.GetMorganFingerprintAsBitVect(
                            mol,
                            radius=2,
                            nBits=self.fingerprint_bits,
                        )
                    onbits = [int(bit) for bit in fp.GetOnBits()]
                    for bit in onbits:
                        fp_rows.append(row_number)
                        fp_cols.append(bit)
                        fp_data.append(1.0)

                if audit_writer is not None:
                    audit_writer.write_row(
                        {
                            "row_number": row_number,
                            "source_index": int(idx),
                            "active": int(row.get("active", 0)),
                            "lig_smiles": smiles,
                            "ex_rec_pdb": row.get("ex_rec_pdb"),
                            "pocket": row.get("pocket"),
                            "murcko_scaffold": scaffold,
                        },
                        continuous_values,
                    )
                    audit_writer.write_onbits(row_number, onbits)
        finally:
            if audit_writer is not None:
                audit_writer.close()

        if not cont_rows:
            raise ValueError("No valid rows remained after SMILES parsing.")

        n_rows = len(cont_rows)
        cont_cols = len(self._continuous_feature_names)
        if cont_cols > 0:
            cont_array = np.asarray(cont_rows, dtype=np.float32)
            cont_array[~np.isfinite(cont_array)] = 0.0
            cont_matrix = sparse.csr_matrix(cont_array, dtype=np.float32)
        else:
            cont_matrix = sparse.csr_matrix((n_rows, 0), dtype=np.float32)

        if self.fingerprint_bits > 0:
            fp_matrix = sparse.csr_matrix(
                (fp_data, (fp_rows, fp_cols)),
                shape=(n_rows, self.fingerprint_bits),
                dtype=np.float32,
            )
        else:
            fp_matrix = sparse.csr_matrix((n_rows, 0), dtype=np.float32)

        total_feature_count = cont_matrix.shape[1] + fp_matrix.shape[1]
        if total_feature_count == 0:
            raise ValueError("Feature config produced zero columns.")

        if cont_matrix.shape[1] == 0:
            X = fp_matrix
        elif fp_matrix.shape[1] == 0:
            X = cont_matrix
        else:
            X = sparse.hstack([cont_matrix, fp_matrix], format="csr", dtype=np.float32)

        y = np.asarray(labels, dtype=np.int32)
        if audit_dir is not None:
            summary_payload = {
                "dropped_invalid_smiles": int(dropped_invalid_smiles),
                "number_of_rows_written": int(n_rows),
                "fingerprint_bits": int(self.fingerprint_bits),
                "continuous_feature_names": list(self._continuous_feature_names),
            }
            summary_path = Path(audit_dir) / f"feature_audit_{audit_tag}_summary.json"
            summary_path.write_text(json.dumps(summary_payload, indent=2), encoding="utf-8")

        return FeaturizedRows(
            X=X,
            y=y,
            source_index=source_index,
            murcko_scaffolds=murcko_scaffolds,
            dropped_invalid_smiles=dropped_invalid_smiles,
        )
