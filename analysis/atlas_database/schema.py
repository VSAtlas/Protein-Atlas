"""SQLite schema for immutable Atlas release snapshots."""

from __future__ import annotations

import sqlite3

SCHEMA_VERSION = 7

SCHEMA_SQL = """
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS releases (
    release_id TEXT PRIMARY KEY,
    schema_version TEXT NOT NULL,
    title TEXT,
    created_at TEXT,
    manifest_sha256 TEXT NOT NULL,
    manifest_json TEXT NOT NULL,
    imported_at TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS scientific_policies (
    release_id TEXT NOT NULL REFERENCES releases(release_id) ON DELETE CASCADE,
    policy_name TEXT NOT NULL,
    policy_json TEXT NOT NULL,
    PRIMARY KEY (release_id, policy_name)
);

CREATE TABLE IF NOT EXISTS runs (
    run_id TEXT PRIMARY KEY,
    release_id TEXT NOT NULL REFERENCES releases(release_id) ON DELETE CASCADE,
    status TEXT,
    manifest_path TEXT NOT NULL,
    manifest_sha256 TEXT NOT NULL,
    command_json TEXT,
    git_json TEXT,
    paths_json TEXT,
    timing_json TEXT,
    resources_json TEXT,
    manifest_json TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS receptor_contexts (
    receptor_context_id INTEGER PRIMARY KEY,
    run_id TEXT NOT NULL REFERENCES runs(run_id) ON DELETE CASCADE,
    pdb_id TEXT NOT NULL,
    variant TEXT NOT NULL DEFAULT '',
    ph_label TEXT NOT NULL DEFAULT '',
    library TEXT,
    status TEXT,
    failure_reason TEXT,
    input_receptor_path TEXT,
    prepared_receptor_path TEXT,
    pocket_method TEXT,
    center_json TEXT,
    box_json TEXT,
    native_redock_status TEXT,
    native_redock_reason TEXT,
    manifest_entry_json TEXT,
    UNIQUE (run_id, pdb_id, variant, ph_label)
);

CREATE TABLE IF NOT EXISTS ligands (
    ligand_id INTEGER PRIMARY KEY,
    canonical_id TEXT NOT NULL UNIQUE,
    display_name TEXT,
    source_path TEXT,
    prepared_state_json TEXT
);

CREATE TABLE IF NOT EXISTS ligand_stereo_evidence (
    ligand_stereo_evidence_id INTEGER PRIMARY KEY,
    ligand_id INTEGER NOT NULL REFERENCES ligands(ligand_id) ON DELETE CASCADE,
    library_name TEXT NOT NULL,
    parent_ligand_id TEXT,
    evidence_schema_version TEXT NOT NULL,
    source_kind TEXT NOT NULL CHECK (source_kind IN ('sdf', 'csv')),
    ligand_source_path TEXT NOT NULL,
    source_file_sha256 TEXT NOT NULL,
    source_record_index INTEGER NOT NULL CHECK (source_record_index >= 1),
    source_record_id TEXT NOT NULL,
    source_record_sha256 TEXT NOT NULL,
    audited_structure_sha256 TEXT NOT NULL,
    structure_column TEXT,
    id_column TEXT,
    structure_format TEXT NOT NULL,
    prepared_artifact_path TEXT,
    prepared_artifact_declared_sha256 TEXT,
    prepared_artifact_computed_sha256 TEXT,
    prepared_artifact_verification_status TEXT NOT NULL,
    source_record_evidence_json TEXT NOT NULL,
    parse_status TEXT NOT NULL CHECK (parse_status IN ('parsed', 'parse_failed')),
    audit_status TEXT NOT NULL CHECK (
        audit_status IN (
            'no_potential_stereo', 'fully_specified',
            'partially_unspecified', 'all_unspecified', 'parse_failed'
        )
    ),
    potential_stereo_count INTEGER NOT NULL CHECK (potential_stereo_count >= 0),
    specified_stereo_count INTEGER NOT NULL CHECK (specified_stereo_count >= 0),
    unspecified_stereo_count INTEGER NOT NULL CHECK (unspecified_stereo_count >= 0),
    unknown_stereo_count INTEGER NOT NULL CHECK (unknown_stereo_count >= 0),
    tetrahedral_count INTEGER NOT NULL CHECK (tetrahedral_count >= 0),
    double_bond_count INTEGER NOT NULL CHECK (double_bond_count >= 0),
    other_stereo_count INTEGER NOT NULL CHECK (other_stereo_count >= 0),
    has_unspecified_potential_stereo INTEGER NOT NULL
        CHECK (has_unspecified_potential_stereo IN (0, 1)),
    has_unresolved_potential_stereo INTEGER NOT NULL
        CHECK (has_unresolved_potential_stereo IN (0, 1)),
    canonical_isomeric_smiles TEXT,
    canonical_smiles TEXT,
    inchikey TEXT,
    stereo_elements_json TEXT NOT NULL,
    rdkit_version TEXT NOT NULL,
    method TEXT NOT NULL,
    error TEXT,
    selected_for_release INTEGER NOT NULL DEFAULT 0
        CHECK (selected_for_release IN (0, 1)),
    selection_method TEXT,
    annotation_source_path TEXT NOT NULL,
    annotation_source_sha256 TEXT NOT NULL,
    annotation_record_index INTEGER NOT NULL CHECK (annotation_record_index >= 1),
    source_record_json TEXT NOT NULL,
    UNIQUE (
        library_name, ligand_id, source_file_sha256, source_record_index,
        source_record_sha256, audited_structure_sha256
    )
);

CREATE TABLE IF NOT EXISTS pair_cells (
    pair_cell_id INTEGER PRIMARY KEY,
    receptor_context_id INTEGER NOT NULL REFERENCES receptor_contexts(receptor_context_id) ON DELETE CASCADE,
    ligand_id INTEGER NOT NULL REFERENCES ligands(ligand_id),
    final_status TEXT NOT NULL,
    failure_code TEXT,
    failure_reason TEXT,
    expected INTEGER NOT NULL DEFAULT 0 CHECK (expected IN (0, 1)),
    has_result INTEGER NOT NULL DEFAULT 0 CHECK (has_result IN (0, 1)),
    pose_valid INTEGER CHECK (pose_valid IN (0, 1) OR pose_valid IS NULL),
    pose_validation_method TEXT,
    pose_validation_scope TEXT,
    pose_validation_thresholds_json TEXT,
    is_control INTEGER NOT NULL DEFAULT 0 CHECK (is_control IN (0, 1)),
    is_decoy INTEGER NOT NULL DEFAULT 0 CHECK (is_decoy IN (0, 1)),
    atlas_score REAL,
    atlas_score_source TEXT,
    selected_docking_score REAL,
    consensus_score REAL,
    final_score REAL,
    final_score_source TEXT,
    final_score_source_reconstructed TEXT,
    final_score_source_classification TEXT,
    final_score_source_evidence_json TEXT,
    final_rank INTEGER,
    source_csv TEXT,
    result_json TEXT,
    UNIQUE (receptor_context_id, ligand_id)
);

CREATE TABLE IF NOT EXISTS completion_records (
    completion_record_id INTEGER PRIMARY KEY,
    run_id TEXT NOT NULL REFERENCES runs(run_id) ON DELETE CASCADE,
    receptor_context_id INTEGER NOT NULL
        REFERENCES receptor_contexts(receptor_context_id) ON DELETE CASCADE,
    engine TEXT NOT NULL DEFAULT '',
    stage TEXT NOT NULL DEFAULT '',
    chunk_id TEXT NOT NULL DEFAULT '',
    completion_path TEXT NOT NULL,
    completion_sha256 TEXT NOT NULL,
    completion_json TEXT NOT NULL,
    UNIQUE (run_id, completion_path)
);

CREATE TABLE IF NOT EXISTS docking_attempts (
    attempt_id INTEGER PRIMARY KEY,
    pair_cell_id INTEGER NOT NULL REFERENCES pair_cells(pair_cell_id) ON DELETE CASCADE,
    completion_record_id INTEGER NOT NULL
        REFERENCES completion_records(completion_record_id) ON DELETE CASCADE,
    engine TEXT NOT NULL DEFAULT '',
    stage TEXT NOT NULL DEFAULT '',
    chunk_id TEXT NOT NULL DEFAULT '',
    status TEXT NOT NULL,
    failure_code TEXT,
    failure_reason TEXT,
    selected_for_release INTEGER NOT NULL DEFAULT 0
        CHECK (selected_for_release IN (0, 1)),
    selection_method TEXT,
    selection_manifest_index INTEGER,
    UNIQUE (pair_cell_id, completion_record_id, engine, stage, chunk_id)
);


CREATE TABLE IF NOT EXISTS result_attempts (
    result_attempt_id INTEGER PRIMARY KEY,
    pair_cell_id INTEGER NOT NULL REFERENCES pair_cells(pair_cell_id) ON DELETE CASCADE,
    completion_record_id INTEGER
        REFERENCES completion_records(completion_record_id) ON DELETE SET NULL,
    completion_link_method TEXT,
    completion_link_evidence_json TEXT,
    input_csv_path TEXT NOT NULL,
    input_csv_sha256 TEXT NOT NULL,
    source_row_number INTEGER NOT NULL,
    result_sha256 TEXT NOT NULL,
    final_status TEXT NOT NULL,
    failure_code TEXT,
    failure_reason TEXT,
    pose_valid INTEGER CHECK (pose_valid IN (0, 1) OR pose_valid IS NULL),
    pose_validation_method TEXT,
    pose_validation_scope TEXT,
    pose_validation_thresholds_json TEXT,
    is_control INTEGER NOT NULL DEFAULT 0 CHECK (is_control IN (0, 1)),
    is_decoy INTEGER NOT NULL DEFAULT 0 CHECK (is_decoy IN (0, 1)),
    atlas_score REAL,
    atlas_score_source TEXT,
    selected_docking_score REAL,
    consensus_score REAL,
    final_score REAL,
    final_score_source TEXT,
    final_score_source_reconstructed TEXT,
    final_score_source_classification TEXT,
    final_score_source_evidence_json TEXT,
    final_rank INTEGER,
    source_csv TEXT,
    result_json TEXT,
    selected_for_release INTEGER NOT NULL DEFAULT 0
        CHECK (selected_for_release IN (0, 1)),
    selection_method TEXT,
    selection_manifest_index INTEGER,
    UNIQUE (pair_cell_id, input_csv_sha256, source_row_number)
);
CREATE TABLE IF NOT EXISTS artifacts (
    artifact_id INTEGER PRIMARY KEY,
    run_id TEXT NOT NULL REFERENCES runs(run_id) ON DELETE CASCADE,
    receptor_context_id INTEGER REFERENCES receptor_contexts(receptor_context_id) ON DELETE SET NULL,
    pair_cell_id INTEGER REFERENCES pair_cells(pair_cell_id) ON DELETE SET NULL,
    stage TEXT,
    mode TEXT,
    original_path TEXT,
    archive_path TEXT,
    ligand_id INTEGER REFERENCES ligands(ligand_id) ON DELETE SET NULL,
    member_name TEXT,
    artifact_role TEXT NOT NULL,
    artifact_scope TEXT NOT NULL
        CHECK (artifact_scope IN ('run', 'receptor', 'ligand', 'pair')),
    sha256 TEXT,
    size_bytes INTEGER,
    file_type TEXT,
    verified INTEGER CHECK (verified IN (0, 1) OR verified IS NULL),
    artifact_json TEXT NOT NULL,
    UNIQUE (run_id, archive_path, member_name)
);

CREATE TABLE IF NOT EXISTS protein_identities (
    protein_identity_id INTEGER PRIMARY KEY,
    protein_key TEXT NOT NULL UNIQUE,
    uniprot_id TEXT,
    gene_symbol TEXT,
    display_name TEXT,
    source_path TEXT NOT NULL,
    source_sha256 TEXT NOT NULL,
    source_record_index INTEGER NOT NULL,
    source_record_json TEXT NOT NULL,
    provenance_json TEXT
);

CREATE TABLE IF NOT EXISTS receptor_annotations (
    receptor_annotation_id INTEGER PRIMARY KEY,
    receptor_context_id INTEGER NOT NULL UNIQUE
        REFERENCES receptor_contexts(receptor_context_id) ON DELETE CASCADE,
    protein_identity_id INTEGER
        REFERENCES protein_identities(protein_identity_id),
    receptor_classification TEXT,
    classification_method TEXT,
    chemistry_evidence_status TEXT,
    prepared_receptor_sha256 TEXT,
    canonical_chemistry_policy_sha256 TEXT,
    retained_metal_atom_count INTEGER
        CHECK (retained_metal_atom_count >= 0 OR retained_metal_atom_count IS NULL),
    retained_cofactor_residue_count INTEGER
        CHECK (
            retained_cofactor_residue_count >= 0
            OR retained_cofactor_residue_count IS NULL
        ),
    requested_observed_conflict INTEGER
        CHECK (requested_observed_conflict IN (0, 1)
               OR requested_observed_conflict IS NULL),
    qualification_status TEXT,
    qualification_reason TEXT,
    native_redock_reason TEXT,
    native_redock_status TEXT,
    native_redock_rmsd REAL,
    source_path TEXT NOT NULL,
    source_sha256 TEXT NOT NULL,
    source_record_index INTEGER NOT NULL,
    source_record_json TEXT NOT NULL,
    provenance_json TEXT
);

CREATE TABLE IF NOT EXISTS receptor_audits (
    receptor_audit_id INTEGER PRIMARY KEY,
    receptor_context_id INTEGER NOT NULL
        REFERENCES receptor_contexts(receptor_context_id) ON DELETE CASCADE,
    audit_kind TEXT NOT NULL,
    source_path TEXT NOT NULL,
    source_sha256 TEXT NOT NULL,
    parse_status TEXT NOT NULL,
    audit_json TEXT NOT NULL,
    UNIQUE (receptor_context_id, audit_kind, source_sha256)
);

CREATE TABLE IF NOT EXISTS known_pair_selections (
    known_pair_selection_id INTEGER PRIMARY KEY,
    receptor_context_id INTEGER NOT NULL UNIQUE
        REFERENCES receptor_contexts(receptor_context_id) ON DELETE CASCADE,
    ligand_canonical_id TEXT NOT NULL,
    pair_cell_id INTEGER REFERENCES pair_cells(pair_cell_id) ON DELETE SET NULL,
    selection_label TEXT,
    evidence_reference TEXT,
    source_path TEXT NOT NULL,
    source_sha256 TEXT NOT NULL,
    source_record_index INTEGER NOT NULL,
    source_record_json TEXT NOT NULL,
    provenance_json TEXT
);

CREATE INDEX IF NOT EXISTS idx_ligand_stereo_ligand
    ON ligand_stereo_evidence(ligand_id, library_name);
CREATE UNIQUE INDEX IF NOT EXISTS idx_ligand_stereo_selected
    ON ligand_stereo_evidence(library_name, ligand_id)
    WHERE selected_for_release = 1;
CREATE INDEX IF NOT EXISTS idx_pair_cells_status ON pair_cells(final_status);
CREATE INDEX IF NOT EXISTS idx_pair_cells_score ON pair_cells(atlas_score DESC);
CREATE INDEX IF NOT EXISTS idx_pair_cells_ligand ON pair_cells(ligand_id, atlas_score DESC);
CREATE INDEX IF NOT EXISTS idx_pair_cells_final_score ON pair_cells(final_score DESC);
CREATE INDEX IF NOT EXISTS idx_pair_cells_ligand_final_score
    ON pair_cells(ligand_id, final_score DESC);
CREATE INDEX IF NOT EXISTS idx_context_pdb ON receptor_contexts(pdb_id, variant, ph_label);
CREATE INDEX IF NOT EXISTS idx_attempt_status ON docking_attempts(status);
CREATE INDEX IF NOT EXISTS idx_attempt_completion
    ON docking_attempts(completion_record_id);
CREATE INDEX IF NOT EXISTS idx_attempt_selected
    ON docking_attempts(pair_cell_id, selected_for_release);
CREATE INDEX IF NOT EXISTS idx_result_attempt_pair
    ON result_attempts(pair_cell_id, selected_for_release);
CREATE INDEX IF NOT EXISTS idx_result_attempt_completion
    ON result_attempts(completion_record_id, selected_for_release);
CREATE INDEX IF NOT EXISTS idx_completion_context
    ON completion_records(receptor_context_id);
CREATE INDEX IF NOT EXISTS idx_artifact_pair ON artifacts(pair_cell_id);
CREATE INDEX IF NOT EXISTS idx_artifact_ligand ON artifacts(ligand_id);
CREATE INDEX IF NOT EXISTS idx_artifact_role ON artifacts(artifact_role);
CREATE INDEX IF NOT EXISTS idx_protein_identity_uniprot
    ON protein_identities(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_protein_identity_gene
    ON protein_identities(gene_symbol);
CREATE INDEX IF NOT EXISTS idx_receptor_annotation_protein
    ON receptor_annotations(protein_identity_id);
CREATE INDEX IF NOT EXISTS idx_receptor_audit_context
    ON receptor_audits(receptor_context_id, audit_kind);
CREATE INDEX IF NOT EXISTS idx_known_pair_ligand
    ON known_pair_selections(ligand_canonical_id);
"""


def create_schema(connection: sqlite3.Connection) -> None:
    """Create or validate the release database schema."""
    current = int(connection.execute("PRAGMA user_version").fetchone()[0])
    if current not in (0, SCHEMA_VERSION):
        raise ValueError(
            f"unsupported Atlas database schema version {current}; expected {SCHEMA_VERSION}"
        )
    connection.executescript(SCHEMA_SQL)
    connection.execute(f"PRAGMA user_version = {SCHEMA_VERSION}")
