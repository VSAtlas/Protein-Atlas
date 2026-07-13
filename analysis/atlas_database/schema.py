"""SQLite schema for immutable Atlas release snapshots."""

from __future__ import annotations

import sqlite3

SCHEMA_VERSION = 1

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
    is_control INTEGER NOT NULL DEFAULT 0 CHECK (is_control IN (0, 1)),
    is_decoy INTEGER NOT NULL DEFAULT 0 CHECK (is_decoy IN (0, 1)),
    atlas_score REAL,
    atlas_score_source TEXT,
    selected_docking_score REAL,
    consensus_score REAL,
    final_score REAL,
    final_rank INTEGER,
    source_csv TEXT,
    result_json TEXT,
    UNIQUE (receptor_context_id, ligand_id)
);

CREATE TABLE IF NOT EXISTS docking_attempts (
    attempt_id INTEGER PRIMARY KEY,
    pair_cell_id INTEGER NOT NULL REFERENCES pair_cells(pair_cell_id) ON DELETE CASCADE,
    engine TEXT NOT NULL DEFAULT '',
    stage TEXT NOT NULL DEFAULT '',
    chunk_id TEXT NOT NULL DEFAULT '',
    status TEXT NOT NULL,
    failure_code TEXT,
    failure_reason TEXT,
    completion_path TEXT,
    completion_json TEXT,
    UNIQUE (pair_cell_id, engine, stage, chunk_id)
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
    member_name TEXT,
    sha256 TEXT,
    size_bytes INTEGER,
    file_type TEXT,
    verified INTEGER CHECK (verified IN (0, 1) OR verified IS NULL),
    artifact_json TEXT NOT NULL,
    UNIQUE (run_id, archive_path, member_name)
);

CREATE INDEX IF NOT EXISTS idx_pair_cells_status ON pair_cells(final_status);
CREATE INDEX IF NOT EXISTS idx_pair_cells_score ON pair_cells(atlas_score DESC);
CREATE INDEX IF NOT EXISTS idx_pair_cells_ligand ON pair_cells(ligand_id, atlas_score DESC);
CREATE INDEX IF NOT EXISTS idx_pair_cells_final_score ON pair_cells(final_score DESC);
CREATE INDEX IF NOT EXISTS idx_pair_cells_ligand_final_score
    ON pair_cells(ligand_id, final_score DESC);
CREATE INDEX IF NOT EXISTS idx_context_pdb ON receptor_contexts(pdb_id, variant, ph_label);
CREATE INDEX IF NOT EXISTS idx_attempt_status ON docking_attempts(status);
CREATE INDEX IF NOT EXISTS idx_artifact_pair ON artifacts(pair_cell_id);
"""


def create_schema(connection: sqlite3.Connection) -> None:
    """Create or validate the release database schema."""
    connection.executescript(SCHEMA_SQL)
    current = int(connection.execute("PRAGMA user_version").fetchone()[0])
    if current not in (0, SCHEMA_VERSION):
        raise ValueError(
            f"unsupported Atlas database schema version {current}; expected {SCHEMA_VERSION}"
        )
    connection.execute(f"PRAGMA user_version = {SCHEMA_VERSION}")
