# config.txt

## DUD test-library override (2026-02-14)

- Updated `TEST_LIBRARY_MAP["3EML"]` from `aa2ar` to `test_library_10` for a requested DUD-class validation run.
- Previous library for `3EML` before this change: `aa2ar`.

## Target Install / Massinstall PDB Query Config (2026-03-12)

Preferred user-facing commands:

```bash
atlas targets panels
atlas targets guide
atlas targets genes --genes EGFR ABL1 --gene-type KINASE
atlas targets search "BRCA DNA repair"
atlas targets install --panel kinases
atlas targets install --uniprot P04637
atlas targets install --genes EGFR ABL1 --dry-run
atlas targets install EGFR --ligand ATP --dedupe sequence-identity 90 --quality publication
atlas targets install EGFR --no-dedupe
atlas targets install --genes EGFR --species "Homo sapiens" --taxonomy-id 9606 --methods X-RAY --resolution-max 2.8 --gene-type KINASE
atlas run-panel kinases --ligands chembl --fast
```

These commands wrap the maintained `tools/massinstall.py --pdb_query` machinery,
write selected-target CSVs, download PDB files into the configured input root,
and emit an `atlas_target_install_manifest.json` provenance file. Direct
`python tools/massinstall.py --pdb_query ...` usage remains available for
developer/backward-compatibility workflows, but new user docs should prefer
`atlas targets install`.
Direct UniProt accessions can be supplied with `--uniprot`, and wetlab
co-crystal use cases can require PDB chemical component IDs with `--ligand` or
`--contains-ligand`. Target queries de-duplicate at 90% sequence identity by
default; use `--dedupe sequence-identity 70` for stricter diversity or
`--no-dedupe` to keep near-identical structures. Use `--quality publication` to
exclude reported low-quality models.

For wetlab-first use, `atlas targets guide` prompts for either a built-in panel
or custom gene symbols, writes a reusable gene CSV, then optionally queries RCSB
or downloads the selected PDB files. `atlas targets genes ...` is the
noninteractive way to write only that reusable `gene,category` CSV.

Supported config keys:

- `PDB_QUERY_GENES_FILE`: optional default genes CSV used when `--genes-file` is omitted.
- `PDB_QUERY_DEFAULT_OUT_DIR`: install destination for downloaded PDB files.
- `PDB_QUERY_SELECTED_OUT`: optional CSV output path for selected IDs/metadata.
- `PDB_QUERY_MISSING_OUT`: optional CSV output path for failed downloads.
- `PDB_QUERY_MAX_INSTALL`: optional cap on downloads in `--pdb_query` mode.
- `MASSINSTALL_MAX_INSTALL`: optional cap in legacy static-list mode.
- `MASSINSTALL_MISSING_OUT`: optional missing-download CSV in legacy mode.

`atlas targets install` and `massinstall --pdb_query` use the built-in
`tools/pdb_query.py` defaults for query behavior unless CLI options override
them. They do not read `PDB_QUERY_*` config for query-selection inputs, filters,
or network/cache knobs.

Default target-query filters are human protein structures (`Homo sapiens`,
taxonomy `9606`, entity type `Protein`), all experimental methods,
ligand-bound entries, up to 10 PDBs per gene, and the pdb_query resolution cutoff
for entries where RCSB reports resolution.
Override those per command with `--species`, `--taxonomy-id`, `--entity-type`,
`--methods`, and `--resolution-max`; use `--allow-apo` to retain structures
without a non-trivial ligand. For ad hoc `--genes` installs, `--gene-type` is an
alias for the selected-target CSV category label used by category limits and
manifests.
Common organism aliases infer NCBI taxonomy IDs, so `--species mouse`,
`--organism zebrafish`, or `--species ecoli` are accepted. Run
`atlas targets species` to list supported aliases, and use `--species any` to
disable organism filtering.

Genes source precedence for `atlas targets install`:
1. `--panel`
2. `--genes`
3. `--genes-file`
4. `PDB_QUERY_GENES_FILE` from config
5. `analysis/gene_list/pilot_pdbs_30.csv` if present

Genes source precedence for direct `massinstall --pdb_query`:
1. `--genes-file`
2. `PDB_QUERY_GENES_FILE` from config
3. `analysis/gene_list/pilot_pdbs_30.csv` if present
Default behavior in direct `massinstall --pdb_query` mode is dry-run unless
`--install` is provided or an interactive prompt is accepted. `atlas targets
install` installs by default; pass `--dry-run` to select without downloading.
