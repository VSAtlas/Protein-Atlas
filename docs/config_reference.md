# Atlas2 Config Reference

Config values are read from `config.txt` with `config.example.txt` / `config.full.example.txt` as templates. Machine-specific paths should stay in local `config.txt`.

## Paths

| Key | Default | Required | Type | Example | Machine-specific |
| --- | --- | --- | --- | --- | --- |
| `OVERALL_DIR` | `.` | yes | path | `/path/to/protein_automation` | yes |
| `CONFIGS_DIR` | `${OVERALL_DIR}/outputs/configs` | derived | path | `outputs/configs` | no |
| `DOCKED_DIR` | `${OVERALL_DIR}/outputs/docked` | derived | path | `outputs/docked` | no |
| `POST_DOCKED_DIR` | `${OVERALL_DIR}/outputs/post_docked` | derived | path | `outputs/post_docked` | no |
| `OUTPUT_DIR` | `${OVERALL_DIR}/outputs/processed_pdbs/<RUN_ID>` at runtime | derived | path | `outputs/processed_pdbs/20260428_210000` | no |
| `DATA_DIR` | `${OVERALL_DIR}/outputs/data` | derived | path | `outputs/data` | no |
| `MANIFESTS_DIR` | `${OVERALL_DIR}/outputs/manifests` | derived | path | `outputs/manifests` | no |
| `LOGS_DIR` | `${OVERALL_DIR}/outputs/logs` | derived | path | `outputs/logs` | no |
| `INPUT_DIR` | `${OVERALL_DIR}/input_pdbs` | derived | path | `input_pdbs` | no |

## External Tools

| Key | Default | Required | Type | Example | Machine-specific |
| --- | --- | --- | --- | --- | --- |
| `ADFRSUITE_BIN` | blank | optional | path | `/opt/ADFRsuite/bin` | yes |
| `P2RANK_PATH` | `prank` | required for P2Rank mode | path/command | `/opt/p2rank/bin/p2rank.jar` | yes |
| `PHENIX_DIR` | blank | optional | path | `/opt/phenix` | yes |
| `PHENIX_LIB_PATH` | blank | optional | path | `/opt/phenix/lib` | yes |
| `REDUCE_EXE` | PATH lookup | optional | command | `reduce` | yes |
| `PDB2PQR_EXE` | PATH lookup | optional | command | `pdb2pqr` | yes |
| `MEEKO_ALLOW_BAD_RES` | `false` | optional/review-only | boolean | `false` | no |
| `VINA_EXE` | PATH lookup or absolute local path | required for Vina | command/path | `vina` | yes |
| `PYMOL_EXE` | PATH lookup | optional screenshots | command/path | `pymol` | yes |

## Docking

| Key | Default | Required | Type | Example | Machine-specific |
| --- | --- | --- | --- | --- | --- |
| `CPU` | detected/default | yes | integer | `32` | no |
| `DOCKING_MODE` | `polypharmacology` | optional | string | `polypharmacology` | no |
| `FORCE_REPROCESS` | `False` | optional | boolean | `True` | no |
| `ALLOW_BOX_EXPAND` | `true` | optional | boolean | `true` | no |
| `CONTROL_CENTER_POLICY` | config default | optional | string | `best_redock` | no |

## Ligand Prep

| Key | Default | Required | Type | Example | Machine-specific |
| --- | --- | --- | --- | --- | --- |
| `TEST_LIBRARY_MAP` | blank | optional | JSON/dict | `{"TEST": "fda_test_library_10"}` | no |
| `LIBRARY_SUBDIR_DEFAULT` | `fda_library` | optional | string | `fda_library` | no |
| `HMDB_LIBRARY_SUBDIR` | `hmdb` | optional | string | `hmdb` | no |
| `FDA_MAPPING_CSV` | `chemdb/data/fda_mapping_from_pdbqt.csv` | optional | path | `chemdb/data/fda_mapping_from_pdbqt.csv` | no |
| `LIGAND_SOURCE_CACHE_DIR` | `outputs/data/ligand_sources` | optional | path | `outputs/data/ligand_sources` | no |
| `TEST_FDA_LIBRARY_SUBDIR` | `fda_test_library_10` | optional | string | `fda_test_library_10` | no |
| `TEST_FIXTURE_PREPPED_LIGANDS_DIR` | `chemdb/tests/fixtures/prepped_ligands` | optional | path | `chemdb/tests/fixtures/prepped_ligands` | no |
| `TOOL_VERIFY_ON_START` | `false` | optional | boolean | `false` | no |

`chemdb/data/fda_mapping_from_pdbqt.csv` is the canonical runtime FDA mapping. It contains the independently verified terminal-v3 mapping; the immutable repair copy and verifier report remain under `outputs/data/fda_repair_20260712/terminal_reconciliation_v3/` for provenance.

Ligand and receptor PDBQT preparation use Meeko in the open publication stack.

Use `atlas ligands sources` to list built-in source downloads, `atlas ligands
install chembl`, `atlas ligands install chebi`, `atlas ligands install coconut`,
`atlas ligands install hmdb`, or `atlas ligands install fda` to download/cache
source SDFs and prepare PDBQT libraries. `atlas chembl --pdb <ID>`,
`atlas chebi --pdb <ID>`, `atlas coconut --pdb <ID>`, `atlas hmdb --pdb <ID>`,
and `atlas fda --pdb <ID>` provide install-if-needed plus docking in one
command.
If HMDB blocks headless downloads, download the HMDB structures SDF or ZIP in a
browser and run `atlas ligands install hmdb --source-sdf <path>`.

## pH / Protonation

| Key | Default | Required | Type | Example | Machine-specific |
| --- | --- | --- | --- | --- | --- |
| `PH_ENSEMBLE` | `false` | optional | boolean | `true` | no |
| `PH_SCOPE` | blank | optional | string | `active_site` | no |
| `PH_RADIUS` | config default | optional | float | `8.0` | no |
| `PH_LIGAND_MODE` | `off` | optional | string | `context_window` | no |

## Rescoring / MMGBSA

| Key | Default | Required | Type | Example | Machine-specific |
| --- | --- | --- | --- | --- | --- |
| `SCORCH` | PATH/config lookup | optional | path | `/path/to/scorch.py` | yes |
| `SCORCH_DEVICE` | `cpu` | optional | enum | `auto`, `amd`, `nvidia` | no |
| `SCORCH_GPU_IDS` | blank | optional | CSV string | `0,1,2,3` | yes |
| `SCORCH_ENV_PREFIX` | blank | optional | path | `/path/to/micromamba/envs/scorch-env` | yes |
| `SCORCH_ENV` | `scorch-env` | optional | string | `scorch-env` | yes |
| `MMGBSA_ENABLED` | `false` | optional | boolean | `false` | no |
| `MMGBSA_AUTO_RUN` | `false` | optional | boolean | `false` | no |
| `MMGBSA_PROTOCOL` | `screening` | optional | enum | `production` | no |
| `MMGBSA_SELECTED_LIGANDS` | blank | optional | CSV string | `lig1,lig2` | no |
| `MMGBSA_INPUT_PH_LABELS` | blank | optional | CSV string | `pH7_0` | no |
| `MMGBSA_AMBERTOOLS_PREFIX` | blank | optional | path | `/opt/ambertools` | yes |
| `AMBERTOOLS_PREFIX` | blank | optional | path | `/opt/ambertools` | yes |
| `MMGBSA_CHEMISTRY_REVIEW_STATUS` | `unreviewed` | optional | enum | `reviewed` | no |
| `MMGBSA_INPUT_CHEMISTRY_AUTHORITATIVE` | `false` | production | boolean | `true` | no |
| `MMGBSA_INPUT_CHEMISTRY_SOURCE` | `unreviewed_sdf` | production provenance | string | `curated_sdf_ligprep_pH7.4` | no |
| `MMGBSA_PROTONATION_REVIEW_STATUS` | `unreviewed` | production | enum | `reviewed` | no |
| `MMGBSA_TAUTOMER_REVIEW_STATUS` | `unreviewed` | production | enum | `reviewed` | no |
| `MMGBSA_STEREOCHEMISTRY_REVIEW_STATUS` | `unreviewed` | production | enum | `reviewed` | no |
| `MMGBSA_NET_CHARGE_REVIEW_STATUS` | `unreviewed` | production | enum | `reviewed` | no |
| `MMGBSA_PARAMETER_REVIEW_STATUS` | `unreviewed` | production | enum | `reviewed` | no |
| `MMGBSA_MD_MIN_STEPS` | `10000` | optional | integer | `10000` | no |
| `MMGBSA_MD_RESCUE_ENABLED` | `true` | optional | boolean | `true` | no |
| `MMGBSA_MD_RESCUE_STEPS` | `50000` | optional | integer | `50000` | no |

`MMGBSA_AUTO_RUN=false` keeps MM/GBSA uncoupled from normal post-run
finalization. Use `atlas-mmgbsa batch-from-report --report <csv> --run-id <id>`
to run selected report/heatmap hits instead.

`MMGBSA_PROTOCOL=production` expands defaults for explicit-solvent MD
replicates, frame-window analysis, QC artifacts, and stricter checks. Production runs require
`MMGBSA_CHEMISTRY_REVIEW_STATUS=reviewed` after ligand charge/protonation,
APO/HOLO/pH, water, and metal choices have been reviewed. They also require
`MMGBSA_INPUT_CHEMISTRY_AUTHORITATIVE=true` and the protonation, tautomer,
stereochemistry, net-charge, and parameter review subgates to be marked
`reviewed`, `approved`, or `ok`.

Successful Vina PDBQT outputs automatically get a docked-pose SDF when ligand
prep has a `*.ligprep_source.json` sidecar. The SDF uses docked coordinates from
the first PDBQT model but preserves bond orders, charges, stereochemistry, and
protonation from the source SDF record. This is preferred over PDBQT-to-SDF
reconstruction, which remains a smoke-test fallback. Report-driven MM/GBSA also
normalizes restored `outputs/post_docked/<run_id>/...` and `outputs/docked/<run_id>/...`
artifact paths during archive extraction so these SDF sidecars are visible to
selection.

For RESP/QM charges, keep using the existing ligand charge-method setting with
`resp`, `qm_resp`, `qm`, or `precharged`. Atlas searches beside the selected SDF
for `<ligand>.resp.mol2`, `<ligand>.qm_resp.mol2`, `<ligand>.qm.mol2`,
`<ligand>.precharged.mol2`, or matching files under nearby `resp/`, `qm/`, or
`mol2/` folders, then runs `parmchk2` on the copied precharged MOL2. If none is
found, Atlas attempts AmberTools internal RESP fitting and records any fallback
limitations in ligand-prep metadata.

Production runs also write structured per-ligand review records under
`mmgbsa/review_records/` for protonation, tautomer, stereochemistry, net charge,
water/metal policy, APO/HOLO, pH, and input chemistry authority.

Explicit-solvent production now performs a pre-heating minimization QC gate.
Systems with non-finite energy markers, extreme VDW/total energy, high RMS, or
high GMAX are blocked before heating. If enabled, Atlas first runs a staged
restrained/unrestrained minimization rescue and records the rescue/QC result in
`explicit_solvent_methods.json`.

QM/RESP charge fitting is exposed as a script workflow instead of config:
`atlas-mmgbsa resp-plan --ligand <sdf|mol2|pdb> --out-dir <dir> --net-charge <n>`
creates a reviewed AmberTools/Gaussian handoff, and
`atlas-mmgbsa resp-finalize --qm-output <log|out> --out-dir <dir> --net-charge <n>`
creates the RESP MOL2/frcmod after the external QM job completes.

Benchmarking is also script-driven:
`atlas-mmgbsa benchmark-wang2016 --out-dir <dir>` downloads/parses the Wang 2016
supplement where reachable and writes reference target metrics plus optional
Atlas-vs-reference comparisons when `--atlas-results` is provided.
`atlas-mmgbsa wang2016-production-package --report <csv> --run-id <id> --out-dir <dir>`
writes a Slurm production script and Wang ligand review-record scaffolds. The
script intentionally requires chemistry review status environment variables to
be set before submission.

## Report / Heatmap

| Key | Default | Required | Type | Example | Machine-specific |
| --- | --- | --- | --- | --- | --- |
| `REPORT_ASSET_MODE` | `auto` | optional | enum | `inline`, `relative`, `cdn`, `auto` | no |
| `REPORT_HTML_WARN_BYTES` | `10485760` | optional | integer | `10000000` | no |
| `REPORT_TOP_N` | code default | optional | integer | `25` | no |
| `TOP_TARGETS_N` | code default | optional | integer | `3` | no |

## Developer / Debug

| Key | Default | Required | Type | Example | Machine-specific |
| --- | --- | --- | --- | --- | --- |
| `ATLAS_RUN_ID` | timestamp | optional env | string | `pilotstudy` | no |
| `ATLAS_DISTRIBUTED_MODE` | blank | optional env | string | `slurm_array` | no |
| `LIBRARY_MANIFEST_BUILD_ON_SCAN` | code default | optional | boolean | `true` | no |

Use `atlas --doctor` to inspect the active Python/tool environment and `atlas --verify-tools` for detailed external-tool checks.
