# Atlas2 Logging Taxonomy

Atlas2 logs should stay searchable without forcing every module into a new logging framework. New or touched maintained modules should prefer compact structured event fragments in the message text.

## Standard Event Keys

Use these keys when the value is known:

| Key | Meaning |
| --- | --- |
| `component` | Broad subsystem, for example `protein_prep`, `ligand_prep`, `docking`, `rescoring`, or `reporting`. |
| `action` | Specific operation or event, for example `prepare_receptor`, `prepare_ligand`, `pose_validate`, or `report_write`. |
| `run_id` | Run identifier when available. |
| `pdb_id` | PDB or target identifier. |
| `ligand_id` | Ligand display name, stem, or stable ID. |
| `stage` | Pipeline stage such as `clean`, `stage1`, `ctrl_redock`, `scorch`, or `html_report`. |
| `status` | `start`, `ok`, `skip`, `warn`, or `fail`. |
| `reason` | Short machine-searchable reason for `skip`, `warn`, or `fail`. |
| `elapsed_sec` | Runtime in seconds for completed operations. |

## Format

Prefer one bracketed event name followed by `key=value` pairs:

```text
[atlas.event] component=protein_prep action=prepare_receptor run_id=pilotstudy pdb_id=1ABC stage=clean status=ok elapsed_sec=3.42
```

Keep values short and avoid absolute paths unless the log is explicitly a local debug log. When a path matters, include the checked path and a short fix hint in the surrounding error message.

## Examples

Protein prep success:

```text
[atlas.event] component=protein_prep action=prepare_receptor run_id=pilotstudy pdb_id=1ABC stage=receptor_pdbqt status=ok elapsed_sec=4.18
```

Protein prep failure:

```text
[atlas.event] component=protein_prep action=prepare_receptor run_id=pilotstudy pdb_id=1ABC stage=receptor_pdbqt status=fail reason=prepare_receptor_missing
```

Ligand prep malformed ligand:

```text
[atlas.event] component=ligand_prep action=prepare_ligand run_id=pilotstudy ligand_id=LIG_A301 stage=mol2_to_pdbqt status=fail reason=malformed_ligand
```

Docking pose invalid:

```text
[atlas.event] component=docking action=pose_validate run_id=pilotstudy pdb_id=1ABC ligand_id=dexamethasone stage=stage1 status=warn reason=pose_outside_box
```

Report generation complete:

```text
[atlas.event] component=reporting action=report_write run_id=pilotstudy stage=html_report status=ok elapsed_sec=1.27
```

## Initial Maintained-Module Targets

Apply this taxonomy opportunistically in stable maintained modules already touched by quality-gate work:

- `src/cli/run_tool_verification.py`
- `analysis/reporting/run_report_core.py`
- `src/protein_prep/receptor_prep.py`
- `src/prep_ligands/prep_ligands_runtime.py`

Do not mass-edit legacy modules just to rename logs. Prefer taxonomy-aligned messages when a module is already being changed for a real behavior, error-message, or reporting fix.
