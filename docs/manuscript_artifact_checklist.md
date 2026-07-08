# Manuscript Artifact Checklist

Use this checklist before submitting an Atlas2 manuscript, preprint, or
reviewer artifact bundle.

## 1. Exact Code Version/Tag

- [ ] Commit hash recorded for every cited run.
- [ ] Release tag or immutable archive created for the manuscript code state.
- [ ] Local diff is empty or attached as a patch artifact.

## 2. Environment File Used

- [ ] Environment file archived: `environment.core.yml`, optional overlay, or
  the internal compatibility `environment.yml`.
- [ ] Active Python executable captured with `atlas --doctor`.
- [ ] Platform and Python version recorded.

## 3. External Tool Versions

- [ ] `atlas --verify-tools` output archived.
- [ ] Vina, Meeko, P2Rank, Reduce/Phenix, receptor-prep tooling, and optional
  conversion/rescoring tools listed when used.
- [ ] Third-party license constraints reviewed.

## 4. Config Snapshot

- [ ] `config.txt` or equivalent config snapshot archived with run artifacts.
- [ ] Machine-specific absolute paths redacted from public report artifacts.
- [ ] Any CLI overrides recorded verbatim.

## 5. Input PDB List

- [ ] Input PDB IDs or local PDB filenames recorded.
- [ ] APO/HOLO mode and water policy recorded.
- [ ] Excluded or failed structures listed with reasons.

## 6. Ligand Library Source and Version

- [ ] Ligand source recorded, including version/date/hash when available.
- [ ] Library preparation command recorded.
- [ ] Any ligand filtering criteria documented.

## 7. Diversity/Background Library Definition

- [ ] Background or diversity library source recorded.
- [ ] Sampling strategy and random seed recorded when applicable.
- [ ] Decoy/control set definitions archived.

## 8. Score Definitions

- [ ] Atlas score / `z_selected` definition cited from report or methods.
- [ ] Percentile-rank reference distribution documented.
- [ ] Coverage, breadth 1%, pose validity, and MM/GBSA usage described when
  present.

## 9. Report Generation Command

- [ ] Report command captured verbatim.
- [ ] Report asset mode recorded (`inline`, `relative`, `cdn`, or `auto`).
- [ ] `run_manifest.yaml`, `report.yaml`, and `heatmap_input.csv` archived.

## 10. Known Limitations

- [ ] Computational docking limitation stated clearly.
- [ ] Rankings described as hypothesis prioritization, not validated binding.
- [ ] Missing external tools, incomplete target coverage, and invalid poses
  reported where relevant.

## 11. Reproducibility Status

- [ ] Fully reproducible from public inputs and open tools.
- [ ] Partially reproducible because private/local data are required.
- [ ] Requires licensed/local external tools.
- [ ] Not included; reason documented.

## 12. Reviewer Quick-start Command

```bash
atlas --doctor
tools/quality_gate.sh
python main.py --test -fast
```
