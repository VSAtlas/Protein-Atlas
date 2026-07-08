# Open Protein Prep Benchmark

`atlas-protein-prep-benchmark` builds an open, reproducible receptor-prep table from RCSB co-crystal structures. It is designed to compare the Atlas open path against the publication-style checks that matter for docking. By default, PDBFixer repairs missing heavy atoms and nonstandard residues but does not build missing residue segments; enable `--add-missing-residues` to test loop-building sensitivity explicitly.

- heavy atom completeness before and after PDBFixer/OpenMM repair
- pH/protonation-stage evidence from PDB2PQR/PROPKA and Reduce/OpenBabel fallback
- metal, cofactor, water, and ligand retention
- HET chemistry assignment viability through CCD/Open Babel plus RDKit sanitization and formal-charge reporting
- HET/cofactor state audit through CCD bond orders, Dimorphite-DL protonation-state enumeration, and RDKit tautomer enumeration
- metal-bound HET classification into competitive ligands, structural cofactors, bridging waters, buffer/salts, and unknown structural ligands
- heuristic metal formal-charge, coordination-number, coordination-geometry, and treatment-route typing
- metal coordination retention by counting nearby N/O/S donors around each retained metal and separating expected competitive-ligand donor loss from unexpected structural loss
- per-chain termini cap plan with explicit ACE/NME detection; audit-only by default
- geometry sanity through nonbonded heavy-atom clash counts and minimum nonbonded distance
- Probe contact validation when `probe` is available
- optional restrained OpenMM protein-only minimization feasibility
- Meeko receptor compatibility
- optional Vina redocking RMSD when `obabel`, `mk_prepare_ligand.py`, and `vina` are available
- water sensitivity through keep-all and remove-all water variants

Run a smoke benchmark:

```bash
atlas-protein-prep-benchmark --limit 3 --out-dir /stor/home/mpg2352/atlas_prep_benchmark_smoke
```

Run the full default 20-target seed set:

```bash
atlas-protein-prep-benchmark --out-dir /stor/home/mpg2352/atlas_prep_benchmark
```

Run optional redocking:

```bash
atlas-protein-prep-benchmark --run-redock --out-dir /stor/home/mpg2352/atlas_prep_benchmark_redock
```

Run strict publication-readiness gates. This enables redocking and requires restrained minimization only for receptors where geometry cleanup actually changed the prepared structure:

```bash
atlas-protein-prep-benchmark --publication-mode --out-dir /stor/home/mpg2352/atlas_prep_benchmark_publication
```

Run the minimization feasibility check when restrained optimization feasibility is needed:

```bash
atlas-protein-prep-benchmark --run-minimize --out-dir /stor/home/mpg2352/atlas_prep_benchmark_minimize
```

Outputs:

- `summary.csv`: one row per target with comparison metrics
- `summary.json`: machine-readable copy of the table
- `<PDB_ID>/`: downloaded input, largest extracted co-crystal ligand, prepared receptor files, water variants, and optional redocking files
- `<PDB_ID>/het_state_audit.json`: per-HET/cofactor state-enumeration details
- `<PDB_ID>/het_state_selection.json`: receptor-contact-scored HET tautomer/protomer/formal-charge selection
- `<PDB_ID>/redock_ligand_prep.json`: extracted-ligand preparation route, including CCD instance-coordinate reconciliation or Open Babel fallback
- `<PDB_ID>/metal_chemistry_audit.json`: heuristic per-metal charge/geometry/treatment-route audit
- `<PDB_ID>/autodock4zn_plan.json`: Zn-specific AD4Zn/TZ map-generation plan, including missing tool checks and runnable command argv
- `<PDB_ID>/metal_parameterization_plan.json`: AmberTools MCPB.py metal-site parameterization handoff templates and required/optional tool checks
- `<PDB_ID>/metal_bound_het_audit.json`: metal-bound HET classification and expected/retained donor-contact accounting
- `<PDB_ID>/termini_audit.json`: per-chain N/C terminus capping plan
- `<PDB_ID>/geometry_clash_audit_pre_fix.json`: exact pre-cleanup close-contact pairs, distances, atom/residue identifiers, atom classes, and ignored/fixable status
- `<PDB_ID>/geometry_cleanup_audit.json`: conservative geometry cleanup actions, including removed waters, dropped repair-added sidechain atoms, and binding-site sidechain rebuild flags
- `<PDB_ID>/geometry_clash_audit.json`: final post-cleanup geometry audit used for the validity gate
- `<PDB_ID>/water_policy_audit.json`: input water retention/removal, occupancy/B-factor metadata, ligand/protein/metal proximity, and bridging-water review flags
- `<PDB_ID>/water_evidence_audit.json`: per-water resolution, occupancy, local protein B-factor, optional CCP4 map sampling, and optional homolog conservation support
- `<PDB_ID>/water_policy_selected.pdb`: conservative receptor variant retaining only supported bridging/metal waters
- `<PDB_ID>/publication_readiness.json`: strict publication-mode gate result, failures, and disclosure-only limitations
- `failures.json`: targets that could not be downloaded or prepared

Validity criteria:

`structure_valid` is true only when the prepared receptor has heavy atoms and hydrogens, Meeko can write a receptor PDBQT, ligand chemistry assignment does not fail when a ligand is present, HET state enumeration does not completely fail for retained HETs, metal count is not reduced relative to the input, no structural metal-bound HET is unexpectedly removed, Probe does not fail when available, and the baseline geometry audit does not find severe nonbonded clashes. `validity_failures` records the exact failed gates.

Publication-readiness criteria:

`publication_ready` is intentionally stricter than `structure_valid`. In `--publication-mode`, Atlas requires the structure-valid gates plus no unresolved Meeko parser-repair fallback, an automatically selected supported-water policy or no input waters, density/conservation/map evidence for any selected retained waters, receptor-context HET state selection, no unresolved metal-physics/protonation review for metal targets, successful seeded redocking with RMSD at or below `--publication-redock-rmsd-max-a` (default 2.5 A), and successful restrained minimization when geometry cleanup changed a receptor. Deterministic Meeko parser normalization that removes waters/hydrogens for receptor PDBQT export is a disclosure, not a hard failure; PDBFixer repair fallback, clash-pruned input, `allow_bad_res`, and legacy parser fallback remain hard review gates. `publication_readiness_failures` is the explicit list of remaining paper-method gaps. `publication_readiness_disclosures` records non-blocking limitations such as local rather than global H-bond optimization, heuristic HET-state selection, heuristic water selection, uncapped termini by default, and optional restrained minimization outside geometry-cleanup cases.

`geometry_clash_audit_pre_fix.json` and `geometry_clash_audit.json` list every heavy-atom close contact below the clash cutoff with atom serial, residue name, chain, residue number, atom name, element, distance, and class flags for water, metal, ligand/cofactor, protein backbone, and protein sidechain. The geometry gate ignores same-chain peptide C-N covalent contacts and chemically plausible metal-donor contacts. Cleanup is intentionally conservative: it removes clashing waters and drops only clashing protein sidechain atoms that were added during repair and were absent from the input PDB. Metals and ligand/cofactor atoms are not removed by this cleanup.

`geometry_cleanup_audit.json` marks dropped repair-added sidechain atoms that are within the binding-site radius of the extracted ligand. Those cases are not publication-complete until the residue is rebuilt or locally minimized successfully; Atlas reports unresolved cases as `binding_site_sidechain_constructive_repair_needed` instead of silently accepting atom deletion as final chemistry. The benchmark now tries a constructive repair before accepting sidechain deletion: it runs a restrained OpenMM/PDBFixer protein repair on the pre-cleanup receptor, restores retained HETs and metal-bound waters, reapplies only conservative water cleanup, and uses the repaired receptor only if the geometry audit and Meeko export pass.

For Meeko receptor export, the benchmark first tries the prepared PDB2PQR/Reduce receptor and ordered-copy retry. If those fail, it tries a Meeko-only normalized input that removes waters/hydrogens, strips Reduce `USER MOD`/`CONECT` records, maps PDB2PQR terminal residue labels such as `NARG`/`CARG` back to standard residue names, and repairs terminal OXT coordinates when the terminal carboxylate geometry would otherwise create false parser bonds. OXT repair writes a sibling `*.oxt_repair_audit.json` and preserves the terminal atom instead of treating parser-driven OXT removal as a clean fix. As parser recovery, Atlas can still produce deterministic sanitized inputs such as salt stripping or clash-pruned copies and rerun Meeko without `allow_bad_res`. Those automatic recovery routes are useful for docking smoke tests, but any Meeko stage with `dropped_atoms>0` is now a publication review gate. The generated sanitized/clash-pruned PDB writes a sibling `*.drop_audit.json` listing the exact dropped atoms and reasons. The `allow_bad_res` path is disabled by default and only runs when `MEEKO_ALLOW_BAD_RES=true`; any such output is reported in `meeko_error`/`meeko_input_stage` and should be treated as a docking-parser rescue, not as solved receptor chemistry.

Production receptor prep writes a `*.meeko_input_stage.json` sidecar next to the receptor PDBQT and in the target `work/` directory. `modern_clean_input` is the direct Meeko path. `normalized_noh_retry` and `normalized_dewatered_noh_retry` are accepted parser-normalization paths because they rerun Meeko after removing receptor hydrogens that can trigger false RDKit valence perception; dewatered export remains a water-policy disclosure. HIS template mapping, atom-dropping sanitized recovery, `allow_bad`, fallback, legacy, existing-receptor, or failed stages remain publication-review cases. The sidecar records both operational `review_required` and stricter `publication_review_required`.

PDB2PQR 3.6.1 does not expose a generic retain-all-HETATM flag in this environment. It supports specific ligand MOL2 handling through `--ligand`/`--reuse-ligand-mol2-files`, but Atlas does not report native PDB2PQR HET retention unless a real generic flag is detected. Metal/cofactor HET preservation after PDB2PQR therefore remains an audited Atlas rescue step, not a solved PDB2PQR-native retention guarantee. Every rescue writes a `*.retained_hets_audit.json` sidecar with source candidates, before/after retained metal/cofactor counts, inserted count, and any retained records still missing after rescue.

`water_policy_audit.json` is an evidence audit, not a conserved-water oracle. A water is flagged as a bridging candidate when it is close to the extracted ligand and a protein polar atom and/or metal. `water_evidence_audit.json` now adds crystallographic support: PDB resolution, water occupancy, water B-factor relative to nearby protein atoms, optional CCP4 `2Fo-Fc`/`Fo-Fc` map sampling when maps are supplied through `ATLAS_WATER_2FOFC_MAP`/`ATLAS_WATER_FOFC_MAP`, and optional homolog conservation using `ATLAS_WATER_CONSERVATION_REFS` or the built-in 3EML/1UYG/1XL2 reference sets. The benchmark also writes `water_policy_selected.pdb`, which keeps only bridging waters without low-occupancy or high-B-factor flags and reruns Meeko on that receptor. `water_redock_sensitivity.json` redocks the same extracted crystal ligand against baseline, selected-water, and baseline-remove-all receptor variants when those PDBQTs are available. In publication mode, a selected-water receptor resolves the water gate only when it passes Meeko, the selected-water redocking RMSD is acceptable and not substantially worse than the remove-all variant, and every selected retained water has density, conservation, or map support. Enrichment evidence is handled by the compact water DUD-E benchmark described below.

`het_state_selection.json` upgrades the older enumeration-only audit. Atlas downloads CCD chemistry, enumerates pH-window protomers with Dimorphite-DL, enumerates tautomers with RDKit, filters invalid states, and scores candidates against an open Protoss-inspired binding-site network optimizer. The optimizer is local and bounded: it builds the HET-connected component within 8 A, includes nearby HIS tautomer/protomer choices, ASN/GLN flips, borderline ASP/GLU/HIS protonation choices, bridging-water orientations, and metal contacts, then scores candidate HET states against that component. The selected state is recorded with candidate count, score, score margin, confidence, formal charge, top alternatives, network mode count, recommended local site modes, and machine-readable review labels such as `multiple_state_candidates`, `low_margin_state_selection`, `metal_adjacent_het_state`, and `open_binding_site_network_optimizer_not_protoss_equivalent`. This is a stronger open heuristic than simple contact counting, but it is still not a full Protoss/Epik-style global coordinate-and-hydrogen-network optimizer. Protoss remains a useful reference method, but Atlas does not vendor it because the available workflow is a ProteinsPlus service rather than a reproducible open-source command-line dependency.

For redocking, Atlas now prepares extracted co-crystal ligands with CCD bond orders and crystal heavy-atom coordinates whenever possible. CCD atom order is reconciled to PDB atom names through the CCD CIF; nucleotide sugar names such as `C5`/`O3` are mapped to PDB prime names such as `C5'`/`O3'` when needed. Open Babel remains a fallback. When the CCD instance route is used, redocking RMSD is computed between the prepared crystal ligand PDBQT and the docked PDBQT so atom-order changes do not inflate RMSD; an RDKit graph/symmetry RMSD is also attempted and the lower valid value is used. Vina redocking is seeded with `--seed 1` for reproducible benchmark gates.

`metal_bound_het_audit.json` determines whether a metal-bound residue should be retained before judging lost contacts. Waters are `bridging_water`, known cofactors are `structural_cofactor`, common salts/ions are `buffer/salt`, the extracted co-crystal ligand or a HET overlapping the docking box is `competitive_ligand`, and metal-bound HETs outside the docking box are `unknown_structural_ligand`. Competitive ligands and salts are expected to be removed for redocking into that site; cofactors, waters, and unknown out-of-box structural ligands are expected to remain. `metal_coord_retained` is policy-aware, while `metal_coord_raw_first_shell_retained` preserves the stricter raw donor-count comparison.

`metal_chemistry_audit.json` is deliberately heuristic. It infers common formal-charge states from the element token, classifies first-shell donor count and approximate geometry from N/O/S donor angles, reports a treatment route, and flags metal-bound histidine/cysteine/carboxylate donors whose protonation state needs curated review. Metal-bound ASP/GLU residues already named as charged carboxylates are now recorded as resolved evidence rather than unresolved protonation failures. `metal_publication_route` separates default docking readiness from force-field claims: retained structural metals with preserved coordination are disclosure-only, while retained competitive ligands, removed ligand-metal donor shells, unresolved protonation, or coordination loss remain hard publication gates. `autodock4zn_candidate` means the site may be suitable for an AutoDock4Zn/TZ map workflow; it does not mean ordinary Vina PDBQT now contains a solved directional zinc force field. Non-Zn metals are reported as preserve/audit/parameterize routes because Atlas has no calibrated generic pseudoatom model for Ca/Mg/Mn/Fe/Cu/Ni/etc.

`autodock4zn_plan.json` follows the open Vina/AutoDock4Zn route: add TZ pseudo-atoms with `zinc_pseudo.py`, generate an AD4Zn GPF with `prepare_gpf4zn.py`, run `autogrid4`, then dock with Vina `--scoring ad4`. The official AutoDock-Vina helper scripts and `AD4Zn.dat` are staged under `tools/autodock4zn/`; `prepare_gpf4zn.py` runs through MGLTools `pythonsh`, and `AD4Zn.dat` is copied into the AutoGrid working directory before execution. The plan is `ready` only when the eligible Zn site, receptor PDBQT, optional ligand PDBQT, and required tools are all present. It intentionally excludes non-Zn metals.

`metal_parameterization_plan.json` stages the open AmberTools/MCPB.py path for non-Zn or non-AD4Zn metal centers. Atlas searches the existing AmberTools prefix, sets `AMBERHOME`, writes an MCPB-only pre-cleaned receptor PDB, remaps `ion_ids`, writes single-metal PDB/MOL2 handoff commands, and records protonation-name suggestions for metal-bound His/Cys/carboxylates. The MCPB staging copy drops hydrogens, preserves terminal `OXT` atoms with the same geometry-repair audit used for Meeko, drops unparameterized non-metal HET residues, and renumbers residues uniquely across chains because MCPB.py indexes model residues by residue number. Metals, waters, and cofactors are now loaded from `aliases.yaml`. Canonical metal-bound cofactors are automatically routed through the existing MMGBSA AmberTools antechamber/parmchk2 fallback when a matching local SDF is already present; successful MOL2/FRCMOD files are staged beside the MCPB input and wired into `naa_mol2files` and `frcmod_files` with MCPB-compatible basenames. For staged cofactor MOL2 files, Atlas now graph-matches the local SDF chemistry against the receptor PDB residue and rewrites atom names to the PDB/CCD residue names when a full mapping is found; heavy atoms are always attempted, and hydrogens/deuteriums are included when the receptor residue contains matching atoms. These edits are only for the MCPB handoff copy, not the production docking receptor. The generated `.mcpb.in` files are heuristic trial inputs, not final chemistry: a publication workflow must still curate metal-site residue selection, oxidation state, structural cofactor/ligand mol2/frcmod files, charge/spin/multiplicity, QM/RESP or empirical stages, and final topology validation.

For MCPB stage 2/3, Atlas records the required QM products and runnable QM command hints when a supported engine is visible. MCPB's default stage 2/3 path consumes Gaussian or GAMESS-US logs (`*_small_fc.log` for Seminario force constants and `*_large_mk.log` for RESP). GAMESS-US is available at no cost but is not open-source/redistributable, so Atlas does not silently install or vendor it. If neither Gaussian nor `rungms` is present, the plan emits `mcpb_qm_required_manual` and keeps stage 2/3 as gated handoff commands rather than pretending the force field is complete. Atlas also writes `ion_info` so MCPB's nonbonded `4n2` route can be executed as an auditable fallback, but this is not equivalent to a bonded MCPB/QM/RESP parameterization.

Other columns are diagnostic rather than hard gates: `termini_status`, `termini_cap_plan`, `restrained_minimization_status`, `metal_coord_*`, `probe_bad_contacts`, water-variant Meeko success, and redocking RMSD.

The default seed set is intentionally RCSB-only rather than redistributing PDBBind or DUD-E assets, and it favors co-crystals with extractable non-water HET ligands. Peptide-only inhibitor complexes are better handled by a separate peptide-docking benchmark because the default ligand extractor intentionally ignores polymer `ATOM` chains. For a formal paper benchmark, replace or extend it with an explicit target list using `--pdb-ids` or `--pdb-id-file`, then archive the exact `summary.csv`, PDB IDs, software versions, config, and output directory manifest.

For stronger reference comparisons:

- Astex Diverse/CCDC validation sets are the most practical prepared-receptor geometry comparator because they include docking-ready protein/ligand files, but download requires accepting CCDC terms.
- DUD-E provides prepared receptors, crystal ligands, actives, and decoys for docking-readiness and enrichment checks.
- PDBBind/CASF and CSAR provide stronger redocking/scoring/protonation benchmarks, but access and redistribution terms need to be documented for any paper artifact.
- PLINDER and PoseX are useful large-scale modern validation options once the smaller prepared-receptor comparisons pass.

Remaining method gaps to track in each benchmark patch:

- HET bond-order and formal-charge assignment now uses CCD/Open Babel/RDKit/Dimorphite-DL plus receptor-contact scoring to select a benchmark HET state. This is still heuristic: Atlas does not yet perform a full coupled ligand/cofactor/water/protein H-bond-network search, and low-confidence selections should be disclosed.
- Termini capping is now planned and audited per chain, but Atlas does not automatically add ACE/NME caps because that changes receptor chemistry and charge.
- Metal handling now audits retention, classifies expected competitive-ligand/salt contact loss, preserves raw first-shell donor counts, heuristically types charge/geometry/treatment route, executes AutoDock4Zn for eligible Zn benchmark sites, emits AMBERHOME-aware MCPB.py heuristic trial inputs for other parameterized sites, pre-cleans MCPB staging copies enough for simple Ca and Mg/ATP stage-1 model generation, attempts aliases.yaml cofactor MOL2/FRCMOD generation when local SDF chemistry is available, stages MCPB-local cofactor files with residue-name and atom-name normalization, records QM log gates for stage 2/3, and distinguishes resolved ASP/GLU carboxylate evidence from unresolved metal-bound His/Cys/Tyr/Lys or protonated-carboxylate cases. Atlas still does not automatically choose the final metal-bound residue protonation state, complete MCPB/QM/RESP parameterization without a supported QM engine and review, fetch or curate missing structural cofactor chemistry automatically, or validate the resulting topology as a production force field. For default docking, structural Ca/Mg-like metals can pass with disclosure when coordination is retained; ligand-metal coordination sites remain blocked unless a metal-aware scoring route is used.
- PDB2PQR/PROPKA, Reduce, and Probe provide pH, flip, hydrogen-network, and contact validation signals, but Atlas does not yet solve a global multi-state H-bond network with ligand/cofactor state penalties.
- Geometry cleanup now audits exact pre/post clash pairs and conservatively removes clashing waters or impossible repair-added sidechain atoms. Benchmark restrained minimization now uses protein heavy atoms, PDBFixer sidechain repair, OpenMM Modeller hydrogens at target pH, and restrained OpenMM minimization. It still does not promote a receptor-wide production minimization stage.
- Meeko-clean receptor export now has normalized and clash-pruned fallback paths for difficult water/metal-heavy receptors; normalized no-hydrogen parser retries are accepted and reported, dewatered normalization is a water-policy disclosure, and any `allow_bad_res`, legacy, open-fallback, existing-output, or failed stage still requires review before publication claims. `allow_bad_res` is off by default and must be explicitly enabled for diagnostic rescue runs.
- Restrained minimization is currently a protein-only OpenMM feasibility check in the benchmark; production promotion should reuse the existing MMGBSA AmberTools/sander minimization machinery rather than adding a second production minimizer.
- Missing residue building is deliberately opt-in because blind loop construction can create docking-parser or chemistry failures; benchmark both settings when missing loops are near the active site.
- Water policy now selects a conservative supported-water receptor, audits density/conservation/map evidence, and redocks baseline/selected/remove-all water variants in publication mode. The optional compact DUD-E water benchmark uses 3EML/AA2AR, 1UYG/HS90A, and 1XL2/HIVPR with <=999 ligands per target. Generate the plan with `python -m protein_prep.water_dude_benchmark plan --out docs/water_dude_benchmark_plan.json`, stage local mini libraries with `python -m protein_prep.water_dude_benchmark stage-libraries`, then run all three policies with `python -m protein_prep.water_dude_benchmark run-policies --run-prefix water_dude_YYYYMMDD --cpu 16 --max-ligands 96`. For individual runs, prefer `atlas analysis dud-eval <RUN_ID>` over direct evaluator module calls. The runner executes `remove_all`, `site_only`, and `keep_all`, writes `docs/water_dude_run_report.json`, and collects metrics in `outputs/data/water_dude_eval/water_policy_summary.tsv`. The water benchmark intentionally truncates docking to the first Vina stage so water-off, selected-water, and keep-all policies are compared on the same compact enrichment pass.
