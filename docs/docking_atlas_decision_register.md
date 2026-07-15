# Docking Atlas decision register

This register separates implemented engineering contracts from scientific choices
that require explicit user approval. Atlas must not infer or auto-fill a pending
scientific decision. Record the reviewer, date, evidence snapshot, and approved
value in the release manifest and/or its hashed external annotation sources before
freezing a public release.

## Already locked

- The public resource is read-only and release-versioned.
- Protein/receptor exploration and drug exploration are bidirectional; receptor
  contexts and pH states remain separate.
- `final_score` is the primary score, higher is better, shared ties are retained,
  and no secondary score silently replaces a missing primary score. A normalized
  value is comparable only when its explicit score-source family is shown, such as
  consensus-versus-decoy Z or SCORCH-versus-decoy Z; an unclassified legacy value
  must remain provenance-pending.
- Invalid and unsuccessful cells remain present, but are excluded from headline
  rankings by default.
- APO contexts remain visible but are excluded from the qualified HOLO headline
  ranking. Ranking uses the controlled observed `receptor_classification`, not the
  requested variant/folder token. A separate, explicitly exploratory APO ranking
  is approved; APO and HOLO values must never be silently pooled into one list.
- Receptor annotations accept only controlled `APO`/`HOLO` classifications or
  an unresolved blank. Readiness schema 6 reports missing/uncontrolled
  classification and exact prepared-receptor evidence gaps rather than treating a
  historical run label as chemistry truth.
- A receptor must pass an explicit receptor-quality outcome and an explicit native
  redock qualification before its drug rows can enter headline rankings.
- The native-redock qualification threshold is 2.5 angstrom on an independently
  score-verified top-ranked generated pose. The current evaluator verifies ranking
  only from exact per-model Vina `REMARK VINA RESULT` records; missing competitor
  scores or unsupported score sources remain evidence-pending and cannot qualify.
  RMSD uses exact heavy-atom molecular-graph identity, symmetry-equivalent atom
  mappings, and a receptor transform independently checked against exact common
  receptor atom keys; it must not truncate atom lists, assume file order, include
  hydrogens, or independently superpose the docked ligand. Best-of-generated-pose
  RMSD is diagnostic and cannot substitute for top-ranked qualification.
- The conference image set contains the native control, up to five top valid
  ligands, the best-scoring invalid result when present, and a representative
  known-pair slot when it is not already selected. The slot is a navigation-only,
  directly curated drug-target annotation with an explicit evidence citation; it
  never changes ranking or validates a docking result.
- For frozen public releases, protein identities, receptor outcomes, and known pairs
  come from explicit annotation rows backed by hashed external sources. Older inline
  known-pair values may remain covered by a manifest hash, but `atlas publish build`
  ignores them; only the advanced callable's explicit parameter remains
  compatibility-only.

## Current evidence snapshot (2026-07-13)

The engineering-only matrix audit of
`spd90_fda_fda_dud_spr10_40h_patchdist_20260525_161503` contains 90 receptor
contexts, 760,878 terminal expected cells, 440,183 calculated-but-unvalidated
cells, 320,695 structured missing-result cells, and 31,950 non-null
`final_score` values. It has zero explicit pair pose validations, zero explicit
native-redock qualifications, and no imported artifact index. The source rows had
no explicit controlled `final_score` source; v0.0.12 now reconstructs all 31,950
finite values from exact row-bound evidence, without resolving their immutable
pose identity. It is not approved as a scientific release.

A non-overwriting full-build rehearsal on 2026-07-13 failed closed during input
audit because the historical run contains duplicate successful completion
attempts without explicit manifest attempt selections. No retry was selected
implicitly. The direct production-equivalent source classifier independently
verified all 31,950 finite scores against input CSV SHA-256
`0dbe4ca37b21d5b61ea5dc4c640a4ea81b60284d676a9d8d4cf5422f365cc742`.
Publication still requires an approved, evidence-bound attempt selection for
every duplicate before rebuilding the frozen release.

The v2 prepared-receptor evidence pass pins those exact 90 context keys with
inventory SHA-256
`f42acf852c878e092cffd6b93f802e405db6f44637478a0a97810811f48c9579`.
Its annotation SHA-256 is
`5f4f115ca7439916dc6b237098b7cb429441d99c64bb44f3627662534cc8298a`.
Exact final prepared chemistry is observed for 75 contexts: 71 classify as APO
and 4 as HOLO under the evidence-only vocabulary. Fifteen remain unresolved,
including cleaned-PDB/final-PDBQT disagreements. All 107 discovered legacy
retained-chemistry audits were rejected for exact receptor-path mismatch; those
rejections are retained as evidence, but no audit is cross-attached or promoted
to qualification.

The canonical status snapshot for `postFDAfixSPDrun_prepcheck_20260713` reported
92 of 93 receptor contexts complete, one pending context, and two structured
control-redock `fallback_missing` examples. Refresh canonical status before review;
do not convert those observations into pass/fail outcomes automatically.

## Scientific constraints approved on 2026-07-13

- The conference inventory is the historical 90-context set identified above.
  Both the count and exact run/PDB/variant/pH key-set hash are frozen; a manifest
  change must produce and explicitly approve a different inventory hash.
- APO/HOLO must describe retained chemistry in the exact prepared receptor content,
  not a folder name, requested run branch, or historical `HOLO` label. The exact
  cleaned preparation and exact docking PDBQT must be separately hashed and any
  disagreement surfaced.
- Metal-retention and metal-site interaction audits are included as receptor
  evidence and must be bound to the exact receptor context. Their presence does not
  by itself qualify receptor quality.
- The approved native-redock cutoff is 2.5 angstrom. The implemented strict v2
  symmetry-aware heavy-atom evaluator can qualify only when every receptor, ligand,
  coordinate-frame, extraction, all-model score, and selected-pose input verifies.
  Unspecified potential stereochemistry is explicitly policy-pending, not invalid;
  it cannot qualify until the user approves a policy. Unprovenanced legacy calls
  remain evidence-pending and cannot qualify.
- Pose-validity thresholds come from the existing validation script. The public
  database must show the validator version/configuration and thresholds rather than
  silently reinterpreting them.
- Normalized primary scores may be compared only when the score-source family is
  explicit and displayed alongside the value.
- The conference image plan includes the native control, up to five top valid
  ligands, the best-scoring invalid result when present, and a representative known
  pair when it is not already among those selections.
- The approved public-artifact baseline is the prepared receptor, prepared ligand,
  and docked-ligand images. The selected docking pose, native ligand and redock pose,
  validation reports, search-box definition, sanitized run configuration, release
  manifest, and content hashes are also approved. Raw stdout/stderr, unselected
  poses, all-pose archives, and unsanitized intermediates remain private by default.
- Pose validity must be recomputed or verified for the exact pose that supplies the
  displayed `final_score`. An any-stage or any-pose pass is retained as legacy
  evidence but cannot qualify a ranked pair. A backfill must record the exact
  validator version, configuration, thresholds, selected-pose hash, and output hash.
- Exact final prepared receptor content is authoritative for what entered docking.
  A cleaned-versus-final chemistry disagreement remains unresolved for release
  qualification until audited. Site-bound canonical cofactors and nucleotides may
  establish HOLO chemistry; free halides and nonspecific ions alone may not.

## Native-redock method basis

The v2 implementation uses the fixed-receptor-frame docking RMSD definition and
graph-isomorphism symmetry rationale described by Bell and Zhang's
[DockRMSD](https://doi.org/10.1186/s13321-019-0362-7); it does not claim to
implement DockRMSD itself. It uses RDKit's documented in-place
[`CalcRMS`](https://www.rdkit.org/docs/source/rdkit.Chem.rdMolAlign.html) with
explicit graph maps and reconstructs pose coordinates from Meeko's documented
[SMILES and index mapping](https://meeko.readthedocs.io/en/develop/export_usage.html).

Before RMSD, every exactly keyed common receptor `ATOM` coordinate must agree
with the declared rigid transform within 0.002 angstrom. That tolerance records
three-decimal PDB/PDBQT coordinate-serialization evidence only; it is not a
receptor-quality criterion or scientific cutoff. Top rank is independently
verified only from exact all-model Vina `REMARK VINA RESULT` records. Unsupported
score sources, missing competitor scores, and unresolved transforms remain
evidence-pending. Policies for unspecified potential stereochemistry and
resonance-equivalent terminal conjugated groups remain pending user approval.

## Fail-closed evidence behavior

The evidence extractor may report observations, hashes, raw retained-chemistry
categories, and audit payloads, but it must not turn them into receptor-quality or
native-redock qualification:

- empty, malformed, unsupported, or protein-coordinate-free receptor content is
  unresolved, never APO;
- cleaned-preparation versus exact-docking-receptor disagreement is explicit and
  leaves the combined classification unresolved;
- salts/free ions, halides, nucleotide-like residues, cofactors, and metals remain
  separately visible so an unapproved conflict rule cannot be applied silently;
- requested run labels and observed chemistry conflicts are recorded, while stale
  or differently scoped metal audits are rejected instead of cross-attached; and
- qualification remains pending until explicit quality and redocking decisions are
  imported from hashed evidence.

## Scientific approvals still required

| Decision | Already constrained | Approval still needed |
|---|---|---|
| Exact 90-context inventory | Historical 90-context key set is frozen at inventory SHA-256 `f42acf852c878e092cffd6b93f802e405db6f44637478a0a97810811f48c9579`. | No further choice for v0.1 unless the manifest changes; any change requires an explicit new inventory approval. |
| Experimental receptor-quality policy | Evidence must remain separate from qualification. | Accepted experimental methods; required source fields; metric definitions and cutoffs for resolution/quality, mutations, missing residues, ligand/chain selection, and pocket completeness; policy identifier; per-context outcome and reason. |
| APO/HOLO conflict rules and ranking use | Exact final prepared content is authoritative; cleaned/final disagreement remains unresolved; canonical site-bound cofactors/nucleotides can establish HOLO; free halides or nonspecific ions alone cannot; APO receives a separate exploratory ranking. | Explicit review of the 15 unresolved contexts and evidence-bound per-context approval. |
| Native-redock RMSD method | The approved cutoff is an inclusive 2.5 angstrom. The strict v2 implementation uses Meeko SMILES-IDX mapping, exact graph/bond identity, heavy atoms, symmetry mappings, in-place RDKit CalcRMS, exact-common-atom receptor-transform verification, and exact all-model Vina score-order verification. Best-of-generated remains diagnostic. Unverified transforms, unsupported score sources, missing competitor scores, and unprovenanced calls remain evidence-pending and cannot qualify. | Approve how unspecified potential stereochemistry should be handled and whether resonance-equivalent terminal conjugated groups should be interchangeable. Generate the required manifests and per-context outcomes only after those choices; no receptor is qualified by the new method yet. |
| Pair pose-validity scope | The exact selected `final_score` pose must pass the existing script predicates with visible validator configuration and thresholds; legacy any-stage evidence cannot qualify. | Backfill exact selected-pose validation for the frozen run. The 31,950 classified historical scores identify SCORCH-stage lineage but no immutable pose, so the contributing-pose/display-pose policy must be approved before this backfill can qualify rows. |
| Legacy `final_score` provenance | All 31,950 finite frozen-run scores are immutably reconstructed as `legacy_scorch_percentile_neutral_cnn_placeholder_vs_decoy_z`: population Z of `0.65 * SCORCH percentile + 0.35 * 0.5`, with CNN absent. The intentional `t_vs_decoys_blend` alias co-matches the same source. | No further source classification decision for scored rows in the frozen v0.1 run. Keep aggregate-score pose linkage separate; do not invent an immutable pose from this classification. |
| Representative known pair and ties | A directly curated, citation-backed drug-target pair is navigation-only and cannot change ranking or validity. | Freeze the accepted evidence-source priority and deterministic tie ordering, then review ambiguous candidates. |
| Public artifact roles | Prepared receptor/ligand, selected docking pose, native ligand/redock pose, images, validation reports, search box, sanitized configuration, release manifest, and hashes are approved; logs and unselected/all-pose artifacts are private. Public projection additionally requires `verified=1` and a valid SHA-256. | Populate and verify the selected artifacts for the frozen release; any later role expansion requires new approval. |

## Engineering approval still required

The small static site remains useful as a local conference backup. Before using one
deployed HTML file per pair for the full matrix, compare its generated file count
and size with the hosting provider's current limits. Cloudflare Pages currently
limits individual assets to 25 MiB; the approximately 439 MB historical public
SQLite projection cannot be uploaded there as a Pages asset. Atlas already
implements an optional, non-deploying bounded architecture prototype: a shared
browser shell plus immutable R2 objects behind a GET/HEAD-only Worker.
The selected delivery architecture is now a provider-neutral static shell plus
immutable object storage; Cloudflare Workers/R2 is the first validated provider
profile, not a data-model dependency. Reconfirm current provider limits and choose
the production domain, resource names, and hosting tier before deployment; none of
these engineering choices changes scientific ranking or qualification policy.

The following engineering contracts are now implemented:

- every master-result row and completion retry remains in a normalized attempt
  ledger; unique successful candidates are selected automatically, while ambiguity
  fails closed unless the release manifest selects one exact path/row plus hash;
- artifact archive members use controlled roles and derived scopes, with explicit
  ligand IDs for ligand/pair evidence; filename guessing no longer establishes pair
  provenance, and only a verified `docking_pose` can satisfy image preflight;
- these artifact roles do not decide scientific validity or public inclusion scope;
- the SQLite exporter streams bounded batches into coarse pair shards and never
  emits one object per receptor-ligand cell; and
- provider-neutral HTTPS download projection and generic/Cloudflare
  credential-free deployment preflight are implemented without uploading or
  deploying;
- public projections omit raw logs, all-pose artifacts, and unselected docking
  poses even from metadata; selected poses require an explicit boolean release
  marker, and every retained public artifact requires `verified=1` plus a valid
  SHA-256; the same filter is applied before the downloadable SQLite is vacuumed;
  and
- HOLO headline ranks and APO exploratory ranks use separate eligibility fields,
  within-receptor ranks, and across-receptor ranks. The browser labels the track
  and the release contract forbids pooling. Exact receptor chemistry, bound score
  evidence, and resolvable score-to-pose lineage are mandatory; the historical
  aggregate scores therefore remain visible but unranked.

Additional engineering choices remain open:

- choose the production download domain and storage provider/tier, create the
  named bucket resources, and authorize a later deployment separately;
- run the full 760,878-cell streaming export only after the frozen public SQLite
  and scientific selection are ready; and
- populate causal completion linkage for selected result rows. Readiness schema 6
  blocks publication when a selected result does not identify its producing
  completion record, link method, structured evidence, and matching selected
  completion attempt;
- extend `atlas stages repeat` from its current planner-only, 25,000-pair-bounded
  checkpoint implementation to a streaming work database and exact-pair
  allowlists/executors for Vina, GNINA, SCORCH, MM/GBSA, and pose validation; and
- keep unknown score-to-pose semantics blocked in every pose-consuming stage plan;
  controlled aggregate sources require complete contributing-pose hashes; and
- choose the partial-stage selection policy. Deterministic hash sampling is
  implemented; top-score sampling remains fail-closed until the exact score
  source, direction, provenance, and tie rule are approved.
