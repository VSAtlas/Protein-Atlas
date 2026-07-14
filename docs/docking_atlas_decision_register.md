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
- APO contexts remain visible but are excluded from the current fail-closed
  qualified headline ranking. Ranking uses the controlled observed
  `receptor_classification`, not the requested variant/folder token. Whether a
  separately labelled APO ranking should also be public remains a scientific
  decision.
- Receptor annotations accept only controlled `APO`/`HOLO` classifications or
  an unresolved blank. Readiness schema 5 reports missing/uncontrolled
  classification and exact prepared-receptor evidence gaps rather than treating a
  historical run label as chemistry truth.
- A receptor must pass an explicit receptor-quality outcome and an explicit native
  redock qualification before its drug rows can enter headline rankings.
- The native-redock RMSD threshold is provisionally 2.5 angstrom. The current RMSD
  implementation must be verified against the approved method before any receptor
  is qualified.
- The conference image set contains the native control, up to five top valid
  ligands, the best-scoring invalid result when present, and a representative
  known-pair slot when it is not already selected. The slot is approved;
  its scientific definition and tie handling are not.
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
native-redock qualifications, no imported artifact index, and no classified
legacy `final_score` source. It is not approved as a scientific release.

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
- The provisional native-redock cutoff is 2.5 angstrom. No pass may be recorded
  until the RMSD calculation method and the current implementation are verified.
- Pose-validity thresholds come from the existing validation script. The public
  database must show the validator version/configuration and thresholds rather than
  silently reinterpreting them.
- Normalized primary scores may be compared only when the score-source family is
  explicit and displayed alongside the value.
- The conference image plan includes the native control, up to five top valid
  ligands, the best-scoring invalid result when present, and a representative known
  pair when it is not already among those selections.
- The approved public-artifact baseline is the prepared receptor, prepared ligand,
  and docked-ligand images. This is not approval to publish every pose, validation
  output, log, or raw intermediate.

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
| APO/HOLO conflict rules and ranking use | Classification uses exact prepared chemistry rather than folder/run labels; disagreement remains unresolved. | Controlled treatment of cofactors, nucleotide/native-ligand cases, free ions/salts, and halides; whether labelled APO contexts receive a separate public ranking or remain outside headline ranks; then explicit review of the 15 unresolved contexts. |
| Native-redock RMSD method | Provisional cutoff is 2.5 angstrom. | Atom mapping, symmetry treatment, alignment frame, hydrogen treatment, generated-pose selection, named method/version, verification of the present implementation, and per-context outcome. |
| Pair pose-validity scope | Existing script thresholds must be shown. | Exact pass predicates and combination rule, whether validation applies to the selected `final_score` pose or any stage/pose, engine-specific exceptions, validator version, and per-pair output mapping. |
| Legacy `final_score` provenance | Comparison requires an explicit displayed source family. | Reconstruct and freeze the source field/family/version for every historical non-null primary score; never infer it from the numeric value alone. |
| Representative known pair and ties | A known pair is a navigation/presentation annotation intended to give each receptor page one recognizable drug-target relationship when it is absent from the top-five valid results. It does not, by itself, validate a docking result or change ranking. | Evidence-source definition, eligibility rules, canonical ligand selection, tie handling, and approval of ambiguous candidates. |
| Additional public artifact roles | Prepared receptor, prepared ligand, and docked-ligand images are approved as the baseline. | Whether downloadable poses, receptor-ligand complexes, validation reports, stdout/stderr, raw logs, all poses/scores, and other intermediates are public or remain only in the private evidence ledger. |

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
  deploying.

Additional engineering choices remain open:

- choose the production download domain and storage provider/tier, create the
  named bucket resources, and authorize a later deployment separately;
- run the full 760,878-cell streaming export only after the frozen public SQLite
  and scientific selection are ready; and
- populate causal completion linkage for selected result rows. Readiness schema 5
  blocks publication when a selected result does not identify its producing
  completion record, link method, structured evidence, and matching selected
  completion attempt.
