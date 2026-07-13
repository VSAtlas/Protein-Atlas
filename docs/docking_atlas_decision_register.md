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
  and no secondary score silently replaces a missing primary score.
- Invalid and unsuccessful cells remain present, but are excluded from headline
  rankings by default.
- APO contexts remain visible but are excluded from qualified headline rankings.
  Only an exact `HOLO` variant can pass the current variant gate; blank, legacy,
  and other variants remain visible but excluded.
- Any non-empty `receptor_classification` value is preserved but excluded as
  `receptor_classification_policy_pending` until the user approves the controlled
  vocabulary and conflict policy; readiness schema 3 reports this as a blocker.
- A receptor must pass an explicit receptor-quality outcome and an explicit native
  redock qualification before its drug rows can enter headline rankings.
- The native-redock RMSD threshold is provisionally 2.5 angstrom, pending the
  calculation details below.
- The conference image set contains the native control, up to five top valid
  ligands, and one frozen known pair when it is not already selected. It never
  selects a best-invalid pose.
- For frozen public releases, protein identities, receptor outcomes, and known pairs
  come from explicit annotation rows backed by hashed external sources. Older inline
  known-pair values may remain covered by a manifest hash, but `atlas publish build`
  ignores them; only the advanced callable's explicit parameter remains
  compatibility-only.

## Current evidence snapshot (2026-07-13)

The engineering-only audit of
`spd90_fda_fda_dud_spr10_40h_patchdist_20260525_161503` contains 90 receptor
contexts, 760,878 terminal expected cells, 440,183 calculated-but-unvalidated
cells, 320,695 structured missing-result cells, and 31,950 non-null
`final_score` values. It has zero explicit pair pose validations, zero explicit
native-redock qualifications, no imported artifact index, and no classified
legacy `final_score` source. It is not approved as a scientific release.

The canonical status snapshot for `postFDAfixSPDrun_prepcheck_20260713` reported
92 of 93 receptor contexts complete, one pending context, and two structured
control-redock `fallback_missing` examples. Refresh canonical status before review;
do not convert those observations into pass/fail outcomes automatically.

## Scientific approvals still required

| Decision | Evidence or input to review | Approval to record |
|---|---|---|
| Conference receptor set | Compare the completed historical 90-context matrix with the current 93-context preparation check and document inclusions/exclusions. | Frozen list of exact run, PDB, variant, and pH keys. |
| Experimental receptor quality policy | Review source structure metadata and the fields that matter for this study, such as experimental method, resolution/quality metrics, bound-ligand identity, chain, mutations, missing residues, and pocket completeness. | Policy identifier, required fields/cutoffs, per-context outcome, and reason. |
| APO/HOLO classification | Freeze a controlled vocabulary and resolve run-label versus source-record conflicts from authoritative structure evidence; free text such as `non-apo` or `apoprotein` must not be interpreted heuristically. | Approved tokens, conflict policy, and explicit classification per exact receptor context. |
| Native-redock RMSD method | The 2.5 angstrom threshold alone does not define atom mapping, symmetry handling, alignment frame, hydrogen treatment, or which generated pose is evaluated. | Named method/version plus per-context RMSD and explicit status. |
| Pair pose-validity policy | Define the checks and acceptable thresholds for geometry, ligand integrity, pocket placement, clashes, and any engine-specific exceptions. | Validator version, structured status vocabulary, and per-pair outputs. |
| Legacy `final_score` provenance | Reconstruct how each of the 31,950 historical values was materialized; do not label a source from the value alone. | Source field/version for every non-null primary score. |
| Known drug-target pair per receptor | Use the frozen ChEMBL mechanism snapshot and surface ambiguous/tied candidates. | One canonical ligand and evidence reference per context; user approval for every tie. |
| Public artifact scope | Decide which prepared inputs, poses, logs, and raw outputs can be public and which remain only in the private evidence ledger. | Artifact inclusion policy and immutable archive index. |

## Engineering approval still required

The small static site remains useful as a local conference backup. Before using one
deployed HTML file per pair for the full matrix, compare its generated file count
and size with the hosting provider's current limits. Cloudflare Pages currently
limits individual assets to 25 MiB; the approximately 439 MB historical public
SQLite projection cannot be uploaded there as a Pages asset. Atlas already
implements an optional, non-deploying bounded architecture prototype: a shared
browser shell plus immutable R2 objects behind a GET/HEAD-only Worker. Approve that
delivery change,
confirm current provider limits, and choose the desired hosting tier before
deployment; none of these engineering choices changes scientific ranking or
qualification policy.

Additional engineering choices are intentionally open:

- keep duplicate canonical master-result rows as a build error for v0.1, or add a
  normalized result-attempt table plus an explicit release selector;
- implement the SQLite-streaming edge exporter before attempting the 760,878-cell
  publication matrix;
- define typed pose/input/log artifact roles so a verified generic pair artifact
  cannot be mistaken for a verified renderable pose.
