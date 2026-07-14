# Optional Cloudflare Worker + R2 Delivery

The conference-sized static explorer remains the default output of
`atlas publish build`. The optional edge bundle replaces thousands of generated
pair HTML files with one dependency-free browser shell and immutable JSON objects
in R2. It does not change, delete, or redeploy the static fallback.

Two build modes intentionally coexist:

- The existing `--site-dir` mode reads `release_browser.json` and preserves the
  conference-sized browser-JSON path. It remains capped at a 128 MiB payload and
  50,000 pair cells and records `publication_scale_supported: false`.
- The `--database` mode reads the compact public SQLite snapshot directly. It
  streams pair rows in bounded batches to a disposable on-disk work database,
  computes rankings with SQLite window functions, and writes coarse pair-detail
  shards instead of one object per cell. It records
  `publication_scale_supported: true`.

The audited 760,878-cell candidate has not yet been run through the streaming mode.
The implementation removes the known memory and object-count blockers, but a full
offline rehearsal, output-size measurement, and deployment preflight are still
required before publication. Do not generate per-pair static HTML solely to feed
the streaming edge path.

```bash
atlas publish build --manifest <RELEASE.yaml> --out-dir outputs/data/<RELEASE_ID>
atlas publish verify --site-dir outputs/data/<RELEASE_ID>/site
```

Then create the optional delivery bundle:

```bash
atlas publish edge-bundle \
  --site-dir outputs/data/<RELEASE_ID>/site \
  --out-dir outputs/data/<RELEASE_ID>/edge \
  --download-base-url https://downloads.example.org
```

For a larger matrix, use the same canonical command against the path-redacted
public SQLite snapshot:

```bash
atlas publish edge-bundle \
  --database outputs/data/<RELEASE_ID>/site/downloads/docking_atlas.sqlite \
  --source-site-dir outputs/data/<RELEASE_ID>/site \
  --out-dir outputs/data/<RELEASE_ID>/edge-streamed \
  --batch-rows 1000 \
  --coarse-shard-rows 2000 \
  --download-base-url https://downloads.example.org
```

`--batch-rows` and `--coarse-shard-rows` accept 100 through 10,000 rows. The
default pair ceiling is 2,000,000 and can be lowered with `--max-pairs`. Source
SQLite schema versions 2 through 5 are accepted; absent legacy validation or
receptor-chemistry fields remain null and therefore cannot silently qualify a row
for headline ranking. The existing browser-JSON command remains unchanged.

`--download-base-url` is a provider-neutral HTTPS object origin. Atlas validates
every download declared by `site_manifest.json`, checks its SHA-256 against the
local bytes, and projects it to an immutable
`releases/<RELEASE_TOKEN>/downloads/<FILE>` key. A Cloudflare R2 custom-domain URL
must be origin-only (no path prefix); other providers may mount an object root below
a path. No upload occurs.

Both modes perform no network requests, create no hosting resources, and do not
deploy. The bounded mode reads only the path-redacted
`site/downloads/release_browser.json`; the streaming mode reads the path-redacted
public SQLite snapshot with a read-only immutable connection and refuses live WAL
or SHM sidecars. It verifies that the source file identity is unchanged before
finalizing output and removes its disposable work database. Use `--overwrite` only
to replace a directory with a valid Atlas edge-bundle marker. Source and output
must be disjoint. The release token includes the public source SHA-256, so changed
data receives a new object prefix and route.

## Bundle layout

```text
edge/
├── .atlas-edge-bundle.json
├── edge_bundle_summary.json
├── object_manifest.json             # normalized JSON-object upload contract
├── download_projection.json         # provider-neutral download upload/URL contract
├── deployment_preflight.json        # written by `atlas publish preflight`
├── wrangler.preflight.toml          # written only after valid Cloudflare names
├── wrangler.toml                    # template; review bucket names
├── src/index.mjs                    # GET/HEAD-only Worker
├── public/
│   ├── _headers                     # static response security headers
│   ├── index.html                   # shared SPA shell
│   ├── release-config.json
│   └── assets/{app.css,app.js}
└── objects/
    └── releases/<RELEASE_TOKEN>/
        ├── manifest.json
        ├── indexes/{targets,drugs,pairs}.json
        ├── indexes/pairs/<SHARD>.json        # bounded JSON mode
        ├── records/{targets,drugs,pairs}/<ROUTE_ID>.json  # bounded JSON mode
        └── records/pair-shards/<SHARD>.json  # streaming SQLite mode
```

`object_manifest.json` is the normalized-record upload contract. Every R2 key has a
byte count, SHA-256, media type, and immutable cache policy.
`download_projection.json` is a separate provider-neutral contract for the frozen
SQLite/CSV/Parquet/JSON downloads. It records only relative source paths, immutable
object keys, HTTPS public URLs, sizes, media types, cache policy, and SHA-256 values;
it contains no credentials or machine-local paths. In bounded JSON mode, pair-index
shards are selected by a deterministic hash. In streaming mode,
`indexes/pairs.json` lists ordinal coarse shards and every target/drug pair summary
includes its `pair_shard_id`; pair details are never emitted as one object per cell.
Route IDs are URL-safe,
one-to-one encodings of release identifiers; the release segment makes all target,
drug, and pair routes release-qualified:

```text
/releases/<RELEASE_TOKEN>
/releases/<RELEASE_TOKEN>/targets/<TARGET_ROUTE_ID>
/releases/<RELEASE_TOKEN>/drugs/<DRUG_ROUTE_ID>
/releases/<RELEASE_TOKEN>/pairs/<PAIR_ROUTE_ID>
```

Streaming pair links append `?shard=<PAIR_SHARD_ID>` so the browser fetches one
coarse object and selects the requested route locally.

The Worker exposes matching `/api/releases/...` reads. It accepts only `GET` and
`HEAD`, validates every path segment against a narrow allowlist, maps only known
route shapes to R2 keys, and uses only R2 `get` and `head`. It contains no bucket
listing, mutation, upload, delete, authentication secret, or administrative route.
Static navigation is handled by Workers Static Assets in single-page-application
mode. The generated `public/_headers` file applies an enforced same-origin content
security policy, frame denial, MIME-sniffing protection, a no-referrer policy, and
a restrictive permissions policy to direct static responses even though
`run_worker_first` is limited to `/api/*`. The HTML meta policy is therefore not
the sole security boundary.

## Credential-free deployment preflight

Run the preflight after choosing names; this command performs no network request,
credential lookup, bucket creation, upload, or deployment:

```bash
atlas publish preflight \
  --edge-dir outputs/data/<RELEASE_ID>/edge \
  --site-dir outputs/data/<RELEASE_ID>/site \
  --provider cloudflare-workers-r2 \
  --cloudflare-plan free \
  --worker-name docking-atlas-edge \
  --production-bucket atlas-docking-releases \
  --preview-bucket atlas-docking-releases-preview
```

The report fails closed on changed object/download bytes, missing or duplicate keys,
symlinks, invalid resource names, absent choices, non-HTTPS download URLs, or current
Cloudflare static-asset count/size limits. A ready Cloudflare preflight writes
`wrangler.preflight.toml` with only non-secret resource names. The report explicitly
records `secret_values_inspected: false`, `network_checks_performed: false`,
`upload_performed: false`, and `deployment_performed: false`. Generic providers use
the same projection and integrity checks but must supply their own host limits.

## Cloudflare account handoff (manual, never automatic)

The preflight reports the exact remaining resources and credential environment names.
For the current Worker + R2 adapter, choose or create:

1. One production Standard-storage R2 bucket matching `--production-bucket`.
2. One distinct preview Standard-storage R2 bucket matching `--preview-bucket`.
3. One Worker service name matching `--worker-name`; the service is created only by
   a future explicitly approved Wrangler deployment, not by Atlas preflight.
4. One production custom domain attached to the download bucket. Use its origin-only
   HTTPS URL as `--download-base-url`. Cloudflare documents `r2.dev` as a
   non-production development URL, so it is not the frozen public-release choice.

Wrangler uses `CLOUDFLARE_ACCOUNT_ID` and `CLOUDFLARE_API_TOKEN`. The preflight lists
those names but never reads or prints their values. Do not place token values in the
manifest, generated TOML, repository, or conference materials.

After explicit upload/deployment approval in a later task, upload `edge/objects/`
using the exact keys and metadata in `object_manifest.json`, and upload the source
downloads using `download_projection.json`. Treat
`releases/<RELEASE_TOKEN>/` as append-never/replace-never. Then test one `HEAD`, one
target, one drug, one successful pair, one failed/invalid pair, and every projected
download hash. Atlas deliberately does none of those remote operations in this patch.

The existing `site/` directory remains a conference-safe local/static backup. Do
not remove it after enabling edge delivery.

## Free plan and $5 Workers plan

The design deliberately keeps only five shared files in Workers Static Assets and
moves record cardinality to R2. As checked on 2026-07-13, Cloudflare documents a
20,000-file limit for Pages Free and up to 100,000 files for paid Pages plans. The
paid limit requires `PAGES_WRANGLER_MAJOR_VERSION=4`. Pages also limits each site
asset to 25 MiB. A generated HTML file for every cell therefore does not scale
safely to a large receptor × FDA matrix, while the shared shell does not grow with
pair count.

The historical candidate's compact public SQLite projection is approximately
439 MB and cannot be uploaded as a Pages or Workers Static Assets file. The edge
Worker serves normalized JSON records, while `--download-base-url` projects verified
SQLite/CSV/Parquet/JSON links to a separate HTTPS object origin. It does not rewrite
the unchanged static fallback. Therefore deploy the shared edge shell for this
architecture; copying downloads to R2 alone still does not make the original Pages
bundle deployable.

For a small public beta, the Workers Free plan can be viable when actual traffic and
stored bytes remain within Cloudflare's current quotas. Cloudflare currently lists
an R2 Standard-storage free tier of 10 GB-month, one million Class A operations, and
ten million Class B operations per month, with no Internet egress charge. Initial
object uploads are Class A operations; Worker `get`/`head` reads are Class B.

The bounded JSON layout still emits approximately
`pairs + targets + drugs + nonempty_pair_shards + 4` R2 objects. Its verified
326-pair checkpoint produced 816 objects totaling 1.35 MB. The streaming layout
emits approximately
`targets + drugs + ceil(pairs / coarse_shard_rows) + 4` objects and no per-cell
pair objects. The same 326-pair public checkpoint produced 302 objects totaling
1,312,013 bytes: 48 target records, 249 drug records, one coarse pair-detail shard,
and four manifest/index objects.

At the default 2,000 rows per coarse shard, the audited inventory of 90 targets,
8,645 drugs, and 760,878 pairs projects to 9,120 JSON objects, including 381 pair
shards. That is a count projection, not a completed full export. Pair summaries are
intentionally present in both target and drug records for fast bidirectional
navigation, so the full offline rehearsal must still measure bytes, largest object,
rendering latency, and upload operations before deployment.

The Workers Paid plan currently has a $5 USD monthly account minimum and larger
included Worker usage. R2 storage and operations beyond its free allocation remain
usage-based; the $5 minimum is not a fixed all-inclusive hosting price. Set billing
alerts and inspect real object counts, storage, and request volume before public
promotion. Neither plan changes the Atlas scientific release contract.

Current primary references:

- [Cloudflare Pages limits](https://developers.cloudflare.com/pages/platform/limits/)
- [Cloudflare Workers limits](https://developers.cloudflare.com/workers/platform/limits/)
- [Cloudflare Workers pricing](https://developers.cloudflare.com/workers/platform/pricing/)
- [Cloudflare R2 pricing](https://developers.cloudflare.com/r2/pricing/)
- [Workers SPA static assets](https://developers.cloudflare.com/workers/static-assets/routing/single-page-application/)
- [R2 Workers API](https://developers.cloudflare.com/r2/api/workers/workers-api-reference/)
- [R2 public buckets](https://developers.cloudflare.com/r2/buckets/public-buckets/)
- [Wrangler R2 commands](https://developers.cloudflare.com/workers/wrangler/commands/r2/)
- [Wrangler environment variables](https://developers.cloudflare.com/workers/wrangler/system-environment-variables/)

Cloudflare quotas and prices can change. Recheck those pages immediately before
deployment; this document is an engineering handoff, not a cost guarantee.

## Remaining engineering decisions


- Run the first full 760,878-cell export only after the public SQLite snapshot and
  scientific release selection are frozen. Measure total bytes, largest target/drug
  record, shard latency, object count, and preflight results before any upload.
- Confirm whether 2,000 pair details per shard is the publication default or tune it
  within the enforced 100-to-10,000 range from the offline measurements.
- Choose production and preview bucket names, a Worker name, and the public download
  custom domain; record them by running `atlas publish preflight`.
- Choose the future upload tool and an independent post-upload checksum audit. The
  bundle provides a manifest but intentionally does not infer credentials or mutate
  R2.
- Decide whether the static fallback retains its local SQLite/CSV/Parquet copies
  after the edge release is deployed. The edge SPA now renders provider-neutral
  projected HTTPS download links, while the current Worker continues to serve only
  normalized JSON records and indexes.
- Decide whether a future release catalog should expose multiple release tokens.
  This v0.1 shell is intentionally pinned to one immutable release.

No scientific policy is selected by this delivery layer. Target identity,
APO/HOLO status, redocking qualification, known-pair selection, score eligibility,
and failure semantics remain exactly as frozen in the source release.
