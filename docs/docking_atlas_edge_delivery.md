# Optional Cloudflare Worker + R2 Delivery

The conference-sized static explorer remains the default output of
`atlas publish build`. The optional edge bundle replaces thousands of generated
pair HTML files with one dependency-free browser shell and immutable JSON objects
in R2. It does not change, delete, or redeploy the static fallback.

> **Bounded builder, not the publication-scale exporter.** The current command
> reads `release_browser.json` into memory and materializes pair/entity summaries.
> It refuses browser payloads above 128 MiB and releases above 50,000 pair cells.
> The audited 760,878-cell publication candidate is intentionally rejected. It
> needs a separate SQLite-streaming edge exporter that bypasses generation of the
> monolithic browser payload and per-pair static HTML. Until that exporter exists,
> this bundle validates the delivery architecture for conference/intermediate
> snapshots only; `edge_bundle_summary.json` records
> `publication_scale_supported: false`.

For a bounded conference site that has already been built, verify the ordinary
release first. Do not run the per-pair static build solely to feed this edge command
when the planned matrix exceeds the bounds above.


```bash
atlas publish build --manifest <RELEASE.yaml> --out-dir outputs/data/<RELEASE_ID>
atlas publish verify --site-dir outputs/data/<RELEASE_ID>/site
```

Then create the optional delivery bundle:

```bash
atlas publish edge-bundle \
  --site-dir outputs/data/<RELEASE_ID>/site \
  --out-dir outputs/data/<RELEASE_ID>/edge
```

The command performs no network requests, creates no Cloudflare resources, and
does not deploy. It reads only the path-redacted public
`site/downloads/release_browser.json`. Use `--overwrite` only to replace a directory
with a valid Atlas edge-bundle marker. Source and output must be disjoint: an
output equal to, inside, or containing the static `site/` is rejected before any
deletion. The release token includes the public payload SHA-256, so changed data
receives a new object prefix and route.

## Bundle layout

```text
edge/
├── .atlas-edge-bundle.json
├── edge_bundle_summary.json
├── object_manifest.json
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
        ├── indexes/pairs/<SHARD>.json
        └── records/{targets,drugs,pairs}/<ROUTE_ID>.json
```

`object_manifest.json` is the upload contract. Every R2 key has a byte count,
SHA-256, media type, and immutable cache policy. Pair-index shards are selected by
a deterministic hash and listed in `indexes/pairs.json`. Route IDs are URL-safe,
one-to-one encodings of release identifiers; the release segment makes all target,
drug, and pair routes release-qualified:

```text
/releases/<RELEASE_TOKEN>
/releases/<RELEASE_TOKEN>/targets/<TARGET_ROUTE_ID>
/releases/<RELEASE_TOKEN>/drugs/<DRUG_ROUTE_ID>
/releases/<RELEASE_TOKEN>/pairs/<PAIR_ROUTE_ID>
```

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

## Cloudflare handoff (manual, never automatic)

1. Review `object_manifest.json` and compare its source hash with the public build.
2. Review the generated `wrangler.toml`. Replace both example bucket names with
   separately chosen production and preview bucket names.
3. Create or select a Standard-storage R2 bucket. Do not use Infrequent Access for
   this read-heavy interactive path without separately reviewing retrieval charges.
4. Upload the contents beneath `edge/objects/` so their paths become the exact R2
   keys listed in `object_manifest.json`. Wrangler v4 requires `--remote` for remote
   object operations. Set `Content-Type: application/json` and
   `Cache-Control: public, max-age=31536000, immutable`.
5. Treat `releases/<RELEASE_TOKEN>/` as append-never/replace-never. Because the token
   is content-qualified, a changed release must be rebuilt under a new token. Never
   upload changed bytes to an existing key.
6. Test locally or against a preview bucket, including one `HEAD`, one target, one
   drug, one successful pair, and one failed/invalid pair.
7. Only after explicit deployment approval, run Wrangler from the edge directory.
   Atlas deliberately does not perform this step.

The existing `site/` directory remains a conference-safe local/static backup. Do
not remove it after enabling edge delivery.

## Free plan and $5 Workers plan

The design deliberately keeps only four shared files in Workers Static Assets and
moves record cardinality to R2. As checked on 2026-07-13, Cloudflare documents a
20,000-file limit for Pages Free and up to 100,000 files for paid Pages plans. The
paid limit requires `PAGES_WRANGLER_MAJOR_VERSION=4`. Pages also limits each site
asset to 25 MiB. A generated HTML file for every cell therefore does not scale
safely to a large receptor × FDA matrix, while the shared shell does not grow with
pair count.

The historical candidate's compact public SQLite projection is approximately
439 MB and cannot be uploaded as a Pages asset. The current edge Worker serves the
normalized JSON record tree only; it does not rewrite the static explorer's local
SQLite/CSV/Parquet download links. A Pages deployment therefore needs an external
download host plus a future supported link-projection build step. Alternatively,
serve the complete verified site from a host that accepts the file. Copying the
download to R2 alone does not make the verified Pages bundle deployable.

For a small public beta, the Workers Free plan can be viable when actual traffic and
stored bytes remain within Cloudflare's current quotas. Cloudflare currently lists
an R2 Standard-storage free tier of 10 GB-month, one million Class A operations, and
ten million Class B operations per month, with no Internet egress charge. Initial
object uploads are Class A operations; Worker `get`/`head` reads are Class B.

The current one-object-per-record layout emits approximately
`pairs + targets + drugs + nonempty_pair_shards + 4` R2 objects. The verified
326-pair checkpoint produced 816 objects totaling 1.35 MB (326 pair records, 48
targets, 249 drugs, 189 nonempty shards, and four manifest/index objects). A
760,878-pair release would require at least 760,882 objects before target, drug, and
shard records. One initial upload would therefore consume more than 76% of the
current monthly one-million Class A free allocation, and rebuilding/re-uploading in
the same month could exceed it. The future streaming exporter must report its exact
object/byte plan before upload and may need a coarser pair-shard design.

The Workers Paid plan currently has a $5 USD monthly account minimum and larger
included Worker usage. R2 storage and operations beyond its free allocation remain
usage-based; the $5 minimum is not a fixed all-inclusive hosting price. Set billing
alerts and inspect real object counts, storage, and request volume before public
promotion. Neither plan changes the Atlas scientific release contract.

Current primary references:

- [Cloudflare Pages limits](https://developers.cloudflare.com/pages/platform/limits/)
- [Cloudflare Workers pricing](https://developers.cloudflare.com/workers/platform/pricing/)
- [Cloudflare R2 pricing](https://developers.cloudflare.com/r2/pricing/)
- [Workers SPA static assets](https://developers.cloudflare.com/workers/static-assets/routing/single-page-application/)
- [R2 Workers API](https://developers.cloudflare.com/r2/api/workers/workers-api-reference/)
- [Wrangler R2 commands](https://developers.cloudflare.com/workers/wrangler/commands/r2/)

Cloudflare quotas and prices can change. Recheck those pages immediately before
deployment; this document is an engineering handoff, not a cost guarantee.

## Remaining engineering decisions


- Implement and verify a bounded-memory SQLite-streaming edge exporter before using
  this architecture for the 760,878-cell publication matrix. It must bypass the
  monolithic browser payload and per-pair HTML generator, preserve ranking/failure
  semantics, and emit incremental integrity metadata.
- Choose production and preview bucket names and the release hostname.
- Choose the upload tool and an independent post-upload checksum audit. The bundle
  provides a manifest but intentionally does not infer credentials or mutate R2.
- Decide whether large downloadable SQLite/CSV/Parquet artifacts remain on the
  static fallback, are mirrored under a separate immutable R2 download prefix, or
  are omitted from the edge hostname. The current Worker serves normalized JSON
  records and indexes only, and Atlas does not yet generate cross-origin download
  URLs for the static explorer.
- Decide whether a future release catalog should expose multiple release tokens.
  This v0.1 shell is intentionally pinned to one immutable release.

No scientific policy is selected by this delivery layer. Target identity,
APO/HOLO status, redocking qualification, known-pair selection, score eligibility,
and failure semantics remain exactly as frozen in the source release.
