"""Embedded, dependency-free assets for the optional Atlas edge explorer."""

from __future__ import annotations

WORKER_SOURCE = r"""const TOKEN = /^[A-Za-z0-9_-]{1,512}$/;
const SHARD = /^[0-9a-f]{2}$/;
const IMMUTABLE = "public, max-age=31536000, immutable";

function objectKey(pathname) {
  const parts = pathname.split("/").filter(Boolean);
  if (parts.length < 4 || parts[0] !== "api" || parts[1] !== "releases") {
    return null;
  }
  const release = parts[2];
  if (!TOKEN.test(release)) return null;
  const prefix = `releases/${release}`;
  if (parts.length === 4 && parts[3] === "manifest") {
    return `${prefix}/manifest.json`;
  }
  if (parts[3] === "indexes") {
    if (parts.length === 5 && ["targets", "drugs", "pairs"].includes(parts[4])) {
      return `${prefix}/indexes/${parts[4]}.json`;
    }
    if (parts.length === 6 && parts[4] === "pairs" && SHARD.test(parts[5])) {
      return `${prefix}/indexes/pairs/${parts[5]}.json`;
    }
    return null;
  }
  if (parts.length === 5 && ["targets", "drugs", "pairs"].includes(parts[3])) {
    if (!TOKEN.test(parts[4])) return null;
    return `${prefix}/records/${parts[3]}/${parts[4]}.json`;
  }
  return null;
}

function headersFor(object) {
  const headers = new Headers();
  object.writeHttpMetadata(headers);
  headers.set("Content-Type", headers.get("Content-Type") || "application/json; charset=utf-8");
  headers.set("Cache-Control", IMMUTABLE);
  headers.set("ETag", object.httpEtag);
  headers.set("X-Content-Type-Options", "nosniff");
  return headers;
}

const STATIC_SECURITY_HEADERS = {
  "Content-Security-Policy": "default-src 'self'; connect-src 'self'; img-src 'self' data:; style-src 'self'; script-src 'self'; object-src 'none'; base-uri 'none'; frame-ancestors 'none'; form-action 'none'",
  "Permissions-Policy": "camera=(), geolocation=(), microphone=(), payment=(), usb=()",
  "Referrer-Policy": "no-referrer",
  "X-Content-Type-Options": "nosniff",
  "X-Frame-Options": "DENY",
};

function secureStaticResponse(response) {
  const headers = new Headers(response.headers);
  Object.entries(STATIC_SECURITY_HEADERS).forEach(([name, value]) => {
    headers.set(name, value);
  });
  return new Response(response.body, {
    status: response.status,
    statusText: response.statusText,
    headers,
  });
}

export default {
  async fetch(request, env) {
    if (request.method !== "GET" && request.method !== "HEAD") {
      return new Response("Method not allowed", {
        status: 405,
        headers: { Allow: "GET, HEAD", "X-Content-Type-Options": "nosniff" },
      });
    }
    const url = new URL(request.url);
    if (!url.pathname.startsWith("/api/")) {
      return secureStaticResponse(await env.ASSETS.fetch(request));
    }
    const key = objectKey(url.pathname);
    if (key === null) {
      return new Response("Not found", {
        status: 404,
        headers: { "Cache-Control": "no-store", "X-Content-Type-Options": "nosniff" },
      });
    }
    const object = request.method === "HEAD"
      ? await env.ATLAS_RELEASES.head(key)
      : await env.ATLAS_RELEASES.get(key);
    if (object === null) {
      return new Response("Not found", {
        status: 404,
        headers: { "Cache-Control": "no-store", "X-Content-Type-Options": "nosniff" },
      });
    }
    const headers = headersFor(object);
    if (request.method === "HEAD") {
      headers.set("Content-Length", String(object.size));
      return new Response(null, { status: 200, headers });
    }
    return new Response(object.body, { status: 200, headers });
  },
};
"""

WRANGLER_TEMPLATE = """name = "docking-atlas-edge"
main = "src/index.mjs"
compatibility_date = "2026-07-13"

[assets]
directory = "./public"
binding = "ASSETS"
not_found_handling = "single-page-application"
run_worker_first = ["/api/*"]

[[r2_buckets]]
binding = "ATLAS_RELEASES"
bucket_name = "atlas-docking-releases"
preview_bucket_name = "atlas-docking-releases-preview"
"""

PUBLIC_HEADERS = """/*
  Content-Security-Policy: default-src 'self'; connect-src 'self'; img-src 'self' data:; style-src 'self'; script-src 'self'; object-src 'none'; base-uri 'none'; frame-ancestors 'none'; form-action 'none'
  Permissions-Policy: camera=(), geolocation=(), microphone=(), payment=(), usb=()
  Referrer-Policy: no-referrer
  X-Content-Type-Options: nosniff
  X-Frame-Options: DENY

/assets/*
  Cache-Control: public, max-age=3600, must-revalidate
"""

INDEX_HTML = r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <meta name="color-scheme" content="light">
  <meta http-equiv="Content-Security-Policy" content="default-src 'self'; connect-src 'self'; img-src 'self' data:; style-src 'self'; script-src 'self'; object-src 'none'; base-uri 'none'; frame-ancestors 'none'">
  <title>Docking Atlas</title>
  <link rel="stylesheet" href="/assets/app.css">
</head>
<body>
  <header><a id="brand" href="/">Docking Atlas</a><span id="release-label"></span></header>
  <main id="app"><p class="loading">Loading release…</p></main>
  <footer>Auditable docking records · invalid and unsuccessful cells remain visible</footer>
  <script src="/assets/app.js" defer></script>
</body>
</html>
"""

APP_CSS = r""":root{--ink:#18231c;--muted:#657068;--paper:#f4f6f1;--card:#fff;--line:#d8dfd6;--green:#24583a;--red:#8c3434;--gold:#87651d}*{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);font:16px/1.5 system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}header{height:64px;padding:0 max(20px,calc((100vw - 1180px)/2));display:flex;align-items:center;justify-content:space-between;background:#183626;color:#fff}header a{color:#fff;text-decoration:none;font-size:1.15rem;font-weight:800}header span{font-size:.82rem;opacity:.8}main{max-width:1180px;margin:auto;padding:36px 22px 80px}footer{padding:22px;text-align:center;color:var(--muted)}a{color:var(--green)}h1{font-size:clamp(2rem,5vw,3.7rem);line-height:1.05;margin:.15em 0}.hero{padding:38px;border-radius:20px;background:linear-gradient(135deg,#e4eee3,#fff);border:1px solid var(--line);margin-bottom:28px}.hero p{max-width:780px;color:#405047}.eyebrow{text-transform:uppercase;letter-spacing:.12em;font-size:.72rem;font-weight:800;color:var(--green)}.metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(140px,1fr));gap:12px;margin:22px 0}.metric,.panel{background:var(--card);border:1px solid var(--line);border-radius:14px;padding:18px}.metric strong{display:block;font-size:1.8rem}.metric span,.muted{color:var(--muted)}.split{display:grid;grid-template-columns:1fr 1fr;gap:18px}.section-head{display:flex;align-items:center;justify-content:space-between;gap:12px}.section-head input{width:min(330px,100%);padding:10px 12px;border:1px solid #aeb9ad;border-radius:9px}.cards{display:grid;gap:10px}.card{display:block;padding:14px;border:1px solid var(--line);border-radius:11px;background:#fff;text-decoration:none;color:inherit}.card:hover{border-color:#8ea292}.card strong,.card span{display:block}.card span{color:var(--muted);font-size:.84rem}.table-wrap{overflow:auto;background:#fff;border:1px solid var(--line);border-radius:12px}table{width:100%;border-collapse:collapse}th,td{text-align:left;padding:11px;border-bottom:1px solid var(--line);vertical-align:top}th{font-size:.72rem;text-transform:uppercase;letter-spacing:.05em;color:var(--muted)}.badge{display:inline-block;padding:3px 8px;border-radius:999px;background:#e8ece7;font-size:.75rem;font-weight:700}.valid,.qualified,.success{background:#d9ecdf;color:#205034}.invalid,.warning{background:#fff0c9;color:#6e500e}.failed,.missing,.error{background:#f3dada;color:#762828}.breadcrumbs{margin-bottom:17px;color:var(--muted);font-size:.84rem}.record{display:grid;grid-template-columns:minmax(140px,.6fr) 1.5fr;margin:0}.record dt,.record dd{padding:9px 0;border-bottom:1px solid var(--line);overflow-wrap:anywhere}.record dt{font-size:.78rem;font-weight:700;color:var(--muted)}.record dd{margin:0}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#edf0ea;border-radius:10px;padding:15px;max-height:620px;overflow:auto}.pager{display:flex;justify-content:center;align-items:center;gap:12px;margin:16px}.pager button{padding:7px 12px;border:1px solid #abb7ab;border-radius:8px;background:#fff;color:var(--green);cursor:pointer}.pager button:disabled{opacity:.4;cursor:default}.empty,.error-box{padding:18px;border:1px dashed #b8c2b8;border-radius:11px;background:#fff}.error-box{border-color:#d1aaaa;color:var(--red)}@media(max-width:780px){main{padding:22px 13px 60px}.hero{padding:26px 20px}.split{grid-template-columns:1fr}.section-head{align-items:stretch;flex-direction:column}.record{grid-template-columns:1fr}.record dd{padding-top:0}}
"""

APP_JS = r"""(() => {
  "use strict";
  const app = document.querySelector("#app");
  const releaseLabel = document.querySelector("#release-label");
  const brand = document.querySelector("#brand");
  const TOKEN = /^[A-Za-z0-9_-]{1,512}$/;
  let config;

  function element(tag, text, className) {
    const node = document.createElement(tag);
    if (text !== undefined) node.textContent = String(text);
    if (className) node.className = className;
    return node;
  }

  function link(text, href, className) {
    const node = element("a", text, className);
    node.href = href;
    return node;
  }

  function route(kind, id) {
    const base = `/releases/${config.release_token}`;
    return kind ? `${base}/${kind}/${id}` : base;
  }

  async function getJson(relative) {
    const response = await fetch(`/api/releases/${config.release_token}/${relative}`, {
      headers: { Accept: "application/json" },
    });
    if (!response.ok) throw new Error(`Request failed (${response.status})`);
    return response.json();
  }

  function statusBadge(value) {
    const text = value || "unknown";
    const className = `badge ${String(text).toLowerCase().replace(/[^a-z]+/g, "-")}`;
    return element("span", text, className);
  }

  function metric(label, value) {
    const node = element("div", undefined, "metric");
    node.append(element("strong", value), element("span", label));
    return node;
  }

  function hero(title, eyebrow, detail) {
    const node = element("section", undefined, "hero");
    node.append(element("p", eyebrow, "eyebrow"), element("h1", title));
    if (detail) node.append(element("p", detail));
    return node;
  }

  function breadcrumbs(items) {
    const node = element("nav", undefined, "breadcrumbs");
    items.forEach((item, index) => {
      if (index) node.append(document.createTextNode(" / "));
      node.append(item.href ? link(item.label, item.href) : document.createTextNode(item.label));
    });
    return node;
  }

  function recordList(values) {
    const list = element("dl", undefined, "record");
    Object.entries(values).forEach(([key, value]) => {
      if (value === null || value === undefined || value === "") return;
      list.append(element("dt", key.replaceAll("_", " ")));
      const display = typeof value === "object" ? JSON.stringify(value) : value;
      list.append(element("dd", display));
    });
    return list;
  }

  function searchableCards(container, records, kind) {
    const input = element("input");
    input.type = "search";
    input.placeholder = `Filter ${kind}…`;
    const cards = element("div", undefined, "cards");
    const count = element("p", undefined, "muted");
    function draw() {
      const query = input.value.trim().toLowerCase();
      const filtered = records.filter((row) =>
        `${row.label} ${row.id}`.toLowerCase().includes(query)
      );
      cards.replaceChildren();
      filtered.slice(0, 250).forEach((row) => {
        const card = link(row.label, route(kind, row.route_id), "card");
        card.append(element("span", `${row.pair_count} pair cells · ${row.primary_score_count} primary scores`));
        cards.append(card);
      });
      count.textContent = `${Math.min(filtered.length, 250)} of ${filtered.length} shown`;
    }
    input.addEventListener("input", draw);
    container.append(input, count, cards);
    draw();
  }

  function pairTable(rows) {
    const host = element("div");
    let page = 0;
    const pageSize = 100;
    function draw() {
      const tableWrap = element("div", undefined, "table-wrap");
      const table = element("table");
      const head = element("tr");
      ["Counterpart", "Status", "Final score", "Rank", "Reason", "Record"].forEach((name) => head.append(element("th", name)));
      const thead = element("thead");
      thead.append(head);
      const body = element("tbody");
      rows.slice(page * pageSize, (page + 1) * pageSize).forEach((row) => {
        const tr = element("tr");
        const counterpart = element("td");
        counterpart.append(link(row.counterpart_label, route(row.counterpart_kind, row.counterpart_route_id)));
        const status = element("td");
        status.append(statusBadge(row.final_status));
        tr.append(
          counterpart,
          status,
          element("td", row.final_score ?? "—"),
          element("td", row.rank ?? "—"),
          element("td", row.failure_code || row.ranking_eligibility_reason || "—"),
        );
        const record = element("td");
        record.append(link("Open", route("pairs", row.pair_route_id)));
        tr.append(record);
        body.append(tr);
      });
      table.append(thead, body);
      tableWrap.append(table);
      const pager = element("div", undefined, "pager");
      const previous = element("button", "Previous");
      const next = element("button", "Next");
      previous.disabled = page === 0;
      next.disabled = (page + 1) * pageSize >= rows.length;
      previous.addEventListener("click", () => { page -= 1; draw(); });
      next.addEventListener("click", () => { page += 1; draw(); });
      pager.append(previous, element("span", `Page ${page + 1} of ${Math.max(1, Math.ceil(rows.length / pageSize))}`), next);
      host.replaceChildren(tableWrap, pager);
    }
    if (!rows.length) return element("p", "No pair cells are available.", "empty");
    draw();
    return host;
  }

  async function renderRelease() {
    const [manifest, targets, drugs] = await Promise.all([
      getJson("manifest"), getJson("indexes/targets"), getJson("indexes/drugs"),
    ]);
    const root = element("div");
    root.append(hero(
      manifest.release.title || "Docking Atlas",
      `Release ${manifest.release.id || config.release_id}`,
      manifest.release.summary || manifest.release.description || "Bidirectional receptor–drug exploration with complete failure visibility.",
    ));
    const metrics = element("section", undefined, "metrics");
    metrics.append(
      metric("receptor contexts", manifest.counts.targets),
      metric("drugs", manifest.counts.drugs),
      metric("pair cells", manifest.counts.pairs),
      metric("rank eligible", manifest.coverage.rank_eligible_count || 0),
    );
    root.append(metrics);
    const split = element("section", undefined, "split");
    const targetPanel = element("div", undefined, "panel");
    const drugPanel = element("div", undefined, "panel");
    targetPanel.append(element("h2", "Targets / receptors"));
    drugPanel.append(element("h2", "Drugs / ligands"));
    searchableCards(targetPanel, targets.records, "targets");
    searchableCards(drugPanel, drugs.records, "drugs");
    split.append(targetPanel, drugPanel);
    root.append(split);
    app.replaceChildren(root);
  }

  async function renderEntity(kind, routeId) {
    const payload = await getJson(`${kind}/${routeId}`);
    const entity = payload.entity;
    const label = entity.display_name || entity.ligand_display_name || entity.id || payload.id;
    const root = element("div");
    root.append(
      breadcrumbs([{ label: "Release", href: route() }, { label }]),
      hero(label, kind === "targets" ? "Receptor context" : "Drug / ligand", `${payload.pairs.length} pair cells`),
    );
    const panel = element("section", undefined, "panel");
    panel.append(element("h2", "Record"), recordList(entity));
    root.append(panel, element("h2", "Pair cells"), pairTable(payload.pairs));
    app.replaceChildren(root);
  }

  async function renderPair(routeId) {
    const payload = await getJson(`pairs/${routeId}`);
    const pair = payload.pair;
    const root = element("div");
    root.append(
      breadcrumbs([
        { label: "Release", href: route() },
        { label: payload.target.label, href: route("targets", payload.target.route_id) },
        { label: payload.drug.label, href: route("drugs", payload.drug.route_id) },
        { label: "Pair" },
      ]),
      hero(`${payload.drug.label} × ${payload.target.label}`, "Pair-level evidence", pair.final_status || "status unavailable"),
    );
    const split = element("section", undefined, "split");
    const summary = element("div", undefined, "panel");
    summary.append(element("h2", "Normalized fields"), recordList(pair));
    const artifacts = element("div", undefined, "panel");
    artifacts.append(element("h2", "Artifacts"));
    if (payload.artifacts.length) {
      const pre = element("pre", JSON.stringify(payload.artifacts, null, 2));
      artifacts.append(pre);
    } else {
      artifacts.append(element("p", "No pair-linked artifacts were exported.", "empty"));
    }
    split.append(summary, artifacts);
    root.append(split);
    app.replaceChildren(root);
  }

  function showError(error) {
    const message = error instanceof Error ? error.message : String(error);
    app.replaceChildren(element("p", `Unable to load this Atlas record: ${message}`, "error-box"));
  }

  async function start() {
    try {
      const response = await fetch("/release-config.json", { cache: "no-store" });
      if (!response.ok) throw new Error("release configuration is unavailable");
      config = await response.json();
      if (!TOKEN.test(config.release_token)) throw new Error("invalid release token");
      releaseLabel.textContent = config.release_id;
      brand.href = route();
      if (location.pathname === "/" || location.pathname === "/index.html") {
        history.replaceState(null, "", route());
      }
      const match = location.pathname.match(/^\/releases\/([^/]+)(?:\/(targets|drugs|pairs)\/([^/]+))?\/?$/);
      if (!match || match[1] !== config.release_token) throw new Error("route does not belong to this release");
      if (!match[2]) await renderRelease();
      else if (!TOKEN.test(match[3])) throw new Error("invalid record identifier");
      else if (match[2] === "pairs") await renderPair(match[3]);
      else await renderEntity(match[2], match[3]);
    } catch (error) {
      showError(error);
    }
  }

  start();
})();
"""


__all__ = [
    "APP_CSS",
    "APP_JS",
    "INDEX_HTML",
    "PUBLIC_HEADERS",
    "WORKER_SOURCE",
    "WRANGLER_TEMPLATE",
]
