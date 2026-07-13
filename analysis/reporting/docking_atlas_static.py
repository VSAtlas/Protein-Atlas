"""Generate a portable, read-only Docking Atlas explorer from release JSON.

This module intentionally performs no scientific computation.  Scores, ranks,
validity decisions, and qualification outcomes are rendered exactly as supplied
by an upstream release builder.
"""

from __future__ import annotations

import html
import json
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any
from urllib.parse import urlsplit


_STATUS_ORDER = {"valid": 0, "success": 0, "invalid": 1, "failed": 2, "missing": 3}


def generate_static_explorer(payload: Mapping[str, Any], output_dir: Path) -> dict[str, Any]:
    """Write overview, protein, drug, and pair pages under ``output_dir``.

    Accepted top-level collections are ``proteins``/``targets``,
    ``drugs``/``ligands``, and ``pairs``/``results``. Collections may be lists
    or ID-keyed objects. Pair records should carry ``protein_id`` or
    ``target_id`` and ``drug_id`` or ``ligand_id``. All scientific values are
    pass-through display fields.
    """

    output_dir = Path(output_dir)
    proteins = _records(payload.get("proteins", payload.get("targets", [])), "protein")
    drugs = _records(payload.get("drugs", payload.get("ligands", [])), "drug")
    pairs = _records(payload.get("pairs", payload.get("results", [])), "pair")
    release = _mapping(payload.get("release"))

    protein_by_id = {_record_id(item, "protein"): item for item in proteins}
    drug_by_id = {_record_id(item, "drug"): item for item in drugs}
    protein_slugs = _unique_slugs(protein_by_id)
    drug_slugs = _unique_slugs(drug_by_id)
    pair_ids, pair_slugs = _pair_identity(pairs)

    pairs_by_protein: dict[str, list[Mapping[str, Any]]] = {}
    pairs_by_drug: dict[str, list[Mapping[str, Any]]] = {}
    for pair in pairs:
        protein_id = _pair_ref(pair, "protein")
        drug_id = _pair_ref(pair, "drug")
        pairs_by_protein.setdefault(protein_id, []).append(pair)
        pairs_by_drug.setdefault(drug_id, []).append(pair)

    (output_dir / "assets").mkdir(parents=True, exist_ok=True)
    (output_dir / "proteins").mkdir(parents=True, exist_ok=True)
    (output_dir / "drugs").mkdir(parents=True, exist_ok=True)
    (output_dir / "pairs").mkdir(parents=True, exist_ok=True)
    (output_dir / "assets" / "atlas.css").write_text(_CSS, encoding="utf-8")
    (output_dir / "assets" / "atlas.js").write_text(_JS, encoding="utf-8")

    title = _text(release.get("title")) or "Docking Atlas"
    index_body = _overview_body(
        title,
        release,
        proteins,
        drugs,
        pairs,
        protein_slugs,
        drug_slugs,
    )
    (output_dir / "index.html").write_text(
        _document(title, index_body, depth=0, release=release), encoding="utf-8"
    )

    for protein_id, protein in protein_by_id.items():
        body = _protein_body(
            protein_id,
            protein,
            _display_pairs(pairs_by_protein.get(protein_id, []), "protein"),
            drug_by_id,
            drug_slugs,
            pair_ids,
            pair_slugs,
        )
        page_title = f"{_entity_name(protein, protein_id)} | {title}"
        (output_dir / "proteins" / f"{protein_slugs[protein_id]}.html").write_text(
            _document(page_title, body, depth=1, release=release), encoding="utf-8"
        )

    for drug_id, drug in drug_by_id.items():
        body = _drug_body(
            drug_id,
            drug,
            _display_pairs(pairs_by_drug.get(drug_id, []), "drug"),
            protein_by_id,
            protein_slugs,
            pair_ids,
            pair_slugs,
        )
        page_title = f"{_entity_name(drug, drug_id)} | {title}"
        (output_dir / "drugs" / f"{drug_slugs[drug_id]}.html").write_text(
            _document(page_title, body, depth=1, release=release), encoding="utf-8"
        )

    for pair, pair_id in zip(pairs, pair_ids):
        body = _pair_body(
            pair,
            pair_id,
            protein_by_id,
            drug_by_id,
            protein_slugs,
            drug_slugs,
        )
        (output_dir / "pairs" / f"{pair_slugs[pair_id]}.html").write_text(
            _document(f"Pair {pair_id} | {title}", body, depth=1, release=release),
            encoding="utf-8",
        )

    manifest = {
        "release_id": _text(release.get("id") or release.get("version")),
        "counts": {"proteins": len(proteins), "drugs": len(drugs), "pairs": len(pairs)},
        "entrypoint": "index.html",
    }
    (output_dir / "site_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def generate_static_explorer_from_json(release_json: Path, output_dir: Path) -> dict[str, Any]:
    """Load a JSON object and pass it to :func:`generate_static_explorer`."""

    with Path(release_json).open(encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, Mapping):
        raise ValueError("Docking Atlas release JSON must contain an object")
    return generate_static_explorer(payload, output_dir)


def _records(value: Any, kind: str) -> list[Mapping[str, Any]]:
    if isinstance(value, Mapping):
        records: list[Mapping[str, Any]] = []
        for key, record in value.items():
            if not isinstance(record, Mapping):
                raise ValueError(f"{kind} record {key!r} must be an object")
            enriched = dict(record)
            enriched.setdefault("id", str(key))
            records.append(enriched)
        return records
    if isinstance(value, list):
        if not all(isinstance(record, Mapping) for record in value):
            raise ValueError(f"{kind} records must be objects")
        return list(value)
    if value in (None, ""):
        return []
    raise ValueError(f"{kind} collection must be a list or object")


def _record_id(record: Mapping[str, Any], kind: str) -> str:
    candidates = {
        "protein": ("id", "protein_id", "target_id", "pdb_id"),
        "drug": ("id", "drug_id", "ligand_id", "name"),
        "pair": ("id", "pair_id"),
    }[kind]
    for key in candidates:
        value = _text(record.get(key))
        if value:
            return value
    raise ValueError(f"{kind} record is missing a stable identifier")


def _pair_ref(pair: Mapping[str, Any], kind: str) -> str:
    keys = ("protein_id", "target_id", "pdb_id") if kind == "protein" else ("drug_id", "ligand_id")
    for key in keys:
        value = _text(pair.get(key))
        if value:
            return value
    return ""


def _pair_identity(pairs: list[Mapping[str, Any]]) -> tuple[list[str], dict[str, str]]:
    ids: list[str] = []
    seen: dict[str, int] = {}
    for pair in pairs:
        explicit = _text(pair.get("id") or pair.get("pair_id"))
        base = explicit or f"{_pair_ref(pair, 'protein')}--{_pair_ref(pair, 'drug')}"
        count = seen.get(base, 0) + 1
        seen[base] = count
        ids.append(base if count == 1 else f"{base}--{count}")
    return ids, _unique_slugs({pair_id: {} for pair_id in ids})


def _unique_slugs(records: Mapping[str, Any]) -> dict[str, str]:
    result: dict[str, str] = {}
    used: set[str] = set()
    for record_id in records:
        base = re.sub(r"[^a-z0-9]+", "-", record_id.lower()).strip("-") or "record"
        slug = base
        counter = 2
        while slug in used:
            slug = f"{base}-{counter}"
            counter += 1
        used.add(slug)
        result[record_id] = slug
    return result


def _document(title: str, body: str, *, depth: int, release: Mapping[str, Any]) -> str:
    root = "../" * depth
    release_label = _text(release.get("id") or release.get("version")) or "unversioned"
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{_h(title)}</title><link rel="stylesheet" href="{root}assets/atlas.css"></head>
<body><header class="site-header"><a class="brand" href="{root}index.html">Docking Atlas</a>
<span class="release">Release {_h(release_label)}</span></header>
<main>{body}</main><footer>Atlas rankings describe computational docking results, not validated binding.</footer>
<script src="{root}assets/atlas.js"></script></body></html>"""


def _overview_body(
    title: str,
    release: Mapping[str, Any],
    proteins: list[Mapping[str, Any]],
    drugs: list[Mapping[str, Any]],
    pairs: list[Mapping[str, Any]],
    protein_slugs: Mapping[str, str],
    drug_slugs: Mapping[str, str],
) -> str:
    summary = _text(release.get("summary") or release.get("description"))
    protein_cards = "".join(
        _entity_card(item, _record_id(item, "protein"), f"proteins/{protein_slugs[_record_id(item, 'protein')]}.html", "protein")
        for item in proteins
    )
    drug_cards = "".join(
        _entity_card(item, _record_id(item, "drug"), f"drugs/{drug_slugs[_record_id(item, 'drug')]}.html", "drug")
        for item in drugs
    )
    valid_count = sum(_status(pair) in {"valid", "success"} for pair in pairs)
    invalid_count = sum(_status(pair) == "invalid" for pair in pairs)
    failed_count = len(pairs) - valid_count - invalid_count
    return f"""<section class="hero"><p class="eyebrow">Auditable docking release</p><h1>{_h(title)}</h1>
<p>{_h(summary or 'Explore receptor-qualified docking results in both directions: protein to drug and drug to protein.')}</p></section>
<section class="metrics"><div><strong>{len(proteins)}</strong><span>proteins</span></div><div><strong>{len(drugs)}</strong><span>drugs</span></div>
<div><strong>{len(pairs)}</strong><span>matrix cells</span></div><div><strong>{valid_count}</strong><span>valid</span></div>
<div><strong>{invalid_count}</strong><span>invalid</span></div><div><strong>{failed_count}</strong><span>failed / missing</span></div></section>
<section><div class="section-head"><h2>Protein explorer</h2><input data-filter="protein-list" type="search" placeholder="Search protein, PDB, gene…"></div>
<div id="protein-list" class="card-grid">{protein_cards or _empty('No proteins in this release.')}</div></section>
<section><div class="section-head"><h2>Drug explorer</h2><input data-filter="drug-list" type="search" placeholder="Search drug or ligand ID…"></div>
<div id="drug-list" class="card-grid">{drug_cards or _empty('No drugs in this release.')}</div></section>"""


def _entity_card(item: Mapping[str, Any], item_id: str, href: str, kind: str) -> str:
    name = _entity_name(item, item_id)
    subtitle_keys = ("gene", "pdb_id", "quality_status") if kind == "protein" else ("generic_name", "library", "approval_status")
    subtitle = " · ".join(_text(item.get(key)) for key in subtitle_keys if _text(item.get(key)))
    searchable = " ".join(_text(value) for value in item.values() if isinstance(value, (str, int, float)))
    return f'<a class="entity-card" href="{_h(href)}" data-search="{_h(searchable.lower())}"><h3>{_h(name)}</h3><p>{_h(subtitle or item_id)}</p></a>'


def _protein_body(
    protein_id: str,
    protein: Mapping[str, Any],
    pairs: list[Mapping[str, Any]],
    drugs: Mapping[str, Mapping[str, Any]],
    drug_slugs: Mapping[str, str],
    pair_ids: list[str],
    pair_slugs: Mapping[str, str],
) -> str:
    heading = _entity_name(protein, protein_id)
    meta = _definition_list(protein, exclude={"id", "name", "display_name", "description", "images"})
    rows = _pair_rows(pairs, pair_ids, pair_slugs, drugs, drug_slugs, "drug")
    return f"""{_breadcrumbs(('Proteins', None), (heading, None))}<section class="hero compact"><p class="eyebrow">Protein</p><h1>{_h(heading)}</h1>
<p>{_h(_text(protein.get('description')) or protein_id)}</p></section>{_qualification(protein)}
<section class="split"><div><h2>Receptor record</h2>{meta}</div>{_image_gallery(protein.get('images'))}</section>
<section><div class="section-head"><h2>Ranked drugs</h2><input data-table-filter="pair-table" type="search" placeholder="Filter drugs or status…"></div>
{_pair_table(rows, 'Drug', 'Protein-normalized rank')}</section>"""


def _drug_body(
    drug_id: str,
    drug: Mapping[str, Any],
    pairs: list[Mapping[str, Any]],
    proteins: Mapping[str, Mapping[str, Any]],
    protein_slugs: Mapping[str, str],
    pair_ids: list[str],
    pair_slugs: Mapping[str, str],
) -> str:
    heading = _entity_name(drug, drug_id)
    meta = _definition_list(drug, exclude={"id", "name", "display_name", "description", "images"})
    rows = _pair_rows(pairs, pair_ids, pair_slugs, proteins, protein_slugs, "protein")
    return f"""{_breadcrumbs(('Drugs', None), (heading, None))}<section class="hero compact"><p class="eyebrow">Drug</p><h1>{_h(heading)}</h1>
<p>{_h(_text(drug.get('description')) or drug_id)}</p></section>
<section class="split"><div><h2>Ligand record</h2>{meta}</div>{_image_gallery(drug.get('images'))}</section>
<section><div class="section-head"><h2>Ranked proteins</h2><input data-table-filter="pair-table" type="search" placeholder="Filter proteins or status…"></div>
{_pair_table(rows, 'Protein', 'Receptor-normalized rank')}</section>"""


def _display_pairs(pairs: list[Mapping[str, Any]], perspective: str) -> list[Mapping[str, Any]]:
    rank_keys = ("protein_rank", "rank_for_protein", "rank") if perspective == "protein" else ("drug_rank", "rank_for_drug", "cross_protein_rank")

    def key(pair: Mapping[str, Any]) -> tuple[Any, ...]:
        rank = next((_number(pair.get(name)) for name in rank_keys if _number(pair.get(name)) is not None), None)
        return (rank is None, rank if rank is not None else 0, _STATUS_ORDER.get(_status(pair), 9))

    return sorted(pairs, key=key)


def _pair_rows(
    pairs: list[Mapping[str, Any]],
    pair_ids: list[str],
    pair_slugs: Mapping[str, str],
    entities: Mapping[str, Mapping[str, Any]],
    entity_slugs: Mapping[str, str],
    entity_kind: str,
) -> list[str]:
    rows = []
    for pair in pairs:
        pair_id = _resolve_pair_id(pair, pair_ids)
        entity_id = _pair_ref(pair, entity_kind)
        entity = entities.get(entity_id, {})
        name = _entity_name(entity, entity_id or "Unresolved entity")
        entity_href = f"../{'drugs' if entity_kind == 'drug' else 'proteins'}/{entity_slugs.get(entity_id, '')}.html"
        rank = _first(pair, "protein_rank", "rank_for_protein", "drug_rank", "rank_for_drug", "cross_protein_rank", "rank")
        score = _score_summary(pair)
        status = _status(pair)
        rows.append(
            f'<tr data-search="{_h((name + " " + entity_id + " " + status).lower())}"><td><a href="{_h(entity_href)}">{_h(name)}</a><small>{_h(entity_id)}</small></td>'
            f'<td>{_h(_display(rank))}</td><td>{score}</td><td>{_badge(status)}</td>'
            f'<td><a class="button small" href="../pairs/{_h(pair_slugs[pair_id])}.html">Inspect</a></td></tr>'
        )
    return rows


def _resolve_pair_id(pair: Mapping[str, Any], pair_ids: list[str]) -> str:
    explicit = _text(pair.get("id") or pair.get("pair_id"))
    base = explicit or f"{_pair_ref(pair, 'protein')}--{_pair_ref(pair, 'drug')}"
    if base in pair_ids:
        return base
    matches = [pair_id for pair_id in pair_ids if pair_id.startswith(f"{base}--")]
    return matches[0] if matches else base


def _pair_table(rows: list[str], entity_label: str, rank_label: str) -> str:
    if not rows:
        return _empty("No pair records are available for this entity.")
    return f"""<div class="table-wrap"><table id="pair-table"><thead><tr><th>{_h(entity_label)}</th><th>{_h(rank_label)}</th>
<th>Atlas score</th><th>Status</th><th></th></tr></thead><tbody>{''.join(rows)}</tbody></table></div>"""


def _pair_body(
    pair: Mapping[str, Any],
    pair_id: str,
    proteins: Mapping[str, Mapping[str, Any]],
    drugs: Mapping[str, Mapping[str, Any]],
    protein_slugs: Mapping[str, str],
    drug_slugs: Mapping[str, str],
) -> str:
    protein_id = _pair_ref(pair, "protein")
    drug_id = _pair_ref(pair, "drug")
    protein_name = _entity_name(proteins.get(protein_id, {}), protein_id or "Unknown protein")
    drug_name = _entity_name(drugs.get(drug_id, {}), drug_id or "Unknown drug")
    status = _status(pair)
    scores = _mapping(pair.get("scores"))
    validity = _mapping(pair.get("validity") or pair.get("validation"))
    provenance = _mapping(pair.get("provenance"))
    failure = _text(pair.get("failure_reason") or validity.get("reason") or pair.get("invalid_reason"))
    status_note = f'<p class="callout {"danger" if status in {"failed", "invalid", "missing"} else "ok"}">{_badge(status)} {_h(failure or "No structured failure reason supplied.")}</p>'
    score_data = dict(scores)
    for key in ("atlas_score", "normalized_score", "z_score", "percentile_rank", "raw_score", "protein_rank", "drug_rank"):
        if key in pair and key not in score_data:
            score_data[key] = pair[key]
    reserved = {"id", "pair_id", "protein_id", "target_id", "pdb_id", "drug_id", "ligand_id", "scores", "validity", "validation", "provenance", "artifacts", "images", "failure_reason", "invalid_reason"}
    extra = {key: value for key, value in pair.items() if key not in reserved and key not in score_data}
    protein_href = f"../proteins/{protein_slugs.get(protein_id, '')}.html"
    drug_href = f"../drugs/{drug_slugs.get(drug_id, '')}.html"
    return f"""{_breadcrumbs(('Protein', protein_href), (protein_name, protein_href), ('Drug', drug_href), (drug_name, drug_href))}
<section class="hero compact"><p class="eyebrow">Protein–ligand result</p><h1>{_h(drug_name)} × {_h(protein_name)}</h1><p>Pair ID: <code>{_h(pair_id)}</code></p></section>
{status_note}{_image_gallery(pair.get('images'))}
<section class="three-col"><div><h2>Scores and ranks</h2>{_definition_list(score_data)}</div>
<div><h2>Validation</h2>{_definition_list(validity) or _empty('No validation fields supplied.')}</div>
<div><h2>Run fields</h2>{_definition_list(extra) or _empty('No additional run fields supplied.')}</div></section>
<section><h2>Pair-level provenance</h2>{_definition_list(provenance) or _empty('No provenance fields supplied.')}</section>
<section><h2>Artifacts and reconstructable inputs</h2>{_artifact_list(pair.get('artifacts'))}</section>"""


def _qualification(protein: Mapping[str, Any]) -> str:
    quality = _text(protein.get("quality_status") or protein.get("qualification_status"))
    redock = _text(protein.get("native_redock_status") or protein.get("redock_status"))
    rmsd = _display(protein.get("native_redock_rmsd") or protein.get("redock_rmsd"))
    if not any((quality, redock, rmsd)):
        return ""
    return f'<section class="qualification"><h2>Receptor qualification</h2><div>{_badge(quality or "not supplied")}<span>Quality policy</span></div><div>{_badge(redock or "not supplied")}<span>Native redock{": " + _h(rmsd) + " Å" if rmsd else ""}</span></div></section>'


def _score_summary(pair: Mapping[str, Any]) -> str:
    scores = _mapping(pair.get("scores"))
    value = _first(pair, "atlas_score", "normalized_score", "z_score")
    if value in (None, ""):
        value = _first(scores, "atlas_score", "normalized_score", "z_score", "selected")
    return _h(_display(value))


def _status(pair: Mapping[str, Any]) -> str:
    value = _text(pair.get("status") or pair.get("final_status")).lower()
    validity = _mapping(pair.get("validity") or pair.get("validation"))
    if not value and validity.get("valid") is not None:
        value = "valid" if bool(validity.get("valid")) else "invalid"
    return value or "missing"


def _definition_list(value: Any, exclude: set[str] | None = None) -> str:
    mapping = _mapping(value)
    rows = []
    for key, item in mapping.items():
        if exclude and key in exclude:
            continue
        if item in (None, "", [], {}):
            continue
        rows.append(f"<dt>{_h(_label(key))}</dt><dd>{_value_html(item)}</dd>")
    return f'<dl class="record">{"".join(rows)}</dl>' if rows else ""


def _value_html(value: Any) -> str:
    if isinstance(value, Mapping):
        return _definition_list(value)
    if isinstance(value, list):
        return "<ul>" + "".join(f"<li>{_value_html(item)}</li>" for item in value) + "</ul>"
    text = _display(value)
    if _safe_url(text):
        return f'<a href="{_h(text)}">{_h(text)}</a>'
    return _h(text)


def _image_gallery(value: Any) -> str:
    images = _iter_objects(value)
    cards = []
    for image in images:
        url = _text(image.get("url") or image.get("path") or image.get("href"))
        if not _safe_url(url):
            continue
        label = _text(image.get("label") or image.get("title")) or "Docked pose"
        alt = _text(image.get("alt")) or label
        cards.append(f'<figure><a href="{_h(url)}"><img loading="lazy" src="{_h(url)}" alt="{_h(alt)}"></a><figcaption>{_h(label)}</figcaption></figure>')
    return f'<div class="gallery">{"".join(cards)}</div>' if cards else ""


def _artifact_list(value: Any) -> str:
    artifacts = _iter_objects(value)
    rows = []
    for artifact in artifacts:
        url = _text(artifact.get("url") or artifact.get("path") or artifact.get("href"))
        if not _safe_url(url):
            continue
        label = _text(artifact.get("label") or artifact.get("name") or artifact.get("kind")) or url
        digest = _text(artifact.get("content_hash") or artifact.get("sha256") or artifact.get("digest"))
        rows.append(f'<li><a href="{_h(url)}">{_h(label)}</a>{f"<code>{_h(digest)}</code>" if digest else ""}</li>')
    return f'<ul class="artifacts">{"".join(rows)}</ul>' if rows else _empty("No downloadable artifacts supplied.")


def _breadcrumbs(*items: tuple[str, str | None]) -> str:
    parts = ['<a href="../index.html">Atlas</a>']
    for label, href in items:
        parts.append(f'<a href="{_h(href)}">{_h(label)}</a>' if href else _h(label))
    return '<nav class="breadcrumbs" aria-label="Breadcrumb">' + " / ".join(parts) + "</nav>"


def _badge(status: str) -> str:
    normalized = re.sub(r"[^a-z0-9]+", "-", status.lower()).strip("-")
    return f'<span class="badge status-{_h(normalized)}">{_h(status)}</span>'


def _entity_name(record: Mapping[str, Any], fallback: str) -> str:
    return _text(record.get("display_name") or record.get("name") or record.get("title")) or fallback


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _iter_objects(value: Any) -> list[Mapping[str, Any]]:
    if isinstance(value, Mapping):
        return [dict(item, label=item.get("label", key)) if isinstance(item, Mapping) else {"label": key, "url": item} for key, item in value.items()]
    if isinstance(value, list):
        return [item if isinstance(item, Mapping) else {"url": item} for item in value]
    return []


def _first(mapping: Mapping[str, Any], *keys: str) -> Any:
    return next((mapping[key] for key in keys if mapping.get(key) not in (None, "")), None)


def _number(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _text(value: Any) -> str:
    return str(value).strip() if value is not None else ""


def _display(value: Any) -> str:
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.6g}"
    return _text(value)


def _label(value: str) -> str:
    return value.replace("_", " ").strip().title()


def _safe_url(value: str) -> bool:
    if not value or any(char in value for char in ('"', "'", "<", ">", "\n", "\r")):
        return False
    parsed = urlsplit(value)
    return parsed.scheme in {"", "http", "https"} and not value.lower().startswith(("javascript:", "data:"))


def _h(value: Any) -> str:
    return html.escape(_text(value), quote=True)


def _empty(message: str) -> str:
    return f'<p class="empty">{_h(message)}</p>'


_CSS = r"""
:root{--ink:#17231b;--muted:#5c695f;--paper:#f6f7f2;--card:#fff;--line:#dce2d9;--green:#285b3b;--gold:#d39a2e;--red:#9f3939;--blue:#315d7c}*{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);font:16px/1.5 system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}a{color:var(--green)}.site-header{height:64px;padding:0 max(24px,calc((100vw - 1180px)/2));display:flex;align-items:center;justify-content:space-between;background:#183626;color:#fff}.brand{color:#fff;text-decoration:none;font-weight:800;font-size:1.15rem}.release{font-size:.85rem;opacity:.8}main{max-width:1180px;margin:auto;padding:36px 24px 80px}footer{padding:24px;text-align:center;background:#e8ece5;color:var(--muted);font-size:.85rem}.hero{padding:48px;border-radius:22px;background:linear-gradient(135deg,#e6efe5,#fff);border:1px solid var(--line);margin-bottom:24px}.hero.compact{padding:30px}.hero h1{font-size:clamp(2rem,5vw,4.2rem);line-height:1.03;margin:.2em 0}.hero.compact h1{font-size:clamp(1.8rem,4vw,3rem)}.hero p{max-width:780px;color:#405047}.eyebrow{text-transform:uppercase;letter-spacing:.12em;font-weight:800;font-size:.75rem;color:var(--green)!important}.metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(130px,1fr));gap:12px;margin:24px 0 50px}.metrics div,.qualification div{background:var(--card);border:1px solid var(--line);border-radius:14px;padding:18px}.metrics strong{display:block;font-size:2rem}.metrics span,.qualification span{color:var(--muted);font-size:.8rem}.section-head{display:flex;gap:18px;align-items:center;justify-content:space-between;margin-top:46px}.section-head input{min-width:280px;padding:11px 14px;border:1px solid #abb8ac;border-radius:10px;background:#fff}.card-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(230px,1fr));gap:14px}.entity-card{display:block;padding:20px;background:var(--card);border:1px solid var(--line);border-radius:14px;text-decoration:none;color:inherit;transition:.15s}.entity-card:hover{transform:translateY(-2px);border-color:#8da393;box-shadow:0 8px 20px #18362614}.entity-card h3{margin:0}.entity-card p{color:var(--muted);margin:.4rem 0 0}.breadcrumbs{color:var(--muted);font-size:.85rem;margin-bottom:18px}.split{display:grid;grid-template-columns:minmax(280px,1fr) minmax(320px,1.4fr);gap:28px}.three-col{display:grid;grid-template-columns:repeat(3,1fr);gap:20px}.record{display:grid;grid-template-columns:minmax(130px,.7fr) 1.5fr;margin:0}.record dt,.record dd{padding:9px 0;border-bottom:1px solid var(--line);overflow-wrap:anywhere}.record dt{color:var(--muted);font-size:.82rem;font-weight:700}.record dd{margin:0}.record .record{grid-column:1/-1}.gallery{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:12px;margin:20px 0}.gallery figure{margin:0;background:#fff;border:1px solid var(--line);border-radius:14px;overflow:hidden}.gallery img{display:block;width:100%;aspect-ratio:4/3;object-fit:cover;background:#e5e8e3}.gallery figcaption{padding:8px 12px;color:var(--muted);font-size:.8rem}.table-wrap{overflow:auto;background:#fff;border:1px solid var(--line);border-radius:14px}table{width:100%;border-collapse:collapse}th,td{text-align:left;padding:13px;border-bottom:1px solid var(--line)}th{font-size:.75rem;text-transform:uppercase;letter-spacing:.05em;color:var(--muted)}td small{display:block;color:var(--muted)}.badge{display:inline-block;border-radius:999px;background:#e8ebe7;padding:3px 9px;font-size:.75rem;font-weight:750}.status-valid,.status-success,.status-qualified,.status-pass{background:#d8eddf;color:#215234}.status-invalid,.status-warning{background:#fff0c9;color:#72500d}.status-failed,.status-missing,.status-fail{background:#f5dada;color:#792727}.button{display:inline-block;padding:8px 13px;border-radius:8px;background:var(--green);color:#fff;text-decoration:none}.button.small{padding:5px 10px;font-size:.8rem}.qualification{display:flex;gap:12px;align-items:stretch;margin:20px 0}.qualification h2{margin-right:auto}.qualification div{min-width:150px}.qualification span{display:block;margin-top:6px}.callout{padding:14px 18px;border-radius:12px}.callout.ok{background:#dff0e4}.callout.danger{background:#f5e1d9}.artifacts{padding:0;list-style:none}.artifacts li{display:flex;justify-content:space-between;gap:20px;padding:10px;border-bottom:1px solid var(--line)}code{font-size:.82em;overflow-wrap:anywhere}.empty{color:var(--muted);font-style:italic;padding:18px;background:#fff;border:1px dashed #bdc7bd;border-radius:10px}@media(max-width:760px){main{padding:20px 14px 60px}.hero{padding:28px 20px}.split,.three-col{grid-template-columns:1fr}.section-head{align-items:stretch;flex-direction:column}.section-head input{min-width:0;width:100%}.qualification{display:grid}.record{grid-template-columns:1fr}.record dd{padding-top:0}}
"""


_JS = r"""
document.querySelectorAll('[data-filter]').forEach(function(input){input.addEventListener('input',function(){var q=input.value.trim().toLowerCase();var root=document.getElementById(input.dataset.filter);if(!root)return;root.querySelectorAll('[data-search]').forEach(function(item){item.hidden=q && !item.dataset.search.includes(q);});});});
document.querySelectorAll('[data-table-filter]').forEach(function(input){input.addEventListener('input',function(){var q=input.value.trim().toLowerCase();var root=document.getElementById(input.dataset.tableFilter);if(!root)return;root.querySelectorAll('tbody tr').forEach(function(row){row.hidden=q && !row.dataset.search.includes(q);});});});
"""
