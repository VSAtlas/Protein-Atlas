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
_PAGE_SIZE = 100


def generate_static_explorer(
    payload: Mapping[str, Any], output_dir: Path
) -> dict[str, Any]:
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
    artifacts = _records(payload.get("artifacts", []), "artifact")
    release = _mapping(payload.get("release"))

    protein_by_id = {_record_id(item, "protein"): item for item in proteins}
    drug_by_id = {_record_id(item, "drug"): item for item in drugs}
    protein_slugs = _unique_slugs(protein_by_id)
    drug_slugs = _unique_slugs(drug_by_id)
    pair_ids, pair_slugs = _pair_identity(pairs)
    pair_id_by_object = {id(pair): pair_id for pair, pair_id in zip(pairs, pair_ids)}
    artifacts_by_pair: dict[str, list[Mapping[str, Any]]] = {}
    for artifact in artifacts:
        pair_reference = _text(artifact.get("pair_cell_id") or artifact.get("pair_id"))
        if pair_reference:
            artifacts_by_pair.setdefault(pair_reference, []).append(artifact)

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
    analysis_payload = _analysis_payload(
        payload, pairs, proteins=protein_by_id, drugs=drug_by_id
    )
    (output_dir / "assets" / "analysis_data.json").write_text(
        json.dumps(analysis_payload, separators=(",", ":"), ensure_ascii=True),
        encoding="utf-8",
    )

    title = _text(release.get("title")) or "Docking Atlas"
    index_body = _overview_body(
        title,
        release,
        proteins,
        drugs,
        pairs,
        protein_slugs,
        drug_slugs,
        _mapping(payload.get("coverage")),
        payload.get("downloads", release.get("downloads")),
    )
    (output_dir / "index.html").write_text(
        _document(title, index_body, depth=0, release=release), encoding="utf-8"
    )
    (output_dir / "analysis.html").write_text(
        _document(
            f"Guided analysis | {title}", _analysis_body(), depth=0, release=release
        ),
        encoding="utf-8",
    )

    for protein_id, protein in protein_by_id.items():
        body = _protein_body(
            protein_id,
            protein,
            _display_pairs(pairs_by_protein.get(protein_id, []), "protein"),
            drug_by_id,
            drug_slugs,
            pair_id_by_object,
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
            pair_id_by_object,
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
            artifacts_by_pair.get(_pair_artifact_ref(pair), []),
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


def generate_static_explorer_from_json(
    release_json: Path, output_dir: Path
) -> dict[str, Any]:
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
    keys = (
        ("protein_id", "target_id", "pdb_id")
        if kind == "protein"
        else ("drug_id", "ligand_id")
    )
    for key in keys:
        value = _text(pair.get(key))
        if value:
            return value

    return ""


def _pair_artifact_ref(pair: Mapping[str, Any]) -> str:
    return _text(pair.get("pair_cell_id") or pair.get("id") or pair.get("pair_id"))


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
<body><header class="site-header"><a class="brand" href="{root}index.html">Docking Atlas</a><nav><a href="{root}analysis.html">Guided analysis</a></nav>
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
    coverage: Mapping[str, Any],
    downloads: Any,
) -> str:
    summary = _text(release.get("summary") or release.get("description"))
    protein_cards = [
        _entity_card(
            item,
            _record_id(item, "protein"),
            f"proteins/{protein_slugs[_record_id(item, 'protein')]}.html",
            "protein",
        )
        for item in proteins
    ]
    drug_cards = [
        _entity_card(
            item,
            _record_id(item, "drug"),
            f"drugs/{drug_slugs[_record_id(item, 'drug')]}.html",
            "drug",
        )
        for item in drugs
    ]
    valid_count = sum(_status(pair) in {"valid", "success"} for pair in pairs)
    invalid_count = sum(_status(pair) == "invalid" for pair in pairs)
    failed_count = len(pairs) - valid_count - invalid_count
    return f"""<section class="hero"><p class="eyebrow">Auditable docking release</p><h1>{_h(title)}</h1>
<p>{_h(summary or "Explore receptor-qualified docking results in both directions: protein to drug and drug to protein.")}</p></section>
<section class="metrics"><div><strong>{len(proteins)}</strong><span>proteins</span></div><div><strong>{len(drugs)}</strong><span>drugs</span></div>
<div><strong>{len(pairs)}</strong><span>matrix cells</span></div><div><strong>{valid_count}</strong><span>valid</span></div>
<div><strong>{invalid_count}</strong><span>invalid</span></div><div><strong>{failed_count}</strong><span>failed / missing</span></div></section>
{_coverage_callout(coverage)}
{_download_section(downloads)}
<section><div class="section-head"><h2>Protein explorer</h2><input data-filter="protein-list" type="search" placeholder="Search protein, PDB, gene…"></div>
{_paged_card_grid(protein_cards, "protein-list", "No proteins in this release.")}</section>
<section><div class="section-head"><h2>Drug explorer</h2><input data-filter="drug-list" type="search" placeholder="Search drug or ligand ID…"></div>
{_paged_card_grid(drug_cards, "drug-list", "No drugs in this release.")}</section>"""


def _paged_card_grid(cards: list[str], element_id: str, empty_message: str) -> str:
    if not cards:
        return f'<div id="{_h(element_id)}" class="card-grid">{_empty(empty_message)}</div>'
    data = [
        {"html": card, "search": html.unescape(_html_data_value(card, "search"))}
        for card in cards
    ]
    initial = "".join(cards[:_PAGE_SIZE])
    payload = _json_script(data)
    return (
        f'<div id="{_h(element_id)}" class="card-grid">{initial}</div>'
        f'<div class="pager" data-card-pager="{_h(element_id)}">'
        '<button type="button" data-page-prev>Previous</button>'
        '<span data-page-status></span><button type="button" data-page-next>Next</button></div>'
        f'<script type="application/json" data-card-source="{_h(element_id)}">{payload}</script>'
    )


def _coverage_callout(coverage: Mapping[str, Any]) -> str:
    if not coverage:
        return ""
    count = _number(coverage.get("primary_score_count", coverage.get("count")))
    pair_count = _number(coverage.get("pair_count"))
    fraction = _number(coverage.get("primary_score_fraction", coverage.get("fraction")))
    if fraction is None and count is not None and pair_count:
        fraction = count / pair_count
    incomplete = (fraction is not None and fraction < 1.0) or (
        count is not None and pair_count is not None and count < pair_count
    )
    count_label = _display(int(count)) if count is not None else "not supplied"
    pair_label = _display(int(pair_count)) if pair_count is not None else "not supplied"
    fraction_label = (
        f"{fraction:.1%}" if fraction is not None else "fraction not supplied"
    )
    note = (
        "Secondary scores are not substituted for missing primary scores."
        if incomplete
        else "Primary-score coverage is complete for the declared release matrix."
    )
    return f'<p class="callout {"danger" if incomplete else "ok"}"><strong>Primary final_score coverage:</strong> {_h(count_label)} / {_h(pair_label)} pairs ({_h(fraction_label)}). {_h(note)}</p>'


def _download_section(downloads: Any) -> str:
    if not downloads:
        return ""
    return (
        '<section class="release-downloads"><h2>Download release</h2>'
        '<p class="muted">Frozen release files and content hashes.</p>'
        f"{_artifact_list(downloads)}</section>"
    )


def _analysis_payload(
    payload: Mapping[str, Any],
    pairs: list[Mapping[str, Any]],
    *,
    proteins: Mapping[str, Mapping[str, Any]],
    drugs: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """Build a browser-safe, allowlisted analysis dataset from release fields."""

    config = _mapping(payload.get("analysis"))
    base_dimensions = [
        "target_id",
        "protein_name",
        "drug_id",
        "drug_name",
        "status",
        "failure_reason",
        "quality_status",
        "native_redock_status",
        "rank_eligible",
        "ranking_eligibility_reason",
        "apo_exploratory_rank_eligible",
        "apo_exploratory_ranking_eligibility_reason",
        "ranking_track",
        "final_score_source",
        "final_score_source_effective",
        "final_score_source_family",
    ]
    base_measures = [
        "final_score",
        "normalized_score",
        "atlas_score",
        "z_score",
        "raw_score",
        "percentile_rank",
    ]
    dimension_keys = _allowed_analysis_keys(
        base_dimensions, config.get("dimensions", [])
    )
    measure_keys = _allowed_analysis_keys(base_measures, config.get("measures", []))
    rows: list[dict[str, Any]] = []
    for pair in pairs:
        protein_id = _pair_ref(pair, "protein")
        drug_id = _pair_ref(pair, "drug")
        protein = proteins.get(protein_id, {})
        drug = drugs.get(drug_id, {})
        scores = _mapping(pair.get("scores"))
        row: dict[str, Any] = {
            "target_id": protein_id,
            "protein_name": _entity_name(protein, protein_id),
            "drug_id": drug_id,
            "drug_name": _entity_name(drug, drug_id),
            "status": _status(pair),
            "failure_reason": _text(
                pair.get("failure_reason") or pair.get("invalid_reason")
            ),
            "quality_status": _text(
                protein.get("quality_status") or protein.get("qualification_status")
            ),
            "native_redock_status": _text(
                protein.get("native_redock_status") or protein.get("redock_status")
            ),
        }
        for key in dimension_keys:
            if key not in row:
                value = pair.get(key)
                if isinstance(value, (str, int, float, bool)):
                    row[key] = value
        for key in measure_keys:
            row[key] = _number(pair.get(key, scores.get(key)))
        rows.append(row)
    dimensions = [
        key
        for key in dimension_keys
        if any(row.get(key) not in (None, "") for row in rows)
    ]
    measures = [
        key for key in measure_keys if any(row.get(key) is not None for row in rows)
    ]
    duckdb = _mapping(config.get("duckdb_wasm"))
    return {
        "schema_version": 1,
        "dimensions": dimensions,
        "measures": measures,
        "rows": rows,
        "optional_engine": {
            "name": "duckdb-wasm",
            "enabled": bool(duckdb.get("enabled", False)),
            "module_url": _text(duckdb.get("module_url")),
            "worker_url": _text(duckdb.get("worker_url")),
            "note": "Optional hook only; this static release never loads CDN code automatically.",
        },
    }


def _allowed_analysis_keys(defaults: list[str], requested: Any) -> list[str]:
    keys = list(defaults)
    if isinstance(requested, list):
        for raw in requested:
            key = _text(raw)
            if re.fullmatch(r"[A-Za-z][A-Za-z0-9_]{0,63}", key) and key not in keys:
                keys.append(key)
    return keys


def _analysis_body() -> str:
    return """<section class="hero compact"><p class="eyebrow">Guided browser analysis</p>
<h1>Explore the released matrix</h1><p>Filter and summarize supplied release values locally in your browser. These controls do not recalculate docking scores or receptor normalization.</p></section>
<noscript><p class="callout danger">JavaScript is required for the guided analysis page. Entity and pair pages remain available without it.</p></noscript>
<section id="analysis-app" class="analysis-app" data-source="assets/analysis_data.json">
<div id="analysis-message" class="callout">Loading release data…</div>
<div id="analysis-controls" class="analysis-controls" hidden>
<label>Filter dimension<select id="analysis-filter-dimension"></select></label>
<label>Filter value<select id="analysis-filter-value"><option value="">All values</option></select></label>
<label>Group by<select id="analysis-group"></select></label>
<label>Measure<select id="analysis-measure"></select></label>
<label>Aggregation<select id="analysis-aggregation"><option value="count">Count</option><option value="mean">Mean</option><option value="min">Minimum</option><option value="max">Maximum</option></select></label>
<button id="analysis-reset" type="button">Reset</button><button id="analysis-download" type="button">Download current CSV</button>
</div>
<div class="analysis-tabs" role="tablist"><button data-analysis-view="table" class="is-active">Table</button><button data-analysis-view="histogram">Histogram</button><button data-analysis-view="matrix">Matrix</button></div>
<p id="analysis-summary" class="muted"></p><div id="analysis-result" class="analysis-result"></div>
</section>
<section><h2>Optional large-release engine</h2><p class="muted">The baseline uses dependency-free browser JavaScript. A release may declare explicit DuckDB-Wasm module and worker URLs; no CDN or remote script is loaded silently.</p></section>"""


def _entity_card(item: Mapping[str, Any], item_id: str, href: str, kind: str) -> str:
    name = _entity_name(item, item_id)
    subtitle_keys = (
        ("gene", "pdb_id", "quality_status")
        if kind == "protein"
        else ("generic_name", "library", "approval_status")
    )
    subtitle = " · ".join(
        _text(item.get(key)) for key in subtitle_keys if _text(item.get(key))
    )
    searchable = " ".join(
        _text(value) for value in item.values() if isinstance(value, (str, int, float))
    )
    return f'<a class="entity-card" href="{_h(href)}" data-search="{_h(searchable.lower())}"><h3>{_h(name)}</h3><p>{_h(subtitle or item_id)}</p></a>'


def _protein_body(
    protein_id: str,
    protein: Mapping[str, Any],
    pairs: list[Mapping[str, Any]],
    drugs: Mapping[str, Mapping[str, Any]],
    drug_slugs: Mapping[str, str],
    pair_id_by_object: Mapping[int, str],
    pair_slugs: Mapping[str, str],
) -> str:
    heading = _entity_name(protein, protein_id)
    meta = _definition_list(
        protein, exclude={"id", "name", "display_name", "description", "images"}
    )
    rankings = _ranking_sections(
        pairs,
        pair_id_by_object,
        pair_slugs,
        drugs,
        drug_slugs,
        "drug",
        "Drug",
        "Protein-normalized rank",
    )
    return f"""{_breadcrumbs(("Proteins", None), (heading, None))}<section class="hero compact"><p class="eyebrow">Protein</p><h1>{_h(heading)}</h1>
<p>{_h(_text(protein.get("description")) or protein_id)}</p></section>{_qualification(protein)}
<section class="split"><div><h2>Receptor record</h2>{meta}</div>{_image_gallery(protein.get("images"))}</section>
{rankings}"""


def _drug_body(
    drug_id: str,
    drug: Mapping[str, Any],
    pairs: list[Mapping[str, Any]],
    proteins: Mapping[str, Mapping[str, Any]],
    protein_slugs: Mapping[str, str],
    pair_id_by_object: Mapping[int, str],
    pair_slugs: Mapping[str, str],
) -> str:
    heading = _entity_name(drug, drug_id)
    meta = _definition_list(
        drug, exclude={"id", "name", "display_name", "description", "images"}
    )
    rankings = _ranking_sections(
        pairs,
        pair_id_by_object,
        pair_slugs,
        proteins,
        protein_slugs,
        "protein",
        "Protein",
        "Receptor-normalized rank",
    )
    return f"""{_breadcrumbs(("Drugs", None), (heading, None))}<section class="hero compact"><p class="eyebrow">Drug</p><h1>{_h(heading)}</h1>
<p>{_h(_text(drug.get("description")) or drug_id)}</p></section>
<section class="split"><div><h2>Ligand record</h2>{meta}</div>{_image_gallery(drug.get("images"))}</section>
{rankings}"""


def _truthy(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return value == 1
    return _text(value).casefold() in {"1", "true", "yes"}


def _eligible_ranking_track(pair: Mapping[str, Any]) -> str | None:
    track = _text(pair.get("ranking_track"))
    if track == "qualified_holo" and _truthy(pair.get("rank_eligible")):
        return track
    if track == "exploratory_apo" and _truthy(
        pair.get("apo_exploratory_rank_eligible")
    ):
        return track
    return None


def _rank_keys(pair: Mapping[str, Any], perspective: str) -> tuple[str, ...]:
    track = _eligible_ranking_track(pair)
    if track == "exploratory_apo":
        return (
            ("apo_rank_within_receptor",)
            if perspective == "protein"
            else ("apo_rank_across_receptors",)
        )
    if track == "qualified_holo":
        return (
            ("rank_within_receptor", "protein_rank", "rank_for_protein", "rank")
            if perspective == "protein"
            else (
                "rank_across_receptors",
                "drug_rank",
                "rank_for_drug",
                "cross_protein_rank",
                "rank",
            )
        )
    return ()


def _display_pairs(
    pairs: list[Mapping[str, Any]], perspective: str
) -> list[Mapping[str, Any]]:
    def key(pair: Mapping[str, Any]) -> tuple[Any, ...]:
        rank = next(
            (
                _number(pair.get(name))
                for name in _rank_keys(pair, perspective)
                if _number(pair.get(name)) is not None
            ),
            None,
        )
        return (
            rank is None,
            rank if rank is not None else 0,
            _STATUS_ORDER.get(_status(pair), 9),
        )

    return sorted(pairs, key=key)


def _pair_rows(
    pairs: list[Mapping[str, Any]],
    pair_id_by_object: Mapping[int, str],
    pair_slugs: Mapping[str, str],
    entities: Mapping[str, Mapping[str, Any]],
    entity_slugs: Mapping[str, str],
    entity_kind: str,
    *,
    mark_ineligible: bool = True,
) -> list[str]:
    rows = []
    for pair in pairs:
        pair_id = pair_id_by_object[id(pair)]
        entity_id = _pair_ref(pair, entity_kind)
        entity = entities.get(entity_id, {})
        name = _entity_name(entity, entity_id or "Unresolved entity")
        entity_href = f"../{'drugs' if entity_kind == 'drug' else 'proteins'}/{entity_slugs.get(entity_id, '')}.html"
        perspective = "protein" if entity_kind == "drug" else "drug"
        rank = _first(pair, *_rank_keys(pair, perspective))
        score = _score_summary(pair)
        status = _status(pair)
        track = _eligible_ranking_track(pair)
        raw_track = _text(pair.get("ranking_track"))
        rank_ineligible = track is None
        if track == "qualified_holo":
            track_label = "HOLO qualified"
        elif track == "exploratory_apo":
            track_label = "APO exploratory"
        elif raw_track == "qualified_holo":
            track_label = "HOLO excluded"
        elif raw_track == "exploratory_apo":
            track_label = "APO exploratory excluded"
        else:
            track_label = "Unqualified"
        visibility = (
            ' data-rank-eligible="false"' if rank_ineligible and mark_ineligible else ""
        )
        rows.append(
            f'<tr data-search="{_h((name + " " + entity_id + " " + status).lower())}"{visibility}><td><a href="{_h(entity_href)}">{_h(name)}</a><small>{_h(entity_id)}</small></td>'
            f"<td>{_h(_display(rank))}<small>{_h(track_label)}</small></td>"
            f"<td>{score}</td><td>{_badge(status)}</td>"
            f'<td><a class="button small" href="../pairs/{_h(pair_slugs[pair_id])}.html">Inspect</a></td></tr>'
        )
    return rows


def _ranking_sections(
    pairs: list[Mapping[str, Any]],
    pair_id_by_object: Mapping[int, str],
    pair_slugs: Mapping[str, str],
    entities: Mapping[str, Mapping[str, Any]],
    entity_slugs: Mapping[str, str],
    entity_kind: str,
    entity_label: str,
    rank_label: str,
) -> str:
    perspective = "protein" if entity_kind == "drug" else "drug"
    ordered = _display_pairs(pairs, perspective)
    holo = [
        pair for pair in ordered if _eligible_ranking_track(pair) == "qualified_holo"
    ]
    apo = [
        pair for pair in ordered if _eligible_ranking_track(pair) == "exploratory_apo"
    ]
    excluded = [pair for pair in ordered if _eligible_ranking_track(pair) is None]

    def rows_for(
        selected: list[Mapping[str, Any]], *, mark_ineligible: bool = True
    ) -> list[str]:
        return _pair_rows(
            selected,
            pair_id_by_object,
            pair_slugs,
            entities,
            entity_slugs,
            entity_kind,
            mark_ineligible=mark_ineligible,
        )

    return f"""
<section><div class="section-head"><div><h2>Qualified HOLO rankings</h2><p class="muted">Receptor-qualified HOLO results only; APO contexts are never pooled into these ranks.</p></div><input data-table-filter="qualified-holo-table" type="search" placeholder="Filter {_h(entity_label.lower())} or status…"></div>
{_pair_table(rows_for(holo), entity_label, rank_label, "qualified-holo-table")}</section>
<section><div class="section-head"><div><h2>Exploratory APO rankings</h2><p class="muted">Separate APO-only exploratory results; these ranks are not interchangeable with qualified HOLO rankings.</p></div><input data-table-filter="exploratory-apo-table" type="search" placeholder="Filter {_h(entity_label.lower())} or status…"></div>
{_pair_table(rows_for(apo), entity_label, rank_label, "exploratory-apo-table")}</section>
<section><div class="section-head"><div><h2>Unranked and excluded pair cells</h2><p class="muted">Failure-complete cells that satisfy neither ranking contract remain visible without a rank.</p></div><input data-table-filter="unranked-pair-table" type="search" placeholder="Filter {_h(entity_label.lower())} or status…"></div>
{_pair_table(rows_for(excluded, mark_ineligible=False), entity_label, "Rank (not assigned)", "unranked-pair-table")}</section>"""


def _pair_table(
    rows: list[str], entity_label: str, rank_label: str, table_id: str
) -> str:
    if not rows:
        return _empty("No pair records are available for this entity.")
    data = [
        {
            "html": row,
            "search": html.unescape(_html_data_value(row, "search")),
            "eligible": 'data-rank-eligible="false"' not in row,
        }
        for row in rows
    ]
    initial_rows = [row for row in rows if 'data-rank-eligible="false"' not in row][
        :_PAGE_SIZE
    ]
    initial = "".join(initial_rows)
    if not initial:
        initial = (
            '<tr><td colspan="5" class="muted">No eligible rows. Use “Show invalid, '
            "failed, and unranked cells” to inspect excluded cells.</td></tr>"
        )
    payload = _json_script(data)
    return f"""<div class="table-wrap"><table id="{_h(table_id)}"><thead><tr><th>{_h(entity_label)}</th><th>{_h(rank_label)}</th>
<th>Final score</th><th>Status</th><th></th></tr></thead><tbody>{initial}</tbody></table></div>
<div class="pager" data-pair-pager="{_h(table_id)}"><button type="button" data-page-prev>Previous</button><span data-page-status></span><button type="button" data-page-next>Next</button></div>
<script type="application/json" data-pair-source="{_h(table_id)}">{payload}</script>"""


def _html_data_value(markup: str, name: str) -> str:
    match = re.search(rf'\bdata-{re.escape(name)}="([^"]*)"', markup)
    return match.group(1) if match else ""


def _json_script(value: Any) -> str:
    return (
        json.dumps(value, separators=(",", ":"), ensure_ascii=True)
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("&", "\\u0026")
    )


def _pair_body(
    pair: Mapping[str, Any],
    pair_id: str,
    proteins: Mapping[str, Mapping[str, Any]],
    drugs: Mapping[str, Mapping[str, Any]],
    protein_slugs: Mapping[str, str],
    drug_slugs: Mapping[str, str],
    linked_artifacts: list[Mapping[str, Any]],
) -> str:
    protein_id = _pair_ref(pair, "protein")
    drug_id = _pair_ref(pair, "drug")
    protein_name = _entity_name(
        proteins.get(protein_id, {}), protein_id or "Unknown protein"
    )
    drug_name = _entity_name(drugs.get(drug_id, {}), drug_id or "Unknown drug")
    status = _status(pair)
    scores = _mapping(pair.get("scores"))
    validity = _mapping(pair.get("validity") or pair.get("validation"))
    artifact_value = pair.get("artifacts") or linked_artifacts
    provenance = _mapping(pair.get("provenance"))
    failure = _text(
        pair.get("failure_reason")
        or validity.get("reason")
        or pair.get("invalid_reason")
    )
    status_note = f'<p class="callout {"danger" if status not in {"valid", "success"} else "ok"}">{_badge(status)} {_h(failure or "No structured failure reason supplied.")}</p>'
    score_data = dict(scores)
    for key in (
        "final_score",
        "final_score_source",
        "final_score_source_effective",
        "final_score_source_family",
        "atlas_score",
        "normalized_score",
        "z_score",
        "percentile_rank",
        "raw_score",
        "protein_rank",
        "drug_rank",
    ):
        if key in pair and key not in score_data:
            score_data[key] = pair[key]
    reserved = {
        "id",
        "pair_id",
        "protein_id",
        "target_id",
        "pdb_id",
        "drug_id",
        "ligand_id",
        "scores",
        "validity",
        "validation",
        "provenance",
        "artifacts",
        "images",
        "failure_reason",
        "invalid_reason",
    }
    extra = {
        key: value
        for key, value in pair.items()
        if key not in reserved and key not in score_data
    }
    protein_href = f"../proteins/{protein_slugs.get(protein_id, '')}.html"
    drug_href = f"../drugs/{drug_slugs.get(drug_id, '')}.html"
    return f"""{_breadcrumbs(("Protein", protein_href), (protein_name, protein_href), ("Drug", drug_href), (drug_name, drug_href))}
<section class="hero compact"><p class="eyebrow">Protein–ligand result</p><h1>{_h(drug_name)} × {_h(protein_name)}</h1><p>Pair ID: <code>{_h(pair_id)}</code></p></section>
{status_note}{_image_gallery(pair.get("images"))}
<section class="three-col"><div><h2>Scores and ranks</h2>{_definition_list(score_data)}</div>
<div><h2>Validation</h2>{_definition_list(validity) or _empty("No validation fields supplied.")}</div>
<div><h2>Run fields</h2>{_definition_list(extra) or _empty("No additional run fields supplied.")}</div></section>
<section><h2>Pair-level provenance</h2>{_definition_list(provenance) or _empty("No provenance fields supplied.")}</section>
<section><h2>Artifacts and reconstructable inputs</h2>{_artifact_list(artifact_value)}</section>"""


def _qualification(protein: Mapping[str, Any]) -> str:
    quality = _text(
        protein.get("quality_status") or protein.get("qualification_status")
    )
    redock = _text(protein.get("native_redock_status") or protein.get("redock_status"))
    rmsd = _display(protein.get("native_redock_rmsd") or protein.get("redock_rmsd"))
    if not any((quality, redock, rmsd)):
        return ""
    return f'<section class="qualification"><h2>Receptor qualification</h2><div>{_badge(quality or "not supplied")}<span>Quality policy</span></div><div>{_badge(redock or "not supplied")}<span>Native redock{": " + _h(rmsd) + " Å" if rmsd else ""}</span></div></section>'


def _score_summary(pair: Mapping[str, Any]) -> str:
    scores = _mapping(pair.get("scores"))
    value = pair.get("final_score")
    if value in (None, ""):
        value = scores.get("final_score")
    source = _text(
        pair.get("final_score_source_effective") or pair.get("final_score_source")
    )
    family = _text(pair.get("final_score_source_family"))
    details = "".join(f"<small>{_h(item)}</small>" for item in (source, family) if item)
    return f"{_h(_display(value))}{details}"


def _status(pair: Mapping[str, Any]) -> str:
    value = _text(pair.get("status") or pair.get("final_status")).lower()
    aliases = {
        "ok": "valid",
        "completed": "valid",
        "unsuccessful": "failed",
        "error": "failed",
    }
    value = aliases.get(value, value)
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
        return (
            "<ul>"
            + "".join(f"<li>{_value_html(item)}</li>" for item in value)
            + "</ul>"
        )
    text = _display(value)
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
        cards.append(
            f'<figure><a href="{_h(url)}"><img loading="lazy" src="{_h(url)}" alt="{_h(alt)}"></a><figcaption>{_h(label)}</figcaption></figure>'
        )
    return f'<div class="gallery">{"".join(cards)}</div>' if cards else ""


def _artifact_list(value: Any) -> str:
    artifacts = _iter_objects(value)
    rows = []
    for artifact in artifacts:
        url = _text(artifact.get("url") or artifact.get("path") or artifact.get("href"))
        label = _text(
            artifact.get("label")
            or artifact.get("name")
            or artifact.get("kind")
            or artifact.get("member_name")
            or artifact.get("file_type")
            or artifact.get("stage")
        )
        digest = _text(
            artifact.get("content_hash")
            or artifact.get("sha256")
            or artifact.get("digest")
        )
        label = (
            label
            or url
            or (
                f"Artifact {artifact['artifact_id']}"
                if artifact.get("artifact_id") is not None
                else ""
            )
        )
        if not label and not digest:
            continue
        content = (
            f'<a href="{_h(url)}">{_h(label or url)}</a>'
            if _safe_url(url)
            else f"<span>{_h(label or 'Unlinked artifact')}</span>"
        )
        rows.append(
            f"<li>{content}{f'<code>{_h(digest)}</code>' if digest else ''}</li>"
        )
    return (
        f'<ul class="artifacts">{"".join(rows)}</ul>'
        if rows
        else _empty("No pair-linked artifacts supplied.")
    )


def _breadcrumbs(*items: tuple[str, str | None]) -> str:
    parts = ['<a href="../index.html">Atlas</a>']
    for label, href in items:
        parts.append(f'<a href="{_h(href)}">{_h(label)}</a>' if href else _h(label))
    return (
        '<nav class="breadcrumbs" aria-label="Breadcrumb">'
        + " / ".join(parts)
        + "</nav>"
    )


def _badge(status: str) -> str:
    normalized = re.sub(r"[^a-z0-9]+", "-", status.lower()).strip("-")
    return f'<span class="badge status-{_h(normalized)}">{_h(status)}</span>'


def _entity_name(record: Mapping[str, Any], fallback: str) -> str:
    return (
        _text(record.get("display_name") or record.get("name") or record.get("title"))
        or fallback
    )


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _iter_objects(value: Any) -> list[Mapping[str, Any]]:
    if isinstance(value, Mapping):
        return [
            dict(item, label=item.get("label", key))
            if isinstance(item, Mapping)
            else {"label": key, "url": item}
            for key, item in value.items()
        ]
    if isinstance(value, list):
        return [item if isinstance(item, Mapping) else {"url": item} for item in value]
    return []


def _first(mapping: Mapping[str, Any], *keys: str) -> Any:
    return next(
        (mapping[key] for key in keys if mapping.get(key) not in (None, "")), None
    )


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
    label = value.replace("_", " ").strip().title()
    return f"{label} (Secondary)" if value == "atlas_score" else label


def _safe_url(value: str) -> bool:
    if not value or any(char in value for char in ('"', "'", "<", ">", "\n", "\r")):
        return False
    parsed = urlsplit(value)
    return parsed.scheme in {"", "http", "https"} and not value.lower().startswith(
        ("javascript:", "data:")
    )


def _h(value: Any) -> str:
    return html.escape(_text(value), quote=True)


def _empty(message: str) -> str:
    return f'<p class="empty">{_h(message)}</p>'


_CSS = r"""
:root{--ink:#17231b;--muted:#5c695f;--paper:#f6f7f2;--card:#fff;--line:#dce2d9;--green:#285b3b;--gold:#d39a2e;--red:#9f3939;--blue:#315d7c}*{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);font:16px/1.5 system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}a{color:var(--green)}.site-header{height:64px;padding:0 max(24px,calc((100vw - 1180px)/2));display:flex;align-items:center;justify-content:space-between;background:#183626;color:#fff}.brand{color:#fff;text-decoration:none;font-weight:800;font-size:1.15rem}.release{font-size:.85rem;opacity:.8}main{max-width:1180px;margin:auto;padding:36px 24px 80px}footer{padding:24px;text-align:center;background:#e8ece5;color:var(--muted);font-size:.85rem}.hero{padding:48px;border-radius:22px;background:linear-gradient(135deg,#e6efe5,#fff);border:1px solid var(--line);margin-bottom:24px}.hero.compact{padding:30px}.hero h1{font-size:clamp(2rem,5vw,4.2rem);line-height:1.03;margin:.2em 0}.hero.compact h1{font-size:clamp(1.8rem,4vw,3rem)}.hero p{max-width:780px;color:#405047}.eyebrow{text-transform:uppercase;letter-spacing:.12em;font-weight:800;font-size:.75rem;color:var(--green)!important}.metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(130px,1fr));gap:12px;margin:24px 0 50px}.metrics div,.qualification div{background:var(--card);border:1px solid var(--line);border-radius:14px;padding:18px}.metrics strong{display:block;font-size:2rem}.metrics span,.qualification span{color:var(--muted);font-size:.8rem}.section-head{display:flex;gap:18px;align-items:center;justify-content:space-between;margin-top:46px}.section-head input{min-width:280px;padding:11px 14px;border:1px solid #abb8ac;border-radius:10px;background:#fff}.card-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(230px,1fr));gap:14px}.entity-card{display:block;padding:20px;background:var(--card);border:1px solid var(--line);border-radius:14px;text-decoration:none;color:inherit;transition:.15s}.entity-card:hover{transform:translateY(-2px);border-color:#8da393;box-shadow:0 8px 20px #18362614}.entity-card h3{margin:0}.entity-card p{color:var(--muted);margin:.4rem 0 0}.breadcrumbs{color:var(--muted);font-size:.85rem;margin-bottom:18px}.split{display:grid;grid-template-columns:minmax(280px,1fr) minmax(320px,1.4fr);gap:28px}.three-col{display:grid;grid-template-columns:repeat(3,1fr);gap:20px}.record{display:grid;grid-template-columns:minmax(130px,.7fr) 1.5fr;margin:0}.record dt,.record dd{padding:9px 0;border-bottom:1px solid var(--line);overflow-wrap:anywhere}.record dt{color:var(--muted);font-size:.82rem;font-weight:700}.record dd{margin:0}.record .record{grid-column:1/-1}.gallery{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:12px;margin:20px 0}.gallery figure{margin:0;background:#fff;border:1px solid var(--line);border-radius:14px;overflow:hidden}.gallery img{display:block;width:100%;aspect-ratio:4/3;object-fit:cover;background:#e5e8e3}.gallery figcaption{padding:8px 12px;color:var(--muted);font-size:.8rem}.table-wrap{overflow:auto;background:#fff;border:1px solid var(--line);border-radius:14px}table{width:100%;border-collapse:collapse}th,td{text-align:left;padding:13px;border-bottom:1px solid var(--line)}th{font-size:.75rem;text-transform:uppercase;letter-spacing:.05em;color:var(--muted)}td small{display:block;color:var(--muted)}.badge{display:inline-block;border-radius:999px;background:#e8ebe7;padding:3px 9px;font-size:.75rem;font-weight:750}.status-valid,.status-success,.status-qualified,.status-pass{background:#d8eddf;color:#215234}.status-invalid,.status-warning{background:#fff0c9;color:#72500d}.status-failed,.status-missing,.status-fail{background:#f5dada;color:#792727}.button{display:inline-block;padding:8px 13px;border-radius:8px;background:var(--green);color:#fff;text-decoration:none}.button.small{padding:5px 10px;font-size:.8rem}.qualification{display:flex;gap:12px;align-items:stretch;margin:20px 0}.qualification h2{margin-right:auto}.qualification div{min-width:150px}.qualification span{display:block;margin-top:6px}.callout{padding:14px 18px;border-radius:12px}.callout.ok{background:#dff0e4}.callout.danger{background:#f5e1d9}.artifacts{padding:0;list-style:none}.artifacts li{display:flex;justify-content:space-between;gap:20px;padding:10px;border-bottom:1px solid var(--line)}code{font-size:.82em;overflow-wrap:anywhere}.empty{color:var(--muted);font-style:italic;padding:18px;background:#fff;border:1px dashed #bdc7bd;border-radius:10px}@media(max-width:760px){main{padding:20px 14px 60px}.hero{padding:28px 20px}.split,.three-col{grid-template-columns:1fr}.section-head{align-items:stretch;flex-direction:column}.section-head input{min-width:0;width:100%}.qualification{display:grid}.record{grid-template-columns:1fr}.record dd{padding-top:0}}.show-all{font-size:.82rem;color:var(--muted);display:flex;align-items:center;gap:6px;white-space:nowrap}.pager{display:flex;align-items:center;justify-content:center;gap:12px;margin:14px 0 28px}.pager button{border:1px solid #aebcaf;border-radius:8px;background:#fff;color:var(--green);padding:7px 12px;cursor:pointer}.pager button:disabled{cursor:default;opacity:.4}.pager span{min-width:150px;text-align:center;color:var(--muted);font-size:.82rem}.analysis-app{background:#fff;border:1px solid var(--line);border-radius:16px;padding:20px}.analysis-controls{display:grid;grid-template-columns:repeat(auto-fit,minmax(155px,1fr));gap:12px;align-items:end}.analysis-controls label{display:grid;gap:5px;color:var(--muted);font-size:.8rem;font-weight:700}.analysis-controls select,.analysis-controls button,.analysis-tabs button{min-height:40px;border:1px solid #aebcaf;border-radius:9px;background:#fff;padding:7px 10px}.analysis-controls button,.analysis-tabs button{cursor:pointer;color:var(--green);font-weight:750}.analysis-tabs{display:flex;gap:8px;margin:22px 0}.analysis-tabs .is-active{background:var(--green);color:#fff}.analysis-result{overflow:auto}.analysis-bars{display:grid;gap:8px}.analysis-bar{display:grid;grid-template-columns:minmax(120px,220px) 1fr 70px;gap:10px;align-items:center}.analysis-bar-track{height:18px;background:#edf0eb;border-radius:5px;overflow:hidden}.analysis-bar-fill{height:100%;background:var(--green)}.matrix td{min-width:58px;text-align:center;font-variant-numeric:tabular-nums}.matrix th{position:sticky;background:#fff}.matrix-cell{background:color-mix(in srgb,var(--blue) calc(var(--heat)*85%),white);color:var(--ink)}.muted{color:var(--muted)}.site-header nav a{color:#fff;font-size:.9rem}
"""


_JS = r"""
document.querySelectorAll('[data-filter]').forEach(function(input){input.addEventListener('input',function(){var q=input.value.trim().toLowerCase();var root=document.getElementById(input.dataset.filter);if(!root)return;root.querySelectorAll('[data-search]').forEach(function(item){item.hidden=q && !item.dataset.search.includes(q);});});});
function updateAtlasPairTable(id){var root=document.getElementById(id);if(!root)return;var input=document.querySelector('[data-table-filter="'+id+'"]'),toggle=document.querySelector('[data-show-all-pairs="'+id+'"]'),q=input?input.value.trim().toLowerCase():'';root.querySelectorAll('tbody tr').forEach(function(row){var excluded=row.dataset.rankEligible==='false'&&!(toggle&&toggle.checked),filtered=q&&!row.dataset.search.includes(q);row.hidden=excluded||filtered;});}document.querySelectorAll('[data-table-filter]').forEach(function(input){input.addEventListener('input',function(){updateAtlasPairTable(input.dataset.tableFilter);});});document.querySelectorAll('[data-show-all-pairs]').forEach(function(toggle){toggle.addEventListener('change',function(){updateAtlasPairTable(toggle.dataset.showAllPairs);});});

(function(){
  var size=100;
  function parse(script){try{return JSON.parse(script.textContent||'[]');}catch(e){return [];}}
  function pager(root,kind){return document.querySelector('[data-'+kind+'-pager="'+root.id+'"]');}
  document.querySelectorAll('script[data-card-source]').forEach(function(script){
    var root=document.getElementById(script.dataset.cardSource),rows=parse(script),page=0,input=document.querySelector('[data-filter="'+script.dataset.cardSource+'"]'),nav=pager(root,'card');
    root._atlasCards=rows;
    function filtered(){var q=input?input.value.trim().toLowerCase():'';return rows.filter(function(row){return !q||String(row.search||'').toLowerCase().indexOf(q)!==-1;});}
    function render(reset){if(reset)page=0;var data=filtered(),pages=Math.max(1,Math.ceil(data.length/size));page=Math.min(page,pages-1);root.innerHTML=data.slice(page*size,(page+1)*size).map(function(row){return row.html;}).join('')||'<p class="empty">No matching records.</p>';if(nav){nav.querySelector('[data-page-status]').textContent=data.length+' records · page '+(page+1)+' of '+pages;nav.querySelector('[data-page-prev]').disabled=page===0;nav.querySelector('[data-page-next]').disabled=page>=pages-1;}}
    if(input)input.addEventListener('input',function(){render(true);});if(nav){nav.querySelector('[data-page-prev]').addEventListener('click',function(){if(page>0){page--;render(false);}});nav.querySelector('[data-page-next]').addEventListener('click',function(){if((page+1)*size<filtered().length){page++;render(false);}});}render(false);
  });
  document.querySelectorAll('script[data-pair-source]').forEach(function(script){
    var table=document.getElementById(script.dataset.pairSource),body=table&&table.querySelector('tbody'),rows=parse(script),page=0,input=document.querySelector('[data-table-filter="'+script.dataset.pairSource+'"]'),toggle=document.querySelector('[data-show-all-pairs="'+script.dataset.pairSource+'"]'),nav=pager(table,'pair');
    if(!table||!body)return;table._atlasRows=rows;
    function filtered(){var q=input?input.value.trim().toLowerCase():'',all=toggle&&toggle.checked;return rows.filter(function(row){return (all||row.eligible)&&(!q||String(row.search||'').toLowerCase().indexOf(q)!==-1);});}
    function render(reset){if(reset)page=0;var data=filtered(),pages=Math.max(1,Math.ceil(data.length/size));page=Math.min(page,pages-1);body.innerHTML=data.slice(page*size,(page+1)*size).map(function(row){return row.html;}).join('')||'<tr><td colspan="5" class="muted">No matching rows. Select “Show invalid, failed, and unranked cells” to include excluded cells.</td></tr>';if(nav){nav.querySelector('[data-page-status]').textContent=data.length+' rows · page '+(page+1)+' of '+pages;nav.querySelector('[data-page-prev]').disabled=page===0;nav.querySelector('[data-page-next]').disabled=page>=pages-1;}}
    if(input)input.addEventListener('input',function(){render(true);});if(toggle)toggle.addEventListener('change',function(){render(true);});if(nav){nav.querySelector('[data-page-prev]').addEventListener('click',function(){if(page>0){page--;render(false);}});nav.querySelector('[data-page-next]').addEventListener('click',function(){if((page+1)*size<filtered().length){page++;render(false);}});}render(false);
  });
})();

(function(){
  var app=document.getElementById('analysis-app'); if(!app)return;
  var message=document.getElementById('analysis-message'), controls=document.getElementById('analysis-controls'), result=document.getElementById('analysis-result'), summary=document.getElementById('analysis-summary');
  var dim=document.getElementById('analysis-filter-dimension'), val=document.getElementById('analysis-filter-value'), group=document.getElementById('analysis-group'), measure=document.getElementById('analysis-measure'), agg=document.getElementById('analysis-aggregation');
  var state={data:null,view:'table',current:[]};
  function esc(x){return String(x==null?'':x).replace(/[&<>\"]/g,function(c){return {'&':'&amp;','<':'&lt;','>':'&gt;','\"':'&quot;'}[c];});}
  function label(x){var text=String(x).replace(/_/g,' ').replace(/\b\w/g,function(c){return c.toUpperCase();});return x==='atlas_score'?text+' (secondary)':text;}
  function option(x){return '<option value="'+esc(x)+'">'+esc(label(x))+'</option>';}
  function csv(rows){if(!rows.length)return '';var keys=Object.keys(rows[0]);function cell(x){var t=String(x==null?'':x);return /[\",\n]/.test(t)?'\"'+t.replace(/\"/g,'\"\"')+'\"':t;}return keys.map(cell).join(',')+'\n'+rows.map(function(r){return keys.map(function(k){return cell(r[k]);}).join(',');}).join('\n')+'\n';}
  function filtered(){var key=dim.value,want=val.value;return state.data.rows.filter(function(r){return !want||String(r[key])===want;});}
  function refreshValues(){var values={};state.data.rows.forEach(function(r){if(r[dim.value]!=null&&r[dim.value]!=='')values[String(r[dim.value])]=1;});if(dim.value==='rank_eligible')values['1']=1;val.innerHTML='<option value="">All values</option>'+Object.keys(values).sort().map(function(x){var shown=dim.value==='rank_eligible'?((x==='1'||x==='true')?'Eligible':'Excluded'):x;return '<option value="'+esc(x)+'">'+esc(shown)+'</option>';}).join('');render();}
  function aggregate(rows){var buckets={};rows.forEach(function(r){var k=String(r[group.value]==null?'Not supplied':r[group.value]);(buckets[k]||(buckets[k]=[])).push(r);});var mode=agg.value,m=measure.value;return Object.keys(buckets).map(function(k){var values=buckets[k].map(function(r){return Number(r[m]);}).filter(Number.isFinite),out;if(mode==='count')out=buckets[k].length;else if(!values.length)out=null;else if(mode==='mean')out=values.reduce(function(a,b){return a+b;},0)/values.length;else if(mode==='min')out=Math.min.apply(Math,values);else out=Math.max.apply(Math,values);return {group:k,value:out,rows:buckets[k].length};}).sort(function(a,b){return (b.value==null?-Infinity:b.value)-(a.value==null?-Infinity:a.value);});}
  function table(data){return '<div class="table-wrap"><table><thead><tr><th>'+esc(label(group.value))+'</th><th>'+esc(label(agg.value))+'</th><th>Rows</th></tr></thead><tbody>'+data.map(function(r){return '<tr><td>'+esc(r.group)+'</td><td>'+esc(r.value==null?'—':Number(r.value).toPrecision(6))+'</td><td>'+r.rows+'</td></tr>';}).join('')+'</tbody></table></div>';}
  function bars(data){var finite=data.map(function(r){return r.value;}).filter(Number.isFinite),max=Math.max.apply(Math,[0].concat(finite.map(Math.abs)));return '<div class="analysis-bars">'+data.slice(0,40).map(function(r){var w=max&&Number.isFinite(r.value)?Math.abs(r.value)/max*100:0;return '<div class="analysis-bar"><span>'+esc(r.group)+'</span><div class="analysis-bar-track"><div class="analysis-bar-fill" style="width:'+w+'%"></div></div><strong>'+esc(r.value==null?'—':Number(r.value).toPrecision(4))+'</strong></div>';}).join('')+'</div>';}
  function matrix(rows){var xs=[],ys=[],cells={};rows.forEach(function(r){var x=String(r.drug_name||r.drug_id||'Unknown drug'),y=String(r.protein_name||r.target_id||'Unknown protein'),k=y+'\u0000'+x;(cells[k]||(cells[k]=[])).push(r);if(xs.indexOf(x)<0)xs.push(x);if(ys.indexOf(y)<0)ys.push(y);});xs=xs.slice(0,25);ys=ys.slice(0,25);var values=[];function value(y,x){var rs=cells[y+'\u0000'+x]||[];if(agg.value==='count')return rs.length||null;var nums=rs.map(function(r){return Number(r[measure.value]);}).filter(Number.isFinite);if(!nums.length)return null;if(agg.value==='min')return Math.min.apply(Math,nums);if(agg.value==='max')return Math.max.apply(Math,nums);return nums.reduce(function(a,b){return a+b;},0)/nums.length;}ys.forEach(function(y){xs.forEach(function(x){var v=value(y,x);if(v!=null)values.push(Math.abs(v));});});var max=Math.max.apply(Math,[0].concat(values));return '<div class="table-wrap"><table class="matrix"><thead><tr><th>Protein \\ Drug</th>'+xs.map(function(x){return '<th>'+esc(x)+'</th>';}).join('')+'</tr></thead><tbody>'+ys.map(function(y){return '<tr><th>'+esc(y)+'</th>'+xs.map(function(x){var v=value(y,x),heat=max&&v!=null?Math.abs(v)/max:0;return '<td class="matrix-cell" style="--heat:'+heat.toFixed(3)+'" title="'+esc(y+' × '+x)+'">'+esc(v==null?'—':Number(v).toPrecision(3))+'</td>';}).join('')+'</tr>';}).join('')+'</tbody></table></div>'+(xs.length===25||ys.length===25?'<p class="muted">Matrix preview limited to 25 × 25 groups; CSV export includes the current grouped result.</p>':'');}
  function render(){if(!state.data)return;var rows=filtered(),data=aggregate(rows);state.current=data.map(function(r){var out={};out[group.value]=r.group;out[agg.value+'_'+(agg.value==='count'?'rows':measure.value)]=r.value;out.source_rows=r.rows;return out;});summary.textContent=rows.length+' of '+state.data.rows.length+' pair rows · '+data.length+' groups';result.innerHTML=state.view==='matrix'?matrix(rows):(state.view==='histogram'?bars(data):table(data));}
  function setup(data){state.data=data;window.AtlasAnalysisDuckDBConfig=data.optional_engine||{enabled:false};if(!Array.isArray(data.rows)||!Array.isArray(data.dimensions))throw new Error('analysis dataset has an unsupported shape');if(!data.rows.length){message.textContent='No pair rows are available for browser analysis.';return;}if(!data.dimensions.length){message.textContent='No allowed analysis dimensions are populated in this release.';return;}dim.innerHTML=data.dimensions.map(option).join('');group.innerHTML=data.dimensions.map(option).join('');measure.innerHTML=(data.measures||[]).map(option).join('');if(!data.measures.length){measure.innerHTML='<option value="">No numeric measure</option>';agg.value='count';agg.querySelectorAll('option:not([value="count"])').forEach(function(o){o.disabled=true;});}message.hidden=true;controls.hidden=false;if(data.dimensions.indexOf('rank_eligible')>=0){dim.value='rank_eligible';refreshValues();var eligible=Array.from(val.options).find(function(o){return o.value==='1'||o.value==='true';});if(eligible){val.value=eligible.value;render();}}else{refreshValues();}}
  [val,group,measure,agg].forEach(function(el){el.addEventListener('change',render);});dim.addEventListener('change',refreshValues);document.getElementById('analysis-reset').addEventListener('click',function(){if(state.data.dimensions.indexOf('rank_eligible')>=0){dim.value='rank_eligible';refreshValues();val.value='1';}else{val.value='';}agg.value='count';render();});document.getElementById('analysis-download').addEventListener('click',function(){var blob=new Blob([csv(state.current)],{type:'text/csv;charset=utf-8'}),a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download='atlas-analysis-current.csv';a.click();setTimeout(function(){URL.revokeObjectURL(a.href);},1000);});document.querySelectorAll('[data-analysis-view]').forEach(function(b){b.addEventListener('click',function(){state.view=b.dataset.analysisView;document.querySelectorAll('[data-analysis-view]').forEach(function(x){x.classList.toggle('is-active',x===b);});render();});});
  fetch(app.dataset.source).then(function(r){if(!r.ok)throw new Error('HTTP '+r.status);return r.json();}).then(setup).catch(function(e){message.className='callout danger';message.innerHTML='Browser analysis data could not be loaded. Serve this directory with a static web server, or use the entity pages and downloadable release files. <small>'+esc(e.message)+'</small>';});
})();

"""
