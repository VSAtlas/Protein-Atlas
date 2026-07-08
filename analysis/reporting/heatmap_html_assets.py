from __future__ import annotations

import base64
import html
import logging
import os
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from analysis.reporting.heatmap_html_runtime import heatmap_call
from config.output_paths import run_output_dir

_LOG = logging.getLogger("heatmap-html")
_CLUSTERGRAMMER_WARNED = False
_FONT_MIME_BY_EXT = {
    ".woff2": "font/woff2",
    ".woff": "font/woff",
    ".ttf": "font/ttf",
    ".eot": "application/vnd.ms-fontobject",
    ".svg": "image/svg+xml",
}

def _warn_clustergrammer_failure(exc: Exception) -> None:
    global _CLUSTERGRAMMER_WARNED
    if _CLUSTERGRAMMER_WARNED:
        return
    _CLUSTERGRAMMER_WARNED = True
    _LOG.warning("heatmap action=clustergrammer_fallback error=%s", exc)


def _ensure_pandas_ix_compat(pd_module: Any) -> None:
    if hasattr(pd_module.DataFrame, "ix"):
        return
    pd_module.DataFrame.ix = property(lambda self: self.loc)  # type: ignore[attr-defined]


def _load_clustergrammer() -> Optional[Tuple[Any, Any]]:
    try:
        import pandas as pd  # type: ignore[import-untyped]
        from clustergrammer import Network  # type: ignore[import-untyped]
    except Exception:
        return None
    _ensure_pandas_ix_compat(pd)
    return Network, pd


def _clustergrammer_key(name: str) -> str:
    return str(name).replace("_", " ")


def _escape_inline_text(text: str) -> str:
    return re.sub(r"</(script|style)", r"<\\/\1", text, flags=re.IGNORECASE)


def _read_asset_text(path: Path) -> str:
    with path.open("r", encoding="utf-8") as handle:
        return handle.read()


def _resolve_clustergrammer_assets_root(repo_root: Path) -> Path:
    direct = repo_root / "report_assets" / "clustergrammer"
    if direct.is_dir():
        return direct

    pointer = repo_root / "report_assets"
    if pointer.is_file():
        try:
            target_root = Path(pointer.read_text(encoding="utf-8").strip())
        except Exception as exc:
            raise FileNotFoundError(
                f"Failed to read report_assets pointer file: {pointer}"
            ) from exc
        candidate = target_root / "clustergrammer"
        if candidate.is_dir():
            return candidate

    raise FileNotFoundError(
        f"Clustergrammer assets root not found under {repo_root / 'report_assets'}"
    )


def _inline_font_urls(css_text: str, fonts_dir: Path) -> str:
    inlined_by_path: Dict[str, str] = {}

    def replacer(match: re.Match[str]) -> str:
        raw_url = match.group(1).strip().strip("'\"")
        if not raw_url or raw_url.startswith("data:") or raw_url.startswith("http"):
            return match.group(0)
        cleaned = raw_url.split("#", 1)[0].split("?", 1)[0]
        filename = Path(cleaned).name
        if not filename:
            return match.group(0)
        font_path = fonts_dir / filename
        if not font_path.exists():
            return match.group(0)
        cache_key = str(font_path.resolve())
        cached = inlined_by_path.get(cache_key)
        if cached is not None:
            return cached
        mime = _FONT_MIME_BY_EXT.get(font_path.suffix.lower())
        if not mime:
            return match.group(0)
        encoded = base64.b64encode(font_path.read_bytes()).decode("ascii")
        replacement = f"url('data:{mime};base64,{encoded}')"
        inlined_by_path[cache_key] = replacement
        return replacement

    return re.sub(r"url\(([^)]+)\)", replacer, css_text)


def _clustergrammer_inline_assets(repo_root: Path) -> str:
    assets_root = _resolve_clustergrammer_assets_root(repo_root)
    fonts_dir = assets_root / "lib" / "fonts"
    css_paths = [
        assets_root / "lib" / "css" / "bootstrap.css",
        assets_root / "lib" / "css" / "font-awesome.min.css",
        assets_root / "css" / "custom.css",
    ]
    js_paths = [
        assets_root / "lib" / "js" / "d3.js",
        assets_root / "lib" / "js" / "jquery-1.11.2.min.js",
        assets_root / "lib" / "js" / "underscore-min.js",
        assets_root / "lib" / "js" / "bootstrap.min.js",
        assets_root / "clustergrammer.js",
    ]
    css_parts = []
    for path in css_paths:
        text = _read_asset_text(path)
        text = _inline_font_urls(text, fonts_dir)
        css_parts.append(text)
    css_blob = _escape_inline_text("\n".join(css_parts))
    blocks = [f"<style>\n{css_blob}\n</style>"]
    for path in js_paths:
        js_text = _escape_inline_text(_read_asset_text(path))
        blocks.append(f"<script>\n{js_text}\n</script>")
    return "\n".join(blocks)


def _clustergrammer_asset_urls(
    repo_root: Path, run_id: str, offline_assets: bool
) -> Dict[str, str]:
    if offline_assets:
        report_dir = run_output_dir(repo_root, "data", run_id)
        assets_root = _resolve_clustergrammer_assets_root(repo_root)
        rel_root = Path(os.path.relpath(assets_root, report_dir)).as_posix()
        return {
            "d3": f"{rel_root}/lib/js/d3.js",
            "jquery": f"{rel_root}/lib/js/jquery-1.11.2.min.js",
            "underscore": f"{rel_root}/lib/js/underscore-min.js",
            "bootstrap_js": f"{rel_root}/lib/js/bootstrap.min.js",
            "bootstrap_css": f"{rel_root}/lib/css/bootstrap.css",
            "font_awesome_css": f"{rel_root}/lib/css/font-awesome.min.css",
            "custom_css": f"{rel_root}/css/custom.css",
            "clustergrammer_js": f"{rel_root}/clustergrammer.js",
        }
    return {
        "d3": "https://d3js.org/d3.v3.min.js",
        "jquery": "https://code.jquery.com/jquery-1.11.2.min.js",
        "underscore": "https://cdnjs.cloudflare.com/ajax/libs/underscore.js/1.8.3/underscore-min.js",
        "bootstrap_js": "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/js/bootstrap.min.js",
        "bootstrap_css": "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css",
        "font_awesome_css": "https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css",
        "custom_css": "https://unpkg.com/clustergrammer@1.19.5/css/custom.css",
        "clustergrammer_js": "https://unpkg.com/clustergrammer@1.19.5/clustergrammer.js",
    }


def _clustergrammer_asset_fallback_urls() -> Dict[str, List[str]]:
    return {
        "d3": [
            "https://d3js.org/d3.v3.min.js",
            "https://cdn.jsdelivr.net/npm/d3@3.5.17/d3.min.js",
            "https://cdnjs.cloudflare.com/ajax/libs/d3/3.5.17/d3.min.js",
        ],
        "jquery": [
            "https://code.jquery.com/jquery-1.11.2.min.js",
            "https://cdn.jsdelivr.net/npm/jquery@1.11.2/dist/jquery.min.js",
            "https://ajax.googleapis.com/ajax/libs/jquery/1.11.2/jquery.min.js",
        ],
        "underscore": [
            "https://cdnjs.cloudflare.com/ajax/libs/underscore.js/1.8.3/underscore-min.js",
            "https://cdn.jsdelivr.net/npm/underscore@1.8.3/underscore-min.js",
            "https://unpkg.com/underscore@1.8.3/underscore-min.js",
        ],
        "bootstrap_js": [
            "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/js/bootstrap.min.js",
            "https://cdn.jsdelivr.net/npm/bootstrap@3.3.6/dist/js/bootstrap.min.js",
            "https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.3.6/js/bootstrap.min.js",
        ],
        "bootstrap_css": [
            "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css",
            "https://cdn.jsdelivr.net/npm/bootstrap@3.3.6/dist/css/bootstrap.min.css",
            "https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.3.6/css/bootstrap.min.css",
        ],
        "font_awesome_css": [
            "https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css",
            "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css",
            "https://cdn.jsdelivr.net/npm/font-awesome@4.7.0/css/font-awesome.min.css",
        ],
        "custom_css": [
            "https://unpkg.com/clustergrammer@1.19.5/css/custom.css",
            "https://cdn.jsdelivr.net/npm/clustergrammer@1.19.5/css/custom.css",
        ],
        "clustergrammer_js": [
            "https://unpkg.com/clustergrammer@1.19.5/clustergrammer.js",
            "https://cdn.jsdelivr.net/npm/clustergrammer@1.19.5/clustergrammer.js",
        ],
    }


def _unique_preserve_order(values: Iterable[str]) -> List[str]:
    seen: Set[str] = set()
    ordered: List[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ordered


def _normalize_asset_mode(raw_mode: Optional[str]) -> str:
    mode = str(raw_mode or "").strip().lower()
    if mode in {"inline", "relative", "cdn", "auto"}:
        return mode
    return "auto"


def _clustergrammer_asset_plan(
    repo_root: Path,
    run_id: str,
    offline_assets: bool,
    inline_assets: bool,
    asset_mode: str = "auto",
) -> Dict[str, Any]:
    normalized_mode = _normalize_asset_mode(asset_mode)
    notices: List[str] = []
    if normalized_mode == "inline":
        try:
            return {
                "mode": "inline",
                "html": heatmap_call("_clustergrammer_inline_assets", repo_root),
                "asset_urls": {},
                "notices": notices,
            }
        except Exception as exc:
            _LOG.warning("heatmap action=inline_assets_failed error=%s", exc)
            notices.append(f"inline_assets_failed: {exc}")
            normalized_mode = "cdn"
    if normalized_mode == "relative":
        assets = heatmap_call("_clustergrammer_asset_urls", repo_root, run_id, True)
        return {
            "mode": "relative",
            "asset_urls": {name: [url] for name, url in assets.items()},
            "notices": notices,
            "html": "\n".join(
                [
                    f"<link rel=\"stylesheet\" href=\"{html.escape(assets['bootstrap_css'], quote=True)}\">",
                    f"<link rel=\"stylesheet\" href=\"{html.escape(assets['font_awesome_css'], quote=True)}\">",
                    f"<link rel=\"stylesheet\" href=\"{html.escape(assets['custom_css'], quote=True)}\">",
                    f"<script src=\"{html.escape(assets['d3'], quote=True)}\"></script>",
                    f"<script src=\"{html.escape(assets['jquery'], quote=True)}\"></script>",
                    f"<script src=\"{html.escape(assets['underscore'], quote=True)}\"></script>",
                    f"<script src=\"{html.escape(assets['bootstrap_js'], quote=True)}\"></script>",
                    f"<script src=\"{html.escape(assets['clustergrammer_js'], quote=True)}\"></script>",
                ]
            ),
        }
    if normalized_mode == "cdn":
        assets = heatmap_call("_clustergrammer_asset_urls", repo_root, run_id, False)
        fallback_urls = _clustergrammer_asset_fallback_urls()
        asset_urls = {
            name: _unique_preserve_order(
                [url for url in [assets[name], *fallback_urls.get(name, [])] if url]
            )
            for name in assets
        }
        return {
            "mode": "cdn",
            "asset_urls": asset_urls,
            "notices": notices,
            "html": "\n".join(
                [
                    f"<link rel=\"stylesheet\" href=\"{html.escape(assets['bootstrap_css'], quote=True)}\">",
                    f"<link rel=\"stylesheet\" href=\"{html.escape(assets['font_awesome_css'], quote=True)}\">",
                    f"<link rel=\"stylesheet\" href=\"{html.escape(assets['custom_css'], quote=True)}\">",
                    f"<script src=\"{html.escape(assets['d3'], quote=True)}\"></script>",
                    f"<script src=\"{html.escape(assets['jquery'], quote=True)}\"></script>",
                    f"<script src=\"{html.escape(assets['underscore'], quote=True)}\"></script>",
                    f"<script src=\"{html.escape(assets['bootstrap_js'], quote=True)}\"></script>",
                    f"<script src=\"{html.escape(assets['clustergrammer_js'], quote=True)}\"></script>",
                ]
            ),
        }

    use_offline_assets = offline_assets
    if inline_assets:
        try:
            return {
                "mode": "inline",
                "html": heatmap_call("_clustergrammer_inline_assets", repo_root),
                "asset_urls": {},
                "notices": notices,
            }
        except Exception as exc:
            _LOG.warning("heatmap action=inline_assets_failed error=%s", exc)
            notices.append(f"inline_assets_failed: {exc}")
            use_offline_assets = False
    try:
        assets = heatmap_call("_clustergrammer_asset_urls", repo_root, run_id, use_offline_assets)
    except Exception as exc:
        if use_offline_assets:
            _LOG.warning("heatmap action=offline_assets_failed error=%s", exc)
            notices.append(f"offline_assets_failed: {exc}")
            assets = heatmap_call("_clustergrammer_asset_urls", repo_root, run_id, False)
        else:
            raise
    mode = "relative" if use_offline_assets else "cdn"
    if mode == "cdn":
        fallback_urls = _clustergrammer_asset_fallback_urls()
        asset_urls = {
            name: _unique_preserve_order(
                [url for url in [assets[name], *fallback_urls.get(name, [])] if url]
            )
            for name in assets
        }
    else:
        asset_urls = {name: [url] for name, url in assets.items()}
    return {
        "mode": mode,
        "asset_urls": asset_urls,
        "notices": notices,
        "html": "\n".join(
        [
            f"<link rel=\"stylesheet\" href=\"{html.escape(assets['bootstrap_css'], quote=True)}\">",
            f"<link rel=\"stylesheet\" href=\"{html.escape(assets['font_awesome_css'], quote=True)}\">",
            f"<link rel=\"stylesheet\" href=\"{html.escape(assets['custom_css'], quote=True)}\">",
            f"<script src=\"{html.escape(assets['d3'], quote=True)}\"></script>",
            f"<script src=\"{html.escape(assets['jquery'], quote=True)}\"></script>",
            f"<script src=\"{html.escape(assets['underscore'], quote=True)}\"></script>",
            f"<script src=\"{html.escape(assets['bootstrap_js'], quote=True)}\"></script>",
            f"<script src=\"{html.escape(assets['clustergrammer_js'], quote=True)}\"></script>",
        ]
        ),
    }
