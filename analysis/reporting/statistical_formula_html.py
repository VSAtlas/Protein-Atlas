"""HTML snippets documenting Atlas statistical calculations."""

from __future__ import annotations

import html

from analysis.reporting.statistical_formula_registry import collect_statistical_formulas


def _formula_card(row: dict[str, str]) -> str:
    title = html.escape(row["title"])
    formula = html.escape(row["formula"])
    notes = html.escape(row["notes"])
    source = html.escape(row["source"])
    return f"""
          <div>
            <h4>{title}</h4>
            <p><code>{formula}</code></p>
            <p>{notes}</p>
            <p class="formula-source">Source: <code>{source}</code></p>
          </div>"""


def render_statistical_formula_block() -> str:
    cards = "\n".join(_formula_card(row) for row in collect_statistical_formulas())
    return f"""
      <details class="formula-panel" open>
        <summary>Statistical formulas used in this report</summary>
        <p class="formula-note">Formula text is collected from metadata exported by the statistical implementation modules. Update the module-level <code>STATISTICAL_FORMULAS</code> next to any changed calculation so future reports stay synchronized with the code.</p>
        <div class="formula-grid">
{cards}
        </div>
      </details>
    """
