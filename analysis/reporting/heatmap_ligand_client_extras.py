from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class HeatmapLigandOrganizationExtras:
    style_block: str
    html_block: str
    script_block: str


def build_heatmap_ligand_organization_extras(
    safe_run_id: str,
) -> HeatmapLigandOrganizationExtras:
    ligand_org_none_id = f"cg-heatmap-ligand-org-none-{safe_run_id}"
    ligand_org_chemotype_id = f"cg-heatmap-ligand-org-chemotype-{safe_run_id}"
    ligand_org_effect_id = f"cg-heatmap-ligand-org-effect-{safe_run_id}"
    ligand_chemotype_panel_id = f"cg-heatmap-ligand-panel-chemotype-{safe_run_id}"
    ligand_effect_panel_id = f"cg-heatmap-ligand-panel-effect-{safe_run_id}"
    motif_search_input_id = f"cg-heatmap-motif-search-{safe_run_id}"
    motif_search_clear_id = f"cg-heatmap-motif-search-clear-{safe_run_id}"

    style_block = """
<style>
  .heatmap-motif-search { margin-top: 8px; display: grid; grid-template-columns: minmax(0, 1fr) auto; gap: 8px; align-items: end; }
  .heatmap-motif-search label { font-size: 0.78em; font-weight: 700; color: #314631; text-transform: uppercase; letter-spacing: 0.03em; display: grid; gap: 4px; }
  .heatmap-motif-search input { height: 34px; border: 1px solid #b8c5af; border-radius: 8px; padding: 0 10px; background: rgba(255,255,255,0.96); color: #243224; }
  .heatmap-motif-search button { height: 34px; border: 1px solid #8ea088; border-radius: 8px; background: #fff; color: #294331; font-weight: 700; cursor: pointer; padding: 0 12px; }
  .cg-heatmap-viz .cg-motif-search-hit { stroke: #295f93 !important; stroke-width: 0.35px !important; }
</style>
""".strip()

    html_block = f"""
  <div class="heatmap-analysis-controls">
    <div class="heatmap-control">
      <label>Ligand Organization</label>
      <div class="heatmap-segmented">
        <button type="button" id="{ligand_org_none_id}" class="is-active">None</button>
        <button type="button" id="{ligand_org_chemotype_id}">Chemotype</button>
        <button type="button" id="{ligand_org_effect_id}">Drug effect</button>
      </div>
    </div>
  </div>
  <div class="heatmap-org-panels">
    <div class="heatmap-org-panel" id="{ligand_chemotype_panel_id}"></div>
    <div class="heatmap-org-panel" id="{ligand_effect_panel_id}"></div>
  </div>
  <div class="heatmap-motif-search">
    <label for="{motif_search_input_id}">Search SAR motifs
      <input type="search" id="{motif_search_input_id}" placeholder="e.g. steroidal core, para-halogen, basic amine">
    </label>
    <button type="button" id="{motif_search_clear_id}">Clear motifs</button>
  </div>
""".strip()

    script_block = f"""
    var ligandOrgNoneEl = document.getElementById({ligand_org_none_id!r});
    var ligandOrgChemotypeEl = document.getElementById({ligand_org_chemotype_id!r});
    var ligandOrgEffectEl = document.getElementById({ligand_org_effect_id!r});
    var ligandChemotypePanelEl = document.getElementById({ligand_chemotype_panel_id!r});
    var ligandEffectPanelEl = document.getElementById({ligand_effect_panel_id!r});
    var motifSearchInputEl = document.getElementById({motif_search_input_id!r});
    var motifSearchClearEl = document.getElementById({motif_search_clear_id!r});
    var ligandOrganizationMode = "none";
    var ligandOrganizationCategory = "";
    var activeMotifSearch = "";
    function ligandOrganizationModeLabel(mode) {{
      if (mode === "chemotype") return "Chemotype";
      if (mode === "drug_effect") return "Drug effect";
      return "None";
    }}
    function ligandOrganizationBlocks(mode) {{
      var meta = analysisMeta && analysisMeta.ligand_organization ? analysisMeta.ligand_organization : {{}};
      if (mode === "chemotype") return Array.isArray(meta.chemotype_blocks) ? meta.chemotype_blocks : [];
      if (mode === "drug_effect") return Array.isArray(meta.drug_effect_blocks) ? meta.drug_effect_blocks : [];
      return [];
    }}
    function ligandOrganizationValue(rowName, mode) {{
      var meta = lookupRowLabelMeta(rowName);
      if (mode === "chemotype") return String(meta.chemotype_primary || "Other / Unassigned");
      if (mode === "drug_effect") return String(meta.drug_effect_primary || "Other / Unassigned");
      return "";
    }}
    function ligandOrganizationColorMap(mode) {{
      var palette = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#117733", "#332288"];
      var blocks = ligandOrganizationBlocks(mode);
      var out = Object.create(null);
      for (var i = 0; i < blocks.length; i++) {{
        var key = String(blocks[i] && blocks[i].key || "");
        out[key] = palette[i % palette.length];
      }}
      return out;
    }}
    function ligandOrganizationColorForValue(value) {{
      if (ligandOrganizationMode === "none") return "";
      var colors = ligandOrganizationColorMap(ligandOrganizationMode);
      return String(colors[String(value || "")] || "#7a6a44");
    }}
    function motifSearchText() {{
      return String(activeMotifSearch || "").trim().toLowerCase();
    }}
    function ligandSearchFields(meta) {{
      if (!meta) return [];
      return []
        .concat(meta.motif_tags || [])
        .concat(meta.chemotype_primary || [])
        .concat(meta.chemotype_series || [])
        .concat(meta.drug_effect_primary || [])
        .concat(meta.drug_effect_memberships || []);
    }}
    function ligandMatchesMotifSearch(meta) {{
      var query = motifSearchText();
      if (!query) return false;
      var haystack = ligandSearchFields(meta).map(function(value) {{ return String(value || "").toLowerCase(); }});
      for (var i = 0; i < haystack.length; i++) {{
        if (haystack[i] && haystack[i].indexOf(query) !== -1) return true;
      }}
      return false;
    }}
    function organizeRowsByMode(rowNames) {{
      if (ligandOrganizationMode === "none") return rowNames.slice();
      var blocks = ligandOrganizationBlocks(ligandOrganizationMode);
      var blockOrder = Object.create(null);
      for (var i = 0; i < blocks.length; i++) blockOrder[String(blocks[i] && blocks[i].key || "")] = i;
      var activeCategory = String(ligandOrganizationCategory || "");
      var tieOrder = stableTieOrder(rowNames);
      return rowNames.slice().sort(function(a, b) {{
        var aKey = ligandOrganizationValue(a, ligandOrganizationMode);
        var bKey = ligandOrganizationValue(b, ligandOrganizationMode);
        if (activeCategory) {{
          var aCat = aKey === activeCategory ? 0 : 1;
          var bCat = bKey === activeCategory ? 0 : 1;
          if (aCat !== bCat) return aCat - bCat;
        }}
        var aRank = Object.prototype.hasOwnProperty.call(blockOrder, aKey) ? blockOrder[aKey] : 99999;
        var bRank = Object.prototype.hasOwnProperty.call(blockOrder, bKey) ? blockOrder[bKey] : 99999;
        if (aRank !== bRank) return aRank - bRank;
        if (aKey !== bKey) return aKey < bKey ? -1 : 1;
        return Number(tieOrder[a] || 0) - Number(tieOrder[b] || 0);
      }});
    }}
    function sortRowsWithinOrganization(rowNames, scorer) {{
      var blocks = ligandOrganizationBlocks(ligandOrganizationMode);
      var blockOrder = Object.create(null);
      for (var i = 0; i < blocks.length; i++) blockOrder[String(blocks[i] && blocks[i].key || "")] = i;
      var activeCategory = String(ligandOrganizationCategory || "");
      var tieOrder = stableTieOrder(rowNames);
      return rowNames.slice().sort(function(a, b) {{
        var aKey = ligandOrganizationValue(a, ligandOrganizationMode);
        var bKey = ligandOrganizationValue(b, ligandOrganizationMode);
        if (activeCategory) {{
          var aCat = aKey === activeCategory ? 0 : 1;
          var bCat = bKey === activeCategory ? 0 : 1;
          if (aCat !== bCat) return aCat - bCat;
        }}
        var aRank = Object.prototype.hasOwnProperty.call(blockOrder, aKey) ? blockOrder[aKey] : 99999;
        var bRank = Object.prototype.hasOwnProperty.call(blockOrder, bKey) ? blockOrder[bKey] : 99999;
        if (aRank !== bRank) return aRank - bRank;
        var av = scorer(a);
        var bv = scorer(b);
        var ax = av == null ? -Infinity : av;
        var bx = bv == null ? -Infinity : bv;
        if (bx !== ax) return bx - ax;
        return Number(tieOrder[a] || 0) - Number(tieOrder[b] || 0);
      }});
    }}
    function renderLigandAnnotationPanel(panelEl, title, note, blocks, activeMode, modeName) {{
      if (!panelEl) return;
      panelEl.classList.toggle("is-active", activeMode === modeName);
      var chips = blocks.map(function(block) {{
        var key = String(block.key || block.label || "");
        var chipColor = ligandOrganizationColorMap(modeName)[key] || "#7a6a44";
        var isActive = activeMode === modeName && key && key === String(ligandOrganizationCategory || "");
        return '<button type="button" class="heatmap-org-chip' + (isActive ? ' is-active' : '') + '" data-org-kind="ligand" data-org-mode="' + escapeHtml(modeName) + '" data-org-category="' + escapeHtml(key) + '" style="border-color:' + chipColor + '55; box-shadow: inset 0 0 0 1px ' + chipColor + '1c;"><span class="heatmap-org-chip-name" style="color:' + chipColor + '">' + escapeHtml(String(block.label || "")) + '</span><span class="heatmap-org-chip-count">' + escapeHtml(String(block.count || 0)) + '</span></button>';
      }}).join("");
      panelEl.innerHTML = '<div class="heatmap-org-panel-head"><div class="heatmap-org-panel-title">' + escapeHtml(title) + '</div><div class="heatmap-org-panel-note">' + escapeHtml(note) + '</div></div><div class="heatmap-org-panel-chips">' + chips + '</div>';
    }}
    function renderLigandAnnotationPanels() {{
      renderLigandAnnotationPanel(ligandChemotypePanelEl, "Chemotype", "Series/scaffold-like row grouping", ligandOrganizationBlocks("chemotype"), ligandOrganizationMode, "chemotype");
      renderLigandAnnotationPanel(ligandEffectPanelEl, "Drug effect", "Pharmacology and therapeutic-class row grouping", ligandOrganizationBlocks("drug_effect"), ligandOrganizationMode, "drug_effect");
    }}
    function updateLigandOrganizationControls() {{
      if (ligandOrgNoneEl) ligandOrgNoneEl.classList.toggle("is-active", ligandOrganizationMode === "none");
      if (ligandOrgChemotypeEl) ligandOrgChemotypeEl.classList.toggle("is-active", ligandOrganizationMode === "chemotype");
      if (ligandOrgEffectEl) ligandOrgEffectEl.classList.toggle("is-active", ligandOrganizationMode === "drug_effect");
      renderLigandAnnotationPanels();
    }}
    function activateLigandOrganizationCategory(mode, category) {{
      var resolvedMode = mode === "chemotype" || mode === "drug_effect" ? mode : "none";
      var resolvedCategory = String(category || "").trim();
      if (resolvedMode === "none") {{
        ligandOrganizationMode = "none";
        ligandOrganizationCategory = "";
      }} else if (ligandOrganizationMode === resolvedMode && ligandOrganizationCategory === resolvedCategory) {{
        ligandOrganizationCategory = "";
      }} else {{
        ligandOrganizationMode = resolvedMode;
        ligandOrganizationCategory = resolvedCategory;
      }}
      if (typeof resetHeatmapSelection === "function") resetHeatmapSelection();
      updateOrganizationControls();
      renderCurrentMode();
    }}
    function appendLigandAnnotationLines(lines, meta) {{
      if (!meta) return lines;
      lines.push("chemotype: " + String(meta.chemotype_primary || "Other / Unassigned"));
      lines.push("drug_effect: " + String(meta.drug_effect_primary || "Other / Unassigned"));
      if (Array.isArray(meta.motif_tags) && meta.motif_tags.length) lines.push("motif_tags: " + String(meta.motif_tags.join("; ")));
      if (meta.ligand_annotation_confidence) lines.push("ligand_annotation_confidence: " + String(meta.ligand_annotation_confidence));
      if (Array.isArray(meta.ligand_annotation_sources) && meta.ligand_annotation_sources.length) lines.push("ligand_annotation_sources: " + String(meta.ligand_annotation_sources.join("; ")));
      return lines;
    }}
    function bindLigandOrganizationControls() {{
      if (ligandOrgNoneEl) ligandOrgNoneEl.addEventListener("click", function() {{ ligandOrganizationMode = "none"; ligandOrganizationCategory = ""; updateOrganizationControls(); renderCurrentMode(); }});
      if (ligandOrgChemotypeEl) ligandOrgChemotypeEl.addEventListener("click", function() {{ ligandOrganizationMode = "chemotype"; ligandOrganizationCategory = ""; updateOrganizationControls(); renderCurrentMode(); }});
      if (ligandOrgEffectEl) ligandOrgEffectEl.addEventListener("click", function() {{ ligandOrganizationMode = "drug_effect"; ligandOrganizationCategory = ""; updateOrganizationControls(); renderCurrentMode(); }});
      var panelClick = function(evt) {{
        var target = evt && evt.target ? evt.target : null;
        if (!target || typeof target.closest !== "function") return;
        var chip = target.closest('[data-org-kind="ligand"][data-org-mode][data-org-category]');
        if (!chip) return;
        evt.preventDefault();
        evt.stopPropagation();
        activateLigandOrganizationCategory(
          String(chip.getAttribute("data-org-mode") || ""),
          String(chip.getAttribute("data-org-category") || ""),
        );
      }};
      if (ligandChemotypePanelEl) ligandChemotypePanelEl.addEventListener("click", panelClick);
      if (ligandEffectPanelEl) ligandEffectPanelEl.addEventListener("click", panelClick);
      if (motifSearchInputEl) motifSearchInputEl.addEventListener("input", function() {{ activeMotifSearch = String(motifSearchInputEl.value || "").trim(); updateActiveLabelStates(); updateViewSummary(baseVizData()); if (typeof renderDiscoveryPanel === "function") renderDiscoveryPanel(); syncReportStateHash(); }});
      if (motifSearchClearEl) motifSearchClearEl.addEventListener("click", function() {{ activeMotifSearch = ""; if (motifSearchInputEl) motifSearchInputEl.value = ""; updateActiveLabelStates(); updateViewSummary(baseVizData()); if (typeof renderDiscoveryPanel === "function") renderDiscoveryPanel(); syncReportStateHash(); }});
    }}
""".strip()

    return HeatmapLigandOrganizationExtras(
        style_block=style_block,
        html_block=html_block,
        script_block=script_block,
    )
