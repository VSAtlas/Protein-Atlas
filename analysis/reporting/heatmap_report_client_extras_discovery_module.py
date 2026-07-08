from __future__ import annotations


def render_discovery_js_module() -> str:
    return """
    function discoveryEntries() {
      return allRowNames().map(function(name) {
        var meta = lookupRowLabelMeta(name);
        var live = currentAnalysis && currentAnalysis.rowMetrics ? (currentAnalysis.rowMetrics[name] || {}) : {};
        return {
          name: name,
          label: name,
          subtitle: String(meta.drug_effect_primary || meta.library || ""),
          group: String(meta.chemotype_primary || ""),
          hits: String(live.significant_hits_text || "0 / 0"),
          prom: asNumber(live.promiscuity_z),
          frac: asNumber(live.hit_fraction),
          meta: meta,
          searchable: [name, meta.chemotype_primary, meta.drug_effect_primary]
            .concat(meta.motif_tags || [], meta.side_effects || [], meta.safety_buckets || [])
            .join(" ").toLowerCase()
        };
      });
    }
    function annotationEntries() {
      var out = Object.create(null);
      function addEntry(kind, label, subtitle, actionMode, extraQuery) {
        var cleanLabel = String(label || "").trim();
        if (!cleanLabel) return;
        var key = kind + "||" + cleanLabel.toLowerCase();
        if (!out[key]) out[key] = { kind: kind, label: cleanLabel, subtitle: subtitle, actionMode: actionMode, query: String(extraQuery || cleanLabel), count: 0 };
        out[key].count += 1;
      }
      allColNames().forEach(function(name) {
        var meta = lookupColLabelMeta(name);
        addEntry("family", meta.protein_family, "Protein family", "family", meta.protein_family);
        addEntry("adme", meta.adme_category, "ADME category", "adme", meta.adme_category);
        addEntry("pathway", meta.primary_display_pathway, "Pathway group", "pathway", meta.primary_display_pathway);
        addEntry("safety_bucket", meta.primary_display_safety, "Safety bucket", "safety", meta.primary_display_safety);
        (meta.secondary_safety_buckets || []).forEach(function(value) {
          addEntry("safety_bucket", value, "Secondary safety bucket", "safety", value);
        });
        (meta.side_effects || []).forEach(function(value) {
          addEntry("side_effect", value, "Side-effect evidence", "safety", value);
        });
      });
      return Object.keys(out).map(function(key) {
        var item = out[key];
        item.searchable = (item.label + " " + item.subtitle).toLowerCase();
        return item;
      });
    }
    function organizationChipHtml(mode, item, label) {
      var blockLabel = String(item.label || "");
      var isActive = activeDiscoveryHighlightMode === mode && activeDiscoveryHighlightValue.toLowerCase() === blockLabel.toLowerCase();
      return '<button type="button" class="heatmap-discovery-chip' + (isActive ? ' is-active' : '') + '" data-org-mode="' + escapeHtml(mode) + '" data-org-value="' + escapeHtml(blockLabel) + '"><strong>' + escapeHtml(blockLabel) + '</strong><span>' + escapeHtml(String(item.count || 0)) + '</span></button>';
    }
    function setDiscoveryAnnotationFilter(mode, value) {
      activeDiscoveryMode = "rows";
      activeDiscoveryHighlightMode = String(mode || "");
      activeDiscoveryHighlightValue = String(value || "");
      activeDiscoveryFilter = "all";
      renderDiscoveryPanel();
      renderActiveStateBar();
      syncReportStateHash();
    }
    function renderOrganizationChips(containerEl, mode, blocks, limit) {
      if (!containerEl) return;
      var shown = (Array.isArray(blocks) ? blocks : []).slice(0, limit);
      containerEl.innerHTML = shown.map(function(item) { return organizationChipHtml(mode, item); }).join("");
      var buttons = containerEl.querySelectorAll("[data-org-mode]");
      for (var i = 0; i < buttons.length; i++) {
        buttons[i].addEventListener("click", function(evt) {
          evt.preventDefault();
          evt.stopPropagation();
          var nextMode = String(this.getAttribute("data-org-mode") || "none");
          var nextValue = String(this.getAttribute("data-org-value") || "");
          var isSame = activeDiscoveryHighlightMode === nextMode && activeDiscoveryHighlightValue.toLowerCase() === nextValue.toLowerCase();
          setDiscoveryAnnotationFilter(isSame ? "" : nextMode, isSame ? "" : nextValue);
        });
      }
    }
    function renderDiscoveryChips() {
      return;
    }
    function entryMatchesDiscoveryFilter(entry) {
      return true;
    }
    function entryMatchesDiscoveryAnnotation(entry) {
      return true;
    }
    function renderDiscoveryPanel() {
      if (!discoveryPanelEl || !discoveryTableEl) return;
      activeDiscoveryMode = "rows";
      if (discoveryRowsEl) discoveryRowsEl.classList.toggle("is-active", true);
      if (discoveryColsEl) discoveryColsEl.classList.toggle("is-active", false);
      if (discoveryOrgFamilyEl) discoveryOrgFamilyEl.classList.toggle("is-active", activeDiscoveryHighlightMode === "family");
      if (discoveryOrgPathwayEl) discoveryOrgPathwayEl.classList.toggle("is-active", activeDiscoveryHighlightMode === "pathway");
      if (discoveryOrgAdmeEl) discoveryOrgAdmeEl.classList.toggle("is-active", activeDiscoveryHighlightMode === "adme");
      if (discoveryOrgSafetyEl) discoveryOrgSafetyEl.classList.toggle("is-active", activeDiscoveryHighlightMode === "safety");
      var query = discoverySearchEl ? String(discoverySearchEl.value || "").trim().toLowerCase() : "";
      var entries = discoveryEntries().filter(function(entry) {
        return (!query || entry.searchable.indexOf(query) !== -1) && entryMatchesDiscoveryFilter(entry) && entryMatchesDiscoveryAnnotation(entry);
      }).sort(function(a, b) {
        var bProm = b.prom == null ? -Infinity : b.prom;
        var aProm = a.prom == null ? -Infinity : a.prom;
        if (bProm !== aProm) return bProm - aProm;
        return a.label < b.label ? -1 : (a.label > b.label ? 1 : 0);
      });
      var limited = entries.slice(0, 40);
      if (discoverySummaryEl) {
        var coverage = analysisMeta.coverage_summary || {};
        discoverySummaryEl.textContent = "Ligands: showing " + limited.length + " of " + entries.length + " matched items | interactive matrix contains "
          + Number(coverage.interactive_ligands || allRowNames().length)
          + " ligands"
          + (activeDiscoveryFilter !== "all" ? (" | filter: " + (activeDiscoveryFilter === "high_safety" ? "high-confidence safety" : "unresolved mappings")) : "")
          + (activeDiscoveryHighlightMode && activeDiscoveryHighlightValue ? (" | annotation: " + activeDiscoveryHighlightMode + " = " + activeDiscoveryHighlightValue) : "")
          + (query ? (" | search: " + query) : "");
      }
      var header = "<thead><tr><th>Ligand</th><th>Chemotype</th><th>Drug effect</th><th>Promiscuity z-score</th><th>Hit fraction</th><th>Significant hits</th><th></th></tr></thead>";
      var body = limited.map(function(entry) {
        return "<tr>"
          + "<td>" + escapeHtml(entry.label) + "</td>"
          + "<td>" + escapeHtml(entry.group || "") + "</td>"
          + "<td>" + escapeHtml(entry.subtitle || "") + "</td>"
          + '<td class="num">' + escapeHtml(formatMetric(entry.prom, 4)) + "</td>"
          + '<td class="num">' + escapeHtml(formatPercent(entry.frac || 0, 2)) + "</td>"
          + "<td>" + escapeHtml(entry.hits) + "</td>"
          + '<td><button type="button" data-discovery-focus-row="' + escapeHtml(entry.name) + '">Sort</button></td>'
          + "</tr>";
      }).join("");
      discoveryTableEl.innerHTML = header + "<tbody>" + body + "</tbody>";
      var rowButtons = discoveryTableEl.querySelectorAll("[data-discovery-focus-row]");
      for (var i = 0; i < rowButtons.length; i++) {
        rowButtons[i].addEventListener("click", function(evt) {
          evt.preventDefault();
          evt.stopPropagation();
          var name = String(this.getAttribute("data-discovery-focus-row") || "");
          activateLigandSort(name);
          renderCurrentMode();
        });
      }
      renderDiscoveryChips();
    }
""".strip()
