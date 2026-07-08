from __future__ import annotations


def render_core_js_module(safe_run_id: str) -> str:
    script = """
    function csvEscape(value) {
      var text = String(value == null ? "" : value);
      if (/[",\\n]/.test(text)) return '"' + text.replace(/"/g, '""') + '"';
      return text;
    }
    function downloadTextFile(filename, text, mimeType) {
      var blob = new Blob([text], { type: mimeType || "text/plain;charset=utf-8" });
      var url = URL.createObjectURL(blob);
      var anchor = document.createElement("a");
      anchor.href = url;
      anchor.download = filename;
      document.body.appendChild(anchor);
      anchor.click();
      setTimeout(function() {
        if (anchor.parentNode) anchor.parentNode.removeChild(anchor);
        URL.revokeObjectURL(url);
      }, 0);
    }
    function copyTextToClipboard(text) {
      if (navigator.clipboard && typeof navigator.clipboard.writeText === "function") {
        navigator.clipboard.writeText(text).catch(function() {
          window.prompt("Copy shareable link", text);
        });
        return;
      }
      window.prompt("Copy shareable link", text);
    }
    function scrollToNode(node) {
      if (!node || typeof node.scrollIntoView !== "function") return;
      node.scrollIntoView({ behavior: "smooth", block: "start" });
    }
    function navigateAnalysisTab(tabName) {
      activeAnalysisTab = tabName;
      renderAnalysisPanels();
      if (panelEl) scrollToNode(panelEl);
    }
    function buildReportState() {
      return {
        scale_mode: activeMode,
        threshold: Number(activeThreshold || 0.01),
        target_organization: String(targetOrganizationMode || "none"),
        ligand_organization: String(ligandOrganizationMode || "none"),
        safety_search: String(activeSafetySearch || ""),
        motif_search: String(activeMotifSearch || ""),
        interaction_mode: activeInteractionMode(),
        selected_pdb_label: String(selectedColLabel || ""),
        selected_pdb: String(selectedPdb || ""),
        selected_ligand: String(selectedLigand || ""),
        analysis_tab: String(activeAnalysisTab || "scatter"),
        scatter_mode: String(activeScatterMode || "rows"),
        rank_mode: String(activeRankMode || "rows"),
        rank_metric: String(activeRankMetric || "promiscuity_z"),
        discovery_mode: String(activeDiscoveryMode || "rows"),
        discovery_search: discoverySearchEl ? String(discoverySearchEl.value || "") : "",
        global_search: heroSearchEl ? String(heroSearchEl.value || "") : "",
        discovery_filter: String(activeDiscoveryFilter || "all"),
        discovery_highlight_mode: String(activeDiscoveryHighlightMode || ""),
        discovery_highlight_value: String(activeDiscoveryHighlightValue || ""),
        side_effect_overlap: !!activeSideEffectOverlap
      };
    }
    function applyReportState(state) {
      if (!state || typeof state !== "object") return;
      if (state.scale_mode === "quantile" && quantileAvailable) activeMode = "quantile";
      else if (state.scale_mode === "fixed") activeMode = "fixed";
      activeThreshold = Math.max(0.001, Math.min(0.1, Number(state.threshold) || activeThreshold));
      syncSideEffectThresholdControl();
      targetOrganizationMode = ["none", "safety", "adme", "family", "pathway"].indexOf(String(state.target_organization || "")) >= 0 ? String(state.target_organization) : targetOrganizationMode;
      ligandOrganizationMode = ["none", "chemotype", "drug_effect"].indexOf(String(state.ligand_organization || "")) >= 0 ? String(state.ligand_organization) : ligandOrganizationMode;
      activeSafetySearch = String(state.safety_search || "");
      if (safetySearchInputEl) safetySearchInputEl.value = activeSafetySearch;
      activeMotifSearch = String(state.motif_search || "");
      if (motifSearchInputEl) motifSearchInputEl.value = activeMotifSearch;
      activeAnalysisTab = ["scatter", "network", "rank"].indexOf(String(state.analysis_tab || "")) >= 0 ? String(state.analysis_tab) : activeAnalysisTab;
      activeScatterMode = state.scatter_mode === "cols" ? "cols" : "rows";
      activeRankMode = state.rank_mode === "cols" ? "cols" : "rows";
      activeRankMetric = ["promiscuity_z", "significant_hits", "hit_fraction", "entropy"].indexOf(String(state.rank_metric || "")) >= 0 ? String(state.rank_metric) : activeRankMetric;
      activeDiscoveryMode = "rows";
      activeDiscoveryFilter = "all";
      activeDiscoveryHighlightMode = "";
      activeDiscoveryHighlightValue = "";
      activeSideEffectOverlap = !!state.side_effect_overlap;
      if (discoverySearchEl) discoverySearchEl.value = String(state.discovery_search || "");
      if (heroSearchEl) heroSearchEl.value = String(state.global_search || "");
      resetHeatmapSelection();
      if (String(state.interaction_mode || "") === "pdb" && state.selected_pdb_label) {
        selectedColLabel = String(state.selected_pdb_label || "");
        selectedPdb = String(state.selected_pdb || "");
      } else if (String(state.interaction_mode || "") === "ligand" && state.selected_ligand) {
        selectedLigand = String(state.selected_ligand || "");
      }
      updateOrganizationControls();
    }
    function serializeReportStateToHash(state) {
      var params = new URLSearchParams();
      var payload = state || buildReportState();
      Object.keys(payload).forEach(function(key) {
        var value = payload[key];
        if (value == null || value === "" || value === false || value === "all") return;
        params.set(key, String(value));
      });
      return params.toString();
    }
    function parseReportStateFromHash() {
      var raw = String(window.location.hash || "").replace(/^#/, "");
      if (!raw) return null;
      var params = new URLSearchParams(raw);
      if (!params.toString()) return null;
      var state = {};
      params.forEach(function(value, key) { state[key] = value; });
      return state;
    }
    function syncReportStateHash() {
      if (suppressHashSync) return;
      var nextHash = serializeReportStateToHash(buildReportState());
      var nextUrl = window.location.pathname + window.location.search + (nextHash ? ("#" + nextHash) : "");
      if (typeof history.replaceState === "function") history.replaceState(null, "", nextUrl);
    }
    function hydrateReportStateFromLocation() {
      var parsed = parseReportStateFromHash();
      if (!parsed) return;
      suppressHashSync = true;
      try {
        applyReportState(parsed);
      } finally {
        suppressHashSync = false;
      }
    }
    function formatNumberText(value, digits) {
      var num = Number(value);
      if (!Number.isFinite(num)) return "NA";
      return num.toFixed(digits);
    }
    function validationItemHtml(label, value, note) {
      return '<div class="heatmap-validation-item"><div class="heatmap-validation-label">' + escapeHtml(label) + '</div><div class="heatmap-validation-value">' + escapeHtml(value) + '</div><div class="heatmap-validation-note">' + escapeHtml(note) + '</div></div>';
    }
    function renderHeroSummary() {
      if (!heroSummaryEl) return;
      var summary = analysisMeta.report_summary || {};
      var coverage = analysisMeta.coverage_summary || {};
      var assetMode = String(summary.asset_mode || "unknown");
      var bits = [
        '<span class="heatmap-report-pill"><strong>Run:</strong> ' + escapeHtml(String(summary.run_id || {safe_run_id!r})) + '</span>',
        '<span class="heatmap-report-pill"><strong>Generated:</strong> ' + escapeHtml(String(summary.generated_at || "unknown")) + '</span>',
        '<span class="heatmap-report-pill"><strong>Asset mode:</strong> ' + escapeHtml(assetMode) + '</span>',
        '<span class="heatmap-report-pill"><strong>Heatmap ligands:</strong> ' + escapeHtml(String(coverage.default_visible_ligands || 0)) + ' / ' + escapeHtml(String(coverage.interactive_ligands || 0)) + '</span>',
        '<span class="heatmap-report-pill"><strong>Full-run ligands:</strong> ' + escapeHtml(String(coverage.full_run_ligands || coverage.interactive_ligands || 0)) + '</span>',
        '<span class="heatmap-report-pill"><strong>Targets:</strong> ' + escapeHtml(String(coverage.interactive_targets || 0)) + ' / ' + escapeHtml(String(coverage.full_run_targets || 0)) + '</span>'
      ];
      heroSummaryEl.innerHTML = bits.join("");
    }
    function activeStateEntries() {
      var entries = [];
      if (selectedPdb) entries.push({ key: "pdb", label: "PDB", value: selectedPdb });
      if (selectedLigand) entries.push({ key: "ligand", label: "Ligand sort", value: selectedLigand });
      if (targetOrganizationMode !== "none") entries.push({ key: "target_org", label: "Target organization", value: organizationModeLabel(targetOrganizationMode) });
      if (activeSafetySearch) entries.push({ key: "ae_search", label: "Side-effect/AE search", value: activeSafetySearch });
      if (activeDiscoveryFilter !== "all") entries.push({ key: "discovery_filter", label: "Discovery filter", value: activeDiscoveryFilter === "high_safety" ? "High-confidence safety" : "Unresolved mappings" });
      if (activeDiscoveryHighlightMode) entries.push({ key: "discovery_org", label: "Discovery annotation", value: activeDiscoveryHighlightMode + (activeDiscoveryHighlightValue ? (" = " + activeDiscoveryHighlightValue) : "") });
      if (activeSideEffectOverlap) entries.push({ key: "side_effect_overlap", label: "Side-effect match", value: "highlighted significant pairs" });
      if (Number(activeThreshold || 0.01) !== Number(analysisMeta.default_threshold_pct_rank || 0.01)) entries.push({ key: "threshold", label: "Threshold", value: formatPercent(activeThreshold, 1) });
      return entries;
    }
    function clearActiveState(key) {
      if (key === "pdb" || key === "ligand") {
        resetHeatmapSelection();
        renderCurrentMode();
        return;
      }
      if (key === "target_org") {
        targetOrganizationMode = "none";
        updateOrganizationControls();
        renderCurrentMode();
        return;
      }
      if (key === "ae_search") {
        activeSafetySearch = "";
        if (safetySearchInputEl) safetySearchInputEl.value = "";
        renderAnnotationPanels();
        updateActiveLabelStates();
        refreshReportExtras(baseVizData(), [], []);
        syncReportStateHash();
        return;
      }
      if (key === "discovery_filter") {
        activeDiscoveryFilter = "all";
        renderDiscoveryPanel();
        renderActiveStateBar();
        syncReportStateHash();
        return;
      }
      if (key === "discovery_org") {
        activeDiscoveryHighlightMode = "";
        activeDiscoveryHighlightValue = "";
        renderDiscoveryPanel();
        renderActiveStateBar();
        syncReportStateHash();
        return;
      }
      if (key === "side_effect_overlap") {
        activeSideEffectOverlap = false;
        renderSideEffectOverlapPanel(baseVizData(), [], []);
        applySideEffectHighlights();
        renderActiveStateBar();
        syncReportStateHash();
        return;
      }
      if (key === "threshold") {
        activeThreshold = Math.max(0.001, Math.min(0.1, Number(analysisMeta.default_threshold_pct_rank) || 0.01));
        renderCurrentMode();
      }
    }
    function renderActiveStateBar() {
      if (!heroStateEl) return;
      var entries = activeStateEntries();
      heroStateEl.classList.toggle("is-empty", entries.length === 0);
      heroStateEl.innerHTML = entries.map(function(entry) {
        return '<span class="heatmap-active-chip"><strong>Active ' + escapeHtml(entry.label) + ':</strong> <span>' + escapeHtml(entry.value) + '</span> <button type="button" data-clear-state="' + escapeHtml(entry.key) + '">Reset</button></span>';
      }).join("");
      var clearButtons = heroStateEl.querySelectorAll("[data-clear-state]");
      for (var i = 0; i < clearButtons.length; i++) {
        clearButtons[i].addEventListener("click", function(evt) {
          evt.preventDefault();
          evt.stopPropagation();
          clearActiveState(String(this.getAttribute("data-clear-state") || ""));
        });
      }
    }
    function renderValidationPanel(payload, visibleRowIndices, visibleColIndices) {
      if (!validationPanelEl) return;
      var body = validationPanelEl.querySelector(".heatmap-report-panel-body");
      if (!body) return;
      var coverage = analysisMeta.coverage_summary || {};
      var validation = analysisMeta.validation_summary || {};
      var summary = analysisMeta.report_summary || {};
      var cache = validation.cache || {};
      var mapping = validation.missing_target_mappings || {};
      var drift = validation.safety_bucket_drift || {};
      var assetMode = String(summary.asset_mode || "unknown");
      var fullLigands = Number(coverage.full_run_ligands || allRowNames().length || 0);
      var interactiveLigands = Number(coverage.interactive_ligands || allRowNames().length || 0);
      var defaultVisibleLigands = Number(coverage.default_visible_ligands || visibleRowIndices.length || 0);
      var fullTargets = Number(coverage.full_run_targets || allColNames().length || 0);
      var interactiveTargets = Number(coverage.interactive_targets || allColNames().length || 0);
      var cacheValue = String(cache.status || "unknown").toUpperCase();
      if (cache.age_days != null) cacheValue += " (" + formatNumberText(cache.age_days, 1) + " days old)";
      var unresolvedTargets = Array.isArray(mapping.targets) ? mapping.targets : [];
      var mappingNote = unresolvedTargets.length ? unresolvedTargets.map(function(item) { return String(item.pdb || item.target_name || ""); }).join(", ") : "All visible targets resolved to UniProt.";
      body.innerHTML = '<div class="heatmap-validation-list">'
        + validationItemHtml("Coverage", "Heatmap default " + defaultVisibleLigands + " / " + interactiveLigands + " ligands; full run " + fullLigands, "Targets: interactive " + interactiveTargets + " / full run " + fullTargets)
        + validationItemHtml("Asset Mode", assetMode, "Clustergrammer asset packaging mode used for this report")
        + validationItemHtml("Safety Cache", cacheValue, String(cache.generated_at || "Cache timestamp unavailable"))
        + validationItemHtml("Target Mapping", String(mapping.count || 0) + " unresolved visible targets", mappingNote)
        + validationItemHtml("Bucket Drift", String(drift.changed_count || 0) + " changed targets on last refresh", String(drift.generated_at || "No drift summary available"))
        + '</div>';
    }
    function flattenSearchValues(value, out, depth) {
      out = out || [];
      depth = depth || 0;
      if (value == null || depth > 4) return out;
      if (Array.isArray(value)) {
        for (var i = 0; i < value.length; i++) flattenSearchValues(value[i], out, depth + 1);
        return out;
      }
      if (typeof value === "object") {
        Object.keys(value).forEach(function(key) {
          if (key === "metadata_json") return;
          flattenSearchValues(value[key], out, depth + 1);
        });
        return out;
      }
      out.push(String(value));
      return out;
    }
    function targetMetaSearchText(meta) {
      return flattenSearchValues(meta || {}, [], 0).join(" ").toLowerCase();
    }
""".strip()
    return script.replace("{safe_run_id!r}", repr(safe_run_id))
