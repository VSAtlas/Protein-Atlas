from __future__ import annotations


def render_side_effect_js_module() -> str:
    return """
    var SIDE_EFFECT_ALIAS_MAP = {
      "qt prolongation": "cardiotoxicity qt",
      "ventricular tachycardia": "cardiotoxicity qt",
      "myocardial infarction": "cardiotoxicity qt",
      "acute myocardial infarction": "cardiotoxicity qt",
      "myocarditis": "cardiotoxicity qt",
      "heart failure": "cardiotoxicity qt",
      "alanine aminotransferase increased": "hepatotoxicity",
      "hepatic function abnormal": "hepatotoxicity",
      "jaundice": "hepatotoxicity",
      "nausea": "gi",
      "diarrhoea": "gi",
      "diarrhea": "gi",
      "colitis ulcerative": "gi",
      "pancreatitis": "gi",
      "rash": "dermatologic",
      "photosensitivity reaction": "dermatologic",
      "anaemia": "hematologic",
      "anemia": "hematologic",
      "haemolytic anaemia": "hematologic",
      "hemolytic anemia": "hematologic",
      "cognitive disorder": "cns"
    };
    function isBlankSideEffectTerm(value) {
      var text = String(value || "").trim().toLowerCase();
      if (!text) return true;
      return text === "unassigned" || text === "other / unassigned" || text === "other unassigned" || text === "non-adme / unassigned";
    }
    function arrayFromMeta(meta, fields) {
      var values = [];
      if (!meta) return values;
      for (var i = 0; i < fields.length; i++) {
        var raw = meta[fields[i]];
        if (Array.isArray(raw)) {
          for (var j = 0; j < raw.length; j++) values.push(raw[j]);
        }
        else if (typeof raw === "string" && raw.trim()) {
          var split = raw.split(/[;|\\n]+/);
          for (var k = 0; k < split.length; k++) values.push(split[k]);
        }
        else if (raw != null && raw !== "") values.push(String(raw));
      }
      return values.map(function(value) { return String(value || "").trim(); }).filter(function(value) { return value && !isBlankSideEffectTerm(value); });
    }
    function canonicalSideEffectTerm(value) {
      var text = String(value || "").toLowerCase().replace(/&/g, " and ").replace(/[^a-z0-9]+/g, " ").trim().replace(/\\s+/g, " ");
      return SIDE_EFFECT_ALIAS_MAP[text] || text;
    }
    function sideEffectTermsForMeta(meta) {
      var fields = [
        "side_effects",
        "safety_evidence_terms",
        "raw_adverse_event_examples",
        "direct_liability_examples",
        "sider_side_effects",
        "openfda_side_effects",
        "drug_ae_buckets",
        "safety_buckets",
        "primary_display_safety",
        "secondary_safety_buckets",
        "direct_safety_buckets"
      ];
      var seen = Object.create(null);
      var out = [];
      arrayFromMeta(meta, fields).forEach(function(label) {
        var key = canonicalSideEffectTerm(label);
        if (!key || isBlankSideEffectTerm(key) || seen[key]) return;
        seen[key] = true;
        out.push({ key: key, label: label });
      });
      return out;
    }
    function numericMetaField(meta, fields) {
      if (!meta) return null;
      for (var i = 0; i < fields.length; i++) {
        var raw = meta[fields[i]];
        if (raw == null || raw === "") continue;
        var value = Number(raw);
        if (Number.isFinite(value)) return value;
      }
      return null;
    }
    function categoricalFactor(value, mapping) {
      var key = String(value || "").trim().toLowerCase();
      return Object.prototype.hasOwnProperty.call(mapping, key) ? Number(mapping[key]) : null;
    }
    function bindingPlausibilityFactor(cell) {
      var pct = asNumber(cell && cell.pct_rank);
      if (pct == null || pct > activeThreshold) return 0;
      if (pct <= 0) return 1;
      return Math.max(0.05, Math.min(1, 1 - (pct / Math.max(activeThreshold, 0.000001))));
    }
    function exposureFactor(rowMetaItem) {
      var freeCmax = numericMetaField(rowMetaItem, ["free_cmax_um", "cmax_free_um", "unbound_cmax_um"]);
      if (freeCmax == null) return { value: 1, label: "not available (neutral)", missing: false };
      var scaled = Math.max(0.05, Math.min(1, Math.log10(1 + Math.max(0, freeCmax)) / 2));
      return { value: scaled, label: formatMetric(freeCmax, 3) + " uM", missing: false };
    }
    function tissueExpressionFactor(colMetaItem) {
      var numeric = numericMetaField(colMetaItem, ["target_tissue_expression_score", "tissue_expression_score", "expression_score", "hpa_expression_score", "gtex_expression_score"]);
      if (numeric != null) {
        var levelLabel = String((colMetaItem && (colMetaItem.target_tissue_expression_label || colMetaItem.target_tissue_expression || colMetaItem.tissue_expression_level || colMetaItem.expression_level)) || "").trim();
        return { value: Math.max(0.05, Math.min(1, numeric)), label: (levelLabel ? levelLabel + " (" + formatMetric(numeric, 3) + ")" : formatMetric(numeric, 3)), missing: false };
      }
      var categorical = categoricalFactor(colMetaItem && (colMetaItem.target_tissue_expression || colMetaItem.tissue_expression_level || colMetaItem.expression_level), { "high": 1, "medium": 0.67, "low": 0.33, "not detected": 0.05 });
      if (categorical != null) return { value: categorical, label: formatMetric(categorical, 3), missing: false };
      return { value: 1, label: "not available (neutral)", missing: false };
    }
    function targetAdrEvidenceFactor(colMetaItem, matches) {
      var base = categoricalFactor(colMetaItem && colMetaItem.safety_confidence, { "high": 1, "medium": 0.7, "low": 0.4, "unassigned": 0.15 });
      if (base == null) base = Array.isArray(colMetaItem && colMetaItem.side_effects) && colMetaItem.side_effects.length ? 0.55 : 0.2;
      if (matches && matches.length > 1) base = Math.min(1, base + 0.1);
      return { value: base, label: String((colMetaItem && colMetaItem.safety_confidence) || "evidence"), missing: false };
    }
    function essentialityPathwayFactor(colMetaItem) {
      var numeric = numericMetaField(colMetaItem, ["essentiality_score", "target_essentiality_score", "pathway_role_score"]);
      if (numeric != null) return { value: Math.max(0.05, Math.min(1, numeric)), label: formatMetric(numeric, 3), missing: false };
      var core = String((colMetaItem && (colMetaItem.core_essentiality || colMetaItem.is_core_essential)) || "").trim().toLowerCase();
      if (["1", "true", "yes", "core essential", "essential"].indexOf(core) >= 0) return { value: 1, label: "essential", missing: false };
      var pathways = Array.isArray(colMetaItem && colMetaItem.pathway_memberships) ? colMetaItem.pathway_memberships : [];
      if (pathways.length) return { value: Math.min(0.85, 0.35 + Math.log(1 + pathways.length) / Math.log(12)), label: String(pathways.length) + " pathways", missing: false };
      return { value: 1, label: "missing", missing: true };
    }
    function visibleNamesFromPayload(payload, indices, axis) {
      var nodes = axis === "row" ? (payload && payload.row_nodes) : (payload && payload.col_nodes);
      var all = (Array.isArray(nodes) ? nodes : []).map(function(node) { return String((node && node.name) || ""); });
      if (!Array.isArray(indices) || !indices.length) return all;
      return indices.map(function(index) { return all[index]; }).filter(function(name) { return !!name; });
    }
    function visibleSideEffectAvailability(payload, visibleRowIndices, visibleColIndices) {
      var rows = visibleNamesFromPayload(payload, visibleRowIndices, "row");
      var cols = visibleNamesFromPayload(payload, visibleColIndices, "col");
      var hasLigand = rows.some(function(rowName) {
        return sideEffectTermsForMeta(lookupRowLabelMeta(rowName)).length > 0;
      });
      var hasTarget = cols.some(function(colName) {
        return sideEffectTermsForMeta(lookupColLabelMeta(colName)).length > 0;
      });
      return { hasLigand: hasLigand, hasTarget: hasTarget };
    }
    function sideEffectNoMatchMessage(payload, visibleRowIndices, visibleColIndices) {
      var availability = visibleSideEffectAvailability(payload, visibleRowIndices, visibleColIndices);
      if (!availability.hasLigand) return "No ligand side-effect terms are available in the visible view.";
      if (!availability.hasTarget) return "No target side-effect/adverse-event terms are available in the visible view.";
      return "No significant shared side-effect/adverse-event terms are present in the visible view at the current threshold.";
    }
    function sideEffectOverlapEntries(payload, visibleRowIndices, visibleColIndices) {
      var rows = visibleNamesFromPayload(payload, visibleRowIndices, "row");
      var cols = visibleNamesFromPayload(payload, visibleColIndices, "col");
      var rowTermsByName = Object.create(null);
      var colTermsByName = Object.create(null);
      var entries = [];
      for (var r = 0; r < rows.length; r++) {
        var rowName = rows[r];
        var rowMetaItem = lookupRowLabelMeta(rowName);
        rowTermsByName[rowName] = sideEffectTermsForMeta(rowMetaItem);
        for (var c = 0; c < cols.length; c++) {
          var colName = cols[c];
          var cell = lookupMeta(rowName, colName);
          if (!metaPassesThreshold(cell)) continue;
          var colMetaItem = lookupColLabelMeta(colName);
          if (!colTermsByName[colName]) colTermsByName[colName] = sideEffectTermsForMeta(colMetaItem);
          var targetTerms = colTermsByName[colName];
          var targetByKey = Object.create(null);
          targetTerms.forEach(function(term) { targetByKey[term.key] = term; });
          var matches = [];
          rowTermsByName[rowName].forEach(function(term) {
            if (!targetByKey[term.key]) return;
            matches.push({ key: term.key, ligand: term.label, target: targetByKey[term.key].label });
          });
          if (!matches.length) continue;
          var binding = bindingPlausibilityFactor(cell);
          var exposure = exposureFactor(rowMetaItem);
          var expression = tissueExpressionFactor(colMetaItem);
          var targetAdr = targetAdrEvidenceFactor(colMetaItem, matches);
          var essentiality = essentialityPathwayFactor(colMetaItem);
          var score = binding * exposure.value * expression.value * targetAdr.value * essentiality.value;
          var missing = [exposure, expression, essentiality].filter(function(item) { return item.missing; }).length;
          entries.push({
            row: rowName,
            col: colName,
            ligand: rowName,
            pdb: String(colMetaItem.pdb || targetShortName(colName, cell)),
            target: String(colMetaItem.target_name || ""),
            matches: matches,
            score: score,
            binding: binding,
            exposure: exposure,
            expression: expression,
            targetAdr: targetAdr,
            essentiality: essentiality,
            missing: missing,
            pct: asNumber(cell && cell.pct_rank)
          });
        }
      }
      return entries.sort(function(a, b) {
        if (b.score !== a.score) return b.score - a.score;
        if (a.missing !== b.missing) return a.missing - b.missing;
        return a.ligand < b.ligand ? -1 : (a.ligand > b.ligand ? 1 : 0);
      });
    }
    function applySideEffectHighlights() {
      if (!container) return;
      var keys = activeSideEffectOverlap ? (sideEffectOverlapKeySet || Object.create(null)) : Object.create(null);
      var tiles = container.querySelectorAll(".tile, .tile_up, .tile_dn");
      for (var i = 0; i < tiles.length; i++) {
        var tile = tiles[i];
        var d = tile.__data__ || {};
        var rowName = d.row_name || d.row || "";
        var colName = d.col_name || d.col || "";
        tile.classList.toggle("cg-side-effect-match", !!keys[tileKey(rowName, colName)]);
      }
    }
    function sideEffectSourceNote() {
      return "Formula scaffold: Atlas binding plausibility x free Cmax/exposure x target tissue expression x known target-ADR evidence x target essentiality/pathway role. Recognized ligand and target side-effect terms are exposed through side_effects, safety_evidence_terms, raw_adverse_event_examples, direct_liability_examples, drug_ae_buckets, and safety_buckets; optional exposure fields include free_cmax_um, cmax_um, and fraction_unbound_plasma. Missing exposure/expression use a neutral prior and are labeled not available.";
    }
    function syncSideEffectThresholdControl() {
      if (!sideEffectThresholdEl) return;
      var current = Math.max(0.001, Math.min(0.1, Number(activeThreshold) || 0.01));
      var options = Array.prototype.slice.call(sideEffectThresholdEl.options || []);
      var best = options[0] ? String(options[0].value) : String(current);
      var bestDelta = Infinity;
      options.forEach(function(option) {
        var value = Number(option.value);
        var delta = Math.abs(value - current);
        if (delta < bestDelta) {
          bestDelta = delta;
          best = String(option.value);
        }
      });
      sideEffectThresholdEl.value = best;
    }
    function renderSideEffectOverlapPanel(payload, visibleRowIndices, visibleColIndices) {
      syncSideEffectThresholdControl();
      var hadSideEffectHighlights = Object.keys(sideEffectOverlapKeySet || {}).length > 0;
      sideEffectOverlapKeySet = Object.create(null);
      if (sideEffectOverlapButtonEl) sideEffectOverlapButtonEl.classList.toggle("is-active", activeSideEffectOverlap);
      if (!activeSideEffectOverlap) {
        if (sideEffectOverlapSummaryEl) sideEffectOverlapSummaryEl.textContent = "";
        if (sideEffectOverlapPanelEl) {
          sideEffectOverlapPanelEl.classList.remove("is-active");
          sideEffectOverlapPanelEl.innerHTML = "";
        }
        if (hadSideEffectHighlights) applySideEffectHighlights();
        return;
      }
      var entries = sideEffectOverlapEntries(payload || baseVizData(), visibleRowIndices, visibleColIndices);
      entries.forEach(function(entry) { sideEffectOverlapKeySet[tileKey(entry.row, entry.col)] = true; });
      if (sideEffectOverlapSummaryEl) {
        sideEffectOverlapSummaryEl.textContent = entries.length
          ? String(entries.length) + " significant protein-drug pairs share at least one known side-effect/adverse-event term at pct_rank <= " + formatPercent(activeThreshold, 1) + "."
          : sideEffectNoMatchMessage(payload || baseVizData(), visibleRowIndices, visibleColIndices);
      }
      if (!sideEffectOverlapPanelEl) return;
      sideEffectOverlapPanelEl.classList.add("is-active");
      if (!entries.length) {
        sideEffectOverlapPanelEl.innerHTML = '<div class="heatmap-side-effect-note">' + escapeHtml(sideEffectNoMatchMessage(payload || baseVizData(), visibleRowIndices, visibleColIndices) + " " + sideEffectSourceNote()) + '</div>';
        return;
      }
      var limited = entries.slice(0, 40);
      var rowsHtml = limited.map(function(entry) {
        var terms = entry.matches.map(function(match) { return match.ligand === match.target ? match.ligand : (match.ligand + " / " + match.target); }).join("; ");
        return "<tr>"
          + "<td>" + escapeHtml(entry.ligand) + "</td>"
          + "<td>" + escapeHtml(entry.pdb) + "</td>"
          + "<td>" + escapeHtml(terms) + "</td>"
          + '<td class="num">' + escapeHtml(formatMetric(entry.score, 4)) + "</td>"
          + '<td class="num">' + escapeHtml(formatMetric(entry.binding, 3)) + "</td>"
          + "<td>" + escapeHtml(entry.exposure.label) + "</td>"
          + "<td>" + escapeHtml(entry.expression.label) + "</td>"
          + "<td>" + escapeHtml(entry.targetAdr.label) + "</td>"
          + "<td>" + escapeHtml(entry.essentiality.label) + "</td>"
          + "</tr>";
      }).join("");
      sideEffectOverlapPanelEl.innerHTML = '<div class="heatmap-side-effect-table-wrap"><table class="heatmap-side-effect-table"><thead><tr><th>Ligand</th><th>PDB</th><th>Matched side effects / AEs</th><th>Hypothesis score</th><th>Binding</th><th>Free Cmax</th><th>Tissue expression</th><th>Target-ADR</th><th>Essentiality / pathway</th></tr></thead><tbody>'
        + rowsHtml
        + '</tbody></table></div><div class="heatmap-side-effect-note">' + escapeHtml(sideEffectSourceNote()) + '</div>';
    }
""".strip()
