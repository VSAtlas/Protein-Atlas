suppressPackageStartupMessages(library(pheatmap))

parse_args <- function(args) {
  opts <- list()
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (startsWith(arg, "--")) {
      key <- sub("^--", "", arg)
      if (grepl("=", key, fixed = TRUE)) {
        parts <- strsplit(key, "=", fixed = TRUE)[[1]]
        key <- parts[1]
        value <- parts[2]
      } else if (i < length(args) && !startsWith(args[i + 1], "--")) {
        value <- args[i + 1]
        i <- i + 1
      } else {
        value <- "true"
      }
      opts[[key]] <- value
    }
    i <- i + 1
  }
  opts
}

args <- commandArgs(trailingOnly = TRUE)
opts <- parse_args(args)

repo_root <- if (!is.null(opts[["repo-root"]])) opts[["repo-root"]] else "."
run_id <- if (!is.null(opts[["run-id"]])) opts[["run-id"]] else ""
if (!nzchar(run_id)) {
  stop("Missing required --run-id")
}

read_config_value <- function(path, key) {
  if (!file.exists(path)) {
    return(NULL)
  }
  lines <- readLines(path, warn = FALSE)
  for (line in lines) {
    line <- trimws(line)
    if (!nzchar(line) || startsWith(line, "#")) {
      next
    }
    line <- sub("\\s+#.*$", "", line)
    if (!nzchar(line)) {
      next
    }
    if (grepl("=", line, fixed = TRUE)) {
      parts <- strsplit(line, "=", fixed = TRUE)[[1]]
      k <- trimws(parts[1])
      v <- trimws(paste(parts[-1], collapse = "="))
      if (toupper(k) == toupper(key)) {
        return(v)
      }
    }
  }
  NULL
}

parse_bool <- function(value, default_value = TRUE) {
  if (is.null(value)) {
    return(default_value)
  }
  v <- tolower(trimws(value))
  if (v %in% c("1", "true", "yes", "on")) {
    return(TRUE)
  }
  if (v %in% c("0", "false", "no", "off")) {
    return(FALSE)
  }
  default_value
}

parse_num <- function(value, default_value) {
  if (is.null(value)) {
    return(default_value)
  }
  v <- suppressWarnings(as.numeric(trimws(value)))
  if (is.na(v)) {
    return(default_value)
  }
  v
}

parse_scale_mode <- function(value) {
  if (is.null(value)) {
    return("fixed")
  }
  v <- tolower(trimws(value))
  if (v %in% c("fixed", "quantile")) {
    return(v)
  }
  "fixed"
}

parse_quantile_value <- function(value, default_value) {
  parsed <- parse_num(value, default_value)
  if (is.na(parsed) || !is.finite(parsed)) {
    return(default_value)
  }
  if (parsed < 0) {
    return(0)
  }
  if (parsed > 1) {
    return(1)
  }
  parsed
}

parse_quantile_params <- function(config_path) {
  q_low <- parse_quantile_value(
    read_config_value(config_path, "HEATMAP_QUANTILE_LOW"),
    0.02
  )
  q_mid2 <- parse_quantile_value(
    read_config_value(config_path, "HEATMAP_QUANTILE_MID2"),
    0.85
  )
  q_high <- parse_quantile_value(
    read_config_value(config_path, "HEATMAP_QUANTILE_HIGH"),
    0.98
  )
  if (!(q_low < q_mid2 && q_mid2 < q_high)) {
    q_low <- 0.02
    q_mid2 <- 0.85
    q_high <- 0.98
  }
  list(q_low = q_low, q_mid2 = q_mid2, q_high = q_high)
}

is_increasing_breaks <- function(breaks) {
  isTRUE(
    breaks$scale_min < breaks$scale_mid &&
      breaks$scale_mid < breaks$scale_mid2 &&
      breaks$scale_mid2 < breaks$scale_max
  )
}

widen_breaks_epsilon <- function(breaks) {
  if (is_increasing_breaks(breaks)) {
    return(list(ok = TRUE, breaks = breaks, widened = FALSE))
  }
  start <- breaks$scale_min
  end <- breaks$scale_max
  if (!is.finite(start) || !is.finite(end) || end <= start) {
    return(list(ok = FALSE, breaks = breaks, widened = FALSE))
  }
  eps <- max(abs(end - start) * 1e-6, 1e-9)
  widened <- breaks
  if (widened$scale_mid <= widened$scale_min) {
    widened$scale_mid <- widened$scale_min + eps
  }
  if (widened$scale_mid2 <= widened$scale_mid) {
    widened$scale_mid2 <- widened$scale_mid + eps
  }
  if (widened$scale_max <= widened$scale_mid2) {
    widened$scale_max <- widened$scale_mid2 + eps
  }
  if (!is_increasing_breaks(widened)) {
    return(list(ok = FALSE, breaks = breaks, widened = FALSE))
  }
  list(ok = TRUE, breaks = widened, widened = TRUE)
}

format_breaks <- function(breaks) {
  sprintf(
    "min=%.6g mid=%.6g mid2=%.6g max=%.6g",
    breaks$scale_min,
    breaks$scale_mid,
    breaks$scale_mid2,
    breaks$scale_max
  )
}

read_color <- function(value, default_value) {
  if (is.null(value)) {
    return(default_value)
  }
  v <- trimws(value)
  if (!nzchar(v)) {
    return(default_value)
  }
  v
}

config_path <- file.path(repo_root, "config.txt")
use_dendrogram <- parse_bool(read_config_value(config_path, "USE_DENDROGRAM"), TRUE)
requested_scale_mode <- parse_scale_mode(read_config_value(config_path, "HEATMAP_SCALE_MODE"))
quantile_params <- parse_quantile_params(config_path)

scale_min <- parse_num(read_config_value(config_path, "HEATMAP_SCALE_MIN"), -2)
scale_mid <- parse_num(read_config_value(config_path, "HEATMAP_SCALE_MID"), 0)
scale_max <- parse_num(read_config_value(config_path, "HEATMAP_SCALE_MAX"), 2)
if (!(scale_min < scale_mid && scale_mid < scale_max)) {
  scale_min <- -2
  scale_mid <- 0
  scale_max <- 2
}

scale_mid2 <- parse_num(
  read_config_value(config_path, "HEATMAP_SCALE_MID2"),
  (scale_mid + scale_max) / 2
)
if (!(scale_mid < scale_mid2 && scale_mid2 < scale_max)) {
  scale_mid2 <- (scale_mid + scale_max) / 2
}
fixed_breaks <- list(
  scale_min = scale_min,
  scale_mid = scale_mid,
  scale_mid2 = scale_mid2,
  scale_max = scale_max
)

color_low <- read_color(read_config_value(config_path, "HEATMAP_COLOR_MIN"), "#2166ac")
color_mid <- read_color(read_config_value(config_path, "HEATMAP_COLOR_MID"), "#f7f7f7")
color_mid2 <- read_color(read_config_value(config_path, "HEATMAP_COLOR_MID2"), "#f4a582")
color_high <- read_color(read_config_value(config_path, "HEATMAP_COLOR_MAX"), "#b2182b")

input_path <- NULL
if (!is.null(opts[["in"]]) && nzchar(opts[["in"]])) {
  input_path <- opts[["in"]]
} else {
  heatmap_path <- file.path(repo_root, "data", run_id, "heatmap_input.csv")
  master_path <- file.path(repo_root, "data", run_id, "master_rows.csv")
  if (file.exists(heatmap_path)) {
    input_path <- heatmap_path
  } else {
    input_path <- master_path
  }
}

out_path <- if (!is.null(opts[["out"]])) {
  opts[["out"]]
} else {
  file.path(repo_root, "data", run_id, "heatmap.png")
}

top_k <- if (!is.null(opts[["top-k"]])) as.integer(opts[["top-k"]]) else 100L
if (is.na(top_k) || top_k <= 0) {
  top_k <- 100L
}

include_decoys <- FALSE
if (!is.null(opts[["include-decoys"]])) {
  flag <- tolower(opts[["include-decoys"]])
  include_decoys <- flag %in% c("1", "true", "yes", "on")
}

if (!file.exists(input_path)) {
  stop(sprintf("Input CSV not found: %s", input_path))
}

df <- read.csv(input_path, stringsAsFactors = FALSE)
if (!("t_selected" %in% names(df))) {
  stop("Missing required column: t_selected")
}
has_target_id <- "target_id" %in% names(df)
if (!has_target_id) {
  required_cols <- c("pdb_id", "variant", "ph_label")
  missing_cols <- required_cols[!(required_cols %in% names(df))]
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }
}

df$t_selected <- suppressWarnings(as.numeric(df$t_selected))
df <- df[!is.na(df$t_selected), ]
if (nrow(df) == 0) {
  stop("No usable t_selected values found in input CSV")
}

if (!include_decoys && "is_decoy" %in% names(df)) {
  decoy_flag <- tolower(as.character(df$is_decoy)) %in% c("1", "true", "yes", "on")
  df <- df[!decoy_flag, ]
}

ligand_display <- if ("ligand_display" %in% names(df)) df$ligand_display else ""
ligand_display <- trimws(as.character(ligand_display))
if ("ligand_base" %in% names(df)) {
  ligand_base <- trimws(as.character(df$ligand_base))
} else {
  ligand_base <- rep("", nrow(df))
}
if (!("ligand_display" %in% names(df)) && !("ligand_base" %in% names(df))) {
  stop("Missing ligand_display/ligand_base columns in input CSV")
}
missing_display <- ligand_display == "" | tolower(ligand_display) == "nan"
ligand_display[missing_display] <- ligand_base[missing_display]
df$ligand_name <- ligand_display

if (has_target_id) {
  df$target_id <- trimws(as.character(df$target_id))
} else {
  df$target_id <- paste(df$pdb_id, df$variant, df$ph_label, sep = "|")
}

# Aggregate duplicate target/ligand pairs with max t_selected.
agg <- aggregate(
  t_selected ~ target_id + ligand_name,
  data = df,
  FUN = function(x) max(x, na.rm = TRUE)
)

# Select top-k ligands by global max t_selected (deterministic tie-breaker by name).
ligand_max <- aggregate(
  t_selected ~ ligand_name,
  data = agg,
  FUN = function(x) max(x, na.rm = TRUE)
)
ligand_max <- ligand_max[order(-ligand_max$t_selected, ligand_max$ligand_name), ]
if (nrow(ligand_max) == 0) {
  stop("No ligands available after filtering")
}
top_k <- min(top_k, nrow(ligand_max))
top_ligands <- ligand_max$ligand_name[seq_len(top_k)]
top_ligands <- sort(top_ligands)

agg <- agg[agg$ligand_name %in% top_ligands, ]
agg$ligand_name <- factor(agg$ligand_name, levels = top_ligands)
target_levels <- sort(unique(agg$target_id))
agg$target_id <- factor(agg$target_id, levels = target_levels)

mat <- tapply(
  agg$t_selected,
  list(agg$ligand_name, agg$target_id),
  function(x) if (length(x)) max(x, na.rm = TRUE) else NA_real_
)
mat <- as.matrix(mat)

row_keep <- apply(mat, 1, function(x) any(is.finite(x)))
col_keep <- apply(mat, 2, function(x) any(is.finite(x)))
mat <- mat[row_keep, col_keep, drop = FALSE]

if (nrow(mat) == 0 || ncol(mat) == 0) {
  stop("Heatmap matrix is empty after filtering")
}

finite_values <- as.numeric(mat[is.finite(mat)])
n_finite <- length(finite_values)
quantile_breaks <- NULL
scale_mode <- "fixed"
fallback_reason <- NULL
quantile_widened <- FALSE

if (requested_scale_mode == "quantile") {
  if (n_finite < 10) {
    fallback_reason <- "n_finite_lt_10"
  } else {
    quantile_breaks <- list(
      scale_min = as.numeric(quantile(finite_values, quantile_params$q_low, names = FALSE, type = 7)),
      scale_mid = as.numeric(quantile(finite_values, 0.5, names = FALSE, type = 7)),
      scale_mid2 = as.numeric(quantile(finite_values, quantile_params$q_mid2, names = FALSE, type = 7)),
      scale_max = as.numeric(quantile(finite_values, quantile_params$q_high, names = FALSE, type = 7))
    )
    if (quantile_breaks$scale_min == quantile_breaks$scale_max) {
      fallback_reason <- "quantile_min_eq_max"
    } else {
      widened <- widen_breaks_epsilon(quantile_breaks)
      if (!isTRUE(widened$ok)) {
        fallback_reason <- "quantile_non_increasing_breaks"
      } else {
        quantile_breaks <- widened$breaks
        quantile_widened <- isTRUE(widened$widened)
        scale_mode <- "quantile"
      }
    }
  }
}

chosen_breaks <- fixed_breaks
if (scale_mode == "quantile" && !is.null(quantile_breaks)) {
  chosen_breaks <- quantile_breaks
}
scale_min <- chosen_breaks$scale_min
scale_mid <- chosen_breaks$scale_mid
scale_mid2 <- chosen_breaks$scale_mid2
scale_max <- chosen_breaks$scale_max

if (!is.null(fallback_reason) && requested_scale_mode == "quantile") {
  warning(
    sprintf(
      "heatmap action=scale_breaks_fallback requested_mode=quantile reason=%s n=%d q_low=%.4f q_mid2=%.4f q_high=%.4f using=fixed breaks=%s",
      fallback_reason,
      n_finite,
      quantile_params$q_low,
      quantile_params$q_mid2,
      quantile_params$q_high,
      format_breaks(fixed_breaks)
    )
  )
} else {
  message(
    sprintf(
      "heatmap action=scale_breaks mode=%s n=%d q_low=%.4f q_mid2=%.4f q_high=%.4f breaks=%s",
      scale_mode,
      n_finite,
      quantile_params$q_low,
      quantile_params$q_mid2,
      quantile_params$q_high,
      format_breaks(chosen_breaks)
    )
  )
}

total_steps <- 100L
span_low <- max(scale_mid - scale_min, 0)
span_mid <- max(scale_mid2 - scale_mid, 0)
span_high <- max(scale_max - scale_mid2, 0)
total_span <- span_low + span_mid + span_high
if (total_span <= 0) {
  span_low <- 1
  span_mid <- 1
  span_high <- 1
  total_span <- 3
}
n_low <- max(1, round(total_steps * span_low / total_span))
n_mid <- max(1, round(total_steps * span_mid / total_span))
n_high <- total_steps - n_low - n_mid
if (n_high < 1) {
  n_high <- 1
  if (n_mid > 1) {
    n_mid <- n_mid - 1
  } else if (n_low > 1) {
    n_low <- n_low - 1
  }
}

breaks_low <- seq(scale_min, scale_mid, length.out = n_low + 1)
breaks_mid <- seq(scale_mid, scale_mid2, length.out = n_mid + 1)
breaks_high <- seq(scale_mid2, scale_max, length.out = n_high + 1)
heatmap_breaks <- c(breaks_low, breaks_mid[-1], breaks_high[-1])
heatmap_colors <- c(
  colorRampPalette(c(color_low, color_mid))(n_low),
  colorRampPalette(c(color_mid, color_mid2))(n_mid),
  colorRampPalette(c(color_mid2, color_high))(n_high)
)

mat_plot <- mat
finite_mask <- is.finite(mat_plot)
mat_plot[finite_mask] <- pmin(pmax(mat_plot[finite_mask], scale_min), scale_max)

cluster_rows <- FALSE
cluster_cols <- FALSE
if (use_dendrogram) {
  mat_for_cluster <- mat
  if (any(!is.finite(mat_for_cluster))) {
    min_val <- min(mat_for_cluster[is.finite(mat_for_cluster)], na.rm = TRUE)
    mat_for_cluster[!is.finite(mat_for_cluster)] <- min_val
  }
  cluster_rows <- if (nrow(mat) > 1) hclust(dist(mat_for_cluster)) else FALSE
  cluster_cols <- if (ncol(mat) > 1) hclust(dist(t(mat_for_cluster))) else FALSE
}
plot_title <- if (use_dendrogram) {
  "t_selected heatmap (rows=ligands, cols=targets; dendrograms show similarity)"
} else {
  "t_selected heatmap (rows=ligands, cols=targets)"
}

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
png(out_path, width = 1600, height = 2000)
pheatmap(
  mat_plot,
  cluster_rows = cluster_rows,
  cluster_cols = cluster_cols,
  scale = "none",
  color = heatmap_colors,
  breaks = heatmap_breaks,
  na_col = "white",
  main = plot_title
)
dev.off()

sidecar_path <- paste0(out_path, ".scale.txt")
sidecar_lines <- c(
  sprintf("scale_mode=%s", scale_mode),
  sprintf("requested_scale_mode=%s", requested_scale_mode),
  sprintf("n_finite=%d", n_finite),
  sprintf("q_low=%.10g", quantile_params$q_low),
  sprintf("q_mid2=%.10g", quantile_params$q_mid2),
  sprintf("q_high=%.10g", quantile_params$q_high),
  sprintf("chosen_scale_min=%.10g", chosen_breaks$scale_min),
  sprintf("chosen_scale_mid=%.10g", chosen_breaks$scale_mid),
  sprintf("chosen_scale_mid2=%.10g", chosen_breaks$scale_mid2),
  sprintf("chosen_scale_max=%.10g", chosen_breaks$scale_max),
  sprintf("fixed_scale_min=%.10g", fixed_breaks$scale_min),
  sprintf("fixed_scale_mid=%.10g", fixed_breaks$scale_mid),
  sprintf("fixed_scale_mid2=%.10g", fixed_breaks$scale_mid2),
  sprintf("fixed_scale_max=%.10g", fixed_breaks$scale_max),
  sprintf("quantile_widened=%s", ifelse(quantile_widened, "true", "false")),
  sprintf("fallback_reason=%s", ifelse(is.null(fallback_reason), "", fallback_reason))
)
if (!is.null(quantile_breaks)) {
  sidecar_lines <- c(
    sidecar_lines,
    sprintf("quantile_scale_min=%.10g", quantile_breaks$scale_min),
    sprintf("quantile_scale_mid=%.10g", quantile_breaks$scale_mid),
    sprintf("quantile_scale_mid2=%.10g", quantile_breaks$scale_mid2),
    sprintf("quantile_scale_max=%.10g", quantile_breaks$scale_max)
  )
}
writeLines(sidecar_lines, con = sidecar_path)
message(sprintf("heatmap action=scale_sidecar path=%s", sidecar_path))
