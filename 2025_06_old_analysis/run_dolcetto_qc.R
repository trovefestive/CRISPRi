#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(htmltools)
})

options(stringsAsFactors = FALSE, scipen = 999)

parse_args <- function(args) {
  out <- list(
    root = getwd(),
    out_dir = "analysis/results",
    sample = "both",
    chunk_reads = 250000L,
    max_reads = Inf,
    skip_md5 = FALSE
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--skip-md5") {
      out$skip_md5 <- TRUE
      i <- i + 1L
      next
    }
    if (i == length(args)) stop("Missing value after ", key)
    value <- args[[i + 1L]]
    if (key == "--root") out$root <- value
    else if (key == "--out-dir") out$out_dir <- value
    else if (key == "--sample") out$sample <- value
    else if (key == "--chunk-reads") out$chunk_reads <- as.integer(value)
    else if (key == "--max-reads") out$max_reads <- as.numeric(value)
    else stop("Unknown argument: ", key)
    i <- i + 2L
  }
  if (!out$sample %in% c("A", "B", "both")) stop("--sample must be A, B, or both")
  if (is.na(out$chunk_reads) || out$chunk_reads < 1L) stop("--chunk-reads must be positive")
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
root <- normalizePath(args$root, mustWork = TRUE)
out_dir <- if (grepl("^/", args$out_dir)) args$out_dir else file.path(root, args$out_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "diagnostics"), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "run.log")
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

sample_config <- data.table(
  sample = c("A", "B"),
  fastq_r1 = file.path(root, c(
    "LIbrary_QC_Novogene/01.RawData/sg_set_A/sg_set_A_CKDL260013754-1A_23KG5MLT4_L3_1.fq.gz",
    "LIbrary_QC_Novogene/01.RawData/sg_set_B/sg_set_B_CKDL260013755-1A_23KG5MLT4_L3_1.fq.gz"
  )),
  fastq_r2 = file.path(root, c(
    "LIbrary_QC_Novogene/01.RawData/sg_set_A/sg_set_A_CKDL260013754-1A_23KG5MLT4_L3_2.fq.gz",
    "LIbrary_QC_Novogene/01.RawData/sg_set_B/sg_set_B_CKDL260013755-1A_23KG5MLT4_L3_2.fq.gz"
  )),
  reference = file.path(root, c(
    "Library target genes/broadgpp-dolcetto-targets-seta.txt",
    "Library target genes/broadgpp-dolcetto-targets-setb.txt"
  )),
  md5_manifest = file.path(root, c(
    "LIbrary_QC_Novogene/01.RawData/sg_set_A/MD5.txt",
    "LIbrary_QC_Novogene/01.RawData/sg_set_B/MD5.txt"
  ))
)

selected <- if (args$sample == "both") c("A", "B") else args$sample
sample_config <- sample_config[sample %in% selected]
all_input_files <- unique(c(sample_config$fastq_r1, sample_config$fastq_r2, sample_config$reference))
missing_files <- all_input_files[!file.exists(all_input_files)]
if (length(missing_files)) stop("Missing input files:\n", paste(missing_files, collapse = "\n"))

read_manifest <- function(path) {
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines)]
  bits <- strsplit(lines, "[[:space:]]+")
  data.table(
    expected_md5 = vapply(bits, `[[`, character(1), 1L),
    basename = basename(vapply(bits, function(x) x[[length(x)]], character(1)))
  )
}

if (!args$skip_md5) {
  log_msg("Calculating input MD5 checksums")
  checksum_rows <- list()
  for (i in seq_len(nrow(sample_config))) {
    cfg <- sample_config[i]
    manifest <- read_manifest(cfg$md5_manifest)
    files <- c(cfg$fastq_r1, cfg$fastq_r2)
    actual <- unname(tools::md5sum(files))
    dt <- data.table(sample = cfg$sample, file = files, basename = basename(files), actual_md5 = actual)
    dt <- manifest[dt, on = "basename"]
    dt[, matches_manifest := actual_md5 == expected_md5]
    checksum_rows[[length(checksum_rows) + 1L]] <- dt
  }
  for (i in seq_len(nrow(sample_config))) {
    cfg <- sample_config[i]
    checksum_rows[[length(checksum_rows) + 1L]] <- data.table(
      sample = cfg$sample,
      file = cfg$reference,
      basename = basename(cfg$reference),
      expected_md5 = NA_character_,
      actual_md5 = unname(tools::md5sum(cfg$reference)),
      matches_manifest = NA
    )
  }
  checksums <- rbindlist(checksum_rows, fill = TRUE)
  fwrite(checksums, file.path(out_dir, "tables", "input_checksums.csv"))
  if (any(checksums$matches_manifest %in% FALSE, na.rm = TRUE)) stop("FASTQ MD5 validation failed")
  log_msg("MD5 validation passed")
}

gini_coefficient <- function(x) {
  x <- as.numeric(x)
  x[is.na(x) | x < 0] <- 0
  if (!length(x) || sum(x) == 0) return(NA_real_)
  x <- sort(x)
  n <- length(x)
  sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x))
}

add_named_counts <- function(target, values) {
  if (!length(values)) return(target)
  tab <- table(values)
  new_names <- setdiff(names(tab), names(target))
  if (length(new_names)) target[new_names] <- 0
  target[names(tab)] <- target[names(tab)] + as.numeric(tab)
  target
}

build_hamming1_lookup <- function(refs) {
  log_msg("Building one-substitution diagnostic lookup for ", length(refs), " reference guides")
  bases <- c("A", "C", "G", "T")
  alt_map <- list(
    A = c("C", "G", "T"),
    C = c("A", "G", "T"),
    G = c("A", "C", "T"),
    T = c("A", "C", "G")
  )
  pieces <- vector("list", 20L)
  for (pos in seq_len(20L)) {
    original <- substring(refs, pos, pos)
    variants <- rep(refs, each = 3L)
    replacements <- unname(unlist(alt_map[original], use.names = FALSE))
    substring(variants, pos, pos) <- replacements
    pieces[[pos]] <- variants
  }
  variants <- unlist(pieces, use.names = FALSE)
  rm(pieces)
  lookup <- data.table(guide = variants)[, .(n_reference_neighbors = .N), by = guide]
  setkey(lookup, guide)
  rm(variants)
  gc(verbose = FALSE)
  log_msg("One-substitution lookup contains ", format(nrow(lookup), big.mark = ","), " unique sequences")
  lookup
}

process_sample <- function(cfg, max_reads, chunk_reads) {
  sample_name <- cfg$sample
  log_msg("Starting Set ", sample_name)
  ref <- fread(cfg$reference, sep = "\t", header = TRUE)
  if (ncol(ref) != 3L) stop("Reference for Set ", sample_name, " must have three columns")
  setnames(ref, c("guide", "gene_symbol", "gene_id"))
  ref[, guide := toupper(guide)]
  if (any(nchar(ref$guide) != 20L) || anyDuplicated(ref$guide)) stop("Invalid or duplicated guides in Set ", sample_name)
  ref[, guide_type := ifelse(gene_symbol == "NO-TARGET", "Non-targeting control", "Gene-targeting")]

  mismatch_lookup <- build_hamming1_lookup(ref$guide)
  exact_counts <- numeric(nrow(ref))
  names(exact_counts) <- ref$guide
  categories <- c(
    total_reads = 0,
    no_CACCG_prefix = 0,
    valid_20nt_insert = 0,
    invalid_base_20nt_insert = 0,
    non20_insert_with_GTTT = 0,
    no_downstream_GTTT = 0,
    exact_reference = 0,
    nonreference_20nt = 0,
    one_mismatch_unique = 0,
    one_mismatch_ambiguous = 0,
    more_than_one_mismatch_or_other = 0
  )
  insert_lengths <- numeric()
  stagger_positions <- numeric()
  read_lengths <- numeric()
  top_nonref_chunks <- list()
  total_processed <- 0
  chunk_number <- 0L
  started <- Sys.time()

  con <- gzfile(cfg$fastq_r1, open = "rt")
  on.exit(close(con), add = TRUE)
  repeat {
    if (is.finite(max_reads) && total_processed >= max_reads) break
    reads_this_chunk <- chunk_reads
    if (is.finite(max_reads)) reads_this_chunk <- min(reads_this_chunk, max_reads - total_processed)
    lines <- readLines(con, n = as.integer(reads_this_chunk * 4), warn = FALSE)
    if (!length(lines)) break
    if (length(lines) %% 4L != 0L) stop("Truncated FASTQ record in Set ", sample_name)
    seqs <- toupper(lines[seq.int(2L, length(lines), by = 4L)])
    headers <- lines[seq.int(1L, length(lines), by = 4L)]
    plus <- lines[seq.int(3L, length(lines), by = 4L)]
    if (any(!startsWith(headers, "@")) || any(!startsWith(plus, "+"))) stop("Malformed FASTQ structure in Set ", sample_name)
    rm(lines, headers, plus)

    n <- length(seqs)
    total_processed <- total_processed + n
    chunk_number <- chunk_number + 1L
    categories["total_reads"] <- categories["total_reads"] + n
    read_lengths <- add_named_counts(read_lengths, nchar(seqs))

    prefix_pos <- regexpr("CACCG", seqs, fixed = TRUE)
    has_prefix <- prefix_pos > 0L
    categories["no_CACCG_prefix"] <- categories["no_CACCG_prefix"] + sum(!has_prefix)
    stagger_positions <- add_named_counts(stagger_positions, prefix_pos[has_prefix])

    candidate_idx <- which(has_prefix)
    guide_raw <- rep(NA_character_, n)
    expected_suffix <- rep(FALSE, n)
    if (length(candidate_idx)) {
      guide_raw[candidate_idx] <- substring(seqs[candidate_idx], prefix_pos[candidate_idx] + 5L, prefix_pos[candidate_idx] + 24L)
      suffix <- substring(seqs[candidate_idx], prefix_pos[candidate_idx] + 25L, prefix_pos[candidate_idx] + 28L)
      expected_suffix[candidate_idx] <- suffix == "GTTT"
    }
    valid_bases <- expected_suffix & grepl("^[ACGT]{20}$", guide_raw)
    invalid_bases <- expected_suffix & !valid_bases
    categories["valid_20nt_insert"] <- categories["valid_20nt_insert"] + sum(valid_bases)
    categories["invalid_base_20nt_insert"] <- categories["invalid_base_20nt_insert"] + sum(invalid_bases)

    unexpected_idx <- which(has_prefix & !expected_suffix)
    if (length(unexpected_idx)) {
      after_prefix <- substring(seqs[unexpected_idx], prefix_pos[unexpected_idx] + 5L)
      downstream <- regexpr("GTTT", after_prefix, fixed = TRUE)
      found <- downstream > 0L
      lengths_found <- downstream[found] - 1L
      insert_lengths <- add_named_counts(insert_lengths, lengths_found)
      categories["non20_insert_with_GTTT"] <- categories["non20_insert_with_GTTT"] + sum(found)
      categories["no_downstream_GTTT"] <- categories["no_downstream_GTTT"] + sum(!found)
    }

    valid_guides <- guide_raw[valid_bases]
    ref_idx <- match(valid_guides, ref$guide)
    is_exact <- !is.na(ref_idx)
    if (any(is_exact)) exact_counts <- exact_counts + tabulate(ref_idx[is_exact], nbins = nrow(ref))
    categories["exact_reference"] <- categories["exact_reference"] + sum(is_exact)
    categories["nonreference_20nt"] <- categories["nonreference_20nt"] + sum(!is_exact)

    nonref <- valid_guides[!is_exact]
    if (length(nonref)) {
      neighbor_count <- mismatch_lookup[.(nonref), n_reference_neighbors]
      categories["one_mismatch_unique"] <- categories["one_mismatch_unique"] + sum(neighbor_count == 1L, na.rm = TRUE)
      categories["one_mismatch_ambiguous"] <- categories["one_mismatch_ambiguous"] + sum(neighbor_count > 1L, na.rm = TRUE)
      categories["more_than_one_mismatch_or_other"] <- categories["more_than_one_mismatch_or_other"] + sum(is.na(neighbor_count))
      chunk_top <- data.table(guide = nonref)[, .N, by = guide][order(-N)][seq_len(min(.N, 200L))]
      top_nonref_chunks[[length(top_nonref_chunks) + 1L]] <- chunk_top
    }

    if (chunk_number %% 20L == 0L || n < reads_this_chunk) {
      elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
      rate <- if (elapsed > 0) total_processed / elapsed else NA_real_
      log_msg("Set ", sample_name, ": ", format(total_processed, big.mark = ","), " reads; ",
              format(round(rate), big.mark = ","), " reads/min")
    }
    rm(seqs, guide_raw, valid_guides, nonref, ref_idx)
    if (n < reads_this_chunk) break
  }
  close(con)
  on.exit(NULL, add = FALSE)
  rm(mismatch_lookup)
  gc(verbose = FALSE)

  median_count <- median(exact_counts)
  total_exact <- sum(exact_counts)
  ref[, exact_count := exact_counts]
  ref[, proportion := if (total_exact > 0) exact_count / total_exact else 0]
  ref[, cpm := proportion * 1e6]
  ref[, fold_vs_median := if (median_count > 0) exact_count / median_count else NA_real_]
  ref[, log2_count_plus1 := log2(exact_count + 1)]
  ref[, zero_count := exact_count == 0]
  ref[, below_0.1x_median := if (median_count > 0) exact_count < 0.1 * median_count else NA]
  ref[, sample := sample_name]
  setcolorder(ref, c("sample", "guide", "gene_symbol", "gene_id", "guide_type"))

  q <- quantile(exact_counts, probs = c(0, .01, .05, .10, .25, .50, .75, .90, .95, .99, 1), names = FALSE)
  q10 <- unname(quantile(exact_counts, .10))
  q90 <- unname(quantile(exact_counts, .90))
  summary_dt <- data.table(
    sample = sample_name,
    total_r1_reads = categories["total_reads"],
    exact_mapped_reads = categories["exact_reference"],
    exact_mapping_percent = 100 * categories["exact_reference"] / categories["total_reads"],
    valid_20nt_percent = 100 * categories["valid_20nt_insert"] / categories["total_reads"],
    unique_one_mismatch_percent = 100 * categories["one_mismatch_unique"] / categories["total_reads"],
    reference_guides = nrow(ref),
    non_targeting_guides = sum(ref$guide_type == "Non-targeting control"),
    mean_exact_reads_per_guide = mean(exact_counts),
    median_exact_reads_per_guide = median_count,
    minimum_exact_reads = min(exact_counts),
    maximum_exact_reads = max(exact_counts),
    zero_guides = sum(exact_counts == 0),
    zero_guides_percent = 100 * mean(exact_counts == 0),
    below_0.1x_median_guides = sum(exact_counts < 0.1 * median_count),
    below_0.1x_median_percent = 100 * mean(exact_counts < 0.1 * median_count),
    q90_q10_ratio = if (q10 > 0) q90 / q10 else Inf,
    coefficient_of_variation = sd(exact_counts) / mean(exact_counts),
    gini_raw_counts = gini_coefficient(exact_counts),
    gini_log2_counts = gini_coefficient(log2(exact_counts + 1)),
    runtime_minutes = as.numeric(difftime(Sys.time(), started, units = "mins"))
  )
  quantiles_dt <- data.table(
    sample = sample_name,
    percentile = c(0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 100),
    exact_count = q
  )
  coverage_dt <- data.table(
    sample = sample_name,
    threshold = c(1, 10, 100, 300),
    guides_at_or_above = vapply(c(1, 10, 100, 300), function(z) sum(exact_counts >= z), numeric(1)),
    percent_at_or_above = 100 * vapply(c(1, 10, 100, 300), function(z) mean(exact_counts >= z), numeric(1))
  )
  classification_dt <- data.table(category = names(categories), reads = as.numeric(categories))
  classification_dt[, sample := sample_name]
  classification_dt[, percent_total := 100 * reads / categories["total_reads"]]
  classification_dt <- classification_dt[, .(sample, category, reads, percent_total)]
  insert_dt <- data.table(insert_length = as.integer(names(insert_lengths)), reads = as.numeric(insert_lengths), sample = sample_name)
  stagger_dt <- data.table(CACCG_start_position = as.integer(names(stagger_positions)), reads = as.numeric(stagger_positions), sample = sample_name)
  read_length_dt <- data.table(read_length = as.integer(names(read_lengths)), reads = as.numeric(read_lengths), sample = sample_name)
  top_nonref <- if (length(top_nonref_chunks)) {
    rbindlist(top_nonref_chunks)[, .(reads = sum(N)), by = guide][order(-reads)][seq_len(min(.N, 1000L))]
  } else data.table(guide = character(), reads = numeric())
  top_nonref[, sample := sample_name]

  gene_dt <- ref[gene_symbol != "NO-TARGET", .(
    exact_count = sum(exact_count),
    guides = .N
  ), by = .(sample, gene_symbol, gene_id)]
  gene_dt[, cpm := exact_count / sum(exact_count) * 1e6]

  fwrite(ref, file.path(out_dir, "tables", paste0("guide_counts_set_", sample_name, ".csv")))
  fwrite(gene_dt, file.path(out_dir, "tables", paste0("gene_counts_set_", sample_name, ".csv")))
  fwrite(classification_dt, file.path(out_dir, "diagnostics", paste0("read_classification_set_", sample_name, ".csv")))
  fwrite(insert_dt, file.path(out_dir, "diagnostics", paste0("insert_lengths_set_", sample_name, ".csv")))
  fwrite(stagger_dt, file.path(out_dir, "diagnostics", paste0("stagger_positions_set_", sample_name, ".csv")))
  fwrite(read_length_dt, file.path(out_dir, "diagnostics", paste0("read_lengths_set_", sample_name, ".csv")))
  fwrite(top_nonref, file.path(out_dir, "diagnostics", paste0("top_nonreference_set_", sample_name, ".csv")))
  fwrite(quantiles_dt, file.path(out_dir, "tables", paste0("count_quantiles_set_", sample_name, ".csv")))
  fwrite(coverage_dt, file.path(out_dir, "tables", paste0("coverage_thresholds_set_", sample_name, ".csv")))
  saveRDS(list(summary = summary_dt, guide = ref, gene = gene_dt, classification = classification_dt,
               quantiles = quantiles_dt, coverage = coverage_dt),
          file.path(out_dir, paste0("set_", sample_name, "_results.rds")))
  log_msg("Finished Set ", sample_name, " in ", round(summary_dt$runtime_minutes, 1), " minutes")
  list(summary = summary_dt, guide = ref, gene = gene_dt, classification = classification_dt,
       quantiles = quantiles_dt, coverage = coverage_dt)
}

results <- lapply(seq_len(nrow(sample_config)), function(i) process_sample(sample_config[i], args$max_reads, args$chunk_reads))
names(results) <- sample_config$sample

summary_all <- rbindlist(lapply(results, `[[`, "summary"), fill = TRUE)
guide_all <- rbindlist(lapply(results, `[[`, "guide"), fill = TRUE)
gene_all <- rbindlist(lapply(results, `[[`, "gene"), fill = TRUE)
classification_all <- rbindlist(lapply(results, `[[`, "classification"), fill = TRUE)
coverage_all <- rbindlist(lapply(results, `[[`, "coverage"), fill = TRUE)

qc_flags <- rbindlist(list(
  summary_all[, .(sample, metric = "Exact mapping (%)", value = exact_mapping_percent, threshold = ">= 65", pass = exact_mapping_percent >= 65)],
  summary_all[, .(sample, metric = "Mean reads per guide", value = mean_exact_reads_per_guide, threshold = ">= 300", pass = mean_exact_reads_per_guide >= 300)],
  summary_all[, .(sample, metric = "Zero-count guides (%)", value = zero_guides_percent, threshold = "<= 1", pass = zero_guides_percent <= 1)],
  summary_all[, .(sample, metric = "Gini of log2(count+1)", value = gini_log2_counts, threshold = "<= 0.1", pass = gini_log2_counts <= 0.1)]
))
fwrite(summary_all, file.path(out_dir, "tables", "sample_summary.csv"))
fwrite(classification_all, file.path(out_dir, "diagnostics", "read_classification_all.csv"))
fwrite(coverage_all, file.path(out_dir, "tables", "coverage_thresholds_all.csv"))
fwrite(qc_flags, file.path(out_dir, "tables", "qc_flags.csv"))

if (all(c("A", "B") %in% names(results))) {
  gene_comparison <- merge(
    results$A$gene[, .(gene_symbol, gene_id, cpm_A = cpm, count_A = exact_count)],
    results$B$gene[, .(gene_symbol, gene_id, cpm_B = cpm, count_B = exact_count)],
    by = c("gene_symbol", "gene_id")
  )
  gene_cor_spearman <- cor(log2(gene_comparison$cpm_A + 1), log2(gene_comparison$cpm_B + 1), method = "spearman")
  gene_cor_pearson <- cor(log2(gene_comparison$cpm_A + 1), log2(gene_comparison$cpm_B + 1), method = "pearson")
  fwrite(gene_comparison, file.path(out_dir, "tables", "set_A_vs_B_gene_abundance.csv"))
} else {
  gene_comparison <- NULL
  gene_cor_spearman <- gene_cor_pearson <- NA_real_
}

palette <- c(A = "#0072B2", B = "#D55E00")
theme_qc <- theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11), legend.position = "top")

rank_dt <- copy(guide_all)
rank_dt[, rank_fraction := frank(-exact_count, ties.method = "average") / .N, by = sample]
p_rank <- ggplot(rank_dt, aes(rank_fraction, exact_count + 1, color = sample)) +
  geom_line(linewidth = 0.7) + scale_y_log10(labels = comma) +
  scale_color_manual(values = palette) +
  labs(title = "Guide rank-abundance", x = "Fraction of guides ranked", y = "Exact reads + 1", color = "Set") + theme_qc

p_dist <- ggplot(guide_all, aes(log2_count_plus1, color = sample, fill = sample)) +
  geom_density(alpha = 0.18, linewidth = 0.7) + scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "Guide count distribution", x = "log2(exact reads + 1)", y = "Density", color = "Set", fill = "Set") + theme_qc

lorenz_dt <- guide_all[order(sample, exact_count), {
  total <- sum(exact_count)
  .(guide_fraction = c(0, seq_len(.N) / .N), read_fraction = c(0, if (total > 0) cumsum(exact_count) / total else rep(0, .N)))
}, by = sample]
p_lorenz <- ggplot(lorenz_dt, aes(guide_fraction, read_fraction, color = sample)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey55") +
  geom_line(linewidth = 0.8) + coord_equal() + scale_color_manual(values = palette) +
  labs(title = "Lorenz curve", x = "Cumulative fraction of guides", y = "Cumulative fraction of exact reads", color = "Set") + theme_qc

p_cov <- ggplot(coverage_all, aes(factor(threshold), percent_at_or_above, fill = sample)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  scale_fill_manual(values = palette) + scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, .03))) +
  labs(title = "Guide coverage thresholds", x = "Minimum exact reads", y = "Guides at or above threshold (%)", fill = "Set") + theme_qc

p_type <- ggplot(guide_all, aes(guide_type, log2_count_plus1, fill = sample)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.55, color = NA) +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.85, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = palette) +
  labs(title = "Targeting and control guides", x = NULL, y = "log2(exact reads + 1)", fill = "Set") + theme_qc +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

if (!is.null(gene_comparison)) {
  p_gene <- ggplot(gene_comparison, aes(log2(cpm_A + 1), log2(cpm_B + 1))) +
    geom_point(size = 0.55, alpha = 0.22, color = "#333333") +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#0072B2") +
    labs(title = "Shared-gene abundance", subtitle = sprintf("Spearman rho = %.3f; Pearson r = %.3f", gene_cor_spearman, gene_cor_pearson),
         x = "Set A log2(CPM + 1)", y = "Set B log2(CPM + 1)") + theme_qc
} else {
  p_gene <- ggplot() + annotate("text", x = 0, y = 0, label = "Both sets required for comparison") + theme_void()
}

overview <- (p_rank | p_dist) / (p_lorenz | p_cov) / (p_type | p_gene) +
  plot_annotation(title = "Dolcetto amplified plasmid library representation QC",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))

figure_list <- list(rank_abundance = p_rank, count_distribution = p_dist, lorenz_curve = p_lorenz,
                    coverage_thresholds = p_cov, guide_type_distribution = p_type,
                    gene_abundance_comparison = p_gene, qc_overview = overview)
for (nm in names(figure_list)) {
  width <- if (nm == "qc_overview") 10 else 5.5
  height <- if (nm == "qc_overview") 11 else 4.2
  ggsave(file.path(out_dir, "figures", paste0(nm, ".png")), figure_list[[nm]], width = width, height = height, dpi = 300, bg = "white")
  ggsave(file.path(out_dir, "figures", paste0(nm, ".pdf")), figure_list[[nm]], width = width, height = height, device = cairo_pdf, bg = "white")
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", ifelse(abs(x) >= 1000, comma(round(x)), format(round(x, digits), nsmall = digits, trim = TRUE)))
}
summary_display <- summary_all[, .(
  Set = sample,
  `R1 reads` = comma(total_r1_reads),
  `Exact mapped (%)` = sprintf("%.2f", exact_mapping_percent),
  `Mean reads/guide` = comma(round(mean_exact_reads_per_guide)),
  `Median reads/guide` = comma(round(median_exact_reads_per_guide)),
  `Zero guides (%)` = sprintf("%d (%.3f%%)", zero_guides, zero_guides_percent),
  `90th/10th ratio` = sprintf("%.2f", q90_q10_ratio),
  `Raw Gini` = sprintf("%.3f", gini_raw_counts),
  `Log-count Gini` = sprintf("%.3f", gini_log2_counts)
)]
flags_display <- copy(qc_flags)
flags_display[, value := ifelse(metric %in% c("Exact mapping (%)", "Zero-count guides (%)"), sprintf("%.2f", value), sprintf("%.3f", value))]
flags_display[, pass := ifelse(pass, "PASS", "FLAG")]
setnames(flags_display, c("sample", "metric", "value", "threshold", "pass"), c("Set", "Metric", "Observed", "Default threshold", "Status"))

html_table <- function(dt) HTML(knitr::kable(dt, format = "html", escape = TRUE, table.attr = 'class="qc-table"'))
report <- tagList(
  tags$html(
    tags$head(
      tags$title("Dolcetto library representation QC"),
      tags$style(HTML("body{font-family:Arial,Helvetica,sans-serif;max-width:1100px;margin:35px auto;padding:0 24px;color:#222;line-height:1.45} h1,h2{color:#15324a} .note{background:#eef6fb;border-left:5px solid #0072B2;padding:12px 16px}.qc-table{border-collapse:collapse;width:100%;margin:12px 0 24px}.qc-table th,.qc-table td{border:1px solid #ddd;padding:7px;text-align:left}.qc-table th{background:#f2f2f2}.figure{max-width:100%;height:auto;border:1px solid #eee}.small{font-size:0.9em;color:#555} code{background:#f5f5f5;padding:2px 4px}"))
    ),
    tags$body(
      tags$h1("Dolcetto amplified plasmid library representation QC"),
      tags$p(class = "small", paste("Generated", format(Sys.time(), "%Y-%m-%d %H:%M %Z"))),
      tags$div(class = "note", tags$b("Counting policy: "),
               "R1 guides were extracted between CACCG and GTTT. Exact 20-nt matches are the authoritative counts. One-substitution matches are diagnostic only and were not added to guide counts."),
      tags$h2("Executive summary"), html_table(summary_display),
      tags$h2("Default QC flags"),
      tags$p("Thresholds are used as transparent flags rather than a single hard pass/fail verdict."), html_table(flags_display),
      tags$h2("Representation overview"), tags$img(src = "figures/qc_overview.png", class = "figure"),
      tags$h2("Read classification"), html_table(classification_all[, .(Set = sample, Category = category, Reads = comma(reads), `Percent of R1` = sprintf("%.3f", percent_total))]),
      tags$h2("Interpretation notes"),
      tags$ul(
        tags$li("Set A and Set B are assessed independently because they contain distinct guide sequences."),
        tags$li("R2 is backbone-dominated and is not used for guide counting."),
        tags$li("The local reference files contain 496 rows annotated NO-TARGET in each set; this is retained without modification."),
        tags$li("Low exact mapping can reflect sequencing substitutions, oligo synthesis errors, cloning variants, or malformed inserts. The diagnostic one-mismatch categories help distinguish these possibilities.")
      ),
      if (!is.null(gene_comparison)) tagList(tags$h2("Set A versus Set B at gene level"),
        tags$p(sprintf("Across %s shared annotated genes, Spearman rho = %.3f and Pearson r = %.3f for log2(CPM + 1).", comma(nrow(gene_comparison)), gene_cor_spearman, gene_cor_pearson)),
        tags$img(src = "figures/gene_abundance_comparison.png", class = "figure")) else NULL,
      tags$h2("Methods and references"),
      tags$p("The workflow streams compressed FASTQ records in chunks, locates the stagger-tolerant CACCG prefix, extracts the following 20 nucleotides, requires the downstream GTTT sequence, and matches exactly against the corresponding Dolcetto reference."),
      tags$ul(
        tags$li(tags$a(href = "https://www.addgene.org/pooled-library/broadgpp-human-crispri-dolcetto/", "Addgene Dolcetto library specification")),
        tags$li(tags$a(href = "https://media.addgene.org/cms/filer_public/61/16/611619f4-0926-4a07-b5c7-e286a8ecf7f5/broadgpp-sequencing-protocol.pdf", "Broad GPP sgRNA sequencing protocol")),
        tags$li(tags$a(href = "https://pmc.ncbi.nlm.nih.gov/articles/PMC4699372/", "MAGeCK-VISPR QC framework"))
      ),
      tags$h2("Output files"),
      tags$p("Guide-level and gene-level counts are in tables/. Read-structure and mismatch diagnostics are in diagnostics/. PNG and vector PDF figures are in figures/."),
      tags$h2("R session"), tags$pre(paste(capture.output(sessionInfo()), collapse = "\n"))
    )
  )
)
save_html(
  report,
  file = file.path(out_dir, "Dolcetto_library_representation_QC.html"),
  libdir = file.path(out_dir, "report_files"),
  background = "white"
)
saveRDS(list(summary = summary_all, qc_flags = qc_flags, gene_cor_spearman = gene_cor_spearman,
             gene_cor_pearson = gene_cor_pearson, arguments = args), file.path(out_dir, "combined_results.rds"))

log_msg("Analysis complete. Report: ", file.path(out_dir, "Dolcetto_library_representation_QC.html"))
