#!/usr/bin/env Rscript
# Candidate Fig. 2: genome-wide CRC-NSCLC microbial-KO concordance.
#
# Question: do KO shifts toward the favorable state agree between
# CRC (Healthy - Cancer) and NSCLC under ICI (Responder - Non-responder)?
#
# Design:
# - Primary unit: biological `Sample name`; repeated HGMT Run IDs are averaged.
# - Sensitivity unit: Run ID, matching the historical project pipelines.
# - Within-cohort prevalence gate: >=10%.
# - Effect: rank-biserial correlation, favorable minus adverse; ties contribute 0.
# - Variance: DeLong AUC placement variance transformed by r_rb = 2*AUC - 1.
# - Disease pooling: DerSimonian-Laird random effects, estimable in >=2 cohorts.
# - Multiplicity: BH FDR separately within CRC and NSCLC.
# - Concordance: Spearman rho and 10,000 joint-prevalence-stratified permutations.
# - Sarcosine KOs are pre-defined annotations, not outcome-selected.

set.seed(42)
options(stringsAsFactors = FALSE)
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(ggplot2)
  library(ggrepel)
  library(here)
  library(jsonlite)
  library(ragg)
})

ROOT <- normalizePath(here::here(), mustWork = TRUE)
OMICS <- file.path(ROOT, "공공_Metabolomics&Metagenomics_분석모음")
NS_ROOT <- file.path(OMICS, "HGMT_NSCLC_ICI_RvsNR_WGS")
CRC_ROOT <- file.path(OMICS, "HGMT_CRC_WGS-Healthy_vs_Cancer")
THEME_FILE <- file.path(OMICS, "_shared", "theme_nc_26.08.18.R")
OUT <- file.path(
  NS_ROOT, "pooled_analysis", "results",
  "Fig2_genomewide_KO_cross_disease_concordance_CANDIDATE_26.08.23"
)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(NS_ROOT), dir.exists(CRC_ROOT), file.exists(THEME_FILE))
source(THEME_FILE)

PREV <- 0.10
MIN_K <- 2L
N_PERM <- 10000L
SEED <- 42L

log_con <- file(file.path(OUT, "run_log.txt"), open = "wt")
on.exit(try(close(log_con), silent = TRUE), add = TRUE)
say <- function(...) {
  x <- sprintf(...)
  cat(x, "\n")
  cat(x, "\n", file = log_con)
  flush(log_con)
}

defs <- data.table(
  disease = c(rep("CRC", 4), rep("NSCLC", 3)),
  cohort = c(
    "PRJEB6070", "PRJEB10878", "PRJEB27928", "PRJNA429097",
    "PRJNA751792", "PRJNA1023797", "PRJEB22863"
  ),
  dir = c(
    file.path(CRC_ROOT, "PRJEB6070_CRC_AdenomatousPolyps"),
    file.path(CRC_ROOT, "PRJEB10878_CRC"),
    file.path(CRC_ROOT, "PRJEB27928_CRC"),
    file.path(CRC_ROOT, "PRJNA429097_CRC"),
    file.path(NS_ROOT, "NSCLC_PRJNA751792"),
    file.path(NS_ROOT, "NSCLC_PRJNA1023797"),
    file.path(NS_ROOT, "NSCLC_RCC_PRJEB22863")
  )
)
defs[, ko_file := file.path(dir, "KO_relative_abundance.tsv")]
stopifnot(all(dir.exists(defs$dir)), all(file.exists(defs$ko_file)))

read_meta <- function(path) {
  fs <- sort(list.files(path, "^selected_project_.*\\.txt$", full.names = TRUE))
  if (!length(fs)) stop("No selected_project file: ", path)
  if (length(fs) > 1L) {
    h <- vapply(fs, digest, character(1), file = TRUE, algo = "sha256")
    if (length(unique(h)) != 1L) stop("Non-identical metadata files: ", path)
  }
  z <- sub("\r$", "", readLines(fs[1], warn = FALSE))
  hdr <- strsplit(z[2], "\t", fixed = TRUE)[[1]]
  body <- z[-c(1, 2)]
  body <- body[nzchar(body)]
  rows <- lapply(strsplit(body, "\t", fixed = TRUE), function(x) {
    length(x) <- length(hdr)
    x
  })
  out <- as.data.table(do.call(rbind, rows))
  setnames(out, hdr)
  out
}

get_response <- function(x) {
  vapply(as.character(x), function(s) {
    m <- regmatches(s, regexpr("response_group:[^;]*", s))
    if (!length(m) || !nzchar(m)) m <- regmatches(s, regexpr("response:[^;]*", s))
    if (!length(m) || !nzchar(m)) return(NA_character_)
    y <- trimws(sub("^response(_group)?:", "", m))
    if (y %in% c("R", "NR")) y else NA_character_
  }, character(1))
}

prepare_meta <- function(d) {
  m <- read_meta(d$dir)
  if (d$disease == "CRC") {
    m <- m[`Assay type` == "WGS" & `Phenotype name` %in% c("Health", "Colorectal Neoplasms")]
    m[, group := fifelse(`Phenotype name` == "Health", "Favorable", "Adverse")]
  } else {
    m <- m[`Assay type` == "WGS" & `Phenotype name` == "Carcinoma, Non-Small-Cell Lung"]
    m[, response := get_response(`Sample description`)]
    m <- m[response %in% c("R", "NR")]
    m[, group := fifelse(response == "R", "Favorable", "Adverse")]
  }
  m <- m[, .(RunID = `Run ID`, SampleID = `Sample name`, group)]
  if (anyNA(m) || any(!nzchar(m$RunID)) || any(!nzchar(m$SampleID))) stop(d$cohort, ": bad metadata")
  if (anyDuplicated(m$RunID)) stop(d$cohort, ": duplicated RunID")
  if (m[, uniqueN(group), by = SampleID][V1 > 1L, .N]) stop(d$cohort, ": SampleID crosses groups")
  m
}

# DeLong placements handle ties explicitly and avoid a normality assumption.
delta_delong <- function(fav, adv) {
  n1 <- length(fav); n0 <- length(adv)
  ys <- sort(adv)
  yle <- findInterval(fav, ys)
  ylt <- findInterval(fav, ys, left.open = TRUE)
  v10 <- (ylt + 0.5 * (yle - ylt)) / n0
  xs <- sort(fav)
  xle <- findInterval(adv, xs)
  xlt <- findInterval(adv, xs, left.open = TRUE)
  v01 <- ((n1 - xle) + 0.5 * (xle - xlt)) / n1
  if (!isTRUE(all.equal(mean(v10), mean(v01), tolerance = 1e-12))) stop("DeLong placement mismatch")
  c(effect = 2 * mean(v10) - 1, variance = 4 * (var(v10) / n1 + var(v01) / n0))
}

matrix_effects <- function(mat, group, d, mode) {
  stopifnot(nrow(mat) == length(group), all(is.finite(mat)), all(mat >= 0))
  pprev <- colMeans(mat > 0)
  mat <- mat[, pprev >= PREV, drop = FALSE]
  pprev <- pprev[pprev >= PREV]
  fav <- group == "Favorable"; adv <- group == "Adverse"
  if (sum(fav) < 3L || sum(adv) < 3L) stop(d$cohort, ": group n < 3")

  ans <- rbindlist(lapply(seq_len(ncol(mat)), function(j) {
    x <- mat[fav, j]; y <- mat[adv, j]
    z <- delta_delong(x, y)
    w <- suppressWarnings(wilcox.test(x, y, exact = FALSE, correct = FALSE))
    data.table(
      KO = colnames(mat)[j], prevalence = unname(pprev[j]),
      mean_favorable = mean(x), mean_adverse = mean(y),
      effect = unname(z["effect"]), variance = unname(z["variance"]),
      wilcox_p = w$p.value
    )
  }))
  ans[, `:=`(
    disease = d$disease, cohort = d$cohort, mode = mode,
    n_favorable = sum(fav), n_adverse = sum(adv)
  )]
  setcolorder(ans, c(
    "disease", "cohort", "mode", "KO", "n_favorable", "n_adverse",
    "prevalence", "mean_favorable", "mean_adverse", "effect", "variance", "wilcox_p"
  ))

  # Independent checks: direct pairwise definition and Wilcoxon W.
  idx <- unique(round(seq(1, ncol(mat), length.out = min(7L, ncol(mat)))))
  for (j in idx) {
    x <- mat[fav, j]; y <- mat[adv, j]
    direct <- mean(outer(x, y, function(a, b) sign(a - b)))
    got <- ans[KO == colnames(mat)[j], effect]
    if (length(got) != 1L || !is.finite(got) || abs(got - direct) > 1e-12) {
      stop(sprintf(
        "%s: direct r_rb check failed for %s (computed=%.17g, direct=%.17g, n_matches=%d)",
        d$cohort, colnames(mat)[j], got[1], direct, length(got)
      ))
    }
    W <- suppressWarnings(wilcox.test(x, y, exact = FALSE, correct = FALSE))$statistic
    from_w <- 2 * unname(W) / (length(x) * length(y)) - 1
    if (length(got) != 1L || !is.finite(got) || abs(got - from_w) > 1e-12) {
      stop(d$cohort, ": Wilcoxon-W r_rb check failed")
    }
  }
  ans
}

analyse_cohort <- function(d) {
  m0 <- prepare_meta(d)
  say("[%s] metadata: runs=%d; samples=%d; duplicate-run rows=%d", d$cohort, nrow(m0), uniqueN(m0$SampleID), nrow(m0) - uniqueN(m0$SampleID))
  say("[%s] reading KO JSON: %.1f MB", d$cohort, file.info(d$ko_file)$size / 1024^2)
  ko <- as.data.table(fromJSON(d$ko_file, simplifyDataFrame = TRUE))
  stopifnot(all(c("ko", "run_id", "abundance") %in% names(ko)))
  ko <- ko[, .(KO = as.character(ko), RunID = as.character(run_id), abundance = as.numeric(abundance))]
  if (anyNA(ko) || any(!is.finite(ko$abundance)) || any(ko$abundance < 0)) stop(d$cohort, ": bad KO values")
  if (ko[, anyDuplicated(paste(RunID, KO, sep = "\r"))]) stop(d$cohort, ": duplicate RunID-KO")

  ko_runs <- unique(ko$RunID)
  missing_ko <- setdiff(m0$RunID, ko_runs)
  outside <- setdiff(ko_runs, m0$RunID)
  m <- m0[RunID %in% ko_runs]
  ko <- ko[RunID %in% m$RunID]
  say("[%s] matched KO runs=%d; metadata without KO=%d; KO runs outside filter=%d", d$cohort, nrow(m), length(missing_ko), length(outside))

  # Run-level sensitivity.
  rl <- merge(ko, m[, .(RunID, group)], by = "RunID", sort = FALSE)
  rw <- dcast(rl, RunID + group ~ KO, value.var = "abundance", fill = 0)
  rg <- rw$group
  run_mat <- as.matrix(rw[, -c("RunID", "group")]); storage.mode(run_mat) <- "double"
  run_eff <- matrix_effects(run_mat, rg, d, "run_level_sensitivity")
  rm(rl, rw, run_mat); invisible(gc(FALSE))

  # Sample-level primary: absent KO-run rows contribute zero to the run average.
  smap <- m[, .(n_runs = .N, group = unique(group)), by = SampleID]
  sl <- merge(ko, m[, .(RunID, SampleID)], by = "RunID", sort = FALSE)
  sl <- sl[, .(sum_abundance = sum(abundance)), by = .(SampleID, KO)]
  sl <- merge(sl, smap, by = "SampleID")
  sl[, abundance := sum_abundance / n_runs]
  sw <- dcast(sl, SampleID + group ~ KO, value.var = "abundance", fill = 0)
  sg <- sw$group
  sample_mat <- as.matrix(sw[, -c("SampleID", "group")]); storage.mode(sample_mat) <- "double"
  sample_eff <- matrix_effects(sample_mat, sg, d, "biological_sample_primary")

  qc <- data.table(
    disease = d$disease, cohort = d$cohort,
    metadata_runs = nrow(m0), matched_KO_runs = nrow(m), biological_samples = uniqueN(m$SampleID),
    duplicate_run_rows = nrow(m) - uniqueN(m$SampleID),
    favorable_runs = sum(m$group == "Favorable"), adverse_runs = sum(m$group == "Adverse"),
    favorable_samples = sum(smap$group == "Favorable"), adverse_samples = sum(smap$group == "Adverse"),
    KOs_run_prev10 = nrow(run_eff), KOs_sample_prev10 = nrow(sample_eff),
    metadata_without_KO = length(missing_ko), KO_runs_outside_filter = length(outside)
  )
  rm(ko, sl, sw, sample_mat); invisible(gc(FALSE))
  list(effects = rbind(sample_eff, run_eff), qc = qc)
}

say("Started: %s", format(Sys.time(), tz = "Asia/Seoul", usetz = TRUE))

# Input hashes document the exact HGMT profiles and metadata.
meta_files <- vapply(defs$dir, function(x) sort(list.files(x, "^selected_project_.*\\.txt$", full.names = TRUE))[1], character(1))
inputs <- c(defs$ko_file, meta_files, file.path(NS_ROOT, "sarcosine_KO_set.csv"), THEME_FILE)
input_hash <- data.table(
  file = vapply(normalizePath(inputs), function(x) substring(x, nchar(ROOT) + 2L), character(1)),
  bytes = file.info(inputs)$size,
  sha256 = vapply(inputs, digest, character(1), file = TRUE, algo = "sha256")
)
fwrite(input_hash, file.path(OUT, "input_sha256.tsv"), sep = "\t")

eff_list <- vector("list", nrow(defs)); qc_list <- vector("list", nrow(defs))
for (i in seq_len(nrow(defs))) {
  z <- analyse_cohort(defs[i])
  eff_list[[i]] <- z$effects; qc_list[[i]] <- z$qc
  fwrite(rbindlist(eff_list[seq_len(i)]), file.path(OUT, "per_cohort_KO_effects_checkpoint.csv"))
  fwrite(rbindlist(qc_list[seq_len(i)]), file.path(OUT, "cohort_sample_QC.csv"))
}
per <- rbindlist(eff_list); qc <- rbindlist(qc_list)
fwrite(per, file.path(OUT, "per_cohort_KO_effects.csv"))
fwrite(qc, file.path(OUT, "cohort_sample_QC.csv"))

dl_meta <- function(x) {
  x <- x[is.finite(effect) & is.finite(variance) & variance > 0]
  if (nrow(x) < MIN_K) return(NULL)
  w <- 1 / x$variance
  mu0 <- sum(w * x$effect) / sum(w)
  Q <- sum(w * (x$effect - mu0)^2)
  C <- sum(w) - sum(w^2) / sum(w)
  tau2 <- if (C > 0) max(0, (Q - (nrow(x) - 1)) / C) else 0
  wr <- 1 / (x$variance + tau2)
  mu <- sum(wr * x$effect) / sum(wr)
  se <- sqrt(1 / sum(wr))
  data.table(
    k = nrow(x), effect = mu, se = se, ci_low = mu - 1.96 * se, ci_high = mu + 1.96 * se,
    p = 2 * pnorm(abs(mu / se), lower.tail = FALSE), tau2 = tau2, Q = Q,
    I2 = if (Q > 0) max(0, (Q - (nrow(x) - 1)) / Q) * 100 else 0,
    mean_prevalence = mean(x$prevalence), cohorts = paste(sort(x$cohort), collapse = ";")
  )
}

make_meta <- function(disease_name, mode_name, exclude = NA_character_) {
  x <- per[disease == disease_name & mode == mode_name]
  if (!is.na(exclude)) x <- x[cohort != exclude]
  out <- rbindlist(lapply(split(x, x$KO), function(y) {
    z <- dl_meta(y)
    if (is.null(z)) return(NULL)
    z[, KO := y$KO[1]]
    z
  }), fill = TRUE)
  if (!nrow(out)) return(out)
  out[, q := p.adjust(p, "BH")]
  out[, `:=`(disease = disease_name, mode = mode_name, excluded_cohort = exclude)]
  setcolorder(out, c("disease", "mode", "excluded_cohort", "KO", "k", "effect", "se", "ci_low", "ci_high", "p", "q", "tau2", "Q", "I2", "mean_prevalence", "cohorts"))
  out
}

modes <- c("biological_sample_primary", "run_level_sensitivity")
meta <- rbindlist(lapply(modes, function(m) rbind(make_meta("CRC", m), make_meta("NSCLC", m))), fill = TRUE)
fwrite(meta, file.path(OUT, "disease_level_random_effects_KO_meta.csv"))

make_cross <- function(meta_dt, mode_name) {
  a <- meta_dt[disease == "CRC" & mode == mode_name]
  b <- meta_dt[disease == "NSCLC" & mode == mode_name]
  merge(
    a[, .(KO, crc_effect = effect, crc_se = se, crc_p = p, crc_q = q, crc_I2 = I2, crc_k = k, crc_prevalence = mean_prevalence)],
    b[, .(KO, nsclc_effect = effect, nsclc_se = se, nsclc_p = p, nsclc_q = q, nsclc_I2 = I2, nsclc_k = k, nsclc_prevalence = mean_prevalence)],
    by = "KO"
  )
}

prev_bin <- function(x, n = 5L) pmin(n, pmax(1L, ceiling(rank(x, ties.method = "average") / length(x) * n)))
concordance <- function(x, label, do_perm = TRUE, prevalence_bins = 5L) {
  rho <- cor(x$crc_effect, x$nsclc_effect, method = "spearman")
  pa <- suppressWarnings(cor.test(x$crc_effect, x$nsclc_effect, method = "spearman", exact = FALSE)$p.value)
  use <- x$crc_effect != 0 & x$nsclc_effect != 0
  agree <- sum(sign(x$crc_effect[use]) == sign(x$nsclc_effect[use]))
  n_dir <- sum(use)
  pb <- binom.test(agree, n_dir, p = 0.5)$p.value
  pp <- NA_real_
  if (do_perm) {
    set.seed(SEED)
    rx <- rank(x$crc_effect); ry <- rank(x$nsclc_effect)
    strata <- interaction(
      prev_bin(x$crc_prevalence, prevalence_bins),
      prev_bin(x$nsclc_prevalence, prevalence_bins),
      drop = TRUE
    )
    ids <- split(seq_along(ry), strata)
    null <- replicate(N_PERM, {
      ord <- seq_along(ry)
      for (ii in ids) if (length(ii) > 1L) ord[ii] <- sample(ii)
      cor(rx, ry[ord])
    })
    pp <- (1 + sum(abs(null) >= abs(rho))) / (N_PERM + 1)
    fwrite(data.table(mode = label, permutation = seq_along(null), rho = null), file.path(OUT, paste0("permutation_null_", label, ".csv")))
  }
  data.table(
    mode = label, n_common_KOs = nrow(x), prevalence_bins = prevalence_bins,
    spearman_rho = rho, analytic_p = pa,
    permutation_p = pp, permutation_reps = if (do_perm) N_PERM else 0L,
    direction_agree_n = agree, direction_test_n = n_dir,
    direction_agreement_percent = 100 * agree / n_dir, direction_binomial_p = pb
  )
}

cross_primary <- make_cross(meta, "biological_sample_primary")
cross_run <- make_cross(meta, "run_level_sensitivity")
fwrite(cross_primary, file.path(OUT, "cross_disease_common_KOs_primary.csv"))
fwrite(cross_run, file.path(OUT, "cross_disease_common_KOs_runlevel_sensitivity.csv"))
stats <- rbind(concordance(cross_primary, "biological_sample_primary"), concordance(cross_run, "run_level_sensitivity"))
fwrite(stats, file.path(OUT, "global_concordance_statistics.csv"))
prevalence_sensitivity <- concordance(
  cross_primary, "biological_sample_primary_decile_prevalence", TRUE, prevalence_bins = 10L
)
fwrite(prevalence_sensitivity, file.path(OUT, "prevalence_stratification_sensitivity.csv"))

# Leave-one-cohort-out sensitivity for the primary biological-sample analysis.
loco_def <- rbind(
  data.table(disease = "CRC", exclude = defs[disease == "CRC", cohort]),
  data.table(disease = "NSCLC", exclude = defs[disease == "NSCLC", cohort])
)
loco <- rbindlist(lapply(seq_len(nrow(loco_def)), function(i) {
  d <- loco_def$disease[i]; ex <- loco_def$exclude[i]
  other <- if (d == "CRC") "NSCLC" else "CRC"
  mm <- rbind(make_meta(d, "biological_sample_primary", ex), meta[disease == other & mode == "biological_sample_primary"], fill = TRUE)
  z <- concordance(make_cross(mm, "biological_sample_primary"), paste0("LOCO_", ex), FALSE)
  z[, `:=`(excluded_disease = d, excluded_cohort = ex)]
  z
}))
fwrite(loco, file.path(OUT, "leave_one_cohort_out_concordance.csv"))

# Verified display names used by the existing individual-gene Fig. 2 panel.
sarc <- data.table(
  KO = c("K00301", "K00302", "K00303", "K00304", "K00305", "K00306", "K00314", "K18897", "K00315", "K08687", "K08688", "K21833", "K21834"),
  display_gene = c("sarcosine oxidase", "soxA", "soxB", "soxD", "soxG", "PIPOX", "sarcosine dehydrogenase", "sdmt", "DMGDH", "N-carbamoylsarcosine amidase", "creatinase", "dgcA/ddhC", "dgcB"),
  role = c(rep("Degradation", 8), rep("Production", 3), rep("Bidirectional", 2))
)
sarc_set <- fread(file.path(NS_ROOT, "sarcosine_KO_set.csv"))
stopifnot(all(sarc$KO %in% sarc_set$KO), !anyDuplicated(sarc$KO))
sarc_summary <- merge(sarc, cross_primary, by = "KO", all.x = TRUE)
sarc_summary <- merge(sarc_summary, cross_run[, .(KO, crc_effect_run = crc_effect, nsclc_effect_run = nsclc_effect)], by = "KO", all.x = TRUE)
sarc_summary[, primary_meta_analysable := !is.na(crc_effect) & !is.na(nsclc_effect)]
fwrite(sarc_summary, file.path(OUT, "sarcosine_KO_cross_disease_summary.csv"))

# Main-panel candidate plot.
pd <- merge(cross_primary, sarc, by = "KO", all.x = TRUE)
pd[, highlight := !is.na(role)]
cols <- c(Degradation = "#2E5F8A", Production = "#C47B3B", Bidirectional = "#7A5195")
lim <- max(0.25, ceiling(max(abs(c(pd$crc_effect, pd$nsclc_effect))) * 20) / 20)
st <- stats[mode == "biological_sample_primary"]
ptext <- if (st$permutation_p < 0.001) "< 0.001" else paste0("= ", signif(st$permutation_p, 2))
lab <- sprintf("Spearman ρ = %.2f\nStratified permutation P %s\nDirection agreement = %.1f%% (%d/%d)", st$spearman_rho, ptext, st$direction_agreement_percent, st$direction_agree_n, st$direction_test_n)

p <- ggplot(pd, aes(crc_effect, nsclc_effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = pt_lw(0.75), color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = pt_lw(0.75), color = "grey70") +
  geom_abline(slope = 1, linetype = "dotted", linewidth = pt_lw(0.6), color = "grey82") +
  geom_point(data = pd[highlight == FALSE], color = "grey55", alpha = 0.28, size = 0.6, stroke = 0) +
  geom_point(data = pd[highlight == TRUE], aes(color = role), size = 2.0) +
  geom_text_repel(
    data = pd[highlight == TRUE], aes(label = paste0(display_gene, " (", KO, ")")),
    family = "Arial", size = mm_text(5.2), color = "black", box.padding = 0.25,
    point.padding = 0.18, min.segment.length = 0, segment.color = "grey45",
    segment.size = pt_lw(0.55), max.overlaps = Inf, seed = 42, show.legend = FALSE
  ) +
  annotate("text", x = -0.97 * lim, y = 0.97 * lim, label = lab, hjust = 0, vjust = 1, family = "Arial", size = mm_text(5.4), lineheight = 0.95) +
  annotate("text", x = -0.94 * lim, y = -0.94 * lim, label = "←  Cancer higher", hjust = 0, vjust = 0, family = "Arial", size = mm_text(5.0), color = "#C47B3B") +
  annotate("text", x = 0.94 * lim, y = -0.94 * lim, label = "Healthy higher  →", hjust = 1, vjust = 0, family = "Arial", size = mm_text(5.0), color = "#2E5F8A") +
  annotate("text", x = 0.025 * lim, y = 0.94 * lim, label = "R higher  ↑", hjust = 0, vjust = 1, family = "Arial", size = mm_text(5.0), color = "#2E5F8A") +
  annotate("text", x = 0.025 * lim, y = -0.94 * lim, label = "NR higher  ↓", hjust = 0, vjust = 0, family = "Arial", size = mm_text(5.0), color = "#C47B3B") +
  scale_color_manual(values = cols, breaks = names(cols), drop = FALSE) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), clip = "off") +
  labs(
    title = "Cross-disease concordance of microbial KO effects",
    x = expression("CRC effect"~("Healthy"-"Cancer,"~meta~r[rb])),
    y = expression("NSCLC effect"~("R"-"NR,"~meta~r[rb])), color = NULL,
    caption = "R, responder; NR, non-responder"
  ) +
  theme_nc(6.0) +
  theme(
    plot.title = element_text(family = "Arial", face = "bold", size = 7.2, hjust = 0, margin = margin(b = 2)),
    legend.position = "top", legend.justification = "left", legend.box.just = "left", aspect.ratio = 1,
    plot.caption = element_text(family = "Arial", size = 5.2, hjust = 0, margin = margin(t = 2))
  )

safe_save <- function(plot, out, type) {
  tmp <- tempfile(fileext = paste0(".", type))
  if (type == "png") ggsave(tmp, plot, width = 4.0, height = 3.65, dpi = 600, bg = "white", device = ragg::agg_png)
  if (type == "pdf") {
    if (!capabilities("aqua")) stop("Quartz PDF output requires macOS Aqua support")
    quartz(file = tmp, type = "pdf", width = 4.0, height = 3.65)
    print(plot)
    grDevices::dev.off()
  }
  stopifnot(file.copy(tmp, out, overwrite = TRUE)); unlink(tmp)
  stopifnot(file.exists(out), file.info(out)$size > 5000)
}
png_file <- file.path(OUT, "Fig2_candidate_genomewide_KO_cross_disease_concordance.png")
pdf_file <- file.path(OUT, "Fig2_candidate_genomewide_KO_cross_disease_concordance.pdf")
safe_save(p, png_file, "png"); safe_save(p, pdf_file, "pdf")

# macOS PowerPoint drag-and-drop compatibility derivative: identical raster pixels,
# explicit sRGB profile and a short ASCII filename. The Desktop convenience copy is
# intentionally not created by this reproducible project script.
ppt_png_file <- file.path(OUT, "PPT_INSERT_Fig2_KO_concordance.png")
sips <- Sys.which("sips")
srgb_profile <- "/System/Library/ColorSync/Profiles/sRGB Profile.icc"
if (!nzchar(sips) || !file.exists(srgb_profile)) stop("macOS sips/sRGB profile unavailable")
sips_out <- system2(
  sips,
  c("-m", shQuote(srgb_profile), shQuote(png_file), "--out", shQuote(ppt_png_file)),
  stdout = TRUE, stderr = TRUE
)
sips_status <- attr(sips_out, "status"); if (is.null(sips_status)) sips_status <- 0L
if (sips_status != 0L || !file.exists(ppt_png_file)) stop("PowerPoint PNG conversion failed")

# Primary-versus-run unit sensitivity.
cm <- merge(cross_primary[, .(KO, crc_primary = crc_effect, ns_primary = nsclc_effect)], cross_run[, .(KO, crc_run = crc_effect, ns_run = nsclc_effect)], by = "KO")
mode_qc <- data.table(
  disease = c("CRC", "NSCLC"),
  spearman_primary_vs_run = c(cor(cm$crc_primary, cm$crc_run, method = "spearman"), cor(cm$ns_primary, cm$ns_run, method = "spearman")),
  max_abs_effect_difference = c(max(abs(cm$crc_primary - cm$crc_run)), max(abs(cm$ns_primary - cm$ns_run)))
)
fwrite(mode_qc, file.path(OUT, "primary_vs_runlevel_QC_statistics.csv"))

# Check the hand-coded DL estimator against metafor on a deterministic subset.
if (requireNamespace("metafor", quietly = TRUE)) {
  kos <- sort(unique(per[mode == "biological_sample_primary", KO]))
  kos <- kos[unique(round(seq(1, length(kos), length.out = min(12L, length(kos)))))]
  max_diff <- 0
  for (dd in c("CRC", "NSCLC")) for (kk in kos) {
    x <- per[disease == dd & mode == "biological_sample_primary" & KO == kk & is.finite(variance) & variance > 0]
    if (nrow(x) < MIN_K) next
    ref <- metafor::rma.uni(yi = x$effect, vi = x$variance, method = "DL")
    got <- meta[disease == dd & mode == "biological_sample_primary" & KO == kk, effect]
    max_diff <- max(max_diff, abs(as.numeric(ref$b) - got))
  }
  if (max_diff > 1e-9) stop("DL estimator differs from metafor")
  fwrite(data.table(max_abs_effect_difference_vs_metafor = max_diff), file.path(OUT, "DL_meta_implementation_check.csv"))
}

plot_hash <- data.table(
  file = basename(c(png_file, pdf_file, ppt_png_file)),
  bytes = file.info(c(png_file, pdf_file, ppt_png_file))$size,
  sha256 = vapply(c(png_file, pdf_file, ppt_png_file), digest, character(1), file = TRUE, algo = "sha256")
)
fwrite(plot_hash, file.path(OUT, "output_plot_sha256.tsv"), sep = "\t")

readme <- c(
  "# Candidate Figure 2: genome-wide KO cross-disease concordance",
  "",
  "Primary analysis averages technical Run IDs within biological Sample name before testing.",
  "CRC effects are Healthy minus Cancer; NSCLC effects are Responder minus Non-responder.",
  "Per-cohort effects are rank-biserial correlations, pooled by DL random effects.",
  "KOs require >=10% prevalence and estimability in >=2 cohorts per disease.",
  "BH FDR is calculated separately within the CRC and NSCLC genome-wide KO meta-analyses.",
  "Cross-disease inference uses all common KOs and 10,000 joint-prevalence-stratified permutations.",
  "The primary five-bin result is repeated with ten prevalence bins as a sensitivity analysis.",
  "Use PPT_INSERT_Fig2_KO_concordance.png for macOS Finder-to-PowerPoint drag-and-drop.",
  "That file has identical raster pixels, an embedded sRGB profile and a short ASCII name.",
  "",
  "The run-level analysis is sensitivity only and mirrors the historical pipeline unit.",
  "Leave-one-cohort-out results test whether one cohort drives the concordance.",
  "",
  "Interpretation: phenotype-associated functional-potential concordance only; no flux,",
  "taxonomic attribution, causation, or independent metabolomics validation.",
  "Any promotion must foreground the robust genome-wide result and transparently retain",
  "the mixed directions and incomplete analysability of the pre-defined sarcosine KOs."
)
writeLines(readme, file.path(OUT, "README.md"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))

say("Primary: nKO=%d rho=%.6f permutation_p=%.6g direction_agreement=%.2f%%", st$n_common_KOs, st$spearman_rho, st$permutation_p, st$direction_agreement_percent)
say("Run-level: rho=%.6f permutation_p=%.6g", stats[mode == "run_level_sensitivity", spearman_rho], stats[mode == "run_level_sensitivity", permutation_p])
unlink(file.path(OUT, "per_cohort_KO_effects_checkpoint.csv"))
say("PNG: %s", substring(png_file, nchar(ROOT) + 2L))
say("PPT PNG: %s", substring(ppt_png_file, nchar(ROOT) + 2L))
say("Completed: %s", format(Sys.time(), tz = "Asia/Seoul", usetz = TRUE))
close(log_con)
