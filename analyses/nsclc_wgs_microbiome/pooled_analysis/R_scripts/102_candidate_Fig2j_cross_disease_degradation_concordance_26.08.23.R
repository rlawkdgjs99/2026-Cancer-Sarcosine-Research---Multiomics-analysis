#!/usr/bin/env Rscript

# Candidate Figure 2 analysis: cross-disease concordance of species associations
# with microbiome sarcosine-degradation potential in CRC and NSCLC ICI cohorts.
#
# Biological question
#   Do the same gut microbial species track community sarcosine-degradation
#   potential across CRC and NSCLC cohorts?
#
# Primary design
#   - CRC independent unit: ENA sample_accession (BioSample). PAIRED runs are
#     averaged within BioSample before analysis; PRJEB10878 and PRJNA429097 have
#     one PAIRED run per BioSample. SINGLE+PAIRED mean and PAIRED median are
#     prespecified sensitivity analyses.
#   - NSCLC independent unit: score-selected baseline sample. Local metadata
#     must show one unique RunID and one unique Sample.name per included row.
#   - Species: species-level HGMT relative abundance, collapsed by species name;
#     >=10% prevalence within each biological cohort.
#   - Association: phenotype-adjusted partial Spearman correlation within cohort
#     (CRC: Healthy/Cancer; NSCLC: R/NR). Unadjusted Spearman is a sensitivity.
#   - Meta-analysis: Fisher-z correlations, DerSimonian-Laird random effects,
#     requiring >=2 cohorts per disease. The custom implementation is checked
#     against metafor::rma.uni(method="DL").
#   - Cross-disease concordance: Spearman correlation of CRC and NSCLC meta-rho
#     over prevalence-defined shared species; inference by 1,000 independent,
#     phenotype-stratified score permutations within every cohort.
#   - Multiplicity: BH-FDR across species within disease.
#
# Interpretation boundary
#   Taxonomy and KO-derived scores come from the same microbiomes. This analysis
#   describes reproducible ecological association, not taxonomic gene carriage,
#   metabolic flux, causal effect, or independent validation of the KO score.

set.seed(42)
options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(metafor)
})

required_packages <- c("ragg", "svglite", "digest")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

get_script_file <- function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(f)) stop("Run this analysis with Rscript.")
  normalizePath(sub("^--file=", "", f[1]), mustWork = TRUE)
}

SCRIPT_FILE <- get_script_file()
SCRIPT_DIR <- dirname(SCRIPT_FILE)
NSCLC_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = TRUE)
COLLECTION_ROOT <- dirname(NSCLC_ROOT)
CRC_ROOT <- file.path(
  COLLECTION_ROOT, "HGMT_CRC_WGS-Healthy_vs_Cancer",
  "Healthy_vs_Cancer_4_CRC_cohorts_integrated"
)
OUT_DIR <- file.path(
  NSCLC_ROOT, "pooled_analysis", "results",
  "Fig2_cross_disease_species_degradation_concordance_CANDIDATE_26.08.23"
)
ENA_DIR <- file.path(
  NSCLC_ROOT, "pooled_analysis", "resources",
  "ENA_CRC_run_sample_mapping_26.08.23"
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(OUT_DIR, "analysis_log.txt")
log_con <- file(LOG_FILE, open = "wt")
say <- function(...) {
  z <- sprintf(...)
  cat(z, "\n")
  cat(z, "\n", file = log_con)
}
on.exit({ if (isOpen(log_con)) close(log_con) }, add = TRUE)

PREV_MIN <- 0.10
N_PERM <- 1000L
DEG_KOS <- c("K00301", "K00302", "K00303", "K00305", "K00306")

CRC_DIRS <- c(
  PRJEB10878 = "PRJEB10878_CRC",
  PRJEB27928 = "PRJEB27928_CRC",
  PRJEB6070 = "PRJEB6070_CRC_AdenomatousPolyps",
  PRJNA429097 = "PRJNA429097_CRC"
)
NSCLC_DISC_DIRS <- c(
  PRJNA751792 = "NSCLC_PRJNA751792",
  PRJNA1023797 = "NSCLC_PRJNA1023797",
  PRJEB22863 = "NSCLC_RCC_PRJEB22863"
)
NSCLC_VALID_DIRS <- c(PRJEB26531 = "NSCLC_PRJEB26531")

if (!dir.exists(CRC_ROOT) || !dir.exists(ENA_DIR)) {
  stop("Required CRC or ENA mapping directory is missing.")
}

parse_hgmt_metadata <- function(path) {
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(x) < 3L) stop("Metadata has fewer than three lines: ", path)
  header <- strsplit(x[2], "\t")[[1]]
  rows <- strsplit(x[-c(1, 2)], "\t")
  rows <- rows[lengths(rows) > 1L]
  rows <- lapply(rows, function(z) {
    if (length(z) >= length(header)) z[seq_along(header)] else
      c(z, rep(NA_character_, length(header) - length(z)))
  })
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(out) <- make.names(header)
  out
}

pick_file <- function(dir, pattern, allow_latest = FALSE) {
  x <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (!length(x)) stop("No file matching ", pattern, " in ", dir)
  if (length(x) > 1L && !allow_latest) {
    stop("Expected one file matching ", pattern, " in ", dir, "; found ", length(x))
  }
  if (allow_latest) sort(x, decreasing = TRUE)[1] else x[1]
}

read_species_run_matrix <- function(path, allowed_runs) {
  b <- fread(path, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
  if (ncol(b) != 3L) stop("Bacteria table must have exactly three columns: ", path)
  setnames(b, c("Taxa", "RunID", "Abundance"))
  if (!is.numeric(b$Abundance) || anyNA(b$Abundance) ||
      any(!is.finite(b$Abundance)) || any(b$Abundance < 0)) {
    stop("Bacteria abundance is not finite, non-negative numeric data: ", path)
  }
  b <- b[
    RunID %in% allowed_runs &
      grepl("|s__", Taxa, fixed = TRUE) &
      !grepl("|t__", Taxa, fixed = TRUE)
  ]
  if (!nrow(b)) stop("No species-level rows after filtering: ", path)
  b[, species := sub(".*\\|s__", "", Taxa)]
  b <- b[, .(Abundance = sum(Abundance)), by = .(RunID, species)]
  w <- dcast(b, RunID ~ species, value.var = "Abundance", fill = 0)
  m <- as.matrix(w[, -1L])
  storage.mode(m) <- "double"
  rownames(m) <- w$RunID
  if (anyDuplicated(rownames(m)) || anyDuplicated(colnames(m))) {
    stop("Duplicated run/species identifiers after species collapse: ", path)
  }
  rs <- rowSums(m)
  attr(m, "run_sum_min") <- min(rs)
  attr(m, "run_sum_median") <- median(rs)
  attr(m, "run_sum_max") <- max(rs)
  m
}

aggregate_rows <- function(m, ids, fun = c("mean", "median")) {
  fun <- match.arg(fun)
  lev <- unique(ids)
  out <- matrix(NA_real_, nrow = length(lev), ncol = ncol(m),
                dimnames = list(lev, colnames(m)))
  for (i in seq_along(lev)) {
    z <- m[ids == lev[i], , drop = FALSE]
    out[i, ] <- if (fun == "mean") colMeans(z) else apply(z, 2, median)
  }
  out
}

load_crc_cohort <- function(cohort, folder, ko_cache) {
  cdir <- file.path(CRC_ROOT, folder)
  meta_file <- pick_file(cdir, "^selected_project_.*\\.txt$", allow_latest = TRUE)
  bact_file <- pick_file(cdir, "^Bacteria_.*\\.txt$")
  ena_file <- file.path(ENA_DIR, paste0(cohort, "_ena_runs.tsv"))
  if (!file.exists(ena_file)) stop("ENA run mapping missing: ", ena_file)

  meta <- parse_hgmt_metadata(meta_file)
  required <- c("Run.ID", "Sample.name", "Assay.type", "Phenotype.name", "Batch")
  if (length(setdiff(required, names(meta)))) stop("CRC metadata missing required columns: ", cohort)
  if (anyDuplicated(meta$Run.ID)) stop("Duplicated CRC Run.ID in metadata: ", cohort)
  meta <- meta[
    meta$Assay.type == "WGS" &
      meta$Phenotype.name %in% c("Health", "Colorectal Neoplasms"),
    , drop = FALSE
  ]
  meta$Group <- ifelse(meta$Phenotype.name == "Health", "Healthy", "Cancer")

  ena <- fread(ena_file, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
  required_ena <- c("run_accession", "sample_accession", "sample_alias", "library_layout")
  if (length(setdiff(required_ena, names(ena)))) stop("ENA mapping missing fields: ", cohort)
  ix <- match(meta$Run.ID, ena$run_accession)
  if (anyNA(ix)) stop("CRC metadata run absent from ENA mapping: ", cohort)
  meta$BioSample <- ena$sample_accession[ix]
  meta$ENA_alias <- ena$sample_alias[ix]
  meta$ENA_layout <- ena$library_layout[ix]
  alias_ok <- meta$Sample.name == meta$BioSample | meta$Sample.name == meta$ENA_alias
  if (!all(alias_ok)) stop("Local Sample.name does not match ENA accession or alias: ", cohort)
  if (any(!nzchar(meta$BioSample))) stop("Invalid BioSample mapping: ", cohort)

  const_group <- aggregate(Group ~ BioSample, meta, function(z) length(unique(z)))
  if (any(const_group$Group != 1L)) stop("Phenotype varies within CRC BioSample: ", cohort)
  bact <- read_species_run_matrix(bact_file, meta$Run.ID)
  ko <- ko_cache[ko_cache$Cohort == cohort, , drop = FALSE]
  if (anyDuplicated(ko$Run.ID)) stop("Duplicated Run.ID in KO cache: ", cohort)
  missing_kos <- setdiff(DEG_KOS, names(ko))
  if (length(missing_kos)) stop("KO cache missing degradation KO(s): ", paste(missing_kos, collapse = ", "))
  ko$degradation <- rowSums(ko[, DEG_KOS, drop = FALSE])

  list(
    cohort = cohort, meta = meta, bact = bact,
    score = setNames(ko$degradation, ko$Run.ID),
    input_files = c(meta_file, bact_file, ena_file),
    bact_sum = unlist(attributes(bact)[c("run_sum_min", "run_sum_median", "run_sum_max")])
  )
}

aggregate_crc_cohort <- function(raw, layout = c("paired", "all"), fun = c("mean", "median")) {
  layout <- match.arg(layout)
  fun <- match.arg(fun)
  map <- raw$meta
  if (layout == "paired") map <- map[map$ENA_layout == "PAIRED", , drop = FALSE]
  common_runs <- Reduce(intersect, list(map$Run.ID, rownames(raw$bact), names(raw$score)))
  map <- map[match(common_runs, map$Run.ID), , drop = FALSE]
  if (!identical(map$Run.ID, common_runs)) stop("CRC run alignment failed: ", raw$cohort)
  if (!nrow(map)) stop("No common CRC runs: ", raw$cohort)

  b <- raw$bact[common_runs, , drop = FALSE]
  s <- raw$score[common_runs]
  ids <- map$BioSample
  bm <- aggregate_rows(b, ids, fun = fun)
  score_m <- vapply(split(as.numeric(s), ids), function(z) {
    if (fun == "mean") mean(z) else median(z)
  }, numeric(1))
  score_m <- score_m[rownames(bm)]
  group_m <- vapply(split(map$Group, ids), function(z) unique(z), character(1))
  group_m <- group_m[rownames(bm)]
  if (anyNA(score_m) || anyNA(group_m)) stop("CRC BioSample aggregation alignment failed: ", raw$cohort)

  unit_qc <- data.frame(
    cohort = raw$cohort, BioSample = rownames(bm), group = group_m,
    n_runs = as.integer(table(ids)[rownames(bm)]), score = as.numeric(score_m),
    stringsAsFactors = FALSE
  )
  list(
    cohort = raw$cohort, species = bm, score = as.numeric(score_m), group = group_m,
    ids = rownames(bm), unit_qc = unit_qc,
    layout = layout, aggregation = fun
  )
}

load_nsclc_cohort <- function(cohort, folder) {
  cdir <- file.path(NSCLC_ROOT, folder)
  meta_file <- pick_file(cdir, "^selected_project_.*\\.txt$", allow_latest = TRUE)
  bact_file <- pick_file(cdir, "^Bacteria_.*\\.txt$")
  score_file <- file.path(cdir, "results", paste0("sarcosine_scores_", cohort, ".csv"))
  if (!file.exists(score_file)) stop("NSCLC score file missing: ", score_file)
  meta <- parse_hgmt_metadata(meta_file)
  sc <- fread(score_file)
  required_sc <- c("RunID", "group", "degradation")
  if (length(setdiff(required_sc, names(sc)))) stop("NSCLC score fields missing: ", cohort)
  if (anyDuplicated(sc$RunID) || anyNA(sc$degradation) || any(sc$degradation < 0)) {
    stop("Invalid NSCLC score table: ", cohort)
  }
  ix <- match(sc$RunID, meta$Run.ID)
  if (anyNA(ix)) stop("NSCLC score run absent from metadata: ", cohort)
  matched <- meta[ix, , drop = FALSE]
  if (anyDuplicated(matched$Sample.name)) stop("NSCLC analysis contains repeated Sample.name: ", cohort)
  bact <- read_species_run_matrix(bact_file, sc$RunID)
  if (!all(sc$RunID %in% rownames(bact))) stop("NSCLC score/taxonomy run mismatch: ", cohort)
  bact_sum <- unlist(attributes(bact)[c("run_sum_min", "run_sum_median", "run_sum_max")])
  bact <- bact[sc$RunID, , drop = FALSE]
  list(
    cohort = cohort, species = bact, score = sc$degradation,
    group = sc$group, ids = matched$Sample.name,
    unit_qc = data.frame(
      cohort = cohort, BioSample = matched$Sample.name, group = sc$group,
      n_runs = 1L, score = sc$degradation, stringsAsFactors = FALSE
    ),
    input_files = c(meta_file, bact_file, score_file),
    bact_sum = bact_sum
  )
}

residualize <- function(x, design) {
  q <- qr.Q(qr(design))
  x - q %*% crossprod(q, x)
}

association_cohort <- function(dat, adjusted = TRUE) {
  if (nrow(dat$species) != length(dat$score) || length(dat$score) != length(dat$group)) {
    stop("Association inputs are not aligned: ", dat$cohort)
  }
  if (length(unique(dat$group)) != 2L) stop("Expected two phenotype groups: ", dat$cohort)
  prev <- colMeans(dat$species > 0)
  keep <- is.finite(prev) & prev >= PREV_MIN
  m <- dat$species[, keep, drop = FALSE]
  if (!ncol(m)) stop("No species pass prevalence: ", dat$cohort)

  rx <- apply(m, 2, rank, ties.method = "average")
  if (is.null(dim(rx))) rx <- matrix(rx, ncol = 1L, dimnames = list(NULL, colnames(m)))
  colnames(rx) <- colnames(m)
  ry <- rank(dat$score, ties.method = "average")
  design <- if (adjusted) model.matrix(~ factor(dat$group)) else matrix(1, nrow(m), 1L)
  xres <- residualize(rx, design)
  yres <- as.numeric(residualize(matrix(ry, ncol = 1L), design))
  xnorm <- sqrt(colSums(xres^2))
  ynorm <- sqrt(sum(yres^2))
  rho <- as.numeric(crossprod(xres, yres)) / (xnorm * ynorm)
  rho[!is.finite(rho)] <- NA_real_
  covariates <- ncol(design) - 1L
  df <- nrow(m) - covariates - 2L
  if (df <= 1L) stop("Insufficient degrees of freedom: ", dat$cohort)
  t_stat <- rho * sqrt(df / pmax(1 - rho^2, .Machine$double.eps))
  p <- 2 * pt(abs(t_stat), df = df, lower.tail = FALSE)
  z_var <- 1 / (nrow(m) - covariates - 3L)

  table <- data.frame(
    cohort = dat$cohort, species = colnames(m), rho = rho,
    p_value = p, p_adj_BH = p.adjust(p, method = "BH"),
    prevalence = prev[colnames(m)], n = nrow(m),
    phenotype_adjusted = adjusted, n_covariates = covariates,
    fisher_z_variance = z_var, stringsAsFactors = FALSE
  )
  table <- table[is.finite(table$rho) & abs(table$rho) < 1, , drop = FALSE]
  state <- list(
    cohort = dat$cohort, species = colnames(m), xres = xres,
    xnorm = xnorm, score = dat$score, group = dat$group,
    design = design, n = nrow(m), covariates = covariates
  )
  list(table = table, state = state)
}

dl_meta_one <- function(rho, vi) {
  keep <- is.finite(rho) & abs(rho) < 1 & is.finite(vi) & vi > 0
  rho <- rho[keep]
  vi <- vi[keep]
  k <- length(rho)
  if (k < 2L) return(NULL)
  z <- atanh(rho)
  w <- 1 / vi
  mu_fixed <- sum(w * z) / sum(w)
  q <- sum(w * (z - mu_fixed)^2)
  cval <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (q - (k - 1)) / cval)
  wr <- 1 / (vi + tau2)
  mu <- sum(wr * z) / sum(wr)
  se <- sqrt(1 / sum(wr))
  p <- 2 * pnorm(abs(mu / se), lower.tail = FALSE)
  data.frame(
    k = k, pooled_z = mu, pooled_rho = tanh(mu), se_z = se,
    ci_lb = tanh(mu - qnorm(0.975) * se),
    ci_ub = tanh(mu + qnorm(0.975) * se),
    p_value = p, tau2 = tau2,
    I2 = ifelse(q > 0, max(0, (q - (k - 1)) / q) * 100, 0)
  )
}

meta_disease <- function(per_cohort, disease, min_k = 2L) {
  all <- rbindlist(lapply(per_cohort, `[[`, "table"), fill = TRUE)
  species <- all[, .N, by = species][N >= min_k, species]
  out <- rbindlist(lapply(species, function(sp) {
    d <- all[species == sp]
    m <- dl_meta_one(d$rho, d$fisher_z_variance)
    if (is.null(m)) return(NULL)
    cbind(data.frame(disease = disease, species = sp, stringsAsFactors = FALSE), m)
  }), fill = TRUE)
  if (!nrow(out)) stop("No meta-analysable species for ", disease)
  out[, q_value_BH := p.adjust(p_value, method = "BH")]
  out[, abs_pooled_rho_for_order := abs(pooled_rho)]
  setorder(out, q_value_BH, -abs_pooled_rho_for_order)
  out[, abs_pooled_rho_for_order := NULL]
  out
}

validate_meta_api <- function(meta_table, per_cohort, disease) {
  all <- rbindlist(lapply(per_cohort, `[[`, "table"), fill = TRUE)
  checks <- rbindlist(lapply(meta_table$species, function(sp) {
    d <- all[species == sp]
    fit <- metafor::rma.uni(yi = atanh(d$rho), vi = d$fisher_z_variance,
                            method = "DL", test = "z")
    ours <- meta_table[species == sp]
    data.table(
      disease = disease, species = sp,
      custom_z = ours$pooled_z, metafor_z = as.numeric(fit$b),
      abs_diff_z = abs(ours$pooled_z - as.numeric(fit$b)),
      custom_se = ours$se_z, metafor_se = fit$se,
      abs_diff_se = abs(ours$se_z - fit$se),
      custom_tau2 = ours$tau2, metafor_tau2 = fit$tau2,
      abs_diff_tau2 = abs(ours$tau2 - fit$tau2)
    )
  }))
  if (max(checks$abs_diff_z, checks$abs_diff_se, checks$abs_diff_tau2) > 1e-8) {
    stop("Custom DL implementation failed metafor validation for ", disease)
  }
  checks
}

concordance_summary <- function(crc_meta, nsc_meta, label) {
  z <- merge(
    crc_meta[, .(species, crc_rho = pooled_rho, crc_q = q_value_BH, crc_k = k)],
    nsc_meta[, .(species, nsclc_rho = pooled_rho, nsclc_q = q_value_BH, nsclc_k = k)],
    by = "species"
  )
  if (nrow(z) < 10L) stop("Too few shared species in scenario: ", label)
  ct <- suppressWarnings(cor.test(z$crc_rho, z$nsclc_rho, method = "spearman", exact = FALSE))
  both <- z$crc_q < 0.05 & z$nsclc_q < 0.05
  data.frame(
    scenario = label, n_shared_species = nrow(z),
    spearman_rho = unname(ct$estimate), naive_species_p = ct$p.value,
    sign_concordance = mean(sign(z$crc_rho) == sign(z$nsclc_rho)),
    both_BH_q_lt_0.05 = sum(both),
    both_FDR_sign_concordance = ifelse(
      sum(both) > 0,
      mean(sign(z$crc_rho[both]) == sign(z$nsclc_rho[both])), NA_real_
    ), stringsAsFactors = FALSE
  )
}

meta_from_rho_matrix <- function(rho, vi) {
  z <- atanh(pmin(pmax(rho, -0.999999), 0.999999))
  ok <- is.finite(z) & is.finite(vi)
  w <- ifelse(ok, 1 / vi, 0)
  z0 <- ifelse(ok, z, 0)
  sw <- colSums(w)
  mu_fixed <- colSums(w * z0) / sw
  mu_mat <- matrix(mu_fixed, nrow(rho), ncol(rho), byrow = TRUE)
  q <- colSums(w * (z0 - mu_mat)^2)
  k <- colSums(ok)
  cval <- sw - colSums(w^2) / sw
  tau2 <- pmax(0, (q - (k - 1)) / cval)
  tau_mat <- matrix(tau2, nrow(rho), ncol(rho), byrow = TRUE)
  wr <- ifelse(ok, 1 / (vi + tau_mat), 0)
  mu <- colSums(wr * z0) / colSums(wr)
  tanh(mu)
}

permuted_rho <- function(state) {
  score_perm <- state$score
  for (g in unique(state$group)) {
    ix <- which(state$group == g)
    score_perm[ix] <- sample(score_perm[ix], length(ix), replace = FALSE)
  }
  ry <- rank(score_perm, ties.method = "average")
  yres <- as.numeric(residualize(matrix(ry, ncol = 1L), state$design))
  as.numeric(crossprod(state$xres, yres)) / (state$xnorm * sqrt(sum(yres^2)))
}

build_meta_matrices <- function(per_cohort, species) {
  k <- length(per_cohort)
  rho <- matrix(NA_real_, k, length(species), dimnames = list(names(per_cohort), species))
  vi <- matrix(NA_real_, k, length(species), dimnames = list(names(per_cohort), species))
  for (i in seq_along(per_cohort)) {
    d <- per_cohort[[i]]$table
    j <- match(d$species, species)
    keep <- !is.na(j)
    rho[i, j[keep]] <- d$rho[keep]
    vi[i, j[keep]] <- d$fisher_z_variance[keep]
  }
  list(rho = rho, vi = vi)
}

permute_concordance <- function(crc_assoc, nsc_assoc, shared_species, observed, reps = N_PERM) {
  crc_template <- build_meta_matrices(crc_assoc, shared_species)
  nsc_template <- build_meta_matrices(nsc_assoc, shared_species)
  out <- numeric(reps)
  for (b in seq_len(reps)) {
    cr <- crc_template$rho
    ns <- nsc_template$rho
    for (i in seq_along(crc_assoc)) {
      rr <- permuted_rho(crc_assoc[[i]]$state)
      j <- match(crc_assoc[[i]]$state$species, shared_species)
      keep <- !is.na(j)
      cr[i, j[keep]] <- rr[keep]
    }
    for (i in seq_along(nsc_assoc)) {
      rr <- permuted_rho(nsc_assoc[[i]]$state)
      j <- match(nsc_assoc[[i]]$state$species, shared_species)
      keep <- !is.na(j)
      ns[i, j[keep]] <- rr[keep]
    }
    cm <- meta_from_rho_matrix(cr, crc_template$vi)
    nm <- meta_from_rho_matrix(ns, nsc_template$vi)
    out[b] <- suppressWarnings(cor(cm, nm, method = "spearman", use = "complete.obs"))
    if (b %% 100L == 0L) say("Permutation %d/%d", b, reps)
  }
  p <- (1 + sum(abs(out) >= abs(observed))) / (reps + 1)
  list(null = out, p = p)
}

save_plot <- function(p, stem, width = 4.2, height = 3.8) {
  png_file <- file.path(OUT_DIR, paste0(stem, ".png"))
  pdf_file <- file.path(OUT_DIR, paste0(stem, ".pdf"))
  svg_file <- file.path(OUT_DIR, paste0(stem, ".svg"))
  tmp_png <- tempfile(fileext = ".png")
  tmp_pdf <- tempfile(fileext = ".pdf")
  tmp_svg <- tempfile(fileext = ".svg")
  ragg::agg_png(tmp_png, width = width, height = height, units = "in", res = 600,
                background = "white")
  print(p); dev.off()
  grDevices::pdf(tmp_pdf, width = width, height = height,
                 family = "Helvetica", useDingbats = FALSE)
  print(p); dev.off()
  svglite::svglite(tmp_svg, width = width, height = height, bg = "white")
  print(p); dev.off()
  stopifnot(file.copy(tmp_png, png_file, overwrite = TRUE))
  stopifnot(file.copy(tmp_pdf, pdf_file, overwrite = TRUE))
  stopifnot(file.copy(tmp_svg, svg_file, overwrite = TRUE))
  unlink(c(tmp_png, tmp_pdf, tmp_svg))
  c(png_file, pdf_file, svg_file)
}

say("[1/9] Loading CRC KO cache and validating value space")
KO_CACHE_FILE <- file.path(
  CRC_ROOT, "results_integrated", "sarcosine", "sarcosine_KO_per_sample_pooled.csv"
)
ko_cache <- read.csv(KO_CACHE_FILE, check.names = FALSE, stringsAsFactors = FALSE)
if (anyDuplicated(ko_cache$Run.ID) || anyNA(ko_cache) ||
    any(as.matrix(ko_cache[, DEG_KOS, drop = FALSE]) < 0)) {
  stop("CRC KO cache failed uniqueness, missingness, or non-negativity checks.")
}
say("CRC KO cache: %d runs x %d fields; range degradation KOs %.6g to %.6g",
    nrow(ko_cache), ncol(ko_cache), min(as.matrix(ko_cache[, DEG_KOS])),
    max(as.matrix(ko_cache[, DEG_KOS])))

say("[2/9] Loading CRC metadata, official ENA mappings, and species profiles")
crc_raw <- Map(load_crc_cohort, names(CRC_DIRS), unname(CRC_DIRS),
               MoreArgs = list(ko_cache = ko_cache))
names(crc_raw) <- names(CRC_DIRS)
crc_primary <- lapply(crc_raw, aggregate_crc_cohort, layout = "paired", fun = "mean")
crc_all_mean <- lapply(crc_raw, aggregate_crc_cohort, layout = "all", fun = "mean")
crc_paired_median <- lapply(crc_raw, aggregate_crc_cohort, layout = "paired", fun = "median")
for (co in names(crc_primary)) {
  d <- crc_primary[[co]]
  say("CRC %s primary: %d BioSamples (%s), %d species, runs/sample %d-%d",
      co, length(d$score), paste(names(table(d$group)), table(d$group), collapse = "/"),
      ncol(d$species), min(d$unit_qc$n_runs), max(d$unit_qc$n_runs))
}

say("[3/9] Loading NSCLC discovery and validation data")
nsc_disc <- Map(load_nsclc_cohort, names(NSCLC_DISC_DIRS), unname(NSCLC_DISC_DIRS))
names(nsc_disc) <- names(NSCLC_DISC_DIRS)
nsc_valid <- Map(load_nsclc_cohort, names(NSCLC_VALID_DIRS), unname(NSCLC_VALID_DIRS))
names(nsc_valid) <- names(NSCLC_VALID_DIRS)
for (co in names(nsc_disc)) {
  d <- nsc_disc[[co]]
  say("NSCLC %s: %d unique samples (%s), %d species",
      co, length(d$score), paste(names(table(d$group)), table(d$group), collapse = "/"), ncol(d$species))
}

run_assoc_set <- function(x, adjusted) {
  out <- lapply(x, association_cohort, adjusted = adjusted)
  names(out) <- names(x)
  out
}

say("[4/9] Computing phenotype-adjusted per-cohort associations and DL meta-analysis")
crc_assoc <- run_assoc_set(crc_primary, adjusted = TRUE)
nsc_assoc <- run_assoc_set(nsc_disc, adjusted = TRUE)
crc_meta <- meta_disease(crc_assoc, "CRC", min_k = 2L)
nsc_meta <- meta_disease(nsc_assoc, "NSCLC", min_k = 2L)
api_checks <- rbind(
  validate_meta_api(crc_meta, crc_assoc, "CRC"),
  validate_meta_api(nsc_meta, nsc_assoc, "NSCLC")
)
say("Meta implementation max absolute difference vs metafor: %.3g",
    max(api_checks$abs_diff_z, api_checks$abs_diff_se, api_checks$abs_diff_tau2))

primary <- concordance_summary(crc_meta, nsc_meta,
                               "Primary: paired-run BioSample mean, phenotype-adjusted")
merged <- merge(
  crc_meta[, .(species, CRC_meta_rho = pooled_rho, CRC_ci_lb = ci_lb, CRC_ci_ub = ci_ub,
               CRC_p = p_value, CRC_q = q_value_BH, CRC_k = k, CRC_I2 = I2)],
  nsc_meta[, .(species, NSCLC_meta_rho = pooled_rho, NSCLC_ci_lb = ci_lb, NSCLC_ci_ub = ci_ub,
               NSCLC_p = p_value, NSCLC_q = q_value_BH, NSCLC_k = k, NSCLC_I2 = I2)],
  by = "species"
)
if (nrow(merged) != primary$n_shared_species) stop("Primary merge count mismatch.")

say("[5/9] Phenotype-stratified permutation test (%d repetitions)", N_PERM)
perm <- permute_concordance(
  crc_assoc, nsc_assoc, merged$species,
  observed = primary$spearman_rho, reps = N_PERM
)
primary$permutation_reps <- N_PERM
primary$stratified_permutation_p <- perm$p
primary$null_rho_mean <- mean(perm$null)
primary$null_rho_sd <- sd(perm$null)
say("Primary concordance: n=%d species, rho=%.6f, permutation P=%.6g, sign concordance=%.3f",
    primary$n_shared_species, primary$spearman_rho,
    primary$stratified_permutation_p, primary$sign_concordance)

say("[6/9] Running sensitivity analyses")
scenario <- function(crc_data, nsc_data, adjusted, label, min_crc = 2L, min_nsc = 2L) {
  ca <- run_assoc_set(crc_data, adjusted)
  na <- run_assoc_set(nsc_data, adjusted)
  concordance_summary(
    meta_disease(ca, "CRC", min_crc), meta_disease(na, "NSCLC", min_nsc), label
  )
}
primary_base <- primary[, names(concordance_summary(crc_meta, nsc_meta, "x")), drop = FALSE]
sensitivity <- rbind(
  primary_base,
  scenario(crc_primary, nsc_disc, FALSE, "Unadjusted Spearman"),
  scenario(crc_all_mean, nsc_disc, TRUE, "All-layout BioSample mean, phenotype-adjusted"),
  scenario(crc_paired_median, nsc_disc, TRUE, "Paired-run BioSample median, phenotype-adjusted"),
  concordance_summary(
    meta_disease(crc_assoc, "CRC", 4L), meta_disease(nsc_assoc, "NSCLC", 3L),
    "Species testable in all 4 CRC and all 3 NSCLC cohorts"
  )
)

loo <- list()
for (co in names(crc_assoc)) {
  loo[[paste0("drop_CRC_", co)]] <- concordance_summary(
    meta_disease(crc_assoc[names(crc_assoc) != co], "CRC", 2L), nsc_meta,
    paste0("Leave out CRC ", co)
  )
}
for (co in names(nsc_assoc)) {
  loo[[paste0("drop_NSCLC_", co)]] <- concordance_summary(
    crc_meta, meta_disease(nsc_assoc[names(nsc_assoc) != co], "NSCLC", 2L),
    paste0("Leave out NSCLC ", co)
  )
}
loo <- rbindlist(loo, fill = TRUE)

valid_assoc <- run_assoc_set(nsc_valid, adjusted = TRUE)[[1]]$table
crc_valid <- merge(crc_meta[, .(species, meta_rho = pooled_rho)],
                   as.data.table(valid_assoc)[, .(species, valid_rho = rho)], by = "species")
nsc_valid_cmp <- merge(nsc_meta[, .(species, meta_rho = pooled_rho)],
                       as.data.table(valid_assoc)[, .(species, valid_rho = rho)], by = "species")
validation_concordance <- rbind(
  data.frame(
    comparison = "CRC meta vs Korean NSCLC validation", n_species = nrow(crc_valid),
    spearman_rho = cor(crc_valid$meta_rho, crc_valid$valid_rho, method = "spearman"),
    p_value = suppressWarnings(cor.test(crc_valid$meta_rho, crc_valid$valid_rho,
                                        method = "spearman", exact = FALSE)$p.value)
  ),
  data.frame(
    comparison = "NSCLC discovery meta vs Korean NSCLC validation", n_species = nrow(nsc_valid_cmp),
    spearman_rho = cor(nsc_valid_cmp$meta_rho, nsc_valid_cmp$valid_rho, method = "spearman"),
    p_value = suppressWarnings(cor.test(nsc_valid_cmp$meta_rho, nsc_valid_cmp$valid_rho,
                                        method = "spearman", exact = FALSE)$p.value)
  )
)

say("[7/9] Writing tables and candidate figure")
merged$significance <- ifelse(
  merged$CRC_q < 0.05 & merged$NSCLC_q < 0.05, "FDR < 0.05 in both",
  ifelse(merged$CRC_q < 0.05, "CRC FDR only",
         ifelse(merged$NSCLC_q < 0.05, "NSCLC FDR only", "Neither FDR"))
)
merged$significance <- factor(
  merged$significance,
  levels = c("Neither FDR", "CRC FDR only", "NSCLC FDR only", "FDR < 0.05 in both")
)
label_map <- c(
  Lachnospira_eligens = "L. eligens",
  Roseburia_faecis = "R. faecis",
  Clostridiales_bacterium_KLE1615 = "KLE1615",
  Clostridium_sp_AF36_4 = "Clostridium sp. AF36_4"
)
merged$plot_label <- unname(label_map[merged$species])
merged$plot_label[is.na(merged$plot_label)] <- ""
merged$plot_significance <- factor(
  as.character(merged$significance), levels = levels(merged$significance),
  labels = c("Neither", "CRC only", "NSCLC only", "Both FDR < 0.05")
)

pal <- c(
  "Neither" = "#BDBDBD", "CRC only" = "#1B9E77",
  "NSCLC only" = "#377EB8", "Both FDR < 0.05" = "#7B3294"
)
lim <- max(abs(c(merged$CRC_meta_rho, merged$NSCLC_meta_rho)))
lim <- min(1, ceiling((lim + 0.03) * 10) / 10)
annot <- sprintf(
  "Shared species = %d\nSpearman rho = %.2f\nStratified permutation P %s\nSign concordance = %.0f%%",
  nrow(merged), primary$spearman_rho,
  ifelse(primary$stratified_permutation_p < 0.001, "< 0.001",
         paste0("= ", formatC(primary$stratified_permutation_p, digits = 3, format = "f"))),
  100 * primary$sign_concordance
)
p <- ggplot(merged, aes(CRC_meta_rho, NSCLC_meta_rho)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.35) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey65", linewidth = 0.35) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "#333333", fill = "grey70", linewidth = 0.5, alpha = 0.20) +
  geom_point(aes(fill = plot_significance), shape = 21, colour = "white",
             stroke = 0.25, size = 2.0, alpha = 0.90) +
  ggrepel::geom_text_repel(
    data = merged[nzchar(merged$plot_label), , drop = FALSE],
    aes(label = plot_label), size = 2.25, fontface = "italic",
    box.padding = 0.28, point.padding = 0.15, min.segment.length = 0,
    segment.colour = "grey55", segment.size = 0.25, max.overlaps = Inf,
    force = 1.4, seed = 42, show.legend = FALSE
  ) +
  annotate("text", x = lim * 0.96, y = -lim * 0.96, label = annot,
           hjust = 1, vjust = 0, size = 2.35, lineheight = 1.06) +
  scale_fill_manual(values = pal, drop = FALSE) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim), clip = "off") +
  labs(
    title = "Cross-cancer species-degradation concordance",
    x = "CRC meta-rho (partial Spearman)",
    y = "NSCLC ICI meta-rho (partial Spearman)",
    fill = NULL,
    caption = paste(
      "Phenotype-adjusted within cohorts; DL random-effects meta-analysis.",
      "Association, not taxonomic KO attribution.", sep = "\n"
    )
  ) +
  theme_classic(base_family = "Helvetica", base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", size = 9.5, hjust = 0.5,
                              margin = margin(b = 4)),
    axis.title = element_text(size = 8.5, colour = "black"),
    axis.text = element_text(size = 7.5, colour = "black"),
    axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
    legend.position = "bottom", legend.direction = "horizontal",
    legend.text = element_text(size = 6.7), legend.key.height = grid::unit(2.8, "mm"),
    legend.spacing.x = grid::unit(1.5, "mm"),
    plot.caption = element_text(size = 5.8, colour = "grey35", hjust = 0, lineheight = 1.02),
    plot.margin = margin(5, 5, 5, 5)
  )
figure_files <- save_plot(
  p, "Proposed_Fig2j_cross_disease_degradation_concordance",
  width = 4.35, height = 4.25
)

null_df <- data.frame(permutation = seq_along(perm$null), null_spearman_rho = perm$null)
null_plot <- ggplot(null_df, aes(null_spearman_rho)) +
  geom_histogram(bins = 40, fill = "grey75", colour = "white") +
  geom_vline(xintercept = primary$spearman_rho, colour = "#B2182B", linewidth = 0.8) +
  labs(x = "Permuted cross-disease Spearman rho", y = "Count",
       title = "Phenotype-stratified permutation null",
       subtitle = sprintf("Observed rho = %.3f; P = %.4g; %d permutations",
                          primary$spearman_rho, primary$stratified_permutation_p, N_PERM)) +
  theme_classic(base_size = 10)
invisible(save_plot(null_plot, "QC_stratified_permutation_null", width = 4.5, height = 3.2))

fwrite(rbindlist(lapply(crc_assoc, `[[`, "table")),
       file.path(OUT_DIR, "CRC_species_degradation_partial_spearman_percohort.csv"))
fwrite(rbindlist(lapply(nsc_assoc, `[[`, "table")),
       file.path(OUT_DIR, "NSCLC_species_degradation_partial_spearman_percohort.csv"))
fwrite(crc_meta, file.path(OUT_DIR, "CRC_species_degradation_DL_meta.csv"))
fwrite(nsc_meta, file.path(OUT_DIR, "NSCLC_species_degradation_DL_meta.csv"))
fwrite(merged, file.path(OUT_DIR, "cross_disease_concordance_plot_source_data.csv"))
fwrite(primary, file.path(OUT_DIR, "cross_disease_concordance_primary_statistics.csv"))
fwrite(sensitivity, file.path(OUT_DIR, "cross_disease_concordance_sensitivity_statistics.csv"))
fwrite(loo, file.path(OUT_DIR, "cross_disease_concordance_leave_one_cohort_out.csv"))
fwrite(validation_concordance, file.path(OUT_DIR, "Korean_validation_concordance.csv"))
fwrite(null_df, file.path(OUT_DIR, "stratified_permutation_null.csv"))
fwrite(api_checks, file.path(OUT_DIR, "DL_meta_custom_vs_metafor_validation.csv"))

unit_qc <- rbindlist(c(
  lapply(crc_primary, `[[`, "unit_qc"),
  lapply(nsc_disc, `[[`, "unit_qc"),
  lapply(nsc_valid, `[[`, "unit_qc")
), fill = TRUE)
fwrite(unit_qc, file.path(OUT_DIR, "independent_unit_source_data.csv"))
unit_summary <- unit_qc[, .(
  n_biological_units = uniqueN(BioSample), n_runs = sum(n_runs),
  min_runs_per_unit = min(n_runs), max_runs_per_unit = max(n_runs),
  degradation_score_min = min(score),
  degradation_score_median = median(score),
  degradation_score_max = max(score),
  degradation_score_zero_n = sum(score == 0)
), by = .(cohort, group)]
fwrite(unit_summary, file.path(OUT_DIR, "independent_unit_QC_summary.csv"))

bact_qc <- rbindlist(c(
  lapply(crc_raw, function(x) data.table(
    disease = "CRC", cohort = x$cohort,
    species_sum_min = x$bact_sum[1], species_sum_median = x$bact_sum[2], species_sum_max = x$bact_sum[3]
  )),
  lapply(c(nsc_disc, nsc_valid), function(x) data.table(
    disease = "NSCLC", cohort = x$cohort,
    species_sum_min = x$bact_sum[1], species_sum_median = x$bact_sum[2], species_sum_max = x$bact_sum[3]
  ))
))
fwrite(bact_qc, file.path(OUT_DIR, "species_abundance_value_space_QC.csv"))

say("[8/9] Writing provenance, checksums, and reproducibility records")
input_files <- unique(c(
  KO_CACHE_FILE, unlist(lapply(crc_raw, `[[`, "input_files")),
  unlist(lapply(nsc_disc, `[[`, "input_files")),
  unlist(lapply(nsc_valid, `[[`, "input_files")), SCRIPT_FILE
))
project_root <- dirname(COLLECTION_ROOT)
input_hash <- data.table(
  file = vapply(input_files, function(x) {
    z <- normalizePath(x)
    if (startsWith(z, paste0(project_root, .Platform$file.sep))) {
      substring(z, nchar(project_root) + 2L)
    } else z
  }, character(1)),
  bytes = file.info(input_files)$size,
  sha256 = vapply(input_files, digest::digest, character(1), algo = "sha256", file = TRUE)
)
fwrite(input_hash, file.path(OUT_DIR, "input_sha256.tsv"), sep = "\t")

readme <- c(
  "# Candidate Figure 2: CRC–NSCLC species–degradation concordance",
  "",
  "Status: candidate only; not installed in the manuscript or PowerPoint.",
  "",
  "## Primary result",
  sprintf("- Shared prevalence-defined species: %d", primary$n_shared_species),
  sprintf("- Cross-disease Spearman rho: %.6f", primary$spearman_rho),
  sprintf("- Phenotype-stratified permutation P: %.6g (%d permutations, seed 42)",
          primary$stratified_permutation_p, N_PERM),
  sprintf("- Sign concordance: %.1f%%", 100 * primary$sign_concordance),
  sprintf("- Species significant at BH q < 0.05 in both diseases: %d", primary$both_BH_q_lt_0.05),
  "",
  "## Independent unit and preprocessing",
  "- CRC: official ENA sample_accession; PAIRED runs averaged within BioSample before score and correlation.",
  "- NSCLC: one unique RunID and local Sample.name per score-selected baseline sample.",
  "- Species prevalence threshold: >=10% within cohort, applied after CRC BioSample aggregation.",
  "- Primary association: phenotype-adjusted partial Spearman; CRC adjusts Healthy/Cancer and NSCLC adjusts R/NR.",
  "- Meta-analysis: Fisher z with DerSimonian-Laird random effects; species present in >=2 cohorts per disease.",
  "- Cross-disease inference: independent score permutations within phenotype strata in every cohort.",
  "",
  "## Interpretation boundary",
  "Species profiles and KO-derived degradation scores originate from the same microbiomes. The result is ecological concordance, not proof that a named taxon carries the KOs, not metabolic flux, not causation, and not an independent validation dataset.",
  "",
  "## Provenance",
  "Official ENA run/sample mapping snapshots are stored under pooled_analysis/resources/ENA_CRC_run_sample_mapping_26.08.23.",
  "The CRC degradation score is reconstructed from the existing sarcosine_KO_per_sample_pooled.csv cache using K00301, K00302, K00303, K00305, and K00306.",
  "See input_sha256.tsv, independent_unit_QC_summary.csv, sensitivity tables, permutation null, and sessionInfo.txt."
)
writeLines(readme, file.path(OUT_DIR, "README.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"), useBytes = TRUE)

say("[9/9] Final sanity checks")
if (!all(file.exists(figure_files))) stop("Candidate figure export missing.")
if (!is.finite(primary$spearman_rho) || !is.finite(primary$stratified_permutation_p)) {
  stop("Non-finite primary result.")
}
if (any(sensitivity$n_shared_species < 10L) || any(loo$n_shared_species < 10L)) {
  stop("A sensitivity analysis has too few shared species.")
}
if (any(!is.finite(sensitivity$spearman_rho)) || any(!is.finite(loo$spearman_rho))) {
  stop("A sensitivity analysis produced a non-finite concordance.")
}

stable_outputs <- list.files(OUT_DIR, full.names = TRUE)
stable_outputs <- stable_outputs[!basename(stable_outputs) %in% c("analysis_log.txt", "output_sha256.tsv")]
output_hash <- data.table(
  file = basename(stable_outputs), bytes = file.info(stable_outputs)$size,
  sha256 = vapply(stable_outputs, digest::digest, character(1), algo = "sha256", file = TRUE)
)
fwrite(output_hash, file.path(OUT_DIR, "output_sha256.tsv"), sep = "\t")

say("DONE: %s", OUT_DIR)
say("Primary rho %.6f | permutation P %.6g | shared species %d",
    primary$spearman_rho, primary$stratified_permutation_p, primary$n_shared_species)
close(log_con)
