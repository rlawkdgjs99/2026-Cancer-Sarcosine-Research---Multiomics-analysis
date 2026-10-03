# Functions copied verbatim from script102, lines106-251 and285-490.
# Source path and SHA256 are recorded in helper_provenance.json.
# No top-level analysis or old cohort definitions copied.
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
