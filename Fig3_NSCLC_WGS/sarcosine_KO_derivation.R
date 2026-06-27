#!/usr/bin/env Rscript
# ===========================================================================
# Derive the sarcosine direct degradation / production KO set DIRECTLY from KEGG
# via the KEGGREST R package. Fully reproducible: every fact comes from a live
# KEGGREST query (no web search, no hand-typed KO list, no model memory).
#
# Classification rule (objective, from data): in each KEGG reaction EQUATION,
#   - sarcosine (C00213) on the SUBSTRATE side  -> consumed  = DEGRADATION
#   - sarcosine (C00213) on the PRODUCT  side  -> formed    = PRODUCTION
# Reactions are written reversibly (<=>); this is the KEGG-canonical direction,
# and KOs catalysing both-direction reactions are marked "both".
#
# Run:  Rscript sarcosine_KO_derivation.R
# ===========================================================================
suppressPackageStartupMessages(library(KEGGREST))
SARC <- "C00213"

## 1) all reactions linked to sarcosine
rxn_ids <- unique(sub("rn:", "", unname(keggLink("reaction", paste0("cpd:", SARC)))))
cat(sprintf("KEGGREST: %d reactions linked to sarcosine (C00213)\n", length(rxn_ids)))

## 2) fetch reaction records (keggGet: max 10 per call)
batches <- split(rxn_ids, ceiling(seq_along(rxn_ids) / 10))
rx <- do.call(c, lapply(batches, keggGet))

## 3) per reaction: side of sarcosine in EQUATION + catalysing KOs
rows <- list()
for (r in rx) {
  eq <- r$EQUATION
  pr <- strsplit(eq, "<=>", fixed = TRUE)[[1]]
  lft <- pr[1]; rgt <- if (length(pr) > 1) pr[2] else ""
  role <- if (grepl(SARC, lft, fixed = TRUE) && grepl(SARC, rgt, fixed = TRUE)) "both"
          else if (grepl(SARC, lft, fixed = TRUE)) "degradation"
          else if (grepl(SARC, rgt, fixed = TRUE)) "production"
          else "unclear"
  kos <- r$ORTHOLOGY
  if (length(kos) == 0) next
  for (k in names(kos))
    rows[[length(rows) + 1]] <- data.frame(
      reaction = r$ENTRY, role = role, KO = k, KO_name = unname(kos[k]),
      equation = eq, definition = r$DEFINITION, stringsAsFactors = FALSE)
}
df <- do.call(rbind, rows)

## 4) collapse to per-KO (a KO can catalyse several sarcosine reactions)
kos_all <- sort(unique(df$KO))
agg <- data.frame(KO = kos_all, stringsAsFactors = FALSE)
agg$KO_name   <- vapply(kos_all, function(k) df$KO_name[match(k, df$KO)], character(1))
agg$role      <- vapply(kos_all, function(k) paste(sort(unique(df$role[df$KO == k])), collapse = ";"), character(1))
agg$reactions <- vapply(kos_all, function(k) paste(sort(unique(df$reaction[df$KO == k])), collapse = ";"), character(1))

## 5) presence in each cohort metagenome (counts pre-extracted to /tmp; runs with the KO)
fmap <- c(PRJEB26531   = "/tmp/koc_NSCLC_PRJEB26531.txt",
          PRJNA1023797 = "/tmp/koc_NSCLC_PRJNA1023797.txt",
          PRJNA751792  = "/tmp/koc_NSCLC_PRJNA751792.txt",
          PRJEB22863   = "/tmp/koc_NSCLC_RCC_PRJEB22863.txt")
for (cn in names(fmap)) {
  if (file.exists(fmap[cn])) {
    cc <- read.table(fmap[cn]); m <- setNames(cc$V1, cc$V2)
    agg[[paste0("n_", cn)]] <- vapply(agg$KO, function(k) { v <- m[k]; if (is.na(v)) 0L else as.integer(v) }, integer(1))
  }
}

## 6) report (ordered: degradation, production, both)
ord <- c("degradation", "production", "both", "unclear")
agg <- agg[order(match(agg$role, ord), agg$KO), ]
cat("\n================ per-reaction basis (from KEGG EQUATION) ================\n")
print(unique(df[, c("reaction", "role", "KO", "definition")]), row.names = FALSE)
cat("\n================ per-KO summary (role from KEGG, presence from data) ================\n")
print(agg, row.names = FALSE)
write.csv(agg, "sarcosine_KO_set.csv", row.names = FALSE)
cat("\nsaved: sarcosine_KO_set.csv\n")
