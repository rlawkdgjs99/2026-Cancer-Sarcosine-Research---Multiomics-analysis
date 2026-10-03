#!/usr/bin/env Rscript
# ===========================================================================
# INDEPENDENT re-verification of the sarcosine KO set (KEGGREST only).
# Route A (sarcosine_KO_derivation.R): C00213 -> reactions -> KOs (+ equation side).
# Route B (here): reverse direction. For EACH KO, query its OWN reactions
#   (keggLink reaction <- ko) and confirm the intersection with sarcosine's
#   reactions is non-empty -> the KO genuinely acts on sarcosine. Also pull the
#   KO NAME from KEGG to confirm the enzyme identity. Agreement of A & B = robust.
# Run:  Rscript sarcosine_KO_verify.R
# ===========================================================================
suppressPackageStartupMessages(library(KEGGREST))
SARC <- "C00213"

stopifnot(file.exists("sarcosine_KO_set.csv"))
ko_set <- read.csv("sarcosine_KO_set.csv", stringsAsFactors = FALSE)

## Route A reference: sarcosine's own reactions
sarc_rxn <- sort(unique(sub("rn:", "", unname(keggLink("reaction", paste0("cpd:", SARC))))))
cat(sprintf("Sarcosine (C00213) reactions [Route A reference]: %d\n  %s\n\n",
            length(sarc_rxn), paste(sarc_rxn, collapse = ", ")))

## Route B: per-KO reverse confirmation
rows <- list()
for (i in seq_len(nrow(ko_set))) {
  k <- ko_set$KO[i]
  kr <- tryCatch(sort(unique(sub("rn:", "", unname(keggLink("reaction", paste0("ko:", k)))))),
                 error = function(e) character(0))
  inter <- intersect(kr, sarc_rxn)
  nm <- tryCatch(paste(keggGet(paste0("ko:", k))[[1]]$NAME, collapse = "; "),
                 error = function(e) NA_character_)
  rows[[i]] <- data.frame(
    KO = k, roleA = ko_set$role[i],
    n_total_rxn_of_KO = length(kr),
    sarc_rxn_hit = paste(inter, collapse = ";"),
    n_sarc_hit = length(inter),
    confirmed = length(inter) > 0,
    KEGG_KO_name = ifelse(is.na(nm), ko_set$KO_name[i], substr(nm, 1, 55)),
    stringsAsFactors = FALSE)
}
v <- do.call(rbind, rows)
print(v[, c("KO", "roleA", "n_sarc_hit", "confirmed", "sarc_rxn_hit", "KEGG_KO_name")], row.names = FALSE)

cat("\n--- AGREEMENT CHECK ---\n")
cat("All KOs in set independently confirmed (Route B intersects sarcosine reactions):",
    all(v$confirmed), "\n")
nc <- v$KO[!v$confirmed]
cat("KOs NOT confirmed by Route B:", if (length(nc)) paste(nc, collapse = ", ") else "(none)", "\n")
write.csv(v, "sarcosine_KO_verification.csv", row.names = FALSE)
cat("saved: sarcosine_KO_verification.csv\n")
