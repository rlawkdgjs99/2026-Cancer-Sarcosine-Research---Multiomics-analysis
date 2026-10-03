#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(png))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not determine this script's path.")
SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", script_arg)))
BASE_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."))
RESULTS_DIR <- file.path(BASE_DIR, "results_integrated")
PROJECT_ROOT <- normalizePath(file.path(BASE_DIR, "..", "..", ".."))
OUT_DIR <- file.path(
  RESULTS_DIR, "CRC_WGS_POOLED_SARCOSINE_SPECIES_PATHWAY_TAGS_26.09.01"
)

input_file <- file.path(
  RESULTS_DIR, "Fig2d_SupFig8_taxonomic_complete_case_26.08.12",
  "sarcosine_bacteria_diff_abundance_pooled_corrected.csv"
)
values_file <- file.path(
  OUT_DIR, "CRC_WGS_pooled_sarcosine_associated_species_with_pathway_tags_values.csv"
)
png_file <- file.path(
  OUT_DIR, "CRC_WGS_pooled_sarcosine_associated_species_with_pathway_tags.png"
)
pdf_file <- file.path(
  OUT_DIR, "CRC_WGS_pooled_sarcosine_associated_species_with_pathway_tags.pdf"
)
transport_file <- file.path(
  PROJECT_ROOT, "Manuscript작업", "FigDesign&Manuscript",
  "PPT_INSERT_CRC_WGS_pooled_sarcosine_species_PATHWAY_TAGS.png"
)
report_file <- file.path(OUT_DIR, "verification_report.csv")
summary_file <- file.path(OUT_DIR, "verification_summary.txt")

checks <- data.frame(Check = character(), Pass = logical(), Detail = character(), stringsAsFactors = FALSE)
record <- function(name, pass, detail) {
  checks <<- rbind(checks, data.frame(Check = name, Pass = isTRUE(pass), Detail = detail, stringsAsFactors = FALSE))
}

record("All expected files exist",
       all(file.exists(c(input_file, values_file, png_file, pdf_file, transport_file))),
       paste(c(input_file, values_file, png_file, pdf_file, transport_file), collapse = " | "))

src <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
val <- read.csv(values_file, check.names = FALSE, stringsAsFactors = FALSE)
src <- src[order(src$log2FC), , drop = FALSE]

record("Exactly 20 unique species", nrow(val) == 20L && !anyDuplicated(val$Species),
       paste("n =", nrow(val), "; duplicates =", anyDuplicated(val$Species)))
record("Species identity/order unchanged", identical(val$Species, src$Species), paste(val$Species, collapse = " | "))
record("log2FC unchanged", identical(val$log2FC, src$log2FC), "exact numeric identity")
record("q values unchanged", identical(val$p_adj, src$p_adj), "exact numeric identity")
record("directions unchanged", identical(val$Direction, src$Direction), paste(val$Direction, collapse = " | "))
record("means unchanged",
       identical(val$Mean_Healthy, src$Mean_Healthy) && identical(val$Mean_Cancer, src$Mean_Cancer),
       "Mean_Healthy and Mean_Cancer exact")
record("max_abs_rho unchanged", identical(val$max_abs_rho, src$max_abs_rho), "exact numeric identity")
record("sample counts unchanged",
       identical(val$N_Healthy, src$N_Healthy) && identical(val$N_Cancer, src$N_Cancer) &&
         all(val$N_Healthy == 745L) && all(val$N_Cancer == 904L),
       "Healthy=745, Cancer=904")

role_map <- c(
  "Degradation" = "Deg",
  "Production" = "Prod",
  "Degradation & Production" = "Deg & Prod"
)
expected_role_label <- unname(role_map[src$Role])
record("Role field unchanged", identical(val$Role, src$Role), paste(unique(val$Role), collapse = " | "))
record("Role labels map exactly", identical(val$Role_label, expected_role_label), paste(table(val$Role_label), collapse = " | "))
record("All three tag classes present",
       setequal(unique(val$Role_label), c("Deg", "Prod", "Deg & Prod")),
       paste(sort(unique(val$Role_label)), collapse = ", "))

img <- png::readPNG(png_file, native = FALSE, info = TRUE)
record("PNG is 3450 x 2640 px", identical(as.integer(dim(img)[1:2]), c(2640L, 3450L)),
       paste("dimensions =", paste(dim(img), collapse = " x ")))
record("PNG is RGB without alpha", length(dim(img)) == 3L && dim(img)[3] == 3L,
       paste("channels =", ifelse(length(dim(img)) >= 3L, dim(img)[3], NA_integer_)))

sips_info <- system2(
  "/usr/bin/sips",
  args = c("-g", "format", "-g", "profile", "-g", "dpiWidth", "-g", "dpiHeight",
           "-g", "pixelWidth", "-g", "pixelHeight", shQuote(png_file)),
  stdout = TRUE, stderr = TRUE
)
record("PNG format reported as png", any(grepl("format: png", sips_info, fixed = TRUE)), paste(sips_info, collapse = " | "))
record("PNG embeds sRGB IEC61966-2.1",
       any(grepl("profile: sRGB IEC61966-2.1", sips_info, fixed = TRUE)),
       paste(sips_info, collapse = " | "))
record("PNG stores 300 dpi",
       any(grepl("dpiWidth: 300.000", sips_info, fixed = TRUE)) &&
         any(grepl("dpiHeight: 300.000", sips_info, fixed = TRUE)),
       paste(sips_info, collapse = " | "))
record("PowerPoint transport is byte-identical",
       identical(unname(tools::md5sum(png_file)), unname(tools::md5sum(transport_file))),
       paste("analysis =", unname(tools::md5sum(png_file)), "; transport =", unname(tools::md5sum(transport_file))))

flags <- system2("/bin/ls", args = c("-lO", shQuote(transport_file)), stdout = TRUE, stderr = TRUE)
record("PowerPoint transport is Finder-visible", !any(grepl("hidden", flags, fixed = TRUE)), paste(flags, collapse = " | "))

write.csv(checks, report_file, row.names = FALSE, quote = TRUE)
writeLines(
  c(
    sprintf("Pooled sarcosine-species pathway-tag verification: %d/%d checks passed", sum(checks$Pass), nrow(checks)),
    if (all(checks$Pass)) "FINAL STATUS: PASS" else "FINAL STATUS: FAIL",
    paste(checks$Check, ifelse(checks$Pass, "PASS", "FAIL"), sep = ": ")
  ),
  summary_file
)
if (!all(checks$Pass)) {
  print(checks[!checks$Pass, , drop = FALSE])
  stop("Verification failed. See: ", report_file)
}
cat(sprintf("Verification passed: %d/%d checks\n", sum(checks$Pass), nrow(checks)))
