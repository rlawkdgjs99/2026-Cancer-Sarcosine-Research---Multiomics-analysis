#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(png)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not determine this script's path.")
SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", script_arg)))
BASE_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."))
RESULTS_DIR <- file.path(BASE_DIR, "results_integrated")
PROJECT_ROOT <- normalizePath(file.path(BASE_DIR, "..", "..", ".."))
OUT_DIR <- file.path(RESULTS_DIR, "CRC_WGS_SARCOSINE_7KO_LOLLIPOP_26.09.01")

input_file <- file.path(RESULTS_DIR, "sarcosine", "sarcosine_KO_comparison_filtered_pooled.csv")
genome_file <- file.path(RESULTS_DIR, "kegg", "diff_KO_abundance_pooled.csv")
values_file <- file.path(OUT_DIR, "CRC_WGS_sarcosine_7KO_lollipop_values.csv")
png_file <- file.path(OUT_DIR, "CRC_WGS_sarcosine_7KO_lollipop.png")
pdf_file <- file.path(OUT_DIR, "CRC_WGS_sarcosine_7KO_lollipop.pdf")
transport_file <- file.path(
  PROJECT_ROOT, "Manuscript작업", "FigDesign&Manuscript",
  "PPT_INSERT_CRC_WGS_sarcosine_7KO_lollipop.png"
)
report_file <- file.path(OUT_DIR, "verification_report.csv")
summary_file <- file.path(OUT_DIR, "verification_summary.txt")

checks <- data.frame(Check = character(), Pass = logical(), Detail = character(), stringsAsFactors = FALSE)
record <- function(name, pass, detail) {
  checks <<- rbind(checks, data.frame(Check = name, Pass = isTRUE(pass), Detail = detail, stringsAsFactors = FALSE))
}

record("All expected files exist", all(file.exists(c(input_file, genome_file, values_file, png_file, pdf_file, transport_file))),
       paste(c(input_file, genome_file, values_file, png_file, pdf_file, transport_file), collapse = " | "))

src <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
genome <- read.csv(genome_file, check.names = FALSE, stringsAsFactors = FALSE)
val <- read.csv(values_file, check.names = FALSE, stringsAsFactors = FALSE)
expected_ko <- c("K00301", "K00302", "K00303", "K00305", "K00306", "K00315", "K08688")

record("Seven KO rows retained", nrow(val) == 7L, paste("n =", nrow(val)))
record("KO identity and order", identical(val$KO, expected_ko), paste(val$KO, collapse = ", "))
record("No duplicated KO", !anyDuplicated(val$KO), paste("duplicates =", anyDuplicated(val$KO)))
record("Role mapping", identical(val$Role, c(rep("Degradation", 5), rep("Production", 2))), paste(val$Role, collapse = ", "))

src <- src[match(expected_ko, src$KO), , drop = FALSE]
record("Frozen means copied exactly",
       identical(val$Mean_Healthy, src$Mean_Healthy) && identical(val$Mean_Cancer, src$Mean_Cancer),
       "Mean_Healthy and Mean_Cancer exact")
record("Frozen Wilcoxon P copied exactly", identical(val$wilcox_p, src$p_value), "p_value exact")
record("Frozen targeted BH q copied exactly", identical(val$targeted_7KO_BH_q, src$p_adj), "p_adj exact")
record("Stored q reproduces across seven KOs",
       max(abs(val$targeted_7KO_BH_q - p.adjust(val$wilcox_p, method = "BH"))) < 1e-14,
       "max absolute difference < 1e-14")

expected_fc <- log2((src$Mean_Cancer + 1e-8) / (src$Mean_Healthy + 1e-8))
record("log2FC reproduces from frozen means",
       max(abs(val$log2FC_Cancer_vs_Healthy - expected_fc)) < 1e-14,
       sprintf("max absolute difference = %.3g", max(abs(val$log2FC_Cancer_vs_Healthy - expected_fc))))
record("Prevalence copied exactly",
       identical(val$Prevalence_percent, src$Prev_Overall) &&
         identical(val$Prevalence_Healthy_percent, src$Prev_Healthy) &&
         identical(val$Prevalence_Cancer_percent, src$Prev_Cancer),
       "overall/Healthy/Cancer prevalence exact")
record("Saved prevalence flag copied exactly", identical(val$Pass_10pct_prevalence, src$Pass_Prevalence),
       paste(val$KO[val$Pass_10pct_prevalence], collapse = ", "))
record("Exactly three KOs pass >=10% prevalence", sum(val$Pass_10pct_prevalence) == 3L,
       paste("pass =", sum(val$Pass_10pct_prevalence), "; fail =", sum(!val$Pass_10pct_prevalence)))

qualified <- val$KO[val$Pass_10pct_prevalence]
g <- genome[match(qualified, genome$KO), , drop = FALSE]
v <- val[match(qualified, val$KO), , drop = FALSE]
record("Qualified KO values match genome-wide table",
       all(!is.na(g$KO)) &&
         max(abs(g$Mean_Healthy - v$Mean_Healthy)) < 1e-15 &&
         max(abs(g$Mean_Cancer - v$Mean_Cancer)) < 1e-15 &&
         max(abs(g$p_value - v$wilcox_p)) < 1e-15 &&
         max(abs(g$log2FC - v$log2FC_Cancer_vs_Healthy)) < 1e-14,
       paste(qualified, collapse = ", "))
record("Low-prevalence KOs absent from genome-wide filtered table",
       !any(val$KO[!val$Pass_10pct_prevalence] %in% genome$KO),
       paste(val$KO[!val$Pass_10pct_prevalence], collapse = ", "))

expected_status <- c("Not significant", "Not significant", "Healthy enriched",
                     "Cancer enriched", "Not significant", "Not significant", "Cancer enriched")
record("Direction/significance status correct", identical(val$Status, expected_status), paste(val$Status, collapse = " | "))

img <- png::readPNG(png_file, native = FALSE, info = TRUE)
record("PNG is 2700 x 2670 px", identical(as.integer(dim(img)[1:2]), c(2670L, 2700L)),
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
record("PNG embeds sRGB IEC61966-2.1", any(grepl("profile: sRGB IEC61966-2.1", sips_info, fixed = TRUE)), paste(sips_info, collapse = " | "))
record("PNG stores 600 dpi", any(grepl("dpiWidth: 600.000", sips_info, fixed = TRUE)) &&
         any(grepl("dpiHeight: 600.000", sips_info, fixed = TRUE)), paste(sips_info, collapse = " | "))
record("PowerPoint transport is byte-identical",
       identical(unname(tools::md5sum(png_file)), unname(tools::md5sum(transport_file))),
       paste("analysis =", unname(tools::md5sum(png_file)), "; transport =", unname(tools::md5sum(transport_file))))

flags <- system2("/bin/ls", args = c("-lO", shQuote(transport_file)), stdout = TRUE, stderr = TRUE)
record("PowerPoint transport is Finder-visible", !any(grepl("hidden", flags, fixed = TRUE)), paste(flags, collapse = " | "))

write.csv(checks, report_file, row.names = FALSE, quote = TRUE)
writeLines(
  c(
    sprintf("CRC WGS seven-KO lollipop verification: %d/%d checks passed", sum(checks$Pass), nrow(checks)),
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
