#!/usr/bin/env Rscript
# ============================================================
#  02_dada2_pe.R  —  PE 16S DADA2 워크플로 (cutadapt 출력 입력)
#  입력:  $WORK/cutadapt/{run}_1.fastq.gz, _2.fastq.gz (프라이머 제거됨)
#  출력:  $WORK/  seqtab_nochim.rds, asv_table_with_taxonomy.tsv,
#               genus_table.tsv, tracking.tsv, quality_profiles.pdf,
#               error_profiles.pdf, taxonomy.rds, sessionInfo.txt
#  참조:  Silva v138.1 trainset ($BASE/ref/...)
# ============================================================
suppressMessages(library(dada2))
set.seed(100)

WORK    <- Sys.getenv("WORK");  BASE <- Sys.getenv("BASE")
THREADS <- as.integer(Sys.getenv("THREADS", "8"))
CUT  <- file.path(WORK, "cutadapt")
FILT <- file.path(WORK, "filtered")
OUT  <- WORK
REF  <- file.path(BASE, "ref", "silva_nr99_v138.1_train_set.fa.gz")
dir.create(FILT, showWarnings = FALSE, recursive = TRUE)

# ---- 파라미터 (조정 가능 — quality_profiles.pdf 와 tracking.tsv 보고 튜닝) ----
# cutadapt 후 read 길이가 가변(스페이서 길이차)이라 고정 절단은 read 손실 위험 →
# truncLen=0(절단 안 함)으로 V3-V4 머지 오버랩을 보존하고, 품질은 maxEE/truncQ로 거른다.
TRUNCLEN <- c(0, 0)
MAXEE    <- c(2, 5)     # R2(역방향)는 품질이 낮아 더 관대하게
TRUNCQ   <- 2
MINLEN   <- 50

cat("== DADA2 시작 ==  WORK=", WORK, " THREADS=", THREADS, "\n", sep = "")
stopifnot(file.exists(REF))

fnFs <- sort(list.files(CUT, pattern = "_1\\.fastq\\.gz$", full.names = TRUE))
fnRs <- sort(list.files(CUT, pattern = "_2\\.fastq\\.gz$", full.names = TRUE))
sample.names <- sub("_1\\.fastq\\.gz$", "", basename(fnFs))
cat("입력 샘플 수:", length(fnFs), "\n")
stopifnot(length(fnFs) > 0, length(fnFs) == length(fnRs))

# 품질 프로파일 (앞 4개 샘플)
k <- min(4, length(fnFs))
pdf(file.path(OUT, "quality_profiles.pdf"), width = 9, height = 6)
print(plotQualityProfile(fnFs[1:k]) + ggplot2::ggtitle("Forward (R1) after cutadapt"))
print(plotQualityProfile(fnRs[1:k]) + ggplot2::ggtitle("Reverse (R2) after cutadapt"))
dev.off()

filtFs <- file.path(FILT, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(FILT, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names; names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = TRUNCLEN, maxEE = MAXEE, truncQ = TRUNCQ,
                     maxN = 0, rm.phix = TRUE, minLen = MINLEN,
                     compress = TRUE, multithread = THREADS)
saveRDS(out, file.path(OUT, "filter_out.rds"))

# 필터 후 read가 남은 샘플만 유지
ok <- file.exists(filtFs) & file.exists(filtRs)
cat("필터 통과 샘플:", sum(ok), "/", length(ok), "\n")
filtFs <- filtFs[ok]; filtRs <- filtRs[ok]; sample.names <- sample.names[ok]

errF <- learnErrors(filtFs, multithread = THREADS)
errR <- learnErrors(filtRs, multithread = THREADS)
pdf(file.path(OUT, "error_profiles.pdf"), width = 9, height = 6)
print(plotErrors(errF, nominalQ = TRUE)); dev.off()

dadaFs  <- dada(filtFs, err = errF, multithread = THREADS)
dadaRs  <- dada(filtRs, err = errR, multithread = THREADS)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
seqtab  <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = THREADS, verbose = TRUE)
saveRDS(seqtab.nochim, file.path(OUT, "seqtab_nochim.rds"))
cat("ASV 수:", ncol(seqtab.nochim), " / 샘플:", nrow(seqtab.nochim), "\n")
cat("머지 후 chimera 제거 비율(read):",
    round(sum(seqtab.nochim) / sum(seqtab), 4), "\n")

# ---- read 추적표 ----
getN <- function(x) sum(getUniques(x))
asV  <- function(z) if (is.list(z)) sapply(z, getN) else getN(z)
track <- cbind(out[ok, , drop = FALSE], asV(dadaFs), asV(dadaRs),
               asV(mergers), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.table(track, file.path(OUT, "tracking.tsv"), sep = "\t", quote = FALSE, col.names = NA)

# ---- 분류 (Silva v138.1) ----
taxa <- assignTaxonomy(seqtab.nochim, REF, multithread = THREADS, tryRC = TRUE)
saveRDS(taxa, file.path(OUT, "taxonomy.rds"))

# ---- ASV 출력 ----
asv.seqs <- colnames(seqtab.nochim)
asv.ids  <- paste0("ASV", seq_along(asv.seqs))
writeLines(paste0(">", asv.ids, "\n", asv.seqs), file.path(OUT, "asv_seqs.fasta"))
at <- t(seqtab.nochim); rownames(at) <- asv.ids
write.table(data.frame(ASV = asv.ids, as.data.frame(taxa), at, check.names = FALSE),
            file.path(OUT, "asv_table_with_taxonomy.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- genus 수준 집계 ----
tax.df <- as.data.frame(taxa, stringsAsFactors = FALSE)
genus.label <- ifelse(is.na(tax.df$Genus),
                      ifelse(is.na(tax.df$Family), "Unclassified",
                             paste0("Unclassified_", tax.df$Family)),
                      tax.df$Genus)
genus.tab <- rowsum(t(seqtab.nochim), group = genus.label)   # genus x sample
write.table(data.frame(genus = rownames(genus.tab), genus.tab, check.names = FALSE),
            file.path(OUT, "genus_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("genus 수:", nrow(genus.tab), "\n")

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("DADA2_PIPELINE_DONE\n")
