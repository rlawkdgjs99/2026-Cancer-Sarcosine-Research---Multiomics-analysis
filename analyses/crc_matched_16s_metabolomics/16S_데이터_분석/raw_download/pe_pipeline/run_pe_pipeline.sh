#!/usr/bin/env bash
# ============================================================
#  run_pe_pipeline.sh — PE-only 16S 파이프라인 (cutadapt -> DADA2)
#  - cutadapt: crc16s 환경의 cutadapt 로 V3-V4 프라이머(319F/806R) 제거
#  - DADA2:    시스템 R(dada2 1.40.0)로 denoise -> ASV -> Silva 분류 -> genus
#  env 변수(선택): SAMPLES, WORK, PARALLEL, THREADS  (테스트/전체 동일 코드)
# ============================================================
set -uo pipefail

BASE="${BASE:-16S_데이터_분석/raw_download}"   # relative to the analysis-folder root (working dir)
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
SAMPLES="${SAMPLES:-$BASE/pe_samples.tsv}"
WORK="${WORK:-$BASE/pe_pipeline}"
PARALLEL="${PARALLEL:-6}"
THREADS="${THREADS:-8}"
CUTADAPT="${CUTADAPT:-cutadapt}"   # path to the cutadapt executable
RSCRIPT="$(command -v Rscript)"
# 정방향: degenerate 341F. reads에서 경험적으로 확인(retention ~97.5%).
# 논문 표기 319F "ACTCCTACGGGAGGCAGCAG"(비축퇴)는 ~5%만 통과해 부적합 → 데이터 기준으로 수정.
FWD="CCTACGGGNGGCWGCAG"      # 341F (degenerate)
REV="GGACTACHVGGGTWTCTAAT"   # 806R (논문 표기와 동일, 이미 degenerate)
IN="$BASE/pe"
CUT="$WORK/cutadapt"
LOG="$WORK/logs/cutadapt"

ts(){ date '+%H:%M:%S'; }
mkdir -p "$CUT" "$LOG"
[ -x "$CUTADAPT" ] || { echo "!! cutadapt 실행불가: $CUTADAPT"; exit 1; }
[ -n "$RSCRIPT" ]  || { echo "!! Rscript 없음"; exit 1; }
[ -s "$SAMPLES" ]  || { echo "!! 샘플목록 없음: $SAMPLES"; exit 1; }

NS=$(wc -l < "$SAMPLES" | tr -d ' ')
echo "[$(ts)] ===== PE 파이프라인 시작 =====  samples=$NS  WORK=$WORK  parallel=$PARALLEL  threads=$THREADS"

# ---------- (1) cutadapt 프라이머 제거 ----------
echo "[$(ts)] (1) cutadapt: V3-V4 프라이머 제거 (319F/806R, --discard-untrimmed)"
export IN CUT LOG CUTADAPT FWD REV
: > "$CUT/.failures"
awk -F'\t' '{print $1}' "$SAMPLES" | xargs -P "$PARALLEL" -I{} bash -c '
  run="$1"
  in1="$IN/${run}_1.fastq.gz"; in2="$IN/${run}_2.fastq.gz"
  if [ ! -s "$in1" ] || [ ! -s "$in2" ]; then echo "MISSING $run" >> "$CUT/.failures"; exit 0; fi
  "$CUTADAPT" -g "$FWD" -G "$REV" --discard-untrimmed -m 50 -j 1 \
    -o "$CUT/${run}_1.fastq.gz" -p "$CUT/${run}_2.fastq.gz" \
    "$in1" "$in2" > "$LOG/${run}.log" 2>&1 || echo "FAIL $run" >> "$CUT/.failures"
' _ {}
NOUT=$(ls "$CUT"/*_1.fastq.gz 2>/dev/null | wc -l | tr -d ' ')
NFAIL=$(grep -c . "$CUT/.failures" 2>/dev/null || echo 0)
echo "[$(ts)] cutadapt 완료: 출력 $NOUT pairs, 실패/누락 $NFAIL"
# 통과율(=프라이머 발견된 read쌍) 요약
awk '/Pairs written \(passing filters\)/{gsub(/[(),%]/,""); print $(NF)}' "$LOG"/*.log 2>/dev/null \
  | awk '{s+=$1;n++} END{if(n)printf "[cutadapt] 평균 프라이머-통과율: %.1f%% (n=%d)\n", s/n, n}'
[ "$NOUT" -gt 0 ] || { echo "!! cutadapt 출력 0 — 중단"; exit 2; }

# ---------- (2) DADA2 ----------
echo "[$(ts)] (2) DADA2 (시스템 R, dada2 1.40.0)"
WORK="$WORK" BASE="$BASE" THREADS="$THREADS" "$RSCRIPT" "$SCRIPTDIR/02_dada2_pe.R"
RC=$?
echo "[$(ts)] ===== 파이프라인 종료 (DADA2 exit=$RC) ====="
exit $RC
