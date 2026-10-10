#!/usr/bin/env bash
# Align paired reads with bwa mem, add a read group, convert to BAM, sort.
#
# Usage:
#   bash bwa.sh R1.(fastq|fq).gz R2.(fastq|fq).gz REF.fasta THREADS OUT.bam TMP_PREFIX LOG
#
# The reference must already be indexed:  bwa index REF.fasta

set -o pipefail

if [ $# -ne 7 ]; then
    echo "Need 7 arguments, got $#." >&2
    echo "Usage: bash bwa.sh R1 R2 REF THREADS OUT_BAM TMP_PREFIX LOG" >&2
    exit 1
fi

R1=$1
R2=$2
REF=$3
THREADS=$4
OUT=$5
TMP=$6
LOG=$7

[ -f "$R1" ]       || { echo "R1 not found: $R1" >&2; exit 1; }
[ -f "$R2" ]       || { echo "R2 not found: $R2" >&2; exit 1; }
[ -f "$REF" ]      || { echo "Reference not found: $REF" >&2; exit 1; }
[ -f "${REF}.bwt" ] || { echo "No bwa index for $REF. Run: bwa index $REF" >&2; exit 1; }

prepare_readgroup() {
    local fq=$1
    case "$fq" in
        *.fastq.gz|*.fq.gz) ;;
        *) echo "Sequence read file extension not recognized: $fq" >&2; return 1 ;;
    esac

    # Sample name: file name without _R1... or extension
    local samplename
    samplename=$(basename "$fq")
    samplename=${samplename%%_R1*}
    samplename=${samplename%.fastq.gz}
    samplename=${samplename%.fq.gz}

    # Read group ID from the first read header (same logic as the original)
    local firstLine id_name
    firstLine=$(gzip -dc "$fq" 2>/dev/null | head -n 1)
    if [[ $firstLine == *":"* ]]; then
        id_name=$(echo "$firstLine" | awk -F':' '{print $2}' | awk '{print $1}')
    elif [[ $firstLine == *"/"* ]]; then
        id_name=$(echo "$firstLine" | awk -F'/' '{print $2}' | awk '{print $1}')
    else
        id_name="1"
    fi
    [ -n "$id_name" ] || id_name="1"

    echo "@RG\\tID:${id_name}\\tSM:${samplename}\\tLB:1\\tPL:Illumina"
}

split_field=$(prepare_readgroup "$R1") || exit 1
echo "split_field: $split_field"
echo "R1: $R1"
echo "R2: $R2"
echo "REF: $REF"
echo "THREADS: $THREADS"
echo "OUT: $OUT"
echo "TMP: $TMP"
echo "LOG: $LOG"

{
    bwa mem -M -R "$split_field" -t "$THREADS" "$REF" "$R1" "$R2" \
        | samtools view -Sb - \
        | samtools sort -m 500M -@ 0 -o "$OUT" -T "$TMP" -
} &> "$LOG"
status=$?

if [ $status -ne 0 ]; then
    echo "Alignment failed (exit $status). See $LOG" >&2
    exit $status
fi
echo "Done: $OUT"