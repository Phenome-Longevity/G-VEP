#!/bin/bash

set -e

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VEP_CACHE="${SCRIPT_DIR}/data/vep_cache"
GPU_INDEX_DIR="${SCRIPT_DIR}/data/gpu_indices"
VEP_IMAGE="ensemblorg/ensembl-vep:release_113.0"
THREADS=10

# Output directory (separate from original pipeline)
OUTPUT_BASE="${SCRIPT_DIR}/outputs"

# Parse arguments
INPUT_VCF="$1"

if [[ -z "$INPUT_VCF" ]]; then
    echo "Usage: $0 <input.vcf.gz>"
    echo ""
    echo "Output will be written to: ${OUTPUT_BASE}/"
    exit 1
fi

# Get paths
INPUT_VCF=$(realpath "$INPUT_VCF")
SAMPLE_NAME=$(basename "${INPUT_VCF}" .vcf.gz)
SAMPLE_NAME="${SAMPLE_NAME%%.hard-filtered}"

# Create output directories (mirrors original structure)
mkdir -p "${OUTPUT_BASE}/logs"
mkdir -p "${OUTPUT_BASE}/vcf_pass"
mkdir -p "${OUTPUT_BASE}/vep"
mkdir -p "${OUTPUT_BASE}/gpu_annotations"

# Output files
LOG_FILE="${OUTPUT_BASE}/logs/${SAMPLE_NAME}_gpu_log.txt"
PASS_VCF="${OUTPUT_BASE}/vcf_pass/${SAMPLE_NAME}_PASS.vcf.gz"
VEP_OUTPUT="${OUTPUT_BASE}/vep/${SAMPLE_NAME}_vep_annotated.vcf.gz"
GPU_TSV="${OUTPUT_BASE}/gpu_annotations/${SAMPLE_NAME}_gpu.tsv"
FINAL_OUTPUT="${OUTPUT_BASE}/vep/${SAMPLE_NAME}_vep_annotated_final.vcf.gz"

# Log everything
exec > >(tee -a "$LOG_FILE") 2>&1

echo "$(date)"
echo "GPU-Accelerated VEP Pipeline v2"
echo "============================================================"
echo "Input:      $INPUT_VCF"
echo "Sample:     $SAMPLE_NAME"
echo "Output dir: $OUTPUT_BASE"
echo "VEP Cache:  $VEP_CACHE"
echo "GPU Index:  $GPU_INDEX_DIR"
echo "============================================================"
echo ""

START_TIME=$(date +%s)

# ============================================================
# Step 1: PASS Filter (same as original)
# ============================================================
echo "[Step 1/4] Filtering PASS variants..."
echo "TIMER|PASS_FILTER_START|$(date +%s)"
STEP1_START=$(date +%s)

zgrep -E "^#|PASS" "$INPUT_VCF" | bgzip > "$PASS_VCF"

STEP1_END=$(date +%s)
PASS_COUNT=$(zgrep -v "^#" "$PASS_VCF" | wc -l)
echo "TIMER|PASS_FILTER_END|$(date +%s)"
echo "VARIANT_COUNT|PASS|$PASS_COUNT"
echo "  PASS variants: ${PASS_COUNT}"
echo "  Done in $((STEP1_END - STEP1_START))s"
echo ""

# ============================================================
# Step 2: VEP Annotation (lightweight - no slow plugins)
# ============================================================
echo "[Step 2/4] Running VEP (lightweight mode)..."
echo "TIMER|VEP_START|$(date +%s)"
STEP2_START=$(date +%s)

docker run --rm \
    -v "${VEP_CACHE}:/data/vep_cache:ro" \
    -v "${OUTPUT_BASE}:/data/output" \
    ${VEP_IMAGE} \
    vep \
    --cache --dir_cache /data/vep_cache \
    --offline \
    --cache_version 113 \
    --fork ${THREADS} \
    --format vcf \
    --vcf \
    --compress_output bgzip \
    --assembly GRCh38 \
    --force_overwrite \
    -i /data/output/vcf_pass/${SAMPLE_NAME}_PASS.vcf.gz \
    -o /data/output/vep/${SAMPLE_NAME}_vep_annotated.vcf.gz \
    --symbol --check_existing --variant_class \
    --sift b --polyphen b \
    --synonyms /data/vep_cache/homo_sapiens_refseq/113_GRCh38/chr_synonyms.txt \
    --hgvs --refseq \
    --fasta /data/vep_cache/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
    --canonical \
    --pick --pick_order mane_select,rank \
    --af --af_gnomade --af_gnomadg --max_af \
    --custom /data/vep_cache/vep_custom/clinvar_20240611_PLPC_new_CPLP.vcf.gz,ClinVar,vcf,exact,0,ID,CLNSIG,CLNDN,CLNHGVS,CLNSIGINCL,CLNVC,GENEINFO,CLNDISDB,CLNSIGCONF,CLNREVSTAT,CLNDNINCL \
    --custom /data/vep_cache/vep_custom/Whole_out_sorted_2024.vcf.gz,Database,vcf,exact,0,Type,SZAID \
    --dir_plugins /data/vep_cache/VEP_plugins \
    --plugin LoF,loftee_path:/data/vep_cache/VEP_plugins/loftee-1.0.4_GRCh38,human_ancestor_fa:/data/vep_cache/VEP_plugins/loftee-1.0.4_GRCh38/human_ancestor.fa.gz,conservation_file:/data/vep_cache/VEP_plugins/loftee-1.0.4_GRCh38/loftee.sql,gerp_bigwig:/data/vep_cache/VEP_plugins/loftee-1.0.4_GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
    --verbose

STEP2_END=$(date +%s)
echo "TIMER|VEP_END|$(date +%s)"
VEP_COUNT=$(zgrep -v "^#" "$VEP_OUTPUT" | wc -l)
echo "VARIANT_COUNT|VEP_OUTPUT|$VEP_COUNT"
echo "  VEP completed in $((STEP2_END - STEP2_START))s"
echo ""

# ============================================================
# Step 3: GPU Annotation Lookup
# ============================================================
echo "[Step 3/4] Running GPU annotation lookup..."
echo "TIMER|GPU_START|$(date +%s)"
STEP3_START=$(date +%s)

python3 "${SCRIPT_DIR}/src/gpu_annotate.py" \
    -i "$PASS_VCF" \
    -o "$GPU_TSV" \
    --index-dir "$GPU_INDEX_DIR"

STEP3_END=$(date +%s)
echo "TIMER|GPU_END|$(date +%s)"
echo "  GPU annotation completed in $((STEP3_END - STEP3_START))s"
echo ""

# ============================================================
# Step 4: Merge into CSQ (exact format match)
# ============================================================
echo "[Step 4/4] Merging annotations into CSQ..."
echo "TIMER|MERGE_START|$(date +%s)"
STEP4_START=$(date +%s)

python3 "${SCRIPT_DIR}/src/merge_annotations.py" \
    --vep-vcf "$VEP_OUTPUT" \
    --gpu-tsv "$GPU_TSV" \
    --output "$FINAL_OUTPUT"

STEP4_END=$(date +%s)
echo "TIMER|MERGE_END|$(date +%s)"
echo "  Merge completed in $((STEP4_END - STEP4_START))s"
echo ""



# ============================================================
# Summary
# ============================================================
END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

echo "============================================================"
echo "Pipeline Complete"
echo "============================================================"
echo "Step 1 (PASS):  $((STEP1_END - STEP1_START))s"
echo "Step 2 (VEP):   $((STEP2_END - STEP2_START))s"
echo "Step 3 (GPU):   $((STEP3_END - STEP3_START))s"
echo "Step 4 (Merge): $((STEP4_END - STEP4_START))s"
echo "------------------------------------------------------------"
echo "TIMER|TOTAL|$START_TIME|$END_TIME|$TOTAL_TIME seconds"