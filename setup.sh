#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${ROOT}/data"
R2="https://data.szags.uk"

info() { echo -e "\033[1;34m[INFO]\033[0m $1"; }
ok() { echo -e "\033[1;32m[OK]\033[0m $1"; }
err() { echo -e "\033[1;31m[ERROR]\033[0m $1"; exit 1; }

download() {
    local url="$1" dest="$2"
    [[ -f "$dest" ]] && return 0
    info "Downloading $(basename "$dest")..."
    curl -fL --progress-bar -o "$dest" "$url"
}

echo ""
echo "G-VEP Setup"
echo "==========="

# GPU indices
info "Checking GPU indices..."
mkdir -p "$DATA/gpu_indices"
for db in alphamissense clinvar dbnsfp dbscsnv revel spliceai; do
    for ext in keys.bin vals.bin meta.json; do
        download "$R2/gpu_indices/${db}.${ext}" "$DATA/gpu_indices/${db}.${ext}"
    done
done
ok "GPU indices ready"

# VEP cache
info "Checking VEP cache..."
mkdir -p "$DATA/vep_cache"
if [[ ! -d "$DATA/vep_cache/homo_sapiens_refseq/113_GRCh38" ]]; then
    download "$R2/vep_cache/homo_sapiens_refseq_vep_113_GRCh38.tar.gz" "$DATA/vep_cache/vep_cache.tar.gz"
    info "Extracting VEP cache..."
    tar -xzf "$DATA/vep_cache/vep_cache.tar.gz" -C "$DATA/vep_cache"
    rm -f "$DATA/vep_cache/vep_cache.tar.gz"
fi
ok "VEP cache ready"

# FASTA
info "Checking FASTA..."
download "$R2/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz" "$DATA/vep_cache/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
download "$R2/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.fai" "$DATA/vep_cache/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.fai"
download "$R2/fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.gzi" "$DATA/vep_cache/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.gzi"
ok "FASTA ready"

# LOFTEE
info "Checking LOFTEE..."
mkdir -p "$DATA/vep_cache/VEP_plugins/loftee-1.0.4_GRCh38"
download "$R2/loftee/LoF.pm" "$DATA/vep_cache/VEP_plugins/loftee-1.0.4_GRCh38/LoF.pm"
download "$R2/loftee/human_ancestor.fa.gz" "$DATA/vep_cache/VEP_plugins/loftee-1.0.4_GRCh38/human_ancestor.fa.gz"
download "$R2/loftee/human_ancestor.fa.gz.fai" "$DATA/vep_cache/VEP_plugins/loftee-1.0.4_GRCh38/human_ancestor.fa.gz.fai"
download "$R2/loftee/human_ancestor.fa.gz.gzi" "$DATA/vep_cache/VEP_plugins/loftee-1.0.4_GRCh38/human_ancestor.fa.gz.gzi"
download "$R2/loftee/gerp_conservation_scores.homo_sapiens.GRCh38.bw" "$DATA/vep_cache/VEP_plugins/loftee-1.0.4_GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
download "$R2/loftee/loftee.sql" "$DATA/vep_cache/VEP_plugins/loftee-1.0.4_GRCh38/loftee.sql"
ok "LOFTEE ready"

# Custom VCFs
info "Checking custom VCFs..."
mkdir -p "$DATA/vep_cache/vep_custom"
download "$R2/vep_custom/clinvar_20240611_PLPC_new_CPLP.vcf.gz" "$DATA/vep_cache/vep_custom/clinvar_20240611_PLPC_new_CPLP.vcf.gz"
download "$R2/vep_custom/clinvar_20240611_PLPC_new_CPLP.vcf.gz.tbi" "$DATA/vep_cache/vep_custom/clinvar_20240611_PLPC_new_CPLP.vcf.gz.tbi"
download "$R2/vep_custom/Whole_out_sorted_2024.vcf.gz" "$DATA/vep_cache/vep_custom/Whole_out_sorted_2024.vcf.gz"
download "$R2/vep_custom/Whole_out_sorted_2024.vcf.gz.tbi" "$DATA/vep_cache/vep_custom/Whole_out_sorted_2024.vcf.gz.tbi"
ok "Custom VCFs ready"

echo ""
ok "All set. Run: ./run_pipeline.sh <input.vcf.gz>"
