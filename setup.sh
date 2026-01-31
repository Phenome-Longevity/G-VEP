#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${ROOT}/data"

# R2 Configuration (read-only)
R2_ENDPOINT="https://cbcc869912f127e202bc081ff4604ff4.r2.cloudflarestorage.com"
R2_ACCESS_KEY="e9d8994654913022e73b4a66ae553fb5"
R2_SECRET_KEY="dae77a24a08761649680ed2adc8425af7f7479908039bb6cec5a0a70c446e65d"
R2_BUCKET="g-vep"

info() { echo -e "\033[1;34m[INFO]\033[0m $1"; }
ok() { echo -e "\033[1;32m[OK]\033[0m $1"; }
err() { echo -e "\033[1;31m[ERROR]\033[0m $1"; exit 1; }

# Install rclone if needed
install_rclone() {
    if command -v rclone &>/dev/null; then
        return 0
    fi
    info "Installing rclone..."
    curl -fsSL https://rclone.org/install.sh | sudo bash
}

# Configure rclone for R2
configure_rclone() {
    mkdir -p ~/.config/rclone
    cat > ~/.config/rclone/rclone.conf << EOF
[g-vep-r2]
type = s3
provider = Cloudflare
access_key_id = ${R2_ACCESS_KEY}
secret_access_key = ${R2_SECRET_KEY}
endpoint = ${R2_ENDPOINT}
no_check_bucket = true
EOF
}

# Download with rclone (handles chunking automatically)
download() {
    local src="$1" dest="$2"
    [[ -f "$dest" ]] && return 0
    info "Downloading $(basename "$dest")..."
    rclone copyto "g-vep-r2:${R2_BUCKET}/${src}" "$dest" --progress
}

download_dir() {
    local src="$1" dest="$2"
    [[ -d "$dest" ]] && [[ "$(ls -A "$dest" 2>/dev/null)" ]] && return 0
    info "Downloading ${src}..."
    mkdir -p "$dest"
    rclone copy "g-vep-r2:${R2_BUCKET}/${src}" "$dest" --progress
}

echo ""
echo "G-VEP Setup"
echo "==========="

# Setup rclone
install_rclone
configure_rclone

# GPU indices
info "Checking GPU indices..."
mkdir -p "$DATA/gpu_indices"
for db in alphamissense clinvar dbnsfp dbscsnv revel spliceai; do
    for ext in keys.bin vals.bin meta.json; do
        download "gpu_indices/${db}.${ext}" "$DATA/gpu_indices/${db}.${ext}"
    done
done
ok "GPU indices ready"

# VEP cache
info "Checking VEP cache..."
mkdir -p "$DATA/vep_cache"
if [[ ! -d "$DATA/vep_cache/homo_sapiens_refseq/113_GRCh38" ]]; then
    download "vep_cache/homo_sapiens_refseq_vep_113_GRCh38.tar.gz" "$DATA/vep_cache/vep_cache.tar.gz"
    info "Extracting VEP cache..."
    tar -xzf "$DATA/vep_cache/vep_cache.tar.gz" -C "$DATA/vep_cache"
    rm -f "$DATA/vep_cache/vep_cache.tar.gz"
fi
ok "VEP cache ready"

# FASTA
info "Checking FASTA..."
download "fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz" "$DATA/vep_cache/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
download "fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.fai" "$DATA/vep_cache/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.fai"
download "fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.gzi" "$DATA/vep_cache/Homo_sapiens.GRCh38.dna.toplevel.fa.gz.gzi"
ok "FASTA ready"

# LOFTEE
info "Checking LOFTEE..."
mkdir -p "$DATA/vep_cache/VEP_plugins/loftee-1.0.4_GRCh38"
download "loftee/LoF.pm" "$DATA/vep_cache/VEP_plugins/loftee-1.0.4_GRCh38/LoF.pm"
download "loftee/human_ancestor.fa.gz" "$DATA/vep_cache/VEP_plugins/loftee-1.0.4_GRCh38/human_ancestor.fa.gz"
download "loftee/human_ancestor.fa.gz.fai" "$DATA/vep_cache/VEP_plugins/loftee-1.0.4_GRCh38/human_ancestor.fa.gz.fai"
download "loftee/human_ancestor.fa.gz.gzi" "$DATA/vep_cache/VEP_plugins/loftee-1.0.4_GRCh38/human_ancestor.fa.gz.gzi"
download "loftee/gerp_conservation_scores.homo_sapiens.GRCh38.bw" "$DATA/vep_cache/VEP_plugins/loftee-1.0.4_GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
download "loftee/loftee.sql" "$DATA/vep_cache/VEP_plugins/loftee-1.0.4_GRCh38/loftee.sql"
ok "LOFTEE ready"

# Custom VCFs
info "Checking custom VCFs..."
mkdir -p "$DATA/vep_cache/vep_custom"
download "vep_custom/clinvar_20240611_PLPC_new_CPLP.vcf.gz" "$DATA/vep_cache/vep_custom/clinvar_20240611_PLPC_new_CPLP.vcf.gz"
download "vep_custom/clinvar_20240611_PLPC_new_CPLP.vcf.gz.tbi" "$DATA/vep_cache/vep_custom/clinvar_20240611_PLPC_new_CPLP.vcf.gz.tbi"
download "vep_custom/Whole_out_sorted_2024.vcf.gz" "$DATA/vep_cache/vep_custom/Whole_out_sorted_2024.vcf.gz"
download "vep_custom/Whole_out_sorted_2024.vcf.gz.tbi" "$DATA/vep_cache/vep_custom/Whole_out_sorted_2024.vcf.gz.tbi"
ok "Custom VCFs ready"

echo ""
ok "Setup complete. Run: ./run_pipeline.sh <input.vcf.gz>"