#!/usr/bin/env bash
# G-VEP setup: env (VEP 113) + compile tools + download all data from R2.
# Total download ~226 GB (VEP cache 48 + plugins 129 + index 49 + FASTA/custom).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GVEP_DATA="${GVEP_DATA:-$HOME}"
VC="$HOME/gvep/data/vep_cache"          # VEP cache + custom (--dir_cache)
DC="$ROOT/data/vep_cache"              # FASTA + plugins (--fasta, --plugin)

# 1. environment
echo "[setup] micromamba env (VEP 113 + bcftools/htslib/pysam/rclone)..."
"$HOME/bin/micromamba" create -y -p "$HOME/mamba/gvep" -c conda-forge -c bioconda \
  ensembl-vep=113 bcftools htslib pysam rclone || true

# 2. compile index tools -> bin/
echo "[setup] compile tools -> bin/"
HTS="$HOME/mamba/gvep"; mkdir -p "$ROOT/bin"
gcc -O3 -o "$ROOT/bin/build_index" "$ROOT/src/build_index.c"
gcc -O3 -I"$HTS/include" -o "$ROOT/bin/lookup" "$ROOT/src/lookup.c" -L"$HTS/lib" -lhts -Wl,-rpath,"$HTS/lib"

# 3. rclone (read-only g-vep bucket)
mkdir -p "$HOME/.config/rclone" "$VC/vep_custom" "$DC" "$GVEP_DATA/gvep_index"
cat > "$HOME/.config/rclone/gvep.conf" <<'RC'
[r2]
type = s3
provider = Cloudflare
access_key_id = e9d8994654913022e73b4a66ae553fb5
secret_access_key = dae77a24a08761649680ed2adc8425af7f7479908039bb6cec5a0a70c446e65d
endpoint = https://cbcc869912f127e202bc081ff4604ff4.r2.cloudflarestorage.com
no_check_bucket = true
RC
RC="$HOME/mamba/gvep/bin/rclone --config $HOME/.config/rclone/gvep.conf -P --transfers 4 --s3-chunk-size 64M --s3-upload-concurrency 8"

# 4. VEP cache (-> VC) + ClinVar/SZA custom
echo "[data] VEP cache (48GB, extracting)..."
$RC copyto r2:g-vep/vep_cache/homo_sapiens_refseq_vep_113_GRCh38.tar.gz /tmp/_vepcache.tar.gz
tar -xaf /tmp/_vepcache.tar.gz -C "$VC" && rm -f /tmp/_vepcache.tar.gz
for f in clinvar_20240611_PLPC_new_CPLP.vcf.gz{,.tbi} Whole_out_sorted_2024.vcf.gz{,.tbi}; do
  $RC copyto "r2:g-vep/vep_custom/$f" "$VC/vep_custom/$f"; done

# 5. FASTA (-> DC and VC) + plugins (-> DC)
echo "[data] FASTA..."
for f in Homo_sapiens.GRCh38.dna.toplevel.fa.gz{,.fai,.gzi}; do
  $RC copyto "r2:g-vep/fasta/$f" "$DC/$f"; cp "$DC/$f" "$VC/$f"; done
echo "[data] plugins (129GB)..."
$RC copy r2:g-vep/plugins/ "$DC/"

# 6. index (-> GVEP_DATA/gvep_index)
echo "[data] index (49GB)..."
for f in cache.bgz cache.bgz.gzi h_sorted.bin v_sorted.bin; do
  $RC copyto "r2:g-vep/gvep_index/$f" "$GVEP_DATA/gvep_index/$f"; done

echo "[setup] DONE. Run:  ./run_gvep.sh <input.vcf.gz> [threads]"
