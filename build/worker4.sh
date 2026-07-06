#!/bin/bash
# REBUILD worker: worker3 + (a) bcftools norm (left-align+trim -> canonical, matches normalized input)
#                          + (b) ALL plugins (dbNSFP/dbscSNV/SpliceAI/AlphaMissense) -> exact plugin scores baked in.
# Cache is now byte-identical to "VEP --plugin ... on normalized input". Resumable per-chrom.
# Usage: worker4.sh "<RUN>" <CACHE> <NPAR> <FORK> <OUT> <chr...>
set -uo pipefail
RUN="$1"; CACHE="$2"; NPAR="$3"; FORK="$4"; OUT="$5"; shift 5; CHROMS="$*"
SRCDIR="$(cd "$(dirname "$0")" && pwd)"
export RUN CACHE FORK SRCDIR
export FASTA="$CACHE/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
export SYN="$CACHE/homo_sapiens_refseq/113_GRCh38/chr_synonyms.txt"
export CLINVAR="$CACHE/vep_custom/clinvar_20240611_PLPC_new_CPLP.vcf.gz,ClinVar,vcf,exact,0,ID,CLNSIG,CLNDN,CLNHGVS,CLNSIGINCL,CLNVC,GENEINFO,CLNDISDB,CLNSIGCONF,CLNREVSTAT,CLNDNINCL"
export SZA="$CACHE/vep_custom/Whole_out_sorted_2024.vcf.gz,Database,vcf,exact,0,Type,SZAID"
# plugin data (staged from Vigil)
export PDIR="$CACHE/VEP_plugins"
export DBNSFP="$CACHE/dbNSFP4.7a/dbNSFP4.7a_grch38.gz"
export DBSCSNV="$CACHE/dbscSNV1.1/dbscSNV1.1_GRCh38.txt.gz"
export SAI_SNV="$CACHE/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz"
export SAI_IND="$CACHE/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz"
export AM="$CACHE/AlphaMissense/AlphaMissense_hg38.tsv.gz"
CHUNK=${CHUNK:-1000000}
mkdir -p "$OUT"
process_chunk(){
  local cf="$1"; local base="${cf%.vcf.gz}"
  $RUN vep --cache --dir_cache "$CACHE" --offline --cache_version 113 --fork "$FORK" --format vcf --vcf \
    --compress_output bgzip --assembly GRCh38 --force_overwrite --refseq --fasta "$FASTA" --synonyms "$SYN" \
    --variant_class --symbol --canonical --pick --pick_order mane_select,rank --hgvs --check_existing \
    --af --af_gnomade --af_gnomadg --max_af --sift b --polyphen b --custom "$CLINVAR" --custom "$SZA" --no_stats \
    --dir_plugins "$PDIR" \
    --plugin dbNSFP,"$DBNSFP",CADD_phred,REVEL_score,REVEL_rankscore,BayesDel_addAF_score,BayesDel_addAF_pred,BayesDel_addAF_rankscore,BayesDel_noAF_score,BayesDel_noAF_pred,BayesDel_noAF_rankscore \
    --plugin dbscSNV,"$DBSCSNV" \
    --plugin SpliceAI,snv="$SAI_SNV",indel="$SAI_IND" \
    --plugin AlphaMissense,file="$AM" \
    -i "$cf" -o "${base}.vep.vcf.gz" 2>"${base}.err"
  $RUN python3 "$SRCDIR/gvep_build_cache.py" --vep-vcf "${base}.vep.vcf.gz" --input-vcf "$cf" --out "${base}.cache.tsv.gz" 2>>"${base}.err"
  rm -f "$cf" "${base}.vep.vcf.gz"
}
export -f process_chunk
for chr in $CHROMS; do
  fa_chr="${chr#chr}"; done_f="$OUT/cache_${chr}.tsv.gz"
  [ -s "$done_f" ] && { echo "[$chr] done, skip"; continue; }
  AV="$CACHE/homo_sapiens_refseq/113_GRCh38/${fa_chr}/all_vars.gz"
  [ -f "$AV" ] || { echo "[$chr] no all_vars, skip"; continue; }
  cd="$OUT/tmp_$chr"; rm -rf "$cd"; mkdir -p "$cd"
  echo "[$chr] $(date) extract"
  $RUN python3 "$SRCDIR/extract_known.py" --all-vars "$AV" --fasta "$FASTA" --fasta-chrom "$fa_chr" --out-chrom "$chr" --out "$cd/body" 2>"$OUT/${chr}.log"
  echo "[$chr] $(date) norm (left-align+trim)"
  { printf '##fileformat=VCFv4.2\n##contig=<ID=%s>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' "$fa_chr"; awk -F'\t' -v OFS='\t' '{sub(/^chr/,"",$1); print}' "$cd/body"; } \
    | $RUN bcftools norm -f "$FASTA" -c w - 2>"$OUT/${chr}.norm.err" \
    | awk -F'\t' -v OFS='\t' '$0!~/^#/{print "chr"$1,$2,$3,$4,$5,$6,$7,$8}' > "$cd/body_norm"
  mv "$cd/body_norm" "$cd/body"
  echo "[$chr] $(date) split+VEP"
  sort -t$'\t' -k2,2n "$cd/body" | split -l "$CHUNK" - "$cd/c"; rm -f "$cd/body"
  for c in "$cd"/c*; do { echo '##fileformat=VCFv4.2'; printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'; cat "$c"; } | $RUN bgzip > "$c.vcf.gz"; rm -f "$c"; done
  nchunks=$(ls "$cd"/c*.vcf.gz 2>/dev/null | wc -l)
  echo "[$chr] $(date) VEP+plugins: $nchunks chunks, $NPAR parallel x fork $FORK"
  ls "$cd"/c*.vcf.gz | xargs -P "$NPAR" -I CF bash -c 'process_chunk "$@"' _ CF
  cat "$cd"/c*.cache.tsv.gz > "$done_f"
  rm -rf "$cd"
  echo "[$chr] $(date) DONE $done_f ($(du -h "$done_f" 2>/dev/null|cut -f1))"
done
echo "[worker4] DONE $(date)"
