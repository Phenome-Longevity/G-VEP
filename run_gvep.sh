#!/bin/bash
# G-VEP Path B: normalize -> plugin-cache lookup (exact base+indels+plugins) -> live-VEP-with-plugins fallback.
# No GPU step, no merge: plugins are baked into the cache. Output is normalized coords = byte-identical to VEP-with-plugins-on-normalized.
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# GVEP_DATA holds the plugin-cache + index (downloaded separately; see setup.sh). Defaults to $HOME.
GVEP_DATA="${GVEP_DATA:-$HOME}"
IN="${1:?usage: run_gvep.sh <input.vcf.gz> [threads]}"
THREADS=${2:-$(nproc)}
SAMPLE=$(basename "$IN" | sed 's/\.vcf\.gz$//;s/\.vcf$//')
OUT="${GVEP_OUT:-$SCRIPT_DIR/outputs}"; mkdir -p "$OUT"/{vcf_pass,norm,out,work}
RUN="$HOME/bin/micromamba run -p $HOME/mamba/gvep"
IDX=$GVEP_DATA/gvep_index
BIN=$SCRIPT_DIR/bin
VC=$HOME/gvep/data/vep_cache
DC=$SCRIPT_DIR/data/vep_cache
FASTA=$DC/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
SYN=$VC/homo_sapiens_refseq/113_GRCh38/chr_synonyms.txt
MEMO=$SCRIPT_DIR/src/gvep_memo.py
CLINVAR="$VC/vep_custom/clinvar_20240611_PLPC_new_CPLP.vcf.gz,ClinVar,vcf,exact,0,ID,CLNSIG,CLNDN,CLNHGVS,CLNSIGINCL,CLNVC,GENEINFO,CLNDISDB,CLNSIGCONF,CLNREVSTAT,CLNDNINCL"
SZA="$VC/vep_custom/Whole_out_sorted_2024.vcf.gz,Database,vcf,exact,0,Type,SZAID"
W=$OUT/work
PASS=$OUT/vcf_pass/${SAMPLE}_PASS.vcf.gz
NORM=$OUT/norm/${SAMPLE}_norm.vcf.gz
FINAL=$OUT/out/${SAMPLE}_pluginB.vcf.gz
# VEP-with-all-plugins flags (for fallback) = exactly the cache-build config
vepf='--cache --dir_cache '"$VC"' --offline --cache_version 113 --format vcf --vcf --compress_output bgzip --assembly GRCh38 --force_overwrite --refseq --fasta '"$FASTA"' --synonyms '"$SYN"' --variant_class --symbol --canonical --pick --pick_order mane_select,rank --hgvs --check_existing --af --af_gnomade --af_gnomadg --max_af --sift b --polyphen b --custom '"$CLINVAR"' --custom '"$SZA"' --no_stats --dir_plugins '"$DC"'/VEP_plugins --plugin dbNSFP,'"$DC"'/dbNSFP4.7a/dbNSFP4.7a_grch38.gz,CADD_phred,REVEL_score,REVEL_rankscore,BayesDel_addAF_score,BayesDel_addAF_pred,BayesDel_addAF_rankscore,BayesDel_noAF_score,BayesDel_noAF_pred,BayesDel_noAF_rankscore --plugin dbscSNV,'"$DC"'/dbscSNV1.1/dbscSNV1.1_GRCh38.txt.gz --plugin SpliceAI,snv='"$DC"'/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel='"$DC"'/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz --plugin AlphaMissense,file='"$DC"'/AlphaMissense/AlphaMissense_hg38.tsv.gz'

echo "===== G-VEP pluginB: $SAMPLE ====="
$RUN bcftools view -f PASS,. "$IN" -Oz -o "$PASS" 2>>$W/err || zcat "$IN" | awk '/^#/||$7=="PASS"||$7=="."' | $RUN bgzip > "$PASS"
NV=$(zcat "$PASS"|grep -vc '^#'); echo "[1/4] PASS: $NV (${SECONDS}s)"; S1=$SECONDS
# 2. normalize MAIN chroms (bcftools norm, parallel (de)compress I/O); non-main pass through -> live-VEP fallback (0 skipped)
BCF=$HOME/mamba/gvep/bin/bcftools
$RUN bgzip -dc -@ $THREADS "$PASS" | awk 'BEGIN{OFS="\t"} /^##contig=<ID=chr/{sub(/ID=chr/,"ID=");print;next} /^#/{print;next} {c=$1;sub(/^chr/,"",c);if(c~/^([1-9]|1[0-9]|2[0-2]|X|Y|MT)$/){$1=c;print}}' \
  | $BCF norm -f "$FASTA" -c w - 2>>$W/err \
  | awk 'BEGIN{OFS="\t"} /^##contig=<ID=/{sub(/ID=/,"ID=chr");print;next} /^#/{print;next} {$1="chr"$1;print}' > $W/norm.vcf
$RUN bgzip -dc -@ $THREADS "$PASS" | awk -F'\t' '$0!~/^#/{c=$1;sub(/^chr/,"",c);if(!(c~/^([1-9]|1[0-9]|2[0-2]|X|Y|MT)$/))print}' >> $W/norm.vcf
$RUN bgzip -@ $THREADS -c $W/norm.vcf > "$NORM"; rm -f $W/norm.vcf
echo "[2/4] norm: $(zcat "$NORM"|grep -vc '^#') ($((SECONDS-S1))s)"; S2=$SECONDS
# launch header-VEP in the BACKGROUND (needs only the NORM header, not the variants) -> overlaps the lookup, free
{ zcat "$NORM"|grep '^#'; zcat "$NORM"|grep -v '^#'|head -20; } | $RUN bgzip > $W/hdr_in.vcf.gz
( $RUN vep $vepf --fork 4 -i $W/hdr_in.vcf.gz -o $W/refhdr.vcf.gz 2>>$W/err ) & RHPID=$!
# 3. INDEX lookup (C: FNV hash -> binary-search sorted index -> bgzf seek+read), parallel across $THREADS cores (shared mmap'd index)
zcat "$NORM" | grep -v '^#' | awk -F'\t' '{print $1":"$2":"$4":"$5}' > $W/inkeys.txt
NK=$(wc -l < $W/inkeys.txt); NP=$(( NK < THREADS*1000 ? 1 : THREADS ))
rm -f $W/kp_* 2>/dev/null
if [ "$NP" -le 1 ]; then
  $BIN/lookup $IDX/h_sorted.bin $IDX/v_sorted.bin $IDX/cache.bgz < $W/inkeys.txt 2>$W/lookup.err | gzip > $W/sample_cache.tsv.gz
else
  split -n l/$NP -d -a 4 $W/inkeys.txt $W/kp_
  ls $W/kp_* | xargs -P $NP -I{} sh -c "$BIN/lookup $IDX/h_sorted.bin $IDX/v_sorted.bin $IDX/cache.bgz < {} > {}.out 2>/dev/null"
  cat $W/kp_*.out 2>/dev/null | gzip > $W/sample_cache.tsv.gz; rm -f $W/kp_*
fi
echo "[3/4] index lookup (${NP}-way): $(zcat $W/sample_cache.tsv.gz|wc -l)/$NV ($((SECONDS-S2))s)"; S3=$SECONDS
# 4. engine: serve NORM from cache + fallback live-VEP-with-plugins
cat > $W/fallback.sh <<FB
#!/bin/bash
$RUN vep $vepf --fork $THREADS -i "\$1/misses.vcf.gz" -o "\$1/misses_vep.vcf.gz" 2>>"\$1/fb.err"
FB
wait $RHPID 2>/dev/null   # header-VEP (was launched during the lookup) is ready by now
$RUN python3 "$MEMO" --input "$NORM" --cache $W/sample_cache.tsv.gz --output "$FINAL" --fallback-script $W/fallback.sh --work $W --ref-header $W/refhdr.vcf.gz --threads $THREADS
echo "[4/4] FINAL: $(zcat "$FINAL"|grep -vc '^#') ($((SECONDS-S3))s)"
echo "===== DONE ${SECONDS}s -> $FINAL ====="
echo "TIMING: pass=${S1}s norm=$((S2-S1))s scan=$((S3-S2))s serve=$((SECONDS-S3))s total=${SECONDS}s"
