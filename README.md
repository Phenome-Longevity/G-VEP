# G-VEP

Byte-identical VEP annotation **with all plugins** in minutes, via memoization.

G-VEP produces output **byte-for-byte identical to Ensembl VEP run with every plugin**
(dbNSFP → CADD/REVEL/BayesDel, dbscSNV, SpliceAI, AlphaMissense) by precomputing the
annotation for every known variant into a cache + index, and falling back to live VEP
only for genuinely novel variants.

## Performance

| Input | G-VEP | Real VEP + all plugins |
|-------|-------|------------------------|
| Whole genome (~5M variants) | **~4 min** | ~100 min |
| Small VCF (gene panel / exome slice) | **~seconds** | minutes |

- **Byte-identical** to VEP-with-all-plugins run on the normalized VCF (verified: 0 diffs).
- **0 variants skipped** — output matches VEP's exact variant set.
- ~97% of a typical genome is served from the cache; the ~3% novel/alt-contig variants
  are run through real live VEP-with-all-plugins, so every variant is exact.

## How it works

1. **Normalize** the input (`bcftools norm`; SNVs unchanged, indels left-aligned/trimmed).
2. **Index lookup** — each variant's annotation is fetched from the precomputed cache via a
   sorted `FNV-hash → bgzf-offset` index (binary search + `bgzf` seek/read, parallel across cores).
3. **Fallback** — cache misses (novel + alt-contig variants) run through live VEP-with-all-plugins.
4. **Assemble** in input order; emit.

The cache stores *real* VEP-with-all-plugins output, so hits are exact; misses are exact by
construction. No GPU step, no plugin approximation.

## Run

```bash
./run_gvep.sh <input.vcf.gz> [threads]      # threads defaults to all cores
# output: outputs/out/<sample>_pluginB.vcf.gz
```

Requires the cache + index in `$GVEP_DATA` (default `$HOME`):

```
$GVEP_DATA/gvep_index/{cache.bgz, cache.bgz.gzi, h_sorted.bin, v_sorted.bin}
```

## Setup

```bash
./setup.sh        # micromamba env (VEP 113) + compile tools (+ fetch cache/index)
```

## Layout

```
run_gvep.sh        main pipeline (normalize → lookup → fallback → assemble)
src/gvep_memo.py   memo engine (cache serve + live-VEP fallback, parallel bgzip I/O)
src/lookup.c       index lookup (FNV hash → binary search → bgzf seek+read)
src/build_index.c  index builder (cache → sorted hash→offset index)
src/sort_index.py  sorts the raw index by hash
build/             cache-build tooling (worker4.sh, extract_known.py, gvep_build_cache.py)
run_pipeline.sh    (legacy) prior GPU-plugin pipeline — superseded by run_gvep.sh
```

## Building the cache (one-time)

The cache is built by running real VEP-with-all-plugins on every known variant (normalized),
sharded by chromosome — see `build/worker4.sh`. ~1.3B variants → `gvep_cache_plugins/`
(per-chromosome `key\tsuffix`). The index is then built with `bin/build_index` + `src/sort_index.py`.
