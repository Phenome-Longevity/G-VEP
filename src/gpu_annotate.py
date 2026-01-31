#!/usr/bin/env python3
"""
GPU Variant Annotator v4
- Matches build_indices v4 (no strings)
- Fast GPU binary search
"""

import sys
import argparse
import numpy as np
from pathlib import Path
import gzip
import json
import time
from typing import List, Tuple

# ============================================================
# GPU DETECTION
# ============================================================

try:
    import cupy as xp
    GPU_AVAILABLE = True
    print("✓ CuPy - GPU enabled")
except ImportError:
    import numpy as xp
    GPU_AVAILABLE = False
    print("✗ CuPy not found - CPU mode")

# ============================================================
# VARIANT KEY ENCODING
# ============================================================

CHROM_MAP = {
    '1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8,
    '9': 9, '10': 10, '11': 11, '12': 12, '13': 13, '14': 14, '15': 15,
    '16': 16, '17': 17, '18': 18, '19': 19, '20': 20, '21': 21, '22': 22,
    'X': 23, 'Y': 24, 'MT': 25, 'M': 25,
    'chr1': 1, 'chr2': 2, 'chr3': 3, 'chr4': 4, 'chr5': 5, 'chr6': 6,
    'chr7': 7, 'chr8': 8, 'chr9': 9, 'chr10': 10, 'chr11': 11, 'chr12': 12,
    'chr13': 13, 'chr14': 14, 'chr15': 15, 'chr16': 16, 'chr17': 17,
    'chr18': 18, 'chr19': 19, 'chr20': 20, 'chr21': 21, 'chr22': 22,
    'chrX': 23, 'chrY': 24, 'chrMT': 25, 'chrM': 25
}

def fnv1a_hash(s: str) -> int:
    h = 0x811c9dc5
    for c in s.encode():
        h ^= c
        h = (h * 0x01000193) & 0xFFFFFFFF
    return h

def normalize_chrom(chrom: str) -> str:
    if chrom.startswith('chr'):
        return chrom[3:]
    return chrom

def encode_variant_key(chrom: str, pos: int, ref: str, alt: str) -> np.uint64:
    chrom_norm = normalize_chrom(chrom)
    chrom_int = CHROM_MAP.get(chrom_norm, 0)
    if chrom_int == 0:
        return np.uint64(0)
    
    ref_hash = fnv1a_hash(ref) & 0x3FFF
    alt_hash = fnv1a_hash(alt) & 0x3FFF
    pos = min(pos, 0x0FFFFFFF)
    
    key = (int(chrom_int) << 56) | (int(pos) << 28) | (int(ref_hash) << 14) | int(alt_hash)
    return np.uint64(key)


# ============================================================
# GPU BINARY SEARCH
# ============================================================

if GPU_AVAILABLE:
    binary_search_kernel = xp.RawKernel(r'''
    extern "C" __global__
    void binary_search(
        const unsigned long long* query_keys,
        const int num_queries,
        const unsigned long long* index_keys,
        const float* index_values,
        const int index_size,
        const int num_fields,
        float* results
    ) {
        int tid = blockDim.x * blockIdx.x + threadIdx.x;
        if (tid >= num_queries) return;
        
        unsigned long long query = query_keys[tid];
        
        int lo = 0;
        int hi = index_size;
        
        while (lo < hi) {
            int mid = lo + (hi - lo) / 2;
            if (index_keys[mid] < query) {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        
        if (lo < index_size && index_keys[lo] == query) {
            for (int f = 0; f < num_fields; f++) {
                results[tid * num_fields + f] = index_values[lo * num_fields + f];
            }
        } else {
            for (int f = 0; f < num_fields; f++) {
                results[tid * num_fields + f] = nanf("");
            }
        }
    }
    ''', 'binary_search')


def gpu_annotate(query_keys, index_keys, index_values, num_fields):
    num_queries = len(query_keys)
    
    query_gpu = xp.asarray(query_keys)
    keys_gpu = xp.asarray(index_keys)
    vals_gpu = xp.asarray(index_values)
    results = xp.full((num_queries, num_fields), np.nan, dtype=xp.float32)
    
    block_size = 256
    grid_size = (num_queries + block_size - 1) // block_size
    
    binary_search_kernel(
        (grid_size,), (block_size,),
        (query_gpu, num_queries, keys_gpu, vals_gpu, len(index_keys), num_fields, results)
    )
    
    return xp.asnumpy(results)


def cpu_annotate(query_keys, index_keys, index_values, num_fields):
    results = np.full((len(query_keys), num_fields), np.nan, dtype=np.float32)
    
    for i, qk in enumerate(query_keys):
        idx = np.searchsorted(index_keys, qk)
        if idx < len(index_keys) and index_keys[idx] == qk:
            results[i] = index_values[idx]
    
    return results


# ============================================================
# INDEX LOADING
# ============================================================

class AnnotationIndex:
    def __init__(self, name: str, index_dir: Path):
        self.name = name
        
        meta_path = index_dir / f"{name}.meta.json"
        with open(meta_path) as f:
            self.metadata = json.load(f)
        
        self.num_variants = self.metadata['num_variants']
        self.num_fields = self.metadata['num_fields']
        self.field_names = self.metadata['field_names']
        
        keys_path = index_dir / f"{name}.keys.bin"
        vals_path = index_dir / f"{name}.vals.bin"
        
        self.keys = np.fromfile(keys_path, dtype=np.uint64)
        self.values = np.fromfile(vals_path, dtype=np.float32).reshape(-1, self.num_fields)
        
        mem_mb = (self.keys.nbytes + self.values.nbytes) / 1e6
        print(f"    Loaded: {mem_mb:.1f} MB")
    
    def annotate(self, query_keys: np.ndarray) -> np.ndarray:
        if GPU_AVAILABLE:
            return gpu_annotate(query_keys, self.keys, self.values, self.num_fields)
        else:
            return cpu_annotate(query_keys, self.keys, self.values, self.num_fields)


# ============================================================
# VCF PARSING
# ============================================================

def parse_vcf(vcf_path: str) -> Tuple[List[Tuple], np.ndarray]:
    variants = []
    keys = []
    
    opener = gzip.open if vcf_path.endswith('.gz') else open
    
    with opener(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            cols = line.strip().split('\t')
            if len(cols) < 5:
                continue
            
            chrom, pos, _, ref, alt = cols[:5]
            
            try:
                pos = int(pos)
            except ValueError:
                continue
            
            for a in alt.split(','):
                key = encode_variant_key(chrom, pos, ref, a)
                if key == 0:
                    continue
                variants.append((chrom, pos, ref, a))
                keys.append(key)
    
    return variants, np.array(keys, dtype=np.uint64)


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(description='GPU Variant Annotator v4')
    parser.add_argument('-i', '--input', required=True, help='Input VCF')
    parser.add_argument('-o', '--output', required=True, help='Output TSV')
    parser.add_argument('--index-dir', required=True, help='Index directory')
    parser.add_argument('--databases', nargs='+', help='Specific databases')
    args = parser.parse_args()
    
    index_dir = Path(args.index_dir)
    
    available_dbs = []
    for meta_file in index_dir.glob('*.meta.json'):
        db_name = meta_file.stem.replace('.meta', '')
        if args.databases is None or db_name in args.databases:
            available_dbs.append(db_name)
    
    if not available_dbs:
        print("ERROR: No indices found")
        sys.exit(1)
    
    print("=" * 60)
    print("GPU Variant Annotator v4")
    print("=" * 60)
    print(f"Input:  {args.input}")
    print(f"Output: {args.output}")
    print(f"DBs:    {available_dbs}")
    print(f"GPU:    {'Yes' if GPU_AVAILABLE else 'No'}")
    
    # Load indices
    print("\nLoading indices...")
    indices = {}
    for db_name in available_dbs:
        print(f"  {db_name}:", end=" ")
        try:
            idx = AnnotationIndex(db_name, index_dir)
            indices[db_name] = idx
            print(f"{idx.num_variants:,} variants")
        except Exception as e:
            print(f"ERROR: {e}")
    
    # Parse VCF
    print("\nParsing VCF...")
    start = time.time()
    variants, query_keys = parse_vcf(args.input)
    print(f"  {len(variants):,} variants in {time.time() - start:.1f}s")
    
    # Annotate
    print("\nAnnotating...")
    all_results = {}
    
    for db_name, index in indices.items():
        start = time.time()
        results = index.annotate(query_keys)
        hits = np.sum(~np.isnan(results[:, 0]))
        elapsed = time.time() - start
        rate = len(query_keys) / elapsed if elapsed > 0 else 0
        pct = 100 * hits / len(query_keys) if len(query_keys) > 0 else 0
        print(f"  {db_name}: {pct:.1f}% hits, {rate:,.0f} var/sec")
        all_results[db_name] = (results, index.field_names)
    
    # Write output
    print("\nWriting output...")
    start = time.time()
    
    with open(args.output, 'w') as f:
        header = ['#CHROM', 'POS', 'REF', 'ALT']
        for db_name, (_, field_names) in all_results.items():
            for fn in field_names:
                header.append(f"{db_name}_{fn}")
        f.write('\t'.join(header) + '\n')
        
        for i, (chrom, pos, ref, alt) in enumerate(variants):
            row = [chrom, str(pos), ref, alt]
            for db_name, (results, _) in all_results.items():
                for val in results[i]:
                    if np.isnan(val):
                        row.append('.')
                    elif val == int(val):
                        row.append(str(int(val)))
                    else:
                        row.append(f"{val:.6g}")
            f.write('\t'.join(row) + '\n')
    
    print(f"  Done in {time.time() - start:.1f}s")
    print("=" * 60)
    print(f"Output: {args.output}")


if __name__ == '__main__':
    main()