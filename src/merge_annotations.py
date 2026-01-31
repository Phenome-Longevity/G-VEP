#!/usr/bin/env python3

import os
import sys
import gzip
import argparse
from typing import Dict, List, Tuple, Optional
import time


# ============================================================
# CSQ FIELD DEFINITIONS
# ============================================================

GPU_FIELDS_ORDER = [
    # From dbscSNV plugin
    'ada_score',
    'rf_score',
    # From REVEL plugin  
    'REVEL',
    # From dbNSFP plugin
    'BayesDel_addAF_pred',
    'BayesDel_addAF_rankscore',
    'BayesDel_addAF_score',
    'BayesDel_noAF_pred',
    'BayesDel_noAF_rankscore',
    'BayesDel_noAF_score',
    'REVEL_rankscore',
    'REVEL_score',
    # From SpliceAI plugin
    'SpliceAI_pred_DP_AG',
    'SpliceAI_pred_DP_AL',
    'SpliceAI_pred_DP_DG',
    'SpliceAI_pred_DP_DL',
    'SpliceAI_pred_DS_AG',
    'SpliceAI_pred_DS_AL',
    'SpliceAI_pred_DS_DG',
    'SpliceAI_pred_DS_DL',
    'SpliceAI_pred_SYMBOL',
    # From AlphaMissense plugin
    'am_class',
    'am_genome',
    'am_pathogenicity',
    'am_protein_variant',
    'am_transcript_id',
    'am_uniprot_id',
]

# Map from GPU TSV column names to CSQ field names
GPU_TO_CSQ_MAP = {
    'dbscsnv_ada_score': 'ada_score',
    'dbscsnv_rf_score': 'rf_score',
    'revel_REVEL': 'REVEL',
    'dbnsfp_REVEL_score': 'REVEL_score',
    'dbnsfp_BayesDel_addAF_score': 'BayesDel_addAF_score',
    'dbnsfp_BayesDel_noAF_score': 'BayesDel_noAF_score',
    'dbnsfp_CADD_phred': 'CADD_phred',
    'spliceai_DS_AG': 'SpliceAI_pred_DS_AG',
    'spliceai_DS_AL': 'SpliceAI_pred_DS_AL',
    'spliceai_DS_DG': 'SpliceAI_pred_DS_DG',
    'spliceai_DS_DL': 'SpliceAI_pred_DS_DL',
    'alphamissense_am_pathogenicity': 'am_pathogenicity',
    'alphamissense_am_class': 'am_class',
    'clinvar_CLNSIG': 'ClinVar_CLNSIG_gpu',
    'clinvar_CLNREVSTAT': 'ClinVar_CLNREVSTAT_gpu',
}

# AlphaMissense class mapping (int to string)
AM_CLASS_MAP = {
    '0': 'likely_benign',
    '0.0': 'likely_benign',
    '1': 'ambiguous',
    '1.0': 'ambiguous', 
    '2': 'likely_pathogenic',
    '2.0': 'likely_pathogenic',
}

# BayesDel prediction thresholds
BAYESDEL_ADDAF_THRESHOLD = 0.0692655
BAYESDEL_NOAF_THRESHOLD = -0.0570105


def load_gpu_annotations(tsv_path: str) -> Tuple[Dict, List[str]]:
    """Load GPU annotations keyed by chr:pos:ref:alt."""
    annotations = {}
    
    opener = gzip.open if tsv_path.endswith('.gz') else open
    with opener(tsv_path, 'rt') as f:
        header = f.readline().strip().split('\t')
        # Skip CHROM, POS, REF, ALT
        gpu_fields = header[4:]
        
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 4:
                continue
            
            chrom, pos, ref, alt = cols[:4]
            # Normalize chromosome
            if chrom.startswith('chr'):
                chrom_norm = chrom[3:]
            else:
                chrom_norm = chrom
            
            # Store with both formats for lookup
            key1 = f"{chrom}:{pos}:{ref}:{alt}"
            key2 = f"chr{chrom_norm}:{pos}:{ref}:{alt}"
            key3 = f"{chrom_norm}:{pos}:{ref}:{alt}"
            
            values = dict(zip(gpu_fields, cols[4:]))
            annotations[key1] = values
            annotations[key2] = values
            annotations[key3] = values
    
    return annotations, gpu_fields


def get_csq_value(gpu_ann: Optional[Dict], csq_field: str, gpu_fields: List[str]) -> str:
    """Get value for a CSQ field from GPU annotations."""
    if gpu_ann is None:
        return ''
    
    # Find the GPU field that maps to this CSQ field
    gpu_field = None
    for gf, cf in GPU_TO_CSQ_MAP.items():
        if cf == csq_field:
            gpu_field = gf
            break
    
    if gpu_field is None or gpu_field not in gpu_ann:
        # Handle derived fields
        if csq_field == 'BayesDel_addAF_pred':
            score = gpu_ann.get('dbnsfp_BayesDel_addAF_score', '.')
            if score and score != '.':
                try:
                    return 'D' if float(score) > BAYESDEL_ADDAF_THRESHOLD else 'T'
                except:
                    pass
            return ''
        elif csq_field == 'BayesDel_noAF_pred':
            score = gpu_ann.get('dbnsfp_BayesDel_noAF_score', '.')
            if score and score != '.':
                try:
                    return 'D' if float(score) > BAYESDEL_NOAF_THRESHOLD else 'T'
                except:
                    pass
            return ''
        elif csq_field == 'BayesDel_addAF_rankscore':
            # We don't have rankscore from GPU, leave empty
            return ''
        elif csq_field == 'BayesDel_noAF_rankscore':
            return ''
        elif csq_field == 'REVEL_rankscore':
            return ''
        elif csq_field in ['SpliceAI_pred_DP_AG', 'SpliceAI_pred_DP_AL', 
                           'SpliceAI_pred_DP_DG', 'SpliceAI_pred_DP_DL',
                           'SpliceAI_pred_SYMBOL']:
            # We don't have DP scores or symbol from GPU indices
            return ''
        elif csq_field in ['am_genome', 'am_protein_variant', 'am_transcript_id', 'am_uniprot_id']:
            # We don't have these from GPU indices
            return ''
        return ''
    
    value = gpu_ann[gpu_field]
    
    # Handle special conversions
    if csq_field == 'am_class' and value and value != '.':
        return AM_CLASS_MAP.get(value, '')
    
    # Convert '.' to empty string for CSQ format
    if value == '.' or value == 'nan':
        return ''
    
    return value


def parse_csq_header(header_line: str) -> List[str]:
    """Parse CSQ field names from VEP header line."""
    if 'Format:' not in header_line:
        return []
    
    format_part = header_line.split('Format:')[1]
    format_part = format_part.rstrip('">')
    return format_part.split('|')


def find_insert_position(csq_fields: List[str]) -> int:
    """Find where to insert GPU fields (after MAX_AF_POPS, before CLIN_SIG)."""
    for i, field in enumerate(csq_fields):
        if field == 'CLIN_SIG':
            return i
    # Fallback: look for LoF
    for i, field in enumerate(csq_fields):
        if field == 'LoF':
            return i
    # Otherwise append at end
    return len(csq_fields)


def merge_vcf(vep_vcf: str, gpu_annotations: Dict, gpu_fields: List[str],
              output_vcf: str, verbose: bool = True):
    """Merge GPU annotations into VEP CSQ field."""
    
    in_opener = gzip.open if vep_vcf.endswith('.gz') else open
    out_opener = gzip.open if output_vcf.endswith('.gz') else open
    out_mode = 'wt' if output_vcf.endswith('.gz') else 'w'
    
    csq_fields = []
    insert_pos = 0
    new_csq_header = None
    
    count = 0
    annotated = 0
    
    with in_opener(vep_vcf, 'rt') as fin, out_opener(output_vcf, out_mode) as fout:
        for line in fin:
            # Handle header lines
            if line.startswith('##INFO=<ID=CSQ'):
                # Parse existing CSQ fields
                csq_fields = parse_csq_header(line)
                insert_pos = find_insert_position(csq_fields)
                
                # Build new CSQ header with GPU fields inserted
                new_fields = csq_fields[:insert_pos] + GPU_FIELDS_ORDER + csq_fields[insert_pos:]
                
                # Reconstruct header line
                new_format = '|'.join(new_fields)
                new_csq_header = f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: {new_format}">\n'
                fout.write(new_csq_header)
                
                if verbose:
                    print(f"  Original CSQ fields: {len(csq_fields)}")
                    print(f"  New CSQ fields: {len(new_fields)}")
                    print(f"  Inserted {len(GPU_FIELDS_ORDER)} GPU fields at position {insert_pos}")
                continue
            
            if line.startswith('#'):
                fout.write(line)
                continue
            
            # Data line
            cols = line.strip().split('\t')
            if len(cols) < 8:
                fout.write(line)
                continue
            
            chrom, pos, id_, ref, alt, qual, filt, info = cols[:8]
            rest = cols[8:] if len(cols) > 8 else []
            
            # Parse INFO field to find CSQ
            info_parts = info.split(';')
            new_info_parts = []
            
            for part in info_parts:
                if part.startswith('CSQ='):
                    # Process CSQ field
                    csq_value = part[4:]  # Remove 'CSQ='
                    transcripts = csq_value.split(',')
                    new_transcripts = []
                    
                    # Handle multi-allelic
                    alts = alt.split(',')
                    
                    for transcript in transcripts:
                        fields = transcript.split('|')
                        
                        # Determine which alt allele this transcript is for
                        # (first field is the allele)
                        allele = fields[0] if fields else ''
                        
                        # Try to match to an alt allele
                        matched_alt = None
                        for a in alts:
                            if allele == a or allele == '-' or a == allele:
                                matched_alt = a
                                break
                        if matched_alt is None and len(alts) == 1:
                            matched_alt = alts[0]
                        
                        # Look up GPU annotations
                        gpu_ann = None
                        if matched_alt:
                            key = f"{chrom}:{pos}:{ref}:{matched_alt}"
                            gpu_ann = gpu_annotations.get(key)
                            if gpu_ann:
                                annotated += 1
                        
                        # Build new field values for GPU columns
                        gpu_values = []
                        for csq_field in GPU_FIELDS_ORDER:
                            val = get_csq_value(gpu_ann, csq_field, gpu_fields)
                            gpu_values.append(val)
                        
                        # Insert GPU values at correct position
                        new_fields = fields[:insert_pos] + gpu_values + fields[insert_pos:]
                        new_transcripts.append('|'.join(new_fields))
                    
                    new_info_parts.append('CSQ=' + ','.join(new_transcripts))
                else:
                    new_info_parts.append(part)
            
            new_info = ';'.join(new_info_parts)
            
            # Reconstruct line
            new_cols = [chrom, pos, id_, ref, alt, qual, filt, new_info] + rest
            fout.write('\t'.join(new_cols) + '\n')
            
            count += 1
            if verbose and count % 500000 == 0:
                print(f"  Processed {count:,} variants...")
    
    return count, annotated


def main():
    parser = argparse.ArgumentParser(description='Merge GPU annotations into VEP CSQ field')
    parser.add_argument('--vep-vcf', required=True, help='VEP annotated VCF (lightweight)')
    parser.add_argument('--gpu-tsv', required=True, help='GPU annotations TSV')
    parser.add_argument('--output', '-o', required=True, help='Output VCF')
    parser.add_argument('--quiet', action='store_true')
    args = parser.parse_args()
    
    verbose = not args.quiet
    
    if verbose:
        print(f"\n{'='*60}")
        print(f"Merging GPU Annotations into CSQ")
        print(f"{'='*60}")
        print(f"VEP VCF:  {args.vep_vcf}")
        print(f"GPU TSV:  {args.gpu_tsv}")
        print(f"Output:   {args.output}")
    
    start = time.time()
    
    # Load GPU annotations
    if verbose:
        print(f"\nLoading GPU annotations...")
    load_start = time.time()
    gpu_annotations, gpu_fields = load_gpu_annotations(args.gpu_tsv)
    if verbose:
        print(f"  Loaded {len(gpu_annotations)//3:,} variants")
        print(f"  GPU fields: {gpu_fields}")
        print(f"  Done in {time.time() - load_start:.1f}s")
    
    # Merge
    if verbose:
        print(f"\nMerging into CSQ...")
    merge_start = time.time()
    count, annotated = merge_vcf(args.vep_vcf, gpu_annotations, gpu_fields,
                                  args.output, verbose)
    
    elapsed = time.time() - start
    
    if verbose:
        print(f"\n{'='*60}")
        print(f"Completed in {elapsed:.1f}s")
        print(f"Variants processed: {count:,}")
        print(f"Transcripts with GPU annotations: {annotated:,}")
        print(f"Output: {args.output}")
        print(f"{'='*60}")


if __name__ == '__main__':
    main()