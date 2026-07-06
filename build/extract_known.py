#!/usr/bin/env python3
"""Extract known variants from a VEP all_vars.gz -> VCF body (chrom pos . ref alt . PASS .).
SNVs are exact (the bulk). Indels are converted from VEP coords to VCF using a reference anchor
base (best-effort; any representation mismatch at runtime simply routes to the VEP backstop, which
is still byte-identical). Streams; emits sorted-within-chrom is not guaranteed (sort downstream).
"""
import sys, gzip, argparse, pysam
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--all-vars', required=True)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--fasta-chrom', required=True)  # contig name in FASTA, e.g. '22'
    ap.add_argument('--out-chrom', required=True)     # emit name, e.g. 'chr22'
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    fa = pysam.FastaFile(a.fasta)
    oc = a.out_chrom
    nsnv = nindel = nskip = 0
    opener = gzip.open if a.all_vars.endswith('.gz') else open
    BASES = ('A', 'C', 'G', 'T')
    def anchor(p):  # 1-based ref base at position p
        b = fa.fetch(a.fasta_chrom, p - 1, p).upper()
        return b if b in BASES else ''
    with opener(a.all_vars, 'rt') as f, open(a.out, 'w') as o:
        w = o.write
        for line in f:
            c = line.split('\t')
            if len(c) < 7 or not c[4].isdigit():
                continue
            start = int(c[4])
            end = int(c[5]) if c[5].isdigit() else start
            parts = c[6].split('/')
            if len(parts) < 2:
                continue
            ref = parts[0]
            for alt in parts[1:]:
                if ref in BASES and alt in BASES:
                    w(f"{oc}\t{start}\t.\t{ref}\t{alt}\t.\tPASS\t.\n"); nsnv += 1
                elif ref == '-':                       # insertion: bases inserted between end and start
                    an = anchor(end)
                    if an: w(f"{oc}\t{end}\t.\t{an}\t{an}{alt}\t.\tPASS\t.\n"); nindel += 1
                    else: nskip += 1
                elif alt == '-':                       # deletion of ref[start..end]
                    an = anchor(start - 1)
                    if an: w(f"{oc}\t{start-1}\t.\t{an}{ref}\t{an}\t.\tPASS\t.\n"); nindel += 1
                    else: nskip += 1
                elif all(b in BASES for b in ref) and all(b in BASES for b in alt):
                    w(f"{oc}\t{start}\t.\t{ref}\t{alt}\t.\tPASS\t.\n"); nindel += 1   # MNV/substitution
                else:
                    nskip += 1
    sys.stderr.write(f"{oc}: {nsnv:,} SNVs, {nindel:,} indel/MNV, {nskip:,} skipped\n")
if __name__ == '__main__':
    main()
