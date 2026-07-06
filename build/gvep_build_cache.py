#!/usr/bin/env python3
"""Build a GVEP cache: variant key -> exact VEP INFO-suffix (the fields VEP appends).
The suffix = INFO fields whose ID exists in the VEP-annotated header but NOT in the original
input header (i.e., CSQ + any --custom fields), preserved in VEP's output order.
By the proven determinism property this suffix is variant-intrinsic and reusable on any sample.
"""
import gzip, argparse, zlib
def op(p, m='rt'): return gzip.open(p, m) if p.endswith('.gz') else open(p, m)
def info_ids(lines):
    s = set()
    for l in lines:
        if l.startswith('##INFO=<ID='):
            s.add(l.split('ID=', 1)[1].split(',', 1)[0])
    return s
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vep-vcf', required=True)    # VEP-annotated VCF = source of truth
    ap.add_argument('--input-vcf', required=True)  # original input (to know pre-existing INFO IDs)
    ap.add_argument('--out', required=True)         # cache: key\tsuffix (gz)
    ap.add_argument('--split-mod', type=int, default=0)   # if >0, KEEP only keys where crc32%mod != split-rem (the "known" set)
    ap.add_argument('--split-rem', type=int, default=0)
    a = ap.parse_args()
    in_hdr = []
    with op(a.input_vcf) as f:
        for l in f:
            if l.startswith('##'): in_hdr.append(l)
            else: break
    in_ids = info_ids(in_hdr)
    vep_hdr = []; added = None; n = 0; kept = 0
    with op(a.vep_vcf) as f, op(a.out, 'wt') as o:
        for l in f:
            if l.startswith('##'): vep_hdr.append(l); continue
            if l.startswith('#'):
                added = info_ids(vep_hdr) - in_ids
                o.write('#vep_added_keys\t' + ','.join(sorted(added)) + '\n')
                continue
            cols = l.rstrip('\n').split('\t')
            chrom, pos, ref, alt = cols[0], cols[1], cols[3], cols[4]
            n += 1
            if a.split_mod and (zlib.crc32(f"{chrom}:{pos}:{ref}:{alt}".encode()) % a.split_mod) == a.split_rem:
                continue  # this key is held out as "novel" -> not in cache
            info = cols[7]
            suffix = ';'.join(fld for fld in info.split(';') if fld.split('=', 1)[0] in added)
            o.write(f"{chrom}:{pos}:{ref}:{alt}\t{suffix}\n")
            kept += 1
    print(f"build_cache: {n:,} variants seen, {kept:,} written to cache, vep_added_keys={sorted(added)}")
if __name__ == '__main__':
    main()
