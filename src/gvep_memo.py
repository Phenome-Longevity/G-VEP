#!/usr/bin/env python3
"""GVEP-memo engine (fast): parallel bgzip I/O + batched assembly. Same output as gvep_memo.py.
Output INFO = original_input_INFO + ';' + cached_suffix (VEP appends at end -> byte-identical)."""
import argparse, subprocess, os, time, gzip

def reader(path, threads):
    """Stream-read a (b)gzip or plain file via parallel bgzip -dc when possible."""
    if path.endswith('.gz'):
        p = subprocess.Popen(['bgzip', '-dc', '-@', str(threads), path], stdout=subprocess.PIPE, bufsize=1<<22)
        return p, p.stdout
    return None, open(path, 'rb')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True)
    ap.add_argument('--cache', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--fallback-script', required=True)
    ap.add_argument('--work', default='outputs/gvep_work')
    ap.add_argument('--ref-header', required=True)
    ap.add_argument('--threads', type=int, default=os.cpu_count())
    a = ap.parse_args()
    T = max(1, a.threads); os.makedirs(a.work, exist_ok=True)
    t0 = time.time()
    # --- load cache (key\tsuffix) ---
    cache = {}
    cp, cf = reader(a.cache, T)
    for l in cf:
        if l[:1] == b'#': continue
        k, _, s = l.rstrip(b'\n').partition(b'\t'); cache[k] = s
    if cp: cp.wait()
    print(f"cache loaded: {len(cache):,} variants ({time.time()-t0:.1f}s)", flush=True)
    # --- read input, partition (keep raw bytes lines) ---
    hdr = []; rows = []  # rows: (rawline_bytes_no_nl, key_bytes)
    ip, inf = reader(a.input, T)
    for l in inf:
        if l[:1] == b'#': hdr.append(l); continue
        line = l.rstrip(b'\n'); c = line.split(b'\t')
        rows.append((line, c, b"%s:%s:%s:%s" % (c[0], c[1], c[3], c[4])))
    if ip: ip.wait()
    miss = [i for i, r in enumerate(rows) if r[2] not in cache]
    print(f"variants: {len(rows):,} | cache hits: {len(rows)-len(miss):,} | misses: {len(miss):,}", flush=True)
    # --- fallback VEP on misses ---
    miss_suffix = {}
    if miss:
        with gzip.open(os.path.join(a.work, 'misses.vcf.gz'), 'wb') as o:
            for l in hdr: o.write(l)
            for i in miss: o.write(rows[i][0] + b'\n')
        s = time.time(); rc = subprocess.call(['bash', a.fallback_script, a.work])
        print(f"fallback VEP on {len(miss):,} misses: rc={rc} ({time.time()-s:.1f}s)", flush=True)
        in_ids = set(l.split(b'ID=',1)[1].split(b',',1)[0] for l in hdr if l.startswith(b'##INFO=<ID='))
        vhdr = []; added = None
        mp, mf = reader(os.path.join(a.work, 'misses_vep.vcf.gz'), T)
        for l in mf:
            if l[:2] == b'##': vhdr.append(l); continue
            if l[:1] == b'#':
                added = set(x.split(b'ID=',1)[1].split(b',',1)[0] for x in vhdr if x.startswith(b'##INFO=<ID=')) - in_ids
                continue
            c = l.rstrip(b'\n').split(b'\t'); k = b"%s:%s:%s:%s" % (c[0], c[1], c[3], c[4])
            miss_suffix[k] = b';'.join(f for f in c[7].split(b';') if f.split(b'=',1)[0] in added)
        if mp: mp.wait()
    # --- assemble in input order, write via parallel bgzip, batched ---
    emitted = dropped = 0
    outf = open(a.output, 'wb')
    op = subprocess.Popen(['bgzip', '-@', str(T), '-c'], stdin=subprocess.PIPE, stdout=outf, bufsize=1<<22)
    ow = op.stdin
    with open(a.ref_header, 'rb') if not a.ref_header.endswith('.gz') else gzip.open(a.ref_header, 'rb') as f:
        for l in f:
            if l[:1] == b'#': ow.write(l)
            else: break
    buf = []; ap_ = buf.append
    for line, c, k in rows:
        suf = cache.get(k)
        if suf is None:
            suf = miss_suffix.get(k)
            if suf is None: dropped += 1; continue
        info = c[7]
        if info and info != b'.':
            c[7] = info + b';' + suf; ap_(b'\t'.join(c))
        else:
            c[7] = suf; ap_(b'\t'.join(c))
        if len(buf) >= 50000:
            ow.write(b'\n'.join(buf)); ow.write(b'\n'); buf.clear()
    if buf: ow.write(b'\n'.join(buf)); ow.write(b'\n')
    emitted = len(rows) - dropped
    ow.close(); op.wait(); outf.close()
    print(f"emitted {emitted:,} | dropped {dropped:,} | TOTAL {time.time()-t0:.1f}s", flush=True)

if __name__ == '__main__':
    main()
