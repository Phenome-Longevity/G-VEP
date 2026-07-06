import numpy as np, sys
raw = np.fromfile(sys.argv[1], dtype=np.uint64).reshape(-1,2)
h = raw[:,0].copy(); v = raw[:,1].copy()
del raw
idx = np.argsort(h, kind='stable')
h[idx].tofile(sys.argv[2]); v[idx].tofile(sys.argv[3])
print(f"sorted {len(h):,} entries -> {sys.argv[2]}, {sys.argv[3]}")
