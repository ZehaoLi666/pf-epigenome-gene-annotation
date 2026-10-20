import sys, gzip
A = 'CTGTAGGCACCATCAAT'; n = k = 0; out = sys.stdout
with sys.stdin as f:
    while True:
        h = f.readline()
        if not h: break
        s = f.readline().rstrip(); f.readline(); q = f.readline().rstrip(); n += 1
        i = s.find(A[:10])
        if i < 0:
            for L in range(9, 4, -1):
                if s.endswith(A[:L]): i = len(s) - L; break
        if i < 25 or i > 35: continue
        out.write(f'@r{n}\n{s[:i]}\n+\n{q[:i]}\n'); k += 1
print('reads', n, 'kept', k, file=sys.stderr)
