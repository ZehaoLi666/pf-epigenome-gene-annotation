"""Per-gene stranded DAFT-seq expression (Chappell 2020, ERP001570): mean sense-strand CPM over CDS exons (v68), per time point (0-48 h),
plus blood-stage summaries. Coverage from daft_cov10.npz (10-bp bins, RNA strand, CPM per million read pairs)."""
from pfconfig import HPCC_BIGDATA, HPCC_LAB  # paths: config/paths.yaml
import numpy as np, collections, re
z0 = np.load(f'{HPCC_BIGDATA}/flank_coverage_2026/daft_cov10.npz'); z = {k: z0[k] for k in z0.files if not k.startswith('pool')}
TPS = sorted({int(k.split('|')[0]) for k in z})
ex = collections.defaultdict(list); info = {}
for ln in open(f'{HPCC_LAB}/PlasmoDB_Genome/PlasmnDB_v68/PlasmoDB-68_Pfalciparum3D7.gff'):
    if ln.startswith('#'): continue
    f = ln.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] != 'CDS' or not re.match(r'Pf3D7_\d\d_v3', f[0]): continue
    g = re.search(r'Parent=([^;.]+)', f[8]).group(1); ex[g].append((int(f[3]) - 1, int(f[4]))); info[g] = (f[0], f[6])
out = open(f'{HPCC_BIGDATA}/flank_coverage_2026/daft_gene_expression_v68.tsv', 'w')
out.write('gene_id\tchrom\tstrand\tcds_len\t' + '\t'.join(f'h{t}' for t in TPS) + '\tantisense_max\n')
for g, iv in ex.items():
    c, s = info[g]; k = 0 if s == '+' else 1; vals = []; L = 0
    bins = np.unique(np.concatenate([np.arange(a // 10, (b + 9) // 10) for a, b in iv]))
    for t in TPS:
        a = z[f'{t}|{c}']; bb = bins[bins < a.shape[1]]; vals.append(float(a[k, bb].mean()) if len(bb) else 0.0)
    anti = max(float(z[f'{t}|{c}'][1 - k, bins[bins < z[f'{t}|{c}'].shape[1]]].mean()) for t in TPS)
    out.write(f'{g}\t{c}\t{s}\t{sum(b - a for a, b in iv)}\t' + '\t'.join(f'{v:.4f}' for v in vals) + f'\t{anti:.4f}\n')
out.close(); print('ok', len(ex), TPS)
