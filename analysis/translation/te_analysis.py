# TE analysis: expects in kernel: cov (per-gene covariates incl. ref_/v71_ census features), pandas/numpy
import statsmodels.formula.api as smf
from scipy import stats as ss
from statsmodels.stats.multitest import multipletests

LIB = {'SRR836068': 'poly_0h', 'SRR836069': 'poly_18h', 'SRR836070': 'poly_36h_r1', 'SRR836071': 'poly_36h_r2',
       'SRR835998': 'ss_0h', 'SRR836066': 'ss_18h', 'SRR836067': 'ss_36h'}

def load_counts(path):
    fc = pd.read_csv(path, sep='\t', comment='#').set_index('Geneid')
    cols = [c for c in fc.columns if c.endswith('.bam')]
    fc = fc[cols]
    fc.columns = [LIB[c.split('/')[-1].split('.')[0]] for c in cols]
    return fc

def te_table(cnt, min_ss=20):
    cpm = cnt / cnt.sum() * 1e6
    te = pd.DataFrame(index=cnt.index)
    for st, p, s in [('0h', ['poly_0h'], 'ss_0h'), ('18h', ['poly_18h'], 'ss_18h'), ('36h', ['poly_36h_r1', 'poly_36h_r2'], 'ss_36h')]:
        pc = cpm[p].mean(axis=1)
        pr = cnt[p].sum(axis=1)
        # pseudocount 0.5 read expressed in CPM units of each library
        ps_p = np.mean([0.5 / cnt[x].sum() * 1e6 for x in p]); ps_s = 0.5 / cnt[s].sum() * 1e6
        lt = np.log2((pc + ps_p) / (cpm[s] + ps_s))
        ok = cnt[s] >= min_ss
        lt = lt.where(ok)
        te['poly_reads_' + st] = pr; te['ss_reads_' + st] = cnt[s]
        te['poly_CPM_' + st] = pc; te['ss_CPM_' + st] = cpm[s]
        te['log2TE_raw_' + st] = lt
        te['log2TE_' + st] = lt - lt.median()   # median-centred per stage
        te['TE_pass_' + st] = ok
    return te

STAGE_DAFT = {'0h': 'h0', '18h': 'h16', '36h': 'daft_36'}

def build_design(cov, te, stage):
    d = cov.copy()
    d['TE'] = te['log2TE_' + stage]
    d['log2_daft'] = np.log2(d[STAGE_DAFT[stage]] + 0.01)
    d['log2_ss'] = np.log2(te['ss_CPM_' + stage] + 0.1)
    d['log2_cds'] = np.log2(d.cds_len)
    for a in ['ref', 'v71']:
        d[a + '_log2_utr5'] = np.log2(d[a + '_utr5_len'])
        d[a + '_log2_utr3'] = np.log2(d[a + '_utr3_len'].where(d[a + '_utr3_len'] > 0))
        d[a + '_log1p_uORF'] = np.log1p(d[a + '_n_uORF'])
        d[a + '_any_uORF'] = (d[a + '_n_uORF'] > 0).astype(float)
        d[a + '_log1p_uAUG'] = np.log1p(d[a + '_n_uAUG'])
        d[a + '_uAUG_per100'] = d[a + '_n_uAUG'] / d[a + '_utr5_len'] * 100
        d[a + '_oORF'] = (d[a + '_n_oORF_overlapping'] > 0).astype(float)
        d[a + '_log1p_strong'] = np.log1p(d[a + '_n_uORF_strongctx'])
    return d

FEATS = ['log2_utr5', 'log1p_uORF', 'any_uORF', 'log1p_uAUG', 'uAUG_per100', 'oORF', 'log1p_strong', 'log2_utr3']
COVS = ['log2_daft', 'log2_cds', 'cds_GC']

def zs(x):
    return (x - x.mean()) / x.std()

def fit(d, y, xs, covs):
    sub = d[[y] + xs + covs].dropna().copy()
    for c in xs + covs:
        if sub[c].nunique() > 2:
            sub[c] = zs(sub[c])
    f = y + ' ~ ' + ' + '.join(xs + covs)
    m = smf.ols(f, data=sub).fit(cov_type='HC3')
    m0 = smf.ols(y + ' ~ ' + ' + '.join(covs), data=sub).fit() if covs else None
    return m, m0, sub

def partial_r2(m, m0):
    return (m0.ssr - m.ssr) / m0.ssr

def run_all(cov, te, te_label, stages=('0h', '18h', '36h'), nboot=1000, seed=1):
    rng = np.random.default_rng(seed)
    rows = []
    strat = []
    for st in stages:
        d = build_design(cov, te, st)
        base = d.TE.notna() & d.log2_daft.notna()
        sets = {
            'common_all': base & d.ref_utr5_len.notna() & d.v71_utr5_len.notna(),
            'common_reference': base & d.ref_utr5_len.notna() & d.v71_utr5_len.notna() & (d.lr_group == 'confirmed'),
            'common_model': base & d.ref_utr5_len.notna() & d.v71_utr5_len.notna() & (d.final_src == 'model'),
            'refined_all': base & d.ref_utr5_len.notna(),
        }
        for sname, mask in sets.items():
            dd = d[mask]
            anns = ['ref'] if sname == 'refined_all' else ['ref', 'v71']
            for feat in FEATS:
                fsets = {}
                for a in anns:
                    x = a + '_' + feat
                    if feat == 'log2_utr3':
                        sub_mask = dd['ref_log2_utr3'].notna() & dd['v71_log2_utr3'].notna() if len(anns) == 2 else dd[x].notna()
                    else:
                        sub_mask = dd[x].notna()
                    fsets[a] = sub_mask
                common_mask = np.logical_and.reduce([fsets[a] for a in anns])
                ddd = dd[common_mask]
                res = {}
                for a in anns:
                    x = a + '_' + feat
                    m, m0, sub = fit(ddd, 'TE', [x], COVS)
                    mfull = smf.ols('TE ~ ' + ' + '.join([x] + COVS), data=sub).fit()
                    pr2 = partial_r2(mfull, m0)
                    ci = m.conf_int().loc[x]
                    rho, prho = ss.spearmanr(sub[x], sub.TE)
                    res[a] = pr2
                    rows.append(dict(te_source=te_label, stage=st, gene_set=sname, annotation=a, model='single_feature+covariates',
                                     feature=feat, n=len(sub), beta_std=m.params[x], ci_low=ci[0], ci_high=ci[1], p=m.pvalues[x],
                                     partial_R2=pr2, spearman_rho_unadj=rho, spearman_p_unadj=prho))
                if len(anns) == 2:
                    # encompassing model: both annotations' feature together
                    xr, xv = 'ref_' + feat, 'v71_' + feat
                    m, m0, sub = fit(ddd, 'TE', [xr, xv], COVS)
                    # bootstrap delta partial R2 (ref - v71)
                    deltas = []
                    arr = sub.reset_index(drop=True)
                    for b in range(nboot):
                        bs = arr.iloc[rng.integers(0, len(arr), len(arr))]
                        m0b = smf.ols('TE ~ ' + ' + '.join(COVS), data=bs).fit()
                        r_ = partial_r2(smf.ols('TE ~ ' + ' + '.join([xr] + COVS), data=bs).fit(), m0b)
                        v_ = partial_r2(smf.ols('TE ~ ' + ' + '.join([xv] + COVS), data=bs).fit(), m0b)
                        deltas.append(r_ - v_)
                    deltas = np.array(deltas)
                    rows.append(dict(te_source=te_label, stage=st, gene_set=sname, annotation='ref_vs_v71', model='encompassing(ref+v71)+covariates',
                                     feature=feat, n=len(sub), beta_std=m.params[xr], ci_low=m.conf_int().loc[xr][0], ci_high=m.conf_int().loc[xr][1], p=m.pvalues[xr],
                                     beta_std_v71=m.params[xv], p_v71=m.pvalues[xv],
                                     delta_partial_R2_ref_minus_v71=res['ref'] - res['v71'],
                                     delta_ci_low=np.percentile(deltas, 2.5), delta_ci_high=np.percentile(deltas, 97.5),
                                     delta_boot_p=2 * min((deltas <= 0).mean(), (deltas >= 0).mean())))
            # joint multivariable model
            for a in anns:
                xs = [a + '_log2_utr5', a + '_log1p_uORF', a + '_oORF', a + '_log2_utr3']
                sub_d = dd.dropna(subset=['ref_log2_utr3', 'v71_log2_utr3'] if len(anns) == 2 else [a + '_log2_utr3'])
                m, m0, sub = fit(sub_d, 'TE', xs, COVS)
                full = smf.ols('TE ~ ' + ' + '.join(xs + COVS), data=sub).fit()
                for x in xs:
                    red = smf.ols('TE ~ ' + ' + '.join([z for z in xs if z != x] + COVS), data=sub).fit()
                    ci = m.conf_int().loc[x]
                    rows.append(dict(te_source=te_label, stage=st, gene_set=sname, annotation=a, model='joint(utr5+uORF+oORF+utr3)+covariates',
                                     feature=x.split('_', 1)[1], n=len(sub), beta_std=m.params[x], ci_low=ci[0], ci_high=ci[1], p=m.pvalues[x],
                                     partial_R2=partial_r2(full, red), model_R2=full.rsquared))
                for cv in COVS:
                    rows.append(dict(te_source=te_label, stage=st, gene_set=sname, annotation=a, model='joint(utr5+uORF+oORF+utr3)+covariates',
                                     feature='covariate:' + cv, n=len(sub), beta_std=m.params[cv], ci_low=m.conf_int().loc[cv][0], ci_high=m.conf_int().loc[cv][1], p=m.pvalues[cv], model_R2=full.rsquared))
            # expression-quintile stratification (Qd = DAFT_max quintile over all genes)
            if sname in ('common_all', 'common_reference'):
                for q, dq in dd.groupby('Qd'):
                    for a in anns:
                        for feat in ['log2_utr5', 'log1p_uORF', 'log1p_uAUG', 'log2_utr3']:
                            x = a + '_' + feat
                            s2 = dq[[x, 'TE']].dropna()
                            if len(s2) < 20: continue
                            rho, p = ss.spearmanr(s2[x], s2.TE)
                            strat.append(dict(te_source=te_label, stage=st, gene_set=sname, Qd=q, annotation=a, feature=feat, n=len(s2), spearman_rho=rho, p=p))
    res = pd.DataFrame(rows); strat = pd.DataFrame(strat)
    # BH within te_source x model family
    for df_, key in [(res, 'p')]:
        df_['q_BH'] = np.nan
        ok = df_[key].notna()
        df_.loc[ok, 'q_BH'] = multipletests(df_.loc[ok, key], method='fdr_bh')[1]
    if len(strat):
        strat['q_BH'] = multipletests(strat.p, method='fdr_bh')[1]
    return res, strat
