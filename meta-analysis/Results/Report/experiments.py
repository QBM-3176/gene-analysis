"""Experiments for the ML report. Same LODO design as the notebooks (Logistic Regression, PyCaret defaults,
mean/mode imputation, one-hot, no scaling; GSE104948 + GSE104954 merged into one ERCB group)."""
import json, warnings
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import OneHotEncoder
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.metrics import roc_auc_score, roc_curve
warnings.filterwarnings('ignore')

ROOT = Path('/Users/striki/Striki/Quant/quant.2/MetaAnalysis')
ML = ROOT / 'Results/ML_training'
OUT = Path(__file__).parent / 'exp'; OUT.mkdir(exist_ok=True)
R = {}

d2 = pd.read_excel(ROOT / 'Data/ML_input/all2.xlsx')
d1 = pd.read_excel(ROOT / 'Data/ML_input/all.xlsx')
assert (d1['GSM'].to_numpy() == d2['GSM'].to_numpy()).all()
META = ['GSM', 'StudyID', 'Batch', 'Compartment', 'Condition', 'Disease_Label', 'Indication', 'Sex']
GENES = [c for c in d2.columns if c not in META]
RRG = ['C4BPA', 'COL13A1', 'HABP2', 'CFH', 'CFI']
y = d2['Disease_Label'].to_numpy()
grp = d2['StudyID'].replace({'GSE104948-GPL22945': 'ERCB', 'GSE104954-GPL22945': 'ERCB'}).to_numpy()
SHORT = {'GSE104948-GPL22945': 'GSE104948', 'GSE104954-GPL22945': 'GSE104954', 'GSE99339-GPL19184': 'GSE99339'}
study = d2['StudyID'].replace(SHORT).to_numpy()
CTYPE = {'GSE104948': 'Living donor + tumour nephrectomy', 'GSE104954': 'Living donor + tumour nephrectomy',
         'GSE108109': 'Living donor', 'GSE200828': 'Living donor', 'GSE116626': 'Living donor',
         'GSE99339': 'Tumour nephrectomy', 'GSE129973': 'Tumour nephrectomy', 'GSE96804': 'Tumour nephrectomy',
         'GSE121233': 'Tumour nephrectomy', 'GSE30122': 'Not stated', 'GSE66494': 'Reference RNA'}
IEG = ['FOS', 'EGR1', 'DUSP1', 'ATF3', 'NR4A2', 'JUNB', 'RGS2', 'ZFP36', 'MAFF', 'NR4A1', 'NFIL3', 'FOSB', 'JUN',
       'PTGS2', 'CCN1', 'HBEGF', 'GDF15']
FASTING = ['ALB', 'G6PC1', 'PCK1', 'PDK4', 'IGF1', 'APOH', 'GSTA1']
MYELOID = ['S100A12', 'MNDA', 'FPR1', 'S100A8']
PROC = IEG + FASTING + MYELOID


def lodo(df, genes, cats=('Compartment', 'Sex')):
    cats = list(cats)
    tr = []
    if genes: tr.append(('num', SimpleImputer(strategy='mean'), genes))
    if cats: tr.append(('cat', Pipeline([('i', SimpleImputer(strategy='most_frequent')),
                                         ('o', OneHotEncoder(handle_unknown='ignore'))]), cats))
    pipe = Pipeline([('pre', ColumnTransformer(tr)), ('m', LogisticRegression(max_iter=1000, random_state=42))])
    X = df[genes + cats]
    score = np.empty(len(df))
    for a, b in LeaveOneGroupOut().split(X, y, grp):
        pipe.fit(X.iloc[a], y[a]); score[b] = pipe.predict_proba(X.iloc[b])[:, 1]
    return score


def summarise(score):
    pred = (score >= 0.5).astype(int)
    rows = []
    for g in pd.unique(grp):
        m = grp == g; yt, p, s = y[m], pred[m], score[m]
        tp = int(((p == 1) & (yt == 1)).sum()); tn = int(((p == 0) & (yt == 0)).sum())
        fp = int(((p == 1) & (yt == 0)).sum()); fn = int(((p == 0) & (yt == 1)).sum())
        fpr, tpr, _ = roc_curve(yt, s)
        rows.append({'Cohort': g, 'n': int(m.sum()), 'Controls': int((yt == 0).sum()), 'Cases': int((yt == 1).sum()),
                     'AUC': roc_auc_score(yt, s), 'Sens': tp / (tp + fn), 'Spec': tn / (tn + fp) if tn + fp else np.nan,
                     'PPV': tp / (tp + fp) if tp + fp else np.nan, 'NPV': tn / (tn + fn) if tn + fn else np.nan,
                     'Spec@Sens0.90': float(1 - fpr[np.argmax(tpr >= 0.9)])})
    t = pd.DataFrame(rows)
    tp = int(((pred == 1) & (y == 1)).sum()); tn = int(((pred == 0) & (y == 0)).sum())
    fp = int(((pred == 1) & (y == 0)).sum()); fn = int(((pred == 0) & (y == 1)).sum())
    pooled = {'TP': tp, 'FN': fn, 'TN': tn, 'FP': fp, 'Sens': tp / (tp + fn), 'Spec': tn / (tn + fp),
              'PPV': tp / (tp + fp), 'NPV': tn / (tn + fn) if tn + fn else None,
              'meanAUC': t['AUC'].mean(), 'meanSens': t['Sens'].mean(), 'meanSpec': t['Spec'].mean(),
              'meanSpec@Sens0.90': t['Spec@Sens0.90'].mean()}
    return t, pooled, pred


# ---------- E1: per-cohort performance of the three best models ----------
top50 = [g for g in pd.read_excel(ML / 'ALL/TOP_50_by_AUC/ML_TOP_50_by_AUC_FEAT_IMP_by_AUC.xlsx')['Source'] if g in GENES]
assert len(top50) == 50
scores = {'RR': lodo(d1, RRG), 'ALL': lodo(d2, GENES), 'TOP50': lodo(d2, top50)}
E1, preds = {}, {}
for k, s in scores.items():
    t, pooled, pred = summarise(s); preds[k] = pred
    E1[k] = {'per_cohort': t.round(3).to_dict('records'), 'pooled': {a: (round(b, 3) if isinstance(b, float) else b) for a, b in pooled.items()}}
R['E1'] = E1
print('E1', {k: v['pooled'] for k, v in E1.items()}, flush=True)

# ---------- E2: detection rate per indication, specificity per control source ----------
ind = d2['Indication'].to_numpy()
rows = []
for i in sorted(set(ind) - {'Control'}):
    m = ind == i
    rows.append({'Indication': i, 'n': int(m.sum()), **{k: round(float(preds[k][m].mean()), 3) for k in preds}})
R['E2_indication'] = rows
rows = []
ct = np.array([CTYPE[s] for s in study])
for c in pd.unique(ct):
    m = (ct == c) & (y == 0)
    rows.append({'Control source': c, 'n': int(m.sum()), **{k: round(float(1 - preds[k][m].mean()), 3) for k in preds}})
R['E2_controls'] = rows
rows = []
for c in ['Glomerular', 'Tubulointerstitial', 'Whole biopsy']:
    m = d2['Compartment'].to_numpy() == c
    rows.append({'Compartment': c, 'n': int(m.sum()), 'Controls': int((y[m] == 0).sum()), 'Cases': int((y[m] == 1).sum()),
                 **{f'{k} sens': round(float(preds[k][m & (y == 1)].mean()), 3) for k in preds},
                 **{f'{k} spec': round(float(1 - preds[k][m & (y == 0)].mean()), 3) for k in preds}})
R['E2_compartment'] = rows

# ---------- E3: what the metadata and the panel contribute ----------
E3 = []
for name, df, genes, cats in [('Compartment + Sex only (no genes)', d1, [], ['Compartment', 'Sex']),
                              ('5 RenoRisk genes only', d1, RRG, []),
                              ('5 RenoRisk genes + Compartment + Sex (RR notebook)', d1, RRG, ['Compartment', 'Sex']),
                              ('All genes without Compartment, Sex', d2, GENES, [])]:
    t, pooled, _ = summarise(lodo(df, genes, cats))
    E3.append({'Features': name, 'AUC': round(pooled['meanAUC'], 3), 'Sens': round(pooled['Sens'], 3), 'Spec': round(pooled['Spec'], 3)})
    print('E3', E3[-1], flush=True)
R['E3'] = E3

# ---------- E4: within-cohort effect sizes (Hedges g) of the panel and marker genes ----------
def hedges(gene):
    out = {}
    for s in pd.unique(study):
        m = study == s; c = d2.loc[m & (y == 0), gene]; x = d2.loc[m & (y == 1), gene]
        sp = np.sqrt(((len(c) - 1) * c.var() + (len(x) - 1) * x.var()) / (len(c) + len(x) - 2))
        out[s] = round(float((x.mean() - c.mean()) / sp), 2) if sp > 0 else None
    return out
R['E4_g'] = {g: hedges(g) for g in RRG + ['FOS', 'EGR1', 'DUSP1', 'ATF3', 'PDK4', 'PCK1', 'G6PC1', 'FN1', 'LUM', 'S100A12', 'HLA-DRB1']}
R['ctype'] = CTYPE

# ---------- E5: stress score and the procurement experiments ----------
z = d2[GENES].groupby(study).transform(lambda v: (v - v.mean()) / (v.std(ddof=0) + 1e-9))
stress = z[IEG].mean(axis=1).to_numpy()
rows = []
for s in pd.unique(study):
    m = study == s
    rows.append({'Series': s, 'Control source': CTYPE[s], 'Controls': round(float(np.mean(stress[m & (y == 0)])), 2),
                 'Cases': round(float(np.mean(stress[m & (y == 1)])), 2),
                 'AUC of stress score alone': round(float(roc_auc_score(y[m], -stress[m])), 3)})
R['E5_stress_by_series'] = rows
# Genes co-varying with the stress score inside each series AND label group (label effect removed)
key = pd.Series(study).astype(str) + '_' + pd.Series(y).astype(str)
zc = d2[GENES].groupby(key.to_numpy()).transform(lambda v: (v - v.mean()) / (v.std(ddof=0) + 1e-9))
sc = pd.Series(stress).groupby(key.to_numpy()).transform(lambda v: (v - v.mean()) / (v.std(ddof=0) + 1e-9)).to_numpy()
r = (zc.mul(sc, axis=0)).mean().fillna(0)
R['E5_n_corr'] = {f'|r|>{t}': int((r.abs() > t).sum()) for t in [0.3, 0.4, 0.5]}
R['E5_top_corr'] = r.sort_values(ascending=False).head(25).round(2).to_dict()
E5 = []
for name, genes in [('17 immediate-early genes only', IEG),
                    ('28 procurement genes only (IEG + fasting + blood)', PROC),
                    ('All genes minus 17 immediate-early genes', [g for g in GENES if g not in IEG]),
                    ('All genes minus 28 procurement genes', [g for g in GENES if g not in PROC]),
                    (f'All genes minus {int((r.abs() > 0.4).sum())} stress-co-varying genes (|r| > 0.4)', [g for g in GENES if abs(r[g]) <= 0.4]),
                    (f'All genes minus {int((r.abs() > 0.3).sum())} stress-co-varying genes (|r| > 0.3)', [g for g in GENES if abs(r[g]) <= 0.3])]:
    t, pooled, _ = summarise(lodo(d2, genes))
    E5.append({'Features': name, 'n genes': len(genes), 'AUC': round(pooled['meanAUC'], 3),
               'Sens': round(pooled['Sens'], 3), 'Spec': round(pooled['Spec'], 3)})
    print('E5', E5[-1], flush=True)
R['E5'] = E5
# Correlation of the RenoRisk genes with the stress score (within series and label)
R['E5_RR_r'] = {g: round(float(r[g]), 2) for g in RRG}

# ---------- E6: are the two TOP_50 folders identical? ----------
same = {}
for f in ['perfs', 'FEAT_IMP_by_AUC', 'FEAT_IMP_by_F1']:
    a = pd.read_excel(ML / f'ALL/TOP_50_by_AUC/ML_TOP_50_by_AUC_{f}.xlsx'); b = pd.read_excel(ML / f'ALL/TOP_50_by_F1/ML_TOP_50_by_F1_{f}.xlsx')
    same[f] = bool(a.equals(b))
R['E6_same'] = same

# ---------- figures ----------
plt.rcParams.update({'font.size': 13})
t = {k: pd.DataFrame(v['per_cohort']).set_index('Cohort') for k, v in E1.items()}
order = t['ALL'].sort_values('AUC').index
fig, ax = plt.subplots(figsize=(10, 5.5)); w = 0.27; xs = np.arange(len(order))
for i, (k, lab, col) in enumerate([('RR', '5 RenoRisk genes', '#9DB4D3'), ('ALL', 'All 10,422 genes', '#1F4E9A'), ('TOP50', 'Top 50 genes', '#E08A2E')]):
    ax.bar(xs + (i - 1) * w, t[k].loc[order, 'AUC'], w, label=lab, color=col)
ax.axhline(0.5, ls='--', color='grey', lw=1); ax.set_xticks(xs); ax.set_xticklabels(order, rotation=35, ha='right')
ax.set_ylabel('ROC-AUC in the left-out cohort'); ax.set_ylim(0, 1.05); ax.legend(loc='lower right', fontsize=11)
ax.set_title('Per-cohort ROC-AUC, Logistic Regression, leave-one-dataset-out')
fig.tight_layout(); fig.savefig(OUT / 'fig_per_cohort_auc.png', dpi=180); plt.close(fig)

fig, ax = plt.subplots(figsize=(10, 5.5)); ser = [r_['Series'] for r_ in sorted(rows, key=lambda q: q['Series'])]
ser = sorted(pd.unique(study), key=lambda s: (CTYPE[s], s))
for i, s in enumerate(ser):
    m = study == s
    for lab, off, col in [(0, -0.18, '#2E7D32'), (1, 0.18, '#C62828')]:
        v = stress[m & (y == lab)]
        ax.boxplot(v, positions=[i + off], widths=0.3, patch_artist=True, showfliers=False,
                   boxprops=dict(facecolor=col, alpha=0.55), medianprops=dict(color='black'))
ax.set_xticks(range(len(ser))); ax.set_xticklabels([f'{s}\n({CTYPE[s].replace("Living donor + tumour nephrectomy", "LD + TN").replace("Tumour nephrectomy", "TN").replace("Living donor", "LD")})' for s in ser], rotation=40, ha='right', fontsize=10)
ax.axhline(0, color='grey', lw=0.8)
ax.set_ylabel('Stress score (mean z of 17 IEGs)'); ax.set_title('Immediate-early-gene stress score: controls (green) vs diseased (red)')
fig.tight_layout(); fig.savefig(OUT / 'fig_stress_score.png', dpi=180); plt.close(fig)

json.dump(R, open(OUT / 'results.json', 'w'), indent=1, default=float)
print('done')
