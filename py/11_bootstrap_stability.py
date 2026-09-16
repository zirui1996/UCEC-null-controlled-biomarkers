#!/usr/bin/env python3
"""11_bootstrap_stability.py -- protocol section 13
ML selection stability conditional on the recurrent consensus feature set
1000 bootstrap iterations, resampling within the six dataset-by-group strata (section 13.1.2)
LASSO and SVM-RFE are rerun on every iteration and the genes entering the intersection are recorded (section 13.1.3)
degenerate iterations are monitored and reported (section 13.1.4)
"""
import numpy as np, pandas as pd, json, os, sys, warnings, time
from collections import Counter
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import RFECV
from sklearn.svm import LinearSVC
warnings.filterwarnings("ignore")

os.chdir(os.environ.get("UCEC_PROJ", os.getcwd()))
FEATSET = sys.argv[1] if len(sys.argv) > 1 else "recurrent"   # recurrent | noCellCycle
NBOOT   = int(sys.argv[2]) if len(sys.argv) > 2 else 1000

d = pd.read_csv("results/ML/ML_input_matrix.tsv", sep="\t")
genes = [c for c in d.columns if c not in ("sample","dataset","group")]
if FEATSET == "noCellCycle":
    keep = set(pd.read_csv("results/consensus/recurrent_noCellCycle.tsv", sep="\t")["gene"])
    genes = [g for g in genes if g in keep]
X_all = d[genes].values
y_all = (d["group"] == "Tumor").astype(int).values
meta  = d[["dataset","group"]].copy()
print(f"feature set = {FEATSET} ({len(genes)} genes) | samples {len(d)} | {NBOOT} bootstrap iterations", flush=True)

# the optimal C from section 12, held fixed across bootstrap iterations rather than retuned each time (otherwise the result is not comparable with the section 12 candidate set)
C_FIXED = json.load(open("results/ML/candidates.json"))["lasso_C"]
rng = np.random.default_rng(42)
strata = [g.index.values for _, g in meta.groupby(["dataset","group"], sort=True)]
strata_names = [f"{k[0]}-{k[1]}" for k, _ in meta.groupby(["dataset","group"], sort=True)]
print("six strata:", dict(zip(strata_names, [len(s) for s in strata])), flush=True)

def make_lasso():
    return Pipeline([("imputer", SimpleImputer(strategy="median")),
                     ("scale", StandardScaler()),
                     ("model", LogisticRegression(C=C_FIXED, penalty="l1", solver="liblinear",
                                                  max_iter=20000, random_state=42))])
def make_svm():
    return Pipeline([("imputer", SimpleImputer(strategy="median")),
                     ("scale", StandardScaler()),
                     ("svm", LinearSVC(random_state=42, max_iter=20000, dual="auto"))])

cnt_lasso, cnt_svm, cnt_inter = Counter(), Counter(), Counter()
n_ok = n_failed = 0
degenerate = Counter()          # section 13.1.4: number of iterations in which a stratum held too few unique samples
t0 = time.time()
for b in range(NBOOT):
    idx = np.concatenate([rng.choice(s, size=len(s), replace=True) for s in strata])
    for nm, s in zip(strata_names, strata):
        u = len(np.unique(idx[np.isin(idx, s)]))
        if u <= 2: degenerate[nm] += 1
    Xb, yb = X_all[idx], y_all[idx]
    if len(np.unique(yb)) < 2:
        n_failed += 1; continue
    try:
        lp = make_lasso().fit(Xb, yb)
        sel_l = {g for g, c in zip(genes, lp.named_steps["model"].coef_.ravel()) if c != 0}
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        sel = RFECV(estimator=make_svm(), step=0.1, cv=cv, scoring="f1_macro",
                    importance_getter="named_steps.svm.coef_", n_jobs=-1).fit(Xb, yb)
        sel_s = {g for g, k in zip(genes, sel.support_) if k}
    except Exception:
        n_failed += 1; continue
    n_ok += 1
    cnt_lasso.update(sel_l); cnt_svm.update(sel_s); cnt_inter.update(sel_l & sel_s)
    if (b+1) % 50 == 0:
        print(f"  ...{b+1}/{NBOOT}  ok {n_ok} failed {n_failed}  elapsed {time.time()-t0:.0f}s", flush=True)

print(f"\nsuccessful iterations {n_ok} / {NBOOT} (failed {n_failed}, {100*n_failed/NBOOT:.1f}%)", flush=True)
print("section 13.1.4 degeneracy monitor (iterations in which a stratum held <=2 unique samples):",
      {k: degenerate.get(k, 0) for k in strata_names}, flush=True)

res = pd.DataFrame({"gene": genes})
res["freq_intersection"] = [cnt_inter.get(g,0)/max(n_ok,1) for g in genes]
res["freq_lasso"]        = [cnt_lasso.get(g,0)/max(n_ok,1) for g in genes]
res["freq_svmrfe"]       = [cnt_svm.get(g,0)/max(n_ok,1) for g in genes]

# the prespecified sort key of section 13.2
rec = pd.read_csv("results/consensus/recurrent_consensus.tsv", sep="\t")
lfc = rec.set_index("gene")[["GSE17025_logFC","GSE63678_logFC","GSE106191_logFC"]]
res["median_abs_logFC"] = [np.median(np.abs(lfc.loc[g].values)) if g in lfc.index else np.nan for g in res.gene]
res = res.sort_values(["freq_intersection","median_abs_logFC","gene"],
                      ascending=[False, False, True]).reset_index(drop=True)
res.to_csv(f"results/ML/bootstrap_stability_{FEATSET}.tsv", sep="\t", index=False)

cand = set(json.load(open("results/ML/candidates.json"))["candidate"])
print(f"\n=== top 15 (in the prespecified order of section 13.2) ===", flush=True)
for i, r in res.head(15).iterrows():
    print(f"{i+1:3d}. {r.gene:11s} intersection frequency {r.freq_intersection:.3f} | LASSO {r.freq_lasso:.3f} "
          f"| SVM-RFE {r.freq_svmrfe:.3f} | median|logFC| {r.median_abs_logFC:.2f}"
          f"{'  ← §12 candidate' if r.gene in cand else ''}", flush=True)

sub = res[res.gene.isin(cand)].reset_index(drop=True)
print(f"\n=== stability ranking of the section 12 candidate set (LASSO n SVM-RFE) ===", flush=True)
for i, r in sub.iterrows():
    print(f"{i+1:3d}. {r.gene:11s} {r.freq_intersection:.3f}   median|logFC| {r.median_abs_logFC:.2f}", flush=True)
gx = sub.iloc[0]
print(f"\n>>> Gene X = {gx.gene} (intersection frequency {gx.freq_intersection:.3f}, median|log2FC| {gx.median_abs_logFC:.2f})", flush=True)
json.dump({"gene_x": gx.gene, "freq": float(gx.freq_intersection), "n_ok": n_ok, "n_failed": n_failed,
           "featset": FEATSET, "degenerate": dict(degenerate)},
          open(f"results/ML/gene_x_{FEATSET}.json","w"), indent=2)
print("===== 11_bootstrap DONE =====", flush=True)
