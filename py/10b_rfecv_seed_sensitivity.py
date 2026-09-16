#!/usr/bin/env python3
"""10b_rfecv_seed_sensitivity.py -- sensitivity of the optimal SVM-RFE (RFECV) feature count to the cross-validation split
Exactly the same Pipeline, step and scoring as py/10_lasso_svmrfe.py; only the random_state of StratifiedKFold changes.
Purpose: to support the statement that a single RFECV feature count is unstable and that bootstrap frequency is therefore used instead (v2.6).
"""
import numpy as np, pandas as pd, os, warnings
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import RFECV
from sklearn.svm import LinearSVC
warnings.filterwarnings("ignore")
os.chdir(os.environ.get("UCEC_PROJ", os.getcwd()))
SEEDS = [1, 2, 3, 7, 11, 21, 42, 99, 123, 2024]

d = pd.read_csv("results/ML/ML_input_matrix.tsv", sep="\t")
genes = [c for c in d.columns if c not in ("sample", "dataset", "group")]
X = d[genes].values; y = (d["group"] == "Tumor").astype(int).values
print(f"X = {X.shape[0]} samples x {X.shape[1]} genes")

rows, curves = [], []
for sd in SEEDS:
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=sd)
    pipe = Pipeline([("imputer", SimpleImputer(strategy="median")), ("scale", StandardScaler()),
                     ("svm", LinearSVC(random_state=42, max_iter=20000, dual="auto"))])
    s = RFECV(estimator=pipe, step=0.05, cv=cv, scoring="f1_macro",
              importance_getter="named_steps.svm.coef_", n_jobs=-1).fit(X, y)
    nf = s.cv_results_["n_features"]; m = s.cv_results_["mean_test_score"]; sdv = s.cv_results_["std_test_score"]
    i = int(np.argmax(m)); se = sdv[i] / np.sqrt(5)
    within = nf[m >= m.max() - se]
    rows.append(dict(cv_seed=sd, optimal_n_features=int(s.n_features_), best_f1=round(float(m.max()), 4),
                     f1_at_all_features=round(float(m[np.argmax(nf)]), 4),
                     n_features_within_1SE_min=int(within.min()), n_features_within_1SE_max=int(within.max())))
    curves += [dict(cv_seed=sd, n_features=int(a), mean_f1=float(b), sd_f1=float(c)) for a, b, c in zip(nf, m, sdv)]
    print(f"  seed {sd:>4}: optimum {s.n_features_:>3} | best f1 {m.max():.4f} | features within 1SE {within.min()}–{within.max()}")

res = pd.DataFrame(rows)
res.to_csv("results/ML/rfecv_seed_sensitivity.tsv", sep="\t", index=False)
pd.DataFrame(curves).to_csv("results/ML/rfecv_seed_curves.tsv", sep="\t", index=False)
vc = res.optimal_n_features.value_counts()
print(f"\noptimal feature count ranges from {res.optimal_n_features.min()} to {res.optimal_n_features.max()}; mode {vc.index[0]} ({vc.iloc[0]}/{len(SEEDS)})")
print("===== 10b DONE =====")
