#!/usr/bin/env python3
"""10_lasso_svmrfe.py -- protocol section 12: feature selection with LASSO and SVM-RFE
Input X = the recurrent consensus genes (v2.6: 182), pooled across the three cohorts after z-scoring within each dataset
the imputer and scaler are inside the Pipeline, so no information leaks across cross-validation folds.
"""
import numpy as np, pandas as pd, json, os
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import RFECV
from sklearn.svm import LinearSVC

os.chdir(os.environ.get("UCEC_PROJ", os.getcwd()))
d = pd.read_csv("results/ML/ML_input_matrix.tsv", sep="\t")
genes = [c for c in d.columns if c not in ("sample","dataset","group")]
X = d[genes].values
y = (d["group"] == "Tumor").astype(int).values
print(f"X = {X.shape[0]} samples x {X.shape[1]} genes | Tumour {y.sum()} / Normal {(1-y).sum()}")
print("per cohort:", d.groupby(["dataset","group"]).size().to_dict())

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# ---------- LASSO ----------
lasso_pipe = Pipeline([
    ("imputer", SimpleImputer(strategy="median")),
    ("scale",   StandardScaler()),
    ("model",   LogisticRegressionCV(
        Cs=np.logspace(-4, 4, 60), cv=cv, penalty="l1",
        solver="liblinear", scoring="neg_log_loss",
        max_iter=20000, random_state=42)),
])
lasso_pipe.fit(X, y)
mdl = lasso_pipe.named_steps["model"]
coef = mdl.coef_.ravel()
sel_lasso = [g for g, c in zip(genes, coef) if c != 0]
print(f"\n[LASSO] optimal C = {mdl.C_[0]:.5g} (lambda = {1/mdl.C_[0]:.5g})")
print(f"[LASSO] genes with non-zero coefficients: {len(sel_lasso)} / {len(genes)}")

pd.DataFrame({"gene": genes, "coef": coef}).query("coef != 0") \
  .assign(abs_coef=lambda t: t.coef.abs()).sort_values("abs_coef", ascending=False) \
  .to_csv("results/ML/lasso_selected_genes.tsv", sep="\t", index=False)

# ---------- SVM-RFE ----------
svm_pipe = Pipeline([
    ("imputer", SimpleImputer(strategy="median")),
    ("scale",   StandardScaler()),
    ("svm",     LinearSVC(random_state=42, max_iter=20000, dual="auto")),
])
selector = RFECV(estimator=svm_pipe, step=0.05, cv=cv, scoring="f1_macro",
                 importance_getter="named_steps.svm.coef_", n_jobs=-1)
selector.fit(X, y)
# section 12.2 self-check: estimator_ must be a Pipeline, not a bare LinearSVC
print(f"\n[self-check] type of selector.estimator_ = {type(selector.estimator_).__name__} "
      f"-> {'pass' if isinstance(selector.estimator_, Pipeline) else 'fail (the v2.0 bug reproduced)'}")
sel_svm = [g for g, s in zip(genes, selector.support_) if s]
print(f"[SVM-RFE] optimal number of features = {selector.n_features_}, best f1_macro = {selector.cv_results_['mean_test_score'].max():.4f}")

pd.DataFrame({"gene": genes, "rank": selector.ranking_, "selected": selector.support_}) \
  .sort_values("rank").to_csv("results/ML/svm_rfe_selected_genes.tsv", sep="\t", index=False)
np.save("results/ML/rfecv_grid_scores.npy", selector.cv_results_["mean_test_score"])

# ---------- intersection ----------
cand = sorted(set(sel_lasso) & set(sel_svm))
print(f"\n=== section 12.3 final candidate set ===")
print(f"LASSO {len(sel_lasso)} ∩ SVM-RFE {len(sel_svm)} = {len(cand)}")
print(cand)

core = ["HAND2-AS1","OSR2","PEG3","RUNX1T1"]
print(f"\nselection status of the four core consensus genes:")
for g in core:
    print(f"  {g:11s} LASSO {'✓' if g in sel_lasso else '✗'} | SVM-RFE {'✓' if g in sel_svm else '✗'} "
          f"| candidate set {'✓' if g in cand else '✗'}")

json.dump({"lasso": sel_lasso, "svm_rfe": sel_svm, "candidate": cand,
           "lasso_C": float(mdl.C_[0]), "rfecv_n": int(selector.n_features_)},
          open("results/ML/candidates.json","w"), indent=2)
pd.DataFrame({"gene": cand}).to_csv("results/ML/candidate_genes.tsv", sep="\t", index=False)
print("\n===== 10_lasso_svmrfe DONE =====")
