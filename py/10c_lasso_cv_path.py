#!/usr/bin/env python3
"""10c_lasso_cv_path.py -- source data for Figure S2 (not plotting): LASSO cross-validation curve and coefficient paths
Same Pipeline, Cs, CV and scoring as py/10_lasso_svmrfe.py; this script only exports numbers, which R plots."""
import numpy as np, pandas as pd, os, warnings
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.model_selection import StratifiedKFold
warnings.filterwarnings("ignore")
os.chdir(os.environ.get("UCEC_PROJ", os.getcwd()))
d = pd.read_csv("results/ML/ML_input_matrix.tsv", sep="\t")
genes = [c for c in d.columns if c not in ("sample", "dataset", "group")]
X = d[genes].values; y = (d["group"] == "Tumor").astype(int).values
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
Cs = np.logspace(-4, 4, 60)
pipe = Pipeline([("imputer", SimpleImputer(strategy="median")), ("scale", StandardScaler()),
                 ("model", LogisticRegressionCV(Cs=Cs, cv=cv, penalty="l1", solver="liblinear",
                                                scoring="neg_log_loss", max_iter=20000, random_state=42))]).fit(X, y)
m = pipe.named_steps["model"]
sc = m.scores_[1]                       # folds × Cs
pd.DataFrame({"C": m.Cs_, "mean_neg_logloss": sc.mean(0), "sd": sc.std(0)}).to_csv("results/ML/lasso_cv_curve.tsv", sep="\t", index=False)
Xs = StandardScaler().fit_transform(SimpleImputer(strategy="median").fit_transform(X))
rows = []
for C in Cs:
    coef = LogisticRegression(C=C, penalty="l1", solver="liblinear", max_iter=20000, random_state=42).fit(Xs, y).coef_.ravel()
    rows += [dict(C=C, gene=g, coef=c) for g, c in zip(genes, coef)]
pd.DataFrame(rows).to_csv("results/ML/lasso_coef_path.tsv", sep="\t", index=False)
print(f"best C = {m.C_[0]:.5g} | non-zero at best C = {int((m.coef_ != 0).sum())}")
