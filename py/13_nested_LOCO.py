#!/usr/bin/env python3
"""13_nested_LOCO.py -- protocol section 14: nested leave-one-dataset-out
the entire feature-selection procedure (fold-internal consensus -> LASSO -> SVM-RFE -> fitting) is carried out inside the outer fold.
The test cohort is only z-scored within itself as described in section 14.3; it takes no part in the DEG analysis, feature selection or tuning.

Key point (section 14.1): limma is run independently per cohort anyway, so reusing the per-cohort DEG tables does not leak information --
the GSE63678 DEG result does not depend on whether GSE17025 is present. What the fold decides is which two tables are intersected.
"""
import numpy as np, pandas as pd, json, os, warnings
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import RFECV
from sklearn.svm import LinearSVC
from sklearn.metrics import roc_auc_score, roc_curve
warnings.filterwarnings("ignore")
os.chdir(os.environ.get("UCEC_PROJ", os.getcwd()))

COH = ["GSE17025", "GSE63678", "GSE106191"]
FDR, LFC = 0.05, 1.0

deg = {c: pd.read_csv(f"results/DEG/{c}_all_genes.tsv", sep="\t").drop_duplicates("gene").set_index("gene")
       for c in COH}
Z = {c: pd.read_csv(f"results/LOCO/z_{c}.tsv", sep="\t", index_col=0) for c in COH}
L = {c: pd.read_csv(f"results/LOCO/lab_{c}.tsv", sep="\t").set_index("sample")["group"] for c in COH}

results, roc_data = [], {}
for test in COH:
    train = [c for c in COH if c != test]
    print(f"\n{'='*62}\nOuter Fold — Test = {test} | Train = {train[0]} + {train[1]}")

    # --- fold-internal consensus: significant in 2/2 training cohorts with the same sign (§14.1) ---
    uni = deg[train[0]].index.intersection(deg[train[1]].index)
    sig = {c: (deg[c].loc[uni, "adj.P.Val"] < FDR) & (deg[c].loc[uni, "logFC"].abs() >= LFC) for c in train}
    both = sig[train[0]] & sig[train[1]]
    same = np.sign(deg[train[0]].loc[uni, "logFC"]) == np.sign(deg[train[1]].loc[uni, "logFC"])
    fold_cons = uni[both & same]
    # v2.6: the GSE56087 direction veto was withdrawn (consistent with the global consensus rule)
    print(f"  fold-internal consensus (2/2): {len(fold_cons)}")

    if len(fold_cons) < 5:
        print("  -> fold-internal consensus < 5; recorded as not evaluable per §14.1.1")
        results.append(dict(fold=test, n_consensus=len(fold_cons), status="not evaluable")); continue

    # --- training matrix: the two training cohorts, each already z-scored within itself ---
    feats = [g for g in fold_cons if all(g in Z[c].index for c in train)]
    Xtr = pd.concat([Z[c].loc[feats].T for c in train])
    ytr = pd.concat([L[c] for c in train]).loc[Xtr.index].eq("Tumor").astype(int).values
    print(f"  training set {Xtr.shape[0]} samples × {Xtr.shape[1]} genes (Tumor {ytr.sum()} / Normal {(1-ytr).sum()})")

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    lasso = Pipeline([("imputer", SimpleImputer(strategy="median")), ("scale", StandardScaler()),
                      ("model", LogisticRegressionCV(Cs=np.logspace(-4,4,60), cv=cv, penalty="l1",
                       solver="liblinear", scoring="neg_log_loss", max_iter=20000, random_state=42))]).fit(Xtr.values, ytr)
    sel_l = {g for g,c in zip(feats, lasso.named_steps["model"].coef_.ravel()) if c != 0}
    svm = Pipeline([("imputer", SimpleImputer(strategy="median")), ("scale", StandardScaler()),
                    ("svm", LinearSVC(random_state=42, max_iter=20000, dual="auto"))])
    rfe = RFECV(estimator=svm, step=0.05, cv=cv, scoring="f1_macro",
                importance_getter="named_steps.svm.coef_", n_jobs=-1).fit(Xtr.values, ytr)
    sel_s = {g for g,k in zip(feats, rfe.support_) if k}
    cand = sorted(sel_l & sel_s)
    print(f"  LASSO {len(sel_l)} ∩ SVM-RFE {len(sel_s)} = {len(cand)} fold-internal candidate genes")

    if not cand:
        print("  -> candidate set empty; recorded as not evaluable per §14.1.1 (thresholds must not be relaxed, genes must not be borrowed from the full data)")
        results.append(dict(fold=test, n_consensus=len(fold_cons), n_candidates=0, status="not evaluable")); continue

    # --- fit the final classifier on the training set (§14.1.1: L2 logistic, C chosen by 5-fold CV within the training set) ---
    clf = Pipeline([("imputer", SimpleImputer(strategy="median")), ("scale", StandardScaler()),
                    ("model", LogisticRegressionCV(Cs=10, cv=cv, penalty="l2", solver="liblinear",
                     scoring="neg_log_loss", max_iter=20000, random_state=42))]).fit(Xtr[cand].values, ytr)

    # --- test: prediction only ---
    avail = [g for g in cand if g in Z[test].index]
    dropped = set(cand) - set(avail)
    Xte = Z[test].loc[avail].T
    yte = L[test].loc[Xte.index].eq("Tumor").astype(int).values
    if dropped:
        print(f"  [note] {len(dropped)} candidate genes are not measurable on the test cohort platform and were dropped: {sorted(dropped)}")
        clf = Pipeline([("imputer", SimpleImputer(strategy="median")), ("scale", StandardScaler()),
                        ("model", LogisticRegressionCV(Cs=10, cv=cv, penalty="l2", solver="liblinear",
                         scoring="neg_log_loss", max_iter=20000, random_state=42))]).fit(Xtr[avail].values, ytr)
    score = clf.predict_proba(Xte.values)[:, 1]
    auc = roc_auc_score(yte, score)
    fpr, tpr, _ = roc_curve(yte, score); roc_data[test] = (fpr.tolist(), tpr.tolist(), auc)
    print(f"  test set {test}: Tumor {yte.sum()} / Normal {(1-yte).sum()}  ->  AUC = {auc:.3f}")
    results.append(dict(fold=test, n_consensus=len(fold_cons), n_candidates=len(avail),
                        n_test_tumor=int(yte.sum()), n_test_normal=int((1-yte).sum()),
                        auc=round(float(auc),4), genes=";".join(avail), status="evaluable"))

print(f"\n{'='*62}\n=== section 14.5 per-fold summary (each AUC must be reported with its test-set size) ===")
df = pd.DataFrame(results)
print(df.drop(columns=["genes"], errors="ignore").to_string(index=False))
df.to_csv("results/LOCO/nested_LOCO_results.tsv", sep="\t", index=False)
json.dump(roc_data, open("results/LOCO/roc_curves.json","w"))
ev = df[df.status=="evaluable"]
if len(ev): print(f"\nevaluable folds {len(ev)}/3 | AUC range {ev.auc.min():.3f}-{ev.auc.max():.3f} | mean {ev.auc.mean():.3f}")
print("\n===== 13_nested_LOCO DONE =====")
