# Null-controlled evaluation of transcriptomic biomarker discovery in endometrial carcinoma

Data-analysis code for a study asking a methodological question: in tumour-versus-normal comparisons
where most of the transcriptome is differentially expressed, does a high AUC or a high machine-learning
selection frequency actually show that a gene carries information? Endometrial carcinoma is the test
case.

**This repository contains the data-analysis code only.** Data, results, tables, figures, plotting
code, logs and the manuscript are not included; they accompany the published article.

## What the pipeline does

1. **Cohort admission.** A candidate cohort is admitted only if its labels produce more differential
   expression than 1,000 label permutations with group sizes preserved.
2. **Consensus.** Differentially expressed genes are combined across admitted cohorts by vote counting
   with a same-sign requirement.
3. **Selection stability.** Genes are ranked by how often LASSO and SVM-RFE both select them across
   1,000 stratified bootstrap iterations, with a seed-sensitivity check on SVM-RFE.
4. **Null distributions.** Every discrimination metric is benchmarked against random gene sets: an
   exhaustive single-gene null, 10,000 random four-gene panels, and size-matched and
   effect-size-matched random signatures.
5. **Interpretation.** Deconvolution with stroma-adjusted partial correlation, single-nucleus RNA
   sequencing, copy number, mutation and immunohistochemistry.

## Layout

| Path | Contents |
|---|---|
| `R/` | Analysis pipeline, numbered in run order: download, preprocessing, QC and admission, differential expression, consensus, enrichment, validation, survival, microenvironment, single-nucleus, mutation and copy number |
| `py/` | LASSO and SVM-RFE selection, RFECV seed sensitivity, bootstrap stability, nested leave-one-cohort-out, Human Protein Atlas retrieval |

Scripts are numbered in execution order; the number prefix of each file matches the stage it belongs
to.

## Data

No data are redistributed here. All inputs are public and are downloaded by the scripts in `R/`:

- **GEO**: GSE17025, GSE63678, GSE106191 (discovery); GSE146889 (held-out replication);
  GSE56087, GSE115810, GSE36389 (screened, not admitted); GSE278879 (single-nucleus)
- **TCGA-UCEC**: GDC Data Portal, via `TCGAbiolinks`
- **Human Protein Atlas**: proteinatlas.org

## Running the scripts

Scripts resolve the project root from the current working directory. Run them from the repository root,
or set `UCEC_PROJ` to point elsewhere:

```bash
git clone https://github.com/zirui1996/UCEC-null-controlled-biomarkers.git
cd UCEC-null-controlled-biomarkers
Rscript R/00_download_GEO.R
```

Each script writes into `results/`, `data_raw/` and `data_processed/`, which are created on first run.

## Licence

MIT — see [LICENSE](LICENSE).
