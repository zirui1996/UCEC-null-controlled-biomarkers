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

### Revision analyses (`R/40a`–`R/40n`)

Reanalyses added during peer review. They run after the main pipeline and read its outputs.

| Script | Analysis |
|---|---|
| `40a` | Four-gene score under alternative null definitions in TCGA-UCEC and GSE146889; patient-clustered bootstrap CI; unpaired, site- and plate-restricted subsets; normalization sensitivity; Monte Carlo CIs |
| `40b` | Selection frequency versus single-gene AUC under Monte Carlo error and co-expression-module resampling |
| `40c` | GSE106191 admission check with unrestricted, within-period and Freedman–Lane permutation |
| `40d` | Superseded per-sample survival analysis (kept for traceability; replaced by `40j`) |
| `40e` | Multivariable gene sets (LASSO, SVM-RFE, bootstrap-stable) evaluated as ridge-logistic classifiers against refitted random sets |
| `40f` | Four-cohort consensus design and GSE17025 downsampling |
| `40g` | Tumour-versus-normal effects of the core genes with and without a fibroblast marker score (genome-wide voom–limma) |
| `40h` | Sample-level inclusion list for the GEO series assessed at sample level |
| `40i` | Normalization diagnostics and scaled Schoenfeld residual figures |
| `40j` | Patient as the unit of analysis: survival with one record per patient; discrimination with one tumour per patient |
| `40k` | Nearest-neighbour effect-size-matched null with distinct genes per panel |
| `40l` | GSE106191 under its prespecified, sensitivity and adopted analysis sets |
| `40m` | Monte Carlo CIs for the empirical P values of the multivariable, admission and four-cohort analyses |
| `40n` | Kaplan–Meier and Cox figure with one record per patient |

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
