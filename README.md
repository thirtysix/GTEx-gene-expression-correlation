

## Simple script for calculating gene expression Spearman correlations between a target gene and all other genes, using the bulk RNA-Seq GTEx dataset.


# 1. Background
The current version of the [GTEx](https://gtexportal.org/home/) dataset is available as open access data and can be downloaded from here:

[V8 GTEx open access data](https://gtexportal.org/home/datasets)

The script accesses the gzipped version of the data and only loads one sample type at a time, thus reducing active memory (RAM) requirements.  This allows the analysis to be run on systems with smaller amounts of installed memory.


# 1.1 GTEx tissues which are analyzed for correlation between your target gene(s) and all other genes.

| Top-level       | Sub-tissue                                |
| --------------- | ----------------------------------------- |
| Adipose Tissue  | Adipose - Subcutaneous                    |
|                 | Adipose - Visceral (Omentum)              |
| Adrenal Gland   | Adrenal Gland                             |
| Bladder         | Bladder                                   |
| Blood           | Whole Blood                               |
|                 | Cells - EBV-transformed lymphocytes       |
| Blood Vessel    | Artery - Aorta                            |
|                 | Artery - Coronary                         |
|                 | Artery - Tibial                           |
| Brain           | Brain - Amygdala                          |
|                 | Brain - Anterior cingulate cortex (BA24)  |
|                 | Brain - Caudate (basal ganglia)           |
|                 | Brain - Cerebellar Hemisphere             |
|                 | Brain - Cerebellum                        |
|                 | Brain - Cortex                            |
|                 | Brain - Frontal Cortex (BA9)              |
|                 | Brain - Hippocampus                       |
|                 | Brain - Hypothalamus                      |
|                 | Brain - Nucleus accumbens (basal ganglia) |
|                 | Brain - Putamen (basal ganglia)           |
|                 | Brain - Spinal cord (cervical c-1)        |
|                 | Brain - Substantia nigra                  |
| Breast          | Breast - Mammary Tissue                   |
| Cervix Uteri    | Cervix - Ectocervix                       |
|                 | Cervix - Endocervix                       |
| Colon           | Colon - Sigmoid                           |
|                 | Colon - Transverse                        |
| Esophagus       | Esophagus - Gastroesophageal Junction     |
|                 | Esophagus - Mucosa                        |
|                 | Esophagus - Muscularis                    |
| Fallopian Tube  | Fallopian Tube                            |
| Heart           | Heart - Atrial Appendage                  |
|                 | Heart - Left Ventricle                    |
| Kidney          | Kidney - Cortex                           |
|                 | Kidney - Medulla                          |
| Liver           | Liver                                     |
| Lung            | Lung                                      |
| Muscle          | Muscle - Skeletal                         |
| Nerve           | Nerve - Tibial                            |
| Ovary           | Ovary                                     |
| Pancreas        | Pancreas                                  |
| Pituitary       | Pituitary                                 |
| Prostate        | Prostate                                  |
| Salivary Gland  | Minor Salivary Gland                      |
| Skin            | Skin - Not Sun Exposed (Suprapubic)       |
|                 | Skin - Sun Exposed (Lower leg)            |
|                 | Cells - Cultured fibroblasts              |
| Small Intestine | Small Intestine - Terminal Ileum          |
| Spleen          | Spleen                                    |
| Stomach         | Stomach                                   |
| Testis          | Testis                                    |
| Thyroid         | Thyroid                                   |
| Uterus          | Uterus                                    |
| Vagina          | Vagina                                    |


# 1.2 Requirements

These two tools will need to be installed to analyze the data:

 - [Pandas](https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html)
 - [SciPy](https://docs.scipy.org/doc/scipy/getting_started.html)


# 2. Input

## 2.1 GTEx data files
Two data files will need to be downloaded from GTEx in order to perform the analysis:

 - [GTEx V8 gene TPMs - 1.6GB](https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz)
 - [GTEx V8 sample attributes - 11MB](https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt)

 These files should be downloaded on your own and then manually placed in the '01.raw_data/V8/' directory.

## 2.2 User-defined list of target genes
You will need to manually edit the '02.data/target_genes.txt' file to contain a list of your target genes; one gene per line.


# 3. Analysis

There are two levels of tissue designation in the GTEx RNA-Seq dataset.  The script will compare the gene expression of your target gene(s) against all other genes in each of the top-level tissues and sub-tissues.



# 4. Output
## 4.1 Results
The output of the script are individual tab-separated (TSV) files, for each combination of target gene and GTEx tissue (top-level tissue or sub-tissue) (e.g., "Breast - Mammary Tissue.BRCA2.correlation.tsv").

Each output file contains the following columns:

| Column name            | Description                                                           |
| ---------------------- | --------------------------------------------------------------------- |
| target_tissue          | tissue in which the correlation occurs.                               |
| tissue_n               | number of GTEx samples comprising this tissue data.                   |
| target_gene_symbol     | gene symbol of your target gene.                                      |
| target_gene_id         | ensembl ID of your target gene.                                       |
| correlated_gene_symbol | gene symbol of the correlated gene.                                   |
| correlated_gene_id     | ensembl ID of the correlated gene.                                    |
| correlation            | Spearman correlation value of the target gene - correlated gene pair. |
| pval                   | raw p-value of the Spearman correlation                               |
| target_gene_expr       | average expression (TPM) of the target gene.                          |
| correlated_gene_expr   | average expression (TPM) of the correlated gene.                      |
| pval corrected         | FDR corrected p-value.                                                |



