# Input data

These are gzip-compressed snapshots of the original project inputs. The
uncompressed input files and original mutation ZIP are also tracked at the
repository root; `project.Rmd` reads the compressed copies. The clinical file
has 512 patients, the mutation file has 356 patients, and the
RNA-seq count matrix has 533 patients with a primary-tumor sample. Matching
the three datasets and retaining recorded overall-survival fields yields 354
patients; two have zero months of follow-up.

The clinical and mutation files match the format of the
[TCGA KIRC PanCancer Atlas study on cBioPortal](https://www.cbioportal.org/study/filesAndLinks?id=kirc_tcga_pan_can_atlas_2018).
The exact source accession for the supplied RNA-seq count matrix was not
recorded. The [NCI PanCancer Atlas page](https://gdc.cancer.gov/about-data/publications/pancanatlas)
provides background on the broader resource.

SHA-256 checksums of the committed archives:

```text
fa428fd74c999610f8c1a4d0c89d621776bf7b99fe153eb5fd256782d78dcd9e  RNAseq_KIRC.csv.gz
1196823532626f1e277f549a1c6d94da9c17c61a9cbea6a87939cd238393d9de  data_clinical_patient.txt.gz
ef7cbc01a437af11cefdad1ed6b1a7773daec4bd32e5cfc0badc101f1d9d9615  data_mutations.txt.gz
```
