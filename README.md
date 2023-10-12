## scRNA-seq

### Unclear Data

For unclear data use the FilterQC.R to filter samples. You will need also to load helperfunctions.R.

### Clear PBMC and B cell Data

- PBMC data are composed by scRNA-seq + CITE-seq data.
- B cells data are just scRNA-seq.

Analyses comprehends: (1) VDJ (TCR,BCR-seq data) Doublet's removal; (2) Normalization; (3) Automatic Annotation;
(4) Integration between scRNA-seq and CITE-seq; (5) Reduce Dimension and Clustering.


