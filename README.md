# MICCA: Skin microbiome composite features in Atopic Dermatitis via integration analysis across cohorts

​                                                                                              [Yihua Wang](https://yihuawww.github.io/)

## Prerequisites

**R** language



## Description

This project mainly includes following five functions for <u>Skin Microbiome Composite Features Analysis</u>:

> (1) lm_meta.R
>
> (2) DA.R
>
> (3) lefse.R
>
> (4) Cross_cohort_Feature_Selection_and_Classification.R
>
> (5) DA_guided_Feature_Selection_Classification.R



Clone this repository and install the corresponding dependent libraries as bellow:

### 1. lm_meta.R

```
Description:
This script performs:
 1. Batch-effect correction using MMUPHin
 2. Meta-analytical differential abundance analysis
 3. Differential genus identification for: ADLS vs HC, ADNLS vs HC, and ADLS vs ADNLS
 4. Visualization of significant genera
```

```
Input:
 data/group.csv # meta data of samples
 data/tax.txt   # abundance of tax
```

```R
Dependent Packages:
install.package ("magrittr")
install.package ("ggplot2")
install.package ("vegan")
install.package ("ggpubr")
install.package ("ggsci")
install.package ("patchwork")
install.package ("MMUPHin")
```

### 2. DA.R

```
Description:
This script identifies differentially abundant genera between two groups using permutation tests and Wilcoxon rank-sum tests.

Note that Genera are filtered according to: 
1. Minimum mean relative abundance
2. Minimum prevalence across samples
```

```R
Input:
 data/tax_L6.txt    # abundance of tax 
 data/meta_all.txt  # meta data of samples
```

```R
Dependent Packages:
install.package ("dplyr")
```

### 3. lefse.R

```
Description:
This script:
1. Performs LEfSe analysis using the microeco package
2. Identifies differentially abundant genera

Note that: It's comparisons among ADLS vs HC, ADNLS vs HC, and ADLS vs ADNLS.
```

```R
Input:
data/
    ├── tax_L6.txt     # abundance of data 
    ├── meta_all.txt   # meta data of samples
    └── lefse/
       ├── otutable.txt   # OTUs abundance
       └── tax_table.txt  # taxonomy information 
```

```R
Dependent Packages:
install.package ("tidyverse")
install.package ("microeco")
install.package ("magrittr")
```

### 4. Cross_cohort_Feature_Selection_and_Classification.R

```
Description:
This script performs a multi-cohort microbiome analysis consisting of three major steps:
1. Batch-effect correction using MMUPHin
2. Cross-cohort meta-analytical differential abundance analysis
3. Random Forest classification with feature selection using RFE and independent-cohort validation

Note that it's main comparison between ADLS vs HC.
Training cohorts: S2, S17, S19, S47, S71, S75
Independent validation cohort: S74
```

```R
Input:
data/
    └── all/
           ├── group.csv   # meta data of samples
           └── tax_L6.txt  #abundance of tax
```

```R
Dependent Packages:
install.package ("MMUPHin")
install.package ("dplyr")
install.package ("ggplot2")
install.package ("vegan")
install.package ("ggpubr")
install.package ("ggsci")
install.package ("patchwork")
install.package ("magrittr")
install.package ("randomForest")
install.package ("caret")
install.package ("pROC")
```

### 5. DA_guided_Feature_Selection_Classification.R

```
Description:
This script performs: 
1. Identification of robust differential genera across cohorts
2. Feature selection based on the intersection of Permutation test, Wilcoxon rank-sum test, and LEfSe
3. Study-stratified 5-fold cross-validation
4. Recursive Feature Elimination (RFE)
5. Random Forest classification
6. Independent validation in a held-out cohort
7. Evaluation of classification performance using Accuracy, Sensitivity, Specificity, Precision, F1-score , Balanced Accuracy, AUC.
8. Saving fold assignments and validation results

Note that it is comparison of ADLS vs HC.

Training cohorts: S2, S17, S19, S47, S71, S75
Independent validation cohort: S74
```

```R
Dependent Packages:
install.package ("MMUPHin")
install.package ("dplyr")
install.package ("ggplot2")
install.package ("vegan")
install.package ("ggpubr")
install.package ("ggsci")
install.package ("patchwork")
install.package ("magrittr")
install.package ("randomForest")
install.package ("caret")
install.package ("pROC")
```
