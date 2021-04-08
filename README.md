# **FS9**

This repo includes all data, code, and relevant results for the study "Route of oxytetracycline administration differentially impacts the growth and gut microbiome of pigs co-infected with *Bordetella bronchiseptica* and *Pasteurella multocida*"

## **Table of contents**
| Chapter | Description |
| -- | -- |
| [data](https://github.com/k39ajdM2/FS9/tree/master/data) | Includes data files needed to carry out R analysis. Also includes `mock_community_16S3.fna` file listing expected mock community sequences for mock sample `Mock5` |
| [scripts](https://github.com/k39ajdM2/FS9/tree/master/scripts) | R scripts for 16S data, qPCR, and other pig sample analyses. AbxConcWeightADGLungLesion.R and AMR_qPCR.R can be run anytime. Only fecal 16S data-related scripts are numbered and must be run in that order. |
| [mothur.md](https://github.com/k39ajdM2/FS9/blob/master/mothur.md) | Commands used in mothur to process FS9 16S data, additional notes |
| [Notes.md](https://github.com/k39ajdM2/FS9/blob/master/Notes.md) | To-do list for FS9, BioProject PRJNA693865 notes |

## Scripts description and the order to run them
| Order | Script file name | Description |
| -- | -- | -- |
| 1a | [mothur.md](https://github.com/k39ajdM2/FS9/tree/master/scripts/mothur.md) | First run `mothur.md` to process 16S sequence data and generate output for R scripts. |
| 1b | [1_OTUtable.R](https://github.com/k39ajdM2/FS9/tree/master/scripts/1_OTUtable.R) | Generate OTU table from mothur output to use for creating phyloseq objects|
| 2a | [2a_phyloseq_NONINFnm_INFnm.R](https://github.com/k39ajdM2/FS9/tree/master/scripts/2a_phyloseq_NONINFnm_INFnm.R) | Generate phyloseq object to use for 3a_alpha_beta_diversity_NONINFnm_INFnm.R. Run [adonis](https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/adonis) function with distance matrices to assess how variation is attributed to different experimental treatments or uncontrolled covariates. |
| 2b | [2b_phyloseq_INFnmfeedinject.R](https://github.com/k39ajdM2/FS9/tree/master/scripts/2b_phyloseq_INFnmfeedinject.R)| Generate phyloseq object to use for 3b_alpha_beta_diversity_INFnmfeedinject.R. Run [adonis](https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/adonis) function with distance matrices to assess how variation is attributed to different experimental treatments or uncontrolled covariates.   |
| 3a | [3a_alpha_beta_diversity_NONINFnm_INFnm.R](https://github.com/k39ajdM2/FS9/tree/master/scripts/3a_alpha_beta_diversity_NONINFnm_INFnm.R) | Run alpha (Shannon, Inverse Simpson) and beta diversity (generating NMDS, pairwise comparisons) analyses, data visualization |
| 3b | [3b_alpha_beta_diversity_INFnmfeedinject.R](https://github.com/k39ajdM2/FS9/tree/master/scripts/3b_alpha_beta_diversity_INFnmfeedinject.R) | Run alpha (Shannon, Inverse Simpson) and beta diversity (generating NMDS, pairwise comparisons) analyses, data visualization |  
| 4a | [4a_DeSeq2_NONINFnm_INFnm.R](https://github.com/k39ajdM2/FS9/tree/master/scripts/4a_DeSeq2_NONINFnm_INFnm.R) | Identify differentially abundant bacterial taxa (order and genus levels) between groups within each day, data visualization |
| 4b | [4b_DeSeq2_INFnmfeedinject.R](https://github.com/k39ajdM2/FS9/tree/master/scripts/4b_DeSeq2_INFnmfeedinject.R) | Identify differentially abundant bacterial taxa (order level) between groups within each day, data visualization |
| 5b | [5b_NMDS_DeSeq2_INFnmfeedinject.R](https://github.com/k39ajdM2/FS9/tree/master/scripts/5b_NMDS_DeSeq2_INFnmfeedinject.R) | Combined steps 2b, 3b, and 4b into one script to generate figure for manuscript |
| Anytime | [AbxConcWeightADGLungLesion.R](https://github.com/k39ajdM2/FS9/tree/master/scripts/AbxConcWeightADGLungLesion.R) | Data visualization of oxytet concentration in various tissue sites, correlation between oxytet concentration and weight for each tissue site, average daily gain, and lung lesion severity |
| Anytime | [AMR_qPCR.R](https://github.com/k39ajdM2/FS9/tree/master/scripts/AMR_qPCR.R) | Data visualization of mean relative abundance (log10) of each AMR gene using line plots; box plots of total abundance of each AMR gene using AULC calculation |

* **Annotation Notes**
  * `(1/2/3/4)a` = NONINFnm vs INFnm
  * `(1/2/3/4/5)b` = INFnm vs INFfeed vs INFinject.
