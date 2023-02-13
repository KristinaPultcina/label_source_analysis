# lmem_source_analysis (based on parcelation)
This scripts made linear mixed effect model for parcelation by S.Khan (aparc_sub)

Step 1. create_stc_for_single_epo.py - Create the stc files for each epoch with extracting needed freq bands. For source localization apply sLoreta

Step 2. df_for_label - create a data frame with information about beta power in each label in each participant

Step 3. Lmem for table. Create lmem model for each label in R

Step 4.Stc&plot with labels. Create the brain vizualization with p-values from lmem (apply space and full FDR)
 
Step 5. brain_with_marginal_means_and_pvalue.py create brain vizualization with pvalue and marginal means on it
