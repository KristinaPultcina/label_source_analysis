# label_source_analysis
This scripts made linear mixed effect model for HCP parcelations 

Step 1. df_for_label - create a data frame with information about beta power in each label. We create stcs with needed freq and apply the parcelation, next avering sources in labels and time points. Put this info at df

Step 2. Lmem for table. Create lmem model for each label in R

Step 3.Stc&plot with labels. Create the brain vizualization with p-values from lmem (apply space and full FDR)
