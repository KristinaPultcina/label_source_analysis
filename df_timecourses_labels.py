#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: pultsinak
"""


import mne
import numpy as np
import pandas as pd
import os
from scipy import stats
from statsmodels.stats import multitest as mul

os.environ['SUBJECTS_DIR'] = '/net/server/data/Archive/prob_learn/freesurfer'
subjects_dir = '/net/server/data/Archive/prob_learn/freesurfer'
mne.viz.set_3d_options(antialias=(False))

labels= mne.read_labels_from_annot("fsaverage", parc = "aparc_sub")
labels = [lab for lab in labels if 'unknown' not in lab.name]
label_names = [label.name for label in labels]

mfc = [labels[155],labels[157]]
donor= mne.read_source_estimate('/net/server/data/Archive/prob_learn/pultsinak/beta_16_30/stc/stc_epo/P001_run2_norisk_fb_cur_positive/0-lh.stc')
donor= donor.crop(tmin=1.100, tmax=1.900, include_tmax=True)
time_int = donor.times

###### fb negative
subjects_nt=['P001','P004','P019','P021','P022','P034','P035',
               'P040','P044','P047','P048','P053','P059','P060',
               'P063','P064','P065','P067']
             
subjects_aut =['P301','P304','P307','P312','P313','P316',
               'P321','P322','P323','P324','P325','P326','P327','P328',
               'P329','P333','P341','P342']
#### fb positive
subjects_nt=['P001','P004','P019','P021','P022','P034','P035',
               'P040','P044','P047','P048','P053','P059','P060',
               'P063','P064','P065','P067']
             
subjects_aut =['P301','P304','P307','P312','P313','P316',
               'P321','P323','P324','P325','P327','P328',
               'P329','P333','P341','P342']
fsaverage = mne.setup_source_space(subject = "fsaverage", spacing="ico5", add_dist=False)

comp1_per_sub = np.zeros(shape=(len(subjects_nt), len(time_int)))


for ind, subj in enumerate(subjects_nt):

        temp1 = mne.read_source_estimate('/net/server/data/Archive/prob_learn/pultsinak/beta_16_30/stc/stc_avg_epo/{0}_risk_fb_cur_positive'.format(subj), subject="fsaverage").crop(tmin=1.100, tmax=1.900, include_tmax=True)
        
        label_ts_nt = mne.extract_label_time_course(temp1, mfc, src=fsaverage, mode='mean_flip')
        label_ts_nt =label_ts_nt.mean(axis=0)
        
        comp1_per_sub[ind, :] = label_ts_nt                                    


pd.DataFrame(comp1_per_sub).to_csv("/net/server/data/Archive/prob_learn/pultsinak/label_stc/nt_risk_positive_mfc.csv")
            
            
      
comp2_per_sub = np.zeros(shape=(len(subjects_aut), len(time_int)))  
for ind, subj in enumerate(subjects_aut):

    temp2 = mne.read_source_estimate('/net/server/data/Archive/prob_learn/pultsinak/beta_16_30/stc/stc_avg_epo/{0}_risk_fb_cur_positive'.format(subj),subject="fsaverage").crop(tmin=1.100, tmax=1.900, include_tmax=True)
    label_ts_nt = mne.extract_label_time_course(temp2, mfc, src=fsaverage, mode='mean_flip')
    label_ts_nt =label_ts_nt.mean(axis=0)
            
    comp2_per_sub[ind, :] = label_ts_nt    

pd.DataFrame(comp2_per_sub).to_csv("/net/server/data/Archive/prob_learn/pultsinak/label_stc/at_risk_positive_mfc.csv")
       
        
####### create brain picture with needed labels ###       
Brain = mne.viz.get_brain_class()

brain = Brain('fsaverage', 'lh', 'inflated', subjects_dir=subjects_dir,
              cortex='low_contrast', background='white', size=(800, 600))


brain.add_label(mfc[0], borders=False)
brain.add_label(mfc[1], borders=False)



