import mne
import os.path as op
import os
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import (make_axes_locatable, ImageGrid,
                                     inset_locator)
from matplotlib  import pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import copy
import statsmodels.stats.multitest as mul
from mne.stats import  fdr_correction
mne.viz.set_3d_options(antialias=False)



os.environ['SUBJECTS_DIR'] = '/Volumes/My Passport for Mac/freesurfer'
subjects_dir = '/Volumes/My Passport for Mac/freesurfer'


labels= mne.read_labels_from_annot("fsaverage", parc = "aparc_sub")


labels = [lab for lab in labels if 'unknown' not in lab.name]
label_names = [label.name for label in labels]


df = pd.read_csv('/Users/kristina/Documents/stc/lmem_label/dfs/tukey_label_1100_1500_betweengroup.csv')


src = mne.setup_source_space(subject = "fsaverage", spacing="ico5", add_dist=False)

###### create plot  with p_value with marginal means ######
data =df['positive-risk-autists - normal']
#df1 = data.drop([448])
data = np.array(data)
data[data>0.05]=1
space_fdr_data = mul.fdrcorrection(data,alpha=0.05)
space_fdr_data=space_fdr_data[1]
space_fdr_data[space_fdr_data>0.05]=1
mean_data =df['lp_pos_aut_norm']
#df1 = data.drop([448])
mean_data= np.array(mean_data)
label_number= np.linspace(0, 447, 448)
arr= np.column_stack((data, mean_data, label_number))


arr[arr[:, 0] == 1, 1] = 0
space_arr= np.column_stack((space_fdr_data, mean_data, label_number))
space_arr[space_arr[:, 0] == 1, 1] = 0

nofdr_stc = mne.labels_to_stc(labels,arr[:,1], subject="fsaverage")

nofdr_stc.save ('/Users/kristina/Documents/stc/lmem_label/stc/beta_diff/1500_1900/group_trial_type_feedback_1500_1900_fdr',overwrite=True)


fdr_stc = mne.labels_to_stc(labels,space_arr[:,1], tmin = -0.800, tstep = 0.1)
fdr_stc.save ('/Users/kristina/Documents/stc/lmem_label/stc/group_trial_type_feedback_fdr', overwrite=True)

####### create brain picture, don't forget change pos_lims
stc=nofdr_stc
            #for t in time_points:
brain = mne.viz.plot_source_estimates(stc, hemi='split', time_viewer=False, background='white', 
                                                  foreground = 'black', cortex='bone', size = (1200, 800),
                                                        views = ['lat', 'med'], clim = dict(kind = 'value', 
                                                                                          pos_lims= [0.25,0.30, 1]), ###### choose scale!!!!
                                                  initial_time = -0.700, spacing ='ico5')
                                



