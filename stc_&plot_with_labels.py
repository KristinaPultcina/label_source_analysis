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




os.environ['SUBJECTS_DIR'] = '/Volumes/My Passport for Mac/freesurfer'
subjects_dir = '/Volumes/My Passport for Mac/freesurfer'


labels= mne.read_labels_from_annot("fsaverage", parc = "aparc_sub")


labels = [lab for lab in labels if 'unknown' not in lab.name]
label_names = [label.name for label in labels]


df = pd.read_csv('/Users/kristina/Documents/stc/lmem_label/lmem_label_1100_1500.csv', sep = ";") 


 
    
src = mne.setup_source_space(subject = "fsaverage", spacing='ico5', add_dist=False)
stc_test = mne.read_source_estimate('/Users/kristina/Documents/stc/lmem_label/P001_run2_norisk_fb_cur_negative_fsaverage/0', "fsaverage").crop(tmin= -0.800, tmax= 2.100, include_tmax=True)
stc_test.resample(10)
    
    

################### CREATE STC ###################
data =df['group:trial_type']
#df1 = data.drop([448])
data = np.array(data)

space_fdr_data = mul.fdrcorrection(data,alpha=0.05)
space_fdr_data=space_fdr_data[1]


####### rescale for vizualization 10p ########
data_pval = []
for p in data:
    pval = 1 - 10*p
    data_pval.append(pval)
data_pval = np.array(data_pval)

fdr_data_pval = []
for f in space_fdr_data:
    fdr_pval = 1 - 10*f
    fdr_data_pval.append(fdr_pval)
fdr_data_pval = np.array(fdr_data_pval)

nofdr_stc = mne.labels_to_stc(label_list, data_pval, tmin = -0.800, tstep = 0.5, src=src)
nofdr_stc.save ('/Users/kristina/Documents/stc/lmem_label/stc_for_article/averaged_stc/hp_lp_avg_10p',overwrite=True)


fdr_stc = mne.labels_to_stc(label_list,fdr_data_pval, tmin = -0.800, tstep = 0.5, src =src)
fdr_stc.save ('/Users/kristina/Documents/stc/lmem_label/stc_for_article/model_subj_epoches/space_fdr_norisk_risk_10p', overwrite=True)


####### plots ############
intervals = [[1.500, 1.900]]

intervals = [[-0.900, -0.300]]
mean_time = [-0.700]
intervals = [[0.700, 0.900]] # time intervals for averaging
intervals = [[1.500, 1.900]] # time intervals for averaging
mean_time = [1.700]
name_int =['before_fb']
mean_time = [0.800]



p_value= [0.45,0.55,1.0]
var_of_ploting = ['p_value']
scale = [p_value]



####### create plot ##########

for idx, inter in enumerate(intervals):
    for ind, v in enumerate(var_of_ploting):
            stc = mne.read_source_estimate('/Users/kristina/Documents/stc/lmem_label/stc_for_article/averaged_stc/hreml', 'fsaverage')
        
            #for t in time_points:
            brain = mne.viz.plot_source_estimates(stc, hemi='split', time_viewer=False, background='white', 
                                                  foreground = 'black', cortex='bone', size = (1200, 600),
                                                        views = ['lat', 'med'], clim = dict(kind = 'value', 
                                                                                          lims= scale[ind]), 
                                                  initial_time = mean_time[idx], time_label=f'{inter[0]} .... {inter[1]} s',
                                                       spacing ='ico5', title = f'nofdr LP vs HP, beta power 16 - 30 Hz, {inter[0]} .... {inter[1]} s')
          
                                                       
            brain.add_text(0.0, 0.9, f'no_fdr_LP_losses_vs_wins_beta 16-30Hz, **{scale[ind]}**',
                               font_size=12, color='black')
            
         
            #brain.add_text(0.0, 0.8, f'{inter[0]} .... inter[1]s',
            #                   font_size=10, color='green')                   
            brain.add_text(0.0, 0.8, f'{v}', font_size=10, color='blue')
            #brain.add_label(label_list[4])
            #brain.add_annotation('HCPMMP1')
            brain.save_image('/Users/kristina/Documents/stc/lmem_label/lp_hp/lp_hp_plus.jpeg')
            brain.close()


