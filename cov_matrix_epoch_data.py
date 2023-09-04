import mne
import os
import os.path as op
import numpy as np
import pandas as pd
from mne.preprocessing import ICA
mne.viz.set_3d_options(antialias=False)

os.environ['SUBJECTS_DIR'] = '/media/kristina/storage/probability/freesurfer'
subjects_dir = '/media/kristina/storage/probability/freesurfer'
subjects= ['P001','P301','P304',"P307",'P312', 'P313', 'P314',
                                       'P316', 'P321','P322','P323','P324',"P325",
                                       'P326','P327','P328', 'P329','P333','P335',
                                       'P334','P341','P342','P004','P019','P021','P022','P034','P035','P039',
                                      'P040','P044','P047','P048','P053','P055','P058', 'P059',"P060", 'P061',
                                      'P065','P063','P064','P067']

rounds = [1, 2, 3, 4, 5, 6]



trial_type = ['norisk', 'risk', 'postrisk','risk']


feedback = ['positive', 'negative']

data_path = '/media/kristina/storage/probability/ICA_cleaned'


for subj in subjects:
    for r in rounds:
        events_list = []
        raw_fname = op.join(data_path, '{0}/run{1}_{0}_raw_ica.fif'.format(subj, r))

        raw_data = mne.io.Raw(raw_fname, preload=True)
        picks = mne.pick_types(raw_data.info, meg =True,exclude='bads')

        print(r)
        for t in trial_type:
            for fb in feedback:

                try:
                    data = np.loadtxt("/media/kristina/storage/probability/fix_cross_mio_corr_lp_hp/{0}_run{1}_{2}_fb_cur_{3}_fix_cross.txt".format(subj, r,t,fb), dtype='int') 
                    
                    if data.size > 0:
                        # Append data to the events_list
                        events_list.append(data)
                        
                except (OSError):
                    print('This file does not exist')
        if len(events_list) > 0:
            events = np.vstack(events_list)
            print(events.shape)
            if events.shape[0] >10:
                epochs_bl = mne.Epochs(raw_data, events, event_id = None, tmin = -1.0, tmax = 1.0, baseline = None, picks = picks, preload = True,event_repeated='merge')

                cov = mne.compute_covariance(epochs=epochs_bl, method=("oas"), tmin=-0.35, tmax = -0.05,n_jobs=-1)
        
                cov = mne.cov.regularize(cov, raw_data.info, mag=0.1, grad=0.1, eeg=0.1, rank='info')
            
                cov.plot(raw_data.info)
                cov.save('/media/kristina/storage/probability/cov_bl_data/{0}_run{1}_er-cov.fif'.format(subj,r),overwrite=True)
              
        
