#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import mne
import os
import os.path as op
import numpy as np
import pandas as pd

os.environ['SUBJECTS_DIR'] = '/net/server/data/Archive/prob_learn/freesurfer'
subjects_dir = '/net/server/data/Archive/prob_learn/freesurfer'
subjects = [
        	'P316','P034',
        	'P320','P035',
        	'P321','P039',
        	'P322','P040',
        	'P323','P044',
        	'P324','P047',
        	'P325','P048',
        	'P326','P053',
        	'P327','P055',
        	'P328','P058',
        	'P329','P059',
        	'P333','P060',
        	'P334','P061',
        	'P335','P063',
        	'P338','P064',
        	'P341','P065',
        	'P342','P067']
L_freq = 16
H_freq = 31
f_step = 2

#time_bandwidth = 4 #(by default = 4)
# if delta (1 - 4 Hz) 
#n_cycles = np.array([1, 1, 1, 2]) # уточнить


period_start = -1.750
period_end = 2.750

baseline = (-0.35, -0.05)

freq_range = 'beta_16_30'
rounds = [1, 2, 3, 4, 5, 6]

trial_type = ['norisk', 'prerisk', 'risk', 'postrisk']


feedback = ['positive', 'negative']

data_path = '/net/server/data/Archive/prob_learn/vtretyakova/ICA_cleaned'

def make_stc_epochs_from_freq_epochs_var2(subj, r, cond, fb, data_path, L_freq, H_freq, f_step, period_start, period_end, baseline, bem, src):

    bands = dict(beta=[L_freq, H_freq])
    #freqs = np.arange(L_freq, H_freq, f_step)
    
    #read events
	#events for baseline
	# download marks of positive feedback
	
    events_pos = np.loadtxt("/net/server/data/Archive/prob_learn/data_processing/fix_cross_mio_corr/{0}_run{1}_norisk_fb_cur_positive_fix_cross.txt".format(subj, r), dtype='int') 
    

        # если только одна метка, т.е. одна эпоха, то выдается ошибка, поэтому приводим shape к виду (N,3)
    if events_pos.shape == (3,):
        events_pos = events_pos.reshape(1,3)
        
    # download marks of negative feedback      
    
    events_neg = np.loadtxt("/net/server/data/Archive/prob_learn/data_processing/fix_cross_mio_corr/{0}_run{1}_norisk_fb_cur_negative_fix_cross.txt".format(subj, r), dtype='int')
    
    
    # если только одна метка, т.е. одна эпоха, то выдается ошибка, поэтому приводим shape к виду (N,3)
    if events_neg.shape == (3,):
        events_neg = events_neg.reshape(1,3) 
    
    #объединяем негативные и позитивные фидбеки для получения общего бейзлайна по ним, и сортируем массив, чтобы времена меток шли в порядке возрастания    
    events = np.vstack([events_pos, events_neg])
    events = np.sort(events, axis = 0) 
    
    #events, which we need
    events_response = np.loadtxt('/net/server/data/Archive/prob_learn/data_processing/events_trained_by_cond_WITH_mio_corrected/{0}_run{1}_{2}_fb_cur_{3}.txt'.format(subj, r, cond, fb), dtype='int')
    # если только одна метка, т.е. одна эпоха, то выдается ошибка, поэтому приводи shape к виду (N,3)
    if events_response.shape == (3,):
        events_response = events_response.reshape(1,3)
    
	           
    raw_fname = op.join(data_path, '{0}/run{1}_{0}_raw_ica.fif'.format(subj, r))

    raw_data = mne.io.Raw(raw_fname, preload=True)
        
    
    picks = mne.pick_types(raw_data.info, meg = True, eog = True)
		    
	# Forward Model
    trans = '/net/server/mnt/Archive/prob_learn/freesurfer/{0}/mri/T1-neuromag/sets/{0}-COR.fif'.format(subj)
        
	   	    
    #epochs for baseline
    # baseline = None, чтобы не вычитался дефолтный бейзлайн
    epochs_bl = mne.Epochs(raw_data, events, event_id = None, tmin = -1.0, tmax = 1.0, baseline = None, picks = picks, preload = True)
    cov = mne.compute_covariance(epochs=epochs_bl, method='auto', tmin=-0.35, tmax = -0.05)
     
    #epochs_bl.resample(100)
    
    ####### ДЛЯ ДАННЫХ ##############
    # baseline = None, чтобы не вычитался дефолтный бейзлайн
    epochs = mne.Epochs(raw_data, events_response, event_id = None, tmin = period_start, 
		                tmax = period_end, baseline = None, picks = picks, preload = True)


                
    fwd = mne.make_forward_solution(info=epochs.info, trans=trans, src=src, bem=bem)	                
    inv = mne.minimum_norm.make_inverse_operator(raw_data.info, fwd, cov, loose=0.2) 	                
		       
    #epochs.resample(100) 
    
    stc_baseline = mne.minimum_norm.source_band_induced_power(epochs_bl.pick('meg'), inv, bands,method='sLORETA', use_fft=False, df = f_step, n_cycles = 8)["beta"].crop(tmin=baseline[0], tmax=baseline[1], include_tmax=True)
    
    
    #усредняем по времени
    b_line  = stc_baseline.data.mean(axis=-1)
    b_line = 10*np.log10(b_line)
    
    stc_epo_list = []
    for ix in range(len(epochs)):
        stc = mne.minimum_norm.source_band_induced_power(epochs[ix].pick('meg'), inv, bands,method='sLORETA', use_fft=False, df = f_step, n_cycles = 8)["beta"]
        data = 10*np.log10(stc.data)
        #stc.data = 10*np.log10(stc.data/b_line[:, np.newaxis])
        stc.data = data - b_line[:, np.newaxis]
        stc_epo_list.append(stc)
    return (stc_epo_list)




for subj in subjects:
    bem = mne.read_bem_solution('/net/server/data/Archive/prob_learn/data_processing/bem/{0}_bem.h5'.format(subj), verbose=None)
    src = mne.setup_source_space(subject =subj, spacing='ico5', add_dist=False ) # by default - spacing='oct6' (4098 sources per hemisphere)
    for r in rounds:
        for cond in trial_type:
            for fb in feedback:

                try:
                
                    
                    stc_epo_list = make_stc_epochs_from_freq_epochs_var2(subj, r, cond, fb, data_path, L_freq, H_freq, f_step, period_start, period_end, baseline, bem, src)
                    print('Количество эпох %s' % len(stc_epo_list))
                    
                    os.makedirs('/net/server/data/Archive/prob_learn/pultsinak/beta_16_30/stc/stc_epo/{0}_run{1}_{2}_fb_cur_{3}'.format( subj, r, cond, fb))
                    
                    for s in range(len(stc_epo_list)):
                        
                        stc_epo_list[s].save('/net/server/data/Archive/prob_learn/pultsinak/beta_16_30/stc/stc_epo/{0}_run{1}_{2}_fb_cur_{3}/{4}'.format(subj, r, cond, fb, s))
                    
                                
         
                except (OSError):
                    print('This file not exist')

############### MORPHING ##################
for subj in subjects:
    for r in rounds:
        for cond in trial_type:
            for fb in feedback:
                try:
                    epochs_num = os.listdir('/net/server/data/Archive/prob_learn/pultsinak/beta_16_30/stc/stc_epo/{0}_run{1}_{2}_fb_cur_{3}'.format( subj, r, cond, fb))
                    #print(subj)
                    #print(r)
                    #print(fb)
                    print (int(len(epochs_num)/2))
                    os.makedirs('/net/server/data/Archive/prob_learn/pultsinak/beta_16_30/stc/stc_epo_fsaverage/{0}_run{1}_{2}_fb_cur_{3}_fsaverage'.format(subj, r, cond, fb))
                
                    for ep in range(int(len(epochs_num)/2)):
                    
                        stc= mne.read_source_estimate("/net/server/data/Archive/prob_learn/pultsinak/beta_16_30/stc/stc_epo/{0}_run{1}_{2}_fb_cur_{3}/{4}".format(subj, r, cond, fb, ep))
                        morph = mne.compute_source_morph(stc, subject_from=subj, subject_to='fsaverage')
                        stc_fsaverage = morph.apply(stc)
                        stc_fsaverage.save('/net/server/data/Archive/prob_learn/pultsinak/beta_16_30/stc/stc_epo_fsaverage/{0}_run{1}_{2}_fb_cur_{3}_fsaverage/{4}'.format(subj, r, cond, fb, ep))
                        
                except (OSError):
                    print('This file not exist')


    
