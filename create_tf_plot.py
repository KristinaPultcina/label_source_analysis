

import mne
import os
import os.path as op
import numpy as np
import pandas as pd

os.environ['SUBJECTS_DIR'] = '/net/server/data/Archive/prob_learn/freesurfer'
subjects_dir = '/net/server/data/Archive/prob_learn/freesurfer'
subjects= ['P301','P304',"P307",'P312', 'P313', 'P314',
                                       'P316', 'P321','P322','P323','P324',"P325",
                                       'P326','P327','P328', 'P329','P333','P335','P334',
                                       'P334','P341','P342', 'P001','P004','P019','P021','P022','P034','P035','P039',
                                      'P040','P044','P047','P048','P053','P055','P058', 'P059',"P060", 'P061',
                                      'P065','P063','P064','P067']
                                      


L_freq = 2
H_freq = 41
f_step = 2
freqs = np.arange(L_freq, H_freq, f_step)


n_cycles = freqs//2


subjects= ['P301']
period_start = -1.750
period_end = 2.750

baseline = (-0.35, -0.05)

freq_range = 'beta_16_30'
rounds = [1, 2, 3, 4, 5, 6]

trial_type = ['norisk', 'risk']


feedback = ['positive', 'negative']

data_path = '/net/server/data/Archive/prob_learn/vtretyakova/ICA_cleaned'

labels =  mne.read_labels_from_annot("fsaverage", "aparc_sub", hemi = "both")
labels = [lab for lab in labels if 'unknown' not in lab.name]
label_names = [label.name for label in labels] 
lattoccip_rh =labels[111]+labels[113]+labels[115]+labels[117]+labels[121]+labels[108]+labels[123]

#def make_tf_array(subj, r, cond, fb, data_path, L_freq, H_freq, f_step, period_start, period_end, baseline, bem, src):

                    
for subj in subjects:
    bem = mne.read_bem_solution('/net/server/data/Archive/prob_learn/data_processing/bem/{0}_bem.h5'.format(subj), verbose=None)
    src = mne.setup_source_space(subject =subj, spacing='ico5', add_dist=False ) # by default - spacing='oct6' (4098 sources per hemisphere)
    for r in rounds:
        for cond in trial_type:
            for fb in feedback:

                try:
                    bands = dict(beta=[L_freq, H_freq])
                   
                 	
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
                    trans = '/net/server/data/Archive/prob_learn/freesurfer/{0}/mri/T1-neuromag/sets/{0}-COR.fif'.format(subj)
                         
                 	   	    
                     #epochs for baseline
                     # baseline = None, чтобы не вычитался дефолтный бейзлайн
                    epochs_bl = mne.Epochs(raw_data, events, event_id = None, tmin = -1.0, tmax = 1.0, baseline = None, picks = picks, preload = True)
                    cov = mne.compute_covariance(epochs=epochs_bl, method='auto', tmin=-0.35, tmax = -0.05)
                    epochs_bl.resample(300)
                    
                     ####### ДЛЯ ДАННЫХ ##############
                     # baseline = None, чтобы не вычитался дефолтный бейзлайн
                    epochs = mne.Epochs(raw_data, events_response, event_id = None, tmin = period_start, 
                 		                tmax = period_end, baseline = None, picks = picks, preload = True)
                    epochs.resample(300)
                                 
                    fwd = mne.make_forward_solution(info=epochs.info, trans=trans, src=src, bem=bem)	                
                    inv = mne.minimum_norm.make_inverse_operator(raw_data.info, fwd, cov, loose=0.2) 	                

                    power_bl, itc_bl = mne.minimum_norm.source_induced_power(epochs_bl, inv, method='sLORETA',freqs=freqs,n_cycles=freqs / 3.0, label=lattoccip_rh,lambda2=1. / 9.,use_fft=False)
                    power_bl = np.mean(power_bl, axis=0)
                    itc_bl = np.mean(itc_bl, axis=0)
                    
                    bl_array=power_bl[:,195:285] ######## take baseline period
                    bl_power=10*np.log10(bl_array.mean(axis=1))
                    
                    itc_array= itc_bl[:,195:285]
                    bl_itc=10*np.log10(itc_array.mean(axis=1))
                    
                    stc_power,stc_itc= mne.minimum_norm.source_induced_power(epochs, inv, method='sLORETA',freqs=freqs,n_cycles=freqs / 3.0, label=lattoccip_rh,lambda2=1. / 9.,use_fft=False)
                    
                    power = 10*np.log10(np.mean(stc_power,axis=0))
                    itc= 10*np.log10(np.mean(stc_itc,axis=0))
                    
                    
                    itc_with_bl = itc - bl_itc.reshape(20,1)
                    power_with_bl = power - bl_power.reshape(20,1)
                    
                    np.save('/net/server/data/Archive/prob_learn/pultsinak/beta_16_30/stc/tfr/lateraloccipital_rh/tf_itc/{0}_run{1}_{2}_fb_cur_{3}'.format(subj, r, cond, fb),itc_with_bl)

                    np.save('/net/server/data/Archive/prob_learn/pultsinak/beta_16_30/stc/tfr/lateraloccipital_rh/tf_power/{0}_run{1}_{2}_fb_cur_{3}'.format(subj, r, cond, fb),power_with_bl)
                    
                                
         
                except (OSError):
                    print('This file not exist')
                                       
                    
     
temp = np.load('/net/server/data/Archive/prob_learn/pultsinak/beta_16_30/stc/tfr/lateraloccipital_rh/tf/P001_run2_norisk_fb_cur_negative.npy')
temp=np.mean(temp,axis=0)
n = temp.shape[1] # количество временных точек (берем у донора, если донор из тех же данных.
sn = temp.shape[0]                    
                    
                    
for subj in subjects:
    for t in trial_type:
        for fb in feedback:
            data_fb = np.empty((0, sn, n))
            for r in rounds:
                
                
                try:
                    power = np.load('/net/server/data/Archive/prob_learn/pultsinak/beta_16_30/stc/tfr/lateraloccipital_rh/tf/{0}_run{1}_{2}_fb_cur_{3}.npy'.format(subj, r, t, fb))                 
                    power = power.mean(axis=0)
                    power = power.reshape(1, sn, n) # добавляем ось блока (run)
                    
                except (OSError):
                    power = np.empty((0, sn, n))
                    print('This file not exist')
                    
                data_fb = np.vstack([data_fb, power])  # собираем все блоки в один массив 
                
            if data_fb.size != 0:
                temp.data = data_fb.mean(axis = 0)    # усредняем между блоками (run)
                np.save('/net/server/data/Archive/prob_learn/pultsinak/beta_16_30/stc/tfr/lateraloccipital_rh/tf_ave_in_subj/{0}_{1}_fb_cur_{2}'.format(subj, t, fb),temp)
            else:
                print('Subject has no feedbacks in this condition')
                pass

temp = np.load('/net/server/data/Archive/prob_learn/pultsinak/beta_16_30/stc/tfr/lateraloccipital_rh/tf_power/P301_run3_norisk_fb_cur_negative.npy')

n = temp.shape[1] # количество временных точек (берем у донора, если донор из тех же данных.
sn = temp.shape[0]                    
                                          
for subj in subjects:
    for t in trial_type:
        ################################ Positive feedback ################################
        positive_fb = np.empty((0,sn,n))
        for r in rounds:
            try:
                               
                epochs_positive = np.load('/net/server/data/Archive/prob_learn/pultsinak/beta_16_30/stc/tfr/lateraloccipital_rh/tf_itc/{0}_run{1}_{2}_fb_cur_positive.npy'.format(subj, r, t))             
                
                positive_fb = np.vstack([positive_fb, epochs_positive.reshape(1,20,1350)])
               
            except (OSError):
                
                print('This file not exist')

        ###### Шаг 1. Усреднили все положительные фидбеки внутри испытуемого (между блоками 1 -6) #################
        if positive_fb.size != 0:
            positive_fb_mean = positive_fb.mean(axis = 0) 
            positive_fb_mean = positive_fb_mean.reshape(1, sn, n) # добавляем ось для фидбека
            
        else:
            positive_fb_mean = np.empty((0, sn, n))
            
            
        ########################## Negative feedback #############################
        negative_fb = np.empty((0, sn, n))
        for r in rounds:
            try:
                               
                epochs_negative = np.load('/net/server/data/Archive/prob_learn/pultsinak/beta_16_30/stc/tfr/lateraloccipital_rh/tf_itc/{0}_run{1}_{2}_fb_cur_negative.npy'.format(subj, r, t))
                negative_fb = np.vstack([negative_fb, epochs_negative.reshape(1,20,1350)])
             
                
            except (OSError):
                print('This file not exist')

        ###### Шаг 1. Усреднили все отрицательные фидбеки внутри испытуемого (между блоками 1 -6) #################
        if negative_fb.size != 0:
            negative_fb_mean = negative_fb.mean(axis = 0) 
            negative_fb_mean = negative_fb_mean.reshape(1, sn, n) # добавляем ось для фидбека
        else:
            negative_fb_mean = np.empty((0, sn, n))
        ####################### Шаг 2 усреднения. Усредняем данные внутри испытуемого #####################################
        if negative_fb_mean.size == 0 and positive_fb_mean.size != 0:
            data_into_subj = positive_fb_mean
                
        elif negative_fb_mean.size != 0 and positive_fb_mean.size == 0:
                
            data_into_subj = negative_fb_mean
                
        elif negative_fb_mean.size != 0 and positive_fb_mean.size != 0:
                                        
            data_into_subj = np.vstack([negative_fb_mean, positive_fb_mean])
            
        else:
            data_into_subj = np.empty((0, sn, n))

        if data_into_subj.size != 0:
            temp.data = data_into_subj.mean(axis = 0)
        
            # сохраняем данные, усредненные внутри испытуемого. Шаг усредения 3, это усреднение между испытуемыми делается при рисовании топомапов
            np.save('/net/server/data/Archive/prob_learn/pultsinak/beta_16_30/stc/tfr/lateraloccipital_rh/tf_ave_itc/{0}_{1}'.format(subj, t),temp)
            
        else:
            print('Subject has no feedbacks in this condition')
            pass



subjects= ['P301','P304',"P307",'P312', 'P313', 'P314',
                                       'P316', 'P321','P322','P323','P324',"P325"]
                                      



autists= ['P301','P304',"P307",'P312', 'P313', 'P314',
                                       'P316', 'P321','P322','P323','P324',"P325",
                                       'P326','P327','P328', 'P329','P333','P334',
                                       'P334','P341','P342']

normal=['P001','P004','P019','P021','P022','P034','P035','P039',
'P040','P044','P047','P048','P053','P055','P058', 'P059',"P060", 'P061',
'P065','P063','P064','P067']
             
data_dir="/net/server/data/Archive/prob_learn/pultsinak/beta_16_30/stc/tfr/lateraloccipital_rh/tf_ave_itc"

temp = np.load('/net/server/data/Archive/prob_learn/pultsinak/beta_16_30/stc/tfr/lateraloccipital_rh/tf_ave_power/P301_norisk.npy')

n = temp.shape[1] # количество временных точек (берем у донора, если донор из тех же данных.
sn = temp.shape[0]                    
                    
                                                      

comp1_per_sub = np.zeros(shape=(len(autists), temp.shape[0], temp.shape[1]))
comp2_per_sub = np.zeros(shape=(len(autists), temp.shape[0], temp.shape[1]))
for ind, subj in enumerate(autists):
    temp1 = np.load(os.path.join(data_dir, "{0}_norisk.npy".format(subj)))
    comp1_per_sub[ind, :, :] = temp1
    
    temp2 = np.load(os.path.join(data_dir, "{0}_risk.npy".format(subj)))
    comp2_per_sub[ind, :, :] = temp2.data
        
    print(comp1_per_sub.shape)
    print(comp2_per_sub.shape)
hp_mean=comp1_per_sub.mean(axis=0)
lp_mean=comp2_per_sub.mean(axis=0)

tfr_diff= lp_mean - hp_mean
import matplotlib.pyplot as plt
times = epochs._raw_times

plt.subplots_adjust(0.1, 0.08, 0.96, 0.94, 0.2, 0.43)
plt.imshow(tfr_diff,extent=[times[0], times[-1], freqs[0], freqs[-1]],
        aspect="auto",
        origin="lower",
        vmin=-1,
        vmax=1,
        cmap="RdBu_r")
plt.xlabel("Time (s)")
plt.ylabel("Frequency (Hz)")

plt.colorbar()







