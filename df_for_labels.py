
import mne
import numpy as np
import pandas as pd
import os


os.environ['SUBJECTS_DIR'] = '/net/server/data/Archive/prob_learn/freesurfer'
subjects_dir = '/net/server/data/Archive/prob_learn/freesurfer'

fsaverage = mne.setup_source_space(subject = "fsaverage", spacing="ico5", add_dist=False)


subjects = pd.read_csv('/home/vtretyakova/Рабочий стол/probability_learning/subj_list.csv')['subj_list'].tolist()
subjects.remove('P062') 
subjects.remove('P052') 
subjects.remove("P032")
subjects.remove('P045') 

rounds = [1, 2, 3, 4, 5, 6]
freq_range = "beta_16_30"
trial_type = ['norisk', 'prerisk', 'risk', 'postrisk']
feedback = ['positive', 'negative']




#parc that we used https://balsa.wustl.edu/WN56
labels =  mne.read_labels_from_annot("fsaverage", "aparc_sub", hemi = "both")
label_names = [label.name for label in labels] 

labels.pop(448)###### delete unknown labels !!!!
label_names = [label.name for label in labels] 

data_path = "/net/server/data/Archive/prob_learn/pultsinak/beta_16_30/stc/stc_epo_fsaverage"

for subj in subjects:
    print(subj)
    df = pd.DataFrame()
    
    for r in rounds:
        for cond in trial_type:
            for fb in feedback:
                    
                try:
                    epochs_num = os.listdir(os.path.join(data_path, '{0}_run{1}_{2}_fb_cur_{3}_fsaverage'.format(subj, r, cond, fb)))
                    epo_n = (int(len(epochs_num) / 2))
                      
                    for ep in range(epo_n):
                        df_epo= pd.DataFrame()
                        print(ep)
                        stc = mne.read_source_estimate(os.path.join(data_path, '{0}_run{1}_{2}_fb_cur_{3}_fsaverage/{4}'.format(subj, r, cond, fb, ep)))
                        stc2 = stc.copy()
                        stc2=stc2.crop(tmin=1.500, tmax=1.900, include_tmax=True) ### crop the time what you want to analyse

                        label_ts = mne.extract_label_time_course(stc2,labels, src=fsaverage, mode='mean')
                        label_ts_avg=label_ts.mean(axis=1)
                        
                        epo = [ep for i in range(448)]
                        subject = [subj for i in range(448)]
                        run = [r for i in range(448)]
                        trial = [cond for i in range(448)]
                        fb_cur=[fb for i in range(448)]
                        
                        df_epo['beta_power'] = label_ts_avg
                        df_epo['label'] = label_names
                        df_epo['epo'] = epo
                        df_epo['subject'] = subject
                        df_epo['round'] = run
                        df_epo['trial_type'] = trial
                        df_epo['feedback_cur'] = fb_cur
                            
                        df = df.append(df_epo)    
                         
                                       
                except (OSError, FileNotFoundError):
                    print('This file not exist')
    df.to_csv('/net/server/data/Archive/prob_learn/pultsinak/label_stc/df_1500_1900_khan/{0}.csv'.format(subj))
                            
                            
