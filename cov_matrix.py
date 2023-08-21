import mne
import os
import os.path as op
import numpy as np
import pandas as pd


os.environ['SUBJECTS_DIR'] = '/net/server/data/Archive/prob_learn/freesurfer'
subjects_dir = '/net/server/data/Archive/prob_learn/freesurfer'
subjects= ['P301','P304',"P307",'P312', 'P313', 'P314',
                                       'P316', 'P321','P322','P323','P324',"P325",
                                       'P326','P327','P328', 'P329','P333','P335',
                                       'P334','P341','P342', 'P001','P004','P019','P021','P022','P034','P035','P039',
                                      'P040','P044','P047','P048','P053','P055','P058', 'P059',"P060", 'P061',
                                      'P065','P063','P064','P067']

date  =['161228','170201', '170330', '170404','170406', '200807','200819','201002','201016', '201106','201120','201120',
         '201130','201204','201205','201205','201209','201218','220204','211203','211222','220706','210226',
         
         '210326','210510','210402','210430','210423','210423','210430','210508','210522','210604',
         '210604','210605','210625','210625','210709','211110','211108','211108','220318','220527']



for subj in subjects:
    for d in date:
        try:
            raw_er = mne.io.Raw('/net/server/data/Archive/prob_learn/experiment/{0}/{1}/RAW/{0}_er_raw.fif'.format(subj,d), preload=True)

            #rank =mne.compute_rank(raw_data)
            cov = mne.compute_raw_covariance(raw_er, method='shrunk')
            
            cov.plot(raw_er.info, exclude="bads", show_svd=False)
            
            
            cov.save('/net/server/data/Archive/prob_learn/pultsinak/cov_empty_room_shrunk/{0}_er-cov.fif'.format(subj))
        except (OSError):
            print('This file not exist')
