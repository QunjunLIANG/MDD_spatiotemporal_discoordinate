#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
This Python script is used to extract the timeseries from functional
images preprocessed by fmriprep.

This pipeline will generate timeseries files, framewise displacement record 
and the labels of the used Atlas.

Schaefer 400 Atlas is used by default.

Note!!!!!!!!!!!!!!!
If the timeseries is extracted for computing time delay, the standard_use
should set as False, using the raw timeseries.
If the timeseries is extracted for computing brain gradient, the standard_use
should set as True, using the z-scored timeseries.

Qunjun Liang 2023/01/10

"""

import os
import pandas as pd
# import numpy as np
from nilearn.interfaces.fmriprep import load_confounds_strategy
from nilearn import datasets
from nilearn.maskers import NiftiLabelsMasker
from multiprocessing.dummy import Pool as ThredPool

# identify the path 
wkDir        = '/home/lqj/MDD_patient/NIfTI_convert/'
sbj_file     = 'MDD_validation_dataset.tsv' # subject file in BIDS

second_dir   = 'derivatives/fmriprep_volumn/xxxx/func/' # xxxx will be replaced by subject name
funcImg_name = 'xxxx_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'
conf_name    = 'xxxx_task-rest_desc-confounds_timeseries.tsv' # confound file

# export timeseries
outDir       = wkDir + 'derivatives/validation_timeseries_for_TD/'
FD_out_file  = outDir + 'xxxx_FD.csv'
TS_out_file  = outDir + 'xxxx_timeseries_schaefer400.csv'

# setting the parameters
para_detrend  = True # detrend the timeseries?
para_lowPass  = 0.1 # low pass in Hz
para_highPass = 0.005 # high pass in Hz
para_TR       = 2 # TR in second
para_smooth   = 8 # full-width at half maximum, in mm
standard_use  = False # extract the raw timeseries or not
para_denosie  = "simple" # one of "simple", "srubbing", "compcor", "ica_aroma"

# obtain parcellations and make it into a labeling masker
atlas_dataset = datasets.fetch_atlas_schaefer_2018()
atlas_map = atlas_dataset.maps

# make sure the output directory is avaliable 
if not os.path.exists(outDir):
    os.makedirs(outDir)

masker = NiftiLabelsMasker(
    labels_img=atlas_map, labels = atlas_dataset.labels,
    detrend=para_detrend, standardize=standard_use, smoothing_fwhm=para_smooth,
    low_pass=para_lowPass, high_pass=para_highPass, t_r=para_TR)

# export the ROI definition for Yeo 7 netw/media/mri/research_data/MDD_patient_kangning/orks
atlas_label = atlas_dataset.labels

fixedList = []
for xind in atlas_label:
    fixedList.append(xind.decode('utf-8'))
    
roiName = pd.DataFrame(fixedList)
roiName.to_csv(outDir+'Yeo_7net_roiAnnotation.csv', index=False, header=False)

# define the function
def extract_tc(sub_name):
    sbj = ''.join(sub_name)
    # replace the string pattern with subject name 
    funcImg_tmp = funcImg_name.replace('xxxx', sbj)
    second_dir_tmp = second_dir.replace('xxxx', sbj)
    funcImg_path_tmp = wkDir + second_dir_tmp + funcImg_tmp
    
    # use ica_aroma strategy
    confounds, sample_mask = load_confounds_strategy(
        funcImg_path_tmp, denoise_strategy=para_denosie)
    
    # extract the timeseries 
    time_series = masker.fit_transform(funcImg_path_tmp, 
                                       confounds=confounds,
                                       sample_mask=None)
    
    # export the time series as dataframe in CSV format
    TS_out_file_tmp = TS_out_file.replace('xxxx', sbj)
    
    time_series_out = pd.DataFrame(time_series)

    time_series_out.to_csv(TS_out_file_tmp, index=False, header=False)

    # obtain the FD for this subject's data
    conf_name_tmp = conf_name.replace('xxxx', sbj)
    conf_sbj = pd.read_table(wkDir+second_dir_tmp+conf_name_tmp).fillna(0) # fill the Nan with 0
    conf_FD = conf_sbj[['framewise_displacement']]
    
    # export the FD 
    FD_out_file_tmp = FD_out_file.replace('xxxx', sbj)
    conf_FD.to_csv(FD_out_file_tmp, index=False, header=False)
    
    print('Work finished in subject %s' % sbj)
    
    
# obtaini subject name 
sbj_tsv = pd.read_table(wkDir+sbj_file)
sbj_name = sbj_tsv[['participant_id']]

# export the subject name 
sbj_name.to_csv(outDir+'subject_name.csv', index=False, header=False)


pool = ThredPool()
pool.map(extract_tc, sbj_name.values) 
pool.close()
pool.join()   




