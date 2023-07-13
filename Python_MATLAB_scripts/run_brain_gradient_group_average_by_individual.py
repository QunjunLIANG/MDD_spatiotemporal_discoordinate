#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script is use to estimate the group level brain gradient bsaed on 
the individual gradient estimated previously

All individual gradients have been aligned to the template gradietn.

@author: lqj 2023/03/01
"""

import pandas as pd
import numpy as np
from brainspace.datasets import load_parcellation, load_conte69
import matplotlib.pyplot as plt
import os
from brainspace.plotting import plot_hemispheres
from brainspace.utils.parcellation import map_to_labels

# indicate the path ----------------------------------------------------------
wkDir = '/home/lqj/MDD_patient/NIfTI_convert/'
grad_path = wkDir + 'derivatives/BrainGradient_alignment/individual_gradient/'
outDir = wkDir + 'derivatives/BrainGradient_alignment/group_mean/'

hc_sbj_file  = wkDir + 'HC_match_data.tsv' 
mdd_sbj_file  = wkDir + 'MDD_match_data.tsv' 

if not os.path.exists(outDir):
    os.makedirs(outDir)

gradient_file  = 'xxxx_gradient_export.csv' # xxxx will be replaced by subject name

# obtain parcellations and mask ----------------------------------------------
labeling = load_parcellation('schaefer', scale=400, join=True)
mask = labeling != 0

## load hemisphere surfaces
surf_lh, surf_rh = load_conte69()

# define function to load the individual gradient and calculate the group-mean
def Calculate_group_mean(sbj_name):
    
    sbj_grad_collect = np.ones((91,400,3))
    
    for i in range(sbj_name.shape[0]):
        
        sbj = ''.join(sbj_name.values[i])
        print('Working in subject %s' % sbj)
        
        # replace the string pattern with subject name 
        gradient_file_tmp = gradient_file.replace('xxxx', sbj)
        gradient_file_tmp = grad_path + gradient_file_tmp
        
        # compute the mean of two arrays
        grad_tmp = pd.read_csv(gradient_file_tmp, delimiter=',')
        sbj_grad_collect[i,:,:] = grad_tmp
        
    group_mean_grad = np.mean(sbj_grad_collect, axis = 0)
    
    return group_mean_grad


# load the subject ----------------------------------------------------------
## MDD group
sbj_tsv = pd.read_table(mdd_sbj_file)
sbj_name = sbj_tsv[['participant_id']]

mdd_group_mean = Calculate_group_mean(sbj_name)

dat_dic = {'gradient1':mdd_group_mean[:,0],
           'gradient2':mdd_group_mean[:,1],
           'gradient3':mdd_group_mean[:,2]}
dat_gradient = pd.DataFrame(dat_dic)
dat_gradient.to_csv(outDir + 'MDD_groupmean_gradient_export.csv', index=None)

## HC group
sbj_tsv = pd.read_table(hc_sbj_file)
sbj_name = sbj_tsv[['participant_id']]

hc_group_mean = Calculate_group_mean(sbj_name)

dat_dic = {'gradient1':hc_group_mean[:,0],
           'gradient2':hc_group_mean[:,1],
           'gradient3':hc_group_mean[:,2]}
dat_gradient = pd.DataFrame(dat_dic)
dat_gradient.to_csv(outDir + 'HC_groupmean_gradient_export.csv', index=None)

## export the scatter plot  
fig = plt.scatter(hc_group_mean[:,0], hc_group_mean[:,1])
plt.savefig(outDir + 'HC_2D_scatter.png')
plt.show()

fig = plt.scatter(mdd_group_mean[:,0], mdd_group_mean[:,1])
plt.savefig(outDir + 'MDD_2D_scatter.png')
plt.show()
# plot first three gradients ----------------------------------------------
## for MDD
gradient_group = [None] * 3

for i in range(3):
    gradient_group[i] = map_to_labels(mdd_group_mean[:,i], labeling, mask=mask,
                                      fill=np.nan)

label_txt = ['gradient1','gradient2', 'gradient3']
plot_hemispheres(surf_lh, surf_rh, array_name=gradient_group, size = (1200,600),
                  cmap='viridis_r', color_bar=True, label_text=label_txt,
                  zoom = 1.2,
                  screenshot=True, transparent_bg = False,
                  filename = outDir + 'MDD_Gradient_visual.jpg')

## for HC
gradient_group = [None] * 3

for i in range(3):
    gradient_group[i] = map_to_labels(hc_group_mean[:,i], labeling, mask=mask,
                                      fill=np.nan)

label_txt = ['gradient1','gradient2', 'gradient3']
plot_hemispheres(surf_lh, surf_rh, array_name=gradient_group, size = (1200,600),
                  cmap='viridis_r', color_bar=True, label_text=label_txt,
                  zoom = 1.2,
                  screenshot=True, transparent_bg = False,
                  filename = outDir + 'HC_Gradient_visual.jpg') 