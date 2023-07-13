#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script is use to estimate the individual-level brain gradient. The template 
gradient is HC-group gradient which estimated by HC-group connectivity.

All individual gradients will be aligned to the template gradietn.

@author: lqj 2023/03/01
"""

import pandas as pd
import numpy as np
from nilearn.connectome import ConnectivityMeasure
from nilearn import plotting
from brainspace.gradient import GradientMaps
from brainspace.datasets import load_parcellation, load_conte69
import matplotlib.pyplot as plt
from sklearn.metrics import pairwise_distances
import os
# from brainspace.plotting import plot_hemispheres
# from brainspace.utils.parcellation import map_to_labels


# identify the inputs path ----------------------------------------------------
wkDir        = '/home/lqj/MDD_patient/NIfTI_convert/'
grad_dir   = wkDir + 'derivatives/REST/power264_timeseries/' 
outDir   = wkDir + 'derivatives/REST/individual_gradient/'
all_sbj_file =  wkDir + 'derivatives/REST_match_data_all.csv' 

timeseries_file  = 'xxxx_Power264.csv' # xxxx will be replaced by subject name

if not os.path.exists(outDir):
    os.makedirs(outDir)

template_conn = wkDir + 'derivatives/BrainGradient_alignment/Template_connectivity.csv'

# connectivity setting
con_method  = "covariance" # {"covariance", "correlation", "partial correlation", "tangent", "precision"}, 
cor_measure = ConnectivityMeasure(kind=con_method)

# gradient setting
grad_kernel = 'normalized_angle' # Possible options: {'pearson', 'spearman', 'cosine', 'normalized_angle', 'gaussian'}.
grad_method = 'dm' # {'dm', 'le', 'pca'} or object, optional
grad_align  = 'procrustes' # {'procrustes', 'joint'}, object or None
grad_Ncomp  = 80 # How much component will be extracted 

# obtain parcellations and mask ----------------------------------------------
labeling = load_parcellation('schaefer', scale=400, join=True)
mask = labeling != 0

## load hemisphere surfaces
surf_lh, surf_rh = load_conte69()

# function definitions -------------------------------------------------------

## function to estimate gradient 
def estimat_gradient(con_mat, method=grad_method, ref = None,
                     con_threshold = 0.9, kernel = grad_kernel,
                     align = grad_align, n_comp = grad_Ncomp):
    
    tmp_conn_z = np.tanh(con_mat)
    tmp_conn_z_pos = tmp_conn_z
    #tmp_conn_z_pos[tmp_conn_z_pos<0] = 0

    aff = tmp_conn_z_pos
    if grad_kernel is None: 
        aff=1-pairwise_distances(aff, metric='cosine')

    gm_tmp = GradientMaps(approach=method, 
                      kernel=kernel,
                      alignment=align,
                      n_components=n_comp)
    gm_tmp.fit(aff,
               sparsity=con_threshold,
               reference = ref)
    return gm_tmp

#########################################
#
# Make template gradient 
#
########################################

# estimate HC brain gradient.gradients_
temp_cor_mat = np.loadtxt(template_conn, delimiter=',')

gm_temp = estimat_gradient(temp_cor_mat, align = None)

print('Template gradient built!')

#########################################
#
# Estimate gradietn in individual
#
########################################

# load the subject ----------------------------------------------------------
sbj_tsv = pd.read_table(all_sbj_file)
sbj_name = sbj_tsv[['ID']]

ind = 1
for i in sbj_name.values:
    sbj = ''.join(i)
    print('Working in subject %s' % sbj)
    
    # replace the string pattern with subject name 
    timeseries_file_tmp = timeseries_file.replace('xxxx', sbj)
    
    if ind > 91 :
        timeseries_file_tmp = grad_dir + timeseries_file_tmp 
    else:
            timeseries_file_tmp = grad_dir + timeseries_file_tmp 
    
    # compute the mean of two arrays
    ts_tmp = np.loadtxt(timeseries_file_tmp, delimiter=',')
    
    # calculate group FC 
    cor_mat_tmp = cor_measure.fit_transform([ts_tmp])[0]
    
    # export the group-average matrix
    np.fill_diagonal(cor_mat_tmp, 0)
    cor_mat_tmp_out = pd.DataFrame(cor_mat_tmp)
    cor_mat_tmp_out.to_csv(outDir + sbj + '_connect_matrix.csv',
                           index=False, header = False)
        
    # plot the correlation matrix and save to disk
    corr_plot = plotting.plot_matrix(cor_mat_tmp_out, 
                                    figure=(15,15),
                                    auto_fit=True, vmax = 1, vmin = -1)
    corr_plot.write_png(outDir + sbj + '_connect_matrix.png')
    
    # estimate gradient ------------------------------------------------
    gm_tmp = estimat_gradient(cor_mat_tmp, ref = gm_temp.gradients_) # align to HC gradient
    
    ## show and save the lambda plot for HC and MDD respectively
    lambdas = gm_tmp.lambdas_
    fig, ax = plt.subplots(1, figsize=(5,4))
    ax.scatter(range(lambdas.size), lambdas/sum(lambdas))
    ax.set_xlabel('Component_Nb')
    ax.set_ylabel('Eigenvalue')
    plt.savefig(outDir+('%s_lambda_kernel-%s_method-%s.png' % (sbj, grad_kernel, grad_method)))
    plt.close()

    
    ## export the scatter plot  
    gradients = gm_tmp.gradients_
    cart = (gradients - np.mean(gradients, axis=0)) / np.max(np.abs(gradients - np.mean(gradients)))
    th, r = np.arctan2(cart[:, 1], cart[:, 0]), np.linalg.norm(cart, axis=1) / 2 + .5
    C = np.hstack([np.cos(.75 * (th + 0 * np.pi))[:, None],
                   np.cos(.75 * th - .5 * np.pi)[:, None], 
                   np.cos(.75 * th + .5 * np.pi)[:, None]])
    C = C * r[:, None] / np.max(C)
    C[C < 0] = 0
    C /= np.max(C)
    fig, ax = plt.subplots(figsize=(9, 9))
    scatter = ax.scatter(gradients[:, 0], gradients[:, 1], s=200, c=C, marker='.')
    ax.set_xlabel('Gradient 1')
    ax.set_ylabel('Gradient 2')
    plt.savefig(outDir + ('%s_scatter_kernel-%s_method-%s.png' % (sbj, grad_kernel, grad_method)))
    plt.close()
    
    ## write the gradient to outer csv
    dat_dic = {'gradient1':gm_tmp.gradients_[:,0],
               'gradient2':gm_tmp.gradients_[:,1],
               'gradient3':gm_tmp.gradients_[:,2]}
    dat_gradient = pd.DataFrame(dat_dic)
    dat_gradient.to_csv(outDir + sbj + '_gradient_export.csv', index=None)
    
    # ## plot first three gradients 
    # gradient_group = [None] * 3

    # for i in range(3):
    #     gradient_group[i] = map_to_labels(gm_tmp.gradients_[:,i], labeling, mask=mask,
    #                                       fill=np.nan)

    # label_txt = ['gradient1','gradient2', 'gradient3']
    # plot_hemispheres(surf_lh, surf_rh, array_name=gradient_group, size = (1200,600),
    #                   cmap='viridis_r', color_bar=True, label_text=label_txt,
    #                   zoom = 1.2,
    #                   screenshot=True, transparent_bg = False,
    #                   filename = outDir + sbj + '_Gradient_visual.jpg')
    
    print('subject %s finished!' % sbj)
    
    ind = ind + 1


