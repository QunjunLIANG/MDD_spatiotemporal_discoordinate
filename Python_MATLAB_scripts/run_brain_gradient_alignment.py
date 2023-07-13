#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script is a update for previous brain gradient estimation pipeline.

The main changes are as follows:
    1. using the alignment method to match the MDD gradient to HC gradient;
    2. using the plotting feature in brainspace instread of nilearn for a better view.

This pipeline includes 4 features:
    1. calculate the connectivity matrix and export to the png;
    2. calculate the gradient map;
    3. export the first three components of the gradients to surface;
    4. export the first three components of the gradients to csv;
    5. plot and export the lambda plot.

Note: 
    1. Please make sure the for loop in subject file is decent before running.
    2. The HC gradient is gm.gradient_[0] while to MDD is gm.gradient_[1].

@author: lqj 2023/02/28
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
# identify the inputs path ----------------------------------------------------
wkDir        = '/home/lqj/MDD_patient/NIfTI_convert/'
grad_dir   = wkDir + 'derivatives/timeseries_for_gradient/' 
outDir   = wkDir + 'derivatives/BrainGradient_alignment/'
hc_sbj_file     = wkDir + 'HC_match_data.tsv' 
mdd_sbj_file  = wkDir + 'MDD_match_data.tsv' 

timeseries_file  = 'xxxx_timeseries_schaefer400.csv' # xxxx will be replaced by subject name

if not os.path.exists(outDir):
    os.makedirs(outDir)
    
# gradient setting
connective_method = "correlation" # {"covariance", "correlation", "partial correlation", "tangent", "precision"}, 
grad_kernel = 'normalized_angle' # Possible options: {'pearson', 'spearman', 'cosine', 'normalized_angle', 'gaussian',None}.
grad_method = 'dm' # {'dm', 'le', 'pca'} or object, optional
grad_align  = 'procrustes' # {'procrustes', 'joint'}, object or None
grad_Ncomp  = 80 # How much component will be extracted 

# obtain parcellations and mask -------------------------------------------
labeling = load_parcellation('schaefer', scale=400, join=True)
mask = labeling != 0

## load hemisphere surfaces
surf_lh, surf_rh = load_conte69()

# function definitions ------------------------------------------------
## function to calculate group-mean connectivity
def obtain_group_conn(sbjname, inputDir,
                      timeseries_file=timeseries_file,
                      con_method=connective_method):
    
    # read the subject's timeseries
    sbj_timeseries_collect = []
    for i in sbjname:
        sbj = ''.join(i)
        print('Working in subject %s' % sbj)
        
        # replace the string pattern with subject name 
        timeseries_file_tmp = timeseries_file.replace('xxxx', sbj)
        timeseries_file_tmp = inputDir + timeseries_file_tmp 
        
        # compute the mean of two arrays
        ts_tmp = np.loadtxt(timeseries_file_tmp, delimiter=',')
        sbj_timeseries_collect.append(ts_tmp)

    # calculate group FC 
    cor_measure = ConnectivityMeasure(kind=con_method)
    cor_mat = cor_measure.fit_transform(sbj_timeseries_collect) # subject-level FC
    grp_cor_mat = np.mean(cor_mat, axis=0)
    return grp_cor_mat

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
                      alignment=align, n_components = n_comp)
    gm_tmp.fit(aff,
               sparsity=con_threshold,
               reference = ref)
    return gm_tmp

#########################################
#
# obtain template gradient from outer sample (Vos de Wael et al., 2018, PNAS)
#
########################################

# load the subject ----------------------------------------------------------
from brainspace.datasets import load_group_fc

con_temp = load_group_fc('schaefer', scale = 400)

#########################################
#
# Make connectivity matrix in Health control 
#
########################################

# load the subject ----------------------------------------------------------
sbj_tsv = pd.read_table(hc_sbj_file)
sbj_name = sbj_tsv[['participant_id']]

# export the group-average matrix
hc_grp_cor_mat = obtain_group_conn(sbjname=sbj_name.values, inputDir=grad_dir)

np.fill_diagonal(hc_grp_cor_mat, 0)
cor_matrix_out = pd.DataFrame(hc_grp_cor_mat)
cor_matrix_out.to_csv(outDir + 'HC_groupmean_connectivity.csv', index=False, header = False)
    
# plot the correlation matrix and save to disk
corr_plot = plotting.plot_matrix(hc_grp_cor_mat, 
                                 figure=(15,15),
                                 auto_fit=True, vmax = 1, vmin = -1)
corr_plot.write_png(outDir + 'HC_groupmean_connectivity.png')


#########################################
#
# Make connectivity matrix in Health control 
#
########################################

# load the subject ----------------------------------------------------------
sbj_tsv = pd.read_table(mdd_sbj_file)
sbj_name = sbj_tsv[['participant_id']]

# export the group-average matrix
mdd_grp_cor_mat = obtain_group_conn(sbjname=sbj_name.values, inputDir=grad_dir)
np.fill_diagonal(mdd_grp_cor_mat, 0)
cor_matrix_out = pd.DataFrame(mdd_grp_cor_mat)
cor_matrix_out.to_csv(outDir + 'MDD_groupmean_connectivity.csv', index=False, header = False)
    
# plot the correlation matrix and save to disk
corr_plot = plotting.plot_matrix(mdd_grp_cor_mat, 
                                 figure=(15,15),
                                 auto_fit=True, vmax = 1, vmin = -1)
corr_plot.write_png(outDir + 'MDD_groupmean_connectivity.png')


#########################################
#
# estimate aligment brain gradient 
#
########################################

from brainspace.plotting import plot_hemispheres
from brainspace.utils.parcellation import map_to_labels

# constructure the template connectivity
# con_temp = (hc_grp_cor_mat + mdd_grp_cor_mat)/2 # template
# np.fill_diagonal(con_temp, 0)
con_temp_out = pd.DataFrame(con_temp)
con_temp_out.to_csv(outDir + 'Template_connectivity.csv', index=False, header = False)
    
# plot the correlation matrix and save to disk
corr_plot = plotting.plot_matrix(con_temp, 
                                 figure=(15,15),
                                 auto_fit=True, vmax = 1, vmin = -1)
corr_plot.write_png(outDir + 'Template_connectivity.png')

# estiamte the brain gradient 
gm_temp = estimat_gradient(con_temp, align=None)

gm_hc = estimat_gradient(hc_grp_cor_mat, ref = gm_temp.gradients_)
gm_mdd = estimat_gradient(mdd_grp_cor_mat, ref = gm_temp.gradients_) # align to HC gradient

## plot the first gradient between HC and MDD ---------------------------------
gradient_group = [None] * 3
gradient_group[0] = map_to_labels(gm_temp.gradients_[:,0],
                                  labeling, mask=mask,
                                  fill=np.nan)
gradient_group[1] = map_to_labels(gm_hc.gradients_[:,0],
                                  labeling, mask=mask,
                                  fill=np.nan)
gradient_group[2] = map_to_labels(gm_mdd.gradients_[:,0],
                                  labeling, mask=mask,
                                  fill=np.nan)

label_txt = ['Template','HC','MDD']
plot_hemispheres(surf_lh, surf_rh, array_name=gradient_group, size = (1200,600),
                 cmap='viridis_r', color_bar=True, label_text=label_txt,
                 zoom = 1.2,
                 screenshot=True, transparent_bg = False,
                 filename = outDir+'Gradient1_visual.jpg')

## plot the second gradient between HC and MDD ---------------------------------
gradient_group = [None] * 3
gradient_group[0] = map_to_labels(gm_temp.gradients_[:,1],
                                  labeling, mask=mask,
                                  fill=np.nan)
gradient_group[1] = map_to_labels(gm_hc.gradients_[:,1],
                                  labeling, mask=mask,
                                  fill=np.nan)
gradient_group[2] = map_to_labels(gm_mdd.gradients_[:,1],
                                  labeling, mask=mask,
                                  fill=np.nan)

label_txt = ['Template','HC','MDD']
plot_hemispheres(surf_lh, surf_rh, array_name=gradient_group, size = (1200,600),
                 cmap='viridis_r', color_bar=True, label_text=label_txt,
                 zoom = 1.2,
                 screenshot=True, transparent_bg = False,
                 filename = outDir+'Gradient2_visual.jpg')

## plot the third gradient between HC and MDD ---------------------------------
gradient_group = [None] * 3
gradient_group[0] = map_to_labels(gm_temp.gradients_[:,2],
                                  labeling, mask=mask,
                                  fill=np.nan)
gradient_group[1] = map_to_labels(gm_hc.gradients_[:,2],
                                  labeling, mask=mask,
                                  fill=np.nan)
gradient_group[2] = map_to_labels(gm_mdd.gradients_[:,2],
                                  labeling, mask=mask,
                                  fill=np.nan)

label_txt = ['Template','HC','MDD']
plot_hemispheres(surf_lh, surf_rh, array_name=gradient_group, size = (1200,600),
                 cmap='viridis_r', color_bar=True, label_text=label_txt,
                 zoom = 1.2,
                 screenshot=True, transparent_bg = False,
                 filename = outDir+'Gradient3_visual.jpg')

# show and save the lambda plot for HC and MDD respectively --------------------
lambdas = gm_temp.lambdas_
fig, ax = plt.subplots(1, figsize=(5,4))
ax.scatter(range(lambdas.size), lambdas/sum(lambdas))
ax.set_xlabel('Component_Nb')
ax.set_ylabel('Eigenvalue')
plt.savefig(outDir+('Template_lambda_kernel-%s_method-%s.png' % (grad_kernel, grad_method)))
plt.close()

lambdas_tmp = gm_hc.lambdas_
fig, ax = plt.subplots(1, figsize=(5,4))
ax.scatter(range(lambdas.size), lambdas/sum(lambdas))
ax.set_xlabel('Component_Nb')
ax.set_ylabel('Eigenvalue')
plt.savefig(outDir+('HC_lambda_kernel-%s_method-%s.png' % (grad_kernel, grad_method)))
plt.close()

lambdas_tmp = gm_mdd.lambdas_
fig, ax = plt.subplots(1, figsize=(5,4))
ax.scatter(range(lambdas.size), lambdas/sum(lambdas))
ax.set_xlabel('Component_Nb')
ax.set_ylabel('Eigenvalue')
plt.savefig(outDir+('MDD_lambda_kernel-%s_method-%s.png' % (grad_kernel, grad_method)))
plt.close()

# export the scatter plot for HC and MDD --------------------------------------
gradients = gm_temp.gradients_
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
plt.savefig(outDir + ('Template_scatter_kernel-%s_method-%s.png' % (grad_kernel, grad_method)))
plt.close()

gradients = gm_hc.gradients_
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
plt.savefig(outDir + ('HC_scatter_kernel-%s_method-%s.png' % (grad_kernel, grad_method)))
plt.close()

gradients = gm_mdd.gradients_
gradients[:,0] *= -1
gradients[:,1] *= -1
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
plt.savefig(outDir + ('MDD_scatter_kernel-%s_method-%s.png' % (grad_kernel, grad_method)))
plt.close()

# write the component to the outer csv ----------------------------------------
## for template
dat_dic = {'gradient1':gm_temp.gradients_[:,0],
           'gradient2':gm_temp.gradients_[:,1],
           'gradient3':gm_temp.gradients_[:,2]}
dat_gradient = pd.DataFrame(dat_dic)
dat_gradient.to_csv(outDir + 'Template_gradient_export.csv', index=None)

## for HC's
dat_dic = {'gradient1':gm_hc.gradients_[:,0],
           'gradient2':gm_hc.gradients_[:,1],
           'gradient3':gm_hc.gradients_[:,2]}
dat_gradient = pd.DataFrame(dat_dic)
dat_gradient.to_csv(outDir + 'HC_gradient_export.csv', index=None)

## for MDD's
dat_dic = {'gradient1':gm_mdd.gradients_[:,0],
           'gradient2':gm_mdd.gradients_[:,1],
           'gradient3':gm_mdd.gradients_[:,2]}
dat_gradient = pd.DataFrame(dat_dic)
dat_gradient.to_csv(outDir + 'MDD_gradient_export.csv', index=None)