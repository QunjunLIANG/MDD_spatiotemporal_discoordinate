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
import glob
from multiprocessing.dummy import Pool as ThredPool
# identify the inputs path ----------------------------------------------------
wkDir        = '/home/lqj/MDD_patient/NIfTI_convert/'
tde_dir      = wkDir + 'derivatives/REST/time_lag_estimation/' 
outDir       = wkDir + 'derivatives/REST/TD_gradient_alignment/'
sbj_file     = '/home/lqj/MDD_patient/R_project/tde_gradient_coupling_in_-mdd/REST_match_data_all.csv' 

tde_file  = 'xxxx_timeDelay.csv' # xxxx will be replaced by subject name

# load the group-mean td matrix for template
grp_tde_mat_hc = tde_dir + 'HC_group_mean_lags.csv'
grp_tde_mat_mdd = tde_dir + 'MDD_group_mean_lags.csv'

if not os.path.exists(outDir):
    os.makedirs(outDir)
    
# gradient setting
grad_kernel = 'normalized_angle' # Possible options: {'pearson', 'spearman', 'cosine', 'normalized_angle', 'gaussian',None}.
grad_method = 'dm' # {'dm', 'le', 'pca'} or object, optional
grad_align  = 'procrustes' # {'procrustes', 'joint'}, object or None
grad_Ncomp  = 80 # How much component will be extracted 

outLambda  = outDir + ('xxxx_lambdas_kernel-%s_method-%s_align-%s.png' % (grad_kernel, grad_method, grad_align))
outScatter = outDir + ('xxxx_scatter_kernel-%s_method-%s_align-%s.png' % (grad_kernel, grad_method, grad_align))
outGrad    = outDir + ('xxxx_gradient_kernel-%s_method-%s_align-%s.csv' % (grad_kernel, grad_method, grad_align))


# function definitions ------------------------------------------------
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

    
####################################
#
# estimate the group template
#
####################################

# step 1: load the MDD and HC group matrix
td_grp_hc = np.loadtxt(grp_tde_mat_hc,delimiter = ',')
td_grp_mdd = np.loadtxt(grp_tde_mat_mdd,delimiter = ',')
grp_tde_mat = (td_grp_hc + td_grp_mdd)/2

# step 2: initalize the Brain gradient function
gm_temp = GradientMaps(approach=grad_method, 
                      kernel=grad_kernel, 
                      n_components = grad_Ncomp)
gm_temp.fit(grp_tde_mat)

colnames = [f'gradint{i}' for i in range(1,grad_Ncomp+1)]
dat_gradient = pd.DataFrame(gm_temp.gradients_, columns=colnames)
dat_gradient.to_csv(outDir + ('Template_gradient_kernel-%s_method-%s.csv' % (grad_kernel, grad_method)), index=None)

# step 3: output the lambda plot
lambdas = gm_temp.lambdas_
fig, ax = plt.subplots(1, figsize=(5,4))
ax.scatter(range(lambdas.size), lambdas/sum(lambdas))
ax.set_xlabel('Component_Nb')
ax.set_ylabel('Eigenvalue')
plt.savefig(outDir+('Template_lambda_kernel-%s_method-%s.png' % (grad_kernel, grad_method)))
plt.close()

# step 4: output the scattter plot for the first two gradients
gradients = gm_temp.gradients_
gradients[:,0] *= -1
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

###########################################
#
# estimate in subjects
#
#########################################

def EstimateIndividualGradient(sbjname):
    
    sbj = ''.join(sbjname)
    print('Working in subject %s' % sbj)
    
    #####################################
    # rename the path  
    
    # rename the input connectivity
    corMat_file_tmp = tde_file.replace('xxxx', sbj)
    
    # rename the output
    outLambda_tmp  = outLambda.replace('xxxx', sbj)
    outScatter_tmp = outScatter.replace('xxxx', sbj)
    outGrad_tmp    = outGrad.replace('xxxx', sbj)

    if not os.path.isfile(outGrad_tmp):
		###################################
		# load the connectivity matrxi 
        conMat_tmp = np.loadtxt(tde_dir + corMat_file_tmp, delimiter = ',')
        conMat_tmp = np.nan_to_num(conMat_tmp)  # replace nan to 0
		##################################
		# estimate the gradietn 
		# step 1: Using fisher-z transform to the connectivity matrix
        conMat_tmp = np.arctan(conMat_tmp)

		# step 2: initalize the Brain gradient function
        gm = GradientMaps(approach=grad_method, 
		                      kernel=grad_kernel, 
		                      alignment=grad_align,
		                      n_components = grad_Ncomp)

		# step 3: calling fit to the group connectivty
        gm.fit(conMat_tmp, reference=gm_temp.gradients_)
		
		# step 4: export to the csv 
        dat_dic = {'gradient1':gm.gradients_[:,0],
		           'gradient2':gm.gradients_[:,1],
		           'gradient3':gm.gradients_[:,2]}
        dat_gradient = pd.DataFrame(dat_dic)
        dat_gradient.to_csv(outGrad_tmp, index=None)
		
		# step 5: output the lambda plot
        lambdas = gm.lambdas_
        fig, ax = plt.subplots(1, figsize=(5,4))
        ax.scatter(range(lambdas.size), lambdas/sum(lambdas))
        ax.set_xlabel('Component number')
        ax.set_ylabel('Scaled eigenvalue')
        fig.savefig(outLambda_tmp)
        plt.close(fig)
		
		# step 6: output the scattter plot for the first two gradients
        gradients = gm.gradients_
        gradients[:,0] *= -1
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
        plt.savefig(outScatter_tmp)
        plt.close()
		
        print('subject %s finished' % sbj)
    else:
        print('subject already finished')

# load the subject ----------------------------------------------------------
sbj_tsv = pd.read_table(sbj_file)
sbj_name = sbj_tsv[['ID']]

pool = ThredPool()
pool.map(EstimateIndividualGradient, sbj_name.values) 
pool.close()
pool.join()   