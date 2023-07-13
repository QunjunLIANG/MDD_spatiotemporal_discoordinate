#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script is used to visualize the time-delay projection maps on a brain template.

For a better view, I chosed to map them on the brain surface.

NOTE: 
    1. nilearn, panda and numpy are the required packages for running.
    2. The TD projection map is supposed to be estimated over the Schaefer's 400-roi Atlas.

Qunjun Liang 2023/01/14
"""

import pandas as pd
import numpy as np
from nilearn import datasets
from nilearn import plotting
from nilearn.maskers import NiftiLabelsMasker

# identiy files
wkDir = '/home/lqj/MDD_patient/NIfTI_convert/derivatives/timeseries_for_TD_MDD/time_lag_estimation/'
grp_projMap_file = wkDir + 'MDD_group_projection_map_weighted.csv' # HealthControl
grp_proflow_file = wkDir + 'MDD_group_probability_flow_projection_map.csv'
grp_projMap_file_unw = wkDir + 'MDD_group_projection_map_unweighted.csv'
grp_meanTD_file = wkDir + 'MDD_group_mean_lags.csv'

# load the group TD project map
projMap = pd.read_csv(grp_projMap_file, header=None)
projMap_array = np.array(projMap)

# load the group TD mean matrix
TDmean = pd.read_csv(grp_meanTD_file, header=None)
TDmean_array = np.array(TDmean)


# load the group TD project map unweighted
projMap_uw = pd.read_csv(grp_projMap_file_unw, header=None)
projMap_uw_array = np.array(projMap_uw)

# laod the group probability flow
prob_flow = pd.read_csv(grp_proflow_file, header=None)
prob_flow_array = np.array(prob_flow)

# obtain the atlas and make it as label masker
atlas_dat = datasets.fetch_atlas_schaefer_2018()
atlas_label = atlas_dat.labels
atlas_label = np.insert(atlas_label, 0, 'Background')

# Instantiate the masker with label image and label values
masker = NiftiLabelsMasker(atlas_dat.maps,
                           labels=atlas_label,
                           standardize=True)

# Note that we need to call fit prior to generating the mask
masker.fit()

# Visualize the TD projection maps

# unweighted TD projection map
td_uw_ima = masker.inverse_transform(projMap_uw_array)

# plot the weighted image 
#plotting.plot_img_on_surf(td_uw_ima, inflate=True, cmap='jet',
#                          title='unweighted Time delay projection map', vmax=0.07)
#plotting.show()

# weighted TD projection map
td_ima = masker.inverse_transform(projMap_array)

# plot the weighted image 
tde_proj_map = plotting.plot_img_on_surf(td_ima, inflate=True, cmap='jet',
                          title='weight Time delay projection map', vmax=0.07)
tde_proj_map_fig = tde_proj_map[0]
tde_proj_map_fig.savefig(wkDir+'weighted_TDE_projection_map.png')

# Visualize the probability flow maps

# inverse value 
pf_ima = masker.inverse_transform(prob_flow_array)

# plot the image 
prbflow_proj_map = plotting.plot_img_on_surf(pf_ima, inflate=True, cmap='jet',
                          title='probability flow projection map', vmax=0.07)
prbflow_proj_map_fig = prbflow_proj_map[0]
prbflow_proj_map_fig.savefig(wkDir+'probability_flow_projection_map.png')

# view the image in browser interactively
#plot_img = plotting.view_img(atlas_dat.maps, cmap='jet')
#plot_img.open_in_browser()

#view_img = plotting.plot_stat_map(td_ima, cmap='jet')