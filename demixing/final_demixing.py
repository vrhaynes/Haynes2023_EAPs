save_plots = True


# for visualization in WSL
import os
import matplotlib.pyplot as plt


# General modules
import sys
import numpy as np
import pandas as pd
from pprint import pprint
import pickle
import seaborn as sns

# Analysis modules
import pywt
from pywt import wavedec, waverec
from scipy.interpolate import CubicSpline
import tensortools as tt
import RegscorePy as regscore
from scipy.signal import find_peaks

import analysis as a
from analysis import *

cwd = os.getcwd()
parent_dir = os.path.dirname(cwd)
sys.path.append(parent_dir)
print('Added\n --> %s \nto path.' %parent_dir)

data_dir = os.path.join(parent_dir,'model_data','bbp')
path_to_figs = os.path.join(parent_dir,'figures')
path_to_save = os.path.join(parent_dir,'analyses')

test_demixing = False
saveModel = True
R = 4

# 12 excitatory types
E_labels =  ['L5_TTPC1_%s'%(i+1) for i in range(5)]+  \
                ['L5_TTPC2_%s'%(i+1) for i in range(5)]+   \
                ['L5_UTPC_%s'%(i+1) for i in range(5)]+  \
                ['L5_STPC_%s'%(i+1) for i in range(5)]+  \
                ['L6_TPC_L4_%s'%(i+1) for i in range(5)]+\
                ['L6_TPC_L1_%s'%(i+1) for i in range(5)]+\
                ['L6_IPC_%s'%(i+1) for i in range(5)]+ \
                ['L6_UTPC_%s'%(i+1) for i in range(5)]+ \
                ['L6_BPC_%s'%(i+1) for i in range(5)]+ \
                ['L4_SP_%s'%(i+1) for i in range(5)]+ \
                ['L4_SS_%s'%(i+1) for i in range(5)]+\
                ['L23_PC_%s'%(i+1) for i in range(5)]

# 9 inhibitory types
I_labels =  ['L23_NBC_%s'%(i+1) for i in range(5)]+ \
                ['L23_MC_%s'%(i+1) for i in range(5)]+ \
                ['L23_BTC_%s'%(i+1) for i in range(5)]+ \
                ['L23_DBC_%s'%(i+1) for i in range(5)]+ \
                ['L23_BP_%s'%(i+1) for i in range(5)]+ \
                ['L23_LBC_%s'%(i+1) for i in range(5)]+ \
                ['L23_SBC_%s'%(i+1) for i in range(5)]+ \
                ['L23_ChC_%s'%(i+1) for i in range(5)]+ \
                ['L23_NGC_%s'%(i+1) for i in range(5)]

cell_labels = E_labels + I_labels

# remove undetectable units
remove_these = ['L23_DBC_1','L23_SBC_4']
for label in remove_these:
    cell_labels.remove(label)

pops = cell_labels




num_noisy_spikes = 200
resample = 128
reduce_dur = 5.6
n_levels = 6



#### Import data
print('>>>> Loading data...')
spike_detect_method = 'psta'
this_dir = '%s_snr_pink_noise_adj_extended_%s'%(num_noisy_spikes,spike_detect_method)
print('from',this_dir)

for i, unit in enumerate(cell_labels):

    # model waveforms
    f = os.path.join(data_dir,this_dir,'opt_model_waves_%s.pkl'%unit)

    waves_df = pd.read_pickle(f)

    if i==0:
        all_waves_df = waves_df.copy()
    else:
        join_frames = [all_waves_df,waves_df]
        all_waves_df = pd.concat(join_frames)

    del waves_df

all_waves_df.set_index('unit_id',inplace=True)




# grab units to use
# -- these units are only canonical units
fname = 'opt_case_train_metadata.pkl'
f = os.path.join(data_dir,fname)

metadata_df = pd.read_pickle(f)
metadata_df.set_index('unit_id',inplace=True)




these_unit_ids = np.unique(metadata_df.index.values)
these_waves_df = all_waves_df[all_waves_df.index.isin(these_unit_ids)]




# get only units with a canonical waveform at the center
these_units = np.unique(these_waves_df.index.values)

true_peak = True
center = 15 # this will be the adjusted center channel for 31 channels
spk_tol = 0.21
detect_scale = 1


columns = ['unit_id','neuron_type','wave_type','center_channel']
waveform_class_df = pd.DataFrame(columns=columns)

which_fams = []
which_ntypes = []
which_wtypes = []
which_dists = []
which_centers= []

for unit_id in these_units:

    # grab the unit and info
    unit_df = these_waves_df[these_waves_df.index==unit_id]

    nt = unit_df.neuron_type.tolist()[0]

    # get its simulated EAPs and take center channel
    eaps, orig_center = get_unit_eaps(unit_id,all_waves_df,norm=True,true_peak=true_peak,return_orig_center=True)
    center_eap = eaps[center,:]

    duration = 5.6
    pre_spk = -1.4
    times = np.linspace(pre_spk,pre_spk+duration,len(center_eap))

    # compute the detection bounds for classification
    detect_ub = np.mean(center_eap) + detect_scale*np.std(center_eap)
    detect_lb = np.mean(center_eap) - detect_scale*np.std(center_eap)

    # define spike voltage tolerances for pre and post spike time
    axon_tol = (times<=0)&(times>=-spk_tol)
    pre_spk = center_eap[axon_tol]

    soma_tol = (times>=0)&(times<=2*spk_tol)
    post_spk = center_eap[soma_tol]

    # classify
    if np.any(pre_spk>detect_ub): # -> has prominent positive deflection during pre-spike interval
        if np.any(post_spk<detect_lb):
            wave_class = 'axp-contaminated'
        else:
            wave_class  = 'dipole-inverted'

    elif times[np.argmin(center_eap)]<=-spk_tol: # -> minimum value comes earlier than pre-spike interval
        wave_class = 'axp-dominant'

    elif np.any(post_spk<detect_lb): # -> has prominent negative deflection during pre-spike interval
        wave_class = 'canonical'

    else:
        wave_class = 'other'

    which_centers.append(orig_center)
    which_ntypes.append(nt)
    which_wtypes.append(wave_class)

waveform_class_df['wave_type'] = which_wtypes
waveform_class_df['neuron_type'] = which_ntypes
waveform_class_df['unit_id'] = these_units
waveform_class_df['center_channel'] = which_centers

waveform_class_df.set_index('unit_id',inplace=True)




canonical_df = waveform_class_df[waveform_class_df['wave_type']=='canonical']
canonical_df = canonical_df[canonical_df.index.isin(these_units)]




these_unit_ids = np.unique(these_waves_df.index.values)
these_features_df = get_multilevel_wavelet_features(unit_ids=these_unit_ids,
                                                  df=these_waves_df,
                                                  wavelet='haar',
                                                  level=n_levels,
                                                  resample=resample,
                                                  reduce_dur=reduce_dur)


algorithm = tt.ncp_bcd # tensor decomposition
this_data = a.reformat_data_for_demixing(data_df=these_features_df,
                                       canonical_df=canonical_df, # gives true "center"
                                        algorithm=algorithm)
if test_demixing:
    a.test_demixing_models(data=this_data,algorithm=algorithm,total_iters=20,save_plots=save_plots)

else:
    source_model = algorithm(this_data, rank=R, verbose=False)
    fname = 'R%s_source_model.pkl'%R
    f = os.path.join(path_to_save,fname)
    if saveModel:
        with open(f, 'wb') as file:
            pickle.dump(source_model, file)
