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
import time

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
# print('Added\n --> %s \nto path.' %parent_dir)

data_dir = os.path.join(parent_dir,'model_data','bbp')
path_to_figs = os.path.join(parent_dir,'figures')
from sklearn.metrics import confusion_matrix


test_mtypes = 'ME'  # 'EI' vs 'M' vs 'ME'

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


# different label groupings
ei_fam_order = ['E','I']
neuron_fam_order = [
    'L5_TTPC2', # 1
    'L5_TTPC1', # 2
    'L23_PC',   # 3
    'L5_STPC',  # 4
    'L4_SS',    # 5
    'L6_BPC',    # 6
    'L6_TPC_L1',# 7
    'L4_SP',    # 8
    'L5_UTPC',  # 9
    'L6_UTPC',  # 10
    'L6_TPC_L4',# 11
    'L6_IPC',   # 12
    'L23_DBC',  # 13
    'L23_LBC',  # 14
    'L23_NBC',  # 15
    'L23_MC',   # 16
    'L23_BTC',  # 17
    'L23_NGC',  # 18
    'L23_ChC',  # 19
    'L23_SBC',  # 20
    'L23_BP'    # 21
]
me_type_order = ['E1','E2','E3','E4','I1','I2']


def get_me_type_from_neuron_fam(nf):
    if nf in ['L5_TTPC2','L5_TTPC1']:
        return 'E1'
    if nf in ['L4_SP','L5_UTPC']:
        return 'E2'
    if nf in ['L5_STPC','L23_PC','L6_BPC','L4_SS','L6_TPC_L1']:
        return 'E3'
    if nf in ['L6_UTPC','L6_IPC','L6_TPC_L4']:
        return 'E4'
    if nf in ['L23_DBC','L23_LBC','L23_NBC','L23_MC', 'L23_BTC','L23_NGC']:
        return 'I1'
    if nf in ['L23_BP','L23_ChC','L23_SBC']:
        return 'I2'


num_noisy_spikes = 200
resample = 128
reduce_dur = 5.6
n_levels = 6



# grab units to use
# -- these units are only canonical units
fname = 'opt_case_train_metadata.pkl'
f = os.path.join(data_dir,fname)

metadata_df = pd.read_pickle(f)
metadata_df.set_index('unit_id',inplace=True)


#### Import data
print('>> Loading data...')
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

#### !--- this should be moved to a function at some point
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


# get labels
neuron_fams = [from_neuron_type_to_fam(nt) for nt in canonical_df.neuron_type.values]
ei_fams = [from_neuron_type_to_EI(nt) for nt in canonical_df.neuron_type.values]
me_types = [get_me_type_from_neuron_fam(nf) for nf in neuron_fams]


comp_keys = [['S1'],
             ['S1','S2'],
             ['S2','S4','S1'],
             ['S1','S3','S4','S2']]

all_accuracy_scores = []
all_importances = []
all_best_params = []


t0_total = time.perf_counter()

# R = 4
test_mtypes = 'M'

for R in range(1,5):
# for test_mtypes in ['EI','ME','M']:

    print('>>>> Tuning R=%s RF Classifier'%R)

    su_prevs_df = pd.read_pickle('r%s_su_demixing_df.pkl'%R)
    spatial_prevs_df = pd.read_pickle('r%s_spatial_demixing_df.pkl'%R)


    these_keys = comp_keys[R-1]

    demixed_df = su_prevs_df
    features = demixed_df[[str(i) for i in range(R)]].values

    # data labels
    if test_mtypes == 'M':
        labels = demixed_df.labels.values
        label_order = neuron_fam_order
    elif test_mtypes == 'ME':
        labels = me_types
        label_order = me_type_order
    else:
        labels = ei_fams
        label_order = ei_fam_order

    # for saving data
    if test_mtypes == 'M':
        test = 'mtype'
        oob_score = True
    elif test_mtypes == 'ME':
        test = 'metype'
        oob_score = True
    else:
        test = 'ei'
        oob_score = False

    t0 = time.perf_counter()

    train_features, test_features, train_labels, test_labels = train_test_split(
    features,labels,test_size=0.2,
    random_state=12345)

    best_params = grid_search_RFC_for_demixing([train_features,train_labels],
                                               plot_results=True,
                                               save_plot=True,
                                               oob_score=oob_score,
                                               figname='r%s_%s_grid_search'%(R,test),
                                               verbose=True)

    all_best_params.append(best_params)
    print(best_params)

    # how long does this take?
    tf = time.perf_counter()
    print('>>>> Grid Search Time (s) = ',(tf - t0)/60.)


    # test parameters:
    # - max_depth
    # - n_estimators
    classifier = RandomForestClassifier(max_depth=best_params['max_depth'],
                                        n_estimators=best_params['n_estimators'],
                                        max_features=None, # this is used to constrain how man features to check, None = total features
                                        oob_score=oob_score,random_state=0)


    classifier.fit(train_features,train_labels)
    predicted_labels = classifier.predict(test_features)

    # mean accuracy (some may be much better than others)
    score = classifier.score(test_features,test_labels)
    all_accuracy_scores.append(score)

    C = confusion_matrix(predicted_labels,test_labels,normalize='true',
                     labels=label_order) # label order should be base on distance matrix of centroids

    fig = plt.figure(figsize=(9,8.5))
    ax = plt.subplot(111)
    cax = ax.imshow(C,interpolation='none',origin='upper',vmin=0,vmax=1)

    # print(100*C)

    avg_acc = 100*np.mean(np.diag(C))
    #     print(avg_acc)
    # title_str +='\nClustering Accuracy: %.2f Percent'%avg_acc
    ax.set_xlabel('Predicted Class',fontsize=16)
    ax.set_ylabel('True Class',fontsize=16)
    ax.set_xticks(range(len(label_order)))
    ax.set_xticklabels(label_order,rotation=90)
    ax.set_yticks(range(len(label_order)))
    ax.set_yticklabels(label_order)#,fontsize=16)
    plt.colorbar(cax)

    fig.savefig('final_r%s_%s_confusion_matrix.png'%(R,test),format='png',bbox_inches='tight')
    fig.savefig('final_r%s_%s_confusion_matrix.pdf'%(R,test),format='pdf',bbox_inches='tight')



    # neuron types correctly predicted above chance
    correct_predictions = np.diag(C)
    ordered_predictions = np.sort(correct_predictions)
    ordered_fams = np.array(label_order)[np.argsort(correct_predictions)]

    fig = plt.figure()
    plt.bar(np.arange(len(ordered_predictions)),ordered_predictions,color='grey')
    plt.xticks(np.arange(len(label_order)),
               ordered_fams,
               rotation=90)

    # more than half correctly identified
    plt.axhline(0.5,linestyle='--',color='k',linewidth=0.8,zorder=1)

    plt.xlim([-0.5,len(label_order)-0.5])
    plt.ylim(bottom=0)
    plt.ylabel('Accuracy (fraction correct)',fontsize=16)

    fig.savefig('final_r%s_%s_true_classes.png'%(R,test),format='png',bbox_inches='tight')
    fig.savefig('final_r%s_%s_true_classes.pdf'%(R,test),format='pdf',bbox_inches='tight')


    # important features
    importances = classifier.feature_importances_
    all_importances.append(importances)

    fig = plt.figure()
    plt.bar(np.arange(R),importances,color='grey')
    plt.axhline(1./R,color='k',linestyle='--',linewidth=0.8)
    plt.xticks(np.arange(R),comp_keys[R-1])

    fig.savefig('final_r%s_%s_importance.png'%(R,test),format='png',bbox_inches='tight')
    fig.savefig('final_r%s_%s_importance.pdf'%(R,test),format='pdf',bbox_inches='tight')



grid_results_df = pd.DataFrame(columns=['R','acc','importance','best_params'])
grid_results_df['acc'] = all_accuracy_scores
grid_results_df['importance'] = all_importances
grid_results_df['best_params'] = all_best_params

if test_mtypes == 'M':
    grid_results_df.to_pickle('final_var_r_mtype_grid_results.pkl',protocol=3)
elif test_mtypes == 'ME':
    grid_results_df.to_pickle('final_metype_grid_results.pkl',protocol=3)
else:
    grid_results_df.to_pickle('final_ei_grid_results.pkl',protocol=3)

tf_total = time.perf_counter()

print('>> Total Performance Time =',(tf_total-t0_total)/60.)
