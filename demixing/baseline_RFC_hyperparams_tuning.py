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


max_features = 4 # used to compare against R=4 sources

for test_mtypes in ['EI','ME','M']:

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
    these_unit_ids = metadata_df.unit_id.values


    features_df = pd.read_csv('final_classic_features.csv',index_col=0)
    # feature_columns = features_df.columns.tolist()[1:]

    # for fair comparison... don't use asymmetric features
    total_prop_vel = features_df['prop_vel_above'].values + features_df['prop_vel_below'].values
    features_df['total_prop_vel'] = total_prop_vel
    feature_columns = ['duration','half_width','pt_ratio','repol_slope','recov_slope','spread','total_prop_vel']

    # constrain to the same units
    these_features_df = features_df[features_df.unit_id.isin(these_unit_ids)].copy()
    constrained_unit_ids = these_features_df.unit_id.values
    these_features_df['neuron_type'] = metadata_df[metadata_df.unit_id.isin(constrained_unit_ids)].neuron_type.values


    # get labels
    neuron_fams = [from_neuron_type_to_fam(nt) for nt in these_features_df.neuron_type.values]
    ei_fams = [from_neuron_type_to_EI(nt) for nt in these_features_df.neuron_type.values]
    me_types = [get_me_type_from_neuron_fam(nf) for nf in neuron_fams]

    all_accuracy_scores = []
    all_importances = []
    all_best_params = []


    t0_total = time.perf_counter()
    print('>>>> Tuning Baseline RF Classifier for Test: %s'%test_mtypes)

    features = these_features_df[feature_columns].values

    # data labels
    if test_mtypes == 'M':
        labels = neuron_fams
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
        oob_score = False # I get errors for this

    t0 = time.perf_counter()

    # split data before hyperparam tuning to prevent test leakage
    train_features, test_features, train_labels, test_labels = train_test_split(
    features,labels,test_size=0.2,
    random_state=12345)


    best_params = grid_search_RFC_for_demixing([train_features,train_labels],
                                               plot_results=True,
                                               save_plot=True,
                                               oob_score=oob_score,
                                               max_features=max_features,
                                               figname='baseline_%s_grid_search'%test,
                                               verbose=True)

    all_best_params.append(best_params)
    print(best_params)

    # how long does this take?
    tf = time.perf_counter()
    print('>>>> Grid Search Time (min) = ',(tf - t0)/60.)


    # test parameters:
    # - max_depth
    # - n_estimators
    classifier = RandomForestClassifier(max_depth=best_params['max_depth'],
                                        n_estimators=best_params['n_estimators'],
                                        max_features=max_features, # this is used to constrain how man features to check, None = total features
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


    # plt.colorbar(cax)
    # if test_mtypes:
    #     test = 'mtype'
    # else:
    #     test = 'ei'
    fig.savefig('baseline_max%s_%s_confusion_matrix.png'%(max_features,test),format='png',bbox_inches='tight')
    fig.savefig('baseline_max%s_%s_confusion_matrix.pdf'%(max_features,test),format='pdf',bbox_inches='tight')



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
    plt.axhline(0.5,linestyle='--',color='k',zorder=1,linewidth=0.8)

    plt.xlim([-0.5,len(label_order)-0.5])
    plt.ylim(bottom=0)
    plt.ylabel('Accuracy (fraction correct)',fontsize=16)


    # if test_mtypes:
    #     test = 'mtype'
    # else:
    #     test = 'ei'
    fig.savefig('baseline_max%s_%s_true_classes.png'%(max_features,test),format='png',bbox_inches='tight')
    fig.savefig('baseline_max%s_%s_true_classes.pdf'%(max_features,test),format='pdf',bbox_inches='tight')


    # important features
    importances = classifier.feature_importances_
    all_importances.append(importances)

    fig = plt.figure()
    plt.bar(np.arange(len(feature_columns)),importances,color='grey')
    plt.axhline(1./len(feature_columns),color='k',linestyle='--',linewidth=0.8) # equal importance line
    plt.xticks(np.arange(len(feature_columns)),feature_columns,rotation=45,fontsize='small')

    # if test_mtypes:
    #     test = 'mtype'
    # else:
    #     test = 'ei'
    fig.savefig('baseline_max%s_%s_importance.png'%(max_features,test),format='png',bbox_inches='tight')
    fig.savefig('baseline_max%s_%s_importance.pdf'%(max_features,test),format='pdf',bbox_inches='tight')



    grid_results_df = pd.DataFrame(columns=['acc','importance','best_params'])
    grid_results_df['acc'] = all_accuracy_scores
    grid_results_df['importance'] = all_importances
    grid_results_df['best_params'] = all_best_params

    if test_mtypes == 'M':
        grid_results_df.to_pickle('baseline_max%s_mtype_grid_results.pkl'%max_features,protocol=3)
    elif test_mtypes == 'ME':
        grid_results_df.to_pickle('baseline_max%s_metype_grid_results.pkl'%max_features,protocol=3)
    else:
        grid_results_df.to_pickle('baseline_max%s_ei_grid_results.pkl'%max_features,protocol=3)

    tf_total = time.perf_counter()

    print('>> Total Performance Time =',(tf_total-t0_total)/60.)
