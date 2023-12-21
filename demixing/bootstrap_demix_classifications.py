'''
    
    Main objectives:
        1. Test robustness to random subsamples using bootstrapping.
        2. Save all bootstrap sampled results for visualizations in J.Notebook
'''
import pandas as pd
import numpy as np
import os 
import sys
import analysis as a
# import tensortools as tt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
import time

cwd = os.getcwd()
parent_dir = os.path.dirname(cwd)
sys.path.append(parent_dir)
# print('Added\n --> %s \nto path.' %parent_dir)

data_dir = os.path.join(parent_dir,'model_data','bbp')
path_to_figs = os.path.join(parent_dir,'figures')


# main parameters
R = 4
n_bootstraps = 1000


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

true_peak = True
center = 15 # this will be the adjusted center channel for 31 channels
spk_tol = 0.21
detect_scale = 1

oob_score = True
which_labels = 'metypes'


#### MAIN CODE ####


# grab units to use
# -- these units are only canonical units
fname = 'opt_case_train_metadata.pkl'
f = os.path.join(data_dir,fname)
metadata_df = pd.read_pickle(f)
metadata_df.set_index('unit_id',inplace=True)

# f = 'mtype_grid_results.pkl'
f = 'metype_grid_results.pkl'
grid_results_df = pd.read_pickle(f)
# these_results_df = grid_results_df[grid_results_df.R==R]
best_params = grid_results_df.best_params.iloc[0]



su_prevs_df = pd.read_pickle('r%s_su_demixing_df.pkl'%R)
spatial_prevs_df = pd.read_pickle('r%s_spatial_demixing_df.pkl'%R)
dynam_prevs_df = pd.read_pickle('r%s_dynam_demixing_df.pkl'%R)

# order of spatial patterns
these_keys = ['prop1','prop2','init','bkgd']


# which_labels = neuron_fam_order
which_labels = me_type_order

#### Main results to save #####

# save importance and correct classes in DataFrames
columns = ['sample','acc'] + [str(i) for i in range(R)]
bs_importance_df = pd.DataFrame(columns=columns)

# columns = ['sample'] + neuron_fam_order 
columns = ['sample'] + me_type_order 
bs_correct_class_df = pd.DataFrame(columns=columns)

# list for collecting these
all_accuracy_scores = []
all_importances = []
all_correct_classes = []

# N = len(neuron_fam_order)
N = len(me_type_order)

# save confusion matrices in numpy.array()
bs_confusion_mat = np.zeros((n_bootstraps,N,N)) 

# only for bs subsamples to get same results
random_states = np.random.randint(low=0,high=12345,size=n_bootstraps)

n_samples = 1000+200 # train+test sizes, from sensitivity analysis??

t0 = time.perf_counter()
for n_iter in range(n_bootstraps):
    
    demixed_df = su_prevs_df
    features = demixed_df[[str(i) for i in range(R)]].values
    orig_labels = demixed_df.labels.values
#     label_order = neuron_fam_order
    labels = np.array([get_me_type_from_neuron_fam(nf) for nf in orig_labels])
    label_order = me_type_order
    
    # get random bootstrap subsample w or w/o replacement
    random_sample = np.random.choice(np.arange(len(labels)),size=n_samples,replace=True)
    bs_features = features[random_sample,:]
    bs_labels = labels[random_sample]
    
    
    train_features, test_features, train_labels, test_labels = train_test_split(bs_features,bs_labels,test_size=0.2,
        random_state=random_states[n_iter]) # use the same train-test random_states on rerun of this script
    
    # test parameters:
    # - max_depth
    # - n_estimators
    classifier = RandomForestClassifier(max_depth=best_params['max_depth'],
                                        n_estimators=best_params['n_estimators'],
                                        max_features=None, # this is used to constrain how man features to check, None = total features
                                        oob_score=oob_score,random_state=0)
    
    
    classifier.fit(train_features,train_labels)
    predicted_labels = classifier.predict(test_features)
    
    score = classifier.score(test_features,test_labels)
    all_accuracy_scores.append(score)
    
    importances = classifier.feature_importances_
    all_importances.append(list(importances)) # in order of these_keys
    
    
    C = confusion_matrix(predicted_labels,test_labels,normalize='true',
                     labels=label_order) # label order should be base on distance matrix of centroids
    
    
    bs_confusion_mat[n_iter] = C
    
    
    
# save confusion matrix
fname = 'r%s_metype_bs_confusion_mat.npy'%R
np.save(fname,bs_confusion_mat)
    
# save importances
bs_importance_df['sample'] = [i for i in np.arange(n_bootstraps)]
bs_importance_df['acc'] = all_accuracy_scores

importance_mat = np.array(all_importances)
for i in range(R):
    bs_importance_df[str(i)] = importance_mat[:,i]
    
fname = 'r%s_metype_bs_importance_df.pkl'%R
bs_importance_df.to_pickle(fname,protocol=3)
    

bs_correct_class_df['sample'] = [i for i in np.arange(n_bootstraps)]

# for i, nt in enumerate(neuron_fam_order):
for i, nt in enumerate(me_type_order):
    bs_correct_class_df[nt] = [val for val in bs_confusion_mat[:,i,i]] # grab the diag
    
fname = 'r%s_metype_bs_correct_class_df.pkl'%R
bs_correct_class_df.to_pickle(fname,protocol=3)

                           
tf = time.perf_counter()                    
print('Analysis took: %.2f min'%((tf-t0)/60.))