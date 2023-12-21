import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import tensortools as tt
import pywt
from pywt import wavedec, waverec
from scipy.interpolate import CubicSpline
from scipy.signal import find_peaks


from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, train_test_split

import warnings


# 12 excitatory types
E_LABELS =  ['L5_TTPC1_%s'%(i+1) for i in range(5)]+  \
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
I_LABELS =  ['L23_NBC_%s'%(i+1) for i in range(5)]+ \
                ['L23_MC_%s'%(i+1) for i in range(5)]+ \
                ['L23_BTC_%s'%(i+1) for i in range(5)]+ \
                ['L23_DBC_%s'%(i+1) for i in range(5)]+ \
                ['L23_BP_%s'%(i+1) for i in range(5)]+ \
                ['L23_LBC_%s'%(i+1) for i in range(5)]+ \
                ['L23_SBC_%s'%(i+1) for i in range(5)]+ \
                ['L23_ChC_%s'%(i+1) for i in range(5)]+ \
                ['L23_NGC_%s'%(i+1) for i in range(5)]


def get_unit_eaps(unit_id,all_waves_df,norm=False,num_channels=31,true_peak=True,return_orig_center=False):
    unit_df = all_waves_df[all_waves_df.index==unit_id]

    if true_peak:
        center = find_true_peak(unit_df)
    else:
        center = unit_df.max_channel_i.iloc[0]

    ub=int(center+np.floor(num_channels/2))+1
    lb=int(center-np.floor(num_channels/2))


    eaps = []
    for _, row_df in unit_df.iterrows():
        wave = row_df.waveform
        if norm:
            amp = row_df.amplitude
            wave = np.divide(wave,amp)
        eaps.append(list(wave))

    eaps = np.array(eaps)
    eaps = eaps[lb:ub,:]

    if return_orig_center:
        return eaps, center
    else:
        return eaps


def find_true_peak(unit_df,perc_amp=0.6,spk_tol=0.2,channel_tol=6):
    amps = unit_df.amplitude.values

    amps = np.array(amps)
    amp_threshold = perc_amp*np.max(amps)
    amp_peaks = find_peaks(amps)[0]

    # use middle point for initial center channel
    center_chan = int(len(amps)/2)
    chan_peaks = np.array([peak for peak in amp_peaks if (np.abs(peak-center_chan)<channel_tol)])

    amp_threshold = perc_amp*np.max(amps[chan_peaks])
    viable_peaks = np.array([peak for peak in chan_peaks if (amps[peak]>amp_threshold)])
    if len(viable_peaks)>2:
        viable_peaks = viable_peaks[np.argsort(amps[viable_peaks])[-2:]]


    spike_times = []
    for peak in viable_peaks:

        wave = unit_df.waveform.iloc[peak]
        duration = 5.6
        pre_spk = -1.4
        times = np.linspace(pre_spk,pre_spk+duration,len(wave))

        ve = wave[times>=0][0] # grab ve at intracellular spike time

        # test if  inverted by looking a center of waveform
        if ve>np.mean(ve):
            min_or_max = np.argmax
        else:
            min_or_max = np.argmin

        spkt = times[min_or_max(wave)]
        spike_times.append(spkt)


    spike_times = np.array(spike_times)



    # just grab the closest spike time
    sorted_peaks = np.argsort(np.abs(spike_times))
    spike_times = spike_times[sorted_peaks]
    viable_peaks = viable_peaks[sorted_peaks]
    true_peak = viable_peaks[0] # the earliest

    return true_peak




def get_multilevel_wavelet_features(unit_ids,df,duration=2.7,level=4,wavelet='haar',resample=False,reduce_dur=False):

    columns = ['unit_id','wavelet_features']
    features_df = pd.DataFrame(columns=columns)

    # reduce waves
    waves_df = df[df.index.isin(unit_ids)]


    for unit_id in unit_ids:

        # get one unit
        these_waves_df = waves_df[waves_df.index==unit_id]

        all_waveform_features = []

        for _, row_df in these_waves_df.iterrows():

            wave = row_df['waveform']
            amp = row_df['amplitude']

            T = len(wave)


            # interpolation grid
            interp_t = np.linspace(0,duration,T)

            if reduce_dur:
                wave = wave[interp_t<=reduce_dur]
                interp_t = interp_t[interp_t<=reduce_dur]
                duration = reduce_dur # ignore previous duration for reduced duration

            if resample:

                # downsample and upsample grids
                re_t = np.linspace(0,duration,resample)

                cs = CubicSpline(interp_t,wave,axis=0)
                wave = cs(re_t)


            norm_wave = np.divide(wave,amp)

            level_coeffs = wavedec(norm_wave,wavelet,level=level)
            wave_features, _, _ = pywt.ravel_coeffs(level_coeffs)

            all_waveform_features.append(wave_features.tolist())

        df = pd.DataFrame(columns=columns)

        df['unit_id'] = [unit_id for _ in range(len(all_waveform_features))]
        df['wavelet_features'] = all_waveform_features


        join_frames = [features_df,df]
        features_df = pd.concat(join_frames)


    # reorganize DataFrame
    features_df.set_index('unit_id',inplace=True)
    these_cols = waves_df.columns.tolist()
    these_cols.remove('waveform')


    for col in these_cols:
        features_df[col] = waves_df[col].tolist()

    return features_df


def from_neuron_type_to_fam(nt):
    if len(nt.split('_'))==3:
        nf = nt.split('_')[0]+'_'+nt.split('_')[1]
    else:
        nf = nt.split('_')[0]+'_'+nt.split('_')[1]+'_'+nt.split('_')[2]
    return nf

def from_neuron_type_to_EI(nt):
    if nt in E_LABELS:
        return 'E'
    elif nt in I_LABELS:
        return 'I'
    else:
        raise Exception("Neuron-type %s is not defined as E or I"%nt)


def test_demixing_models(data,algorithm,total_iters=10,show_plots=True,save_plots=False):


    all_sims = []
    all_errs = []

    print('Testing models (total = %s)'%total_iters)
    for R in np.arange(1,total_iters):
        R = int(R)

        W = algorithm(data, rank=R, verbose=False)
        X = algorithm(data, rank=R, verbose=False)
        Y = algorithm(data, rank=R, verbose=False)
        Z = algorithm(data, rank=R, verbose=False)


        # models
        models = [W,X,Y,Z] #[A,B,W,X,Y,Z]
        errors = [W.obj,X.obj,Y.obj,Z.obj] #[A.obj,B.obj,W.obj,X.obj,Y.obj,Z.obj]

        all_errs.append(errors)

        # find lowest recon error
        print('Lowest error for R=%s is %s' %(R,np.min(errors)))
        i = np.argmin(errors)
        U = models[i]

        sims = []

        for j, V in enumerate(models):

            if not j == i:

                sim = tt.kruskal_align(U.factors, V.factors, permute_U=True, permute_V=True)
                sims.append(sim)

        all_sims.append(sims)


    if show_plots:
        fig = plt.figure(figsize=(6,4))
        mean_errs = []

        for i, errs in enumerate(all_errs):
            plt.scatter([i,i,i,i],errs,color='k',s=80,alpha=0.6)

            mean_errs.append(np.mean(errs))

        plt.plot(mean_errs,'r',zorder=0)


        plt.axhspan(0.1,0.5,color='g',alpha=0.3,zorder=0)
        plt.ylim(bottom=np.min(mean_errs)-0.05)

        plt.xticks(range(total_iters-1),[i+1 for i in range(total_iters-1)])
        plt.xlabel('# components')
        plt.ylabel('error');

        if save_plot:
            fname = 'demix_test_recon_err.pdf'
            fig.savefig(fname,bbox_inches='tight',format='pdf')
            fname = 'demix_test_recon_err.png'
            fig.savefig(fname,bbox_inches='tight',format='png')

        fig, ax = plt.subplots(figsize=(6,6))

        means = []
        varis = []

        for i, sims in enumerate(all_sims):
            ax.scatter([i,i,i],sims,color='k',s=80,alpha=0.6)

            means.append(np.mean(sims))
            varis.append(np.var(sims))

        ax.plot(means,'r',zorder=0,linewidth=2)

        ax.axhline(0.98,color='gray',linestyle='--',zorder=0)

        ax.set_xticks(range(total_iters-1))
        ax.set_xticklabels([i+1 for i in range(total_iters-1)])
        ax.set_xlabel('# components')
        ax.set_ylabel('similarity (mean)')

        if save_plot:
            fname = 'demix_test_similarity.pdf'
            fig.savefig(fname,bbox_inches='tight',format='pdf')
            fname = 'demix_test_similarity.png'
            fig.savefig(fname,bbox_inches='tight',format='png')




def get_max_amps(data_df):

    unit_ids = np.unique(data_df.index.values)

    max_amps_i = []
    max_amps = []

    for unit_id in unit_ids:
        df = data_df[data_df.index==unit_id]
        amplitudes = df.amplitude.values

        max_amps_i.append(np.argmax(amplitudes))
        max_amps.append(np.max(amplitudes))

    return max_amps, max_amps_i



def reformat_data_for_demixing(data_df,canonical_df=None,tensor_feature='wavelet_features',algorithm=tt.ncp_bcd,force_positive=False,D=31):


    unit_ids = np.unique(data_df.index.values)

    # data structure dimensions
    T = len(data_df[tensor_feature].iloc[0])
    N = len(unit_ids)

    if canonical_df is None:
        max_amps, max_amps_i = get_max_amps(data_df)
        these_chans_i = max_amps_i
    else:
        which_peak_i = []
        for unit_id in unit_ids:
            unit_df = canonical_df[canonical_df.index==unit_id]
            true_peak = unit_df.center_channel.iloc[0]
            which_peak_i.append(true_peak)

        these_chans_i = which_peak_i





    data = np.zeros((N,D,T))

    print('>>>> Organizing data')
    # fix probe alignment and grab only the inner channels
    for i, exid in enumerate(unit_ids):
        df = data_df[data_df.index==exid]
        max_chan = len(df)

        # grab only the inner channels
        center_chan_i = these_chans_i[i]
        lb, ub = int(center_chan_i-np.floor(D/2)), int(center_chan_i+np.floor(D/2)+1)

        # bad units, but will include regardless
        if ub>=max_chan:
            print('Max channel (%s) too high for'%ub,df.neuron_type.iloc[0])
            ub = max_chan-1
            lb = ub-D


        these_channels = np.arange(lb,ub)
        for j, idx in enumerate(these_channels):
            these_features = df[tensor_feature].iloc[idx]

            # non-negative CP tensor decompositions need data>=0
            if algorithm in [tt.ncp_bcd,tt.ncp_hals]:
                these_features = np.abs(these_features)


            data[i,j,:] = these_features

    if force_positive:
        min_ = np.abs(np.min(data)) # should be non-positive, if 0, will not do anything
        data = np.add(data,min_) # add the absolute value of mininum, performed element-wise



    return data


def get_constrained_neuron_type_reps(model):

    unit_factors = model.factors[0].T
    _, N = np.shape(unit_factors)

    I_bar = np.ones(N)
    N_bar = unit_factors
    N_inv = np.linalg.pinv(N_bar)

    R_bar = np.dot(I_bar,N_inv)

    error = np.dot(R_bar,N_bar)-I_bar
    rel_error = np.linalg.norm(error/I_bar)
    print('Relative Error :',rel_error)

    constrained_units = np.transpose(N_bar)*R_bar

    # R_bar is needed for other factors
    return constrained_units, R_bar


def compute_centroids(data,labels,test_stat=np.median):

    N,R = np.shape(data)
    centroids = np.zeros((len(np.unique(labels)),R))

    # get centroids for labeled groups
    for l in np.unique(labels):
        grouped_data = data[labels==l]
        for r in range(R):
            x = grouped_data[:,r]
            center_x = test_stat(x)

            centroids[l,r] = center_x

    return centroids


def bootstrap_centroids(data,group_labels,sample_frac,test_stat=np.median,n_bootstraps=1000,with_replacement=True,within_group=True):

    N,R = np.shape(data)
    L = len(np.unique(group_labels))
    total_labels = len(group_labels)

    labels = np.zeros(total_labels,dtype='int64')

    for i,l in enumerate(np.unique(group_labels)):
        labels[np.argwhere(group_labels==l)] = i


    all_centroids = np.zeros((n_bootstraps,L,R))

    if not within_group:
        n_samples = int(sample_frac*N)
        print('>>>> Bootstrapping n=%s samples x %s'%(n_samples,n_bootstraps))
    else:
        print('>>>> Bootstrapping x %s within groups'%n_bootstraps)

    for i in range(n_bootstraps):

        if within_group:


            # gather all the ids
            for l in np.unique(labels):

                group_ids = np.where(labels==l)[0]

                G = len(group_ids)
                g_samples = int(sample_frac*G)

                if l==0:
                    sample_ids = np.random.choice(group_ids,size=g_samples,
                                              replace=with_replacement)
                else:
                    temp_ids = np.random.choice(group_ids,size=g_samples,
                                              replace=with_replacement)

                    sample_ids = np.concatenate([sample_ids,temp_ids])



        else:
            sample_ids = np.random.choice(np.arange(N),size=n_samples,
                                          replace=with_replacement)


        sample_data = data[sample_ids,:]
        sample_labels = labels[sample_ids]

#         sample_size, _ = np.shape(sample_data)
#         if (sample_size!=total_samples):
#             raise Exception("Sample size doesn't match desired size")


        centroids = compute_centroids(sample_data,sample_labels,test_stat=test_stat)
        all_centroids[i,:,:] = centroids


    return all_centroids

def grid_search_RFC_for_demixing(demixed_df,plot_results=False,save_plot=False,oob_score=False,max_features=False,figname=None,verbose=False):
    '''
        Used to find hyperparameters for neuron-type identification limits.
    '''

    # warnings associated with lack of Out-of-Bag estimates in bootstrap aggregation
    # - UserWarning is due to no OOB samples in which case the sum of their predictions is 0
    # - RuntimeWarning is due to prediction=0 as denominator for OOB score
    warnings.filterwarnings("ignore",category=UserWarning)
    warnings.filterwarnings("ignore",category=RuntimeWarning)

    if type(demixed_df) is list:
        features = demixed_df[0]
        labels = demixed_df[1]

    else:
        comp_ids = demixed_df.columns.tolist()[:-1] # except for label column

        # features should have been ordered for dataframe for consistency
        features = demixed_df[comp_ids].values
        labels = demixed_df.labels.values


    # for tuning hyperparameters
    depth_range = np.arange(2,22,2)
    est_range = np.arange(5,100,5)

    total_depths = len(depth_range)
    total_ests = len(est_range)

    N_samples, N_features = np.shape(features)
    if not max_features:
        if N_features>10:
            max_features = 'auto'
        else:
            max_features = None



    # tuning hyperparameters
    params_to_tune = {'max_depth' : np.arange(2,22,2),
                 'n_estimators' : np.arange(5,100,5)}

    def score_accuracy(classifier,features,labels):
        return classifier.score(features,labels)

    if verbose:
        print('Performing grid search for:',list(params_to_tune.keys()))

    grid_classifier = GridSearchCV(RandomForestClassifier(max_features=max_features,
                                                          oob_score=oob_score, # lower n_estimators cause an issue with this
                                                          random_state=0),
                                    params_to_tune,
                                    scoring=score_accuracy)

    grid_classifier.fit(features,labels) # based on full dataset

    if verbose:
        print('Results found')
    # Final results of tuning RFC
    grid_results = grid_classifier.cv_results_


    if plot_results:
        print('Plotting results...')
        fig = plt.figure(figsize=(8,5))
        means = grid_results['mean_test_score']
        stds = grid_results['std_test_score']

        plot_means = np.reshape(means,(total_depths,total_ests))
        plot_stds = np.reshape(stds,(total_depths,total_ests))

        colors = sns.color_palette('Oranges',total_depths)
        for i, (plot_mean, plot_std) in enumerate(zip(plot_means,plot_stds)):
            plt.errorbar(est_range,plot_mean,
                    yerr=plot_std,
                     color=colors[i],
                     label='max_depth=%s'%depth_range[i])

        plt.xticks(np.arange(0,110,10))
        plt.ylabel('Accurancy Score')
        plt.xlabel('n_estimators')
        plt.title('Grid Search Results for Hyperparameters')
        plt.legend(loc='lower right');

        if save_plot:
            if figname is None:
                figname = 'grid_search'
            fig.savefig(figname+'.png',format='png',bbox_inches='tight')
            fig.savefig(figname+'.pdf',format='pdf',bbox_inches='tight')

    # return this Dict object
    best_params = grid_classifier.best_params_
    return best_params
