import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.patches import Patch

import os
from shutil import copyfile


from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split
from sklearn.metrics import precision_recall_fscore_support
from sklearn.svm import SVC
from sklearn.metrics.pairwise import laplacian_kernel
from sklearn.mixture import BayesianGaussianMixture as bgmm

from sklearn.preprocessing import OneHotEncoder
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import dendrogram
from scipy.spatial.distance import pdist
    
from  MUA_utility import fetch, slice_data, compute_si
import MUA_constants as const

def onset_offset_response(plots_dest_dir_appdx='', csv_dest_dir_appdx='', 
                          generate_plots=True, single_channels=True, 
                          draw_gmm_fit=True):
    """This function presents the start of the onset-offset analysis pipeline.
    This functions iterates all mice, all paradigms, all stimulus types as usual
    then fetches the dataframe with all the negative time stamps. Out of this,
    a histogram of the 0-20ms trace is generated (200 bins, 0.1ms timeframe 
    bins). The bin counts serve as the input feature for the classifier. The 
    final CSV with the computed bin counts will be saved at 
    P["outputPath"]/csv_dest_dir_appdx/onset_offset_spikebins_*_.csv'. If 
    `csv_dest_dir_appdx` is None, no data is saved. In addition, 
    `generate_plots` gives the option to draw the histograms.`draw_gmm_fit` will
    fit a Gaussian Mixture Model to the histogram with 10 components. When 
    `single_channels` is True (default), the 32 channels are not collpassed into
    the previously defined region mapping. The classifier is designed to 
    classify single channels, not collapsed regions. Setting this to False is 
    implemented for plotting region histograms. The plots are saved at 
    P["outputPath"]/plots_dest_dir_appdx/onset_offset_*_*_.png'. The produced
    histogram data is also returned besides being being saved."""

    # get all the available data from the output dir
    if single_channels:
        data = fetch()
        which_region = 'channels'
        plt_spacers = {'hspace':0, 'right':.97, 'top':.96, 'left':.1, 'bottom':.07}
    else:
        data = fetch(collapse_ctx_chnls=True, collapse_th_chnls=True, 
                     drop_not_assigned_chnls=True)
        which_region = 'regions'
        plt_spacers = {'hspace':0, 'right':.97, 'top':.85, 'left':.1, 'bottom':.22}

    # due to different base level actitvity, set heatmap vmax seperately 
    vmaxs =  {'mGE82': 20 ,
              'mGE83': 50 ,
              'mGE84': 50,
              'mGE85': 30,}

    # iter the usual dimensions
    all_spike_bins = []
    for m_id in const.ALL_MICE:
        for parad in const.ALL_PARADIGMS:
            for stim_t in const.ALL_STIMTYPES:
                key = '-'.join([m_id, parad, stim_t])
                # not all combinations of paradigm/stimtype exist
                if key not in data.keys():
                    continue

                # get the negative sptike time stamps (a MultiIndex DataFrame) 
                spikes = slice_data(data, m_id, parad, stim_t, neg_spikes=True, 
                                    drop_labels=True)[0].sort_index(axis=1)
                if not single_channels:
                    spikes = spikes.reindex(reversed(const.REGIONS.keys()), 
                                            axis=1, level=0)

                if generate_plots:
                    # init plot with fitting size (height depends on n channels/ regions)
                    nregions = len(spikes.columns.unique(0))
                    fig, axes = plt.subplots(nrows=nregions, figsize=(6,.4*nregions))
                    fig.subplots_adjust(**plt_spacers)
                    [ax.tick_params(bottom=False, left=False, labelbottom=False, 
                    labelleft=False) for ax in axes.flatten()]

                    axes[0].set_title(key)
                else:
                    # dummy
                    axes = range(len(spikes.columns.unique(0)))
                
                for region, ax in zip(spikes.columns.unique(0), axes):
                    # get the spike timestamp data sliced to the channel/ region
                    region_spikes = spikes[region].values.flatten()
                    # set 0 to nan _> ignored by np.hist
                    region_spikes[region_spikes == 0] = np.nan
                    hist = np.histogram(region_spikes, bins=2000, range=(-50,200))

                    # slice to 0-20ms bins
                    start, stop = 400, 600  # relative to 2000 bins from -50-200
                    spike_bins = hist[0][np.newaxis, start:stop]
                    
                    if spikes.shape[0] != 200:
                        # norm spike counts to 200 trails
                        spike_bins = spike_bins/ spikes.shape[0]
                        spike_bins = (spike_bins*200).astype(int)
                    
                    lbl = key+f'-{region:0>2d}' if single_channels else key+f'-{region}'
                    all_spike_bins.append(pd.Series(spike_bins[0], name=lbl))

                    if not generate_plots:
                        continue
                    
                    # draw the heatmap, setup axis
                    ax.imshow(spike_bins, aspect='auto', extent=(hist[1][start], 
                              hist[1][stop], 0, 1), vmin=0, vmax=vmaxs[m_id])
                    ax.set_ylabel(region, rotation=0, va='center', labelpad=20)
                    ax.set_ylim(0, .4)

                    xt = [0,2,4,6,8,10,12,14,16,18,20]
                    ax.set_xticks(xt)
                    ax.set_xlim(xt[0], xt[-1])
                    
                    if draw_gmm_fit and (spike_bins>1).sum() > 15:
                        model = bgmm(10, 
                                     covariance_type='diag', 
                                     random_state=1, 
                                     mean_prior=(8,), 
                                     covariance_prior=(.1,),
                                     degrees_of_freedom_prior=10)
                        
                        post_stim = np.logical_and(region_spikes>0, region_spikes<20)
                        model.fit(region_spikes[post_stim, np.newaxis])

                        x = np.linspace(-50, 200, 2000).reshape(2000,1)
                        logprob = model.score_samples(x)
                        pdf = np.exp(logprob)
                        ax.plot(x, pdf, '-w', linewidth=.8)

                    # the last plot gets a red stimulus indication bar
                    if region == spikes.columns.unique(0)[-1]:
                        ax.tick_params(bottom=True, labelbottom=True)
                        ax.set_xlabel('[ms]')
                        ax.hlines(0, 0, 8, clip_on=False, linewidth=6, color='r')
                        ax.annotate('Stimulus', (2.3,-.6), color='r', 
                                    fontsize=15, annotation_clip=False)
                
                if generate_plots:
                    f = (f'{const.P["outputPath"]}/{plots_dest_dir_appdx}/'
                         f'onset_offset_{key}_{which_region}.png')
                    fig.savefig(f)
                    plt.close()

    on_off_scores = pd.concat(all_spike_bins, axis=1).T
    if csv_dest_dir_appdx is not None:
        f = f'{const.P["outputPath"]}/{csv_dest_dir_appdx}/onset_offset_spikebins_{which_region}.csv'
        on_off_scores.to_csv(f)
    return on_off_scores

def onset_offset_labels(dest_dir_appdx):
    """Second function in the onset-offset pipeline. This script constructs
    an empty labels-tsv file with all combinations of 
    mouse-paradigm-stimulustype-channel in the index (examples). The column
    `label` describes the label (=0), the column `file` has the full path to the 
    raster plot corresbonding to that example. Also, a .txt file is generated 
    that contains gwenview commands (image display program on linux) that open
    a batch of 32 rasters (mouse-paradigm-stimulustype). The idea is to execute 
    one such command, go through the 32 raster plots, then change labels from 0
    to 1, 2, or 3. if a positive was found. I labeled perfect examples with 1,
    less clear examples with 2, and on-off example with shorter intervals with 
    3. The empty labels tsv is saved at 
    P["outputPath]/dest_dir_appdx/onset_offset_labels_empty.tsv, the commands in
    P["outputPath]/dest_dir_appdx/gwenview_command.txt. Note that the next 
    function reads in onset_offset_labels.tsv instead of 
    onset_offset_labels_empty.tsv
    """
    rasterfiles = []
    gwenview_cmd = ''
    # iter data, built gwenview command and rasterfiles dataframe + label=0
    for parad in const.ALL_PARADIGMS:
        for m_id in const.ALL_MICE:
            for stim_t in const.PARADIGMS_STIMTYPES[parad]:
                gwenview_cmd += f'{m_id}-{parad}-{stim_t}\ngwenview '
                for channel in range(1,33):
                    key = f'{m_id}-{parad}-{stim_t}-{channel:0>2d}'
                    
                    rasterf = (f'{const.P["outputPath"]}/{const.MICE_DATES[m_id]}'
                               f'_{parad}.mcd/Raster_NegativeSpikes_Triggers_'
                               f'{stim_t}_ElectrodeChannel_{channel:0>2d}.png')
                    rasterfiles.append(pd.Series([0, rasterf], ['label', 'file'], 
                                                 name=key))

                    gwenview_cmd += f'{rasterf} '
                gwenview_cmd += '\n\n'
    rasterfiles = pd.concat(rasterfiles, axis=1).T

    labels_tsv = f'{const.P["outputPath"]}/{dest_dir_appdx}/onset_offset_labels_empty.tsv'
    rasterfiles.to_csv(labels_tsv, sep='\t')

    openimages_txt = f'{const.P["outputPath"]}/{dest_dir_appdx}/gwenview_command.txt'
    with open(openimages_txt, 'w') as file:
        file.write(gwenview_cmd)
    
    print(labels_tsv)
    print(openimages_txt)


def lapl_kernel_SVM(dest_dir_appdx='', training_data_dir='', parameter_search=False, 
                    plot_cv_result=False, analyize_confusions=False, pred_X=None):
    """The third function in the onset-offset pipeline. This function combines
    the definition, CV-optimization, result analysis, and prediction on new data.
    The model being used a laplace-kernel SVM trained on weighted examples.
    In all use cases, `training_data_dir` must be passed. This dir should contain
    onset_offset_labels.tsv and onset_offset_spikebins_channels.csv generated in
    the previous 2 functions. The 4 sub-functions mentioned above are all set to
    False by default. The natural order would be to 1, pass `parameter_search` 
    as True to do CV on a set of hyperparameters hardcoded in the function. This
    will produce 2 output files with the CV results (1 trained on weighted data,
    one on unweighted) saved at 
    P["outputPath"]/dest_dir_appdx/cv_laplace_kernel_SVM_*.csv. Second, if 
    `plot_cv_result` is True, the file CV result file mentinoed above will be
    read in and plotted. 3 specific handpicked parametersets are annotated (this
    needs adjustment for new data of course). The plot is saved at 
    P["outputPath"]/dest_dir_appdx/laplce_kernel_SVM_cv.png. These two usees of 
    this function use the `dest_dir_appdx` argument, the next two dont't because
    no files are produced. When passing analyze_confusions as True, the 
    classifier (explicit optimal hyperparamters in code) is evaluated at a 
    random test/train split and the FP and FN classifications (raster plots) 
    are printed out (the command to open them). This should help to get a feel 
    on the features the classifier (doesn't) consider. The last option, `pred_X`
    is passed when we want to use the model to predict labels on new data. 
    pred_x is a pd.DataFrame with examples as rows and features (200 bin counts)
    as columns. The predicted labels are returned.
    """

    def predict(X_train, X_test, y_train, y_test, train_sweights, 
                weight_samples, gamma, C, dec_b):

        # standardize the bin counts
        scaler = StandardScaler().fit(X_train)
        X_train = scaler.transform(X_train)
        X_test = scaler.transform(X_test)

        # define the model
        comp_gram = lambda X, X_: laplacian_kernel(X, X_, gamma/X.shape[1])
        classifier = SVC(class_weight = 'balanced', 
                         random_state = 1, 
                         kernel = comp_gram, 
                         probability = True,
                         C = C)
        
        if not weight_samples:
            classifier.fit(X_train, y_train)
        else:
            classifier.fit(X_train, y_train, train_sweights)

        y_pred_p = classifier.predict_proba(X_test)
        # convert to binary prediction using the specified desicion boundry
        y_pred = np.where(y_pred_p[:,1] <dec_b, 0, 1)
        # for known label compute score, else return prediction
        if y_test is not None:
            pres, recall, f1, _ = precision_recall_fscore_support(y_test, y_pred)
            return pres[1], recall[1], f1[1]
        # this is entered when pred_X is passed
        else:
            return y_pred_p, y_pred

    def cv(weighted, gamma, C, dec_b):
        kf = KFold(6, shuffle=True, random_state=1)
        
        scores = []
        for i, (train_idx, test_idx) in enumerate(kf.split(X.index)):
            # slice X, y and the example weight vector to the fold index
            X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
            train_sweights = sample_weights[train_idx]
            y_train, y_test = y[train_idx], y[test_idx]

            # run the classifier
            score = predict(X_train, X_test, y_train, y_test, train_sweights,
                            weighted, gamma, C, dec_b)

            # save the scores
            scores.append(pd.Series(score, name=i,
                          index=['presicion', 'recall', 'f1']))
        scores = pd.concat(scores, axis=1).T
        return scores.mean()
    
    def find_hyperparamters(weighted):
        # the hyperparamters to test (all combinations)
        gammas = [.3, .4, .5, .7, .8, .9, 1, 1.1, 1.2, 1.5]
        Cs = [.3,.5,.7,.8, 1, 1.2, 1.75, 2, 2.5 ]
        boundries = [.3, .25, .2, .15, .1,.08,.05]
        
        i = 0
        scores = []
        n_params = len(gammas)*len(Cs)*len(boundries)
        for gamma in gammas:
            for C in Cs:
                for dec_b in boundries:
                    # log progress
                    print(f'{i+1}/{n_params}', end='  ')
                    
                    # save recall, presicion and f1 score for parameter set
                    # parameter values are indicated in the label
                    params = f'weighted:{weighted}-gamma:{gamma:.2f}-C:{C:.2f}-bound:{dec_b:.2f}'
                    score = cv(weighted, gamma, C, dec_b).rename(params)
                    print(score.name)
                    
                    scores.append(score)
                    i += 1
        
        scores = pd.concat(scores, axis=1).T
        weighted = 'weighted' if weighted else 'unweighted'
        # save the cv scores
        f = f'{const.P["outputPath"]}/{dest_dir_appdx}/cv_laplace_kernel_SVM_{weighted}.csv'
        scores.to_csv(f)
        print(f, end='\n\n\n')
        return scores

    def plot_hyperparamter_search_result(scores_w, scores_unw):
        plt.subplots(figsize=(8,7))
        plt.xlim(0,1.01)
        plt.ylim(0,1.01)
        plt.xlabel('Recall')
        plt.ylabel('Presicion')
        plt.grid()
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.title('Cross Validated Kernel SVM\nred: weighted examples, blue: unweighted', pad=10)
        
        scores_w = scores_w.reindex(scores_w.f1.sort_values(ascending=False).index)
        # print(scores_w)
        for at, i in enumerate([0,2,17]):
            score = scores_w.iloc[i]
            score_str = f'presic.: {score.presicion:.3f}\nrecall: {score.recall:.3f}\nF1: {score.f1:.3f}'
            plt.annotate(score_str, (score.recall, score.presicion), (.5+at*.17, .9-.02*at),
                        fontsize=8.5, arrowprops={'arrowstyle': '->'},)
            print(score, end='\n\n')
        plt.scatter(scores_unw.recall, scores_unw.presicion, s=3, alpha=.7, color='b')
        plt.scatter(scores_w.recall, scores_w.presicion, s=3, alpha=.7, color='r')
        
        plt.savefig(f'{const.P["outputPath"]}/{dest_dir_appdx}/laplce_kernel_SVM_cv.png')
        print(f'{const.P["outputPath"]}/{dest_dir_appdx}/laplce_kernel_SVM_cv.png')

    
    print('\n\n\n')
    # get the rasterplot labels
    labels_file = f'{training_data_dir}/onset_offset_labels.tsv'
    labels = pd.read_csv(labels_file, sep='\t', index_col=0)
    y = labels.label.copy()
        
    # set the sample weights based on the label, then make label binary
    sample_weights = np.ones(len(y))
    sample_weights[y==3] = 2
    sample_weights[y==2] = 4
    sample_weights[y==1] = 8
    y[y==3] = 1
    y[y==2] = 1

    # get the features, ie the histogram counts
    hist_bins_file = f'{training_data_dir}/onset_offset_spikebins_channels.csv'
    X = pd.read_csv(hist_bins_file, index_col=0)
    
    # some hist bins might be all 0 and not in the hist_bins_file, add here
    missing_bins = pd.DataFrame(0, index=y.index.difference(X.index),
                                columns=X.columns)
    X = pd.concat((X, missing_bins)).reindex(y.index)
    print('Training data loaded')

    # found optimal classifier
    weighted, gamma, C, bound = True, 1.1, 2.5, .1
    # Recall 0.953, Presicion 0.625 F1: 0.725

    if parameter_search:
        # calls all of the functions above
        print('Starting CV hyperparamter search')
        find_hyperparamters(weighted=True)
        find_hyperparamters(weighted=False)
    
    if plot_cv_result:
        # read in the csv outputfile prodcued by find_hyperparameter and plot it
        scores_w = pd.read_csv(f'{const.P["outputPath"]}/{dest_dir_appdx}/cv_laplace_kernel_SVM_weighted.csv', index_col=0)
        scores_unw = pd.read_csv(f'{const.P["outputPath"]}/{dest_dir_appdx}/cv_laplace_kernel_SVM_unweighted.csv', index_col=0)
        plot_hyperparamter_search_result(scores_w, scores_unw)
    
    if analyize_confusions:
        print('Analyzing misclassifications')
        split = train_test_split(X, y, sample_weights)
        X_train, X_test, y_train, y_test, sample_weights_train, _ = split
        y_pred_p, y_pred = predict(X_train, X_test, y_train, None, 
                                   sample_weights_train, weighted, gamma, C, 
                                   bound)

        y_pred_p = pd.Series(y_pred_p[:,1], index=X_test.index, name='prob')
        y_pred = pd.Series(y_pred, index=X_test.index, name='label')
        pres, recall, f1, _ = precision_recall_fscore_support(y_test, y_pred)
        print(f'presic.: {pres[1]:.3f}  recall: {recall[1]:.3f}  F1: {f1[1]:.3f}\n')

        false = (y_pred != y_test)
        false_pos = y_pred[false.values & (y_pred==1).values]
        false_neg = y_pred[false.values & (y_pred==0).values]

        rasters = ' '.join(labels.reindex(false_pos.index).file)
        print(f'False Positives ({len(false_pos)}): \ngwenview {rasters}\n')

        rasters = ' '.join(labels.reindex(false_neg.index).file)
        print(f'False Negatives ({len(false_neg)}): \ngwenview {rasters}\n')

    if pred_X is not None:
        print('Predicting unseen data')
        y_pred_p, y_pred = predict(X, pred_X, y, None, sample_weights, 
                                   weighted, gamma, C, bound)

        y_pred_p = pd.Series(y_pred_p[:,1], index=pred_X.index, name='prob')
        y_pred = pd.Series(y_pred, index=pred_X.index, name='label')
        return pd.concat((y_pred_p, y_pred), axis=1)









def get_onset_offset_classification(which_data, training_data_dir, dest_dir_appdx,
                                    cached_prediction=True, training_data_chnl_map_file='', 
                                    MUA_output_data_chnl_map_file='', keep_labels=[1,2,3]):
    """The fourth function in the onset-offset pipeline. THis is used to fetch 
    the onset-offset data and slice it/ format it for plotting. The core 
    parameter passed is `which_data` which should be passed as `training`, 
    `MUA_output` or `both`. `training` fetches the onset-offset examples the 
    classifier was trained on. `MUA_output` is intended to be used to classify
    unseen data. Thus the process of using the pipeline invloves training the
    classifier by setting the P["outputPath"] to the training data and also 
    changing the MUA_constants.ALL_MICE to the mice serving as training
    examples. Then, the user changes both of the above to the unseen data in 
    which case this function (together with passing `which_data='MUA_output') 
    will generate the feature vectors and have the trained SVM classify the 
    unseen data. Importantly, the rasterplot classifications should be verified 
    by eye. The function will create a directory at
    P["outputPath"]}/dest_dir_appdx/predicted_rasters and copy all positive 
    rasterplot examples into it. Then, the user should go through each plot and
    label it by prepanding a 1,2,3 or 0 to the front of the rasterplots filename.
    In the next run of this function the new labels indicated by the rasterplots
    filenames will be used. Otherwise, the original classifier labels will used
    which probably include 40% false positives.
    The last option, `both` reads in both training examples and
    unseen examples to merge them. Finally, the labeled examples are returned
    as a pd.DataFrame where rows refer to examples, the columns label, 
    propability, rasterplot_filename, mouse, paradigm, stimulus_type and channel 
    indicate the positive examples identity. Again, the `which_data` argument 
    fundamentally controls which examples will bin included in this returned 
    DataFrame.`training_data_dir` is again essential no matter what is passed in 
    `which_data`. Same as described in the previous function of the pipeline.  
    Since prediction and feature generatation might be computationally intense,
    `cached_prediction` provides the option to use the previous prediction 
    instead of computing a new one. When this is False or no cached prediction
    exists, a new prediction is computed and cached. The caching involves simply
    saving the prediction output at 
    P["outputPath"]/dest_dir_appdx/cached_prediction.csv. Channels may be mapped
    to regions using the `training_data_chnl_map_file` and 
    `MUA_output_data_chnl_map_file`. If None are passed, channels are not 
    mapped. For the format of these files, check the working examples of this
    pipeline (training one and MUA_data one differ slightly!). `keep_labels` is 
    the final option to control which examples are included in the output. 
    Pass a list of labels to be included by default [1,2,3] (all positives)
    
    """

    # training data
    if which_data in ('training', 'both'):
        labels_file = f'{training_data_dir}/onset_offset_labels.tsv'
        y_train = pd.read_csv(labels_file, sep='\t', index_col=0)
            
        # split the index indicating mouse-paradigm-stimtype-channel into the 
        # single components and add them as columns
        split_labels = pd.DataFrame(np.stack([lbl.split('-') for lbl in y_train.index], axis=0), 
                                    columns=['mouse', 'paradigm', 'stimulus_type', 'channel'],
                                    index=y_train.index)
        y_train = pd.concat((y_train, split_labels), axis=1)
        # reorder to label importance
        order = [idx for i in (1,3,2,0) for idx in y_train.index[y_train.label==i]]
        y_train = y_train.reindex(order)

        # Use the previously defined channelmapping to convert from channel to region
        if training_data_chnl_map_file:
            chnl_map = pd.read_csv(training_data_chnl_map_file, index_col=0)
            # the training data was processed using mapped channels, ie channel=1
            # refers to the first channel of the shank. Therefore, when indexing
            # the channel map csv, we use iloc instead of loc (the channelmap is 
            # ordered as the physical electrode)
            y_train.channel = [chnl_map[value.mouse+'-'+value.paradigm].iloc[int(value.channel)-1] 
                               for idx, value in y_train.iterrows()]
        

    # MUA_output data, ie. all the mice defined in MUA_constants, intended to be
    # the new data the SVM should classify
    if which_data in ('MUA_output', 'both'):

        # pepare a tmp folder in the base dir to save the prediction.csv
        file = f'{const.P["outputPath"]}/{dest_dir_appdx}/cached_prediction.csv'

        # Run the prediction, including generating features
        if not cached_prediction or not os.path.exists(file):
            print('Computing features...\n')
            X = onset_offset_response(None, generate_plots=False)
            prediction = lapl_kernel_SVM(training_data_dir=training_data_dir, pred_X=X)

            split_labels = pd.DataFrame(np.stack([lbl.split('-') for lbl in prediction.index], axis=0),
                                        columns=['mouse', 'paradigm', 'stimulus_type', 'channel'],
                                        index=prediction.index)
            
            # get the corresbonding rasterfiles, reconstruct path from the index
            #  and const.MICE_DATES
            get_raster = lambda m_id, parad, stim_t, channel: (f'{const.P["outputPath"]}'
                                f'/{const.MICE_DATES[m_id]}_{parad}.mcd/'
                                f'Raster_NegativeSpikes_Triggers_{stim_t}'
                                f'_ElectrodeChannel_{channel}.png')
            prediction['file'] = [get_raster(*lbl) for lbl in split_labels.values]
            prediction = pd.concat((prediction, split_labels), axis=1)
            prediction = prediction.sort_values('prob', ascending=False)
            prediction.to_csv(file)
            print(f'New cached prediction saved at {file}')
        else:
            prediction = pd.read_csv(file, index_col=0)

        # same as for training data, but here, loc is used to get the region 
        # from the mapping csv. This is because the prediction channel index
        # actually refers to the unmapped channel, so it matches the mapping.csv
        if MUA_output_data_chnl_map_file:
            chnl_map = pd.read_csv(MUA_output_data_chnl_map_file, index_col=0)
            prediction.channel = [chnl_map.loc[value.channel, value.mouse] 
                                    for idx, value in prediction.iterrows()]
    
        # here the predicted labels are preparred for human double checking.      
        # get the positively predicted rasterplots
        pos_data = prediction[prediction.label==1]
        rasters = pos_data.file

        # make a new dictionary and copy all the rasters in it with a mofified
        # name. The prependix _000_*filename, _001_*filename ... is attachted to
        # beginning of the filename the number maintains the order, and before 
        # the underscore the human should put in the true label of this raster
        relabeled_dir = f'{const.P["outputPath"]}/{dest_dir_appdx}/predicted_rasters'
        os.makedirs(relabeled_dir, exist_ok=True)
        relabeled_rasters = [f'{relabeled_dir}/_{i+1:0>3d}_{os.path.basename(raster)}'
                             for i, raster in enumerate(rasters)]

        relabeled_raster_files = os.listdir(relabeled_dir)
        # if there is not files in this dir yet, copy them in with the modified name
        if not relabeled_raster_files:
            print('Copying raster plots into directory, please label')
            [copyfile(raster, raster_ranked) for raster, raster_ranked in zip(rasters, relabeled_rasters)]
        # if there are already copied files in there, read them in and sort 
        else:
            print('Copied raster plots found in directory')
            # reorder to original order, crucial
            relabeled_rasters = sorted(relabeled_raster_files, key= lambda str: str[1:])

        # if the rasters have been lablled (no files starts with _ but  1,2,3,0)
        if not any([True if raster[0] == '_' else False for raster in relabeled_rasters]):
            print('Lableled raster plots found in directory, setting new labels according to it')
            corrected_lbls = [int(raster[0]) for raster in relabeled_rasters]
            corrected_lbls = pd.DataFrame({'label': corrected_lbls, 'file': rasters}, 
                                          index=rasters.index)
            prediction.loc[corrected_lbls.index, 'label'] = corrected_lbls
            corrected_lbls.to_csv(f'{const.P["outputPath"]}/{dest_dir_appdx}/pred_corrected_labels.tsv', sep='\t')
        else:
            print('Raster plots haven`t all been labeled, using orignial SVM labels.')
        
    # final data is simply the training data 
    if which_data == 'training':
        data = y_train
    # just the prediction of the data
    elif which_data == 'MUA_output':
        data = prediction
    # merges both training and prediction data to one chunk (y_train does not 
    # have probabilities associated with it, in contrast to `prediction`)
    elif which_data == 'both':
        data = prediction.append(y_train)
    
    # slice the data to the passed labels, and to the mouse, paradigm, stimtype, channel
    data = data[[True if lbl in keep_labels else False for lbl in data.label]]
    print(data)
    return data

def onoff_heatmap(data, rank='', plot_labeled_data=False, print_pos_rasters_pred=False,
                  print_labeled_data=False, split_mice=False, cluster=False):

    def do_plot(dat, y):
        if rank:
            sort_data = dat[rank]
            groups = [sort_data[sort_data==group].index for group in np.unique(sort_data.values)]
            group_sizes = [len(group) for group in groups]
            groups = [group for _, group in sorted(zip(group_sizes, groups), key=lambda pair: pair[0])]
            dat = dat.reindex(np.concatenate(groups)[::-1])
            if plot_labeled_data:
                y = y.reindex(dat.index)
        elif cluster:
            print(dat)
            enc = OneHotEncoder()
            dat_enc = enc.fit_transform(dat).toarray()

            d = dat_enc
            Y = pdist(d, metric='euclidean')
            Z = linkage(Y, method='complete', metric='euclidean')
            den = dendrogram(Z,
                             count_sort = True,
                             no_labels = True,
                             orientation = 'right',
                             labels = dat.index.values)
            order = den['ivl']
            # print(order)
            dat = dat.reindex(order)


        nrows, ncols = dat.shape
        height = nrows*.08 if nrows*.08 > 9 else 9
        fig, ax = plt.subplots(figsize=(ncols*1.8, height))
        [sp.set_visible(False) for sp in ax.spines.values()]
        ax.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False,
                        top=True, labeltop=True)
        fig.subplots_adjust(left=.03, right=.7, bottom=.01, top=.92)
        fig.suptitle('Positive Onset-Offset examples', y=.99, fontsize=13)

        widths = [1,1,1.5,1]
        for row in range(nrows):
            left = 0
            for col, width in zip(range(ncols), widths):
                value = dat.iloc[row, col]
                ax.barh(y=row, width=width, height=1, left=left, align='edge', 
                        color=const.GENERAL_CMAP[value])
                left += width
                
        if plot_labeled_data:
            ax.set_yticklabels(y.values, fontsize=7)
            ax.set_yticks(np.arange(.5, nrows+.5))
            ax.tick_params(labelleft=True)

        ax.set_ylim(nrows, 0)
        ax.set_xlim(0, sum(widths))
        ax.set_xticks((.5,1.5,2.75,4))
        lbls = 'MOUSE', 'PARADIGM', 'STIMULUS TYPE', 'REGION'
        ax.set_xticklabels((lbls), fontsize=12.5)

        mice_legend = [(key, const.GENERAL_CMAP[key]) for key in const.ALL_MICE+['mGE82', 'mGE83', 'mGE84', 'mGE85']]
        parad_legend = [(key, const.GENERAL_CMAP[key]) for key in const.ALL_PARADIGMS]
        region_legend = [(key, const.GENERAL_CMAP[key]) for key in const.REGIONS.keys()]
        stimt_legend = (('Standard', const.GENERAL_CMAP['Standard']),
                        ('Deviant', const.GENERAL_CMAP['Deviant']),
                        ('(Unqiue)Predeviant', const.GENERAL_CMAP['Predeviant']),
                        ('(Unqiue)Postdeviant', const.GENERAL_CMAP['Postdeviant']),
                        ('MS (C1, C2, D1, B1)', const.GENERAL_CMAP['C1']))
        legends = [mice_legend, parad_legend, stimt_legend, region_legend]

        for at_y, which_legend in zip((.92,.52, .3, .17), range(4)):
            legend = legends[which_legend]
            handles = [Patch(color=legend[j][1], label=legend[j][0]) 
                    for j in range(len(legend))]
            fig.legend(handles=handles, loc='upper left', ncol=1,
                    bbox_to_anchor=(.7, at_y))
            ax.annotate(lbls[which_legend], (.7+.02, at_y), ha='left', va='bottom', 
                        xycoords='figure fraction', fontsize=12.5)
        return fig

    def do_barplot(dat, feature, scnd_feature):
        realizs = np.unique(prediction[feature])
        scnd_realizs = np.unique(prediction[scnd_feature])

        props = [(dat[feature]==realiz).sum()/ (prediction[feature]==realiz).sum() for realiz in realizs]
        realizs = [realz for _, realz in sorted(zip(props, realizs), key=lambda pair: pair[0], reverse=True)]
        props.sort(reverse=True)
        
        
        scnd_props = []
        scnd_realizs_sorted = []
        for realiz in realizs:

            scnd_prop = [((dat[scnd_feature]==scnd_realiz).values & (dat[feature]==realiz).values).sum() / (prediction[feature]==realiz).sum() for scnd_realiz in scnd_realizs]
            scnd_realizs = [sncd_realz for _, sncd_realz in sorted(zip(scnd_prop, scnd_realizs), key=lambda pair: pair[0])]
            scnd_prop.sort()
            
            scnd_realizs_sorted.append(scnd_realizs)
            scnd_props.append(np.array(scnd_prop))

        # scnd_props = [np.sort(arr) for arr in scnd_props]


        # print(realizs)
        # print(scnd_realizs)
        # print(props)
        # print([sum(l) for l in scnd_props])
        # exit()

        fig, ax = plt.subplots(figsize=(11,6))
        fig.subplots_adjust(bottom=.12)


        cols = [const.GENERAL_CMAP[key] for key in realizs]
        if feature == scnd_feature:
            ax.bar(np.arange(len(props)), props, color=cols, edgecolor='grey')
            ax.set_title(feature+' - proportion of onset-offset observations')
        else:
            ax.set_title(f'{scnd_feature} in {feature} - proportion of onset-offset observations')
            
            for which_bar in np.arange(len(props)):
                bottom = 0
                colors = [const.GENERAL_CMAP[key] for key in scnd_realizs_sorted[which_bar]]
                for col, scnd_bar_x in zip(colors, scnd_props[which_bar]):
                    ax.bar(which_bar, scnd_bar_x, bottom=bottom, edgecolor='grey', color=col)
                    bottom += scnd_bar_x
            
                legend = [(key, const.GENERAL_CMAP[key]) for key in scnd_realizs_sorted[which_bar]]
                handles = [Patch(color=legend[j][1], label=legend[j][0]) 
                           for j in range(len(legend))]
                fig.legend(handles=handles, loc='upper right', ncol=1,
                           bbox_to_anchor=(.9, .84))
                ax.annotate(scnd_feature.upper(), (.89, .84), ha='right', va='center', 
                            xycoords='figure fraction', fontsize=12.5)


        ax.set_xticks(np.arange(len(props)))
        ax.set_xticklabels(realizs, rotation=30, rotation_mode='anchor', ha='right')
        ax.set_ylabel('Proportions')
        ax.set_xlabel(feature.upper())
        return fig


    fname_labled = 'true_labels' if plot_labeled_data else ''
    fname_rank = f'{rank}_sorted' if rank else ''
    if not split_mice:
        # fig = do_plot(data, y)
        # f = f'{const.P["outputPath"]}/{dest_dir_appdx}/positives_{fname_labled}_{fname_rank}.png'
        # fig.savefig(f)

        for feature in data.columns:
            for scnd_feature in data.columns[:]:
                data = data.reindex(data[scnd_feature].sort_values().index)
                fig = do_barplot(data, feature, scnd_feature)
                f = f'{const.P["outputPath"]}/{dest_dir_appdx}/proprtions_{feature}-{scnd_feature}.png'
                print(f)
                fig.savefig(f)
    else:
        for mouse in np.unique(data.mouse.values):
            mouse_dat = data[data.mouse==mouse]
            fig = do_plot(mouse_dat, y.reindex(mouse_dat.index))
            f = f'{const.P["outputPath"]}/{dest_dir_appdx}/{mouse}_positives_{fname_labled}_{fname_rank}.png'
            fig.savefig(f)





# _200_Raster_NegativeSpikes_Triggers_Deviant_ElectrodeChannel_21.png
# _231_Raster_NegativeSpikes_Triggers_C1_ElectrodeChannel_21.png
# _342_Raster_NegativeSpikes_Triggers_Deviant_ElectrodeChannel_12.png