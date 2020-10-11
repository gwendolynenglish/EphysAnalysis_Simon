import numpy as np
import pandas as pd  
from collections import OrderedDict
from glob import glob
import os
from shutil import copyfile, move

import warnings
import concurrent.futures

import MUA_constants as const

def process_data(how={}):
    """
    Runs the processing of the .mcd dirs. If `how` is not passed/ empty, simply
    iterate through all dirs. If `nbatches` and `batch` passed, the dirs are 
    split into nbatches and only the subset batch is run. batch is simply the 
    index of the split dirs. If `threads` is passed, the value of threads will
    be passed the the concurrent.futures.ProcessPoolExecutor to process many dirs
    in parallel. 
    """
    from MUA_cycle_dirs import MUA_analyzeMouseParadigm
    warnings.filterwarnings('ignore')
    dirs = os.listdir(const.P['inputPath'])

    if 'threads' in how:
        with concurrent.futures.ProcessPoolExecutor(max_workers=how['threads']) as executer:
            [executer.submit(MUA_analyzeMouseParadigm, folder) for folder in dirs]
    elif 'nbatches' in how and 'batch' in how:

        batches = np.array_split(dirs, how['nbatches'])
        batch = batches[how['batch']]
        str_batch = '\n\t'.join(batch)
        s = f'====== Processing batch {how["batch"]}, n={len(batch)} =======\n\t{str_batch}'
        print(s)

        [MUA_analyzeMouseParadigm(folder) for folder in batch]
    else:
        dirs_string = "\n".join(dirs)
        print(f'Processing the following files:\n{dirs_string}')
        [MUA_analyzeMouseParadigm(folder) for folder in dirs]
    warnings.filterwarnings('default')
    

def compress_CSVs():
    """Reads in all the CSVs produced and saves them as gzip binary. The 32(pos)
    + 32(neg) spike timestamp CSVs are merged into a pos and neg gzip`ed frame.
    """
    
    # iterate over passed mouse id's
    print('Compressing CSVs and merging channel wise spike time stamps...')
    for m_id in const.ALL_MICE:
        # get all dirs containing the mouse id
        mouse_files = glob(f'{const.P["outputPath"]}/*{m_id}*') 
        print(f'{m_id} paradigms found: {len(mouse_files)}')

        # iterate paradims
        for parad in const.ALL_PARADIGMS:
            # get the one dir matching the mouse_id and paradigm
            parad_dir = [f for f in mouse_files if parad in f]
            if not parad_dir:
                print(f'\t{parad} not found. Skipping {parad}.',)
                continue
            else:
                parad_dir = parad_dir[0]

            # iterate the stimulus types, eg `Deviant`, `Predeviant`... for MS `C1`...
            for stim_t in const.ALL_STIMTYPES:
                # get a list of all CSVs for mouse-paradigm-stim_type
                stim_files = glob(f'{parad_dir}/*{stim_t}*.csv')

                # only enter when the stim_type applies in the current paradigm
                if stim_files:
                    # get firingrates.csv, save as gzip`ed pandas frame
                    fr_file = [f for f in stim_files if 'FiringRate' in f][0]
                    compr_fn = fr_file[:-4] + '.gzip'
                    pd.read_csv(fr_file, index_col=0).to_pickle(compr_fn)

                    # get summary.csv, save as gzip`ed pandas frame
                    summary_file = [f for f in stim_files if 'Summary' in f][0]
                    compr_fn = summary_file[:-4] + '.gzip'
                    pd.read_csv(summary_file, index_col=0).to_pickle(compr_fn)

                    # aggregate channel specific spike.csv`s in one gzip`ed frame
                    for which in ('pos', 'neg'):
                        # get 32 (pos or neg) spike.csv`s 
                        spike_files = [f for f in stim_files if which in f]
                        
                        # read in csv, delete default 0-column, add Multiindex
                        # with channel at level 0, and spike# at level 1
                        def format_spikes(f):
                            channel = int(f[-19:-17])
                            df = pd.read_csv(f).iloc(1)[2:]
                            if df.empty:
                                return
                            multi_idx = [[channel], df.columns]
                            df.columns = pd.MultiIndex.from_product(multi_idx,
                                                     names=('channel', 'spike'))
                            return df
                        all_spikes  = [format_spikes(f) for f in spike_files]
                        # merge into one df, save gzip
                        compr_fn = spike_files[0][:-36] + f'{which}Spikes.gzip'
                        pd.concat(all_spikes, axis=1).to_pickle(compr_fn)
    print('Done.')




def fetch(mouseids=const.ALL_MICE, paradigms=const.ALL_PARADIGMS, 
          stim_types=const.ALL_STIMTYPES, collapse_ctx_chnls=False, 
          collapse_th_chnls=False, drop_not_assigned_chnls=False):
    """
    Get the processed data by passing the mice-, paradigms-, and stimulus
    types of interst from the saved .gzip`s. Returns a dictionary with key: 
    mouseid-paradigm-stimulus_type and value: (firingrate_df, summary_df, 
    pos_spike_df, neg_spike_df, lfp, lfp_summary). To have the last two elements 
    in the list not set to None, MUA_const.LFP_OUTPUT needs to direct to the 
    output of the LFP analysis pipeline. The DataFrames may be converted from 
    having the 32 channels as rows, to the mapped region indicated by a 
    mapping file. This mapping file should look like this: 
    P["outputPath"]/../chnls_map.csv. When channels are collased to a region,
    the arithmetic mean is taken for firingrate and summary stats. Spike time
    stamps lists of different channels are simply merged together. The mapping
    is controlled by the parameters `collapse_ctx_chnls`, `collapse_th_chnls`, 
    and `drop_not_assigned_chnls`, which should be self explanatory.
    """
    
    if collapse_ctx_chnls or collapse_th_chnls:
        chnls_map = pd.read_csv(f'{const.P["outputPath"]}/../chnls_map.csv', index_col=0)
    data = OrderedDict()

    invalid_mid = [mouse for mouse in mouseids if mouse not in const.ALL_MICE]
    invalid_pard = [parad for parad in paradigms if parad not in const.ALL_PARADIGMS]
    invalid_stimt = [stimt for stimt in stim_types if stimt not in const.ALL_STIMTYPES]
    if any(invalid_mid+invalid_pard+invalid_stimt):
        err = (f'Invalid data request: mice: {invalid_mid}\nparadigms: '
               f'{invalid_pard}\nstim_types: {invalid_stimt}')
        print(err)
        exit()
    
    # iterate over passed mouse id's
    for m_id in mouseids:
        # get all dirs containing the mouse id
        mouse_files = glob(f'{const.P["outputPath"]}/*{m_id}*') 

        # iterate paradims
        for parad in paradigms:
            # get the one dir matching the mouse_id and paradigm
            parad_dir = [f for f in mouse_files if parad in f]
            if not parad_dir:
                continue
            else:
                parad_dir = parad_dir[0]

            # iterate the stimulus types, eg `Deviant`, `Predeviant`... for MS `C1`...
            for stim_t in stim_types:
                # get a list of all CSVs for mouse-paradigm-stim_type
                stim_files = glob(f'{parad_dir}/*_{stim_t}_*.gzip')

                # only enter when the stim_type applies in the current paradigm
                if stim_files:
                    # get firingrates.csv, save as gzip`ed pandas frame
                    fr_file = [f for f in stim_files if 'FiringRate' in f][0]
                    frate = pd.read_pickle(fr_file)

                    # get summary.csv, save as gzip`ed pandas frame
                    summary_file = [f for f in stim_files if 'Summary' in f][0]
                    mua_summary = pd.read_pickle(summary_file)

                    # get 32 (pos or neg) spike.csv`s 
                    pos_spike_file = [f for f in stim_files if 'pos' in f][0]
                    neg_spike_file = [f for f in stim_files if 'neg' in f][0]
                    pos_spikes = pd.read_pickle(pos_spike_file)
                    neg_spikes = pd.read_pickle(neg_spike_file)

                    # lfp summary and lfp time series
                    # for oddball 25, only UniquePredeviant/Postd. exists 
                    if 'O25' in parad and stim_t in ('Predeviant', 'Postdeviant'):
                        lfp_summary, lfp = None, None
                    else:
                        lfp_avg_csv = f'{const.LFP_OUTPUT}/{os.path.basename(parad_dir)}/Triggers_{stim_t}_LFPAverages.csv'
                        if os.path.exists(lfp_avg_csv):
                            lfp_summary, lfp = np.split(pd.read_csv(lfp_avg_csv, index_col=0), [7], axis=1)
                        else:
                            # in the case where LFP anaylsis wasn't run at all
                            lfp_summary, lfp = None, None

                    key = f'{m_id}-{parad}-{stim_t}'
                    data[key] = []
                    if not collapse_ctx_chnls and not collapse_th_chnls:
                        data[key] = [frate, mua_summary, pos_spikes, neg_spikes, lfp, lfp_summary]
                    
                    # collapse channels index to region using chnls_map.csv
                    else:
                        this_map = chnls_map.reset_index(drop=True)[f'{m_id}-{parad}']
                        regions = ['SG', 'G', 'IG', 'dIG']
                        if collapse_th_chnls:
                            regions.append('Th')
                        region_map = {region: this_map.index[this_map==region]
                                      for region in regions}
                        assigned = [chnl for chnls in region_map.values() for chnl in chnls]
                        for df in [frate, mua_summary, pos_spikes, neg_spikes, lfp, lfp_summary]:
                            if df is not pos_spikes and df is not neg_spikes:
                                if df is None:  # O25 pre, postdeviant (unique exists)
                                    data[key].append(None)
                                    continue
                                region_collps = [pd.Series(df.iloc[chnls].mean(), name=region) 
                                                for region, chnls in region_map.items() 
                                                if any(chnls)]
                                df_region_idx = pd.concat(region_collps, axis=1).T
                                if not drop_not_assigned_chnls:
                                    df_not_ass = df.drop(assigned)
                                    df_region_idx = pd.concat((df_region_idx, df_not_ass))

                            # spikestamp data 
                            else:
                                def chnl_to_region(df, reg):
                                    df.columns = pd.MultiIndex.from_product([[reg], np.arange(df.shape[1])+1])
                                    return df
                                df_region_idx = [chnl_to_region(df.loc[:, chnls], region) 
                                                 for region, chnls in region_map.items() if any(chnls)]
                                df_region_idx = pd.concat(df_region_idx, axis=1)
                                if not drop_not_assigned_chnls:
                                    df_not_ass = df.drop(assigned, axis=1, level=0)
                                    df_region_idx = pd.concat((df_region_idx, df_not_ass), axis=1)

                            data[key].append(df_region_idx)
    return data

def slice_data(data, mouseids=const.ALL_MICE, paradigms=const.ALL_PARADIGMS, 
               stim_types=const.ALL_STIMTYPES, firingrate=False, mua_summary=False, 
               pos_spikes=False, neg_spikes=False, lfp=False, lfp_summary=False,
               drop_labels=False, frate_noise_subtraction='paradigm_wise'):
    """Convenient`s data selection function. Takes in the data (obtained from 
    fetch()) and returns the subset of interst, eg. a specific mouse/stimulus 
    type combination and one of the 6 datatypes eg. by passing firingrate=True.
    Only one data modality can be returned per function call. 
    By default returns the sliced data in the original input format (dict). If 
    `drop_labels` is passed only the values of the dictionary are returned. 
    `frate_noise_subtraction` is only relevent for the case where 
    `firingrate`=True. The string you pass here is passed to the 
    `subtract_noise()` function as the `method` parameter. Valid methods right 
    now are `paradigm_wise` and `deviant_alone`. Check the `subtract_noise()` 
    function for details.
    """
    mask = [firingrate, mua_summary, pos_spikes, neg_spikes, lfp, lfp_summary]
    # check that exactly one type of data is set to True
    if not any(mask):
        print('Failed to slice data: at least one of `firing_rate`, `mua_summary`, '
              '`pos_spikes`, `neg_spikes`, must be retrieved.')
        exit()
    elif sum(mask) == 2:
        print('Failed to slice data: more then one datatype was requested. '
              'Can only retrieve one type of data per call.')
        exit()

    # iterate all keys in data, select the ones matching the slice
    new_data = dict()
    for key in data.keys():
        m_id, parad, stim_t = key.split('-')
        if m_id in mouseids and parad in paradigms and stim_t in stim_types:
            # built new dict up with values of type pd.DataFrame
            df = [data[key][i] for i in range(len(mask)) if mask[i]][0]
            if firingrate and frate_noise_subtraction:
                df = subtract_noise(df, frate_noise_subtraction, m_id, parad)
            # delete the odd replicated channel, differnt location for neg_spikes
            if firingrate and 14 in df.index:
                df.iloc[14,:] = 0
            if neg_spikes and 15 in df.columns.unique(0):
                df.drop(15, axis=1, inplace=True)
            new_data[key] = df
    
    if not drop_labels:
        return new_data
    else:
        return list(new_data.values())


def subtract_noise(firingrate, method, mouse_id, paradigm):
    """Takes in a firingrate dataframe and subtracts the baseline activity which
    is defined by the `method` passed. Firingrates that fall below 0 after base
    line subtraction are set to 0.  
    method = `deviant_alone` will read in the DA paradigm for the respective 
    mouse and calculate a mean firingrate based on the 50ms pre stimulus. 
    This is done seperately for C1 and C2 paradigms. MS is processed using the 
    mean between C1 and C2 baselines.
    method = `paradigm_wise` will get the standard stimulus type of the paradigm
    the input belongs to. For paradigms without a standard stimulus: deviant 
    alone will use the deviant as a baseline, MS uses the mean prestim 
    firingrate across the four standards."""
    pre_stim = [str(float(time_bin)) for time_bin in range(-50, 1, 5)]
    ctx_idx = any([reg in firingrate.index for reg in ['SG', 'G', 'IG', 'dIG']])
    th_idx = 'Th' in firingrate.index
    only_region_idx = len(firingrate) == len(['SG', 'G', 'IG', 'dIG', 'Th'])
    
    def compute_baseline(paradigms, stim_types):
        data = fetch([mouse_id], paradigms, stim_types, collapse_ctx_chnls=ctx_idx, collapse_th_chnls=th_idx, drop_not_assigned_chnls=only_region_idx)
        base_fr = [dat[0][pre_stim] for dat in data.values()] # dat[0] indexes the firing rate
        # for example MS is based on the average across C1, C2, B1, D1, collapse here
        base_fr = base_fr[0] if len(base_fr) == 1 else sum(base_fr) /len(base_fr)

        # compute ther average over time- and channel domain, subtract from input
        parad_base = base_fr.mean(1).astype(int)
        corr_firingrate = firingrate.apply(lambda time_bin: time_bin-parad_base)
        return corr_firingrate.mask(corr_firingrate < 0, 0)
        # return corr_firingrate

         
    if method == 'deviant_alone':
        if not (paradigm == 'MS'):
            return compute_baseline(['DAC1', 'DAC2'], ['Deviant'])
        elif paradigm == 'MS':
            return compute_baseline(['DAC1', 'DAC2'], ['Deviant'])

    elif method == 'paradigm_wise':
        if paradigm not in ('DAC1', 'DAC2', 'MS'):
            return compute_baseline([paradigm], ['Standard'])
        elif paradigm in ('DAC1', 'DAC2'):
            return compute_baseline([paradigm], ['Deviant'])
        elif paradigm == 'MS':
            return compute_baseline([paradigm], ['C1', 'C2', 'D1', 'B1'])


def compute_si(data, MS=False, start=5, stop=20):
    parads = [key[key.find('-')+1:key.rfind('-')] for key in data.keys()]
    parads = list(dict().fromkeys(parads))
    mice = [key[:key.find('-')] for key in data.keys()]
    mice = list(dict().fromkeys(mice))

    post_stim = [str(float(time_bin)) for time_bin in range(start, stop, 5)]
    parads_paris = [[c1, c2] for c1, c2 in const.PARAD_PAIRS if c1 in parads and c2 in parads]

    SI_values = []
    frates = []
    for parad_pair in parads_paris:
        for m_id in mice:
            print()
            print()
            print()
            print()
            print(m_id)
            dat = slice_data(data, mouseids=m_id, paradigms=parad_pair+['MS'], 
                             firingrate=True, frate_noise_subtraction='paradigm_wise')

            if not MS:         # C1 Standard                C2 Standard
                compare_with = parad_pair[1]+'-Predeviant', parad_pair[0]+'-Predeviant'
            else:
                compare_with = 'MS-C1', 'MS-C2'
            
            c1_dev = dat[f'{m_id}-{parad_pair[0]}-Deviant'][post_stim].mean(1)
            c1_stnd = dat[f'{m_id}-{compare_with[0]}'][post_stim].mean(1)
            # print('c1_dev: ', f'{m_id}-{parad_pair[0]}-Deviant')
            # print(dat[f'{m_id}-{parad_pair[0]}-Deviant'][post_stim])
            # print(c1_dev)
            # print('c1_stnd: ', f'{m_id}-{compare_with[0]}')
            # print(c1_stnd)
            
            c2_stnd = dat[f'{m_id}-{compare_with[1]}'][post_stim].mean(1)
            c2_dev = dat[f'{m_id}-{parad_pair[1]}-Deviant'][post_stim].mean(1)
            # print('c2_stnd: ',f'{m_id}-{compare_with[1]}')
            # print('c2_dev: ', f'{m_id}-{parad_pair[1]}-Deviant')
            
            c1_dev[c1_dev < const.SI_MIN_FRATE_5MS] = 0
            c1_stnd[c1_stnd < const.SI_MIN_FRATE_5MS] = 0
            c2_dev[c2_dev < const.SI_MIN_FRATE_5MS] = 0
            c2_stnd[c2_stnd < const.SI_MIN_FRATE_5MS] = 0

            c1_SI = (c1_dev - c1_stnd) / (c1_dev + c1_stnd)
            c2_SI = (c2_dev - c2_stnd) / (c2_dev + c2_stnd)

            c1_SI.index = pd.MultiIndex.from_product([[m_id], ['C1'], c1_SI.index])
            c2_SI.index = pd.MultiIndex.from_product([[m_id], ['C2'], c2_SI.index])

            SI_values.append(c1_SI)
            SI_values.append(c2_SI)

            c1_SI_str = [f'Dev fr: {dev:.2f}\nStd fr: {std:.2f}' for dev, std in zip(c1_dev.values, c1_stnd.values)]
            c2_SI_str = [f'Dev fr: {dev:.2f}\nStd fr: {std:.2f}' for dev, std in zip(c2_dev.values, c2_stnd.values)]
            frates.append(pd.Series(c1_SI_str, index=c1_SI.index))
            frates.append(pd.Series(c2_SI_str, index=c2_SI.index))
    frates = pd.concat(frates).unstack(level=0).T.swaplevel(axis=1)
    SI_values = pd.concat(SI_values).unstack(level=0).T.swaplevel(axis=1)
    return SI_values.reindex(const.REGIONS.keys(), axis=1, level=0), frates