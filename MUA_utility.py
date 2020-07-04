import numpy as np
import pandas as pd  
from glob import glob
from collections import OrderedDict
import os


import MUA_constants as const

def fetch(mouseids=const.ALL_MICE, paradigms=const.ALL_PARADIGMS, stim_types=const.ALL_STIMTYPES, 
          collapse_ctx_chnls=True, collapse_th_chnls=False, drop_not_assigned_chnls=False):
    """Get the processed data by passing the mice-, paradigms-, and stimulus
    types of interst from the saved .gzip`s. Returns a dictionary with key: 
    mouseid-paradigm-stimulus_type and value: (firingrate_df, summary_df, 
    pos_spike_df, neg_spike_df)."""
    
    path = const.P['outputPath']
    if collapse_ctx_chnls or collapse_th_chnls:
        chnls_map = pd.read_csv(path+'/../chnls_map.csv', index_col=0)
    data = OrderedDict()

    invalid_mid = [mouse for mouse in mouseids if mouse not in const.ALL_MICE]
    invalid_pard = [parad for parad in paradigms if parad not in const.ALL_PARADIGMS]
    invalid_stimt = [stimt for stimt in stim_types if stimt not in const.ALL_STIMTYPES]
    if any(invalid_mid+invalid_pard+invalid_stimt):
        err = (f'Invalid data request: mice: {invalid_mid}\nparadigms: '
               f'{invalid_pard}\nstim_types: {invalid_stimt}')
    
    # iterate over passed mouse id's
    for m_id in mouseids:
        # get all dirs containing the mouse id
        mouse_files = glob(f'{path}/*{m_id}*') 

        # iterate paradims
        for parad in paradigms:
            # get the one dir matching the mouse_id and paradigm
            parad_dir = [f for f in mouse_files if parad in f][0]

            # iterate the stimulus types, eg `Deviant`, `Predeviant`... for MS `C1`...
            for stim_t in stim_types:
                # get a list of all CSVs for mouse-paradigm-stim_type
                stim_files = glob(f'{parad_dir}/*{stim_t}*.gzip')

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
                        lfp_summary, lfp = np.split(pd.read_csv(lfp_avg_csv, index_col=0), [7], axis=1)

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
                        for df in [frate, mua_summary, pos_spikes, neg_spikes, lfp, lfp_summary]:
                            if df is not pos_spikes and df is not neg_spikes:
                                if df is None:  # O25 pre, postdeviant
                                    data[key].append(None)
                                    continue
                                region_collps = [pd.Series(df.iloc[chnls].mean(), name=region) 
                                                for region, chnls in region_map.items() 
                                                if any(chnls)]
                                df_region_idx = pd.concat(region_collps, axis=1).T
                                if not drop_not_assigned_chnls:
                                    assigned = [chnl for chnls in region_map.values() for chnl in chnls]
                                    df_not_ass = df.drop(assigned)
                                    df_region_idx = pd.concat((df_region_idx, df_not_ass))
                                data[key].append(df_region_idx)
                            else:
                                # spacers for pos_spikes and neg_spikes, bc not implemented
                                data[key].append(None)
                                pass
    return data

def slice_data(data, mouseids=const.ALL_MICE, paradigms=const.ALL_PARADIGMS, 
               stim_types=const.ALL_STIMTYPES, firingrate=False, mua_summary=False, 
               pos_spikes=False, neg_spikes=False, lfp=False, lfp_summary=False,
               drop_labels=False, frate_noise_subtraction=True):
    """Convenient`s data selection function. Takes in the data (obtained from 
    fetch()) and returns the subset of interst, eg. a specific mouse/stimulus 
    type combination and one of the 4 datatypes (eg. the firing rate). Returns
    the sliced data in the original input format (dict)."""
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
        # return corr_firingrate.mask(corr_firingrate < 0, 0)
        return corr_firingrate

         
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


def compute_si(data):
    parads = [key[key.find('-')+1:key.rfind('-')] for key in data.keys()]
    parads = list(dict().fromkeys(parads))
    mice = [key[:key.find('-')+1] for key in data.keys()]
    mice = list(dict().fromkeys(mice))

    post_stim = [str(float(time_bin)) for time_bin in range(0, 20, 5)]

    parads_paris = [(c1, c2) for c1, c2 in const.PARAD_PAIRS if c1 in parads and c2 in parads]

    SI_values = dict()
    for parad_pair in parads_paris:
        for m_id in mice:
            dat = slice_data(data, mouseids=m_id, paradigms=parad_pair, 
                             firingrate=True, frate_noise_subtraction='deviant_alone')
            stim_ts = [key[key.rfind('-')+1:] for key in data.keys() if 'Deviant' not in key]
            stim_ts = list(dict().fromkeys(stim_ts))
            
            c1_dev, c2_dev = [fr[post_stim] for key, fr in dat.items() if 'Deviant' in key]
            c1c2_dev_resp = ((c1_dev+c2_dev) /2).mean().mean()
            c1_dev_resp = c1_dev.mean().mean()
            c2_dev_resp = c2_dev.mean().mean()

            for stim_t in stim_ts:
                c1, c2 = [fr[post_stim] for key, fr in dat.items() if stim_t in key]
                c1c2_resp = ((c1+c2) /2).mean().mean()
                c1_resp = c1.mean().mean()
                c2_resp = c2.mean().mean()

                SI = (c1c2_dev_resp - c1c2_resp) / (c1c2_dev_resp + c1c2_resp)

                SI_values[f'{m_id}-{parad_pair[0][:-2]}-{stim_t}'] = SI

    return SI_values
