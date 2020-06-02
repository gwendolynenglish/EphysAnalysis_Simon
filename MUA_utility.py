import numpy as np
import pandas as pd  
from glob import glob

from MUA_constants import ALL_MICE, ALL_PARADIGMS, ALL_STIMTYPES

def fetch(p, mouseids=ALL_MICE, paradigms=ALL_PARADIGMS, stim_types=ALL_STIMTYPES):
    """Get the processed data by passing the mice-, paradigms-, and stimulus
    types of interst from the saved .gzip`s. Returns a dictionary with key: 
    mouseid-paradigm-stimulus_type and value: (firingrate_df, summary_df, 
    pos_spike_df, neg_spike_df)."""
    path = p['outputPath']
    data = dict()
    
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
                    summary = pd.read_pickle(summary_file)

                    # get 32 (pos or neg) spike.csv`s 
                    pos_spike_file = [f for f in stim_files if 'pos' in f][0]
                    neg_spike_file = [f for f in stim_files if 'neg' in f][0]
                    pos_spikes = pd.read_pickle(pos_spike_file)
                    neg_spikes = pd.read_pickle(neg_spike_file)

                    key = f'{m_id}-{parad}-{stim_t}'
                    data[key] = [frate, summary, pos_spikes, neg_spikes]
    return data

def slice_data(data, mouseids=ALL_MICE, paradigms=ALL_PARADIGMS, 
               stim_types=ALL_STIMTYPES, firingrate=False, summary=False, 
               pos_spikes=False, neg_spikes=False):
    """Convenient`s data selection function. Takes in the data (obtained from 
    fetch()) and returns the subset of interst, eg. a specific mouse/stimulus 
    type combination and one of the 4 datatypes (eg. the firing rate). Returns
    the sliced data in the original input format (dict)."""
    mask = [firingrate, summary, pos_spikes, neg_spikes]

    # check that exactly one type of data is set to True
    if not any(mask):
        print('Failed to slice data: at least one of `firing_rate`, `summary`, '
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
            new_data[key] = [data[key][i] for i in range(len(mask)) if mask[i]][0]
    return new_data