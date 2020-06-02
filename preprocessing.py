"""
This file contains functions for filtering, downsampling, and aligning the data. 

Created on 02.11.18. 
"""

from scipy import signal
import numpy as np
import pandas as pd
from scipy.signal import decimate   
from glob import glob

from MUA_constants import all_mice, all_parad, all_stimtypes

################################################################################
#Filtering functions

def lowpass_filter(data, samplerate, cutoffFreq, order):
    """Butterworth lowpass filter."""
    b, a = signal.butter(order, cutoffFreq / (samplerate / 2), btype = 'low')
    filtered_signal = signal.filtfilt(b, a, data)
    
    return filtered_signal 


def highpass_filter(data, samplerate, cutoffFreq, order):
    """Butterworth lowpass filter."""
    b, a = signal.butter(order, cutoffFreq / (samplerate / 2), btype = 'high')
    filtered_signal = signal.filtfilt(b, a, data)
    
    return filtered_signal 

################################################################################
#Alignment Functions

def extract_stim_timestamps(stim_array, p):
    #Transform trigger times from seconds to samples 
    stim_timestamps = np.round(stim_array * p['sample_rate']) 

    return stim_timestamps 

def align_to_trigger(filtered_data, stim_timestamps, p):
    #Create array with dimension rows: # of stimulations, columns: window length 
    ncols = int(p['sample_rate']*(p['evoked_pre']+p['evoked_post']))
    evoked_array = np.zeros((len(stim_timestamps), ncols))
    for stim in range(len(stim_timestamps)):
        frm = int(stim_timestamps[stim] - p['evoked_pre']*p['sample_rate'])
        to = int(stim_timestamps[stim] + p['evoked_post']*p['sample_rate'])
        evoked_array[stim,:] = np.array(filtered_data[frm:to])
    return evoked_array 
                             
################################################################################
#Function for making list of uneven-lengthed arrays into matrix (0-filled values)     

def numpy_fillna(array):
    lens = np.array([len(item) for item in array])
    mask = lens[:,None] > np.arange(lens.max())
    out = np.zeros(mask.shape,dtype=np.float64)
    out[mask] = np.concatenate(array)
    
    return out

################################################################################


def compress_CSVs(p, mouseids=all_mice, paradigms=all_parad, stim_types=all_stimtypes):
    """Reads in all the CSVs produced and saves them as gzip binary. The 32(pos)
    + 32(neg) spike timestamp CSVs are merged into a pos and neg gzip`ed frame.
    """
    path = p['outputPath']
    
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