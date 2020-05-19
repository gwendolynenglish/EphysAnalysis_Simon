"""
This file contains functions for filtering, downsampling, and aligning the data. 

Created on 02.11.18. 
"""

from scipy import signal
import numpy as np
from scipy.signal import decimate 

#################################################################################################
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

#################################################################################################
#Alignment Functions

def extract_stim_timestamps(stim_array, p):
    
    #Transform trigger times from seconds to samples 
    stim_timestamps = np.round(stim_array * p['sample_rate']) 

    return stim_timestamps 

def align_to_trigger(filtered_data, stim_timestamps, p):
    
    #Create holder array with dimension rows: # of stimulations, columns: window length 
    evoked_array = np.zeros((len(stim_timestamps), int(p['sample_rate']*(p['evoked_pre']+p['evoked_post']))))
    for stim in range(len(stim_timestamps)):
        evoked_array[stim,:] = np.array(filtered_data[int(stim_timestamps[stim] - p['evoked_pre']*p['sample_rate']):int(stim_timestamps[stim] + p['evoked_post']*p['sample_rate'])])
    return evoked_array 
                             

#################################################################################################
#Function for making list of uneven-lengthed arrays into matrix (0-filled values)     

def numpy_fillna(array):
    lens = np.array([len(item) for item in array])
    mask = lens[:,None] > np.arange(lens.max())
    out = np.zeros(mask.shape,dtype=np.float64)
    out[mask] = np.concatenate(array)
    
    return out
    
    
    
    
    
    
    
    
    