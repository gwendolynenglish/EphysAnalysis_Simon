##########################################################################################################
# Gwendolyn English 15.04.2020
#
# Functions created April 2020 for better modularizing code.   
##########################################################################################################

##########################################################################################################
#Import required packages & functions
import sys
import os
import numpy as np
import pandas as pd  
from loadFiles import * 
from preprocessing import *
from plotting import * 
##########################################################################################################

##########################################################################################################
#Preprocess MUA data 
#Inputs: data array, array of stimulus triggers, parameter dictionary
#Returns: data array of size #stim x window - 1 for positive and negative threshold crossings and aligned data 

def preprocessMUA(data, trigger_array,  p):
    
    #high-pass filter data
    hp_data = highpass_filter(data, p['sample_rate'], p['high_pass_freq'], order = 4)  
    
    #extract stimulus time stamps and align data within specified windows 
    stim_ts = extract_stim_timestamps(trigger_array, p)
    aligned_hp_data = align_to_trigger(hp_data, stim_ts, p)
    
    #determine the negative and positive thresholds
    neg_thresh = -p['threshold'] * np.std(hp_data) 
    pos_thresh = p['threshold'] * np.std(hp_data)
    
    #artifact detection - eliminiate 
    
    #identify all instances of threshold responses
    neg_abovethresh = np.diff((aligned_hp_data < neg_thresh).astype(np.int))
    pos_abovethresh = np.diff((aligned_hp_data > pos_thresh).astype(np.int))
    
    #Identify all threshold crossings
    neg_x = np.where(neg_abovethresh == 1, neg_abovethresh, 0)
    pos_x = np.where(pos_abovethresh == 1, pos_abovethresh, 0)
    
    return neg_x, pos_x, aligned_hp_data  

##########################################################################################################
#Extract Time Stamps
#Inputs: Matrix of threshold Crossing, parameter dictionary  
#Returns: set of arrays, each array contains the spike timestamps of one stimluation  

def extract_ts(data_x, p): 
    
    #compute timing array in ms, where the first element is removed to align with threshold-crossings array 
    time = np.linspace(-p['sample_rate'] * p['evoked_pre'], p['sample_rate'] * p['evoked_post'],\
                                           p['sample_rate'] * p['evoked_pre'] + p['sample_rate'] * p['evoked_post'] + 1)
    time = time / p['sample_rate'] * 1000
    time = np.delete(time, 0)
    
    #compute threshold-crossing time stamps
    ts_byStim = []
    
    for stimulus in range(np.shape(data_x)[0]):
        indices = np.where(data_x[stimulus]==1)
        ts = np.take(time, indices)
        
        if ts.size != 0:
            timestamps = np.insert(ts, 0, stimulus)
        if ts.size == 0:
            timestamps= ([stimulus])
            
        ts_byStim.append(timestamps)
    
    #Reformat timestamp array 
    ts_byStim = np.asarray(ts_byStim)
    ts_byStim = numpy_fillna(ts_byStim)
    
    return ts_byStim 

##########################################################################################################
#Extract and plot waveforms
#Inputs: Array of binned highpass-filtered data, matrix of threshold crossings, parameter dictionary
#Returns: 2D Array containing all detected waveforms  

def extract_wf(data, data_x,  p):
    
    #Identify samples to collect before and after threshold crossing (2ms)
    wf_window = p['sample_rate']/1000 * 2 
    
    #Initialize waveform holder
    waveforms = [] 
    for stimulus in range(np.shape(data_x)[0]):
        indices = np.where(data_x[stimulus]==1)
        indices = np.asarray(indices).flatten()
        for crossing in range(len(indices)):
            wf_before = int(indices[crossing] - wf_window)
            wf_after = int(indices[crossing] + wf_window)
            waveform = data[stimulus, wf_before:wf_after]
            
            waveforms.append(np.asarray(waveform))
    
    return waveforms 
    
##########################################################################################################
#Calculate smoothed firing rate
#Input: Matrix of threshold crossings, parameter dictionary 
#Returns: 1D array containing the average firing rate per bin

def firing_rate(data_x, p):
    #Compensate for removal of first column of data crossings during preprocessMUA 
    zeros = np.zeros((np.shape(data_x)[0], 1))
    data_x_fulllen = np.hstack((zeros, data_x))
    
    #Calculate time bins 
    bins = np.arange(-p['evoked_pre'] *1000, p['evoked_post'] * 1000, p['psth_binsize'])

    #Bin data into corresponding time bins, sum activity within bins, then average activity across trials 
    binned_data_x = data_x_fulllen.reshape((np.shape(data_x_fulllen)[0], len(bins), int(np.shape(data_x_fulllen)[1]/len(bins))))
    summed_binned_data_x =  binned_data_x.sum(axis = 2)
    avg_summed_binned_data_x = summed_binned_data_x.mean(axis = 0) 
        
    #Calculate firing rate according to bin size, dividing by the bin size in milliseconds 
    binned_firing_rate = avg_summed_binned_data_x / (p['psth_binsize']/1000) 
        
    return binned_firing_rate 
        

