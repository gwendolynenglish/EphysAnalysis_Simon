"""
Gwendolyn English 20.04.2020

Functions direct MUA analysis called from cycleDirectoryMUA.py 
"""

################################################################################
#Import required packages & functions
import pickle 
import pandas as pd
import scipy.stats as ss
import matplotlib 
# matplotlib.use('agg')

from loadFiles import * 
from preprocessing import *
from plotting import * 
from MUA_functions import * 


################################################################################
    
def triggers(trigger_array, channel_array, outputpathFolder, trigger_filename, 
             channel_filename):
    """
    Inputs: Trigger array, Channel array, Output folder path, Parameter dictionary 
    Outputs: Plots: PSTH, PSTH with Gamma Fit, Raster Plots, Firing Rate Plots 
    """
        
    # #Prepare channel array data and extract relevant information 
    # Preprocess Data 
    neg_crossings, pos_crossings, aligned_hpf_data = preprocessMUA(channel_array, 
                                                               trigger_array) 

    #Extract all spike timestamps 
    #Extracts all pre- and post- stimulus neg spike timestamps
    neg_timestamps_withLabel = extract_ts(neg_crossings)    
    #Remove stimulus label column     
    neg_timestamps = neg_timestamps_withLabel[:,1:]            
  
    #Extracts all pre- and post- stimulus pos spike timestamps
    pos_timestamps_withLabel = extract_ts(pos_crossings)    
    # Remove stimulus label column 
    # pos_timestamps = pos_timestamps_withLabel[:,1:]
    
    #Extract all spike waveforms 
    neg_waveforms = extract_wf(aligned_hpf_data, neg_crossings)
    # pos_waveforms = extract_wf(aligned_hpf_data, pos_crossings)
    
    #Extract firing rates
    # neg_firingrate = firing_rate(neg_crossings)    #deprecated
    # pos_firingrate = firing_rate(pos_crossings)    #deprecated
    neg_firingrate = firing_rate(neg_timestamps)
    # pos_firingrate = firing_rate(pos_timestamps)
    

    #Plotting    
    # #Plot Peri-Stimulus-Time-Histograms  
    outputpath = outputpathFolder + '/' + 'PSTH_NegativeSpikes_' + \
                 trigger_filename[:-4] + '_' + channel_filename[:-4] + '.png'
    plot_PSTH(neg_timestamps, outputpath)
    
    #Plot Raster plots
    outputpath = outputpathFolder + '/' + 'Raster_NegativeSpikes_' \
                 + trigger_filename[:-4] + '_' + channel_filename[:-4] + '.png'
    plot_raster(neg_timestamps, outputpath)
    
    #Plot waveforms
    outputpath = outputpathFolder + '/' + 'Waveforms_NegativeSpikes_' + \
                 trigger_filename[:-4] + '_' + channel_filename[:-4] + '.png'
    plot_waveforms(neg_waveforms, outputpath)
 
    #Plot firing rates
    outputpath = outputpathFolder + '/' + 'FiringRate_NegativeSpikes_' + \
                 trigger_filename[:-4] + '_' + channel_filename[:-4] + '.png'
    plot_firing_rate(neg_firingrate, outputpath) 
    # return None, None


    #Collect Summary Data -- UPDATE HEADERS TO REFLECT NEW ENTRIES
    summary_data = []
    
    #Create another matrix instance to remove pre-stim spikes 
    neg_timestamps_postStim = neg_timestamps
    #Exchange negative times for 0
    neg_timestamps_postStim[neg_timestamps_postStim < 0] = 0   
    
    totTrials = len(trigger_array)
    #Sum all trials with no spike timestamps after stimulus
    totUnresponsiveTrials = np.sum(~neg_timestamps_postStim.any(1))  
    totSpikespostTrigger = np.count_nonzero(neg_timestamps_postStim)
    avgSpikespostTrigger = totSpikespostTrigger/totTrials 
    
    #Average time to first spike 
    if np.sum(neg_timestamps_postStim) == 0: 
        avgTimetoFirstSpike = 0 
    if np.sum(neg_timestamps_postStim) != 0:
        mins = np.amin(neg_timestamps_postStim, axis = 1)
        avgTimetoFirstSpike = mins[np.nonzero(mins)].mean()
    
    #Linear Regression over spikes evoked by trial throughout paradigm 
    y = np.count_nonzero(neg_timestamps_postStim, axis = 1)
    x = np.arange(len(y))
    r = ss.linregress(x,y)
    slope_over_trials = r.slope 

    #Gamma distribution fit - accounts for only post-stimulus spikes
    neg_timestamps_postStim_clipped = neg_timestamps[neg_timestamps>0] 
    k, alpha, theta = ss.gamma.fit(neg_timestamps_postStim_clipped, floc = 0)
    gamma_mu = k * theta
    gamma_std = np.sqrt(k * theta**2)
    gamma_skew =  2 / np.sqrt(k)

    #Gaussian distribution fit - accounts for only post-stimulus spikes
    neg_timestamps_postStim_clipped = neg_timestamps[neg_timestamps>0]
    mu, sigma = ss.norm.fit(neg_timestamps_postStim_clipped)
    normal_mu = mu
    normal_std = sigma
    
    #Create another matrix instance to isolate spikes in 0-50ms post stim window 
    neg_timestamps_0to50ms = neg_timestamps_postStim          
    #Exchange spike times >50ms with 0  
    neg_timestamps_0to50ms[neg_timestamps_0to50ms > 50] = 0    
    
    #Gamma distribution fit - accounts for only post-stimulus spikes
    neg_timestamps_0to50ms_clipped = neg_timestamps_postStim_clipped[ \
                                     neg_timestamps_postStim_clipped<=50] 
    k_50, alpha_50, theta_50 = ss.gamma.fit(neg_timestamps_0to50ms_clipped, 
                                            floc = 0)
    gamma_mu_50 = k_50 * theta_50
    gamma_std_50 = np.sqrt(k_50 * theta_50**2)
    gamma_skew_50 =  2 / np.sqrt(k_50)
    
    #Gaussian distribution fit - accounts for only post-stimulus spikes
    neg_timestamps_0to50ms_clipped = neg_timestamps_postStim_clipped[ \
                                     neg_timestamps_postStim_clipped<=50] 
    mu_50, sigma_50 = ss.norm.fit(neg_timestamps_0to50ms_clipped)
    normal_mu_50 = mu_50
    normal_std_50 = sigma_50   

    summary_data.append([channel_filename, totTrials, totUnresponsiveTrials, \
                         avgTimetoFirstSpike, totSpikespostTrigger, \
                         avgSpikespostTrigger, slope_over_trials, \
                         gamma_mu, gamma_std, gamma_skew, normal_mu, \
                         normal_std, gamma_mu_50, gamma_std_50, gamma_skew_50, \
                         normal_mu_50, normal_std_50])
        

    #Save Data to File, return summary data  
    avg_firing_rate_over_time = neg_firingrate
    fname = outputpathFolder + '/' + trigger_filename[:-4] + '_' + \
            channel_filename[:-4] + '_TS_negSpikes.csv'
    pd.DataFrame(neg_timestamps_withLabel).to_csv(fname)
    fname = outputpathFolder + '/' + trigger_filename[:-4] + '_' + \
            channel_filename[:-4] + '_TS_posSpikes.csv'
    pd.DataFrame(pos_timestamps_withLabel).to_csv(fname)
    
    return summary_data, avg_firing_rate_over_time  