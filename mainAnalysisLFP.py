"""
Gwendolyn English 03.09.2019

Functions return LFP analysis called from cycleDirectory.py 
"""

##########################################################################################################
#Import required packages & functions
import pickle 
import pandas as pd 
import matplotlib 
# matplotlib.use('agg')
from matplotlib.pyplot import *

from loadFiles import * 
from preprocessing import *
from plotting import * 
from analysisFunctionsLFP import *
##########################################################################################################
def allStimuli(p, dataPath, trigger_array, channel_file, folderPath, outputpathFolder):
    """
    Inputs: Path to stimulation file.
    Output: CSV & Pickle file containing evoked responses & statistics info for each stimuli for each individual electrode channel.
    """

    #Load electrode channel raw file 
    channel_array = readChannelRawData(dataPath)
    
    #Low pass filter signal
    filtered_channel_array = lowpass_filter(data = channel_array, samplerate = p['sample_rate'], \
                                            cutoffFreq = p['low_pass_freq'], order = 3)

    #Trigger align signals  
    stim_timestamps = extract_stim_timestamps(trigger_array, p)
    aligned_filtered_channel_array =align_to_trigger(filtered_channel_array,stim_timestamps, p)
    
    #Downsample signals
    ds_aligned_filtered_channel_array = down_sample_2D(aligned_filtered_channel_array, p['sample_rate'])
    
    #Norm all evoked responses to 0V at stimulus onset time
    normed_ds_aligned_filtered_channel_array = norm_to_zero(ds_aligned_filtered_channel_array, p)

    #Downsampled time course
    time = np.linspace(-p['evoked_pre'], p['evoked_post'], int(p['sample_rate']*(p['evoked_pre'] +p['evoked_post'])))
    ds_time = down_sample_1D(time, p['sample_rate'])
   
    #Extract peak negative values & time stamps
    peak_negs = np.amin(normed_ds_aligned_filtered_channel_array, axis = 1)
    peak_negs_ts = ds_time[np.argmin(normed_ds_aligned_filtered_channel_array, axis = 1)] * 1000
    
    #Read Labview output
    whisker_stim = loadLabviewOutput(folderPath, outputpathFolder, stim_timestamps)
    
    #Write data to file 
    datatofile = np.vstack((peak_negs * 1000000, peak_negs_ts, whisker_stim, np.transpose(normed_ds_aligned_filtered_channel_array)))
    datatofile = np.transpose(datatofile)
    headers = np.hstack((['Peak_negs_uV', 'Peak_negs_ts_ms', 'Whisker'], ds_time))
    pd.DataFrame(datatofile).to_csv(outputpathFolder + '/' + channel_file[:-4] + '.csv', header = headers)
  
    pickle.dump({'Peak_negs_uV':datatofile[:,0], 'Peak_negs_ts_ms': datatofile[:,1], 'Whisker': datatofile[:,4], \
                 'ERP':datatofile[:,5:]}, 
                 open(outputpathFolder + '/' + channel_file[:-4] + '.pickle' , 'wb'), protocol = -1)
    
    #Plot all evoked responses onto one figure 
    #fig = matplotlib.pyplot.figure()
    #for stim in range(len(trigger_array)):
    #    plot(ds_time,normed_ds_aligned_filtered_channel_array[stim, :] *1000000)
    #savefig(outputpathFolder + '/' + channel_file[:-4] + '.png', format = 'png')
    #close(fig)
    
    #Clear variables 
    del peak_negs, peak_negs_ts

##########################################################################################################    
def avgStimuli(p, avg_dataPath, avg_trigger_array, outputpathFolder, triggerFile, channelFile): 
    """
    Input: Path to stimulation file. 
    Output: CSV & Pickle file containing averaged evoked responses & statistics for all electrode channels 
    """   
    
    #Load electrode channel raw file
    avg_channel_array = readChannelRawData(avg_dataPath)
    
    #Low pass filter signal
    avg_filtered_channel_array = lowpass_filter(data = avg_channel_array, samplerate = p['sample_rate'], \
                                                cutoffFreq = p['low_pass_freq'], order = 3)

    #Trigger align signals
    avg_stim_timestamps = extract_stim_timestamps(avg_trigger_array, p)
    avg_aligned_filtered_channel_array =align_to_trigger(avg_filtered_channel_array,avg_stim_timestamps, p)
    
    #Downsample signals
    avg_ds_aligned_filtered_channel_array = down_sample_2D(avg_aligned_filtered_channel_array, p['sample_rate']) 
    
    #Norm all evoked responses to 0mV at stimulus onset time
    avg_normed_ds_aligned_filtered_channel_array = norm_to_zero(avg_ds_aligned_filtered_channel_array, p)
    
    #Downsampled time course
    time = np.linspace(-p['evoked_pre'], p['evoked_post'], int(p['sample_rate']*(p['evoked_pre'] +p['evoked_post'])))
    ds_time = down_sample_1D(time, p['sample_rate'])
    
    #Extract peak negative values & time stamps - identify average and standard deviation of each. Add to output array. 
    peak_negs = np.amin(avg_normed_ds_aligned_filtered_channel_array, axis = 1)
    peak_negs_ts = ds_time[np.argmin(avg_normed_ds_aligned_filtered_channel_array, axis = 1)] * 1000
    peak_negs_avg = np.mean(peak_negs) 
    peak_negs_sd = np.std(peak_negs)
    peak_negs_ts_avg = np.mean(peak_negs_ts)
    peak_negs_ts_sd = np.std(peak_negs_ts) 
    
    #Average all evoked responses
    avg_evoked = []
    avg_evoked = np.mean(avg_normed_ds_aligned_filtered_channel_array, axis = 0)
    peak_neg = np.amin(avg_evoked)
    peak_neg_ts = ds_time[np.argmin(avg_evoked)] * 1000
    
    summary_data = []
    summary_data = np.hstack((peak_neg * 1000000, peak_neg_ts, peak_negs_avg * 1000000, peak_negs_sd * 1000000, \
                              peak_negs_ts_avg, peak_negs_ts_sd, avg_evoked))
    
    #Plot average evoked response
    plot_evoked_channel(avg_evoked, outputpathFolder, triggerFile, channelFile)
    
    #Compute wavelet transform for each trial and average coefficients
    avg_coeff, freq = wavelet_cycle_trials(avg_normed_ds_aligned_filtered_channel_array, p)
    
    #Plot wavelet transform
    plot_wavelet_heatmap(avg_coeff, freq, outputpathFolder, triggerFile, channelFile, 'trial', p)
      
    #Save wavelet info 
    headers = freq
    pd.DataFrame(avg_coeff.T).to_csv(outputpathFolder + '/WaveletCoefficients_' + \
                                     triggerFile[:-4] + '_' + channelFile[:-4] + '.csv', header = headers)
    pickle.dump({'Wavelet_coeff': avg_coeff, 'Wavelet_freq' : freq}, 
                 open(outputpathFolder + '/WaveletCoefficients_' + \
                                     triggerFile[:-4] + '_' + channelFile[:-4] +  '.pickle' , 'wb'), protocol = -1)

    #Clear variables
    del peak_negs, peak_negs_ts, peak_negs_avg, peak_negs_sd, peak_negs_ts_avg, peak_negs_ts_sd, \
        avg_evoked, avg_channel_array, avg_stim_timestamps, avg_aligned_filtered_channel_array, \
        avg_ds_aligned_filtered_channel_array, avg_normed_ds_aligned_filtered_channel_array 
    
    return summary_data
    