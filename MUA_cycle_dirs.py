################################################################################
# Gwendolyn English 20.04.2020
#
# Function cycles through all files for analysis and outputs two primary 
# analysis files:
#    1. Evoked Responses for All Stimuli in paradigm 
#           a. CSV & pickle file, and PNG image produced for each electrode 
#              channel  
#    2. Average Evoked Responses for specific Trigger subset
#           a. CSV & pickle file, containing all channels, produced for each 
#              trigger subset  
################################################################################

#Import required packages & functions
import sys
import os
import numpy as np
import scipy.stats as ss
import pandas as pd  

import concurrent.futures

import MUA_constants as const
from MUA_core import * 
from plotting import covariance_artifact_heatmap 


# redo_dict = {
#     'mGE84_30.07.2019_O25C1.mcd': ['Deviant'],
    # 'mGE82_24.07.2019_DOC1.mcd': ['Deviant', 'Predeviant'],
    # 'mGE85_31.07.2019_DAC2.mcd': ['Deviant', 'Standard'],
    # 'mGE84_30.07.2019_DAC2.mcd': ['Deviant', 'Standard'],

# }


################################################################################

#Input to primary function: dictionary pickle file created in 
# LocalFieldPotentialEvaluation jupyter notebook 

#####identify and cycle all folders in directory
def MUA_analyzeMouseParadigm(folder):
    print('\n\n\n\n', folder, '\n')

    # Create folder in output path for each raw file 
    outputpathFolder = const.P['outputPath'] + '/' + folder
    if not os.path.exists(outputpathFolder): 
        os.mkdir(outputpathFolder)
        
    # Cycle through all files in folder 
    for file in os.listdir(const.P['inputPath'] + '/' + folder):
    
################################################################################                
        # Identify all additional trigger files   
        if 'Triggers' in file and 'AllStimuli' not in file: 
            stim_t = file[file.rfind('_')+1:-4]
            print(file, '  -  Processing 32 channels now:')  
            
            # if stim_t not in ['Deviant']:
            #     continue
        
            
            #load trigger file
            fname = const.P['inputPath'] + '/' + folder + '/' + file
            trigger_array = readTrigger(fname)
            
            def iterate_channels(artifact_trials=None):
                #Initialize array for trigger summary data                 
                summary_data_to_file = []
                firing_rates_over_time = [] 
                all_trial_frates = []
                
                # Cycle through channels
                for channel_file in os.listdir(const.P['inputPath']+'/'+folder):
                    if 'ElectrodeChannel' in channel_file:
                        electrodechannel = int(channel_file[-6:-4])
                        print(electrodechannel, end='..')
                        # if not electrodechannel in [1,2]:
                        #     continue

                        #Load channel array 
                        fname = const.P['inputPath'] +'/'+ folder+'/'+channel_file
                        channel_array = readChannelRawData(fname)
                        
                        #Complete Analysis
                        summary_data, firing_rates, trial_frates = triggers(trigger_array, 
                                                                channel_array, 
                                                                outputpathFolder, 
                                                                file, 
                                                                channel_file,
                                                                artifact_trials,
                                                                folder) 
                                                            
    
                        #Append summary data        
                        summary_data_to_file.append(summary_data)
                        firing_rates_over_time.append(firing_rates)
                        all_trial_frates.append(trial_frates)
                print('. Done.')
                return summary_data_to_file, firing_rates_over_time, all_trial_frates
            
            # get all firingrates of all channels
            out = iterate_channels()
            summary_data_to_file, firing_rates_over_time, all_trial_frates = out
            
            # compute faulty trials, possibly rerun iter_channels passing them
            cov = np.cov(pd.concat(all_trial_frates, axis=1))
            fault_trials = (cov > const.ARTIFACT_TRIAL_COV_THR).any(axis=0)
            covariance_artifact_heatmap(cov, fault_trials, f'{outputpathFolder}/trial_covariance_{stim_t}.png')
            if fault_trials.any():
                print(f'Trial covariance analysis found {fault_trials.sum()}'
                        f'/{len(fault_trials)} trials to be artifacts. Rerunning '
                        'processing with these trials set to NaN.')
                out = iterate_channels(fault_trials)
                summary_data_to_file, firing_rates_over_time, all_trial_frates = out

################################################################################                
            # Restructure summary data and write to file 
            #Summary Data 
            headers = (['ElectrodeChannel', 'NumTrials', 'NumNonresponsiveTrials', 
                        'AvgTimetoFirstSpike', 'TotalSpikesPostTrigger', 
                        'AvgSpikesPostTrigger', 'NumSpikesSlopeOverTrials', 
                        'Gamma mu (0-end ms)', 'Gamma std (0-end ms)', 
                        'Gamma skew (0-end ms)', 'Normal mu (0-end ms)', 
                        'Normal std (0-end ms)', 'Gamma mu (0-50ms)', 
                        'Gamma std (0-50ms)', 'Gamma skew (0-50ms)', 
                        'Normal mu (0-50ms)', 'Normal std (0-50ms)']) 
            
            channelmap = np.asarray(const.P['id'] - 1).flatten()
            datatofile_summary = np.asarray(summary_data_to_file)
            datatofile_summary = datatofile_summary[channelmap, :]
            datatofile = datatofile_summary.reshape(const.P['nr_of_electrodes'], 17)
            fname = outputpathFolder + '/' + file[:-4] + '_SpikeSummary.csv'
            pd.DataFrame(datatofile).to_csv(fname, header=headers)
            
            #Firing Rates 
            frm, to = -const.P['evoked_pre']*1000, const.P['evoked_post']*1000
            headers = np.arange(frm, to, const.P['psth_binsize'])   
            datatofile_firingrates = np.asarray(firing_rates_over_time)
            datatofile_firingrates = datatofile_firingrates[channelmap, :]
            fname = outputpathFolder + '/' + file[:-4] + '_FiringRates.csv'
            pd.DataFrame(datatofile_firingrates).to_csv(fname, 
                                                        header=headers)
            print('FiringRates and Summary stats saved.\n')