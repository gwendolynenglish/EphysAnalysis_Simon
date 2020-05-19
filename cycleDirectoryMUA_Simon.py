##########################################################################################################
# Gwendolyn English 20.04.2020
#
# Function cycles through all files for analysis and outputs two primary analysis files:
#    1. Evoked Responses for All Stimuli in paradigm 
#           a. CSV & pickle file, and PNG image produced for each electrode channel  
#    2. Average Evoked Responses for specific Trigger subset
#           a. CSV & pickle file, containing all channels, produced for each trigger subset  
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
from mainAnalysisMUA_Simon import *  
##########################################################################################################

############################################################################################################
#Input to primary function: dictionary pickle file created in LocalFieldPotentialEvaluation jupyter notebook 

def MUA_analyzeAllFiles(p):
#####identify and cycle all folders in directory
    dirs = os.listdir(p['inputPath'])
        
    for folder in dirs:
        print(folder)
#########Create folder in output path for each raw file 
        outputpathFolder = p['outputPath'] + '/' + folder
        if not os.path.exists(outputpathFolder): 
            os.mkdir(outputpathFolder)
            
#########Cycle through all files in folder 
        for file in os.listdir(p['inputPath'] + '/' + folder):
        
############################################################################################################                
#############Identify all additional trigger files   
            if 'Triggers' in file and 'AllStimuli' not in file: 
                print(file)  
            
                #Initialize array for trigger summary data                 
                summary_data_to_file = []
                firing_rates_over_time = [] 
            
                #load trigger file
                trigger_array = readTrigger(p['inputPath'] + '/' + folder + '/' + file)

#################Cycle through channels
                for channel_file in os.listdir(p['inputPath'] + '/' + folder):
                    if 'ElectrodeChannel' in channel_file:

                        #Load channel array 
                        electrodechannel = channel_file[:-4]
                        channel_array = readChannelRawData(p['inputPath'] + '/' + folder + '/' + channel_file)
                        
                        #Complete Analysis
                        summary_data, firing_rates = triggers(trigger_array, channel_array, outputpathFolder, file, channel_file, p) 
    
                        #Append summary data        
                        summary_data_to_file.append(summary_data)
                        firing_rates_over_time.append(firing_rates)
            
#################Restructure summary data and write to file 
                #Summary Data 
                headers = (['ElectrodeChannel', 'NumTrials', 'NumNonresponsiveTrials', 'AvgTimetoFirstSpike', \
                            'TotalSpikesPostTrigger', 'AvgSpikesPostTrigger', 'NumSpikesSlopeOverTrials', \
                            'Gamma mu (0-end ms)', 'Gamma std (0-end ms)', 'Gamma skew (0-end ms)', \
                            'Normal mu (0-end ms)', 'Normal std (0-end ms)',
                            'Gamma mu (0-50ms)', 'Gamma std (0-50ms)', 'Gamma skew (0-50ms)', \
                            'Normal mu (0-50ms)', 'Normal std (0-50ms)']) 
                

                channelmap = np.asarray(p['id'] - 1).flatten()
                datatofile_summary = np.asarray(summary_data_to_file)                
                datatofile_summary = datatofile_summary[channelmap, :]
                datatofile = datatofile_summary.reshape(p['nr_of_electrodes'], 17)
                pd.DataFrame(datatofile).to_csv(outputpathFolder + '/' + file[:-4] + '_SpikeSummary.csv', header = headers)
                
                #Firing Rates 
                headers = np.arange(-p['evoked_pre'] *1000, p['evoked_post'] * 1000, p['psth_binsize'])
                datatofile_firingrates = np.asarray(firing_rates_over_time)
                datatofile_firingrates = datatofile_firingrates[channelmap, :]
                pd.DataFrame(datatofile_firingrates).to_csv(outputpathFolder + '/' + file[:-4] + '_FiringRates.csv', header = headers)
                
                
                