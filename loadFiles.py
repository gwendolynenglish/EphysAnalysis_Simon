"""
Loads different file types for use in analysis 

Created on 03.09.2019
"""

#Import required packages 
import pickle
import os
import sys
import numpy as np
import pandas as pd

def readChannelRawData(filepath):
    """
    This function reads the data from a .dat file in binary float64 format.
    
    Inputs:
        filepath: The path to the .dat file to be read
        
    Outputs:
        amplifier_file: numpy array containing the data from the .dat file (in uV)
    """
    
    with open(filepath, 'rb') as fid:
        channel_array = np.fromfile(fid, np.float64)
    return channel_array 

def readTrigger(filepath):
    """
    This function reads the data from a .dat file in binary float64 format.
    
    Inputs:
        filepath: The path to the .dat file to be read
        
    Outputs:
        amplifier_file: numpy array containing the data from the .dat file (in uV)
    """
    
    with open(filepath, 'rb') as fid:
        trigger_array = np.fromfile(fid, np.float64)
    return trigger_array

def loadLabviewOutput(folderPath, outputpathFolder):
    LVfilepath = ''
    for LVfile in os.listdir(folderPath):
        if 'StimSeq' in LVfile:
            LVfilepath = folderPath + '/' + LVfile 
    if LVfilepath != '':   
        whisker_stim = readLabviewOuput(LVfilepath, outputpathFolder)
    else: whisker_stim = np.zeros(len(stim_timestamps))
    
    return whisker_stim 

def readLabviewOuput(filepath, paradigm):
    """ 
    This function reads the Labview output folder that identifies stimulated whiskers.
  
    Inputs:
        filepath: The path to the .txt file to be read
        
    Outputs:
        whisker_stim_array: array containing order of whisker stimulations 
    """    

    with open(filepath) as stimseq:
        whiskers = np.loadtxt(stimseq)
        
    #if 'PSC' in paradigm:
    #    whiskers = np.asarray(whiskers)
    #    whiskers = whiskers.reshape((40,100))
    #    whiskers = whiskers[:,0:72]
    #    whiskers = whiskers.flatten()
    
    #Identify C2 stimulations
    C2 = ['C2' if x==1 else x for x in whiskers]#PSC1
    #C2 = ['C2' if x==0 else x for x in whiskers]#PSC2
    
    #Identify C1 stimulations
    C2C1 = ['C1' if x==0 else x for x in C2]#PSC1
    #C2C1 = ['C1' if x==1 else x for x in C2]#PSC2
    
    #Identify D1 stimulations
    C2C1D1 = ['D1' if x==2 else x for x in C2C1]
    
    #Identify B1 stimulations
    C2C1D1B1 = ['B1' if x==3 else x for x in C2C1D1] 

    return C2C1D1B1 

def loadCSDAssignment(p, filepath):
    """ 
    This function reads the csv file created by the user which assigns electrode channels to different cortical regions. 
    
    Inputs:
        filepath: The path to the .csv file to be read
        
    Outputs:
        electrode_assignment: array containing order of electrodes along cortical depth
    """ 
   
    datafromfile = pd.read_csv(filepath)
    if p['nr_of_shanks'] == 1:
        shank1electrodes = datafromfile.shank1electrodes
        shank1assignments = datafromfile.shank1assignments 
        shank2electrodes = [None]
        shank2assignments = [None]
         
    if p['nr_of_shanks'] == 2:
        shank1electrodes = datafromfile.shank1electrodes
        shank1assignments = datafromfile.shank1assignments
        shank2electrodes = datafromfile.shank2electrodes
        shank2assignments = datafromfile.shank2assignments 
    
    if p['nr_of_shanks'] ==2: return shank1electrodes, shank1assignments, shank2electrodes, shank2assignments
    if p['nr_of_shanks'] ==1: return shank1electrodes, shank1assignments
    
    
         
def headerID(p, paradigm):
    
    if p['evoked_pre'] == 0.05 and p['evoked_post'] == 0.25 and p['binsize'] == 0.05 and paradigm == 'PSC' :
        headers = np.vstack((['Electrode', 'C1_PSC_-50-0ms', 'C1_PSC_0-50ms', 'C1_PSC_50-100ms', 'C1_PSC_100-150ms', \
                                  'C1_PSC_150-200ms', 'C1_PSC_200-250ms', 'C2_PSC_-50-0ms', 'C2_PSC_0-50ms', 'C2_PSC_50-100ms', \
                                  'C2_PSC_100-150ms', 'C2_PSC_150-200ms', 'C2_PSC_200-250ms', 'B1_PSC_-50-0ms', 'B1_PSC_0-50ms', \
                                  'B1_PSC_50-100ms', 'B1_PSC_100-150ms', 'B1_PSC_150-200ms', 'B1_PSC_200-250ms', \
                                  'D1_PSC_-50-0ms', 'D1_PSC_0-50ms', 'D1_PSC_50-100ms', 'D1_PSC_100-150ms', 'D1_PSC_150-200ms', \
                                  'D1_PSC_200-250ms', 'C1_PSCC_-50-0ms', 'C1_PSCC_0-50ms', 'C1_PSCC_50-100ms', 'C1_PSCC_100-150ms', \
                                  'C1_PSCC_150-200ms', 'C1_PSCC_200-250ms', 'C2_PSCC_-50-0ms','C2_PSCC_0-50ms', 'C2_PSCC_50-100ms', \
                                  'C2_PSCC_100-150ms', 'C2_PSCC_150-200ms', 'C2_PSC_200-250ms', 'B1_PSCC_-50-0ms', 'B1_PSCC_0-50ms', \
                                  'B1_PSCC_50-100ms', 'B1_PSCC_100-150ms', 'B1_PSCC_150-200ms', 'B1_PSC_200-250ms', \
                                  'D1_PSCC_-50-0ms', 'D1_PSCC_0-50ms', 'D1_PSCC_50-100ms', 'D1_PSCC_100-150ms', \
                                  'D1_PSCC_150-200ms', 'D1_PSC_200-250ms']))
                
    if p['evoked_pre'] == 0.05 and p['evoked_post'] == 0.25 and p['binsize'] == 0.05 and paradigm == 'CompPSC':
        headers = np.vstack((['Electrode', 'LayerAssignment', 'C1_-50-0ms', 'C1_0-50ms', 'C1_50-100ms', 'C1_100-150ms', \
                             'C1_150-200ms', 'C1_200-250ms',\
                             'C2_-50-0ms', 'C2_0-50ms', 'C2_50-100ms', 'C2_100-150ms', 'C2_150-200ms', 'C2_200-250ms',\
                             'B1_-50-0ms', 'B1_0-50ms', 'B1_50-100ms', 'B1_100-150ms', 'B1_150-200ms', 'B1_200-250ms',\
                             'D1-50-0ms', 'D1_0-50ms', 'D1_50-100ms', 'D1_100-150ms', 'D1_150-200ms', 'D1_200-250ms']))
                
    if p['evoked_pre'] == 0.05 and p['evoked_post'] == 0.20 and p['binsize'] == 0.05 and paradigm == 'Spikes':
        headers = (['ElectrodeChannel', 'NegSpikes_-50-0ms', 'NegSpikes_0-50ms', 'NegSpikes_50-100ms', \
                                'NegSpikes_100-150ms', 'NegSpikes_150-200ms', 'PosSpikes_-50-0ms', \
                                'PosSpikes_0-50ms', 'PosSpikes_50-100ms', 'PosSpikes_100-150ms', \
                                'PosSpikes_150-200ms'])  
    if p['evoked_pre'] == 0.01 and p['evoked_post'] == 0.05 and p['binsize'] == 0.01 and paradigm == 'Spikes':    
        headers = (['ElectrodeChannel', 'NegSpikes_-10-0ms', 'NegSpikes_0-10ms', 'NegSpikes_10-20ms', \
                                'NegSpikes_20-30ms', 'NegSpikes_30-40ms', 'NegSpikes_40-50ms', 'PosSpikes_-10-0ms', \
                                'PosSpikes_0-10ms', 'PosSpikes_10-20ms', 'PosSpikes_20-30ms', \
                                'PosSpikes_30-40ms', 'PosSpikes_40-50ms']).flatten() 
         
    return headers 

