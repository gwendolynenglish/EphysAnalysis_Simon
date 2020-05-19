##########################################################################################################
#Import required packages & functions
import sys
import os
import numpy as np

from loadFiles import * 
from mainAnalysisLFP import * 
from analysisFunctionsLFP import * 
from plotting import * 
##########################################################################################################

def analyzeAllFiles(p):
    
#####Identify and cycle all folders in directory
    dirs = os.listdir(p['inputPath'])    
    for folder in dirs:
        
#########Create folder in output path for each raw file 
        outputpathFolder = p['outputPath'] + '/' + folder
        if not os.path.exists(outputpathFolder): 
            os.mkdir(outputpathFolder)
            
#########Cycle through all files in folder 
        for file in os.listdir(p['inputPath'] + '/' + folder):
        
#############Identify AllStimuli file 
            if 'AllStimuli' in file:
        
                print(p['inputPath'] + '/' + folder + '/' + file)  
                triggerArray = readTrigger(p['inputPath'] + '/' + folder + '/' + file)
                
                for channelFile in os.listdir(p['inputPath'] + '/' + folder):    
                    if 'ElectrodeChannel' in channelFile:
                        
                        dataPath = p['inputPath'] + '/' + folder + '/' + channelFile
                        folderPath = p['inputPath'] + '/' + folder
                
                        #Analyze the response to each stimulus for entire paradigm 
                        allStimuli(p, dataPath, triggerArray, channelFile, folderPath, outputpathFolder) 
       
 ############Identify all other trigger files
            elif 'Triggers' in file:
                
                print(p['inputPath'] + '/' + folder + '/' + file) 
                triggerArray = readTrigger(p['inputPath'] + '/' + folder + '/' + file)
                
                #Initialize data holder 
                avgdatatofile = []
                
                for channelFile in os.listdir(p['inputPath'] + '/' + folder):    
                    if 'ElectrodeChannel' in channelFile:
                        
                        dataPath = p['inputPath'] + '/' + folder + '/' + channelFile
                        
                        summarydata = avgStimuli(p, dataPath, triggerArray, outputpathFolder, file, channelFile)
                        avgdatatofile.append(summarydata)
                        
                
                #Write summary data to pickle file 
                avgdatatofile = np.asarray(avgdatatofile)
                channelorder = p['id']
                channelorder = channelorder.flatten() - 1 
                avgdata = avgdatatofile[channelorder, :]
                pickle.dump({'Electrode': p['id'], 'Peak_neg_uV': avgdata[:,0], 'Peak_neg_ts': avgdata[:,1], \
                             'Avg_Peak_neg_uV': avgdata[:,2], 'SD_Peak_negs_uV': avgdata[:,3], \
                             'Avg_peak_neg_TS_ms': avgdata[:,4], 'SD_peak_neg_TS_ms': avgdata[:,5], \
                             'ERP':avgdata[:,10:]}, \
                              open(outputpathFolder + '/' + file[:-4] + '.pickle' , 'wb'), protocol = -1)
                                          
                #Write summary data to csv
                #Downsampled time course
                time = np.linspace(-p['evoked_pre'], p['evoked_post'], p['sample_rate']*(p['evoked_pre'] +p['evoked_post']))
                ds_time = down_sample_1D(time, p['sample_rate']) 
                avgdata = np.hstack((p['id'].reshape((p['nr_of_electrodes'], 1)), avgdata))
                headers = np.hstack((['Electrode', 'Peak_neg_uV', 'Peak_neg_ts', 'Avg_Peak_neg_uV', 'SD_Peak_negs_uV', \
                                      'Avg_peak_neg_TS_ms',  'SD_peak_neg_TS_ms'], ds_time))
                pd.DataFrame(avgdata).to_csv(outputpathFolder + '/' + file[:-4] + '_LFPAverages.csv', header = headers)                     
                
                #Plot evoked shank
                avg_evoked = avgdatatofile[:,6:] #remove all summary data
                channelmap = p['id']
                
                for shank in range(p['nr_of_shanks']):
                    shankmap = channelmap[shank] - 1 #returns indices, not channel numbers
                    avg_evoked_shank = avg_evoked[shankmap, :]
                    plot_evoked_shank(avg_evoked_shank, outputpathFolder, file, ds_time, shank, shankmap)
                
                
#################CSD Analysis for each shank, if dictated
                    if p['csd'] == True:
                    
                        CSD = delta_iCSD(avg_evoked_shank, p)
                        pd.DataFrame(CSD).to_csv(outputpathFolder + '/CSD_shank_' +  str(shank) + '_' + file[:-4] + '.csv')
                    
                        plot_CSD(CSD, outputpathFolder, file, ds_time, shank, shankmap)
                        plot_CSD_heatmap(CSD, outputpathFolder, file, ds_time, shank, shankmap)
                    

                    
                
                
                
                
                