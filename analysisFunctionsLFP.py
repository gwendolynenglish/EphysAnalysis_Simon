"""
Gwendolyn English 03.09.2019

Functions return LFP analysis called from mainAnalysisLFP.py 
"""

##########################################################################################################
# Gwendolyn English 27.04.2020
#
# Functions created April 2020 for better modularizing code.   
##########################################################################################################

##########################################################################################################
#Import required packages & functions
import sys
import os
import numpy as np
import pandas as pd 
import pywt

from loadFiles import * 
from preprocessing import *
from plotting import * 
##########################################################################################################

##########################################################################################################
#Downsampling Function for 2D arrays
#Inputs: Data array & sample rate
#Returns: Data array downsampled yo 1kHz 
        
def down_sample_2D(data, sample_rate):
    
    #Downsample signal to 1kHZ
    ds_rate = int(sample_rate / 1000)
    ds_signal_array = data[:,::ds_rate]
    
    return ds_signal_array 

##########################################################################################################
#Downsampling Function for 1D arrays
#Inputs: Data array & sample rate
#Returns: Data array downsampled yo 1kHz 
        
def down_sample_1D(data, sample_rate):
    
    #Downsample signal to 1kHZ
    ds_rate = int(sample_rate / 1000)
    ds_signal_array = data[::ds_rate]
    
    return ds_signal_array 


##########################################################################################################
#Normalize downsampled data to zero
#Inputs: Stimulus-aligned data matrix (rows: trials, columns: ERP across analysis window), parameter dictionary
#Returns: Stimulus-aligned data matrix, where value of ERP are shifted to 0V at stimulus onset 

def norm_to_zero(data, p):
    
    stim_onset_time = int(p['sample_rate']*p['evoked_pre'] / int(p['sample_rate'] / 1000)) 
    stim_onset_value = data[:,stim_onset_time]        
    
    subtract_mat = np.tile(stim_onset_value, (int((int(p['sample_rate']*p['evoked_pre']) + int(p['sample_rate']*p['evoked_post']))/int(p['sample_rate'] / 1000)), 1)).transpose() 
    normed_data = data - subtract_mat
    
    return normed_data

##########################################################################################################
#Delta inverted-Current Source Density Analysis described in Pettersen et al 2006. 
#Inputs: Matrix of averaged LFPs (Rows: Channels organized as on shank, Columns: Averaged evoked response), parameter dictionary
#Returns: CSD Matrix 

#Function to compute the delta-iCSD
def delta_iCSD(avg_LFPs, p):

    #Set Parameters
    radius = 0.5						#radius of assumed constant planar current source density (in mm) 
    conductivity = 0.3						#extracellular electrical conductivity in S/m
    numElec = np.shape(avg_LFPs)[0]				#number of electrode sites on shank
    siteSpacing = p['sitespacing'] / 1000000		#distance in between electrode sites in meters 
    tres = 	1/ p['sample_rate']			#time resolution of samples (1/sampling frequency)

    #Array containing Spatial Position of Electrodes 
    z = np.linspace(siteSpacing, siteSpacing * numElec, numElec) 

    for column in avg_LFPs.T:
 
        #Generates/clears the F matrix for estimating the CSD 
        F_matrix = np.zeros((numElec, numElec)) # holder matrix

        for i in np.arange(numElec):
            for j in np.arange(numElec): 
                F_matrix[j, i] = (siteSpacing/(2 * conductivity)) * (np.sqrt(np.square(z[j] - z[i]) + np.square(radius)) - abs(z[j] - z[i]))


    #Calculate CSD using the inverse of the F matrix 
    CSD = np.matmul(np.linalg.inv(F_matrix), avg_LFPs)
    CSD = CSD / (1000000e3)   #converts CSD data units to um 
    CSD = CSD * 1000000       #converts CSD data units to uA 
    return CSD

##########################################################################################################
#Continuous 'Morlet' Wavelet Transformation 
#Inputs: Signal
#Returns: Frequencies (Scales) and Coefficients Matrix 

def wavelet_transform(signal, p):
    
    #Assumes that the down-sampled signal at 1kHz is analyzed
    sampling_period = 1 / 1000
   
    #Use built in pywt function to determine physical frequencies 
    scales = pywt.scale2frequency( 'morl', np.arange(1,120)) / sampling_period
   
    #Complete wavelet anlaysis 
    coef, freq = pywt.cwt(signal, scales, 'morl', sampling_period) 
    
    return coef, freq 

##########################################################################################################
#Compute trial by trial wavelet 
#Inputs: 2D matrix of trials and evoked responses 
#Returns: averaged wavelet 

def wavelet_cycle_trials(data, p):

    #Holder array for wavelet coefficiants
    wavelet_coeff = [] 
    
    #Compute wavelet coefficients for each trial 
    for trial in np.arange(np.shape(data)[0]):
        coef, freq = wavelet_transform(data[trial], p)
        wavelet_coeff.append(coef) 
    
    #Average wavelet coefficients across trials 
    avg_coeff = np.mean(wavelet_coeff, axis=0)    

    return avg_coeff, freq

##########################################################################################################
#Compute trial by trial wavelet 
#Inputs: 2D matrix of trials and evoked responses 
#Returns: averaged wavelet 

def wavelet_cycle_trials_remove_avg(avg, data, p):

    #Holder array for wavelet coefficiants
    wavelet_coeff = [] 
    
    #Compute wavelet coefficients for averaged evoked response over all trials 
    avg_coef, avg_freq = wavelet_transform(avg, p)
    
    #Compute wavelet coefficients for each trial and subtract average 
    for trial in np.arange(np.shape(data)[0]):
        coef, freq = wavelet_transform(data[trial], p)
        coef_minusEvoked = coef - avg_coef
        wavelet_coeff.append(coef_minusEvoked) 
    
    #Average wavelet coefficients across trials 
    avg_coeff = np.mean(wavelet_coeff, axis=0)    

    return avg_coeff, freq
    
    