"""
Uploaded on GitHub on Tuesday, Aug 1st, 2017

author: Tansel Baran Yasar

Contains the information and the maps for the Neuronexus probes that are used in the lab.
"""

import numpy as np
import pickle

def load_probe_info(probe):
    print('loadding')
    """
    This function generates a dictionary containing the information about the probe used in the experiment.

    Inputs:
        probe: String indicating the model of the probe.

    Outputs:
        probe_info: Dictionary including the information about the probe. The dictionary must include the following keys:
            'numShanks' : Number of shanks
            'type': Configuration of the electrodes (tetrode, linear, polytrode)
            'numTetrodes': Total number of tetrodes on the probe (only for tetrode configuration)
            'numTetrodesPerShank': Number of tetrodes on each shank (only for tetrode configuration)
            'numTrodesPerShank': Number of electrodes on each shank (only for linear configuration)
            'numTrodes': Total number of electrodes on the probe (only for linear configuration)
            'id': For tetrode configuration, (numTetrodesPerShank) x (numShanks) x 4; for linear   configuration, (numShanks) x (numTrodesPershank) Numpy array
        containing how the physical mapping of the probe corresponds to the mapping of the electrodes on Intan. When generating this list, please use the following
        convention: Shanks are numbered from 0, starting from the left. Tetrodes or the electrodes of a linear probe are numbered from 0, starting from top of the
        shank. Electrodes in a tetrode are numbered from 0, starting from the left-most electrode and continuing counter clockwise.
    """
    
    probe_info = {}
    probe_info['name'] = probe

      
    if probe == 'a2x16_10mm_100_500_177':
        probe_info['numShanks'] = 2
        probe_info['type'] = 'linear'
        probe_info['numTrodesPerShank'] = 16
        probe_info['numTrodes'] = 32

        probe_info['bottom_ycoord'] = 50
        probe_info['top_ycoord'] = 1550    
        id = np.array(([18,27,29,19,17,25,31,20,26,24,28,21,30,23,32,22],[8,9,11,4,10,6,13,2,12,7,15,5,14,3,16,1]))
        probe_info['id'] = id 
        probe_info['sitespacing'] = 100
    
    elif probe == 'a4x8_5mm_200_400_177':
        probe_info['numShanks'] = 4
        probe_info['type'] = 'linear'
        probe_info['numTrodesPerShank'] = 8
        probe_info['numTrodes'] = 32
        
        probe_info['bottom_ycoord'] = 50
        probe_info['top_ycoord'] = 1450
        probe_info['sitespacing'] = 200
        id = np.array(([20,24,25,21,19,23,27,22],[26,31,28,17,30,29,32,18],[2,7,6,5,4,3,9,1],[12,13,15,10,14,11,16,8]))
        probe_info['id'] = id 
    
    elif probe == 'a1x32_6mm_100_177':
        probe_info['numShanks'] = 1
        probe_info['type'] = 'linear'
        probe_info['numTrodesPerShank'] = 32
        probe_info['numTrodes'] = 32
        
        probe_info['bottom_ycoord'] = 0
        probe_info['top_ycoord'] = 3200
        probe_info['sitespacing'] = 100
        id = np.array(([1,32,3,30,5,28,7,26,2,31,6,17,4,29,9,18,8,27,11,19,10,25,13,20,12,24,15,21,14,23,16,22],))
        probe_info['id'] = id
        
    elif probe == 'a2x16_10mm_50_500_177':
        probe_info['numShanks'] = 2
        probe_info['type'] = 'linear'
        probe_info['numTrodesPerShank'] = 16
        probe_info['numTrodes'] = 32
        
        probe_info['bottom_ycoord'] = 0
        probe_info['top_ycoord'] = 800
        probe_info['sitespacing'] = 50
        id = np.array(([18,27,29,19,17,25,31,20,26,24,28,21,30,23,32,22],[8,9,11,4,10,6,13,2,12,7,15,5,14,3,16,1]))
        probe_info['id'] = id

    return probe_info

