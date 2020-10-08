import numpy as np
from MUA_init import initialize

# get the user input 
P = initialize()

# all the conditions of the experiment
ALL_MICE = 'mGE82', 'mGE83', 'mGE84', 'mGE85'
ALL_PARADIGMS = 'DAC1', 'DAC2', 'O10C1', 'O10C2', 'O25C1', 'O25C2', 'O25UC1', 'O25UC2', 'MS', 'DOC1', 'DOC2'
ALL_STIMTYPES = 'Standard', 'Predeviant', 'Deviant', 'Postdeviant', 'D1', 'C2', 'C1', 'B1', 'UniquePredeviant', 'UniquePostdeviant'
PARADIGMS_STIMTYPES = {'DAC1': ['Deviant',],
                       'DAC2': ['Deviant',],
                       'O10C1': ['Standard', 'Predeviant', 'Deviant', 'Postdeviant'],
                       'O10C2': ['Standard', 'Predeviant', 'Deviant', 'Postdeviant'],
                       'O25C1': ['Standard', 'UniquePredeviant', 'Deviant', 'UniquePostdeviant'],
                       'O25C2': ['Standard', 'UniquePredeviant', 'Deviant', 'UniquePostdeviant'],
                       'O25UC1': ['Standard', 'Predeviant', 'Deviant', 'Postdeviant'],
                       'O25UC2': ['Standard', 'Predeviant', 'Deviant', 'Postdeviant'],
                       'MS': ['D1', 'C2', 'C1', 'B1'],
                       'DOC1': ['Standard', 'Deviant',],
                       'DOC2': ['Standard', 'Deviant',]}
MICE_DATES = {
    'mGE82': 'mGE82_24.07.2019',
    'mGE83': 'mGE83_29.07.2019',
    'mGE84': 'mGE84_30.07.2019',
    'mGE85': 'mGE85_31.07.2019',
}
# order of paradigms in experiment - not applicable paradigms were removed
PARAD_ORDER =  {'mGE82': ('O10C1',	'DAC1',	    'DAC2',	    'O10C2',	'O25UC1',	'O25C1',	'MS',	    'O25C2',	'O25UC2',	'DOC1',	'DOC2'),
                'mGE83': ('O25UC2',	'O25C2',	'MS',	    'O25C1',	'O25UC1',	'O10C2',	'DAC2',	    'DAC1',	    'O10C1',	'DOC2',	'DOC1'),
                'mGE84': ('O10C2',	'DAC2',	    'DAC1',	    'O10C1',	'O25UC2',	'O25C2',	'MS',	    'O25C1',	'O25UC1',	'DOC2',	'DOC1'),
                'mGE85': ('O10C1',	'DAC1',	    'DAC2',	    'O10C2',	'O25UC1',	'O25C1',	'MS',	    'O25C2',	'O25UC2',	'DOC1',	'DOC2')
}

PARAD_PAIRS = (('DAC1', 'DAC2'),
               ('DOC1', 'DOC2'),
               ('O10C1', 'O10C2'),
               ('O25C1', 'O25C2'),
               ('O25UC1', 'O25UC2'))


# specify the filetype used to save plots
PLOT_FORMAT = 'svg'
# PLOT_FORMAT = 'png'

# ------------------GWENDOLYN NEW DATA---------------
ALL_MICE = ['mGE33', 'mGE35', 'mGE36', 'mGE47', 'mGE48', 'mGE49', 'mGE50', 
            'mGE51', 'mGE52', 'mGE53', 'mGE54', 'mGE57', 'mGE58', 'mGE71', 
            'mGE73', 'mGE74', 'mGE76', 'mGE77', 'mGE79', 'mGE80', 
            # ]
            'mGE82', 'mGE83', 'mGE84', 'mGE85']
MICE_DATES = {
 'mGE33': 'mGE33_04.02.2019',
 'mGE35': 'mGE35_12.02.2019',
 'mGE36': 'mGE36_28.02.2019',
 'mGE47': 'mGE47_02.04.2019',
 'mGE48': 'mGE48_04.04.2019',
 'mGE49': 'mGE49_08.04.2019',
 'mGE50': 'mGE50_09.04.2019',
 'mGE51': 'mGE51_12.04.2019',
 'mGE52': 'mGE52_17.04.2019',
 'mGE53': 'mGE53_18.04.2019',
 'mGE54': 'mGE54_23.04.2019',
 'mGE57': 'mGE57_09.05.2019',
 'mGE58': 'mGE58_10.05.2019',
 'mGE71': 'mGE71_27.06.2019',
 'mGE73': 'mGE73_01.07.2019',
 'mGE74': 'mGE74_02.07.2019',
 'mGE76': 'mGE76_04.07.2019',
 'mGE77': 'mGE77_05.07.2019',
 'mGE79': 'mGE79_09.07.2019',
 'mGE80': 'mGE80_10.07.2019',
 
 'mGE82': 'mGE82_24.07.2019',
 'mGE83': 'mGE83_29.07.2019',
 'mGE84': 'mGE84_30.07.2019',
 'mGE85': 'mGE85_31.07.2019',
}


# More readable paradigm string for plot annoation
PARAD_FULL = {'DAC1': 'Deviant alone C1',
              'DAC2': 'Deviant alone C2',
              'DOC1': 'Deviant omission C1',
              'DOC2': 'Deviant omission C2',
              'MS': 'Many Standards',
              'O10C1': 'Oddball 10% C1',
              'O10C2': 'Oddball 10% C2',
              'O25C1': 'Oddball 25% C1',
              'O25C2': 'Oddball 25% C2',
              'O25UC1': 'Oddball Unif. 25% C1',
              'O25UC2': 'Oddball Unif. 25% C2'}

PROJ_DIR = '/media/loaloa/Samsung_T5/gdrive/projects/ephys'
LFP_OUTPUT = P['outputPath'] + '/../../output/LFP_output'

# set to something like a 100000 to never delete artifact channels
# importantly, this value is used to compare against the negative firingrate but 
# then slice both the negative AND positive trials. The positive firingrate
# generally follows the pattern of the negative firingrate artifacts (checked visually)
# but is usually a bit lower. If you work with the positive firingrates consider 
# classifying positive and negative seperately.
ARTIFACT_TRIAL_COV_THR = 100
ARTIFACT_TRIAL_COV_HM_MIN = 0
ARTIFACT_TRIAL_COV_HM_MAX = 200
SI_MIN_FRATE_5MS = .5

# OPtion for processing the raw data. If True, all plots are renamed to the 
# physical channel location based on P["id"]. So 1 .. 32 would go from dorsal
# to ventral for standard 32 channel electrode. 
CHNL_TO_PHYSICAL_ORDER = False

# predefined colors to use for labeling 
COLORS = {'red':       '#e6194B',
          'deep_red':  '#8a0b25',
          'green':     '#3cb44b',
          'yellow':    '#ffe119',
          'orange':    '#f58231',
          'deep_blue': '#0a3b70',
          'cyan':      '#42d4f4',
          'purple':    '#911eb4',
          'pink':      '#fabebe',
          'lavender':  '#e6beff',
          'beige':     '#fffac8',
          'blue':      '#4363d8',
          'apricot':   '#ffd8b1',
          'lime':      '#bfef45',
          'magenta':   '#f032e6',
          'teal':      '#469990',
          'mint':      '#aaffc3',
          'olive':     '#808000',
          'brown':     '#9A6324',
          'grey':      '#a9a9a9',
          'white':     '#ffffff',
          'black':     '#000000',
          'light yellow':   '#FFFF80',
          'deep green':     '#005C31',
          'neon yellow':    '#FFFF00',
          'neon red':       '#FF0010',
          'deep orange':    '#FF5005',
}

REGION_CMAP  = {'not_assigned': COLORS['white'], 
                'SG': COLORS['green'], 
                'G': COLORS['deep_blue'], 
                'IG': COLORS['orange'], 
                'dIG': COLORS['red'],
                'Th': COLORS['teal'],
                'SG_lateSI': COLORS['green'], 
                'G_lateSI': COLORS['deep_blue'], 
                'IG_lateSI': COLORS['orange'], 
                'dIG_lateSI': COLORS['red'],
                'Th_lateSI': COLORS['teal'],
                }

REGIONS = {
           'VPM': 'VPM',
           'POM': 'POM',
           'Th': 'VPM',
           'G': 'granular',
           'G_mid': 'granular',
           'SG': 'supra-granular',
           'SG_mid': 'supra-granular',
           'IG': 'infra-granular',
           'IG_mid': 'infra-granular',
           'dIG': 'deep infra-granular',
           'dIG_mid': 'deep infra-granular',
            1: 'Channel 1',
            2: 'Channel 2',
            3: 'Channel 3',
            3: 'Channel 3',
            4: 'Channel 4',
            5: 'Channel 5',
            6: 'Channel 6',
            7: 'Channel 7',
            8: 'Channel 8',
            9: 'Channel 9',
            10: 'Channel 10',
            11: 'Channel 11',
            12: 'Channel 12',
            13: 'Channel 13',
            14: 'Channel 14',
            15: 'Channel 15',
            16: 'Channel 16',
            17: 'Channel 17',
            18: 'Channel 18',
            19: 'Channel 19',
            20: 'Channel 20',
            21: 'Channel 21',
            22: 'Channel 22',
            23: 'Channel 23',
            24: 'Channel 24',
            25: 'Channel 25',
            26: 'Channel 26',
            27: 'Channel 27',
            28: 'Channel 28',
            29: 'Channel 29',
            30: 'Channel 30',
            31: 'Channel 31',
            32: 'Channel 32',
            }


REGIONS_EXT = {
               'Th': 'VPM',
               'G': 'granular',
               'SG': 'supra-granular',
               'IG': 'infra-granular',
               'dIG': 'deep infra-granular',
               'Th_lateSI': 'VPM_lateSI',
               'G_lateSI': 'granular_lateSI',
               'SG_lateSI': 'supra-granular_lateSI',
               'IG_lateSI': 'infra-granular_lateSI',
               'dIG_lateSI': 'deep infra-granular_lateSI',
}

GENERAL_CMAP = {

    'mGE82': COLORS['deep green'],
    'mGE83': COLORS['neon yellow'],
    'mGE84': COLORS['neon red'],
    'mGE85': COLORS['deep orange'],

    'mGE33': COLORS['blue'],
    'mGE35': COLORS['green'],
    'mGE36': COLORS['yellow'],
    'mGE47': COLORS['orange'],
    'mGE48': COLORS['deep_blue'],
    'mGE49': COLORS['cyan'],
    'mGE50': COLORS['purple'],
    'mGE51': COLORS['pink'],
    'mGE52': COLORS['lavender'],
    'mGE53': COLORS['beige'],
    'mGE54': COLORS['white'],
    'mGE57': COLORS['apricot'],
    'mGE58': COLORS['lime'],
    'mGE71': COLORS['magenta'],
    'mGE73': COLORS['teal'],
    'mGE74': COLORS['mint'],
    'mGE76': COLORS['olive'],
    'mGE77': COLORS['brown'],
    'mGE79': COLORS['grey'],
    'mGE80': COLORS['black'],

    'DAC1': COLORS['deep_red'], 
    'DAC2': COLORS['deep_red']+'44', 
    'O10C1': COLORS['deep_blue'], 
    'O10C2': COLORS['deep_blue']+'44', 
    'O25C1': COLORS['yellow'], 
    'O25C2': COLORS['yellow']+'44', 
    'O25UC1': COLORS['green'], 
    'O25UC2': COLORS['green']+'44', 
    'MS': COLORS['orange'], 
    'DOC1': COLORS['purple'], 
    'DOC2': COLORS['purple']+'44', 

    'Standard': COLORS['teal'],
    'Predeviant': COLORS['mint'],
    'UniquePredeviant': COLORS['mint'],
    'Deviant': COLORS['magenta'],
    'Postdeviant': COLORS['lime'],
    'UniquePostdeviant': COLORS['lime'],
    'C1': COLORS['grey'],
    'C2': COLORS['grey'],
    'D1': COLORS['white'],
    'B1': COLORS['white'],
    
    'SG': COLORS['green'], 
    'SG_mid': COLORS['green']+'44', 
    'G': COLORS['deep_blue'], 
    'G_mid': COLORS['deep_blue']+'44', 
    'IG': COLORS['orange'], 
    'IG_mid': COLORS['orange']+'44', 
    'dIG': COLORS['red'],
    'dIG_mid': COLORS['red']+'44',
    'Th': COLORS['teal'],
    'VPM': COLORS['teal'],
    'POM': COLORS['mint'],
    'not_assigned': COLORS['white'],

    1: COLORS['grey'],
    2: COLORS['grey'],
    3: COLORS['grey'],
    3: COLORS['grey'],
    4: COLORS['grey'],
    5: COLORS['grey'],
    6: COLORS['grey'],
    7: COLORS['grey'],
    8: COLORS['grey'],
    9: COLORS['grey'],
    10: COLORS['grey'],
    11: COLORS['grey'],
    12: COLORS['grey'],
    13: COLORS['grey'],
    14: COLORS['grey'],
    15: COLORS['grey'],
    16: COLORS['grey'],
    17: COLORS['grey'],
    18: COLORS['grey'],
    19: COLORS['grey'],
    20: COLORS['grey'],
    21: COLORS['grey'],
    22: COLORS['grey'],
    23: COLORS['grey'],
    24: COLORS['grey'],
    25: COLORS['grey'],
    26: COLORS['grey'],
    27: COLORS['grey'],
    28: COLORS['grey'],
    29: COLORS['grey'],
    30: COLORS['grey'],
    31: COLORS['grey'],
    32: COLORS['grey'],
}