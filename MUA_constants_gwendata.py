import numpy as np
from MUA_init import initialize

# get the user input 
P = initialize()

MUA_output = P["outputPath"]
PROJ_DIR = '/media/loaloa/Samsung_T5/gdrive/projects/ephys'

# all the conditions of the experiment
ALL_MICE = 'mGE22', 'mGE70'
ALL_PARADIGMS = 'DOC1',
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

LFP_OUTPUT = P['outputPath'] + '/../../output/LFP_output'

# set to something like a 100000 to never delete artifact channels
# importantly, this value is used to compare against the negative firingrate but 
# then slice both the negative AND positive trials. The positive firingrate
# generally follows the pattern of the negative firingrate artifacts (checked visually)
# but is usually a bit lower. If you work with the positive firingrates consider 
# classifying positive and negative seperately.
ARTIFACT_TRIAL_COV_THR = 75
ARTIFACT_TRIAL_COV_HM_MIN = 0
ARTIFACT_TRIAL_COV_HM_MAX = 200
SI_MIN_FRATE_5MS = .5

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
           'Th': 'VPM',
           'G': 'granular',
           'SG': 'supra-granular',
           'IG': 'infra-granular',
           'dIG': 'deep infra-granular'}

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
    'mGE82': COLORS['red'],
    'mGE83': COLORS['green'],
    'mGE84': COLORS['yellow'],
    'mGE85': COLORS['blue'],

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
    'G': COLORS['deep_blue'], 
    'IG': COLORS['orange'], 
    'dIG': COLORS['red'],
    'Th': COLORS['teal'],
}