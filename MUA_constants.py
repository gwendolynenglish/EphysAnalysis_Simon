import numpy as np
from MUA_init import initialize

# get the user input 
P = initialize()

MUA_output = P["outputPath"]
PROJ_DIR = '/home/loaloa/gdrive/projects/ephys'

# all the conditions of the experiment
ALL_MICE = 'mGE82', 'mGE83', 'mGE84', 'mGE85'
ALL_PARADIGMS = 'DAC1', 'DAC2', 'DOC1', 'DOC2', 'MS', 'O10C1', 'O10C2', 'O25C1', 'O25C2', 'O25UC1', 'O25UC2'
ALL_STIMTYPES = 'Standard', 'Predeviant', 'Deviant', 'Postdeviant', 'D1', 'C2', 'C1', 'B1', 'UniquePredeviant', 'UniquePostdeviant'

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
              'O25UC1': 'Oddball Uniform 25% C1',
              'O25UC2': 'Oddball Uniform 25% C2'}

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

LFP_OUTPUT = P['outputPath'] + '/../LFP_output'

# set to something like a 100000 to never delete artifact channels
# importantly, this value is used to compare against the negative firingrate but 
# then slice both the negative AND positive trials. The positive firingrate
# generally follows the pattern of the negative firingrate artifacts (checked visually)
# but is usually a bit lower. If you work with the positive firingrates consider 
# classifying positive and negative seperately.
ARTIFACT_TRIAL_COV_THR = 25
ARTIFACT_TRIAL_COV_HM_MIN = 0
ARTIFACT_TRIAL_COV_HM_MAX = 200

# predefined colors to use for labeling 
COLORS = {'red':       '#e6194B',
          'green':     '#3cb44b',
          'yellow':    '#ffe119',
          'blue':      '#4363d8',
          'orange':    '#f58231',
          'purple':    '#911eb4',
          'cyan':      '#42d4f4',
          'magenta':   '#f032e6',
          'lime':      '#bfef45',
          'pink':      '#fabebe',
          'teal':      '#469990',
          'lavender':  '#e6beff',
          'brown':     '#9A6324',
          'beige':     '#fffac8',
          'deep_red':  '#8a0b25',
          'mint':      '#aaffc3',
          'olive':     '#808000',
          'apricot':   '#ffd8b1',
          'deep_blue': '#0a3b70',
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
                }