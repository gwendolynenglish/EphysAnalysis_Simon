import numpy as np
from MUA_init import initialize

# get the user input 
P = initialize()

# all the conditions of the experiment
ALL_MICE = 'mGE82', 'mGE83', 'mGE84', 'mGE85'
ALL_PARADIGMS = 'DAC1', 'DAC2', 'DOC1', 'DOC2', 'MS', 'O10C1', 'O10C2', 'O25C1', 'O25C2', 'O25UC1', 'O25UC2'
ALL_STIMTYPES = 'Standard', 'Predeviant', 'Deviant', 'Postdeviant', 'D1', 'C2', 'C1', 'B1', # 'UniquePredeviant', 'UniquePostdeviant'

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

# set to something like a 100000 to never delete artifact channels
# importantly, this value is used to compare against the negative firingrate but 
# then slice both the negative AND positive trials. The positive firingrate
# generally follows the pattern of the negative firingrate artifacts (checked visually)
# but is usually a bit lower. If you work with the positive firingrates consider 
# classifying positive and negative seperately.
ARTIFACT_TRIAL_COV_THR = 9

# predefined colors to use for labeling 
colors = ['#e6194B', #    0 = red
          '#3cb44b', #    1 = green
          '#ffe119', #    2 = yellow
          '#4363d8', #    3 = blue
          '#f58231', #    4 = orange
          '#911eb4', #    5 = purple
          '#42d4f4', #    6 = cyan
          '#f032e6', #    7 = magenta
          '#bfef45', #    8 = lime
          '#fabebe', #    9 = pink
          '#469990', #    10 = teal
          '#e6beff', #    11 = lavender
          '#9A6324', #    12 = brown
          '#fffac8', #    13 = beige
          '#8a0b25', #    14 = deep red
          '#aaffc3', #    15 = mint
          '#808000', #    16 = olive
          '#ffd8b1', #    17 = apricot
          '#0a3b70', #    18 = deep blue
          '#a9a9a9', #    19 = grey
          '#ffffff', #    20 = white
          '#000000'  #    21 = black
]