from MUA_init import initialize

ALL_MICE = 'mGE82', 'mGE83', 'mGE84', 'mGE85'
ALL_PARADIGMS = 'DAC1', 'DAC2', 'DOC1', 'DOC2', 'MS', 'O10C1', 'O10C2', 'O25C1', 'O25C2', 'O25UC1', 'O25UC2'
ALL_STIMTYPES = 'Standard', 'Predeviant', 'Deviant', 'Postdeviant', 'D1', 'C2', 'C1', 'B1'

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

P = initialize()