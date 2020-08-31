# coding: utf-8
import MUA_constants as const
from preprocessing import compress_CSVs, process_data
import plotting

# for this to work, change the input and output path in MUA_init.py
# then, change ALL_MICE, MICE_DATES, GENERAL_CMAP, (PARD_ORDER) in 
# MUA_constants.py to the new mice. This should do it.


"""Process data """
# process_data()
compress_CSVs()

"""Onset offset"""
train_lbls_tsv = f'{const.P["outputPath"]}/../../output_lowthr/onset_offset_labels.tsv'
train_features_csv = f'{const.P["outputPath"]}/../../output_lowthr/onset_offset_spikebins_channels.csv'
plotting.classify_onset_offset(dest_dir_appdx='../onset_offset_classification/',
                               train_lbls_tsv=train_lbls_tsv, 
                               train_features_csv=train_features_csv)


