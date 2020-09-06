# coding: utf-8
import MUA_constants as const
from preprocessing import compress_CSVs, process_data, fix_channelmap1x32_to2x16, adjust_channelmap
import plotting

# for this to work, change the input and output path in MUA_init.py
# then, change ALL_MICE, MICE_DATES, GENERAL_CMAP, (PARD_ORDER) in 
# MUA_constants.py to the new mice. This should do it.


"""Process data """
# process_data(how={'nbatches': 8, 'batch': 1},)
# process_data()
# fix_channelmap1x32_to2x16()
# compress_CSVs()
# adjust_channelmap()

"""Onset offset"""
train_lbls_tsv = f'{const.P["outputPath"]}/../../output_lowthr/onset_offset_labels.tsv'
train_features_csv = f'{const.P["outputPath"]}/../../output_lowthr/onset_offset_spikebins_channels.csv'
# for rank in ['mouse', 'paradigm', 'stimulus_type', 'channel']:
#     plotting.classify_onset_offset(dest_dir_appdx='../onset_offset_classification',
#                                 train_lbls_tsv=train_lbls_tsv, print_pos_rasters_pred=False,
#                                 train_features_csv=train_features_csv, plot_labeled_data=True,
#                                 rank=rank)

plotting.classify_onset_offset(dest_dir_appdx='../onset_offset_classification',
                            train_lbls_tsv=train_lbls_tsv, print_pos_rasters_pred=False,
                            train_features_csv=train_features_csv, plot_labeled_data=True,
                            cluster=True
)