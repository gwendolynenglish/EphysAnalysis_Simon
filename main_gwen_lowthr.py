# coding: utf-8
import MUA_constants as const
from preprocessing import compress_CSVs, process_data
import plotting
import os
from MUA_constants import P


# change MUA_constans ALL_MICE as well!
if not P["outputPath"].endswith('output_gwen_data_lowthr/MUA_output'):
    print('Outputpath wrong for this script? Check MUA_init.')
    exit()


"""Process data """
# process_data()
# compress_CSVs()


"""Onset offset"""
from onset_offset_classif import onset_offset_response
from onset_offset_classif import onset_offset_labels
from onset_offset_classif import lapl_kernel_SVM
from onset_offset_classif import get_onset_offset_classification
from onset_offset_classif import onoff_heatmap
from onset_offset_classif import onoff_barplot
from onset_offset_classif import idk_yet

training_data_dir = '/mnt/Samsung_T5/output_lowthr/onset_offset_model'
training_data_chnl_map_file = '/mnt/Samsung_T5/output_lowthr/chnls_map.csv'
MUA_output_data_chnl_map_file = '/mnt/Samsung_T5/output_gwen_data_lowthr/channel_mappings/S1Th_LayerAssignment_22.10.19.csv'


# dest_dir_appdx = '../pred_gwendata_lowthr'
# dest_dir_plot_appdx = '../pred_gwendata_lowthr/gwendata_plots'
# os.makedirs(f'{P["outputPath"]}/{dest_dir_plot_appdx}', exist_ok=True)

# new_data = get_onset_offset_classification(which_data='MUA_output', training_data_dir=training_data_dir, dest_dir_appdx=dest_dir_appdx, cached_prediction=True,
#                                        training_data_chnl_map_file=training_data_chnl_map_file, MUA_output_data_chnl_map_file=MUA_output_data_chnl_map_file)
# onoff_heatmap(new_data, dest_dir_appdx=dest_dir_plot_appdx, fig_height=14)
# new_data_all = get_onset_offset_classification(which_data='MUA_output', training_data_dir=training_data_dir, dest_dir_appdx=dest_dir_appdx, keep_labels=[0,1,2,3],
#                                                  training_data_chnl_map_file=training_data_chnl_map_file, MUA_output_data_chnl_map_file=MUA_output_data_chnl_map_file)
# onoff_barplot(new_data_all, dest_dir_appdx=dest_dir_plot_appdx)


# both training and new data
# dest_dir_appdx = '../pred_gwendata_lowthr'
# dest_dir_plot_appdx = '../pred_gwendata_lowthr/gwendata_plus_train_plots'
# os.makedirs(f'{P["outputPath"]}/{dest_dir_plot_appdx}', exist_ok=True)

# both_data = get_onset_offset_classification(which_data='both', training_data_dir=training_data_dir, dest_dir_appdx=dest_dir_appdx, cached_prediction=True,
#                                        training_data_chnl_map_file=training_data_chnl_map_file, MUA_output_data_chnl_map_file=MUA_output_data_chnl_map_file)
# onoff_heatmap(both_data, dest_dir_appdx=dest_dir_plot_appdx, fig_height=15)

# both_data = get_onset_offset_classification(which_data='both', training_data_dir=training_data_dir, dest_dir_appdx=dest_dir_appdx, keep_labels=[0,1,2,3],
#                                                  training_data_chnl_map_file=training_data_chnl_map_file, MUA_output_data_chnl_map_file=MUA_output_data_chnl_map_file)
# onoff_barplot(both_data, dest_dir_appdx=dest_dir_plot_appdx)

# dest_dir_plot_appdx = '../pred_gwendata_lowthrgwendata_plus_train_plots_fastonoff'
# os.makedirs(f'{P["outputPath"]}/{dest_dir_plot_appdx}', exist_ok=True)
# onoff_barplot(both_data, dest_dir_appdx=dest_dir_plot_appdx, keep_labels=[3])









# _200_Raster_NegativeSpikes_Triggers_Deviant_ElectrodeChannel_21.png
# _231_Raster_NegativeSpikes_Triggers_C1_ElectrodeChannel_21.png
# _342_Raster_NegativeSpikes_Triggers_Deviant_ElectrodeChannel_12.png

# low threshold nspikes 12,252,875
            # n onset-offsets: 294
# high threshold nspikes 1,986,462
            # n onset-offsets: 108

