import os
from MUA_constants import P

# change MUA_constans ALL_MICE as well!
if not P["outputPath"].endswith('output_all_data_highthr/MUA_output'):
    print('Outputpath wrong for this script? Check MUA_init.')
    exit()


"""Process data """
# from MUA_utility import compress_CSVs, process_data
# process_data(how={'nbatches': 3, 'batch':0})
# compress_CSVs()


"""Onset offset"""
# from onset_offset_classif import onset_offset_response
# from onset_offset_classif import onset_offset_labels
# from onset_offset_classif import lapl_kernel_SVM
# from onset_offset_classif import get_onset_offset_classification
# from onset_offset_classif import onoff_heatmap
# from onset_offset_classif import onoff_barplot

# # this is a merge of the training data mapping and gwen data
# MUA_output_data_chnl_map_file = '/mnt/Samsung_T5/output_all_data_highthr/S1th_LayerAssignment_22.10.19.csv'
# training_data_dir = '/mnt/Samsung_T5/output_lowthr/onset_offset_model'

# # both training (but highthreshold, trained on low threshold) and new data
# dest_dir_appdx = '../pred_alldata_highthr'
# dest_dir_plot_appdx = '../pred_alldata_highthr/alldata_plots'
# os.makedirs(f'{P["outputPath"]}/{dest_dir_plot_appdx}', exist_ok=True)

# # Careful here, the classification is from the raw output of the SVM,
# # I accidentally deleted my filtered classification ...
# both_data = get_onset_offset_classification(which_data='MUA_output', training_data_dir=training_data_dir, dest_dir_appdx=dest_dir_appdx, cached_prediction=True, 
#                                             MUA_output_data_chnl_map_file=MUA_output_data_chnl_map_file)
# onoff_heatmap(both_data, dest_dir_appdx=dest_dir_plot_appdx, fig_height=15)

# both_data = get_onset_offset_classification(which_data='MUA_output', training_data_dir=training_data_dir, dest_dir_appdx=dest_dir_appdx, keep_labels=[0,1,2,3],
#                                             MUA_output_data_chnl_map_file=MUA_output_data_chnl_map_file)
# onoff_barplot(both_data, dest_dir_appdx=dest_dir_plot_appdx)