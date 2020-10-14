import os
import MUA_constants as const

if not const.P["outputPath"].endswith('output_lowthr/MUA_output'):
    print('Outputpath wrong for this script? Check MUA_init.')
    exit()

"""Process data """
from MUA_utility import compress_CSVs, process_data
# process_data()
# compress_CSVs()



"""Mapping cortical and thalamic channels"""
# from map_chanels import cortical_mapping_panel
# from map_channels import plot_first_response
# from map_channels import thalamic_mapping_panel
# cortical_mapping_panel(dest_dir_appdx='../cortical_mapping')
# first_response_dir = '../first_response'
# plot_first_response(dest_dir_appdx=first_response_dir)
# anatomy_dir = '/media/loaloa/gdrive/projects/ephys/metadata/Histology_pngs'
# thalamic_mapping_panel(dest_dir_appdx = '../thalamus_mapping', 
#                        anatomy_dir = anatomy_dir,
#                        ts_plots_dir = f'{const.P["outputPath"]}/{first_response_dir}')



"""Noise subtraction over experiment timeline"""
# from plotting import firingrate_noise_timeline
# firingrate_noise_timeline(dest_dir_appdx='../plots/firingrates_heatmaps')
# firingrate_noise_timeline(dest_dir_appdx='../plots/firingrates_heatmaps', fname_postfix='_parad_subtr', subtr_noise='paradigm_wise')
# firingrate_noise_timeline(dest_dir_appdx='../plots/firingrates_heatmaps', fname_postfix='_deviant_subtr', subtr_noise='deviant_alone')



"""General explorative / summarizing  plots"""
# from plotting import firingrate_heatmaps
# firingrate_heatmaps(dest_dir_appdx='../plots/firingrates_heatmaps', subtr_noise='paradigm_wise')
# firingrate_heatmaps(dest_dir_appdx='../plots/frates_whiskerwise_regionwise', chnls_to_regions=True, subtr_noise='paradigm_wise', grouping='whisker_wise')
# firingrate_heatmaps(dest_dir_appdx='../plots/frates_whiskerwise', subtr_noise='paradigm_wise', grouping='whisker_wise')
# firingrate_heatmaps(dest_dir_appdx='../plots/frates_whiskerwise_reduced', subtr_noise='paradigm_wise', grouping='whisker_wise_reduced')


"""SSA and SSA correlation"""
# start, stop = 5, 20
# plotting.oddball10_si(which='O10', dest_dir_appdx='../plots/oddball_si', fname_appdx=f'{start}-{stop}ms_Oddball_10')
# plotting.oddball10_si(which='O25', dest_dir_appdx='../plots/oddball_si', fname_appdx=f'{start}-{stop}ms_Oddball_25')
# plotting.oddball10_si(which='MS', dest_dir_appdx='../plots/oddball_si', fname_appdx=f'{start}-{stop}ms_Oddball_25vsMS')
# plotting.ssa_correlation(which='O10', dest_dir_appdx='../plots/si_corr_5-20ms_010_late', fname_appdx='', post_stim=True)
# plotting.ssa_correlation(which='O10', dest_dir_appdx='../plots/si_corr_5-20ms_010_check', fname_appdx='', post_stim=False)


"""Onset offset"""
from onset_offset_classif import onset_offset_response
from onset_offset_classif import onset_offset_labels
from onset_offset_classif import lapl_kernel_SVM
from onset_offset_classif import get_onset_offset_classification
from onset_offset_classif import onoff_heatmap
from onset_offset_classif import onoff_barplot
from onset_offset_classif import idk_yet

# onset_offset_response(plots_dest_dir_appdx='../plots/onset_offset_chnls', csv_dest_dir_appdx='../onset_offset_model')
# onset_offset_response(plots_dest_dir_appdx='../plots/onset_offset_regions', single_channels=False, csv_dest_dir_appdx=None)
# onset_offset_labels(dest_dir_appdx='../onset_offset_model')

training_data_dir = '/mnt/Samsung_T5/output_lowthr/onset_offset_model'
# lapl_kernel_SVM(training_data_dir=training_data_dir, dest_dir_appdx='../onset_offset_model', parameter_search=True)
# lapl_kernel_SVM(training_data_dir=training_data_dir, dest_dir_appdx='../onset_offset_model', plot_cv_result=True)
# lapl_kernel_SVM(training_data_dir=training_data_dir, analyize_confusions=True)u

training_data_chnl_map_file = '/mnt/Samsung_T5/output_lowthr/chnls_map.csv'

# training data
# train_data = get_onset_offset_classification(which_data='training', training_data_dir=training_data_dir, 
#                                              training_data_chnl_map_file=training_data_chnl_map_file)
# onoff_heatmap(train_data, dest_dir_appdx='../onset_offset_model/plots', fig_height=11)
# train_data_all = get_onset_offset_classification(which_data='training', training_data_dir=training_data_dir, keep_labels=[0,1,2,3],
#                                                  training_data_chnl_map_file=training_data_chnl_map_file)
# onoff_barplot(train_data_all, dest_dir_appdx='../onset_offset_model/plots')


# idk_yet('mGE84', chnl_map_file=training_data_chnl_map_file, dest_dir_appdx='../plots/rastercanvas')
# idk_yet('mGE83', chnl_map_file=training_data_chnl_map_file, dest_dir_appdx='../plots/rastercanvas')
