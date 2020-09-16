# coding: utf-8
import MUA_constants as const
from preprocessing import compress_CSVs, process_data
import plotting


"""Process data """
# process_data()
# compress_CSVs()


"""General explorative / summarizing  plots"""
# plotting.firingrate_heatmaps(dest_dir_appdx='../plots/frates_whiskerwise', subtr_noise='paradigm_wise', grouping='whisker_wise')
# plotting.firingrate_heatmaps(dest_dir_appdx='../plots/frates_whiskerwise_reduced', subtr_noise='paradigm_wise', grouping='whisker_wise_reduced')
# plotting.firingrate_heatmaps('noisy', False)
# plotting.firingrate_noise_timeline('noisyy')
# plotting.firingrate_noise_timeline(fname_prefix='parad_subtr', subtr_noise='paradigm_wise')



"""Mapping cortical and thalamic channels"""
# plotting.make_CSD_summary_plots(lfp_output_appdx='/../LFP_output', dest_dir_appdx='/../CSD_lowthr')
# plotting.plot_time_to_first(dest_dir=f'{const.PROJ_DIR}/output/time_to_first')
# plotting.plot_evoked_lfp(dest_dir_appdx='../thalamus_mapping_lowthr', 
#                          anatomy_dir=f'{const.PROJ_DIR}/metadata/Histology_pngs',
                        #  ts_dir=f'{const.PROJ_DIR}/output/time_to_first')



"""SSA and SSA correlation"""
# start, stop = 5, 20
# plotting.oddball10_si(which='O10', dest_dir_appdx='../plots/oddball_si', fname_appdx=f'{start}-{stop}ms_Oddball_10')
# plotting.oddball10_si(which='O25', dest_dir_appdx='../plots/oddball_si', fname_appdx=f'{start}-{stop}ms_Oddball_25')
# plotting.oddball10_si(which='MS', dest_dir_appdx='../plots/oddball_si', fname_appdx=f'{start}-{stop}ms_Oddball_25vsMS')
# plotting.ssa_correlation(which='O10', dest_dir_appdx='../plots/si_corr_5-20ms_010_late', fname_appdx='', post_stim=True)
# plotting.ssa_correlation(which='O10', dest_dir_appdx='../plots/si_corr_5-20ms_010_check', fname_appdx='', post_stim=False)


"""Onset offset"""
from onset_offset_classif import onset_offset_response, onset_offset_labels, lapl_kernel_SVM, get_onset_offset_classification
# onset_offset_response(plots_dest_dir_appdx='../plots/onset_offset_chnls', csv_dest_dir_appdx='../onset_offset_model')
# onset_offset_response(plots_dest_dir_appdx='../plots/onset_offset_regions', single_channels=False, csv_dest_dir_appdx=None)
# onset_offset_labels(dest_dir_appdx='../onset_offset_model')

training_data_dir = '/mnt/Samsung_T5/output_lowthr/onset_offset_model'
# lapl_kernel_SVM(training_data_dir=training_data_dir, dest_dir_appdx='../onset_offset_model', parameter_search=True)
# lapl_kernel_SVM(training_data_dir=training_data_dir, dest_dir_appdx='../onset_offset_model', plot_cv_result=True)
# lapl_kernel_SVM(training_data_dir=training_data_dir, analyize_confusions=True)

training_data_chnl_map_file = f'{const.P["outputPath"]}/../../output_lowthr/chnls_map.csv'
MUA_output_data_chnl_map_file = f'{const.P["outputPath"]}/../S1Th_LayerAssignment_22.10.19.csv'

# get_onset_offset_classification(training_data_dir=training_data_dir, dest_dir_appdx='../plots/classifier_train_perf', rank='mouse',
#                       training_data_chnl_map_file=training_data_chnl_map_file, MUA_output_data_chnl_map_file=MUA_output_data_chnl_map_file)
get_onset_offset_classification(which_data='both', training_data_dir=training_data_dir, dest_dir_appdx='../pred_gwendata_lowthr',
                      training_data_chnl_map_file=training_data_chnl_map_file, MUA_output_data_chnl_map_file=MUA_output_data_chnl_map_file)
# get_onset_offset_classification(dest_dir_appdx='../plots/classifier_train_perf', rank='paradigm')
# get_onset_offset_classification(dest_dir_appdx='../plots/classifier_train_perf', rank='stimulus_type')
# get_onset_offset_classification(dest_dir_appdx='../plots/classifier_train_perf', rank='channel')

# get_onset_offset_classification(dest_dir_appdx='../plots/classifier_train_perf/splitmice/', 
#                                rank='', plot_labeled_data=True, print_labeled_data=False, 
#                                split_mice=True, )
# get_onset_offset_classification(dest_dir_appdx='../plots/classifier_train_perf/splitmice/', rank='paradigm',
#                                plot_labeled_data=True, print_labeled_data=False, split_mice=True)
