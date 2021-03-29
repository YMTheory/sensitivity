# import sys
# import os
# import pandas as pd
# import histlite as hl
# import numpy as np
# from matplotlib import pyplot as plt
# from scipy import stats
#
# sys.path.append('../../modules')
#
# # Import the nEXO sensitivity classes
# import nEXOFitWorkspace
# import nEXOFitModel
# import nEXOFitLikelihood
#
# path_home = '/Users/czyz1/lc-home/nexo_sensitivity/work'
# path_result = '/Users/czyz1/lc-nexouser'
#
# config_loc = "{}/config/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1.yaml".format(path_home)
# date = '21_01_13'
# database_num = '023'
# # Set some switches
# lt_years = 10
#
#
#
# for resolution in ['0.008']:#, '0.009', '0.01', '0.011', '0.012', '0.013', '0.014', '0.015', '0.016', '0.017', '0.018']:
#
#     hist_loc = "{}/workdir/histogram_files/{}/Baseline2019_Histograms_Energy_Res={}.h5".\
#         format(path_result, date, resolution)
#     df_mc_histograms = pd.read_hdf(hist_loc, key='SimulationHistograms')
#     print('hi')
#     # # Create the workspace
#     # workspace = nEXOFitWorkspace.nEXOFitWorkspace(config_loc)
#     # workspace.SetHandlingOfRadioassayData(fluctuate=True)
#     #
#     # workspace.LoadComponentsTableFromFile(comp_loc)
#     # workspace.livetime = lt_years * 365.25 * 24. * 60. * 60.
#     # workspace.CreateGroupedPDFs()
#     #
#     # # Create the likelihood object
#     # likelihood = nEXOFitLikelihood.nEXOFitLikelihood()
#     # likelihood.AddPDFDataframeToModel(workspace.df_group_pdfs, workspace.histogram_axis_names)

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# import sys
# sys.path.append('../../../modules')
#
# import os
# import pandas as pd
# import histlite as hl
# import numpy as np
# from matplotlib import pyplot as plt
# import nEXOFitWorkspace
# import nEXOFitModel
# import nEXOFitLikelihood
# import nEXOMaterialsDBInterface
# import importlib
# importlib.reload( nEXOMaterialsDBInterface )
# importlib.reload( nEXOFitWorkspace )
# importlib.reload( nEXOFitModel )
#
# import matplotlib
# from matplotlib.backends.backend_pdf import PdfPages
#
# plt.rcParams.update({'font.size': 18})
# plt.rcParams['figure.figsize'] = (10,8)
#
# # path_home = '/Users/czyz1/lc-home/nexo_sensitivity/work'
# # path_result = '/Users/czyz1/lc-nexouser'
# path_home = '/p/lustre2/czyz1/nexo_sensitivity/work'
# path_result = '/p/lustre2/nexouser/czyz1'
# config_loc = "{}/config/Sensitivity2020_Optimized_DNN_Standoff_Binning_version1.yaml".format(path_home)
# date = '21_03_01'
# database_num = '025'
# res = ['0.008']#, '0.009', '0.01', '0.011', '0.012', '0.013', '0.014', '0.015', '0.016', '0.017', '0.018']
#
# optimized_binning_yaml = config_loc
#
# for resolution in res:
#     workspace2020 = nEXOFitWorkspace.nEXOFitWorkspace(config=optimized_binning_yaml)
#     workspace2020_ind = nEXOFitWorkspace.nEXOFitWorkspace(config=optimized_binning_yaml)
#
#
#     optimized_binning_components_table = "{}/workdir/components_tables/{}/ComponentsTable_D-{}_Energy_Res={}.h5".\
#             format(path_result, date, database_num, resolution)
#     # optimized_binning_components_table = '/p/vast1/nexo/sensitivity2020/pdfs/component_tables/ComponentsTable_D-{}_merged-v11_Optimized_DNN_Standoff_Binning_version1.h5'.format(database_num)
#     workspace2020.LoadComponentsTableFromFile(optimized_binning_components_table)
#     workspace2020_ind.LoadComponentsTableFromFile(optimized_binning_components_table)
#
#
#     groupdict2020 = {}
#     for index, row in workspace2020_ind.df_components.iterrows():
#         if row['Histogram'] is None:
#             continue
#         groupdict2020[row['PDFName']] = row['Group']
#         workspace2020_ind.df_components.loc[index,'Group'] = row['PDFName']
#
#
#     groups = set(groupdict2020.values())
#     group_component_dict = dict()
#     for group in groups:
#         group_component_dict[group] = [k for k, v in groupdict2020.items() if v == group]
#
#
#     workspace2020.CreateGroupedPDFs()
#     model2020 = nEXOFitModel.nEXOFitModel()
#     model2020.AddPDFsFromDataframe(workspace2020.df_group_pdfs,\
#                                    workspace2020.histogram_axis_names)
#     model2020.GenerateModelDistribution()
#
#
#     workspace2020_ind.CreateGroupedPDFs()
#     model2020_ind = nEXOFitModel.nEXOFitModel()
#     model2020_ind.AddPDFsFromDataframe(workspace2020_ind.df_group_pdfs,\
#                                    workspace2020_ind.histogram_axis_names)
#     model2020_ind.GenerateModelDistribution()
#
#     # Relative contributions to ROI
#     roi_cut_dict_2020 = {'DNN': (0.86,1.),
#                         'Energy (keV)': (2434., 2480.),
#                         'Standoff (mm)': (104.5, 650.)
#                         }
#
#     workspace2020.DefineROI( roi_cut_dict_2020 )
#     workspace2020_ind.DefineROI( roi_cut_dict_2020 )
#
#     ##
#     total_in_ROI = model2020.GetIntegralInBinRange(workspace2020.GetROIBinIndices())
#     print('{:4.4} total events in ROI'.format(total_in_ROI))
#
#     print('{:<20} {:>13} {:>13}'.format('Name','Counts','Percent'))
#
#     group_bkg_dict = dict()
#
#     for row in model2020.variable_list:
#         this_group_in_roi = model2020.GetComponentIntegralInBinRange(row['Name'],workspace2020.GetROIBinIndices())
#         print('{:<20} {:>13.3} {:>13.3}%'.format(row['Name'],\
#                                                 this_group_in_roi,\
#                                                 this_group_in_roi/total_in_ROI*100.))
#         group_bkg_dict[row['Name'][4:]] = this_group_in_roi
#
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# ##
# plt.rcParams.update({'font.size': 18})
# plt.rcParams['figure.figsize'] = (10, 8)
#
# spectrum_cut_dict_2020 = {'DNN': (0.86, 1.),
#                           'Energy (keV)': (1000., 3500.),
#                           'Standoff (mm)': (100., 650.)
#                           }
#
# counter = 0
# for group in groups:
#
#     plt.figure(counter)
#
#     if group == 'Off': continue
#
#     component_list = group_component_dict[group]
#
#     # Get first component histogram, so we can get the shape and create an empty sum_hist
#     idx = model2020_ind.GetVariableIndexByName(component_list[0])
#     test_hist = model2020_ind.GetSlicedDistribution(cut_dict=spectrum_cut_dict_2020, \
#                                                     renormalize=False, \
#                                                     var_name=component_list[0], \
#                                                     verbose=False)
#     sum_hist = hl.Hist(bins=test_hist.bins, values=np.zeros(test_hist.values.shape))
#
#     print('\n********** {} ***********'.format(group))
#     print('{:<33} {:>13} {:>16}'.format('Name', 'Counts (ROI)', 'Percent of grp'))
#
#     for component in component_list:
#         this_hist = model2020_ind.GetSlicedDistribution(cut_dict=spectrum_cut_dict_2020, \
#                                                         renormalize=False, \
#                                                         var_name=component, \
#                                                         verbose=False)
#         idx = model2020_ind.GetVariableIndexByName(component)
#         this_val = model2020_ind.variable_list[idx]['Value']
#         hl.plot1d(this_hist.project([1]) * this_val, label=component)
#
#         sum_hist = sum_hist + this_hist * this_val
#
#         this_comp_in_roi = model2020_ind.GetComponentIntegralInBinRange(component, \
#                                                                         workspace2020_ind.GetROIBinIndices())
#         print('{:<33} {:>13.3} {:>16.4}%'.format(component, \
#                                                  this_comp_in_roi, \
#                                                  this_comp_in_roi / group_bkg_dict[group] * 100.))
#
#     hl.plot1d(sum_hist.project([1]), color='k', linewidth=2, label="Sum")
#     plt.yscale('log')
#     plt.legend(fontsize=8)
#     plt.ylabel('Single-site counts in 10yr in inner ~2T')
#     plt.xlabel('Energy (keV)')
#     plt.title(group)
#     plt.xlim(1000., 3500.)
#
#     counter += 1
#
# ##
# plt.rcParams.update({'font.size': 18})
# plt.rcParams['figure.figsize'] = (10, 8)
#
# spectrum_cut_dict_2020 = {'DNN': (0.86, 1.),
#                           'Energy (keV)': (1000., 3500.),
#                           'Standoff (mm)': (100., 650.)
#                           }
#
# counter = 0
# for group in groups:
#
#     if 'Rn222' not in group: continue
#
#     plt.figure(counter)
#
#     # if group == 'Off': continue
#
#     component_list = group_component_dict[group]
#
#     # Get first component histogram, so we can get the shape and create an empty sum_hist
#     idx = model2020_ind.GetVariableIndexByName(component_list[0])
#     test_hist = model2020_ind.GetSlicedDistribution(cut_dict=spectrum_cut_dict_2020, \
#                                                     renormalize=False, \
#                                                     var_name=component_list[0], \
#                                                     verbose=False)
#     sum_hist = hl.Hist(bins=test_hist.bins, values=np.zeros(test_hist.values.shape))
#
#     print('\n********** {} ***********'.format(group))
#     print('{:<33} {:>13} {:>16}'.format('Name', 'Counts (ROI)', 'Percent of grp'))
#
#     for component in component_list:
#         this_hist = model2020_ind.GetSlicedDistribution(cut_dict=spectrum_cut_dict_2020, \
#                                                         renormalize=False, \
#                                                         var_name=component, \
#                                                         verbose=False)
#         idx = model2020_ind.GetVariableIndexByName(component)
#         this_val = model2020_ind.variable_list[idx]['Value']
#         hl.plot1d(this_hist.project([1]) * this_val, label=component)
#
#         sum_hist = sum_hist + this_hist * this_val
#
#         this_comp_in_roi = model2020_ind.GetComponentIntegralInBinRange(component, \
#                                                                         workspace2020_ind.GetROIBinIndices())
#         print('{:<33} {:>13.3} {:>16.4}%'.format(component, \
#                                                  this_comp_in_roi, \
#                                                  this_comp_in_roi / group_bkg_dict[group] * 100.))
#
#     # Now look at all "off" components
#     component_list = group_component_dict['Off']
#
#     for component in component_list:
#         if 'Rn222' not in component: continue
#         this_hist = model2020_ind.GetSlicedDistribution(cut_dict=spectrum_cut_dict_2020, \
#                                                         renormalize=False, \
#                                                         var_name=component, \
#                                                         verbose=False)
#         idx = model2020_ind.GetVariableIndexByName(component)
#         this_val = model2020_ind.variable_list[idx]['Value']
#         hl.plot1d(this_hist.project([1]) * this_val, label=component)
#
#         # sum_hist = sum_hist + this_hist*this_val
#
#         this_comp_in_roi = model2020_ind.GetComponentIntegralInBinRange(component, \
#                                                                         workspace2020_ind.GetROIBinIndices())
#         print('{:<33} {:>13.3} {:>16.4}%'.format(component, \
#                                                  this_comp_in_roi, \
#                                                  this_comp_in_roi / group_bkg_dict[group] * 100.))
#
#     hl.plot1d(sum_hist.project([1]), color='k', linewidth=2, label="Sum")
#     plt.yscale('log')
#     plt.legend(fontsize=10)
#     plt.ylabel('Single-site counts in 10yr in inner ~2T')
#     plt.xlabel('Energy (keV)')
#     plt.title(group)
#     plt.xlim(1000., 3500.)
#
#     counter += 1
#
# ##
# int_idx = model2020.GetVariableIndexByName('Internals_U238')
# rn222_idx = model2020.GetVariableIndexByName('Rn222')
#
# ##
# plt.rcParams.update({'font.size': 18})
# plt.rcParams['figure.figsize'] = (8,6)
#
# histint = model2020.pdfs[int_idx]
#
#
# hl.plot2d(histint.project([0,1]).log10())
# plt.xlabel('DNN')
# plt.ylabel('Energy (keV)')
#
#
# plt.figure(2)
# hl.plot2d(histint.project([1,2]).log10())
# plt.xlabel('Energy (keV)')
# plt.ylabel('Standoff (mm)')
#
# plt.figure(3)
# hl.plot2d(histint.project([0,2]).log10())
# plt.xlabel('DNN')
# plt.ylabel('Standoff (mm)')
#
# ##
# plt.rcParams.update({'font.size': 18})
# plt.rcParams['figure.figsize'] = (8,6)
#
# histrn222 = model2020.pdfs[rn222_idx]
#
#
# hl.plot2d(histrn222.project([0,1]).log10())
# plt.xlabel('DNN')
# plt.ylabel('Energy (keV)')
#
#
# plt.figure(2)
# hl.plot2d(histrn222.project([1,2]).log10())
# plt.xlabel('Energy (keV)')
# plt.ylabel('Standoff (mm)')
#
# plt.figure(3)
# hl.plot2d(histrn222.project([0,2]).log10())
# plt.xlabel('DNN')
# plt.ylabel('Standoff (mm)')
#
#
# ##
# eproj_cut_dict = {'DNN': (0.86,1.),
#                     'Energy (keV)': (1000., 3500.),
#                     'Standoff (mm)': (104.5, 650.)
#                     }
# stproj_cut_dict = {'DNN': (0.86,1.),
#                     'Energy (keV)': (2434., 2480.),
#                     'Standoff (mm)': (0., 650.)
#                     }
# dnn_cut_dict = {'DNN': (0.,1.),
#                     'Energy (keV)': (2434., 2480.),
#                     'Standoff (mm)': (104.5, 650.)
#                     }
#
# rn222_in_roi = model2020.GetComponentIntegralInBinRange('Rn222',\
#                                                         workspace2020.GetROIBinIndices())
# rn222_idx = model2020.GetVariableIndexByName('Rn222')
# int_in_roi = model2020.GetComponentIntegralInBinRange('Internals_U238',\
#                                                         workspace2020.GetROIBinIndices())
# int_idx = model2020.GetVariableIndexByName('Internals_U238')
# bb0n_in_roi = model2020.GetComponentIntegralInBinRange('Bb0n',\
#                                                         workspace2020.GetROIBinIndices())
# bb0n_idx = model2020.GetVariableIndexByName('Bb0n')
#
# rn222_scale = model2020.variable_list[rn222_idx]['Value']/rn222_in_roi
# int_scale = model2020.variable_list[int_idx]['Value']/int_in_roi
# bb0n_scale = model2020.variable_list[bb0n_idx]['Value']/bb0n_in_roi
#
#
# # Get the energy projections
# rn222eproj = model2020.GetSlicedDistribution(cut_dict=eproj_cut_dict,\
#                                                 renormalize=False,\
#                                                 var_name='Rn222',\
#                                                 verbose=False)
# inteproj = model2020.GetSlicedDistribution(cut_dict=eproj_cut_dict,\
#                                                 renormalize=False,\
#                                                 var_name='Internals_U238',\
#                                                 verbose=False)
# bb0neproj = model2020.GetSlicedDistribution(cut_dict=eproj_cut_dict,\
#                                                 renormalize=False,\
#                                                 var_name='Bb0n',\
#                                                 verbose=False)
# #model2020.GetComponentIntegralInBinRange()
# hl.plot1d( rn222eproj.project([1])*rn222_scale, color='b', label='Rn222')
# hl.plot1d( inteproj.project([1])*int_scale, color='r', label='Internals_U238')
# hl.plot1d( bb0neproj.project([1])*bb0n_scale,color='g', label='Bb0n')
#
# plt.xlabel('Energy (keV)')
# plt.ylabel('Normalized to counts in ROI')
# plt.title('Energy spectrum for DNN>0.86 in inner 2T',fontsize=16)
# plt.legend(fontsize=16)
# plt.xlim(2000.,3000.)
# #plt.yscale('log')
# plt.ylim(0.,1.)
#
# # plt.savefig('plots/rn222_internalsU238_energy_spectrum_comparison.png',\
# #             dpi=200,bbox_inches='tight')
#
# plt.figure(2)
#
# rn222stproj = model2020.GetSlicedDistribution(cut_dict=stproj_cut_dict,\
#                                                 renormalize=False,\
#                                                 var_name='Rn222',\
#                                                 verbose=False)
# intstproj = model2020.GetSlicedDistribution(cut_dict=stproj_cut_dict,\
#                                                 renormalize=False,\
#                                                 var_name='Internals_U238',\
#                                                 verbose=False)
# bb0nstproj = model2020.GetSlicedDistribution(cut_dict=stproj_cut_dict,\
#                                                 renormalize=False,\
#                                                 var_name='Bb0n',\
#                                                 verbose=False)
# hl.plot1d( rn222stproj.project([2])*rn222_scale, color='b', label='Rn222')
# hl.plot1d( intstproj.project([2])*int_scale, color='r', label='Internals_U238')
# hl.plot1d( bb0nstproj.project([2])*bb0n_scale,color='g', label='Bb0n')
#
# plt.xlabel('Standoff (mm)')
# plt.ylabel('Normalized to counts in ROI')
# plt.title('Standoff for DNN>0.86 and energy between 2430-2480 keV',fontsize=16)
# plt.legend(fontsize=16)
# #plt.yscale('log')
# plt.ylim(0.,1.)
# # plt.savefig('plots/rn222_internalsU238_standoff_comparison.png',\
# #             dpi=200,bbox_inches='tight')
#
#
# plt.figure(3)
#
# rn222dnnproj = model2020.GetSlicedDistribution(cut_dict=dnn_cut_dict,\
#                                                 renormalize=False,\
#                                                 var_name='Rn222',\
#                                                 verbose=False)
# intdnnproj = model2020.GetSlicedDistribution(cut_dict=dnn_cut_dict,\
#                                                 renormalize=False,\
#                                                 var_name='Internals_U238',\
#                                                 verbose=False)
# bb0ndnnproj = model2020.GetSlicedDistribution(cut_dict=dnn_cut_dict,\
#                                                 renormalize=False,\
#                                                 var_name='Bb0n',\
#                                                 verbose=False)
# hl.plot1d( rn222dnnproj.project([0])*rn222_scale, color='b', label='Rn222')
# hl.plot1d( intdnnproj.project([0])*int_scale, color='r', label='Internals_U238')
# hl.plot1d( bb0ndnnproj.project([0])*bb0n_scale,color='g', label='Bb0n')
#
# plt.xlabel('DNN')
# plt.ylabel('Normalized to counts in ROI')
# plt.title('DNN for standoff > 100mm and energy between 2430-2480 keV',fontsize=16)
# plt.legend(fontsize=16)
# #plt.yscale('log')
# plt.ylim(0.,1.)
# # plt.savefig('plots/rn222_internalsU238_DNN_comparison.png',\
# #             dpi=200,bbox_inches='tight')


# Scratch work
# import os
#
# loc = '/p/lustre2/nexouser/czyz1/workdir/lambda/21_03_01_DNN1_023'
# files = os.listdir(loc)
# a = []
# for file in files:
#     # print(file[-3:])
#     if file[-3:] == '.h5':
#         a.append(int(file.split('_')[3]))
#     else:
#         print(file[-3:] == '.h5')
# a.sort()
# print(a)
# missing = []
# for num in range(121):
#     if num != a[num-len(missing)]:
#         missing.append(num)
# print(missing)
#
# lam_list{
# '023' : [23, 31, 32, 33, 54, 58, 71, 72, 95, 98, 102, 103]
# '024' : [12, 13, 16, 18, 20, 23, 24, 28, 30, 35, 36, 39, 40, 43, 45, 46, 47, 48, 49, 50, 52, 53, 54, 55, 56, 58, 59,  60, 62, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 84, 85, 86, 87, 88, 89, 90, 91, 92, 94, 96, 97, 101, 102, 103, 107, 109, 112, 113, 114, 115, 116, 118, 120]
# '025' : [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 16, 17, 18, 19, 22, 23, 24, 25, 27, 28, 30, 32, 34, 37, 38, 39, 40,  41, 42, 43, 44, 45, 46, 47, 48, 49, 52, 55, 56, 57, 58, 60, 62, 63, 64, 65, 66, 67, 70, 71, 72, 73, 75, 82, 83, 84, 85, 86, 87, 88, 89, 92, 93, 95, 99, 101, 102, 103, 105, 108, 109, 116, 117, 119, 121]
# }