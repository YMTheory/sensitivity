import pandas
import os
from matplotlib import pyplot as plt
import numpy as np
import histlite as hl

dnn = 'DNN1'
materialdb = "024"
# # dates = ['Brian']
dates = ["21_03_01_{}_{}".format(dnn, "023"),
         "21_03_01_{}_{}".format(dnn, materialdb),
         "21_03_01_{}_{}".format(dnn, "025")]
# materialdb = "024"
# dates = ["21_01_21_{}_{}".format(dnn, materialdb)]#, "21_01_21_{}_023".format(dnn)]
# dates = ["20_10_22", "20_10_21", "20_10_19"]
# root_dir = '/p/lustre2/nexouser/czyz1/output/'
# path_home = '/p/lustre2/czyz1/nexo_sensitivity/work'
# path_result = '/p/lustre2/nexouser/czyz1'
path_home = '/Users/czyz1/lc-home/nexo_sensitivity/work'
path_result = '/Users/czyz1/lc-nexouser'
root_dir = '{}/output/'.format(path_result)


num_dataset = 100
start_it = 0
end_it = 100
lres = '0.008'
# eres = ['1']
# eres = ['0.011', '0.015', '0.018']
resols = ['0.008', '0.009', '0.01', '0.011', '0.012', '0.013', '0.014', '0.015', '0.016', '0.017', '0.018']
# resols = ['0.008', '0.01', '0.012', '0.014', '0.016', '0.018']

def calc_atoms_136():
    """ Number of Xe136 atoms in nEXO fiducial volume """
    mmass134 = 0.133905395  # kg/mol 134
    mmass136 = 0.135907219  # kg/mol 136
    at_frac = 0.9           # atomic fraction 136 / (136 + 134)
    avog_num = 6.022141E23  # Avogadro's number
    fid_mass = 3281         # mass of fiducial volume [kg]

    atoms136 = (fid_mass * avog_num * at_frac) / ((mmass136 * at_frac) + ((1 - at_frac) * mmass134))

    return atoms136


def sensitivity_calc(atoms136, lt_years, cross_median):
    """Calculate the sensitivty of nEXO in terms of half-life (years)"""
    eff = 0.963  # hit efficiency
    sensitivity = eff * atoms136 * lt_years * np.log(2) / cross_median

    return sensitivity


if __name__ == "__main__":

    lt_years = 10
    atoms136 = calc_atoms_136()
    fig2, ax2 = plt.subplots()

    colors = ['b', 'r', 'g', 'k', 'm']
    coppers = ['Electroform Copper', 'Aurubis Copper', 'D025']

    for index, date in enumerate(dates):
        if date == '20_11_04_DNN1_024':     # special case where the code was ran differently
            eres = ['0.008']
            end_its = 450
            num_datasets = 50
        elif date == '20_12_03_DNN1_023' or date == '21_03_01_DNN1_023' or date == '21_03_01_DNN1_025':
            eres = ['0.008']
            end_its = end_it
            num_datasets = num_dataset
        else:
            eres = resols
            end_its = end_it
            num_datasets = num_dataset

        num_its = end_its - start_it
        num_toys = num_its * len(eres) * num_datasets #* 2  # 2 = different lam curves created at different resolutions
                                                           # (0.008 and the resolution)
        print('num_toys =' + str(num_toys))

        fig, ax = plt.subplots()

        failed = 0
        sensitivity = []

        # input_dir = '/p/lustre2/lenardo1/sensitivity_output/Jan19_Rn222Study_merged-v10b_OptimizedV1Binning_D024_REANALYZED_WITH_STEVENS_CURVE/'
        input_dir = root_dir + "h5/reanalyzed/" + date

        output_dir = root_dir + "plots/"# + date

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        plot_index = 0
        converged = 0
        for res in eres:
            for lres in ['0.008']:
                print('lres = {}'.format(lres))
                num_runs = 0
                crossing_masked = []

                # if (res == '0.008' and date == '21_01_21_DNN1_024'):
                #     end_its = 100
                #     num_toys += num_datasets * end_it
                # else:
                #     end_its = end_it

                for iter in range(start_it, end_its):
                    num_runs += 1
                    # if iter < 10:
                    #     filename = '{}/reanalyzed_sens_output_file_rn222study_01.0x_90CL_00{}_D024.h5'.format(input_dir,
                    #                                                                                      iter)
                    # else:
                    #     filename = '{}/reanalyzed_sens_output_file_rn222study_01.0x_90CL_0{}_D024.h5'.format(input_dir, iter)

                    filename = '{}/reanalyzed_lres=0.008_sens_output_file_90CL_{}_resolution_{}.h5'.format(input_dir, iter, res)
                    # filename = '{}/sens_output_file_90CL_{}_resolution_{}.h5'.format(input_dir, iter, res)
                    if not os.path.exists(filename):
                        num_toys -= num_datasets
                        failed += 1
                        print('resolution = {}, iteration = {}, failed = {}'.format(res, iter, failed))
                        print(filename)
                        continue
                    df = pandas.read_hdf(filename)

                    fix_fit_converged = []
                    for i in df['fixed_fit_converged']:
                        numfalse = 0
                        for j in i:
                            if not j:
                                numfalse += 1
                        if numfalse > 2:
                            fix_fit_converged.append(False)
                        else:
                            fix_fit_converged.append(True)

                    fix_acc_covar = []
                    for i in df['fixed_fit_acc_covar']:
                        numfalse = 0
                        for j in i:
                            if not j:
                                numfalse += 1
                        if numfalse > 2:
                            fix_acc_covar.append(False)
                        else:
                            fix_acc_covar.append(True)

                    if 'best_fit_nll' in df.columns:
                        print('placeholder')


                    crossing_masked = crossing_masked + [f for a, b, c, d, e, f in zip(df['best_fit_converged'],
                                                 df['best_fit_covar'], fix_fit_converged, fix_acc_covar, df['lambda'],
                                                 df['90CL_crossing_EXACT']) if (a and b and c and d and min(e) > -0.1 and f > 0)]
                                                 # df['90CL_crossing']) if (a and b and c and d and min(e) > -0.1 and f > 0)]

                    # print('Length of Crossing Masked = ' + str(len(crossing_masked)))
                converged += len(crossing_masked)
                sensitivity.append(sensitivity_calc(atoms136, lt_years, np.median(crossing_masked)))
                ax.axvline(np.median(crossing_masked), color=colors[plot_index % len(colors)], linestyle='--')
                hteststats = hl.hist(crossing_masked, bins=np.linspace(0, 30, 31))
                print(hteststats.values)
                hl.plot1d(ax, hteststats,  color=colors[plot_index % len(colors)], label='E_res = {} Y, Lam_curve = {}'
                                                                                         '\n T_1/2 = {:.3e} Y,Median = {:.3}'
                          .format(res, lres, sensitivity[-1], np.median(crossing_masked)))
                ax.set_xlabel('Upper limit on 0nuBB counts at 90% confidence limit', fontsize=16)
                ax.set_ylabel('Counts ({} toys, {} converged)'.format(num_toys, converged), fontsize=16)
                ax.set_xlim([0, 30])
                ax.legend()

                plot_index += 1

            fig.savefig('{}/sens_res={}_lres=_{}_hist_{}.png'.format(output_dir,res, lres, date), dpi=800)

        eres_num = [float(i) * 100 for i in eres]
        if date == '20_11_30_DNN1_023' or date == '20_12_03_DNN1_023':
            eres_num = [0.8]
        ax2.plot(eres_num, sensitivity, '-o', label=coppers[index])
        ax2.set_xlim([.7, 1.9])
        ax2.set_ylim([5E27, 1.7E28])
        ax2.set_xlabel('Resolution $\sigma$Q$_{\\beta\\beta}$ [%]', fontsize=16)
        ax2.set_ylabel('$^{136}$Xe 0$\\nu\\beta\\beta$ T$_{1/2}$ [yr]', fontsize=16)
        ax2.set_yscale('linear')
        ax2.legend()

    fig2.savefig('{}/sens_vs_res_{}_lres={}.png'.format(output_dir, date, lres), dpi=800, bbox_inches='tight')

    fig.show()
    fig2.show()
