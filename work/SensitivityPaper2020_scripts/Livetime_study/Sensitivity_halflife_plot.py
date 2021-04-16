import pandas
import os
from matplotlib import pyplot as plt
import numpy as np
import histlite as hl

# dates = ["20_11_12_Binary_023", "20_11_12_DNN1_023", "20_11_12_DNN1_024"]
# dates = ["20_10_22", "20_10_21", "20_10_19"]
dates = ["21_03_11_DNN1_024"]
root_dir = '/p/lustre2/nexouser/czyz1/output/'
# root_dir = '/Users/czyz1/lc-nexouser/output'
# root_dir = '/Users/czyz1/OneDrive - LLNL/Projects/nEXO/testoutput/'

num_datasets = 100
start_it = 0
end_its = 100
lts = ['0.25', '0.5', '1.0', '2.0', '5.0', '10.0']

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
    eff = 0.9598  # hit efficiency
    sensitivity = eff * atoms136 * lt_years * np.log(2) / cross_median

    return sensitivity


if __name__ == "__main__":

    atoms136 = calc_atoms_136()
    fig2, ax2 = plt.subplots()

    colors = ['b', 'r', 'g', 'k', 'm', 'c']

    for date in dates:

        num_its = end_its - start_it
        num_toys = num_its * len(lts) * num_datasets

        fig, ax = plt.subplots()
        dnn = date.split('_')[-2]
        materialdb = date.split('_')[-1]
        failed = 0
        converged = 0
        sensitivity = []

        # input_dir = root_dir + "h5/" + date
        input_dir = root_dir + "h5/reanalyzed/" + date
        output_dir = root_dir + "plots/"# + date

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        for plot_index, lt_years in enumerate(lts):
            num_runs = 0
            crossing_masked = []


            for iter in range(start_it, end_its):
                num_runs += 1
                filename = '{}/sens_output_file_90CL_{}_livetime_{}.h5'.format(input_dir, iter, lt_years)
                if not os.path.exists(filename):
                    num_toys -= num_datasets
                    failed += 1
                    print('livetime = {}, iteration = {}, failed = {}'.format(float(lt_years), iter, failed))
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

                # if 'best_fit_nll' in df.columns:
                #     print('placeholder')

                crossing_masked = crossing_masked + [f for a, b, c, d, e, f in zip(df['best_fit_converged'],
                                     df['best_fit_covar'], fix_fit_converged, fix_acc_covar, df['lambda'],
                                     # df['90CL_crossing']) if (a and b and c and d and min(e) > -0.1 and f > 0)]
                                     df['90CL_crossing_EXACT']) if (a and b and c and d and min(e) > -0.1 and f > 0)]

            converged += len(crossing_masked)
            sensitivity.append(sensitivity_calc(atoms136, float(lt_years), np.median(crossing_masked)))
            ax.axvline(np.median(crossing_masked), color=colors[plot_index % len(colors)], linestyle='--')
            hteststats = hl.hist(crossing_masked, bins=np.linspace(0, 30, 121))
            hl.plot1d(ax, hteststats,  color=colors[plot_index % len(colors)], label='LT = {} Y, T_1/2 = {:.3e} Y,'
                                                                              '\n Median = {:.3}'
                      .format(float(lt_years), sensitivity[-1], np.median(crossing_masked)))
            ax.set_xlabel('Upper limit on 0nuBB counts at 90% confidence limit', fontsize=16)
            ax.set_ylabel('Counts ({} toys, {} converged)'.format(num_toys, converged), fontsize=16)
            ax.set_xlim([0, 30])
            ax.legend()

        fig.savefig('{}/sens_hist_{}.png'.format(output_dir, date, dnn, materialdb), dpi=800)

        ax2.plot([float(lt) for lt in lts], sensitivity, '-o', label='MatDB = {}, DNN = {}, t_1/2,10 = {:.3e}'.format(
            materialdb, dnn, sensitivity[-1]))
        ax2.set_xlim([0, 12])
        ax2.set_ylim([1E26, 1.7E28])
        ax2.set_xlabel('Livetime, years')
        ax2.set_ylabel('Sensitivity (T_1/2)')
        ax2.set_yscale('linear')
        ax2.legend()

    fig2.savefig('{}/sens_vs_hl_{}.png'.format(output_dir, date), dpi=800)

    # fig.show()
    # fig2.show()
