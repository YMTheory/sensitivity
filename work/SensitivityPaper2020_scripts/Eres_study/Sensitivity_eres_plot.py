import pandas
import os
from matplotlib import pyplot as plt
import numpy as np
import histlite as hl

dnn = 'DNN1'
materialdb = "024"
dates = ["21_01_21_{}_{}".format(dnn, materialdb), "21_01_21_{}_023".format(dnn)]
# dates = ["20_10_22", "20_10_21", "20_10_19"]
# root_dir = '/p/lustre2/nexouser/czyz1/output/'
path_home = '/p/lustre2/czyz1/nexo_sensitivity/work'
path_result = '/p/lustre2/nexouser/czyz1'
# path_home = '/Users/czyz1/lc-home/nexo_sensitivity/work'
# path_result = '/Users/czyz1/lc-nexouser'
root_dir = '{}/output/'.format(path_result)

num_dataset = 100
start_it = 0
end_it = 50
eres = ['0.008', '0.009', '0.01', '0.011', '0.012', '0.013', '0.014', '0.015', '0.016', '0.017', '0.018']

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
    coppers = ['Electroform Copper', 'Aurubis Copper']

    for index, date in enumerate(dates):
        if date == '20_11_04_DNN1_024':     # special case where the code was ran differently
            end_its = 450
            num_datasets = 50
        elif date == '20_12_03_DNN1_023':
            eres = ['0.008']
        else:
            end_its = end_it
            num_datasets = num_dataset

        num_its = end_its - start_it
        num_toys = num_its * len(eres) * num_datasets

        fig, ax = plt.subplots()

        failed = 0
        converged = 0
        sensitivity = []

        input_dir = root_dir + "h5/" + date
        output_dir = root_dir + "plots/"# + date

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        for plot_index, res in enumerate(eres):
            num_runs = 0
            crossing_masked = []


            for iter in range(start_it, end_its):
                num_runs += 1
                filename = '{}/sens_output_file_90CL_{}_resolution_{}.h5'.format(input_dir, iter, res)
                if not os.path.exists(filename):
                    num_toys -= num_datasets
                    failed += 1
                    print('resolution = {}, iteration = {}, failed = {}'.format(res, iter, failed))
                    print(filename)
                    continue
                df = pandas.read_hdf(filename)
                crossing_masked = crossing_masked + [b for a, b in zip(df['best_fit_converged'],
                                                                       df['90CL_crossing']) if (a and b > 0)]
            converged += len(crossing_masked)
            sensitivity.append(sensitivity_calc(atoms136, lt_years, np.median(crossing_masked)))
            ax.axvline(np.median(crossing_masked), color=colors[plot_index % len(colors)], linestyle='--')
            hteststats = hl.hist(crossing_masked, bins=np.linspace(0, 30, 31))
            hl.plot1d(ax, hteststats,  color=colors[plot_index % len(colors)], label='E_res = {} Y, T_1/2 = {:.3e} Y,'
                                                                                     '\n Median = {:.3}'
                      .format(res, sensitivity[-1], np.median(crossing_masked)))
            ax.set_xlabel('Upper limit on 0nuBB counts at 90% confidence limit', fontsize=16)
            ax.set_ylabel('Counts ({} toys, {} converged)'.format(num_toys, converged), fontsize=16)
            ax.set_xlim([0, 30])
            ax.legend()

        fig.savefig('{}/sens_hist_{}_{}.png'.format(output_dir, dnn, materialdb), dpi=800)

        eres_num = [0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8]
        if date == '20_11_30_DNN1_023' or date == '20_12_03_DNN1_023':
            eres_num = [0.8]
        ax2.plot(eres_num, sensitivity, '-o', label=coppers[index])
        ax2.set_xlim([.7, 1.9])
        ax2.set_ylim([5E27, 1.3E28])
        ax2.set_xlabel('Resolution $\sigma$Q$_{\\beta\\beta}$ [%]', fontsize=16)
        ax2.set_ylabel('$^{136}$Xe 0$\\nu\\beta\\beta$ T$_{1/2}$ [yr]', fontsize=16)
        ax2.set_yscale('linear')
        ax2.legend()

    fig2.savefig('{}/sens_vs_res_{}.png'.format(output_dir, date), dpi=800, bbox_inches='tight')

    fig.show()
    fig2.show()
