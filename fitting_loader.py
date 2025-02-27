import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import sys
sys.path.append('/fs/ddn/sdf/group/nexo/users/miaoyu/Sterile_nu/sensitivity/')
from contour_tool import *

def load_other_experiment_sensitivities():
    path = '/fs/ddn/sdf/group/nexo/users/miaoyu/Sterile_nu/sensitivity/data/'
    arr = np.loadtxt(path+'reactor_global_allowed_region_0.csv')
    x_reactor0, y_reactor0 = arr[:,0], arr[:, 1]
    arr = np.loadtxt(path+'reactor_global_allowed_region_1.csv')
    x_reactor1, y_reactor1= arr[:,0], arr[:, 1]
    arr = np.loadtxt(path+'reactor_global_allowed_region_2.csv')
    x_reactor2, y_reactor2 = arr[:,0], arr[:, 1]
    arr = np.loadtxt(path+'reactor_global_allowed_region_3.csv')
    x_reactor3, y_reactor3 = arr[:,0], arr[:, 1]

    arr = np.loadtxt(path+'LZ_natural_LXe_shapeonly.csv')
    x_LZ_natural, y_LZ_natural = arr[:, 0], arr[:, 1]

    #arr = np.loadtxt(path+'gallium_allowed_90CL.csv')
    #x_gallium, y_gallium = arr[:, 0], arr[:, 1]

    arr = np.loadtxt(path+'gallium_anomaly.csv')
    x_gallium, y_gallium = arr[:, 0], arr[:, 1]

    arr = np.loadtxt(path+'carbon12.csv')
    x_c12, y_c12 = arr[:, 0], arr[:, 1]

    data_dict = {'Reactor': {0: [x_reactor0, y_reactor0], \
                             1: [x_reactor1, y_reactor1], \
                             2: [x_reactor2, y_reactor2], \
                             3: [x_reactor3, y_reactor3], \
                             }, \
                  'LZ_natural' : [x_LZ_natural, y_LZ_natural],\
                  'Gallium': [x_gallium, y_gallium], \
                  'Carbon12': [x_c12, y_c12]}

    return data_dict


def get_fitting_filename(exp, channel, dist, Ethr):
    fitfilepath = '/fs/ddn/sdf/group/nexo/users/miaoyu/Sterile_nu/Fits/csv'
    if channel == 'combine':
        fitfile = f'{fitfilepath}/spectralFit_{exp}_source14cm_dist{int(dist*100)}cm_combine_Ethr{int(Ethr*1000)}keV_smeared.csv'
    elif channel == 'ES':
        fitfile = f'{fitfilepath}/spectralFit_{exp}_source14cm_dist{int(dist*100)}cm_{channel}only_Ethr{int(Ethr*1000)}keV_smeared.csv'
    else:
        fitfile = f'{fitfilepath}/spectralFit_{exp}_source14cm_dist{int(dist*100)}cm_{channel}only_smeared.csv'
    label = f'{exp}_{channel}_{int(dist*100)}cm_Ethr={int(Ethr*1000)}keV'
    print(label, fitfile)
    return fitfile, label


def split_label(lb):
    parts = lb.split('_')
    return parts


def load_my_fitting(filename):
    try:
        res = pd.read_csv(filename)
        delta_chisquare = res['delta_chisquare'].to_numpy()
        delta_chisquare = np.reshape(delta_chisquare, (100, 100))
        x, y = get_contour(delta_chisquare, level=4.605)
        return delta_chisquare, x, y
    except Exception as e:
        print(e)
        sys.exit(-1)

def load_my_fitting_parameter(filename, param):
    try:
        res = pd.read_csv(filename)
        data = res[param].to_numpy()
        return data
    except Exception as e:
        print(e)
        sys.exit(-1)


def plot_fitting(filelist, labels, colors, other=False):
    expcls = {'nEXO': "royalblue", 'LZ': 'black', 'XLZD': 'chocolate'}
    chals = {'CC': '-.', 'ES': '--', 'combine': '-'}
    fig, ax = plt.subplots(figsize=(8, 6))
    for filename, lb, c in zip(filelist, labels, colors):
        _, x, y = load_my_fitting(filename)
        exp, cha, _, _ = split_label(lb)
        if c is None:
            ax.plot(x, y, lw=3, linestyle=chals[cha], color=expcls[exp], label=lb)
        else:
            ax.plot(x, y, lw=3, linestyle=chals[cha], color=c, label=lb)
    if other:
        data = load_other_experiment_sensitivities()
        ax.fill(data['Reactor'][0][0], data['Reactor'][0][1], alpha=0.2, color='blueviolet', label='Reactor, JHEP05(2013)050')
        for i in range(1, 4, 1):
            ax.fill(data['Reactor'][i][0], data['Reactor'][i][1], alpha=0.2, color='blueviolet')
        ax.plot(data['LZ_natural'][0], data['LZ_natural'][1], linestyle=chals['ES'], lw=3, color='coral', label='LZ, JHEP11(2014)042')
        
        # For Ga and C12 data, the x-axis is Ue4^2:
        Ue4_square = data['Gallium'][0]
        sin2theta_square = 4 * Ue4_square * (1 - Ue4_square)
        ax.fill(sin2theta_square, data['Gallium'][1], alpha=0.1, color='forestgreen', label='Gallium, JHEP05(2013)050')
        Ue4_square = data['Carbon12'][0]
        sin2theta_square = 4 * Ue4_square * (1 - Ue4_square)
        ax.plot(data['Carbon12'][0], data['Carbon12'][1], ':', lw=3, color='darkred', label='C-12, JHEP05(2013)050')
    
    #ax.hlines(0.55, 1e-2, 1e0, color='gray')
    #ax.vlines(0.2, 1e-2, 1e1, colors='gray')
    #ax.hlines(0.4, 1e-2, 1e0, color='gray')
    
    ax.set_ylabel(r'$\Delta m^2 [\mathrm{eV}^2]$', fontsize=14)
    ax.set_xlabel(r'$\sin^2(2\theta)$', fontsize=14)
    ax.set_xlim(1e-2, 1e0)
    ax.set_ylim(1e-2, 1e1)
    ax.tick_params(labelsize=13)
    ax.loglog()
    ax.legend(fontsize=13)
    plt.tight_layout()
    plt.show()
    return fig
        

def plot_param_histograms(parameters, labels, ylog=True):
    '''
    Input arguments:
    1. parameters (list of 1D parameter arrays)
    2. labels (list of strings)
    3. ylog (bool, default=True, if log scale is used for y-axis)
    Return: fig
    '''
    fig, ax = plt.subplots(figsize=(8, 6))
    for lb, par in zip(labels, parameters):
        ax.hist(par, bins=100, range=(np.min(np.array(parameters)), np.max(np.array(parameters))), label=lb, histtype='step')
    ax.set_xlabel('Fitted values', fontsize=13)
    ax.set_ylabel('', fontsize=0)
    if ylog:
        ax.semilogy()
    ax.legend(fontsize=14)
    ax.tick_params(labelsize=13)
    plt.tight_layout()
    plt.show()
    return fig


