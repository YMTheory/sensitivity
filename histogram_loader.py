import pickle
import histlite as hl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os


def get_oscillation_parameters(dm_square, sin2theta_square):
    dm_square_arr           = np.logspace(-2, 1, 100)
    sin2theta_square_arr    = np.logspace(-2, 0, 100)
    
    if dm_square < dm_square_arr[0]:
        dm_square = dm_square_arr[0]
    elif dm_square > dm_square_arr[-1]:
        dm_square = dm_square_arr[-1]
    else:
        for i in range(len(dm_square_arr)):
            if dm_square_arr[i] <= dm_square < dm_square_arr[i+1]:
                dm_square = dm_square_arr[i]

    if sin2theta_square < sin2theta_square_arr[0]:
        sin2theta_square = sin2theta_square_arr[0]
    elif sin2theta_square > sin2theta_square_arr[-1]:
        sin2theta_square = sin2theta_square_arr[-1]
    else:
        for i in range(len(sin2theta_square_arr)):
            if sin2theta_square_arr[i] <= sin2theta_square < sin2theta_square_arr[i+1]:
                sin2theta_square = sin2theta_square_arr[i]
    
    return dm_square, sin2theta_square



def parse_filename(exp, dist, channel, dm_square, sin2theta_square):
    if channel == "CC" or channel == "ES":
        path = f"/fs/ddn/sdf/group/nexo/users/miaoyu/Sterile_nu/PDFs/{exp}/dist{int(dist*100)}cm/{channel}/histograms/dmsquare{dm_square:.5f}eV2/"
        filename = path + f"MChist_dmsquare{dm_square:.5f}eV2sin2thetasquare{sin2theta_square:.5f}_source14cm_{exp}_dist{int(dist*100)}cm_{channel}_smeared_scaled.p"
    elif channel == "ESbkg" :
        dm_square, sin2theta_square = 0., 0.
        path = f"/fs/ddn/sdf/group/nexo/users/miaoyu/Sterile_nu/PDFs/{exp}/dist{int(dist*100)}cm/{channel}/histograms/dmsquare{dm_square:.5f}eV2/"
        filename = path + f"MChist_dmsquare{dm_square:.5f}eV2sin2thetasquare{sin2theta_square:.5f}_source14cm_{exp}_dist{int(dist*100)}cm_{channel}_smeared_scaled.p"
        
    elif channel == 'EStotal':
        path = f"/fs/ddn/sdf/group/nexo/users/miaoyu/Sterile_nu/PDFs/{exp}/dist{int(dist*100)}cm/ES/histograms/dmsquare{dm_square:.5f}eV2/"
        filename = path + f"MChist_dmsquare{dm_square:.5f}eV2sin2thetasquare{sin2theta_square:.5f}_source14cm_{exp}_dist{int(dist*100)}cm_ES_smeared_total.p"
    
    if not os.path.exists(filename):
        print(f"ERROR: {filename} does not exist!!!")
        return None
    else:
        return filename


def load_histogram(filename):
    try:
        with open(filename, 'rb') as f:
            h = pickle.load(f)
        return h
    except Exception as e:
        print(e)
        return None
    
def parse_oscillation_label(dm_square, sin2theta_square):
    lb = r'$\Delta m^2=$' + f'{dm_square:.5f} eV' + r'$^2$' + r', $\sin^2(2\theta)=$' + f'{sin2theta_square:.5f}'
    return lb
    
        
def plot1D(histograms, labels, projection=-1, logx=False, logy=False):
    fig, ax = plt.subplots(figsize=(8, 6))
    for h, l in zip(histograms, labels):
        if projection == -1:
            hl.plot1d(ax, h, label=l)
        else:
            hl.plot1d(ax, h.project(axes=[projection]), label=l)
    if projection == 0:
        ax.set_xlabel('Baseline [m]', fontsize=14)
        ax.set_ylabel('Count per 3 cm', fontsize=14) 
    elif projection == 1:
        ax.set_xlabel('Electron kinetic energy [MeV]', fontsize=14)
        ax.set_ylabel('Count per 10 keV', fontsize=14)
    ax.tick_params(labelsize=13)
    ax.legend(fontsize=14)
    if logx:
        ax.semilogx()
    if logy:
        ax.semilogy()
    plt.tight_layout()
    plt.show()
    return fig


def plot2D(histogram):
    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(2, 2, width_ratios=[2, 1], height_ratios=[1, 1])

    ax1 = fig.add_subplot(gs[:, 0])
    im = hl.plot2d(ax1, histogram, cmap='viridis',  )
    ax1.set_xlabel('Baseline [m]', fontsize=13)
    ax1.set_ylabel('Electron kinetic energy [MeV]', fontsize=13)
    ax1.tick_params(labelsize=12)
    cb = plt.colorbar(im, ax=ax1)
    cb.set_label('Signal count', fontsize=12)

    ax2 = fig.add_subplot(gs[0, 1])
    hl.plot1d(ax2, histogram.project(axes=[0]), color='blue' )
    ax2.set_xlabel('Baseline [m]', fontsize=11)
    ax2.set_ylabel('Count per 3 cm', fontsize=11)

    ax3 = fig.add_subplot(gs[1, 1])
    hl.plot1d(ax3, histogram.project(axes=[1]), color='darkorange' )
    ax3.set_xlabel('Electron kinetic energy [MeV]', fontsize=11)
    ax3.set_ylabel('Count per 10 keV', fontsize=11)

    plt.tight_layout()
    plt.show()
    #plt.savefig('./plots/BaselineEnergy_distribution_14cmsource_1cmsmear_nonoscillation_ES.eps')    
    return fig
       
    





