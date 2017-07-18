
import ROOT

import cPickle as pickle
import glob, os
import matplotlib.pyplot as plt
import numpy as np

import matplotlib
import matplotlib.patches as patches

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'

def MakePlots(outname,logscale=False):
    markerList = {'fwhm':'bo-','1-sigma':'gD-','2-sigma':'rs-','3-sigma':'m^-'}

    axis_label_fontsize = 15
    legend_fontsize = 14
    labels_fontsize = 14

    min_x, max_x = 0.25, 2.25
    if logscale:
        min_x, max_x = 0.4, 3.0
    min_y, max_y = 1e-7, 1e-2
    # exo-like resolution function, Gaussian assumption
    results= {0.005: {'1-sigma': 0.00012156702423737897, '3-sigma': 0.0033944465261427914, 'fwhm': 0.0001692518237632612, 'klz': 6.297484902155599, '2-sigma': 0.0007526356585967342}, 0.025: {'1-sigma': 1.3823112294776365, '3-sigma': 40.42842561928248, 'fwhm': 1.990805720575736, 'klz': 21.588934716815682, '2-sigma': 8.782721038129239}, 0.015: {'1-sigma': 0.07094488432994694, '3-sigma': 2.174251918499273, 'fwhm': 0.10469339584233239, 'klz': 10.598819211820723, '2-sigma': 0.4668243877624718}, 0.01: {'1-sigma': 0.0067560696852524416, '3-sigma': 0.20346954650063775, 'fwhm': 0.009434754255835287, 'klz': 7.8150010856290635, '2-sigma': 0.04540726057558686}, 0.047: {'1-sigma': 46.97428731061518, '3-sigma': 1285.7508966185196, 'fwhm': 67.73201160365716, 'klz': 96.1696785929671, '2-sigma': 295.85924841600354}, 0.029: {'1-sigma': 3.227244202804286, '3-sigma': 91.78470542957535, 'fwhm': 4.545140399015509, 'klz': 28.837731118066934, '2-sigma': 20.473612865109317}, 0.02: {'1-sigma': 0.3891544728685403, '3-sigma': 11.404611063743879, 'fwhm': 0.5427439321938436, 'klz': 15.018502268045767, '2-sigma': 2.5010184130460402}}
    # roiList = [roi for roi in results[0.005]]
    roiList = ['fwhm']

    xList = [resol*100 for resol in sorted(results)]
    yList = {}
    for roi in roiList:
        yList[roi] = [results[resol][roi]/1000 for resol in sorted(results)]

    fig, ax = plt.subplots()

    # plt.plot( xList, yList['1-sigma'], markerList['1-sigma'], label=r'$\pm$ $1\sigma$')
    plt.plot( xList, yList['fwhm'], markerList['fwhm'], label='FWHM')
    # plt.plot( xList, yList['2-sigma'], markerList['2-sigma'], label=r'$\pm$ $2\sigma$')
    # plt.plot( xList, yList['3-sigma'], markerList['3-sigma'], label=r'$\pm$ $3\sigma$')
    
    plt.xlabel(r"Resolution $\sigma(E)/E$ at the Q-value [%]", fontsize=axis_label_fontsize)
    # plt.ylabel(r"$^{136}$Xe  $2\nu\beta\beta$ Background [cts/(tonne$\cdot$yr)]", fontsize=axis_label_fontsize)
    plt.ylabel(r"$^{136}$Xe  $2\nu\beta\beta$ [cts/(FWHM$\cdot$kg$\cdot$yr)]", fontsize=axis_label_fontsize)
    plt.xlim([min_x, max_x])
    plt.ylim([min_y, max_y])
    
    if logscale:
        plt.xscale('log')

        ax.xaxis.set_minor_formatter(matplotlib.ticker.FormatStrFormatter("%.1f"))
        ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.1f"))
        ax.xaxis.set_ticklabels(['','','','','','','','','','','','0.5','','','','','2.0','3.0','4.0','5.0','',''],minor=True)
        

    plt.yscale('log')

    # plt.legend(loc="lower right",title="ROI",prop={'size':legend_fontsize}) #,frameon=False
    # ax.get_legend().get_title().set_fontsize(legend_fontsize)


    ax.grid(True, which='both', linestyle='dotted')
    fig.set_size_inches(7.5,5.5)

    ax.add_patch(patches.Polygon([[0.9, min_y], [0.9, max_y], [1.1, max_y], [1.1, min_y]], alpha=0.3, color='grey', lw=0, fill=True))

    plt.savefig("pdf/" + outname + '.pdf', bbox_inches='tight' )
    plt.savefig("png/" + outname + '.png', bbox_inches='tight' )
    # plt.show()


if __name__ == "__main__":

    MakePlots("bb2n_bkg_vs_resol_v3_log", logscale=True)
