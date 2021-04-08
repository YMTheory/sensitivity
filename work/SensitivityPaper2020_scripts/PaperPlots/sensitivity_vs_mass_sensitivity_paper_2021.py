import glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import cPickle as pickle

matplotlib.rcParams['font.family'] = 'Times New Roman'
#matplotlib.rcParams['text.usetex'] = True

plot_vs_total_mass = False  ## True to plot vs m_tot, False to plot vs m_min
five_or_ten_yr = False ## True for 5 yr, False for 10 yr

allowed_regions = pickle.load(open("allowed_regions.pkl","rb"))

if( plot_vs_total_mass):
    xlim = [1e-2,1]
    idx_to_use = [2,3]
    xticklist = [1e-2, 1e-1]
    suffix = '_mtot'
    xlab = "$m_{\mathrm{tot}}$ [eV]"
else:
    xlim = [1e-4,1]
    idx_to_use = [0,1]
    xticklist = [1e-4, 1e-3, 1e-2, 1e-1]
    suffix = '_mmin'
    xlab = "$m_{\mathrm{min}}$ [eV]"
    
fig,(ax1,ax2) = plt.subplots(1,2,sharey=True)
for f in allowed_regions[idx_to_use[0]]:
    ax1.loglog( f[0], f[1], 'r'+f[2], linewidth=1.5 )

ax1.set_xlim(xlim)
ax1.set_ylim([1e-3,1])
ax1.set_xticks(xticklist)

for f in allowed_regions[idx_to_use[1]]:
    ax2.loglog( f[0], f[1], 'b'+f[2], linewidth=1.5 )
    
fig.subplots_adjust(wspace=0,bottom=0.12,top=0.95,right=0.95)

ax2.set_xlim(xlim)
ax2.set_ylim([1e-3,1])
ax2.set_xticks(xticklist + [1,])
ax2.tick_params(axis='both',which='major',labelsize='10')

ax1.set_xlabel( xlab, fontsize=18)
ax2.set_xlabel( xlab, fontsize=18)
ax1.xaxis.set_label_coords(0.5,-0.08)
ax2.xaxis.set_label_coords(0.5,-0.08)
ax1.set_ylabel( r"$<m_{\beta\beta}>$ [eV]", fontsize=18)
ax1.yaxis.set_label_coords(-0.18,0.5)
ax1.tick_params(axis='both',which='major',labelsize='14')
ax2.tick_params(axis='both',which='major',labelsize='14')

ax1.set_xlabel( xlab ) 
ax2.set_xlabel( xlab ) 
ax1.set_ylabel( r"$<m_{\beta\beta}>$ [eV]" ) 

fig.set_size_inches(6,4.5)

plt.savefig("sens_allowed_only"+suffix+".pdf")

if( plot_vs_total_mass ):
    ax1.fill_between([0.23, 1],[1e-3,1e-3],[1,1], color='k', alpha= 0.15, edgecolor='none')
    ax2.fill_between([0.23, 1],[1e-3,1e-3],[1,1], color='k', alpha= 0.15, edgecolor='none')
#    ax2.text(0.4,1.2e-3, "Planck,\narXiv:1502.01589",rotation="vertical",ha='left',va='bottom')
    ax2.text(0.4,1.2e-3, "Planck 2015 Results,\nA&A 594, A13 (2016)",rotation="vertical",ha='left',va='bottom')

    ax1.annotate('N.O.', xy=(0.25,0.055), color='red', size=22, fontweight='bold')
    ax2.annotate('I.O.', xy=(0.25,0.055), color='blue', size=22, fontweight='bold')
else:
    ax1.annotate('N.O.', xy=(0.0008,0.002), color='red', size=22, fontweight='bold')
    ax2.annotate('I.O.', xy=(0.001,0.025), color='blue', size=22, fontweight='bold')
            
    

exo200_final = 3.5e25 ## half-life sensitivity from final paper arXiv: 1906.02723 - Phys. Rev. Lett. 123, 161802 (2019)
## old value - 2019 sensitivity paper
#	nexo_10yrs_old = 9.06e27 ## half-life limit in 10 yrs (yr)
## new value - 2021 sensitivity paper
nexo_10yrs = 1.35e28 ## half-life limit in 10 yrs (yr)

mmin = 1.11 ## minimum matrix element, Phys. Rev. C 97, 045503 (2018)
mmax = 4.77 ## maximum matrix element, PhysRevC.91.024316 (2015)
Gxe = 1/6.88e24 ## from Benato paper arXiv:1705.02996 and PhysRevC.85.034316
lmin_exo200, lmax_exo200 = np.sqrt(1./(mmin**2 * Gxe * exo200_final)), np.sqrt(1./(mmax**2 * Gxe * exo200_final))
lmin_10yrs, lmax_10yrs = np.sqrt(1./(mmin**2 * Gxe * nexo_10yrs)), np.sqrt(1./(mmax**2 * Gxe * nexo_10yrs))

## "one-sigma" band
mmin = 1.55 ## minimum matrix element, 1 sigma band
mmax = 4.32 ## maximum matrix element, 1 sigma band 
lmin_10yrs_1sigma, lmax_10yrs_1sigma = np.sqrt(1./(mmin**2 * Gxe * nexo_10yrs)), np.sqrt(1./(mmax**2 * Gxe * nexo_10yrs))
lmin_exo200_1sigma, lmax_exo200_1sigma = np.sqrt(1./(mmin**2 * Gxe * exo200_final)), np.sqrt(1./(mmax**2 * Gxe * exo200_final))

print('Mass range(meV): ', lmin_10yrs*1000, lmax_10yrs*1000)

## print sensitivity for all values
#for val_ in [3.6, 1.55, 2.46, 2.54, 1.11, 4.77, 4.32, 2.28, 2.45]:
#	print(val_, np.sqrt(1./(val_**2 * Gxe * nexo_10yrs))*1000)

## exo200 current
col=np.array([102.,178.,255.])/256.
ax1.fill_between(xlim,[lmax_exo200, lmax_exo200],[lmin_exo200, lmin_exo200],edgecolor=col,facecolor=col,alpha=0.5)
ax2.fill_between(xlim,[lmax_exo200, lmax_exo200],[lmin_exo200, lmin_exo200],edgecolor=col,facecolor=col,alpha=0.5)

ax1.fill_between(xlim,[lmax_exo200_1sigma, lmax_exo200_1sigma],[lmin_exo200_1sigma, lmin_exo200_1sigma],edgecolor='blue',facecolor='blue',alpha=0.3)
ax2.fill_between(xlim,[lmax_exo200_1sigma, lmax_exo200_1sigma],[lmin_exo200_1sigma, lmin_exo200_1sigma],edgecolor='blue',facecolor='blue',alpha=0.3)


ax1.text( xlim[0]*1.15, 0.75*np.mean([lmin_exo200, lmax_exo200]), r"EXO-200", ha="left", va="center", fontsize=20)
ax1.text((xlim[0]*1.2,xlim[0]*6.0)[plot_vs_total_mass], 0.5*np.mean([lmin_exo200, lmax_exo200]), r"PRL 123, 161802 (2019)", ha="left", va="center", fontsize=13, fontstyle='italic', fontweight='bold')


## nexo
col=np.array([180.,50.,150.])/256.
ax1.fill_between(xlim,[lmax_10yrs, lmax_10yrs],[lmin_10yrs, lmin_10yrs],edgecolor='lightgreen',facecolor='lightgreen',alpha=0.7)
ax2.fill_between(xlim,[lmax_10yrs, lmax_10yrs],[lmin_10yrs, lmin_10yrs],edgecolor='lightgreen',facecolor='lightgreen',alpha=0.7)

ax1.fill_between(xlim,[lmax_10yrs_1sigma, lmax_10yrs_1sigma],[lmin_10yrs_1sigma, lmin_10yrs_1sigma],edgecolor='green',facecolor='green',alpha=0.7)
ax2.fill_between(xlim,[lmax_10yrs_1sigma, lmax_10yrs_1sigma],[lmin_10yrs_1sigma, lmin_10yrs_1sigma],edgecolor='green',facecolor='green',alpha=0.7)

ax1.text( 4e-2, 8.2e-3, r"nEXO", ha="left", va="center", fontsize=22, fontweight='bold')
ax2.text( 1.1e-4, 8.2e-3, r"10 Years", ha="left", va="center", fontsize=22, fontweight='bold')

plt.savefig("sens_nexo"+suffix+"_%dyrs_v20.pdf"%((10,5)[five_or_ten_yr]))
plt.savefig("sens_nexo"+suffix+"_%dyrs_v20.png"%((10,5)[five_or_ten_yr]))

plt.show()
