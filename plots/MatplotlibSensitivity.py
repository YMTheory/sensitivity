import numpy as np
import matplotlib
#if matplotlib.get_backend() != "MacOSX": matplotlib.use('QT4Agg')
import matplotlib.pyplot as plt
import math

matplotlib.rcParams['font.family'] = 'Times New Roman'
build_KKDC = True
show_individual_NMEs = False ## plot NME lines separately, False for a single filled band
show_limits_sens = False ## plot the limits and sensitivities as arrows


ylims = [1e24, 2e26]
xlims = ylims

def AxisToWorld(x):
    upper_limit = math.log10(xlims[1])
    lower_limit = math.log10(xlims[0])
    return (math.log10(x)-lower_limit)/(upper_limit - lower_limit)

DistTextToLine = 0.01
DistTextFromEnd = 0.015
ArrowHeadSize = 0.03

def DrawArrowOnLogLog(Vertical, StartPos, RealEnd, FracPos, color, ax_arr, label):
    # Vertical is true or false.
    # StartPos is the axis coordinate we wish to indicate.
    # All other coordinates are relative to the plot size: FracLength, and FracPos.
    def_args = {"length_includes_head" : True, "width" : .005,
                "head_width" : ArrowHeadSize, "head_length" : ArrowHeadSize}
    awidth = 0.015
    thepos = AxisToWorld(StartPos)
    FracLength = AxisToWorld(RealEnd) - thepos
    if Vertical:
        ax_arr.arrow(FracPos, thepos,
                        0, FracLength,
                        color = color, **def_args)
        ax_arr.plot([FracPos - awidth, FracPos + awidth], [thepos, thepos], '-', color=color, lw=3)
    else:
        ax_arr.arrow(thepos, FracPos,
                        FracLength, 0,
                        color = color, **def_args)
        ax_arr.plot([thepos]*2, [FracPos - awidth, FracPos + awidth], color=color, lw=3)
    ax_arr.set_ylim([0, 1])
    ax_arr.set_xlim([0, 1])

    # Place label so that it is uniformly offset from the arrow head.
    if Vertical:
        ax_arr.text(FracPos - DistTextToLine, thepos + FracLength - ArrowHeadSize - DistTextFromEnd, label,
                    rotation = "vertical", ha = 'right', va = 'top',
                    color = color, backgroundcolor="white")
    else:
        #ax_arr.text(thepos + FracLength - ArrowHeadSize - DistTextFromEnd, FracPos + DistTextToLine, label, rotation = "horizontal", ha = 'right', va = 'bottom', color = color)
        ax_arr.text(0.99, FracPos + DistTextToLine*2.2, label, rotation = "horizontal", ha = 'right', va = 'bottom', color = color, backgroundcolor="white")

def draw_axis(ax, x, y, z, fmt=dict(), txt_fmt=dict(), ignore=[],
              label = "", label_pos = 0.28, label_above = True, return_band=False):
    xlow, xhi = x
    ylow, yhi = y
    zlow = z
    allx = np.arange(xlow, xhi, xlow/10000.)
    ally = (yhi - ylow)/(xhi - xlow) * (allx - xlow) + ylow
    allz = zlow*math.sqrt(xlow)/np.sqrt(allx) 
    if( return_band ):
        return allx, ally
    
    f = lambda xx: np.power(10, np.log10(xx) - np.floor(np.log10(xx)))
    integers = f(allz) % 1 < 1e-3

    ax.plot(allx, ally, '-', **fmt)
    already_seen_values = []
    for (X, Y, Z) in zip(allx[integers], ally[integers], allz[integers]):
        if len(already_seen_values) > 0 and np.abs(np.array(already_seen_values) - Z).min() < 1e-3: continue 
        tick_label = '{:.1f}'.format(Z)
        if Z < 0.1:
            tick_label = '{:.2f}'.format(Z)
        #if int(Z) != 0: tick_label = '%i ' % (int(Z))
        already_seen_values.append(Z)
        ax.annotate('', xy=(X,Y), xytext=(5, -5), ha='right',
                    textcoords='offset points', 
                    arrowprops=dict(arrowstyle='-', shrinkA=0, shrinkB=0, **fmt), **txt_fmt)
        if tick_label in ignore: continue
        ax.annotate(tick_label, xy=(X,Y), xytext=(-1, 1), ha='right',
                    textcoords='offset points', **txt_fmt) 

    # Place the label.
    x_for_label = pow(zlow*math.sqrt(xlow)/label_pos, 2)
    y_for_label = (yhi - ylow)/(xhi - xlow) * (x_for_label - xlow) + ylow
    if label_above: y_for_label *= 1.15
    else: y_for_label *= 0.88
    ax.text(x_for_label, y_for_label, label, ha = 'left', va = 'bottom', rotation = 45, color = fmt['color'])

fig, ax = plt.subplots(figsize = (7,7))
label_size = 20
ax.set_ylabel(r"T$_{1/2}$ $^{76}$Ge (yr)", fontsize=label_size)
ax.set_xlabel(r"T$_{1/2}$ $^{136}$Xe (yr)", fontsize=label_size)
plt.setp(ax.get_xticklabels(), fontsize=label_size) 
plt.setp(ax.get_yticklabels(), fontsize=label_size) 
plt.tight_layout()

GERDAColor = "k"
KZColor = "red"
EXOColor = "blue"

ax.set_xlim(xlims)
ax.set_ylim(ylims)
ax.set_yscale("log")
ax.set_xscale("log")

ax_world = fig.add_axes(ax.get_position()) # For linear coordinates, needed to make distances reasonable.
ax_world.set_axis_off()

# Klapdor claim
if build_KKDC and show_limits_sens:
    ax.fill_between(xlims, 1.92e25,2.67e25, color='0.875')
    ax_world.text(0.02, AxisToWorld(1.97e25), "KK&K 68% CL", color='0.3', verticalalignment="bottom",
      horizontalalignment="left")

if show_limits_sens:
    ax_world.plot([AxisToWorld(1.9e25)]*2, [0,1], '-', lw=2, color=EXOColor)
    ax_world.text(AxisToWorld(1.9e25) - DistTextToLine, DistTextFromEnd, "EXO-200 Sensitivity",
                  rotation="vertical", verticalalignment="bottom", horizontalalignment="right",
                  color=EXOColor)

    ax_world.plot([AxisToWorld(5.7e25)]*2, [0,1], '--', lw=2, color=EXOColor)
    ax_world.text(AxisToWorld(5.7e25) + DistTextToLine, DistTextFromEnd, "EXO-200 Projected\nSensitivity (Phase I+II)",
                  rotation="vertical", verticalalignment="bottom", horizontalalignment="left",
                  color=EXOColor)

    ax_world.plot([AxisToWorld(4.9e25)]*2, [0,1], '-', lw=2, color=KZColor)
    ax_world.text(AxisToWorld(4.9e25) - 7*DistTextToLine, DistTextFromEnd, "KamLAND-Zen Sensitivity\n(Phase II only)",
                  rotation="vertical", verticalalignment="bottom", horizontalalignment="left",
                  color=KZColor)

    ax_world.plot([0,1], [AxisToWorld(2.4e25)]*2, '-', color=GERDAColor , lw=2)

    ax_world.text(AxisToWorld(0.7e26) , AxisToWorld(2.5e25), "GERDA\nSensitivity",
                  rotation="horizontal", verticalalignment="bottom", horizontalalignment="left",
                  color=GERDAColor)

    # arrows
    exo_vertical = 0.38
    DrawArrowOnLogLog(True, 2.1e25, 0.7*ylims[1], AxisToWorld(0.9*xlims[1])-ArrowHeadSize/2, GERDAColor, ax_world, "GERDA Limit")
    DrawArrowOnLogLog(False, 1.1e25, 0.9*xlims[1], 0.36, EXOColor, ax_world, "EXO-200 Limit")
    DrawArrowOnLogLog(False, 1.1e26, 0.9*xlims[1], 0.46, KZColor, ax_world, "KamLAND-Zen\nLimit (Phase I+II)")

# OLD VALUES
#(1.618948e+25,1e+26,1e+24,6.176851e+24,0.5593971,"GCM", "darkorange", [], 0.503, True),
#(2.22198e+25,1e+26,1e+24,4.50049e+24,1.072816,"NSM", "firebrick", ["1.0"], 0.668, False),
#(3.585184e+25,1e+26,1e+24,2.789257e+24,0.7013337,"IBM-2", "purple",["0.7"], 0.39, True),
#(6.00272e+25,1e+26,1e+24,1.665911e+24,0.9141899,"RQRPA-1", "darkolivegreen", ["0.9"], 0.408, True),
# NEW VALUES, uses new 
# (1.94e25,1e26,1e24,5.144e24,0.625,"GCM", "darkorange", ["0.6", "0.1"], 0.503, True),
# (2.74e25,1e26,1e24,3.654e24,1.482,"NSM", "firebrick", ["0.2"], 0.668, False),
# (4.29e25,1e26,1e24,2.33e24,0.788,"IBM-2", "purple",["0.09"], 0.39, True),
# (8.31e25,1e26,1e24,1.203e24,1.067,"RQRPA", "darkolivegreen", ["1.0"], 0.408, True),
# (9.48767e+25,1e+26,1e+24,1.054e+24,0.1526218,1.486607,"QRPA-2", "rosybrown", [], 1.3, False),

# NEW VALUES (from document by Petr Vogel, 2015/02/17)   
# (1.95e25,1e26,1e24,5.14e24,0.625,"GCM", "darkorange", ["0.6", "0.1"], 0.503, True),
# (2.67e25,1e26,1e24,3.74e24,1.2,"NSM", "firebrick", ["0.2"], 0.668, False),
# (4.29e25,1e26,1e24,2.33e24,0.703,"IBM-2", "purple",["0.09"], 0.39, True),
# (8.55e25,1e26,1e24,1.17e24,1.29,"QRPA", "darkolivegreen", ["1.0"], 0.408, True),

# NEW VALUES as of 2016/05/29.  Uses the following references (same as KZ paper,
# arXiv:1605.02889v1):
# EDF: Mxe = 4.20, Mge = 4.60 PRL 105, 252503 (2010)
# ISM: Mxe = 2.19, Mge = 2.81 Nucl Phys A 818, 139 (2009)
# IBM-2: Mxe = 3.05, Mge = 4.68 PRC 91, 034304 (2015)
# QRPA: Mxe = 2.02, Mge = 4.64 PRC 89, 064308 (2014)
# Skyrme QRPA: Mxe = 1.55, Mge = 5.09 PRC 87, 064302 (2013)
Gxe = 1./6.88e24 ## from document by Petr Vogel, 2015/02/17
Gge = 1./4.24e25 ## from document by Petr Vogel, 2015/02/17
minx, maxy = xlims[0], ylims[1]

## order is M_Xe, M_Ge, Color, Tick labels to ignore, Label position, Label Above?
NME_dict = {"EDF": [4.20,4.60,"darkorange",["0.6","0.3","0.1"],0.503, True],
            "ISM": [2.19,2.81,"firebrick",["0.9","0.7","0.5"],0.68, False],
            "IBM-2": [3.05,4.68,"purple",["0.8","0.4","0.3","0.1"],0.39, False],
            "QRPA": [2.02,4.64,"darkolivegreen",["0.4","0.3","0.1"],0.498, True],
            "SkyrmeQRPA": [1.55,5.09,"dodgerblue",["1.0","0.5","0.3"],1.38, True],
}
NME_list = []
for nme in ["EDF", "ISM","IBM-2","QRPA","SkyrmeQRPA"]:
    if( nme not in NME_dict ):
        print "Warning could not find %s, skipping" % nme
        continue
    
    prop_fac = Gxe/Gge * (NME_dict[nme][0]/NME_dict[nme][1])**2
    masslo = np.sqrt( 1./(NME_dict[nme][0]**2 * Gxe * minx) )
    NME_list.append( (maxy/prop_fac, maxy, minx, minx*prop_fac, masslo, nme, NME_dict[nme][2], NME_dict[nme][3], NME_dict[nme][4], NME_dict[nme][5] ) )
print NME_list

if show_individual_NMEs:

    for xhi, yhi, xlo, ylo, zlo,label,c,ignore,label_pos, label_above in NME_list:
        print xhi, yhi, xlo, ylo, zlo
        if xhi < xlims[1]:
            xhi = xlims[1] 
            yhi = xlims[1]*ylo/xlo 
        draw_axis(ax, (xlo, xhi), (ylo, yhi), (zlo), 
          ignore=ignore, fmt={"color" : c, "linewidth" : 2}, txt_fmt={"color" : c},
          label = label, label_pos = label_pos, label_above = label_above) 

## plot just a band at the min max instead        
else:
    ## find min and max nme:
    min_nme = ["",1e10]
    max_nme = ["",0]
    for k in NME_dict:
        if( NME_dict[k][0]/NME_dict[k][1] > max_nme[1] ):
            max_nme = [k,NME_dict[k][0]/NME_dict[k][1]]
        if( NME_dict[k][0]/NME_dict[k][1] < min_nme[1] ):
            min_nme = [k,NME_dict[k][0]/NME_dict[k][1]]    

    print "Min NME:", min_nme
    print "Max NME:", max_nme
    
    for xhi, yhi, xlo, ylo, zlo,label,c,ignore,label_pos, label_above in NME_list:
        if xhi < xlims[1]:
            xhi = xlims[1] 
            yhi = xlims[1]*ylo/xlo
        if( label == min_nme[0] ): 
            xmin,ymin = draw_axis(ax, (xlo, xhi), (ylo, yhi), (zlo), return_band=True) 
        if( label == max_nme[0] ):
            xmax,ymax = draw_axis(ax, (xlo, xhi), (ylo, yhi), (zlo), return_band=True) 
    xboth = np.linspace(np.min(np.hstack((xmin,xmax))),np.max(np.hstack((xmin,xmax))),1e3)
    ymin_both = np.interp( xboth, xmax, ymax )
    ymax_both = np.interp( xboth, xmin, ymin )
    ymin_both[ymin_both < ylims[0]] = ylims[0]
    ymax_both[ymax_both > ylims[1]] = ylims[1]

    ax.fill_between(xboth, ymin_both, ymax_both, edgecolor='none', color='darkolivegreen', alpha=0.2)
    
    for xhi, yhi, xlo, ylo, zlo,label,c,ignore,label_pos, label_above in NME_list:
        if xhi < xlims[1]:
            xhi = xlims[1] 
            yhi = xlims[1]*ylo/xlo
        if( label == min_nme[0] ):
            draw_axis(ax, (xlo, xhi), (ylo, yhi), (zlo), 
                      ignore=ignore, fmt={"color" : 'darkolivegreen', "linewidth" : 2}, txt_fmt={"color" : 'darkolivegreen'},
                      label = label, label_pos = label_pos, label_above = True) 
        if( label == max_nme[0] ):
            draw_axis(ax, (xlo, xhi), (ylo, yhi), (zlo), 
                      ignore=ignore, fmt={"color" : 'darkolivegreen', "linewidth" : 2}, txt_fmt={"color" : 'darkolivegreen'},
                      label = label, label_pos = label_pos-0.03, label_above = False) 
                
NMEsuffix = "" if show_individual_NMEs else "_band"
limssuffix = "" if show_limits_sens else "_nolims"
plt.savefig('Sensitivity_GeXe'+NMEsuffix+limssuffix+'.pdf')
plt.show()

