#!/usr/bin/python
#-*- coding: UTF-8 -*-
import os
import sys
sys.path.append("/nashome2/hejajama/lib/")
sys.path.append("/home/hejajama/lib/")
sys.path.append("/home/heikki/lib")
sys.path.append("/Users/heikki/lib")
import math
from math import pow
from matplotlibhelper import *
import pylab
import scipy.integrate
import numpy as np

#import scipy.integrate

from matplotlib import rcParams
#rcParams.update({'figure.autolayout': True})

slides=True

ShowStatErrs = True
ShowBand = False  # show largest/smallest xs at given t

rc('text',usetex=True)
rc('text.latex',  preamble='\usepackage{amsmath},\usepackage{amssymb},\usepackage{mathtools}')
textsize=19

rc("xtick", labelsize=textsize)
rc("ytick", labelsize=textsize)

minx=0.0
maxx=2.50
miny=0.02
maxy=1.5e3

#widths
lw_coh = 2.5
lw_incoh = 0.9

if slides:
    lw_incoh += 0.7

markersize_coh = 5.0
markersize_incoh = 6.0

PI=3.141592

GEVSQRTONB = 1.0e7/(5.068*5.068)



fig = figure()
p1=fig.add_subplot(111)
#xlabel(r"$r$ $[\mathrm{GeV}^{-1}]$ ", fontsize=textsize+2)
xlabel(r"$|t|$ $[\mathrm{GeV}^2]$", fontsize=textsize+2)
#xlabel(r"$rQ_s$")
ylabel(r"$\mathrm{d}\sigma/\mathrm{d}t$ $[\mathrm{nb}/\mathrm{GeV}^2]$ ", fontsize=textsize+2)

    
#fig.suptitle(r"$W = 100 \, \mathrm{GeV}, Q^2=0\,\mathrm{GeV}^2$" , fontsize=textsize )

# fname title style normalization,facecolor, hatch, linecolor
files = [
         #["final/ipsat2012_w_100_q2_0", r"IPsat, $B_p=4.0\, \mathrm{GeV}^{-2}$", Linestyle(0), "", "", 1.0, "black"],
         
         
         # B_p = 3.5
         #["final/ipsat2012_bp_3.5_bq_0.5_w_100_q2_0", r"$B_{qc}=3.5\,\mathrm{GeV}^{-2}, B_q=0.5\,\mathrm{GeV}^{-2}$", Linestyle(1), 'blue', "", 1.0, "blue"],
         #["final/ipsat2012_bp_3.5_bq_1.0_w_100_q2_0", r"$B_{qc}=3.5\,\mathrm{GeV}^{-2}, B_q=1.0\,\mathrm{GeV}^{-2}$", Linestyle(1), 'red', "", 1.0, "red"],
         #["paper_2/ipsat_bp_3.5_bq_1.0_new", r"$B_{qc}=3.5\,\mathrm{GeV}^{-2}, B_q=1.0\,\mathrm{GeV}^{-2}$ new", Linestyle(1), 'blue', "", 1.0, "blue"],
         #["paper_2/ipsat_bp_3.5_bq_1.0_new2", r"$B_{qc}=3.5\,\mathrm{GeV}^{-2}, B_q=1.0\,\mathrm{GeV}^{-2}$ new2", Linestyle(3), 'red', "", 1.0, "red"],
         
         ##["paper_2/ipsat_bp_3.5_bq_1.0_origin2", r"$B_{qc}=3.5\,\mathrm{GeV}^{-2}, B_q=1.0\,\mathrm{GeV}^{-2}$ origin", Linestyle(1), 'blue', "", 1.0, "blue"],
         
         ## B_p = 1.0
         #["final/ipsat2012_bp_1.0_bq_3.0_w_100_q2_0", r"$B_{qc}=1.0\,\mathrm{GeV}^{-2}, B_q=3.0\,\mathrm{GeV}^{-2}$", Linestyle(3), "blue", "", 1.0, "blue"],  # paper: grey black

         
         

         
         #fluct 3.5 1.0
         #["final/ipsat2012_bp_3.5_bq_1.0_w_100_q2_0_gridfluct_05", r"IPsat, $a=0.5\,\mathrm{fm}$", Linestyle(2),1.0],
         #["final/ipsat2012_bp_3.5_bq_1.0_w_100_q2_0_gridfluct_04", r"IPsat, $a=0.4\,\mathrm{fm}$", Linestyle(3), "blue", "", 1.0, Color(0)],
         #["final/ipsat2012_bp_3.5_bq_1.0_w_100_q2_0_quarkfluct_0.5", r"$B_{qc}=3.5\,\mathrm{GeV}^{-2}, B_q=1.0\,\mathrm{GeV}^{-2}, \sigma=0.5$", Linestyle(2), "blue", "", 1.0, "blue"],
         #["paper_2/ipsat_bp_4.0_gridfluct_04", r"IPsat $B_p=4\,\mathrm{GeV}^{-2}, \sigma=0.5$", Linestyle(2), "blue", "", 1.0, "blue"],
         
         #
         
         # round ipsat, fluctuations
         #["paper_2/ipsat_bp_4.0_gridfluct_04", r"$B_{p}=4.0\,\mathrm{GeV}^{-2}, a=0.4\,\mathrm{fm}, \sigma=0.5$", Linestyle(3), "grey", "", 1.0, Color(1)],
         
         # ipglasma
         #["final/ipglasma_bp_2.0_bq_0.3_m04_n07_ncf_416_w_100_q2_0", r"$B_{qc}=2.0\,\mathrm{GeV}^{-2}, B_{q}=0.3\,\mathrm{GeV}^{-2}$", Linestyle(0), 'grey', "", 1.0, "black"],
         #["final/ipglasma_bp_4.0_m04_n065_ncf_288_w_100_q2_0", r"$B_{p}=4\,\mathrm{GeV}^{-2}, m=0.4\,\mathrm{GeV}$", Linestyle(0), Color(0), "", 1.0, "black" ],
         #["paper_2/ipglasma_bp_2.0_bq_0.3_m04_n075_qsfluct_ncf_96", r"$B_p=2.0\,\mathrm{GeV}^{-2}, B_q=0.3\,\mathrm{GeV}^{-2}$, fluct", Linestyle(2), Color(2), "", 1.0, Color(2) ],
         #["paper_2/ipglasma_bp_2.0_bq_0.5_m08_n075_qsfluct_ncf_96", r"$B_p=2.0\,\mathrm{GeV}^{-2}, B_q=0.5\,\mathrm{GeV}^{-2}, m=0.8$", Linestyle(2), Color(2), "", 1.0, Color(2) ],
         #["paper_2/ipglasma_bp_2.0_bq_0.5_m04_n075_qsfluct_ncf_384", r"$B_p=2.0\,\mathrm{GeV}^{-2}, B_q=0.5\,\mathrm{GeV}^{-2}, m=0.4$", Linestyle(3), Color(3), "", 1.0, Color(3) ],
         #["paper_2/ipglasma_bp_2.0_bq_0.5_m02_n075_qsfluct_ncf_96", r"$B_p=2.0\,\mathrm{GeV}^{-2}, B_q=0.5\,\mathrm{GeV}^{-2}, m=0.2$", Linestyle(1), Color(1), "", 1.0, Color(2) ],
         #["paper_2/ipglasma_bp_4.0_ny_1", r"$B_p=4.0\,\mathrm{GeV}^{-2}, N_y=1$", Linestyle(1), Color(1), "", 1.4, Color(2) ],
         
         #["paper_2/ipglasma_bp_4.0_mdep_m_0.2", r"$B_p=4.0\,\mathrm{GeV}^{-2}, m=0.2\,\mathrm{GeV}$", Linestyle(1), Color(0), "", 1.4, Color(1) ],
         #["paper_2/ipglasma_bp_4.0_mdep_m_0.4", r"$B_p=4.0\,\mathrm{GeV}^{-2}, m=0.4$", Linestyle(0), Color(2), "", 1.4, Color(2) ],
         #["paper_2/ipglasma_bp_4.0_mdep_m_0.6", r"$B_p=4.0\,\mathrm{GeV}^{-2}, m=0.6$", Linestyle(3), Color(3), "", 1.4, Color(3) ],
         #["paper_2/ipglasma_fluct_m_0.2", r"$m=0.2$", Linestyle(3), Color(3), "", 1.4, Color(3) ],
         #["paper_2/ipglasma_fluct_m_0.6", r"$m=0.6$", Linestyle(1), Color(1), "", 1.4, Color(1) ],
         
         # #fluxtube
         #["paper_2/ipsat_bp_3.5_bq_0.5_fluxtube_norm_0.11", r"Fluxtube $B_{qc}=3.5, B_q=0.5$", Linestyle(0), Color(1), "", 1.0, Color(1)],
         #["paper_2/ipsat_bp_3.5_bq_1.0_fluxtube_norm_0.11", r"Fluxtube $B_{qc}=3.5, B_q=1.0$", Linestyle(1), Color(2), "", 1.0, Color(2)],
         #["paper_2/ipsat_bp_3.5_bq_1.5_fluxtube_norm_0.11", r"flix $B_p=3.5, B_q=1.5$", Linestyle(3), Color(3), "", 1.0, Color(3)],
 
         
         # w = 75
         ["paper_2/ipsat_bp_1.0_bq_3.0_w_75_q2_0", r"$B_{qc}=1.0\,\mathrm{GeV}^{-2}, B_q=3.0\,\mathrm{GeV}^{-2}$", Linestyle(2), 'red', "", 1.0, "red"],
         ["paper_2/ipsat_bp_3.5_bq_1.0_w_75_q2_0", r"$B_{qc}=3.5\,\mathrm{GeV}^{-2}, B_q=1.0\,\mathrm{GeV}^{-2}$", Linestyle(1), 'blue', "", 1.0, "blue"],

]



#p1.plot(np.NaN, np.NaN, '-', color="white", label=r"IPsat")

style=-1
color=-1
for f in files:
    color = color + 1
    style = style + 1
    xdata=[]
    ydata=[]
    staterrs=[]
    lower=[]
    upper=[]
    tmp=[]
    fname = "proton/coherent/" +f[0]
    print fname
    try:
        readfile_xy(fname, xdata, ydata)
        scale_list(ydata, GEVSQRTONB*f[5])
    except Exception, e:
        print "Error with file " + fname + ": " + str(e)
        continue

    if ShowStatErrs:
        try:
            readfile_xy(fname, tmp, staterrs, ycol=4)
            scale_list(staterrs, GEVSQRTONB*f[5])
            readfile_xy(fname, tmp, upper, ycol=2)
            readfile_xy(fname, tmp, lower, ycol=3)
            scale_list(lower, GEVSQRTONB*f[5])
            scale_list(upper, GEVSQRTONB*f[5])
        except:
            print "Cant read staterrs for file " + fname
            staterrs=[]

    # band
    if ShowBand and len(upper)>0:
        p1.fill_between(xdata, lower, upper, alpha=0.2, edgecolor='#1B2ACC', facecolor='#089FFF')

    if ShowStatErrs and len(staterrs)>0:
        p1.plot(xdata, ydata, linestyle=f[2], color=f[6], label=f[1], linewidth=lw_coh)
        p1.fill_between(xdata, np.array(ydata)-np.array(staterrs), np.array(ydata)+np.array(staterrs), alpha=0.2, facecolor=f[3], edgecolor='none', hatch=f[4])

        #p1.errorbar(xdata, ydata, marker=datadashes[1], yerr=staterrs, label=f[1], linewidth=1.8, markersize=2, linestyle=f[2], color=Color(color))
        totxs = scipy.integrate.simps(np.array(ydata),np.array(xdata))
        print f[1] + " totxs " + str(totxs)
    else:
        #plot
        p1.plot(xdata, ydata, linestyle=f[2], color=f[6], label=f[1], linewidth=lw_coh)

        totxs = scipy.integrate.simps(np.array(ydata),np.array(xdata))
        print f[1] + " totxs " + str(totxs)

# incoh
color = 0
for f in files:
    
    style = style + 1
    xdata=[]
    ydata=[]
    tmp=[]
    staterr=[]
    fname = "proton/incoherent/" + f[0]
    print fname
    try:
        readfile_xy(fname, xdata, ydata)
        scale_list(ydata, GEVSQRTONB*f[5])
        readfile_xy(fname, tmp, staterr, ycol=2)
        scale_list(staterr, GEVSQRTONB*f[5])
    
    except Exception, e:
        print "Error " + str(e) + " with file " + fname
        continue
    
    #p1.plot(xdata, ydata, linestyle=f[2], color=Color(color), label=r"", linewidth=1.0)
    p1.plot(xdata, ydata, label="", linestyle=f[2], color=f[6], linewidth=lw_incoh)
    p1.fill_between(xdata, np.array(ydata)-np.array(staterr), np.array(ydata)+np.array(staterr), alpha=0.2, edgecolor='none', facecolor=f[3], hatch=f[4])
    color = color + 1

#p1.plot(np.NaN, np.NaN, '-', color='white', label=r"$\mathrm{}$")
#p1.plot(np.NaN, np.NaN, '-', color='white', label=r"$\mathrm{}$")


# h1 data
expx=[]
expy=[]
experr=[]
tmp=[]
readfile_xy("proton/coherent/exp/h1_q2_0_w_100", expx, expy)
readfile_xy("proton/coherent/exp/h1_q2_0_w_100", tmp, experr, ycol=2)
#p1.errorbar(expx, expy, yerr=experr, marker=datadashes[2], linestyle='None', linewidth=1, markersize=markersize_coh, label=r"Coherent H1", color=Color(1), markeredgecolor=Color(1))


# h1 data
expx=[]
expy=[]
experr=[]
tmp=[]
readfile_xy("proton/coherent/exp/h1_q2_3.2_w_100", expx, expy)
readfile_xy("proton/coherent/exp/h1_q2_3.2_w_100", tmp, experr, ycol=2)
##p1.errorbar(expx, expy, yerr=experr, marker=datadashes[2], linestyle='None', linewidth=1, label=r"H1 $Q^2=3.2\,\mathrm{GeV}^2$")


# h1 data
expx=[]
expy=[]
experr=[]
tmp=[]
readfile_xy("proton/coherent/exp/h1_q2_22.4_w_100", expx, expy)
readfile_xy("proton/coherent/exp/h1_q2_22.4_w_100", tmp, experr, ycol=2)
##p1.errorbar(expx, expy, yerr=experr, marker=datadashes[2], linestyle='None', linewidth=1, label=r"H1 $Q^2=22.4\mathrm{GeV}^2$")

expx=[]
expy=[]
experr=[]
tmp=[]
readfile_xy("proton/coherent/exp/zeus_q2_0_w_100", expx, expy)
readfile_xy("proton/coherent/exp/zeus_q2_0_w_100", tmp, experr, ycol=2)
#p1.errorbar(expx, expy, yerr=experr, marker=datadashes[1], linestyle='None', linewidth=0.7, markersize=markersize_coh, color=Color(2), markeredgecolor=Color(2),label=r"Coherent ZEUS")

expx=[]
expy=[]
experr=[]
tmp=[]
readfile_xy("proton/incoherent/exp/h1_jpsi_w_100", expx, expy)
readfile_xy("proton/incoherent/exp/h1_jpsi_w_100", tmp, experr, ycol=2)
#p1.errorbar(expx, expy, yerr=experr, marker=datadashes[2], linestyle='None', linewidth=0.7, markersize=markersize_incoh, fillstyle='none', color=Color(1), markeredgecolor=Color(1), label=r"Total H1")


expx=[]
expy=[]
pluserr=[]
minuserr=[]
readfile_xyerrors("proton/incoherent/exp/zeus_jpsi_w_100", expx, expy, pluserr, minuserr)
# units
#scale_list(expy, 1000)
#scale_list(pluserr, 1000)
#scale_list(minuserr, 1000)
#p1.errorbar(expx, expy, yerr=[minuserr,pluserr], marker=datadashes[1], linestyle='None', linewidth=0.7, markersize=markersize_incoh, fillstyle='none', color=Color(2), markeredgecolor=Color(2), label=r"Incoherent ZEUS")

#old zeus
expx=[]
expy=[]
pluserr=[]
minuserr=[]
readfile_xyerrors("proton/incoherent/exp/zeus_jpsi_w_94", expx, expy, pluserr, minuserr)
# units
scale_list(expy, 1000)
scale_list(pluserr, 1000)
scale_list(minuserr, 1000)
##p1.errorbar(expx, expy, yerr=[minuserr,pluserr], marker=datadashes[4], linestyle='None', linewidth=0.7, markersize=4.4, label=r"incoh old ZEUS $Q^2=0, W=95$")

### W = 75 data
expx=[]
expy=[]
experr=[]
tmp=[]
readfile_xy("proton/coherent/exp/h1_jpsi_w_75", expx, expy)
readfile_xy("proton/coherent/exp/h1_jpsi_w_75", tmp, experr, ycol=2)
p1.errorbar(expx, expy, yerr=experr, marker=datadashes[2], linestyle='None', linewidth=0.7, markersize=markersize_coh, color=Color(2), markeredgecolor=Color(2), label=r"Coherent H1")
expx=[]
expy=[]
experr=[]
tmp=[]
readfile_xy("proton/incoherent/exp/h1_jpsi_w_75", expx, expy)
readfile_xy("proton/incoherent/exp/h1_jpsi_w_75", tmp, experr, ycol=2)
p1.errorbar(expx, expy, yerr=experr, marker=datadashes[2], linestyle='None', linewidth=0.7, markersize=markersize_incoh, fillstyle='none', color=Color(1), markeredgecolor=Color(1), label=r"Incoherent H1")



yscale("log",nonposy='clip')
#xscale("log")
axis([minx,maxx,miny,maxy])

legfont = textsize-8.5
leg=legend(prop=dict(size=legfont+3),labelspacing=0.001,ncol=2,numpoints=1, loc=1)
leg.draw_frame(False)

#plt.gcf().subplots_adjust(top=0.9)

plt.gcf().subplots_adjust(bottom=0.13)

file = "./proton_spectra.pdf"
print file
pp = PdfPages(file)
savefig(pp, format='pdf')
pp.close()


