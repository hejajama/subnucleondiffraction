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

#import scipy.integrate

slides=False

rc('text',usetex=True)
rc('text.latex',  preamble='\usepackage{amsmath},\usepackage{amssymb},\usepackage{mathtools}')
textsize=18
if slides:
    textsize=26
rc("xtick", labelsize=textsize)
rc("ytick", labelsize=textsize)

minx=0.0
maxx=2
miny=0.01
maxy=1e3


PI=3.141592

GEVSQRTONB = 1.0e7/(5.068*5.068)

mc="1e6"
nconf="1000"

# filename title normalization style
files = [
         [ "old_code/ipsat_nofactor_proton_w_100", "IPsat", 1, Linestyle(0), Color(0)],
         
         ]


fig = figure()
p1=fig.add_subplot(111)
#xlabel(r"$r$ $[\mathrm{GeV}^{-1}]$ ", fontsize=textsize+2)
xlabel(r"$t$", fontsize=textsize+2)
#xlabel(r"$rQ_s$")
ylabel(r"$\mathrm{d}\sigma/\mathrm{d}t$ ", fontsize=textsize+2)

col=-1
for f in files:
    col=col+1
    xdata=[]
    ydata=[]
    try:
        readfile_xy(f[0], xdata, ydata)
        scale_list(ydata, f[2])
    except:
        print "Error with file " + f[0]
        continue
    p1.plot(xdata, ydata, label=f[1], linestyle=f[3], linewidth=1, color=f[4])
    

    
    fig.suptitle(r"$Q^2=0, x=0.001$, $W \sim 100 \mathrm{GeV}$"   )

# fname title style marker
files = [
         #["gaussian/ipsat_miser_1e7_bp_1.0_bq_2.75", r"$B_p=1.0, B_q=2.75$", Linestyle(1),""],
         #["gaussian/ipsat_miser_1e7_bp_1.0_bq_3.0", r"$B_p=1.0, B_q=3.0$", Linestyle(2),""],
         #["gaussian/ipsat_miser_1e7_bp_1.0_bq_3.25", r"$B_p=1.0, B_q=3.25$", Linestyle(3),""],
         #["gaussian/ipsat_miser_1e7_bp_1.5_bq_2.25", r"$B_p=1.5, B_q=2.25$", Linestyle(1),""],
         #["gaussian/ipsat_miser_1e7_bp_1.5_bq_2.5", r"$B_p=1.5, B_q=2.5$", Linestyle(2),""],
         #["gaussian/ipsat_miser_1e7_bp_1.5_bq_2.75", r"$B_p=1.5, B_q=2.75$", Linestyle(3),""],
         ["gaussian/ipsat_miser_1e7_bp_2.0_bq_1.75", r"$B_p=2.0, B_q=1.75$", Linestyle(1),""],
         ["gaussian/ipsat_miser_1e7_bp_2.0_bq_2.0", r"$B_p=2.0, B_q=2.0$", Linestyle(2),""],
         ["gaussian/ipsat_miser_1e7_bp_2.0_bq_2.25", r"$B_p=2.0, B_q=2.25$", Linestyle(3),""],
         #["exponential/ipsat_miser_1e7_a_1.3_bq_2.0", r"$a=1.3, B_q=2.0$", Linestyle(1),""],
         #["exponential/ipsat_miser_1e7_a_1.5_bq_1.5", r"$a=1.5, B_q=1.5$", Linestyle(2),""],
         #["exponential/ipsat_miser_1e7_a_1.1_bq_2.0", r"$a=1.1, B_q=2.0$", Linestyle(3),""],
]

style=0
color=0
for f in files:
    color = color + 1
    style = style + 1
    xdata=[]
    ydata=[]
    fname = "proton/coherent/" +f[0]
    print fname
    try:
        readfile_xy(fname, xdata, ydata)
        scale_list(ydata, GEVSQRTONB)
    except:
        print "Error with file " + fname
    p1.plot(xdata, ydata, linestyle=f[2], marker=f[3], color=Color(color), label=f[1], linewidth=1,markersize=5)


# incoh
color = 0
for f in files:
    color = color + 1
    style = style + 1
    xdata=[]
    ydata=[]
    fname = "proton/coherent/" + f[0]
    print fname
    try:
        readfile_xy(fname, xdata, ydata)
        scale_list(ydata, GEVSQRTONB)
    except:
        print "Error with file " + fname
    
    # read total
    xdata_tot=[]
    ydata_tot=[]
    fname = "proton/total/" + f[0]
    try:
        readfile_xy(fname, xdata_tot, ydata_tot)
        scale_list(ydata_tot, GEVSQRTONB)
    except:
        print "Error with file " + fname

    ydata_incoh=[]
    tdata=[]
    for t,tot,coh in zip(xdata,ydata_tot, ydata):
        ydata_incoh.append(tot-coh)
        tdata.append(t)

    p1.plot(tdata, ydata_incoh, linestyle=f[2], marker=f[3], color=Color(color), label=r"", linewidth=0.5, markersize=2.5)

# h1 data
expx=[]
expy=[]
experr=[]
tmp=[]
readfile_xy("proton/coherent/exp/h1_q_0_w_100", expx, expy)
readfile_xy("proton/coherent/exp/h1_q_0_w_100", tmp, experr, ycol=2)
p1.errorbar(expx, expy, yerr=experr, marker=datadashes[0], linestyle='None', linewidth=1, label=r"H1 $Q^2=0.05\mathrm{GeV}^2, 40 < W  < 160 \mathrm{GeV}$")

expx=[]
expy=[]
experr=[]
tmp=[]
#readfile_xy("proton/coherent/exp/zeus_qsqr_0", expx, expy)
#readfile_xy("proton/coherent/exp/zeus_qsqr_0", tmp, experr, ycol=2)
#p1.errorbar(expx, expy, yerr=experr, marker=datadashes[1], linestyle='None', linewidth=0.7, markersize=0.8, label=r"ZEUS $Q^2=0\mathrm{GeV}^2, 90 < W < 110 \mathrm{GeV}$")

expx=[]
expy=[]
pluserr=[]
minuserr=[]
#readfile_xyerrors("proton/incoherent/exp/zeus", expx, expy, pluserr, minuserr)
# units
scale_list(expy, 1000)
scale_list(pluserr, 1000)
scale_list(minuserr, 1000)
#p1.errorbar(expx, expy, yerr=[minuserr,pluserr], marker=datadashes[2], linestyle='None', linewidth=0.7, markersize=0.8, label=r"Incoh ZEUS $Q^2=0\mathrm{GeV}^2$")

expx=[]
expy=[]
experr=[]
tmp=[]
#readfile_xy("proton/incoherent/exp/h1_thesis", expx, expy)
#readfile_xy("proton/incoherent/exp/h1_thesis", tmp, experr, ycol=2)
#p1.errorbar(expx, expy, yerr=experr, marker=datadashes[2], linestyle='None', linewidth=0.7, markersize=2.1, label=r"incoh H1 $Q^2 \le 2.5\mathrm{GeV}^2, 40 < W < 110 \mathrm{GeV}$")

yscale("log")
#xscale("log")
axis([minx,maxx,miny,maxy])

legfont = textsize-10
leg=legend(prop=dict(size=legfont),labelspacing=0.001,ncol=2,numpoints=1, loc=1)
leg.draw_frame(False)
        
file = "./proton_incoh_spectra.pdf"
print file
pp = PdfPages(file)
savefig(pp, format='pdf')
pp.close()


