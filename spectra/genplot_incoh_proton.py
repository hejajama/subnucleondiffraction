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
         [ "old_code/proton_ipsat", "IPsat", 1, Linestyle(0), Color(0)],
         
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

bpvals = ["5.0","6.0"]
bpvals=["4.5","7.0"]
bqvals = ["2.0","2.5"] #["2.0","3.0"]

style=0
color=0
for bp in bpvals:
    color=color+1
    #for
    #color = color + 1
    for bq in  bqvals:
        style = style + 1
        xdata=[]
        ydata=[]
        fname = "proton/coherent/ipsat_miser_1e7_bp_" + bp + "_bq_" + bq
        print fname
        try:
            readfile_xy(fname, xdata, ydata)
            scale_list(ydata, GEVSQRTONB)
        except:
            print "Error with file " + fname
        p1.plot(xdata, ydata, linestyle=Linestyle(color), marker=Datastyle(style), color=Color(color), label=r"$B_p=" + bp +", B_q=" + bq + r"$", linewidth=1,markersize=5)


# systematical
style=0
color=0
#params=[ ["5.0","2.5"],["6.0","2.5"],["7.0","2.5"],["5.0","3.0"], ["5.5","3.0"]]
#for p in params:
#        bp = p[0]
#        bq=p[1]
for bp in bpvals:
    color=color+1
    for bq in bqvals:
        style = style + 1
        xdata=[]
        ydata=[]
        fname = "proton/coherent/ipsat_miser_1e7_bp_" + bp + "_bq_" + bq
        print fname
        try:
            readfile_xy(fname, xdata, ydata)
            scale_list(ydata, GEVSQRTONB)
        except:
            print "Error with file " + fname
        
        # read total
        xdata_tot=[]
        ydata_tot=[]
        fname = "proton/total/ipsat_miser_1e7_bp_" + bp + "_bq_" + bq
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

        p1.plot(tdata, ydata_incoh, linestyle=Linestyle(3), marker=Datastyle(style), color=Color(color), label=r"", linewidth=0.5, markersize=2.5)

# h1 data
expx=[]
expy=[]
experr=[]
tmp=[]
readfile_xy("proton/coherent/exp/h1_qsqr_005", expx, expy)
readfile_xy("proton/coherent/exp/h1_qsqr_005", tmp, experr, ycol=2)
p1.errorbar(expx, expy, yerr=experr, marker=datadashes[0], linestyle='None', linewidth=1, label=r"H1 $Q^2=0.05\mathrm{GeV}^2, 40 < W  < 160 \mathrm{GeV}$")

expx=[]
expy=[]
experr=[]
tmp=[]
readfile_xy("proton/coherent/exp/zeus_qsqr_0", expx, expy)
readfile_xy("proton/coherent/exp/zeus_qsqr_0", tmp, experr, ycol=2)
p1.errorbar(expx, expy, yerr=experr, marker=datadashes[1], linestyle='None', linewidth=0.7, markersize=0.8, label=r"ZEUS $Q^2=0\mathrm{GeV}^2, 90 < W < 110 \mathrm{GeV}$")

expx=[]
expy=[]
pluserr=[]
minuserr=[]
readfile_xyerrors("proton/incoherent/exp/zeus", expx, expy, pluserr, minuserr)
# units
scale_list(expy, 1000)
scale_list(pluserr, 1000)
scale_list(minuserr, 1000)
#p1.errorbar(expx, expy, yerr=[minuserr,pluserr], marker=datadashes[2], linestyle='None', linewidth=0.7, markersize=0.8, label=r"Incoh ZEUS $Q^2=0\mathrm{GeV}^2$")

expx=[]
expy=[]
experr=[]
tmp=[]
readfile_xy("proton/incoherent/exp/h1_thesis", expx, expy)
readfile_xy("proton/incoherent/exp/h1_thesis", tmp, experr, ycol=2)
p1.errorbar(expx, expy, yerr=experr, marker=datadashes[2], linestyle='None', linewidth=0.7, markersize=2.1, label=r"incoh H1 $Q^2 \le 2.5\mathrm{GeV}^2, 40 < W < 110 \mathrm{GeV}$")

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


