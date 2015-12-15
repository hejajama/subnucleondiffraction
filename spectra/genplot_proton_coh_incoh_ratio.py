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
miny=0.1
maxy=14


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
ylabel(r"coherent/incoherent", fontsize=textsize+2)



fig.suptitle(r"$Q^2=0, x=0.001$, $W \sim 100 \mathrm{GeV}$"   )


files = [
         ["gaussian/ipsat_miser_1e7_bp_1.0_bq_2.75", r"$B_p=1.0, B_q=2.75$", Linestyle(1),""],
         ["gaussian/ipsat_miser_1e7_bp_1.5_bq_2.25", r"$B_p=1.5, B_q=2.25$", Linestyle(2),""],
         ["gaussian/ipsat_miser_1e7_bp_2.0_bq_1.75", r"$B_p=2.0, B_q=1.75$", Linestyle(3),""],
         ]




# systematical
style=0
color=0
for f in files:
    style = style + 1
    xdata_coh=[]
    ydata_coh=[]
    fname = "proton/coherent/" + f[0]
    print fname
    try:
        readfile_xy(fname, xdata_coh, ydata_coh)
        scale_list(ydata_coh, GEVSQRTONB)
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
    for t,tot,coh in zip(xdata_coh,ydata_tot, ydata_coh):
        ydata_incoh.append(coh/(tot-coh))
        tdata.append(t)
    
    
    p1.plot(tdata, ydata_incoh, linestyle=f[2], marker=f[3], color=Color(color), label=f[1], linewidth=0.5, markersize=2.5)


expx=[]
expy=[]
pluserr=[]
minuserr=[]
readfile_xyerrors("proton//exp/coh_incoh_ratio_zeus", expx, expy, pluserr, minuserr)
p1.errorbar(expx, expy, yerr=[minuserr,pluserr], marker=datadashes[2], linestyle='None', linewidth=0.7, markersize=0.8, label=r"ZEUS $Q^2=0\mathrm{GeV}^2, <W>=94$ GeV")

yscale("log", nonposy='clip')
#xscale("log")
axis([minx,maxx,miny,maxy])

legfont = textsize-10
leg=legend(prop=dict(size=legfont),labelspacing=0.001,ncol=2,numpoints=1, loc=1)
leg.draw_frame(False)

file = "./proton_cohincoh_ratio.pdf"
print file
pp = PdfPages(file)
savefig(pp, format='pdf')
pp.close()


