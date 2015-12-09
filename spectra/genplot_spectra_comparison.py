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
maxx=0.4
miny=1e0
maxy=1e7


PI=3.141592

GEVSQRTONB = 1.0e7/(5.068*5.068)

mc="1e6"
nconf="1000"

# filename title normalization style
files = [
         [ "old_code/ipsat_incoh", "fIPsat incoh", 1, Linestyle(0)],
         ["old_code/ipnonsat_incoh_step_0005", "Total ipnonsat", 1, Linestyle(0) ],
         #[ "old_code/ipsatnofactor_coh", "coherent IPsat", 1, Linestyle(0)],
         ["old_code/ipnonsat_coh_step_0005", "coherent IPnonsat", 1, Linestyle(0)],
         #["coherent/ipnonsat_vegas_" + mc + "_fluct_0_nconf_" + nconf, "coherent ipnonsat vegas " + mc + " nc " + nconf, GEVSQRTONB, Linestyle(2)],
         #["coherent/ipnonsat_vegas_" + mc + "_fluct_0.2_nconf_" + nconf, "coherent ipnonsat vegas " + mc + " nc " + nconf + "fl 0.2", GEVSQRTONB, Linestyle(3)],
         ["total/ipnonsat_vegas_" + mc + "_fluct_0_nconf_" + nconf, " ipnonsat vegas " + mc + " nc " + nconf, GEVSQRTONB, Linestyle(2)],
         ["total/ipnonsat_vegas_" + mc + "_fluct_0.2_nconf_" + nconf, "ipnonsat vegas " + mc + " nc " + nconf + " fl 0.2", GEVSQRTONB, Linestyle(3)],
         ["coherent/ipnonsat_miser_5e6_fluct_0_nconf_500", "miser 5e6 nc 500", GEVSQRTONB, Linestyle(1)],
         ["coherent/ipnonsat_miser_5e6_fluct_0_nconf_1200", "miser 5e6 nc 1200", GEVSQRTONB, Linestyle(2)],
         ["coherent/ipnonsat_miser_5e6_fluct_0_nconf_1500", "miser 5e6 nc 1500", GEVSQRTONB, Linestyle(3)],
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
    p1.plot(xdata, ydata, label=f[1], linestyle=f[3], linewidth=1, color=Color(col))
    

    
    fig.suptitle(r"$Q^2=0, x=0.001$"   )
    
    yscale("log")
    #xscale("log")
    axis([minx,maxx,miny,maxy])

    legfont = textsize-10
    leg=legend(prop=dict(size=legfont),labelspacing=0.001,ncol=2,numpoints=1, loc=1)
    leg.draw_frame(False)

            
        
file = "./spectra.pdf"
print file
pp = PdfPages(file)
savefig(pp, format='pdf')
pp.close()


