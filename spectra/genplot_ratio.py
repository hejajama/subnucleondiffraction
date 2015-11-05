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

slides=True

rc('text',usetex=True)
rc('text.latex',  preamble='\usepackage{amsmath},\usepackage{amssymb},\usepackage{mathtools}')
textsize=18
if slides:
    textsize=26
rc("xtick", labelsize=textsize)
rc("ytick", labelsize=textsize)

minx=0.0
maxx=0.4
miny=0
maxy=3


PI=3.141592

GEVSQRTONB = 1.0e7/(5.068*5.068)

# filename normalization

# nominator fname nominator norm  denominator name denom norm  lable style
data=[
      ["mc_1e8/ipnonsat_1e8_nconf_100", GEVSQRTONB, "old_code/ipnonsat_incoh_step_001", 1, r"IPnonsat $10^8$", Linestyle(0)],
      ["mc_1e7/ipnonsat_1e7_nconf_100", GEVSQRTONB, "old_code/ipnonsat_incoh_step_0005", 1, r"IPnonsat $10^7$ nc 100", Linestyle(1)],
      ["total/ipnonsat_miser_1e7_nconf_200", GEVSQRTONB, "old_code/ipnonsat_incoh_step_0005", 1, r"IPnonsat $10^7$ nc 200", Linestyle(1)],
      #["coherent/coherent_ipnonsat_vegas_1e6", GEVSQRTONB, "old_code/ipnonsat_coh_step_0005", 1, r"Coherent ipnonsat", Linestyle(3)]
      ]

fig = figure()
p1=fig.add_subplot(111)
#xlabel(r"$r$ $[\mathrm{GeV}^{-1}]$ ", fontsize=textsize+2)
xlabel(r"$t$", fontsize=textsize+2)
#xlabel(r"$rQ_s$")
ylabel(r"$\mathrm{ratio}$ ", fontsize=textsize+2)
for d in data:
    nominator = [d[0], d[1]]
    denominator=[d[2],d[3]]
    x1=[]
    y1=[]
    x2=[]
    y2=[]

    readfile_xy(nominator[0], x1, y1)
    scale_list(y1, nominator[1])
        
    readfile_xy(denominator[0], x2, y2)
    scale_list(y1, denominator[1])

    ratio=[]
    for nom,denom,x,xx in zip(y1,y2,x1,x2):
        if x != xx:
            print "t values do not match: " + str(x) + " and " + str(xx)
        ratio.append(nom/denom)
        print ratio[-1]

    p1.plot(x1, ratio, label=d[4], linestyle=d[5], linewidth=2)
    

    
fig.suptitle(r"$Q^2=0, x=0.001$"   )
    
#yscale("log")
#xscale("log")
axis([minx,maxx,miny,maxy])

legfont = textsize-10
leg=legend(prop=dict(size=legfont),labelspacing=0.001,ncol=2,numpoints=1, loc=1)
leg.draw_frame(False)

            
        
file = "./ratio.pdf"
print file
pp = PdfPages(file)
savefig(pp, format='pdf')
pp.close()


