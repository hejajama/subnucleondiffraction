#!/usr/bin/python
#-*- coding: UTF-8 -*-
import os
import sys
sys.path.append("/nashome2/hejajama/lib/")
sys.path.append("/home/hejajama/lib/")
sys.path.append("/Users/heikkimantysaari/lib/")

import math
from matplotlibhelper import *
#from virtualphoton import *
import scipy.stats
#from fit import *
import pylab
import sys
import pdb
rc('text',usetex=True)
rc('text.latex',  preamble='\usepackage{amsmath},\usepackage{amssymb},\usepackage{mathtools}')

PLOT_PDF = True

textsize=15
rc("xtick", labelsize=textsize)
rc("ytick", labelsize=textsize)




minx=-0.002
maxx=0.08
miny=0
maxy=9

ms=4

fig = figure()
p1=fig.add_subplot(111)

#xlabel(r"$r$ $[\mathrm{GeV}^{-1}]$ ", fontsize=textsize+2)
xlabel(r"$v$ $[\mathrm{fm}]$", fontsize=textsize+2)
#xlabel(r"$rQ_s$")
ylabel(r"$\chi^2/N$", fontsize=textsize+2)


xdata=[]
ydata=[]
readfile_xy("chisqr.txt", xdata, ydata, ycol=2)

p1.plot(xdata, ydata, linestyle=Linestyle(3), linewidth=1, marker=Datastyle(2), markersize=ms, color=Color(1), label="Fit quality, $Q^2>50\,\mathrm{GeV}^2$")

xdata=[]
ydata=[]
readfile_xy("chisqr.txt", xdata, ydata)


p1.plot(xdata, ydata, linestyle=Linestyle(3), linewidth=1, marker=Datastyle(1), markersize=ms, color=Color(0), label="Fit quality, all $Q^2$")




axis([minx,maxx,miny,maxy])
leg=legend(prop=dict(size=textsize-2),labelspacing=0.001,ncol=1,numpoints=1, loc=1)
leg.draw_frame(False)


#fig.suptitle(r"Fit $A*Exp(-B*|t|)$ to calculation, $N=" + str(fitpoints) + r", n=2$", fontsize=textsize+2)

#plt.gcf().subplots_adjust(top=0.9)

#plt.gcf().subplots_adjust(bottom=0.13)

file = "./chisqr_uv.pdf"
print file
pp = PdfPages(file)
savefig(pp, format='pdf')
pp.close()

