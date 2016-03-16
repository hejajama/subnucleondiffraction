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
textsize=19
if slides:
    textsize=26
rc("xtick", labelsize=textsize)
rc("ytick", labelsize=textsize)

def sqr(x):
    return x*x

minx=0.0
maxx=2.0
miny=0.01
maxy=40


PI=3.141592

GEVSQRTONB = 1.0e7/(5.068*5.068)

# filename normalization

# fname title style marker facecolor hatch, linecolor
data=[
      #["final/ipsat2012_bp_1.0_bq_3.0_w_100_q2_0", r"IPsat $B_{qc}=1.0\,\mathrm{GeV}^{-2}, B_q=3.0\,\mathrm{GeV}^{-2}$", Linestyle(3),"grey", "", "black"],
      #["final/ipsat2012_bp_3.5_bq_1.0_w_100_q2_0", r"IPsat, $B_{qc}=3.5\,\mathrm{GeV}^{-2}, B_q=1.0\,\mathrm{GeV}^{-2}$", Linestyle(2), 'red', "", "red"], # use linestyle 2
      #["final/ipsat2012_bp_3.5_bq_1.0_w_100_q2_0_gridfluct_04", r"IPsat$B_p=3.5\,\mathrm{GeV}^{-2}, B_q=1.0\,\mathrm{GeV}^{-2}, a=0.4\,\mathrm{fm}$", Linestyle(1),"grey", "", "black"],
      #["final/ipsat2012_bp_3.5_bq_1.0_w_100_q2_0_gridfluct_01", r"$B_p=3.5\,\mathrm{GeV}^{-2}, B_q=1.0\,\mathrm{GeV}^{-2}, a=0.1\,\mathrm{fm}$", Linestyle(3),""],
      #["final/ipsat2012_bp_3.5_bq_1.0_w_100_q2_0_quarkfluct_0.5", r"IPsat, $B_{qc}=3.5\,\mathrm{GeV}^{-2}, B_q=1.0\,\mathrm{GeV}^{-2}, \sigma=0.5$", Linestyle(1),"green", "", "green"],
      #["test/ipglasma_bp_2.0_bq_0.3_norm07_nc_96", r"IP-Glasma, $B_p=2.0\,\mathrm{GeV}^{-2}, B_q=0.3\,\mathrm{GeV}^{-2}$", Linestyle(3), 1.4],
      #["final/ipglasma_bp_2.0_bq_0.3_m04_n07_ncf_416_w_100_q2_0", r"IP-Glasma, $B_{qc}=2.0\,\mathrm{GeV}^{-2}, B_q=0.3\,\mathrm{GeV}^{-2}$", Linestyle(0),"grey", "", Color(0)],
      # ["paper_2/ipglasma_bp_2.0_bq_0.3_m04_n075_qsfluct_ncf_96", r"IP-Glasma, $B_{qc}=2.0\,\mathrm{GeV}^{-2}, B_q=0.3\,\mathrm{GeV}^{-2}$ fluct", Linestyle(2),Color(2), "", Color(2)],
      ["paper_2/ipglasma_bp_2.0_bq_0.5_m08_n075_qsfluct_ncf_96", r"IP-Glasma, $B_{qc}=2.0\,\mathrm{GeV}^{-2}, B_q=0.5\,\mathrm{GeV}^{-2} m=0.8$ ", Linestyle(2),Color(2), "", Color(2)],
      ["paper_2/ipglasma_bp_2.0_bq_0.5_m04_n075_qsfluct_ncf_384", r"IP-Glasma, $B_{qc}=2.0\,\mathrm{GeV}^{-2}, B_q=0.5\,\mathrm{GeV}^{-2} m=0.4$ ", Linestyle(3),Color(3), "", Color(3)],
      ["paper_2/ipglasma_bp_2.0_bq_0.5_m02_n075_qsfluct_ncf_96", r"IP-Glasma, $B_{qc}=2.0\,\mathrm{GeV}^{-2}, B_q=0.5\,\mathrm{GeV}^{-2}, m=0.2$", Linestyle(1),Color(1), "", Color(1)]
]

fig = figure()
p1=fig.add_subplot(111)
#xlabel(r"$r$ $[\mathrm{GeV}^{-1}]$ ", fontsize=textsize+2)
xlabel(r"$|t|$ $[\mathrm{GeV}^2]$", fontsize=textsize+2)
ylabel(r"$\mathrm{coherent/incoherent}$ ", fontsize=textsize+2)
cohdir = "proton/coherent/"
incohdir = "proton/incoherent/"
color=-1
for d in data:
    color=color+1
    x1=[]
    y1=[]
    x2=[]
    y2=[]
    
    tmp=[]
    err_1=[]
    err_2=[]

    readfile_xy(cohdir + d[0], x1, y1)
    readfile_xy(cohdir + d[0], tmp, err_1, ycol=4)
    readfile_xy(incohdir + d[0], x2, y2)
    readfile_xy(incohdir + d[0], tmp, err_2, ycol=2)

    ratio=[]
    err=[]
    
    for nom,denom,x,xx,err1,err2 in zip(y1,y2,x1,x2,err_1,err_2):
        if x != xx:
            print "t values do not match: " + str(x) + " and " + str(xx)
        ratio.append(nom/denom)

        #uncertainty
        err.append( sqrt( sqr(1.0/denom * err1) + sqr(nom/sqr(denom)*err2) ) )
        #print ratio[-1]

    p1.plot(x1, ratio, label=d[1], linestyle=d[2], linewidth=2, color=d[5])
    p1.fill_between(x1, np.array(ratio)-np.array(err), np.array(ratio)+np.array(err), alpha=0.2, edgecolor="none", facecolor=d[3], hatch=d[4])
    

plt.gcf().subplots_adjust(bottom=0.12)
plt.gcf().subplots_adjust(left=0.12)

#fig.suptitle(r"$Q^2=0, W=100\,\mathrm{GeV}$"   )
    
yscale("log", nonposy='clip')
#xscale("log")
axis([minx,maxx,miny,maxy])



# data
expx=[]
expy=[]
tmp=[]
#err_p=[]
#err_m = []
staterr=[]
syst_p = []
syst_m = []
readfile_xy("proton/ratio/zeus_coh_incoh_w_94_q2_0", expx, expy)
readfile_xy("proton/ratio/zeus_coh_incoh_w_94_q2_0", tmp, staterr, ycol=2)
readfile_xy("proton/ratio/zeus_coh_incoh_w_94_q2_0", tmp, syst_p, ycol=3)
readfile_xy("proton/ratio/zeus_coh_incoh_w_94_q2_0", tmp, syst_m, ycol=4)
#readfile_xy("proton/ratio/zeus_coh_incoh_w_94_q2_0_tote", tmp, err_p, ycol=2)
#readfile_xy("proton/ratio/zeus_coh_incoh_w_94_q2_0_tote", tmp, err_m, ycol=3)
#p1.errorbar(expx, expy, yerr=[err_m, err_p], marker=datadashes[1], linestyle='None', linewidth=0.7, markersize=4.4, label=r"ZEUS", color=Color(0))


#statistical error as a box
width = 0.05
ex = np.ones(np.array(expx).shape)*width
ey=[]  # total yerr
for p,m in zip(syst_p, syst_m):
    ey.append(p+m)
bottoms=[] # shift the center of the box
for i in range(len(ey)):
    bottoms.append( expy[i] - syst_m[i])
    print ey[i],bottoms[i]
p1.bar(expx, ey, bottom=bottoms, width=width, align='center', alpha=0.4, linewidth=0)
p1.errorbar(expx, expy, yerr=staterr, marker=datadashes[1], linestyle='None', linewidth=0.7, markersize=5.0, label=r"ZEUS", color=Color(0), zorder=13)

legfont = textsize-8.5
leg=legend(prop=dict(size=legfont+3),labelspacing=0.001,ncol=1,numpoints=1, loc=1)
leg.draw_frame(False)

file = "./ratio.pdf"
print file
pp = PdfPages(file)
savefig(pp, format='pdf')
pp.close()


