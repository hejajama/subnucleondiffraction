import sys
sys.path.append("/nashome2/hejajama/lib/")
sys.path.append("/home/hejajama/lib/")
sys.path.append("/home/heikki/lib")
sys.path.append("/Users/heikki/lib")
import math
from math import pow
from matplotlibhelper import *

FMGEV = 5.08
R_p = 0.557


minx=-2*R_p
maxx=2*R_p
miny=minx
maxy=maxx


fname = sys.argv[1]

xcoords=[]
ycoords=[]
radius=[]
tmp=[]

readfile_xy(fname, xcoords, ycoords)
readfile_xy(fname, tmp, radius, ycol=2)


fig=figure()

p1=fig.add_subplot(111, aspect='equal')

xlabel(r"$x$ $[\mathrm{fm}]$",fontsize=15)
ylabel(r"$x$ $[\mathrm{fm}]$",fontsize=15)

col=-1
for (x,y,r) in zip(xcoords,ycoords,radius):
    col=col+1
    x=x / FMGEV
    y = y / FMGEV
    r = r / FMGEV

    circle = Circle((x,y),r)
    p1.add_artist(circle)

# draw proton
circle = Circle((0,0),R_p, alpha=0.1, color="pink")
p1.add_artist(circle)

axis([minx,maxx,miny,maxy])

pp=PdfPages("proton.pdf")
savefig(pp, format='pdf')
pp.close()
