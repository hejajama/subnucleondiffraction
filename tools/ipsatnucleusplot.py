import sys
sys.path.append("/nashome2/hejajama/lib/")
sys.path.append("/home/hejajama/lib/")
sys.path.append("/home/heikki/lib")
sys.path.append("/Users/heikki/lib")
import math
from math import pow
from matplotlibhelper import *

minx=-60
maxx=60
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
p1=fig.add_subplot(111)
xlabel(r"$x$ $[\mathrm{GeV}^{-1}]$")
ylabel(r"$x$ $[\mathrm{GeV}^{-1}]$")

col=-1
for (x,y,r) in zip(xcoords,ycoords,radius):
    col=col+1

    circle = Circle((x,y),r)
    p1.add_artist(circle)


axis([minx,maxx,miny,maxy])

pp=PdfPages("nuke.pdf")
savefig(pp, format='pdf')
pp.close()
