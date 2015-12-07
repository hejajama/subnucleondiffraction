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
import numpy


def ComparisonFunction(x):
    #a = 0.55
    #return x*x*numpy.exp(-a*x)/( 2/math.pow(a,3.0) )
    #width = 0.337556/100
    #return 2.0/(math.sqrt(width*2*math.pi))*numpy.exp(-x*x/(2.0*width))
    #return 1.0/(2*math.exp*(2)/math.pow(a,3.0)) * x*x*numpy.exp(-a*x)
    return 0.5*x*numpy.exp(-x*x/4)

if len(sys.argv) != 3:
    print "Syntax: fname number_of_bins"
    sys.exit(0)

fname=sys.argv[1]
nbins=int(sys.argv[2])

f = open(fname)
lines = f.readlines()
data=[]
for n in lines:
    try:
        data.append(float(n))
    except:
        continue

n,bins,patches = hist(data, nbins, normed=1)

xlabel("Value")
ylabel("Frequency")
plt.grid(True)


# Compare with a test function
t = np.arange( min(data), max(data), (max(data)-min(data))/(nbins*10))
plot(t, ComparisonFunction(t), label='Distribution')

plt.show()