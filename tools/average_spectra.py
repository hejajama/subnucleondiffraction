# utf-8-

# Average the calculated spectra

import os
import sys
import math



sys.path.append("/Users/heikki/lib/")
sys.path.append("/nashome2/hejajama/lib/")
from matplotlibhelper import *

dir = ""
coherent = False

if sys.argv[1] == "-coherent":
    coherent= True
    dir = sys.argv[2]
    print "# Calculating coherent cross section"
else:
    dir=sys.argv[1]

maxnconfs = 99999999 # can limit number of configs

if coherent == False:
    if len(sys.argv)>2:
        maxnconfs = int(sys.argv[2])

tdata=[]
ydata=[]

first=1

tmpxdatas=[]
tmpydatas=[]

for f in os.listdir(dir):
    tmpydata=[]
    tmpxdata=[]
    
    fname = dir + f

    readfile_xy(fname, tmpxdata, tmpydata)

    tmpxdatas.append(tmpxdata)
    tmpydatas.append(tmpydata)

    if len(tmpydatas) >= maxnconfs:
        break

# average
nconf = len(tmpydatas)

print "# Dir: " + dir
print "# Number of files: " + str(len(tmpydatas))

for i in range(len(tmpydatas[0])):
    t = tmpxdatas[0][i]
    try:    # if some file is not calcualted to high enough t, this we just skip that t value
        sum=0
        for j in range(len(tmpydatas)):
            sum+=tmpydatas[j][i]
        avg = sum/len(tmpydatas)
        
        if coherent:
            avg = avg*avg / (16.0*pi)

        tdata.append(t)
        ydata.append(avg)
    except:
        print "# Skipping t index " + str(i) + " t=" + str(t)

for x,y in zip(tdata, ydata):
    print x,y
