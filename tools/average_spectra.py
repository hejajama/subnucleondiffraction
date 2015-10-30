# utf-8-

# Average the calculated spectra

import os
import sys

sys.path.append("/Users/heikki/lib/")
sys.path.append("/nashome2/hejajama/lib/")
from matplotlibhelper import *

dir=sys.argv[1]

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

# average
nconf = len(tmpydatas)

print "# Number of files: " + str(len(tmpydatas))

for i in range(len(tmpydatas[0])):
    t = tmpxdatas[0][i]
    try:    # if some file is not calcualted to high enough t, this we just skip that t value
        sum=0
        for j in range(len(tmpydatas)):
            sum+=tmpydatas[j][i]
        avg = sum/len(tmpydatas)

        tdata.append(t)
        ydata.append(avg)
    except:
        print "# Skipping t index " + str(i) + " t=" + str(t)

for x,y in zip(tdata, ydata):
    print x,y
