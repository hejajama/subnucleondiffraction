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
maxnconfs = 99999999 # can limit number of configs

for i in range(len(sys.argv)):
    if sys.argv[i]=="-coherent":
        coherent= True
    elif sys.argv[i]=="-incoherent":
        coherent = False
    elif sys.argv[i]=="-dir":
        dir = sys.argv[i+1]
    elif sys.argv[i]=="-maxconf":
        maxnconfs = int(sys.argv[i+1])
    elif sys.argv[i][0]=="-":
        print "Unknown argument " + sys.argv[i]
        sys.exit(1)


tdata=[]
ydata=[]

first=1

tmpxdatas=[]
tmpydatas=[]

if coherent == False:
    print "# Total diffractive cross section"
else:
    print "# Coherent cross section"

files = os.listdir(dir)
for f in files:
    tmpydata=[]
    tmpxdata=[]
    
    fname = dir + f
    
    parse = fname.split("_")
    if int(parse[-1]) > maxnconfs:
        continue

    readfile_xy(fname, tmpxdata, tmpydata)

    # If we calculate incoherent scattering, we average the squared amplitude
    if coherent == False:
        for i in range(len(tmpydata)):
            tmpydata[i] = tmpydata[i]*tmpydata[i]

    tmpxdatas.append(tmpxdata)
    tmpydatas.append(tmpydata)


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
        
        if coherent:    # Coherent scattering: now we have averated, then we squared
            avg = avg*avg / (16.0*pi)
        else:   # Incoherent
            avg = avg / (16.0*pi)

        tdata.append(t)
        ydata.append(avg)
    except:
        print "# Skipping t index " + str(i) + " t=" + str(t)

for x,y in zip(tdata, ydata):
    print x,y
