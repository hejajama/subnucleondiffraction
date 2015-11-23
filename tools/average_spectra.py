# utf-8-

# Average the calculated spectra

import os
import sys
import math



sys.path.append("/Users/heikki/lib/")
sys.path.append("/nashome2/hejajama/lib/")
from matplotlibhelper import *

dir = ""
imagdir = ""
coherent = False
maxnconfs = 99999999 # can limit number of configs

for i in range(len(sys.argv)):
    if sys.argv[i]=="-coherent":
        coherent= True
    elif sys.argv[i]=="-total":
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

tdata_imag=[]
ydata_imag=[]

first=1

tmpxdatas=[]
tmpydatas=[]
tmpxdatas_imag=[]
tmpydatas_imag=[]

imagdir = dir + "/imag/"
dir = dir + "/real/"

fnames=[]

if coherent == False:
    print "# Total diffractive cross section"
else:
    print "# Coherent cross section"

files = os.listdir(dir)
# Note: we assume that imaginary part is in the file with a same name but in dir imagdir
for f in files:
    tmpydata=[]
    tmpxdata=[]
    
    tmpxdata_imag=[]
    tmpydata_imag=[]
    
    fname_real = dir + f
    fname_imag = imagdir + f
    
    parse = fname_real.split("_")
    try:
        if int(parse[-1]) > maxnconfs:
            continue
    except ValueError:
            print "WTF file " + fname
            continue

    try:
        readfile_xy(fname_real, tmpxdata, tmpydata)
        readfile_xy( fname_imag, tmpxdata_imag, tmpydata_imag)
    except IOError:
        print "#File not found: " + fname_real + " or " + fname_imag
        continue



    # If we calculate incoherent scattering, we average the squared amplitude
    if coherent == False:
        for i in range(len(tmpydata)):
            tmpydata[i] = tmpydata[i]*tmpydata[i]
        for i in range(len(tmpydata_imag)):
            if imagdir != "":
                tmpydata_imag[i] = tmpydata_imag[i]*tmpydata_imag[i]

    tmpxdatas.append(tmpxdata)
    tmpydatas.append(tmpydata)
    fnames.append(fname_real)

    if imagdir != "":
        tmpxdatas_imag.append(tmpxdata_imag)
        tmpydatas_imag.append(tmpydata_imag)


# average
nconf = len(tmpydatas)

print "# Dir: " + dir
print "# Number of files: " + str(len(tmpydatas))

for i in range(len(tmpydatas[0])):
    t = tmpxdatas[0][i]
    fileind=0
    try:    # if some file is not calcualted to high enough t, this we just skip that t value
        sum=0
        sum_imag=0
        for j in range(len(tmpydatas)):
            fileind = j
            sum+=tmpydatas[j][i]
            if imagdir != "":
                sum_imag+=tmpydatas_imag[j][i]

        avg = sum/len(tmpydatas)
        avg_imag=0
        if imagdir != "":
            avg_imag = sum_imag/len(tmpydatas_imag)
        
        if coherent:    # Coherent scattering: now we have averated, then we squared
            avg = (avg*avg + avg_imag*avg_imag) / (16.0*pi)
        else:   # Incoherent
            avg = (avg + avg_imag) / (16.0*pi)

        tdata.append(t)
        ydata.append(avg)
    except:
        print "# Skipping t index " + str(i) + " t=" + str(t) + " because of file " +fnames[fileind]

for x,y in zip(tdata, ydata):
    print x,y
