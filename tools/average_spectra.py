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
corrections_file=""  # if set, read corrections from separated file
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
    elif sys.argv[i]=="-corrections":
        corrections_file = sys.argv[i+1]
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

include_imag = True
if os.path.exists(imagdir):
    include_imag = False
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
        tmplist=[]
        corrections=[]
        readfile_xy(fname_real, tmpxdata, tmpydata)
        if include_imag:
            readfile_xy( fname_imag, tmpxdata_imag, tmpydata_imag)
    except IOError:
        print "#File not found: " + fname_real + " or " + fname_imag
        continue



    # If we calculate incoherent scattering, we average the squared amplitude
    if coherent == False:
        for i in range(len(tmpydata)):
            tmpydata[i] = tmpydata[i]*tmpydata[i]
        
        for i in range(len(tmpydata_imag)):
            tmpydata_imag[i] = tmpydata_imag[i]*tmpydata_imag[i]

    tmpxdatas.append(tmpxdata)
    tmpydatas.append(tmpydata)
    fnames.append(fname_real)

    if imagdir != "":
        tmpxdatas_imag.append(tmpxdata_imag)
        tmpydatas_imag.append(tmpydata_imag)

# Read corrections if asked
if corrections_file != "":
    tmptdata = []
    corrections = []
    readfile_xy(corrections_file, tmptdata, corrections)


# average
nconf = len(tmpydatas)

print "# Dir: " + dir
if not include_imag:
    print "# Not including imaginary part, dir does not exist"
print "# Number of files: " + str(len(tmpydatas))

correction_maxt_error = False  # display this error only once
for i in range(len(tmpydatas[0])):
    t = tmpxdatas[0][i]
    correction = 1.0
    if corrections_file != "" and i >= len(corrections):
        # We were not able to calculate corrections at high t
        # as corrections depend very weakly on t, we could just use the highest t value
        correction = corrections[len(corrections)-1]
        if correction_maxt_error == False:
            print "# No correction calculated at t=" + str(t) + ", from now on use correction(t=" + str(tmpxdatas[0][i-1]) + "=" + str(correction) + " (note: correction(t=0)=" + str(corrections[0]) + ")"
            correction_maxt_error = True
    elif corrections_file != "":
        correction = corrections[i]
    
    fileind=0
    try:    # if some file is not calcualted to high enough t, this we just skip that t value
        sum=0
        sum_imag=0
        for j in range(len(tmpydatas)):
            fileind = j
            sum+=tmpydatas[j][i]
            if include_imag:
                sum_imag+=tmpydatas_imag[j][i]

        avg = sum/len(tmpydatas)
        avg_imag=0
        if include_imag:
            avg_imag = sum_imag/len(tmpydatas_imag)
        
        
        if coherent:    # Coherent scattering: now we have averated, then we squared
            avg = (avg*avg + avg_imag*avg_imag) / (16.0*pi)
        else:   # Incoherent
            avg = (avg + avg_imag) / (16.0*pi)

        tdata.append(t)
        ydata.append(avg*correction)
    except:
        print "# Skipping t index " + str(i) + " t=" + str(t) + " because of file " +fnames[fileind]

for x,y in zip(tdata, ydata):
    print x,y
