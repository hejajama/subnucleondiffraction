# utf-8-

# Average the calculated spectra

import os
import sys
import math

def sqr(x):
    return x*x


sys.path.append("/Users/heikki/lib/")
sys.path.append("/nashome2/hejajama/lib/")
from matplotlibhelper import *
import numpy as np

dir = ""
imagdir = ""
corrections_file=""  # if set, read corrections from separated file
maxnconfs = 99999999 # can limit number of configs

coherent = True

for i in range(len(sys.argv)):
    if sys.argv[i]=="-dir":
        dir = sys.argv[i+1]
    elif sys.argv[i]=="-incoherent":
        coherent=False
    elif sys.argv[i]=="-coherent":
        coherent=True
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

fnames=[]

imagdir = dir + "/imag/"
dir = dir + "/real/"

include_imag = True
if not os.path.exists(imagdir):
    include_imag = False
fnames=[]


# Read amplitudes, imag and real part, and transverse and longitudinal polarizations in each case
realparts = []
imagparts = []
tvals=[]
# structure: realparts[fileind][tind][transverse, longitudinal]


files = os.listdir(dir)
# Note: we assume that imaginary part is in the file with a same name but in dir imagdir
for f in files:
    longitudinal_real = []   # amplitudes for longitudinal scattering at each t
    transverse_real = []
    tvals_real = []
    
    longitudinal_imag = []   # amplitudes for longitudinal scattering at each t
    transverse_imag = []
    tvals_imag = []
    
    fname_real = dir + f
    fname_imag = imagdir + f
    
    
    
    
    parse = fname_real.split("_")
    try:
        if int(parse[-1]) > maxnconfs:
            continue
    except ValueError:
        print "# WTF file " + fname_real
        continue

    try:
        tmplist=[]
        readfile_xy(fname_real,tvals_real, transverse_real)
        readfile_xy(fname_real, tmplist, longitudinal_real, ycol=2)
        if include_imag:
            readfile_xy(fname_imag,tvals_imag, transverse_imag)
            readfile_xy(fname_imag, tmplist, longitudinal_imag, ycol=2)
    except IOError:
        print "#File not found: " + fname_real + " or " + fname_imag
        continue
    
    # assume that all files go to largest t
    if tvals == []:
        tvals = tvals_real
    
    # Save ampltiude values as a function of t
    tmp_realparts = []
    tmp_imagparts=[]
    for t,l in zip(transverse_real, longitudinal_real):
        tmp_realparts.append([t,l])
    if include_imag:
        for t,l in zip(transverse_imag, longitudinal_imag):
            tmp_imagparts.append([t,l])
    else:
        for i in range(len(transverse_real)):
            tmp_imagparts.append([0,0])

    realparts.append(tmp_realparts)
    imagparts.append(tmp_imagparts)

    fnames.append(fname_real)

print "# Read " + str(len(realparts)) + " amplitudes"


# Read corrections
corrections_t = []
corrections_l = []
tmp=[]
if corrections_file != "":
    readfile_xy(corrections_file, tmp, corrections_t)
    readfile_xy(corrections_file, tmp, corrections_l)

# Ok, calculate coheret xs
if coherent:
    # <A_T>^2 + <A_L>^2
    for t in range(len(realparts[0])):
        sum_real_t=0
        sum_real_l=0
        sum_imag_t=0
        sum_imag_l=0
        nconfs=0
        for conf in range(len(realparts)):
            # check if the calculations is done at high t
            if len(realparts[conf]) <= t or len(imagparts[conf])<=t:
                print "# Skip file " + fnames[conf] + " at t " + str(tvals[t])
                continue
            try:
                sum_real_t += realparts[conf][t][0]
                sum_real_l += realparts[conf][t][1]
                sum_imag_t += imagparts[conf][t][0]
                sum_imag_l += imagparts[conf][t][1]
            #print conf,t,realparts[conf][t][0], sum_real_t
            except IndexError:  # Some calculation was not done at high enough t, skip
                print "# WHY INDEX ERROR at " + str(conf) + " at t " + str(tvals[t])
                continue
            nconfs += 1

        transverse = sqr(sum_real_t/nconfs) + sqr(sum_imag_t/nconfs)
        longitudinal = sqr(sum_real_l/nconfs) + sqr(sum_imag_l/nconfs)

        # Get correction
        correction_t=1.0
        correction_l=1.0
        if corrections_file != "":
            if len(corrections_t) <= t: # Use highest t calculated
                correction_t = corrections_t[len(corrections_t)-1]
                correction_l = corrections_l[len(corrections_l)-1]
            else:
                correction_t = corrections_t[t]
                correction_l = corrections_l[t]

        print tvals[t], (correction_t*transverse+correction_l*longitudinal)/(16.0*pi)





# Incoherent xs
# Variance, <|A|^2> - |<A>|^2
# Variance is calculated using scipy and more stable alogirthm
if coherent == False:
    for t in range(len(realparts[0])):
        # create lists
        real_l=[]
        real_t=[]
        imag_l=[]
        imag_t=[]
        for rc,ic in zip(realparts,imagparts):
            if len(rc) <= t or len(ic)<= t:
                print "# Skip file at t=" + str(tvals[t])
                continue
            real_l.append(rc[t][1])
            real_t.append(rc[t][0])
            imag_l.append(ic[t][1])
            imag_t.append(ic[t][0])
        
        var_real_t = np.var(real_t)
        var_real_l = np.var(real_l)
        var_imag_t = np.var(imag_t)
        var_imag_l = np.var(imag_l)

        xs_t = (var_real_t + var_real_l)/(16.0*pi)
        xs_l = ( var_imag_t + var_imag_l)/(16.0*pi)
        
        # Get correction
        correction_t=1.0
        correction_l=1.0
        if corrections_file != "":
            if len(corrections_t) <= t: # Use highest t calculated
                correction_t = corrections_t[len(corrections_t)-1]
                correction_l = corrections_l[len(corrections_l)-1]
            else:
                correction_t = corrections_t[t]
                correction_l = corrections_l[t]

        print tvals[t], correction_t * xs_t + correction_l*xs_l

