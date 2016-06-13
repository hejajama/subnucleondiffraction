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
        print >> sys.stderr, "#File not found: " + fname_real + " or " + fname_imag
        continue
    except ValueError, e:
        print >> sys.stderr, "Error while reading file " + fname_real + " or " + fname_imag +", error: " + str(e)
        sys.exit(1)

    # assume that all files go to largest t
    if tvals == []:
        tvals = tvals_real
    
    # Save amplitude values as a function of t
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

if len(realparts) != len(imagparts):
    print "ERROR! different number of real and imaginary parts!"
    sys.exit(-1)

print "# Read " + str(len(realparts)) + " amplitudes from dir " + dir


# Read corrections
corrections_t = []
corrections_l = []
tmp=[]
tvals_corrections=[]
if corrections_file != "":
    readfile_xy(corrections_file, tvals_corrections, corrections_t)
    readfile_xy(corrections_file, tmp, corrections_l, ycol=2)
    # Check that tvals match
    for i in range(len(tvals_corrections)):
        if tvals_corrections[i] != tvals[i]:
            print "tvals do not agree at index " + str(i) + ", priting tvals and tvals_corrections"
            print str(tvals)
            print str(tvals_corrections)
            sys.exit(-1)

# Ok, calculate coheret xs
if coherent:
    # <A_T>^2 + <A_L>^2
    for t in range(len(realparts[0])):
        sum_real_t=0
        sum_real_l=0
        sum_imag_t=0
        sum_imag_l=0
        nconfs=0
        amp_t_real = []    # to calculate stdev
        amp_t_imag = []
        amp_l_real = []
        amp_l_imag = []
        xs_t=[]   # cross sections from individual files
        xs_l=[]
        for conf in range(len(realparts)):
            # check if the calculations is done at high t
            if len(realparts[conf]) <= t or len(imagparts[conf])<=t:
                print "# Skip file " + fnames[conf] + " at t " + str(tvals[t])
                continue
            try:
                #sum_real_t += realparts[conf][t][0]
                #sum_real_l += realparts[conf][t][1]
                #sum_imag_t += imagparts[conf][t][0]
                #sum_imag_l += imagparts[conf][t][1]
                amp_t_real.append(realparts[conf][t][0])
                amp_t_imag.append(imagparts[conf][t][0])
                amp_l_real.append(realparts[conf][t][1])
                amp_l_imag.append(imagparts[conf][t][1])
                xs_t.append(sqr(amp_t_real[-1]) + sqr(amp_t_imag[-1]))
                xs_l.append(sqr(amp_l_real[-1])+ sqr(amp_l_imag[-1]))
            #print conf,t,realparts[conf][t][0], sum_real_t
            except IndexError:  # Some calculation was not done at high enough t, skip
                print "# WHY INDEX ERROR at " + str(conf) + " at t " + str(tvals[t])
                continue
            nconfs += 1
        
        transverse =  sqr(mean(amp_t_real)) + sqr(mean(amp_t_imag))
        longitudinal =sqr(mean(amp_l_real)) + sqr(mean(amp_l_imag))
        
        

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
        
        # staterr
        #staterr_t = sqr(np.std( amp_t_real ) / sqrt(len(amp_t_real))) + sqr(np.std( amp_t_imag ) / sqrt(len(amp_t_imag)))
        #staterr_t *= correction_t
        #staterr_l =sqr(np.std( amp_l_real ) / sqrt(len(amp_l_real))) + sqr(np.std( amp_l_imag ) / sqrt(len(amp_l_imag)))
        #staterr_l *= correction_l
        
        #staterr_t = correction_t * np.std(xs_t)/sqrt(len(xs_t))
        # staterr_l =correction_l * np.std(xs_l)/sqrt(len(xs_l))

        staterr_t = correction_t*sqrt( sqr(2.0*mean(amp_t_real)*np.std(amp_t_real)/sqrt(len(amp_t_real)) ) + sqr(2.0*mean(amp_t_imag)*np.std(amp_t_imag)/sqrt(len(amp_t_imag)) ) )
        staterr_l = correction_l*sqrt( sqr(2.0*mean(amp_l_real)*np.std(amp_l_real)/sqrt(len(amp_l_real)) ) + sqr(2.0*mean(amp_l_imag)*np.std(amp_l_imag)/sqrt(len(amp_l_imag)) ) )

#print "old, new staterr: ", staterr_t, err_t
        
        
        # err = max/min
        max_xs_t =correction_t*(max(xs_t))
        max_xs_l =correction_l*(max(xs_l))
        min_xs_t =correction_t*(min(xs_t))
        min_xs_l =correction_l*(min(xs_l))
        
        max_xs = (max_xs_t + max_xs_l)/(16.0*pi)
        min_xs = (min_xs_t + min_xs_l)/(16.0*pi)
        
        xs =(correction_t*transverse+correction_l*longitudinal)/(16.0*pi)
        
        #staterr = 2.0*transverse*staterr_t + 2.0*longitudinal*staterr_l
        staterr= staterr_t + staterr_l

        print tvals[t], xs, max_xs, min_xs, staterr/(16.0*pi)





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

        xs_t = (var_real_t + var_imag_t)/(16.0*pi)
        xs_l = ( var_real_l + var_imag_l)/(16.0*pi)
        
        # Jackknife resampling to get error estimates
        jackknife_xs_t_real=[]
        jackknife_xs_t_imag=[]
        jackknife_xs_l_real=[]
        jackknife_xs_l_imag=[]
        for skip_i in range(len(realparts)):
            i=0
            tmplist_t_r = []
            tmplist_t_i = []
            tmplist_l_r = []
            tmplist_l_i = []
            for rc,ic in zip(realparts, imagparts):
                if i != skip_i:
                    try:
                        tmplist_t_r.append(rc[t][0])
                        tmplist_t_i.append(ic[t][0])
                        tmplist_l_r.append(rc[t][1])
                        tmplist_l_i.append(ic[t][1])
                    except:
                        print >> sys.stderr, "Error at t=" + str(t) #+ ", rc " + str(rc) +", ic " + str(ic)
                i=i+1
            
            jackknife_xs_t_real.append( np.var(tmplist_t_r))
            jackknife_xs_t_imag.append( np.var(tmplist_t_i))
            jackknife_xs_l_real.append( np.var(tmplist_l_r))
            jackknife_xs_l_imag.append( np.var(tmplist_l_i))
            
        # Error estimate
        #err_t = np.std(jackknife_xs_t)/sqrt(len(jackknife_xs_t))/(16.0*pi)
        #err_l = np.std(jackknife_xs_l)/sqrt(len(jackknife_xs_l))/(16.0*pi)
        
        # properly
        var_t_real=0
        var_t_imag=0
        var_l_real=0
        var_l_imag=0
        for i in range(len(jackknife_xs_t_real)):
            var_t_real += sqr( jackknife_xs_t_real[i] - var_real_t)
            var_t_imag += sqr( jackknife_xs_t_imag[i] - var_imag_t)
            var_l_real += sqr( jackknife_xs_l_real[i] - var_real_l)
            var_l_imag += sqr( jackknife_xs_l_imag[i] - var_imag_l)

        


        coef = (len(jackknife_xs_t_real)-1.0)/len(jackknife_xs_t_real)
        var_t =  coef * var_t_real + coef * var_t_imag
        var_l = coef * var_l_real + coef * var_l_imag

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

# TODO: ASSUMING NOW ONLY TRANSVERSE!!!!
        #err = correction_t * sqrt( var_t / len(jackknife_xs_t_real)) / (16.0*pi)
        err = (correction_t * sqrt(var_t_real/len(jackknife_xs_t_real)) + correction_t*sqrt(var_t_imag/len(jackknife_xs_t_real)))/(16.0*pi)

        print tvals[t], correction_t * xs_t + correction_l*xs_l, err

