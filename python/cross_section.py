'''
Calculate cross sections reading in the data generated using the subnucleondiffraction code

This code supports both the new and old output format for the subnucleondiffraction code.
The format changed in 11/2024. By default the new format is used.

Heikki Mäntysaari <heikki.mantysaari@jyu.fi>
'''

import numpy as np
import scipy
from typing import Tuple
import argparse

########## Cross section calculation

GEVSQRTONB = 1.0e7/(5.068*5.068)
GEVSQRTOMB=GEVSQRTONB*1e-6


def CoherentCrossSection(dirname:str,minconf:int=0,maxconf:int=250, amplitude:bool=False, file_format:str="new", 
                         polarization:str="T", show_warnings=True) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    '''Coherent cross section in mb, imaginary part neglected
    Returns tvals [GeV^2], dsigma/dt [mb/GeV^2], staterrr

    polarization: T or L (transverse, longitudinal)
    
    file_format: new or old
        new: both real and imaginary part in the same file
        old: real and imaginary part in separate files
    '''

    xsvals=[]
    tvals=[]
    for i in range(minconf,maxconf+1):
        #try:
        if file_format == "old":
            fn=dirname+"/real/spectra_"+str(i)
            if polarization=="T":
                column=1
            elif polarization=="L":
                column=2
        else:
            fn=dirname+"/spectra_"+str(i)
            if polarization=="T":
                column=1
            elif polarization=="L":
                column=3

        try:
            dat=np.loadtxt(fn)
        except IOError:
            if show_warnings:
                print("# File ",fn," does not exist, skip")
            continue
        except ValueError:
            if show_warnings:
                print("# Error with file ",fn)
            continue
        except Exception as e:
            if show_warnings:
                print("# Error with file ",fn, " ",e)
            continue
            
        if len(dat)==0 or np.isnan(dat).any() or np.isinf(dat).any():
            if show_warnings:
                print("# Skip file with NaN or Inf entries " , i)
                print("# Skip emtpy file " , i)
            continue
    
        if len(tvals)==0:
            tvals=dat[:,0]
        else:
            if dat.ndim != 2:
                if show_warnings:
                    print("# Wrong shape ", dat.shape , " with file " , fn)
                continue
            if len(tvals) != len(dat[:,0]):
                if show_warnings:
                    print("# Error with file - t values do not match " + fn)
                continue
            if not np.allclose(tvals, dat[:,0], atol=1e-5):
                if show_warnings:
                    print("# Error with file - t values do not match " + fn)
                continue
                    
        # Potentially empty file?
        if np.abs(scipy.integrate.simpson(dat[:,column]**2, x=dat[:,0]))<1e-25:
            print("# Zero cross section, dir ", dirname, " config ", i)
            continue
            
        xsvals.append(dat[:,column])
            
        
    if amplitude:
        return tvals,np.mean(xsvals, axis=0), scipy.stats.sem(xsvals, axis=0)
    
    xs = np.square(np.mean(xsvals,axis=0))/(16.0*np.pi)*GEVSQRTONB/1000/1000
    #staterr=scipy.stats.sem(np.square(xsvals),axis=0)/(16.0*np.pi)*GEVSQRTONB/1000/1000
    staterr=np.abs(2.0*np.array(np.mean(xsvals,axis=0))*scipy.stats.sem(xsvals,axis=0)/(16.0*np.pi)*GEVSQRTONB/1000/1000)
    
    #print("staterr shape",staterr.shape)
    
    return tvals,xs,staterr
    

def IncoherentCrossSection(dirname: str, minconf: int = 0, maxconf: int = 400, 
                           file_format: str = "new", polarization: str = "T",
                           show_warnings=True) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    '''Coherent cross section in mb
    Returns tvals [GeV^2], dsigma/dt [mb/GeV^2], staterrr
    
    polarization: T or L (transverse, longitudinal)
    
    file_format: new or old
        new: both real and imaginary part in the same file
        old: real and imaginary part in separate files
    '''
    tvals = []
    data_r = []
    data_i = []
    for i in range(minconf, maxconf + 1):  
        if file_format == "old":
            fn = dirname + "/real/spectra_" + str(i)
            if polarization == "T":
                column = 1
            elif polarization == "L":
                column = 2
        else:
            fn = dirname + "/spectra_" + str(i)
            if polarization == "T":
                column = 1
            elif polarization == "L":
                column = 3

        try:  
            tmpdata = np.loadtxt(fn)
        except IOError:
            if show_warnings:
                print("# File ", fn, " does not exist, skip")
            continue
        if len(tmpdata) == 0 or np.isnan(tmpdata).any() or np.isinf(tmpdata).any():
            if show_warnings:
                print("# Skip file with NaN or Inf entries " + fn)
            continue

        tmptvals = tmpdata[:, 0]
        
        if len(tvals) == 0:
            tvals = tmptvals
        else:
            if tmpdata.ndim != 2:
                if show_warnings:
                    print("# Wrong shape ", tmpdata.shape , " with file " , fn)
                continue
            if len(tmptvals) != len(tvals):
                if show_warnings:
                    print("# Error with file - t values do not match " + fn)
                continue
            if not np.allclose(tmptvals, tvals, atol=1e-5):
                if show_warnings:
                    print("# Error with file - t values do not match " + fn)
                continue
                    
        # Potentially empty file?
        if np.abs(scipy.integrate.simpson(tmpdata[:, column]**2, x=tmpdata[:, 0])) < 1e-25:
            if show_warnings:
                print("Inoch: zero cross section, dir ", dirname, " config ", i)
            continue
        
        if file_format == "old":
            ## imag
            try:
                tmpdata_i = np.loadtxt(dirname + "/imag/spectra_" + str(i))
            except IOError:
                print("File ", dirname + "/imag/spectra_" + str(i), " does not exist, but real part was there, RESULTS MAY NOT MAKE SENSE")
                return np.array([]), np.array([]), np.array([])
        
            if len(tmpdata_i) == 0:
                continue
            
            # Potentially empty file?
            if np.abs(scipy.integrate.simpson(tmpdata_i[:, 1]**2, x=tmpdata_i[:, 0])) < 1e-25:
                print("# Inoch: zero cross section, dir ", dirname, " config ", i)
                continue
            
            tmptvals_i = tmpdata[:, 0]

            if not np.allclose(tvals, tmpdata[:, 0], atol=1e-5):
                print("# Error with file - t values do not match " + fn)
                continue
            
            data_i.append(tmpdata_i)

        data_r.append(tmpdata) # new file format: this contains everything

    data_r = np.array(data_r)
    data_i = np.array(data_i)
    
    if file_format == "old":
        var = np.var(data_r[:, :, column], axis=0) + np.var(data_i[:, :, column], axis=0)
        halfvar1 = np.var(data_r[0:maxconf // 2, :, column], axis=0) + np.var(data_i[0:maxconf // 2, :, column], axis=0)
        halfvar2 = np.var(data_r[maxconf // 2:, :, column], axis=0) + np.var(data_i[maxconf // 2:, :, column], axis=0)
    else:
        var = np.var(data_r[:, :, column], axis=0) + np.var(data_r[:, :, column + 1], axis=0)
        halfvar1 = np.var(data_r[0:maxconf // 2, :, column], axis=0) + np.var(data_r[0:maxconf // 2, :, column + 1], axis=0)
        halfvar2 = np.var(data_r[maxconf // 2:, :, column], axis=0) + np.var(data_r[maxconf // 2:, :, column + 1], axis=0)

    err1 = np.abs(halfvar1 - var)
    err2 = np.abs(halfvar2 - var)
    errs = (err1 + err2) * 0.5
    
    return tvals, var / (16.0 * np.pi) * GEVSQRTONB / 1000 / 1000, errs / (16.0 * np.pi) * GEVSQRTONB / 1000 / 1000
    

######### Kinematics

def MesonMass(meson):
    if meson=="rho": 
        return 0.776
    elif meson=="jpsi":
        return 3.097
    elif meson=="psi_2s":
        return 3.68609
    elif meson=="upsilon_1s" or meson=="upsilon":
        return 9.460
    elif meson=="upsilon_2s":
        return 10.023
    elif meson=="upsilon_3s":
        return 10.355
    else:
        print("Unknown meson " + meson)
        return 0
    

def xpom(meson: str, Q2: float, W: float) -> float:
    '''xp for the given meson at fixed Q^2 and W^2
    Assumes t=0
    '''
    M2 = MesonMass(meson)**2
    
    return (M2 + Q2) / (W**2 + Q2)

def JIMWLK_steps(meson: str, Q2: float, W: float, alphas: float = -1, ds: float=-1, x0: float = 0.01) -> int:
    '''Number of JIMWLK evolution steps (when using our JIMWLK code), epending on initial x0 and 
    fixed/running coupling
    
    alphas: -1 for running coupling, otherwise value of the fixed coupling
    x0: initial x for the JIMWLK evolution
    W: center of mass energy [GeV]
    Q2: photon virtuality [GeV^2]
    ds: step size specified in the JIMWLK code. If <0, default is used
    meson: meson name
    '''
    xp = xpom(meson, Q2, W)
    y = np.log(x0 / xp)   
    
    if alphas > 0:
        if ds < 0:
            ds = 0.0004
        s = alphas * y / (ds * np.pi * np.pi)
    else:
        # RC
        if ds < 0:
            ds = 0.004
        s = y / (ds * np.pi * np.pi)
    return int(s+0.5)

def WFromSteps(steps: int, meson: str, alphas: float = -1, ds: float = -1, x0: float = 0.01) -> float:
    '''Compute center-of-mass energy W from JIMWLK evolution steps
    
    steps: number of JIMWLK evolution steps
    meson: meson name
    alphas: -1 for running coupling, otherwise value of the fixed coupling
    x0: initial x for the JIMWLK evolution
    ds: step size specified in the JIMWLK code. If <0, default is used
    '''
    if alphas > 0:
        if ds < 0:
            ds = 0.0004
        y = steps * ds * np.pi * np.pi / alphas
    else:
        if ds < 0:
            ds = 0.004
        y = steps * ds * np.pi * np.pi
        
    xp = x0 * np.exp(-y)
    W = np.sqrt(MesonMass(meson)**2 / xp)
    return W


def W_rapidity(y, meson="jpsi", sqrts=5020):
    Mv=MesonMass(meson)
        
    return np.sqrt(sqrts*Mv*np.exp(y))



########## Visual aspects
def GetStyle(experiment):
    if experiment=="alice":
        return "."
    elif experiment=="cms":
        return "^"
    elif experiment=="lhcb":
        return "p"
    elif experiment=="star":
        return "x"
    else:    
        print("Error: experiment not recognized")
        return "o"
    
def GetColor(experiment):
    if experiment=="alice":
        return "black"
    elif experiment=="cms":
        return "blue"
    elif experiment=="lhcb":
        return "red"
    elif experiment=="star":
        return "red"
    else:    
        print("Error: experiment not recognized")
        return "black"

def GetLabel(experiment):
    if experiment=="alice":
        return "ALICE"
    elif experiment=="cms":
        return "CMS"
    elif experiment=="lhcb":
        return "LHCb"
    elif experiment=="star":
        return "STAR"
    else:    
        print("Error: experiment not recognized")
        return "Unknown"
    


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Calculate coherent and incoherent cross sections.")
    parser.add_argument("-dir", type=str, required=True, help="Directory containing the data files.")
    parser.add_argument("-maxconf", type=int, required=True, help="Maximum number of configurations.")
    parser.add_argument("-ds", type=float, default=-1, help="Step size specified in the JIMWLK code.")
    parser.add_argument("-steps", type=int, default=0, help="Number of JIMWLK evolution steps.")
    args = parser.parse_args()

    dirname = args.dir
    maxconf = args.maxconf

    tvals_coh, xs_coh, staterr_coh = CoherentCrossSection(dirname, maxconf=maxconf)
    tvals_incoh, xs_incoh, staterr_incoh = IncoherentCrossSection(dirname, maxconf=maxconf)

    print("# Cross sections at W="+str(WFromSteps(steps=args.steps, meson="jpsi", ds=args.ds))+" GeV")
    print("# t [GeV^2] cohxs [mb/GeV^2] staterr_coh [mb/GeV^2] incohxs [mb/GeV^2] staterr_incoh [mb/GeV^2]")
    for t, cohxs, cohxserr, incohxs, incohxserr in zip(tvals_coh, xs_coh, staterr_coh, xs_incoh, staterr_incoh):
        print(t, cohxs, cohxserr, incohxs, incohxserr)
