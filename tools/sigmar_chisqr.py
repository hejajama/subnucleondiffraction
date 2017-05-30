# utf-8-

# Average the calculated reduced cross section
import sys
sys.path.append("/Users/heikki/lib/")
from matplotlibhelper import *

maxconf = 63
minq2=0
maxq2=500

# prefix
theory = [
          ["jimwlkfit/log/cutoff_exp_0_qsmu_1.2_bp_4.0_N_512/fixed_as_m_0.2/sigmar_cc_as_0.30_", r"m=0.2 b=0, as=0.30"],
          ["jimwlkfit/log/cutoff_exp_0_qsmu_1.2_bp_4.0_N_512/fixed_as_m_0.2/sigmar_cc_as_0.32_", r"m=0.2 b=0, as=0.32"],
          ["jimwlkfit/log/cutoff_exp_0_qsmu_1.2_bp_4.0_N_512/fixed_as_m_0.2/sigmar_cc_as_0.34_", r"m=0.2 b=0, as=0.34"],
          ["jimwlkfit/log/cutoff_exp_0_qsmu_1.2_bp_4.0_N_512/fixed_as_m_0.2/sigmar_cc_as_0.36_", r"m=0.2 b=0, as=0.36"],
          ["jimwlkfit/log/cutoff_exp_0_qsmu_1.2_bp_4.0_N_512/fixed_as_m_0.2/sigmar_cc_as_0.38_", r"m=0.2 b=0, as=0.38"],
          ["jimwlkfit/log/cutoff_exp_0_qsmu_1.2_bp_4.0_N_512/fixed_as_m_0.2/sigmar_cc_as_0.40_", r"m=0.2 b=0, as=0.40"],
          
          ["jimwlkfit/log/cutoff_exp_0_qsmu_1.2_bp_4.0_N_512/fixed_as_m_0.1/sigmar_cc_as_0.12_", r"m=0.1 b=0, as=0.12"],
          ["jimwlkfit/log/cutoff_exp_0_qsmu_1.2_bp_4.0_N_512/fixed_as_m_0.1/sigmar_cc_as_0.14_", r"m=0.1 b=0, as=0.14"],
          ["jimwlkfit/log/cutoff_exp_0_qsmu_1.2_bp_4.0_N_512/fixed_as_m_0.1/sigmar_cc_as_0.16_", r"m=0.1 b=0, as=0.16"],
          ["jimwlkfit/log/cutoff_exp_0_qsmu_1.2_bp_4.0_N_512/fixed_as_m_0.1/sigmar_cc_as_0.18_", r"m=0.1 b=0, as=0.18"],
          ["jimwlkfit/log/cutoff_exp_0_qsmu_1.2_bp_4.0_N_512/fixed_as_m_0.1/sigmar_cc_as_0.20_", r"m=0.1 b=0, as=0.20"],
          ["jimwlkfit/log/cutoff_exp_0_qsmu_1.2_bp_4.0_N_512/fixed_as_m_0.1/sigmar_cc_as_0.22_", r"m=0.1 b=0, as=0.22"],
          ["jimwlkfit/log/cutoff_exp_0_qsmu_1.2_bp_4.0_N_512/fixed_as_m_0.1/sigmar_cc_as_0.24_", r"m=0.1 b=0, as=0.24"],
          
          ["jimwlkfit/log/cutoff_exp_0_qsmu_1.2_bp_2.0_N_512/fixed_as_m_0.2/sigmar_cc_as_0.28_", r"B=2 m=0.2 b=0, as=0.28"],
          ["jimwlkfit/log/cutoff_exp_0_qsmu_1.2_bp_2.0_N_512/fixed_as_m_0.2/sigmar_cc_as_0.30_", r"B=2 m=0.2 b=0, as=0.30"],
          ["jimwlkfit/log/cutoff_exp_0_qsmu_1.2_bp_2.0_N_512/fixed_as_m_0.2/sigmar_cc_as_0.32_", r"B=2 m=0.2 b=0, as=0.32"],
          ["jimwlkfit/log/cutoff_exp_0_qsmu_1.2_bp_2.0_N_512/fixed_as_m_0.2/sigmar_cc_as_0.34_", r"B=2 m=0.2 b=0, as=0.34"],
          ["jimwlkfit/log/cutoff_exp_0_qsmu_1.2_bp_2.0_N_512/fixed_as_m_0.2/sigmar_cc_as_0.36_", r"B=2 m=0.2 b=0, as=0.36"],
          ["jimwlkfit/log/cutoff_exp_0_qsmu_1.2_bp_2.0_N_512/fixed_as_m_0.2/sigmar_cc_as_0.38_", r"B=2 m=0.2 b=0, as=0.38"],
          
     
          
          
          
          #["jimwlkfit/log/cutoff_exp_0_qsmu_1.3_bp_4.0_N_512/fixed_as_m_0.1/sigmar_cc_as_0.14_", r"b=0, m=0.1, as=0.14"],
          #["jimwlkfit/log/cutoff_exp_0_qsmu_1.3_bp_4.0_N_512/fixed_as_m_0.1/sigmar_cc_as_0.16_", r"b=0, m=0.1 , as=0.16"],

          
	 #["jimwlkfit/log/cutoff_exp_0_qsmu_1.2_bp_4.0_N_512/fixed_as_m_0.2/sigmar_cc_as_0.28_", r"b=0, as=0.28"],
          #["jimwlkfit/log/cutoff_exp_5_bp_4.0_qsmu_0.75_N_512_m_0.4/fixed_as_m_0.2/sigmar_cc_as_0.3_", r"b=5, as=0.28"],
          ["jimwlkfit/log/cutoff_exp_5_bp_4.0_qsmu_0.75_N_512_m_0.4/fixed_as_m_0.2/sigmar_cc_as_0.26_", r"b=5, m=0.2, as=0.26"],
          ["jimwlkfit/log/cutoff_exp_5_bp_4.0_qsmu_0.75_N_512_m_0.4/fixed_as_m_0.2/sigmar_cc_as_0.28_", r"b=5, m=0.2, as=0.28"],
          ["jimwlkfit/log/cutoff_exp_5_bp_4.0_qsmu_0.75_N_512_m_0.4/fixed_as_m_0.2/sigmar_cc_as_0.30_", r"b=5, m=0.2, as=0.30"],
          ["jimwlkfit/log/cutoff_exp_5_bp_4.0_qsmu_0.75_N_512_m_0.4/fixed_as_m_0.2/sigmar_cc_as_0.32_", r"b=5, m=0.2, as=0.32"],
           ["jimwlkfit/log/cutoff_exp_5_bp_4.0_qsmu_0.75_N_512_m_0.4/fixed_as_m_0.2/sigmar_cc_as_0.34_", r"b=5, m=0.2, as=0.34"],
          ["jimwlkfit/log/cutoff_exp_5_bp_4.0_qsmu_0.75_N_512_m_0.4/fixed_as_m_0.2/sigmar_cc_as_0.36_", r"b=5, m=0.2, as=0.36"],
          ["jimwlkfit/log/cutoff_exp_5_bp_4.0_qsmu_0.75_N_512_m_0.4/fixed_as_m_0.2/sigmar_cc_as_0.38_", r"b=5, m=0.2, as=0.38"],
          
          ["jimwlkfit/log/cutoff_exp_5_bp_2.0_qsmu_0.75_N_512_m_0.4/fixed_as_m_0.2/sigmar_cc_as_0.26_", r"B=2 b=5, m=0.2, as=0.26"],
          ["jimwlkfit/log/cutoff_exp_5_bp_2.0_qsmu_0.75_N_512_m_0.4/fixed_as_m_0.2/sigmar_cc_as_0.28_", r"B=2 b=5, m=0.2, as=0.28"],
          ["jimwlkfit/log/cutoff_exp_5_bp_2.0_qsmu_0.75_N_512_m_0.4/fixed_as_m_0.2/sigmar_cc_as_0.30_", r"B=2 b=5, m=0.2, as=0.30"],
          ["jimwlkfit/log/cutoff_exp_5_bp_2.0_qsmu_0.75_N_512_m_0.4/fixed_as_m_0.2/sigmar_cc_as_0.32_", r"B=2 b=5, m=0.2, as=0.32"],
          
           ["jimwlkfit/log/cutoff_exp_5_bp_2.0_qsmu_0.75_N_512_m_0.4/fixed_as_m_0.25/sigmar_cc_as_0.32_", r"B=2 b=5, m=0.25, as=0.32"],
          ["jimwlkfit/log/cutoff_exp_5_bp_2.0_qsmu_0.75_N_512_m_0.4/fixed_as_m_0.25/sigmar_cc_as_0.34_", r"B=2 b=5, m=0.25, as=0.34"],
          ["jimwlkfit/log/cutoff_exp_5_bp_2.0_qsmu_0.75_N_512_m_0.4/fixed_as_m_0.25/sigmar_cc_as_0.36_", r"B=2 b=5, m=0.25, as=0.36"],
          ["jimwlkfit/log/cutoff_exp_5_bp_2.0_qsmu_0.75_N_512_m_0.4/fixed_as_m_0.25/sigmar_cc_as_0.38_", r"B=2 b=5, m=0.25, as=0.38"],
          ["jimwlkfit/log/cutoff_exp_5_bp_2.0_qsmu_0.75_N_512_m_0.4/fixed_as_m_0.25/sigmar_cc_as_0.40_", r"B=2 b=5, m=0.25, as=0.40"],
          ["jimwlkfit/log/cutoff_exp_5_bp_2.0_qsmu_0.75_N_512_m_0.4/fixed_as_m_0.25/sigmar_cc_as_0.42_", r"B=2 b=5, m=0.25, as=0.42"],
          
          ["jimwlkfit/log/cutoff_exp_5_bp_2.0_qsmu_0.75_N_512_m_0.4/fixed_as_m_0.1/sigmar_cc_as_0.16_", r"B=2 b=5, m=0.1, as=0.16"],
          ["jimwlkfit/log/cutoff_exp_5_bp_2.0_qsmu_0.75_N_512_m_0.4/fixed_as_m_0.1/sigmar_cc_as_0.18_", r"B=2 b=5, m=0.1, as=0.18"],
          ["jimwlkfit/log/cutoff_exp_5_bp_2.0_qsmu_0.75_N_512_m_0.4/fixed_as_m_0.1/sigmar_cc_as_0.20_", r"B=2 b=5, m=0.1, as=0.20"],
          ["jimwlkfit/log/cutoff_exp_5_bp_2.0_qsmu_0.75_N_512_m_0.4/fixed_as_m_0.1/sigmar_cc_as_0.22_", r"B=2 b=5, m=0.1, as=0.22"],
          ]


for th in theory:
    xvals=[]
    q2vals=[]
    hera_sigmarvals=[]
    hera_sigmarerrs = []
    theory_sigmar = []

    # read all files and average
    for c in range(maxconf):
        fname = th[0] + str(c)
        tmpxlist=[]
        tmpexplist=[]
        tmpexperrlist=[]
        tmptheorylist=[]
        tmpq2list=[]
        tmplist=[]
        try:
            readfile_xy(fname, tmpxlist, tmpq2list, xcol=1, ycol=0 )
            readfile_xy(fname, tmpexplist, tmpexperrlist, xcol=3, ycol=4)
            readfile_xy(fname, tmptheorylist, tmplist, xcol=5, ycol=0)
        except IOError:
            print >> sys.stderr, "IOError, Skipping file " + fname
            continue

        if xvals == []:
            xvals = tmpxlist
            q2vals = tmpq2list
            hera_sigmarerrs = tmpexperrlist
            hera_sigmarvals = tmpexplist
            for i in range(len(tmptheorylist)):
                theory_sigmar.append([tmptheorylist[i]])
        else:
            for i in range(len(tmptheorylist)):
                try:
                    if xvals[i] != tmpxlist[i] or q2vals[i] != tmpq2list[i]:
                        print >> sys.stderr, "ERROR! x,Q2 vals don't match, fname " + fname
                    theory_sigmar[i].append(tmptheorylist[i])
                except:
                    print "Error with file " + fname
    if len(hera_sigmarvals)==0:
        print "No data with files " + th[0]
        continue

    # Calculate chi^2
    chisqr=0
    points=0
    average_theory_sigmar = []
    for i in range(len(theory_sigmar)):
        average_theory_sigmar.append(mean(theory_sigmar[i]))
        if q2vals[i] <= maxq2 and q2vals[i] >= minq2:
            chisqr += pow( (hera_sigmarvals[i] - average_theory_sigmar[-1]) / hera_sigmarerrs[i], 2.0)
            points = points+1

    print th[1] + ", Chi^2/N=" + str(chisqr/points) + ", N=" + str(points)
