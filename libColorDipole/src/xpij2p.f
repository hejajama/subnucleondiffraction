*
* ..File: xpij2p.f 
*
*                           __
* ..The parametrized 3-loop MS singlet splitting functions P^(2) for 
*    the evolution of unpol. singlet parton densities at mu_r = mu_f.
*    The expansion parameter is alpha_s/(4 pi).
*
* ..The distributions (in the mathematical sense) are given as in eq.
*   (B.27) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
*   The name-endings A, B, and C of the functions below correspond to 
*   the kernel superscripts [2], [3], and [1] in that equation.
*
* ..The relative accuracy of these parametrisations, as well as of
*    the convolution results, is better than one part in thousand.

* ..The coefficients of 1/(1-x)_+, (ln x)/x and 1/x are exact (up
*    to a truncation of irrational coefficients).  Furthermore all
*    coefficients written as fractions (e.g., 160./27.D0) are exact.
*    The other terms at x < 1 have fitted to the exact results for x 
*    between 10^-6 and 1 - 10^-6.  The coefficient of delta(1-x) of
*    P_gg^(2) have been slightly adjusted using the second moments.
*
* ..References: S. Moch, J. Vermaseren and A. Vogt,
*               hep-ph/0403192 = Nucl. Phys. B688 (2004) 101
*               A. Vogt, S. Moch and J. Vermaseren,
*               hep-ph/0404111 = Nucl. Phys. B691 (2004) 129
* 
* =====================================================================
*
*
* ..The (regular) pure-singlet splitting functions P_ps^(2).
*    P_qq^(2) is obtained by adding the non-singlet quantity P_NS^(2)+.
*    A parametrization of the latter is provided in the file  xpns2p.f.

       FUNCTION P2PSA (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       P2PS1 = - 3584./(27.D0*Y) * DL - 506.0/ Y + 160./27.D0 * DL**4
     ,         - 400./9.D0 * DL**3 + 131.4 * DL**2 - 661.6 * DL
     ,         - 5.926  * DL1**3 - 9.751 * DL1**2 - 72.11 * DL1
     ,         + 177.4 + 392.9 * Y - 101.4 * Y**2 - 57.04 * DL*DL1
       P2PS2 =   256./(81.*Y) + 32./27.D0 * DL**3 + 17.89 * DL**2
     ,         + 61.75 * DL + 1.778 * DL1**2 + 5.944 * DL1 + 100.1
     ,         - 125.2 * Y + 49.26 * Y**2 - 12.59 * Y**3 
     ,         - 1.889 * DL*DL1 
*
       P2PSA = (1.-Y) * NF * ( P2PS1 + NF * P2PS2 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The gluon->quark splitting functions P_qg^(2).
*
       FUNCTION P2QGA (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER  NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       P2QG1 = - 896./(3.D0*Y) * DL - 1268.3 / Y + 536./27.D0 * DL**4 
     ,         - 44./3.D0 * DL**3 + 881.5 * DL**2 + 424.9 * DL 
     ,         + 100./27.D0 * DL1**4 - 70./9.D0 * DL1**3 
     ,         - 120.5 * DL1**2 + 104.42 * DL1
     ,         + 2522. - 3316.* Y + 2126.* Y**2
     ,         + DL*DL1 * (1823. - 25.22 * DL) - 252.5 * Y*DL**3  
       P2QG2 =   1112./(243.D0*Y) - 16./9.D0 * DL**4 
     ,         - 376./27.D0 * DL**3 - 90.8 * DL**2 - 254.0 * DL   
     ,         + 20./27.D0 * DL1**3 + 200./27.D0 * DL1**2 - 5.496 * DL1
     ,         - 252.0  + 158.0 * Y + 145.4 * Y**2 - 139.28 * Y**3
     ,         - DL*DL1 * ( 53.09  + 80.616 * DL) - 98.07 * Y*DL**2
     ,         + 11.70 * Y*DL**3
* 
       P2QGA = NF * ( P2QG1 + NF * P2QG2 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The quark->gluon splitting functions P_gq^(2).  P2GQ2 is exact.
*
       FUNCTION P2GQA (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       P2GQ0 =   1189.3 * DL/Y + 6163.1 / Y - 4288./81.D0 * DL**4
     ,         + 1568./9.D0 * DL**3 - 1794. * DL**2 + 4033. * DL
     ,         + 400./81.D0 * DL1**4 + 2200./27.D0 * DL1**3
     ,         + 606.3 * DL1**2 + 2193.* DL1 
     ,         - 4307. + 489.3 * Y + 1452.* Y**2 + 146.0 * Y**3
     ,         - 447.3 * DL**2*DL1 - 972.9 * Y*DL**2
       P2GQ1 =   71.082 * DL/Y  - 46.41 / Y + 128./27.D0 * DL**4
     ,         + 704/81.D0 * DL**3 + 20.39  * DL**2 + 174.8 * DL
     ,         - 400./81.D0 * DL1**3 - 68.069 * DL1**2 - 296.7 * DL1
     ,         - 183.8 + 33.35 * Y - 277.9 * Y**2 + 108.6 * Y*DL**2
     ,         - 49.68 * DL*DL1
       P2GQ2 = ( 64. * ( - 1./Y + 1. + 2.* Y)
     ,         + 320.* DL1 * ( 1./Y - 1. + 0.8 * Y)
     ,         + 96.* DL1**2 * ( 1./Y - 1. + 0.5 * Y) ) / 27.D0
*
       P2GQA = ( P2GQ0 + NF * (P2GQ1 + NF * P2GQ2) )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The regular piece of the gluon-gluon splitting function P_gg^(2).
*
       FUNCTION P2GGA (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       P2GGA0 = 2675.8 * DL/Y + 14214./ Y - 144. * DL**4 + 72. * DL**3
     1          - 7471. * DL**2 + 274.4 * DL + 3589. * DL1 - 20852. 
     2          + 3968.* Y - 3363. * Y**2 + 4848. * Y**3 
     3          + DL*DL1 * ( 7305. + 8757. * DL )
       P2GGA1 = 157.27 * DL/Y + 182.96 / Y + 512./27.D0 * DL**4
     1          + 832./9.D0 * DL**3 + 491.3 * DL**2 + 1541. * DL
     2          - 320.0 * DL1 - 350.2 + 755.7 * Y - 713.8 * Y**2 
     3          + 559.3 * Y**3 + DL*DL1 * ( 26.15 - 808.7 * DL )
       P2GGA2 = - 680./(243.D0 * Y) - 32./27.D0 * DL**3 + 9.680 * DL**2
     1          - 3.422 * DL - 13.878 + 153.4 * Y - 187.7 * Y**2 
     2          + 52.75 * Y**3 - DL*DL1 * (115.6 - 85.25* Y + 63.23* DL)
*
       P2GGA = P2GGA0 + NF * ( P2GGA1 + NF * P2GGA2 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The singular piece of the gluon-gluon splitting function P_gg^(2).
*
       FUNCTION P2GGB (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       P2GGB = ( 2643.521 - NF * 412.172 - NF**2 * 16./9.D0 ) / ( 1.-Y) 
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The 'local' piece of the gluon-gluon splitting function P_gg^(2).  
*
       FUNCTION P2GGC (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       DL1 = LOG (1.-Y)
*
       P2GGC =       2643.521 * DL1 + 4425.448 + 0.446
     ,       - NF * ( 412.172 * DL1 +  528.720 + 0.003 )
     ,       + NF**2 * ( - 16./9.D0 * DL1 + 6.4630)
*
       RETURN
       END
*
* =================================================================av==
