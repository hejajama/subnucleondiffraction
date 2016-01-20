*
* ..File: xpns2p.f 
*
*                           __
* ..The parametrized 3-loop MS non-singlet splitting functions P^(2)
*    for the evolution of unpolarized partons densities, mu_r = mu_f.
*    The expansion parameter is alpha_s/(4 pi).
*
* ..The distributions (in the mathematical sense) are given as in eq.
*    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
*    The name-endings A, B, and C of the functions below correspond to
*    the kernel superscripts [2], [3], and [1] in that equation.
*
*  ..The relative accuracy of these parametrizations, as well as of
*    the convolution results, is better than one part in thousand.
*
* ..References: S. Moch, J. Vermaseren and A. Vogt,
*               hep-ph/0209100 = Nucl. Phys. B646 (2002) 181,
*               hep-ph/0403192 = Nucl. Phys. B688 (2004) 101
*
* =====================================================================
*
*
* ..This is the regular piece of P2_NS+.  The rational coefficients are
*    exact, the rest has been fitted for x between 10^-6 and 1 - 10^-6.
*    The N_f^2 part is exact and was first determined in N-space by 
*    J.A. Gracey in Phys. Lett. B322 (1994) 141.
*
       FUNCTION P2NSPA (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
       D81 = 1./81.D0
*
       P2NSPA =   1641.1 - 3135.* Y + 243.6 * Y**2 - 522.1 * Y**3
     ,            + 128.*D81 * DL**4 + 2400.*D81 * DL**3
     ,            + 294.9 * DL**2 + 1258.* DL
     ,            + 714.1 * DL1 + DL*DL1 * (563.9 + 256.8 * DL)
     ,        + NF * ( -197.0 + 381.1 * Y + 72.94 * Y**2 + 44.79 * Y**3
     ,            - 192.*D81 * DL**3  - 2608.*D81 * DL**2 - 152.6 * DL
     ,            - 5120.*D81 * DL1 - 56.66 * DL*DL1 - 1.497 * Y*DL**3 )
     ,        + NF**2 * ( 32.* Y*DL/(1.-Y) * (3.* DL + 10.) + 64.
     ,            + (48.* DL**2 + 352.* DL + 384.) * (1.-Y) ) * D81
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the regular piece of P2_NS-.  The rational coefficients are 
*    exact, the rest has been fitted for x between 10^-6 and 1 - 10^-6.
*    The N_f^2 part is exact (and identical to that of P2_NS+). 
*
       FUNCTION P2NSMA (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
       D81 = 1./81.D0
*
       P2NSMA =   1860.2 - 3505.* Y + 297.0 * Y**2 - 433.2 * Y**3
     ,            + 116.*D81 * DL**4 + 2880.*D81 * DL**3 
     ,            + 399.2 * DL**2 + 1465.2 * DL
     ,            + 714.1 * DL1 + DL*DL1 * (684.0 + 251.2 * DL)
     ,        + NF * ( -216.62 + 406.5 * Y + 77.89 * Y**2 + 34.76 * Y**3
     ,            - 256.*D81 * DL**3  - 3216.*D81 * DL**2 - 172.69 * DL 
     ,            - 5120.*D81 * DL1 - 65.43 * DL*DL1 - 1.136 * Y*DL**3 )
     ,        + NF**2 * ( 32.* Y*DL/(1.-Y) * (3.* DL + 10.) + 64.
     ,            + (48.* DL**2 + 352.* DL + 384.) * (1.-Y) ) * D81
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the singular piece of both P2_NS+ and P2_NS-. It is exact 
*    up to the truncation of the irrational coefficients.
*
       FUNCTION P2NSB (Y, NF)
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       P2NSB = ( 1174.898 - NF * 183.187 - NF**2 * 64./81.D0 ) / (1.-Y)
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' piece of P2_NS+. The coefficients of delta(1-x)
*    have been partly shifted relative to the exact (truncated) values.
*
       FUNCTION P2NSPC (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       DL1 = LOG (1.-Y)
*
       P2NSPC =       1174.898 * DL1 + 1295.624 - 0.24
     ,        - NF * ( 183.187 * DL1 + 173.938 - 0.011 )
     ,        + NF**2 * ( - 64./81.D0 * DL1 + 1.13067 )
*
       RETURN
       END
*
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' piece of P2_NS-. The coefficients of delta(1-x) 
*    have been partly shifted relative to the exact (truncated) values.
*
       FUNCTION P2NSMC (Y, NF)
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
       DL1 = LOG (1.-Y)
*
       P2NSMC =       1174.898 * DL1 + 1295.624 - 0.154
     ,        - NF * ( 183.187 * DL1 + 173.938  - 0.005 )
     ,        + NF**2 * ( - 64./81.D0 * DL1 + 1.13067 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..This is P2_NSS, the difference of P2_NSV and P2_NS-.
*
       FUNCTION P2NSSA (Y, NF)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       D27 = 1./27.D0
       DL  = LOG (Y)
       Y1  = 1.- Y
       DL1 = LOG (Y1)
*
       P2NSSA = Y1* ( 151.49 + 44.51 * Y - 43.12 * Y**2 + 4.820 * Y**3 )
     1          + 40.*D27 * DL**4 - 80.*D27 * DL**3 + 6.892 * DL**2 
     2          + 178.04 * DL + DL*DL1 * ( - 173.1 + 46.18 * DL )
     4          + Y1*DL1 * ( - 163.9 / Y - 7.208 * Y ) 
*
       P2NSSA  = NF * P2NSSA
*
       RETURN
       END
*
* =================================================================av==
