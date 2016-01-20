
*-MB--------------------------------------------------------------------
*    This is the parameterization by Andreas Vogt of the coefficient
*    \~{A}^{S,(2)}_{Hg} eq. (B.3) in Buza et al., (BSMN, hep-ph/9612398).
*    The parameterization is made for the case m^2 = mu^2.  
*-MB--------------------------------------------------------------------

*
* ..File: xa2hgp.f  
*
*
* ..Calculation of the alpha_s^2 heavy-quark singlet operator matrix 
*    element (OME) a^2_Hg(x) in the MS(bar) scheme for mu_f = m_H via 
*    a compact parametrization involving only logarithms.
*    The coupling constant is normalized as  a_s = alpha_s/(4*pi).
*
* ..This quantity, presented in Appendix B of M. Buza, Y. Matiounine,
*    J. Smith and W.L. van Neerven, Eur. Phys. J. C1 (1998) 301 (BSMN), 
*    is required for the N_f matching of the NNLO parton densities.
*
*  ..The distributions (in the mathematical sense) are given as in eq.
*    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
*    The name-endings A, B, and C of the functions below correspond to
*    the kernel superscripts [2], [3], and [1] in that equation.
*
*  ..The relative accuracy of the OME, as well as of its convolutions
*    with the gluon distribution, amounts to a few thousandth.
*
*  ..Reference: A. Vogt, 
*               hep-ph/0408244 = Comput. Phys. Commun. 170 (2005) 65
*             
*  ..The user should also cite the original calculation by BSMN.
*
*
* =====================================================================
*
* ..This is the regular piece.
*
       FUNCTION A2HGA (Y)
       IMPLICIT REAL*8 (A-Z)
*
       DL  = LOG (Y)
       DL1 = LOG (1.-Y)
*
       A2HGA = - 24.89 / Y - 187.8 + 249.6 * Y - 146.8 * DL**2 * DL1  
     1         - 1.556 * DL**3  - 3.292 * DL**2  - 93.68 * DL
     2         - 1.111 * DL1**3 - 0.400 * DL1**2 - 2.770 * DL1
*
       RETURN
       END
* 
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' piece, which has no counterpart in WvN's
*    program, as it does not exist in the exact expressions. 
*
       FUNCTION A2HGC (Y)
       IMPLICIT REAL*8 (A-Z)
*
       DUM   = Y          !avoid compiler warning
       A2HGC = - 0.006  
*
       RETURN
       END
*
* =================================================================av==
