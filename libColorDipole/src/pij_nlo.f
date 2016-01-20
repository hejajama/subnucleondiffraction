
C     ========================================
      DOUBLE PRECISION FUNCTION XF1TFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (10)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      XF1TFUNC = 0!X * ( FF1TFUNC(X,NF) + FG1TFUNC(X,NF) )

      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION XG1TFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (10)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      XG1TFUNC = X * ( GG1TFUNC(X,NF) + GF1TFUNC(X,NF) )

      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION PP1SFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio NPB175(1980)27

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
!       include 'qconst.inc'
!  
!       CX2    = X**2
!       C1PX   = 1.+X
!       C1MX   = 1.-X
!       CLX    = LOG(X)
!       CLX2   = CLX**2
!       CL1MX  = LOG(C1MX)
!       CL1PX  = LOG(C1PX)
!       CPFFX  = (1.+CX2) / C1MX
!       CPFFMX = (1.+CX2) / C1PX
!       CS3X   = -DMB_DILOG(-X)
!       CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)
!  
! C--   P_F(x) in FP eq. (4.52)
!       AAA = - CPFFX*2.*CLX*CL1MX - (2.*X+3./C1MX)*CLX - .5*C1PX*CLX2
!      +      - 5.*C1MX
! C--   0.5*P_G(x) in FP eq. (4.53)     
!       BBB =   CPFFX*(.5*CLX2+C11S6*CLX+CPIA) + C1PX*CLX + C20S3*C1MX
! C--   P_N(x) in FP eq. (4.54)   
!       CCC = - CPFFX*C2S3*(C5S3+CLX) - C4S3*C1MX
! C--   P_QQ(x) in FP eq. (4.50) 
!       PQQ  = C16S9*AAA + 4.*BBB + C2S3*NF*CCC
! C--   0.5*P_A(x) in FP eq. (4.55) and P_QQBAR(x) eq. (4.51)      
!       PQQB = - C4S9 * ( CPFFMX*CS2X + C1PX*CLX + 2.*C1MX )
!  
!       PP1SFUNC = PQQ + PQQB
      PP1SFUNC = 0
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION PM1SFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio NPB175(1980)27

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
!       include 'qconst.inc'
!  
!       CX2    = X**2
!       C1PX   = 1.+X
!       C1MX   = 1.-X
!       CLX    = LOG(X)
!       CLX2   = CLX**2
!       CL1MX  = LOG(C1MX)
!       CL1PX  = LOG(C1PX)
!       CPFFX  = (1.+CX2) / C1MX
!       CPFFMX = (1.+CX2) / C1PX
!       CS3X   = -DMB_DILOG(-X)
!       CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)
! 
! C--   P_F(x) in FP eq. (4.52) 
!       AAA = - CPFFX*2.*CLX*CL1MX - (2.*X+3./C1MX)*CLX - .5*C1PX*CLX2
!      +      - 5.*C1MX
! C--   0.5*P_G(x) in FP eq. (4.53)     
!       BBB =   CPFFX*(.5*CLX2+C11S6*CLX+CPIA) + C1PX*CLX + C20S3*C1MX
! C--   P_N(x) in FP eq. (4.54)      
!       CCC = - CPFFX*C2S3*(C5S3+CLX) - C4S3*C1MX
! C--   P_QQ(x) in FP eq. (4.50) 
!       PQQ  = C16S9*AAA + 4.*BBB + C2S3*NF*CCC
! C--   0.5*P_A(x) in FP eq. (4.55) and P_QQBAR(x) eq. (4.51)      
!       PQQB = - C4S9 * ( CPFFMX*CS2X + C1PX*CLX + 2.*C1MX )
!  
!       PM1SFUNC = PQQ - PQQB
      PM1SFUNC = 0
      RETURN
      END
      
C     ========================================
      DOUBLE PRECISION FUNCTION PP1TFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio NPB175(1980)27

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
!       include 'qconst.inc'
!  
!       CX2    = X**2
!       C1PX   = 1.+X
!       C1MX   = 1.-X
!       CLX    = LOG(X)
!       CLX2   = CLX**2
!       CL1MX  = LOG(C1MX)
!       CL1PX  = LOG(C1PX)
!       CPFFX  = (1.+CX2) / C1MX
!       CPFFMX = (1.+CX2) / C1PX
!       CS3X   = -DMB_DILOG(-X)
!       CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)
!  
! C--   P_F(x) in FP eq. (4.52)
!       AAA = - CPFFX*2.*CLX*CL1MX - (2.*X+3./C1MX)*CLX - .5*C1PX*CLX2
!      +      - 5.*C1MX
! C--   0.5*P_G(x) in FP eq. (4.53)     
!       BBB =   CPFFX*(.5*CLX2+C11S6*CLX+CPIA) + C1PX*CLX + C20S3*C1MX
! C--   P_N(x) in FP eq. (4.54)   
!       CCC = - CPFFX*C2S3*(C5S3+CLX) - C4S3*C1MX
! C--   DEL(x) in FP eq. (6.40)      
!       DDD = 4.D0*CPFFX*CLX*CL1MX + (6.D0/C1MX-5.D0-X)*CLX
!      +      + (C1PX-2.D0*CPFFX)*CLX2
! C--   P_QQ(x) in FP eq. (4.50) 
!       PQQ  = C16S9*AAA + 4.*BBB + C2S3*NF*CCC
! C--   0.5*P_A(x) in FP eq. (4.55) and P_QQBAR(x) eq. (4.51)      
!       PQQB = - C4S9 * ( CPFFMX*CS2X + C1PX*CLX + 2.*C1MX )
!  
!       PP1TFUNC = PQQ + PQQB + C16S9*DDD
      PP1TFUNC = 0
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION PM1TFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio NPB175(1980)27

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
!       include 'qconst.inc'
!  
!       CX2    = X**2
!       C1PX   = 1.+X
!       C1MX   = 1.-X
!       CLX    = LOG(X)
!       CLX2   = CLX**2
!       CL1MX  = LOG(C1MX)
!       CL1PX  = LOG(C1PX)
!       CPFFX  = (1.+CX2) / C1MX
!       CPFFMX = (1.+CX2) / C1PX
!       CS3X   = -DMB_DILOG(-X)
!       CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)
! 
! C--   P_F(x) in FP eq. (4.52) 
!       AAA = - CPFFX*2.*CLX*CL1MX - (2.*X+3./C1MX)*CLX - .5*C1PX*CLX2
!      +      - 5.*C1MX
! C--   0.5*P_G(x) in FP eq. (4.53)     
!       BBB =   CPFFX*(.5*CLX2+C11S6*CLX+CPIA) + C1PX*CLX + C20S3*C1MX
! C--   P_N(x) in FP eq. (4.54)      
!       CCC = - CPFFX*C2S3*(C5S3+CLX) - C4S3*C1MX
! C--   DEL(x) in FP eq. (6.40)      
!       DDD = 4.D0*CPFFX*CLX*CL1MX + (6.D0/C1MX-5.D0-X)*CLX
!      +      + (C1PX-2.D0*CPFFX)*CLX2      
! C--   P_QQ(x) in FP eq. (4.50) 
!       PQQ  = C16S9*AAA + 4.*BBB + C2S3*NF*CCC
! C--   0.5*P_A(x) in FP eq. (4.55) and P_QQBAR(x) eq. (4.51)      
!       PQQB = - C4S9 * ( CPFFMX*CS2X + C1PX*CLX + 2.*C1MX )
!  
!       PM1TFUNC = PQQ - PQQB + C16S9*DDD
      PM1TFUNC = 0
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION GG1SFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (11)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'qconst.inc'

      CX2    = X**2
      C1PX   = 1.+X
      C1MX   = 1.-X
      CLX    = LOG(X)
      CLX2   = CLX**2
      CL1MX  = LOG(C1MX)
      CL1PX  = LOG(C1PX)
      CL1MX2 = CL1MX**2
      CS3X   = -DMB_DILOG(-X)
      CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)

      CPGG   = 1./C1MX + 1./X -2. + X - CX2
      CMPGG  = 1./C1PX - 1./X -2. - X - CX2

      AAA   = -16.+ 8.*X+ C20S3*CX2 + C4S3/X + (-6.-10.*X)*CLX +
     +        (-2.)*C1PX*CLX2
      BBB   = 2.* C1MX +  26./9.*(CX2-1./X) - C4S3*C1PX*CLX -
     +        20./9.*CPGG
      CCC   = 27./2.*C1MX + 67./9.*(CX2-1./X)+(-25./3.+11./3.*x-
     +        44./3.*CX2)*CLX+4.*C1PX*CLX2+(67./9.-4.*CLX*CL1MX +
     +        CLX2-CPI2S3)*CPGG + 2.*CMPGG*CS2X

      GG1SFUNC = C2S3*NF*AAA + 1.5*NF*BBB + 9.* CCC

      RETURN
      END
      

C     ========================================
      DOUBLE PRECISION FUNCTION GG1TFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (12)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'qconst.inc'
      
      CX2    = X**2
      C1PX   = 1.+X
      C1MX   = 1.-X
      CLX    = LOG(X)
      CLX2   = CLX**2
      CL1MX  = LOG(C1MX)
      CL1PX  = LOG(C1PX)
      CL1MX2 = CL1MX**2
      CPGFX  = CX2 + C1MX**2
      CPGFMX = CX2 + C1PX**2
      CS1X   = -DMB_DILOG(1.D0-X)
      CS3X   = -DMB_DILOG(-X)
      CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)

      CPGG  = 1./C1MX + 1./X -2. + X - CX2
      CMPGG = 1./C1PX - 1./X -2. - X - CX2

      AAA     = -4.+12.*x-164./9.*CX2+92./9./X+(10.+14.*X+C16S3*CX2+
     +          C16S3/X)*CLX + 2.*C1PX*CLX2
      BBB     = 2.-2.*X+26./9.*(CX2-1./X)-C4S3*C1PX*CLX-
     +          (20./9.+8./3.*CLX)*CPGG
      CCC     = 27./2.*(C1MX)+67./9.*(CX2-1./X)+(11./3.-25./3.*X-
     +          44./3./X)*CLX -4.*(C1PX) * CLX2 + (4.*CLX*CL1MX -
     +          3.*CLX2+22./3.*CLX-CPI2S3+67./9.)*CPGG+
     +          2.*CMPGG*CS2X
      GG1TFUNC = 2./3.*NF*AAA+3./2.*NF*BBB+9.*CCC
      
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION GF1SFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (11)

*     Warning: 'Furmanski-Petronzio' notation where
*                          
*              |F'| = |FF GF| |F| 
*              |G'| = |FG GG| |G|
*                  
*     In ../src/qcdwfun.f we will swap to math notation
*
*              |F'| = |FF FG| |F| 
*              |G'| = |GF GG| |G|                   

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!       include 'qconst.inc'
! 
!       CX2    = X**2
!       C1PX   = 1.+X
!       C1MX   = 1.-X
!       CLX    = LOG(X)
!       CLX2   = CLX**2
!       CL1MX  = LOG(C1MX)
!       CL1PX  = LOG(C1PX)
!       CL1MX2 = CL1MX**2
!       CPGFX  = CX2 + C1MX**2
!       CPGFMX = CX2 + C1PX**2
!       CS3X   = -DMB_DILOG(-X)
!       CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)
! 
!       AAA =   4. - 9.*X + (4.*X-1.)*CLX + (2.*X-1.)*CLX2
!      +      + 4.*CL1MX
!      +      + (2.*CLX-2.*CLX*CL1MX+CLX2-2.*CL1MX+CL1MX2+CPIE)
!      +      * 2. * CPGFX
!       DDD =   C182S9 + C14S9*X + C40S9/X + (C136S3*X-C38S3)*CLX
!      +      - 4.*CL1MX - (2.+8.*X)*CLX2 + 2.*CS2X*CPGFMX
!      +      + (C44S3*CLX-CLX2-2.*CL1MX2+4.*CL1MX+CPIF) * CPGFX
! 
!       GF1SFUNC = C2S3*NF*AAA + 1.5*NF*DDD
	GF1SFUNC=0
      RETURN
      END
      
C     ========================================
      DOUBLE PRECISION FUNCTION GF1TFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (12)

*     Warning: 'Furmanski-Petronzio' notation where
*                          
*              |F'| = |FF GF| |F| 
*              |G'| = |FG GG| |G|
*                  
*     In ../src/qcdwfun.f we will swap to math notation
*
*              |F'| = |FF FG| |F| 
*              |G'| = |GF GG| |G|

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!       include 'qconst.inc'
!       
!       CX2    = X**2
!       C1PX   = 1.+X
!       C1MX   = 1.-X
!       CLX    = LOG(X)
!       CLX2   = CLX**2
!       CL1MX  = LOG(C1MX)
!       CL1PX  = LOG(C1PX)
!       CL1MX2 = CL1MX**2
!       CPGFX  = CX2 + C1MX**2
!       CPGFMX = CX2 + C1PX**2
!       CS1X   = -DMB_DILOG(1.D0-X)
!       CS3X   = -DMB_DILOG(-X)
!       CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)
! 
!       CPGG   = 1./C1MX + 1./X -2. + X - CX2
!       CMPGG  = 1./C1PX - 1./X -2. - X - CX2
! 
!       AAA     = -8./3.-(16./9.+8./3.*CLX+8./3.*CL1MX)*CPGFX
!       BBB     = -2.+3.*X+(-7.+8.*X)*CLX-4.*CL1MX + (1.-2.*X)*CLX2
!      +          +(-4.*CLX*CL1MX-2.*CLX2-2.*CL1MX+2.*CLX-2.*CL1MX2
!      +          +16.*CS1X+ 2.*PI*PI - 10.)*CPGFX
!       CCC     = -152./9.+166./9.*X-40./9./X+ (-C4S3-76./3.*X)*CLX+
!      +          4.*CL1MX + (2.+8.*X)*CLX2+ (8.*CLX*CL1MX-CLX2-
!      +          C4S3*CLX+10./3.*CL1MX+2.*CL1MX2-16.*CS1X-7.*CPI2S3+
!      +          178./9.)*CPGFX+2.*CPGFMX*CS2X
!       GF1TFUNC = (0.5*NF)**2*AAA+2./3.*NF*BBB+3./2.*NF*CCC
      GF1TFUNC = 0
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION FG1SFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (11)

*     Warning: 'Furmanski-Petronzio' notation where
*                          
*              |F'| = |FF GF| |F| 
*              |G'| = |FG GG| |G|
*                  
*     In ../src/qcdwfun.f we will swap to math notation
*
*              |F'| = |FF FG| |F| 
*              |G'| = |GF GG| |G|

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!       include 'qconst.inc'
! 
!       CX2    = X**2
!       C1PX   = 1.+X
!       C1MX   = 1.-X
!       CLX    = LOG(X)
!       CLX2   = CLX**2
!       CL1MX  = LOG(C1MX)
!       CL1PX  = LOG(C1PX)
!       CL1MX2 = CL1MX**2
!       CPFGX  = (1.+C1MX**2) / X
!       CPFGMX = - (1.+C1PX**2) / X
!       CS3X   = -DMB_DILOG(-X)
!       CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)
! 
!       AAA   = -5./2.- 7./2.*X+(2.+7./2.*X)*CLX+(-1.+0.5*X)*CLX2 
!      +        -2.*X*CL1MX+ (-3.*CL1MX-CL1MX2)*CPFGX
!       BBB   = 28./9.+65./18.*X+44./9.*CX2+(-12.-5.*X-8./3.*CX2)*CLX+
!      +        (4.+X)*CLX2+2.*X*CL1MX+ (-2.*CLX*CL1MX+0.5*CLX2+
!      +        11./3.*CL1MX+CL1MX2-0.5*CPI2S3+0.5)*CPFGX+CPFGMX*CS2X
!       CCC   = -C4S3*X- (20./9.+C4S3*CL1MX)*CPFGX
! 
!       FG1SFUNC = C16S9*AAA+4.*BBB+2./3.*NF*CCC
      FG1SFUNC = 0
      RETURN
      END
      
C     ========================================
      DOUBLE PRECISION FUNCTION FG1TFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (12)

*     Warning: 'Furmanski-Petronzio' notation where
*                          
*              |F'| = |FF GF| |F| 
*              |G'| = |FG GG| |G|
*                  
*     In ../src/qcdwfun.f we will swap to math notation
*
*              |F'| = |FF FG| |F| 
*              |G'| = |GF GG| |G|

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!       include 'qconst.inc'
!       
!       idum   = nf      !avoid compiler warning
!       
!       CX2    = X**2
!       C1PX   = 1.+X
!       C1MX   = 1.-X
!       CLX    = LOG(X)
!       CLX2   = CLX**2
!       CL1MX  = LOG(C1MX)
!       CL1PX  = LOG(C1PX)
!       CL1MX2 = CL1MX**2
!       CPFFX  = (1.+CX2) / C1MX
!       CPFFMX = (1.+CX2) / C1PX
!       CPFGX  = (1.+C1MX**2) / X
!       CPFGMX = - (1.+C1PX**2) / X
!       CS1X   = -DMB_DILOG(1.D0-X)
!       CS3X   = -DMB_DILOG(-X)
!       CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)
! 
!       AAA =   -.5 + 4.5*X + (-8.+.5*X)*CLX + 2.*X*CL1MX
!      +        + (1.-.5*X)*CLX2
!      +        + (CL1MX2+4.*CLX*CL1MX-8.*CS1X-CPIB) * CPFGX
!       BBB =     C62S9 - C35S18*X - C44S9*CX2
!      +        + (2.+12.*X+C8S3*CX2) * CLX
!      +        - 2.*X*CL1MX - (4.+X)*CLX2 + CPFGMX*CS2X
!      +        + ( - 2.*CLX*CL1MX - 3.*CLX - 1.5*CLX2
!      +        - CL1MX2 + 8.*CS1X + CPIC ) * CPFGX
!       FG1TFUNC = C16S9*AAA + 4.*BBB

      FG1TFUNC = 0
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION FF1SFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (11)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!       include 'qconst.inc'
! 
!       CX2    = X**2
!       C1PX   = 1.+X
!       C1MX   = 1.-X
!       CLX    = LOG(X)
!       CLX2   = CLX**2
!       CL1MX  = LOG(C1MX)
!       CL1PX  = LOG(C1PX)
!       CPFFX  = (1.+CX2) / C1MX
!       CPFFMX = (1.+CX2) / C1PX
!       CS3X   = -DMB_DILOG(-X)
!       CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)
! 
!       AAA = - CPFFX*CLX*(1.5+2.*CL1MX) + 2.*CPFFMX*CS2X
!      +      - 1. + X + (.5-1.5*X)*CLX - .5*C1PX*CLX2
!       BBB =   CPFFX*(C11S6*CLX+.5*CLX2+CPIA) - CPFFMX*CS2X
!      +      + C14S3*C1MX
!       CCC = - CPFFX*(C10S9+C2S3*CLX) + C40S9/X - 2.*C1PX*CLX2
!      +      - C16S3 + C40S3*X + (10.*X+C16S3*CX2+2.)*CLX
!      +      - C112S9*CX2
! 
!       FF1SFUNC = C16S9*AAA + 4.*BBB + C2S3*NF*CCC
      FF1SFUNC = 0
      RETURN
      END
      
C     ========================================
      DOUBLE PRECISION FUNCTION FF1TFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (12)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!       include 'qconst.inc'
!       
!       CX2    = X**2
!       C1PX   = 1.+X
!       C1MX   = 1.-X
!       CLX    = LOG(X)
!       CLX2   = CLX**2
!       CL1MX  = LOG(C1MX)
!       CL1PX  = LOG(C1PX)
!       CL1MX2 = CL1MX**2
!       CPFFX  = (1.+CX2) / C1MX
!       CPFFMX = (1.+CX2) / C1PX
!       CPFGX  = (1.+C1MX**2) / X
!       CPFGMX = - (1.+C1PX**2) / X
!       CS1X   = -DMB_DILOG(1.D0-X)
!       CS3X   = -DMB_DILOG(-X)
!       CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)
! 
!       AAA =     CPFFX*(1.5*CLX-2.*CLX2+2.*CLX*CL1MX) + 2.*CPFFMX*CS2X
!      +        - 1. + X + (-1.5+.5*X)*CLX + .5*C1PX*CLX2
!       BBB =     CPFFX*(C11S6*CLX+.5*CLX2+CPIA) - CPFFMX*CS2X
!      +        + C14S3*C1MX
!       CCC =   - CPFFX*(C2S3*CLX+C10S9)
!      +        - C52S3 + C28S3*X + C112S9*CX2 - C40S9/X
!      +        - (10.+18.*X+C16S3*CX2)*CLX + 2.*C1PX*CLX2
!       FF1TFUNC = C16S9*AAA + 4.*BBB + C2S3*NF*CCC
      FF1TFUNC = 0
      RETURN
      END            
