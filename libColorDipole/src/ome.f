*-MB----------------------------------------------------------------------
*-MB  Coefficient functions from Buza et al., (BMSN, hep-ph/9612398). The
*-MB  code is passed to me by Andreas Vogt but originally comes from 
*-MB  a program written by Buza and van Neerven. Not all functions in this
*-MB  file are used by QCDNUM
*-MB---------------------------------------------------------------------- 

*-MB----------------------------------------------------------------------
*-MB  Not used

C OPERATOR MATRIX ELEMENT A^(NS,1)_QQ

      REAL*8 FUNCTION A1QQNS(Z,FS2,HM2,OR)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CHOICE/NOP,NL,MODE
      L=IDINT(OR)
      DLZ=DLOG(Z)
      DLMZ=DLOG(1.0D0-Z)
      DL=DLOG(HM2/FS2)
C SINGLE LOGARITHMIC TERMS LN^2(M^2/MU^2)
C C_F - PART
      A1=2.0D0*(1.0D0+Z)
      A0=0.0D0
      IF(L.LE.1) GO TO 1
C NON LOGARITHMIC TERMS
C C_F - PART
      A0=2.0D0*(1.0D0+Z)*(2.0D0*DLMZ+1.0D0)
    1 A1QQNS=4.0D0*(A1*DL+A0)/3.0D0
      RETURN
      END

*-MB----------------------------------------------------------------------
*-MB  Not used

C SOFT GLUON PART OF A^(NS,1)_QQ

      REAL*8 FUNCTION SOFTQ1(Z,FS2,HM2,OR)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CHOICE/NOP,NL,MODE
      L=IDINT(OR)
      DLMZ=DLOG(1.0D0-Z)
      DM=1.0D0/(1.0D0-Z)
      DL=DLOG(HM2/FS2)
C SINGLE LOGARITHMIC TERMS LN^2(M^2/MU^2)
C C_F - PART
      A1=-4.0D0
      A0=0.0D0
      IF(L.LE.1) GO TO 1
C NON LOGARITHMIC TERMS
C C_F - PART
      A0=-8.0D0*DLMZ-4.0D0
    1 SOFTQ1=4.0D0*(A1*DL+A0)*DM/3.0D0
      RETURN
      END

*-MB----------------------------------------------------------------------
*-MB  Not used

C SOFT PLUS VIRTUAL GLUON PART OF A^(NS,1)_QQ

      REAL*8 FUNCTION CORQ1(Z,FS2,HM2,OR)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CHOICE/NOP,NL,MODE
      L=IDINT(OR)
      ZETA2=1.6449340668D0
      DLMZ=DLOG(1.0D0-Z)
      DL=DLOG(HM2/FS2)
C SINGLE LOGARITHMIC TERMS LN^2(M^2/MU^2)
C C_F - PART
      A1=-4.0D0*DLMZ-3.0D0
      A0=0.0D0
      IF(L.LE.1) GO TO 1
C NON LOGARITHMIC TERMS
C C_F - PART
      A0=-4.0D0*DLMZ*DLMZ-4.0D0*DLMZ+4.0D0
    1 CORQ1=4.0D0*(A1*DL+A0)/3.0D0
      RETURN
      END

*-MB----------------------------------------------------------------------
*-MB  Not used

C OPERATOR MATRIX ELEMENT A^(1)_Qg

      REAL*8 FUNCTION A1QG(Z,FS2,HM2)
      IMPLICIT REAL*8(A-H,O-Z)
      DL=DLOG(HM2/FS2)
      A1QG=-2.0D0*(1.0D0-2.0D0*Z+2.0D0*Z*Z)*DL
      RETURN
      END

*-MB----------------------------------------------------------------------
*-MB  \~{A}^{S,(2)}_{Hg} eq. (B.3) in BMSN. QCDNUM does not use this 
*-MB  function but, instead, the parameterization by Andreas Vogt which can
*-MB  be found in the file xa2hgp.f

C OPERATOR MATRIX ELEMENT A^(2)_Qg

      REAL*8 FUNCTION A2QG(Z,FS2,HM2,OR)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CHOICE/NOP,NL,MODE
      COMPLEX*16 WGPLG
      L=IDINT(OR)
*-MB  Convert to dble to avoid compiler warnings (01-01-2012)       
      S121MZ=dble(WGPLG(1,2,1.0D0-Z))
      S12MZ=dble(WGPLG(1,2,-Z))
      S211MZ=dble(WGPLG(2,1,1.0D0-Z))
      S21MZ=dble(WGPLG(2,1,-Z))
      S111MZ=dble(WGPLG(1,1,1.0D0-Z))
      S11MZ=dble(WGPLG(1,1,-Z))
      ZETA2=1.6449340668D0
      ZETA3=1.2020569031D0
      DLZ=DLOG(Z)
      DLZ2=DLZ*DLZ
      DLZ3=DLZ2*DLZ
      DLM=DLOG(1.0D0-Z)
      DLM2=DLM*DLM
      DLM3=DLM2*DLM
      DLP=DLOG(1.0D0+Z)
      DLP2=DLP*DLP
      DL=DLOG(HM2/FS2)
      DL2=DL*DL
C DOUBLE LOGARITHMIC TERMS LN^2(M^2/MU^2)
C C_F.T_f - PART
      A2=(8.0D0-16.0D0*Z+16.0D0*Z*Z)*DLM-(4.0D0-8.0D0*Z
     1+16.0D0*Z*Z)*DLZ-2.0D0+8.0D0*Z
C C_A.T_f - PART
      B2=-(8.0D0-16.0D0*Z+16.0D0*Z*Z)*DLM-(8.0D0+32.0D0*Z)*DLZ
     1-16.0D0/Z/3.0D0-4.0D0-32.0D0*Z+124.0D0*Z*Z/3.0D0
C T_f^2 - PART
C      C2=-16.0D0*(1.0D0-2.0D0*Z+2.0D0*Z*Z)/3.0D0
      C2=0.0D0
      A1=0.0D0
      A0=0.0D0
      B1=0.0D0
      B0=0.0D0
      IF(L.LE.1) GO TO 1
C SINGLE LOGARITHMIC TERMS LN(M^2/MU^2)
C C_F.T_f - PART
      A1=(8.0D0-16.0D0*Z+16.0D0*Z*Z)*(2.0D0*DLZ*DLM-DLM2+2.0D0
     1*ZETA2)-(4.0D0-8.0D0*Z+16.0D0*Z*Z)*DLZ2-32.0D0*Z*(1.0D0-Z)
     2*DLM-(12.0D0-16.0D0*Z+32.0D0*Z*Z)*DLZ-56.0D0+116.0D0*Z
     3-80.0D0*Z*Z
C C_A.T_f - PART
      B1=(16.0D0+32.0D0*Z+32.0D0*Z*Z)*(S11MZ+DLZ*DLP)+(8.0D0-16.0D0
     1*Z+16.0D0*Z*Z)*DLM2+(8.0D0+16.0D0*Z)*DLZ2+32.0D0*Z*ZETA2
     2+32.0D0*Z*(1.0D0-Z)*DLM-(8.0D0+64.0D0*Z+352.0D0*Z*Z/3.0D0)
     3*DLZ-160.0D0/Z/9.0D0+16.0D0-200.0D0*Z+1744.0D0*Z*Z/9.0D0
      A0=0.0D0
      B0=0.0D0
      IF(L.LE.2) GO TO 1
C NON LOGARITHMIC TERMS
C C_F.T_f - PART
      A01=(1.0D0-2.0D0*Z+2.0D0*Z*Z)*(8.0D0*ZETA3+4.0D0*DLM3/3.0D0
     1-8.0D0*DLM*S111MZ+8.0D0*ZETA2*DLZ-4.0D0*DLZ*DLM2+2.0D0*DLZ3
     2/3.0D0-8.0D0*DLZ*S111MZ+8.0D0*S211MZ-24.0D0*S121MZ)
      A02=-(4.0D0+96.0D0*Z-64.0D0*Z*Z)*S111MZ-(4.0D0-48.0D0*Z
     1+40.0D0*Z*Z)*ZETA2-(8.0D0+48.0D0*Z-24.0D0*Z*Z)*DLZ*DLM
     2+(4.0D0+8.0D0*Z-12.0D0*Z*Z)*DLM2-(1.0D0+12.0D0*Z-20.0D0*Z*Z)
     3*DLZ2-(52.0D0*Z-48.0D0*Z*Z)*DLM-(16.0D0+18.0D0*Z+48.0D0*Z*Z)
     4*DLZ+26.0D0-82.0D0*Z+80.0D0*Z*Z+Z*Z*(-16.0D0*ZETA2*DLZ
     5+4.0D0*DLZ3/3.0D0+ 16.0D0*DLZ*S111MZ+ 32.0D0*S121MZ)
      A0=A01+A02
C C_A.T_f - PART
      B01=(1.0D0-2.0D0*Z+2.0D0*Z*Z)*(-4.0D0*DLM3/3.0D0+8.0D0*DLM
     1*S111MZ-8.0D0*S211MZ)+(1.0D0+2.0D0*Z+2.0D0*Z*Z)*(-8.0D0*ZETA2
     2*DLP-16.0D0*DLP*S11MZ-8.0D0*DLZ*DLP2+4.0D0*DLZ2*DLP+8.0D0*DLZ
     3*S11MZ-8.0D0*S21MZ-16.0D0*S12MZ)+(16.0D0+64.0D0*Z)*(2.0D0*S121MZ
     4+DLZ*S111MZ)-(4.0D0+8.0D0*Z)*DLZ3/3.0D0+(8.0D0-32.0D0*Z
     5+16.0D0*Z*Z)*ZETA3-(16.0D0+64.0D0*Z)*ZETA2*DLZ
      B02=(16.0D0*Z+16.0D0*Z*Z)*(S11MZ+DLZ*DLP)+(32.0D0/Z/3.0D0+12.0D0
     1+64.0D0*Z-272.0D0*Z*Z/3.0D0)*S111MZ-(12.0D0+48.0D0*Z
     2-260.0D0*Z*Z/3.0D0+32.0D0/Z/3.0D0)*ZETA2-4.0D0*Z*Z*DLZ*DLM
     3-(2.0D0+8.0D0*Z-10.0D0*Z*Z)*DLM2+(2.0D0+8.0D0*Z+46.0D0*Z*Z/3.0D0)
     4*DLZ2+(4.0D0+16.0D0*Z-16.0D0*Z*Z)*DLM-(56.0D0/3.0D0+172.0D0*Z
     5/3.0D0+1600.0D0*Z*Z/9.0D0)*DLZ-448.0D0/Z/27.0D0-4.0D0/3.0D0
     6-628.0D0*Z/3.0D0+6352.0D0*Z*Z/27.0D0
      B0=B01+B02
    1 A2QG=(2.0D0*A2/3.0D0+1.5D0*B2+0.25D0*C2)*DL2+(2.0D0*A1/3.0D0
     1+1.5D0*B1)*DL+(2.0D0*A0/3.0D0+1.5D0*B0)
      RETURN
      END

*-MB----------------------------------------------------------------------
*-MB  \~{A}^{PS,(2)}_{Hq} eq. (B.1) in BMSN

C OPERATOR MATRIX ELEMENT A^(PS,2)_Qq

      REAL*8 FUNCTION A2QQ(Z,FS2,HM2,OR)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CHOICE/NOP,NL,MODE
      COMPLEX*16 WGPLG
      L=IDINT(OR)
*-MB  Convert to dble to avoid compiler warnings (01-01-2012)      
      S121MZ=dble(WGPLG(1,2,1.0D0-Z))
      S111MZ=dble(WGPLG(1,1,1.0D0-Z))
      ZETA2=1.6449340668D0
      DLZ=DLOG(Z)
      DLZ2=DLZ*DLZ
      DLZ3=DLZ2*DLZ
      DL=DLOG(HM2/FS2)
      DL2=DL*DL
C DOUBLE LOGARITHMIC TERMS LN^2(M^2/MU^2)
C C_F.T_f - PART
      A2=-8.0D0*(1.0D0+Z)*DLZ-16.0D0/Z/3.0D0-4.0D0+4.0D0*Z
     1+16.0D0*Z*Z/3.0D0
      A1=0.0D0
      A0=0.0D0
      IF(L.LE.1) GO TO 1
C SINGLE LOGARITHMIC TERMS LN(M^2/MU^2)
C C_F.T_f - PART
      A1=8.0D0*(1.0D0+Z)*DLZ2-(8.0D0+40.0D0*Z+64.0D0*Z*Z/3.0D0)
     1*DLZ-160.0D0/Z/9.0D0+16.0D0-48.0D0*Z+448.0D0*Z*Z/9.0D0
      A0=0.0D0
      IF(L.LE.2) GO TO 1
C NON LOGARITHMIC TERMS
C C_F.T_f - PART
      A0=(1.0D0+Z)*(32.0D0*S121MZ+16.0D0*DLZ*S111MZ-16.0D0*ZETA2
     1*DLZ-4.0D0*DLZ3/3.0D0)+(32.0D0/Z/3.0D0+8.0D0-8.0D0*Z-32.0D0
     2*Z*Z/3.0D0)*S111MZ+(-32.0D0/Z/3.0D0-8.0D0+8.0D0*Z
     3+32.0D0*Z*Z/3.0D0)*ZETA2+(2.0D0+10.0D0*Z+16.0D0*Z*Z/3.0D0)
     4*DLZ2-(56.0D0/3.0D0+88.0D0*Z/3.0D0+448.0D0*Z*Z/9.0D0)*DLZ
     5-448.0D0/Z/27.0D0-4.0D0/3.0D0-124.0D0*Z/3.0D0+1600.0D0*Z*Z
     6/27.0D0
    1 A2QQ=2.0D0*(A2*DL2+A1*DL+A0)/3.0D0
c     A2QQ= 2./3.D0* ( 
c    1      (1.0D0+Z)*(32.0D0*S121MZ+16.0D0*DLZ*S111MZ 
c    2      - 16.0D0*ZETA2*DLZ-4.0D0*DLZ3/3.0D0)
c    2      + (32.0D0/Z/3.0D0+8.0D0-8.0D0*Z-32.0D0
c    3         *Z*Z/3.0D0)*S111MZ 
c    5 -448.0D0/Z/27.0D0-4.0D0/3.0D0-124.0D0*Z/3.0D0+1600.0D0*Z*Z
c    6 /27.0D0 ) ! + 
c    7 (-32.0D0/Z/3.0D0-8.0D0+8.0D0*Z+32.0D0*Z*Z/3.0D0)*ZETA2) 
      RETURN
      END

*-MB----------------------------------------------------------------------
*-MB  Regular piece of \~{A}^{NS,(2)}_{qq,H} eq. (B.4) in BMSN (A-piece)

C OPERATOR MATRIX ELEMENT A^(NS,2)_qq,Q

      REAL*8 FUNCTION A2QQNS(Z,FS2,HM2,OR)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CHOICE/NOP,NL,MODE
      L=IDINT(OR)
      DLZ=DLOG(Z)
      DLZ2=DLZ*DLZ
      DL=DLOG(HM2/FS2)
      DL2=DL*DL
C DOUBLE LOGARITHMIC TERMS LN^2(M^2/MU^2)
C C_F.T_f - PART
      A2=-4.0D0*(1.0D0+Z)/3.0D0
      A1=0.0D0
      A0=0.0D0
      IF(L.LE.1) GO TO 1
C SINGLE LOGARITHMIC TERMS LN(M^2/MU^2)
C C_F.T_f - PART
      A1=8.0D0*(1.0D0+Z*Z)*DLZ/(1.0D0-Z)/3.0D0+8.0D0/9.0D0
     1-88.0D0*Z/9.0D0
      A0=0.0D0
      IF(L.LE.2) GO TO 1
C NON LOGARITHMIC TERMS
C C_F.T_f - PART
      A0=(1.0D0+Z*Z)*(2.0D0*DLZ2/3.0D0+20.0D0*DLZ/9.0D0)
     1/(1.0D0-Z)+8.0D0*(1.0D0-Z)*DLZ/3.0D0+44.0D0/27.0D0
     2-268.0D0*Z/27.0D0
    1 A2QQNS=2.0D0*(A2*DL2+A1*DL+A0)/3.0D0
      RETURN
      END

*-MB----------------------------------------------------------------------
*-MB  Singular piece of \~{A}^{NS,(2)}_{qq,H} eq. (B.4) in BMSN (B-piece)

C SOFT GLUON PART OF A^(NS,2)_qq,Q

      REAL*8 FUNCTION SOFTQ2(Z,FS2,HM2,OR)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CHOICE/NOP,NL,MODE
      L=IDINT(OR)
      DM=1.0D0/(1.0D0-Z)
      DL=DLOG(HM2/FS2)
      DL2=DL*DL
C DOUBLE LOGARITHMIC TERMS LN^2(M^2/MU^2)
C C_F.T_f - PART
      A2=8.0D0/3.0D0
      A1=0.0D0
      A0=0.0D0
      IF(L.LE.1) GO TO 1
C SINGLE LOGARITHMIC TERMS LN(M^2/MU^2)
C C_F.T_f - PART
      A1=80.0D0/9.0D0
      A0=0.0D0
      IF(L.LE.2) GO TO 1
C NON LOGARITHMIC TERMS
C C_F.T_f - PART
      A0=224.0D0/27.0D0
    1 SOFTQ2=2.0D0*(A2*DL2+A1*DL+A0)*DM/3.0D0
      RETURN
      END

*-MB----------------------------------------------------------------------
*-MB  Local piece of \~{A}^{NS,(2)}_{qq,H} eq. (B.4) in BMSN (C-piece)

C SOFT PLUS VIRTUAL GLUON PART OF A^(NS,2)_qq,Q

      REAL*8 FUNCTION CORQ2(Z,FS2,HM2,OR)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CHOICE/NOP,NL,MODE
      L=IDINT(OR)
      ZETA3=1.2020569031D0
      ZETA2=1.6449340668D0
      DLM=DLOG(1.0D0-Z)
      DL=DLOG(HM2/FS2)
      DL2=DL*DL
C DOUBLE LOGARITHMIC TERMS LN^2(M^2/MU^2)
C C_F.T_f - PART
      A2=8.0D0*DLM/3.0D0+2.0D0
      A1=0.0D0
      A0=0.0D0
      IF(L.LE.1) GO TO 1
C SINGLE LOGARITHMIC TERMS LN(M^2/MU^2)
C C_F.T_f - PART
      A1=80.0D0*DLM/9.0D0+16.0D0*ZETA2/3.0D0+2.0D0/3.0D0
      A0=0.0D0
      IF(L.LE.2) GO TO 1
C NON LOGARITHMIC TERMS
C C_F.T_f - PART
      A0=224.0D0*DLM/27.0D0-8.0D0*ZETA3/3.0D0+40.0D0*ZETA2/9.0D0
     1+73.0D0/18.0D0
    1 CORQ2=2.0D0*(A2*DL2+A1*DL+A0)/3.0D0
      RETURN
      END

*-MB----------------------------------------------------------------------
*-MB  Regular piece of \~{A}^{S,(2)}_{gg,H} eq. (B.7) in BMSN (A-piece)

C OPERATOR MATRIX ELEMENT A^(2)_gg,Q

      REAL*8 FUNCTION A2GG(Z,FS2,HM2,OR)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CHOICE/NOP,NL,MODE
      L=IDINT(OR)
      DLZ=DLOG(Z)
      DLZ2=DLZ*DLZ
      DLZ3=DLZ2*DLZ
      DLM=DLOG(1.0D0-Z)
      DL=DLOG(HM2/FS2)
      DL2=DL*DL
C DOUBLE LOGARITHMIC TERMS LN^2(M^2/MU^2)
C C_F.T_f - PART
      A2=8.0D0*(1.0D0+Z)*DLZ+16.0D0/Z/3.0D0+4.0D0-4.0D0*Z
     1-16.0D0*Z*Z/3.0D0
C C_A.T_f - PART
      B2=8.0D0/Z/3.0D0-16.0D0/3.0D0+8.0D0*Z/3.0D0-8.0D0*Z*Z/3.0D0
      A1=0.0D0
      A0=0.0D0
      B1=0.0D0
      B0=0.0D0
      IF(L.LE.1) GO TO 1
C SINGLE LOGARITHMIC TERMS LN(M^2/MU^2)
C C_F.T_f - PART
      A1=8.0D0*(1.0D0+Z)*DLZ2+(24.0D0+40.0D0*Z)*DLZ-16.0D0/Z/3.0D0
     1+64.0D0-32.0D0*Z-80.0D0*Z*Z/3.0D0
C C_A.T_f - PART
      B1=16.0D0*(1.0D0+Z)*DLZ/3.0D0+184.0D0/Z/9.0D0-232.0D0/9.0D0
     1+152.0D0*Z/9.0D0-184.0D0*Z*Z/9.0D0
      A0=0.0D0
      B0=0.0D0
      IF(L.LE.2) GO TO 1
C NON LOGARITHMIC TERMS
C C_F.T_f - PART
      A0=4.0D0*(1.0D0+Z)*DLZ3/3.0D0+(6.0D0+10.0D0*Z)*DLZ2
     1+(32.0D0+48.0D0*Z)*DLZ
     2-8.0D0/Z+80.0D0-24.0D0*Z*Z-48.0D0*Z
C C_A.T_f - PART
      B0=4.0D0*(1.0D0+Z)*DLZ2/3.0D0+(52.0D0+88.0D0*Z)*DLZ/9.0D0
     1-4.0D0*Z*DLM/3.0D0
     2+(556.0D0/Z-628.0D0+548.0D0*Z-700.0D0*Z*Z)/27.0D0
    1 A2GG=(2.0D0*A2/3.0D0+1.5D0*B2)*DL2+(2.0D0*A1/3.0D0+1.5D0*B1)
     1*DL+(2.0D0*A0/3.0D0+1.5D0*B0)
      RETURN
      END

*-MB----------------------------------------------------------------------
*-MB  Singular piece of \~{A}^{S,(2)}_{gg,H} eq. (B.7) in BMSN (B-piece)

C SOFT GLUON PART OF A^(2)_gg,Q

      REAL*8 FUNCTION SOFTG(Z,FS2,HM2,OR)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CHOICE/NOP,NL,MODE
      L=IDINT(OR)
      DM=1.0D0/(1.0D0-Z)
      DL=DLOG(HM2/FS2)
      DL2=DL*DL
C DOUBLE LOGARITHMIC TERMS LN^2(M^2/MU^2)
C C_A.T_f - PART
      B2=8.0D0/3.0D0
      B1=0.0D0
      B0=0.0D0
      IF(L.LE.1) GO TO 1
C SINGLE LOGARITHMIC TERMS LN(M^2/MU^2)
C C_A.T_f - PART
      B1=80.0D0/9.0D0
      B0=0.0D0
      IF(L.LE.2) GO TO 1
C NON LOGARITHMIC TERMS
C C_A.T_f - PART
      B0=224.0D0/27.0D0
    1 SOFTG=1.5D0*(B2*DL2+B1*DL+B0)*DM
      RETURN
      END

*-MB----------------------------------------------------------------------
*-MB  Not used

C SOFT PLUS VIRTUAL GLUON PART OF A^(1)_gg,Q

      REAL*8 FUNCTION CORG1(FS2,HM2)
      IMPLICIT REAL*8(A-H,O-Z)
      DL=DLOG(HM2/FS2)
      CORG1=2.0D0*DL/3.0D0
      RETURN
      END

*-MB----------------------------------------------------------------------
*-MB  Local piece of \~{A}^{S,(2)}_{gg,H} eq. (B.7) in BMSN (C-piece)

C SOFT PLUS VIRTUAL GLUON PART OF A^(2)_gg,Q

      REAL*8 FUNCTION CORG2(Z,FS2,HM2,OR)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CHOICE/NOP,NL,MODE
      L=IDINT(OR)
      DLM=DLOG(1.0D0-Z)
      DL=DLOG(HM2/FS2)
      DL2=DL*DL
C DOUBLE LOGARITHMIC TERMS LN^2(M^2/MU^2)
C C_F.T_f - PART
      A2=0.0D0
C C_A.T_f - PART
      B2=8.0D0*DLM/3.0D0
C T_f^2 - PART
      C2=16.0D0/9.0D0
      A1=0.0D0
      A0=0.0D0
      B1=0.0D0
      B0=0.0D0
      IF(L.LE.1) GO TO 1
C SINGLE LOGARITHMIC TERMS LN(M^2/MU^2)
C C_F.T_f - PART
      A1=4.0D0
C C_A.T_f - PART
      B1=80.0D0*DLM/9.0D0+16.0D0/3.0D0
      B0=0.0D0
      A0=0.0D0
      IF(L.LE.2) GO TO 1
C NON LOGARITHMIC TERMS
C C_F.T_f - PART
      A0=-15.0D0
C C_A.T_f - PART
      B0=224.0D0/27.0D0*DLM+10.0D0/9.0D0
    1 CORG2=(2.0D0*A2/3.0D0+1.5D0*B2+0.25D0*C2)*DL2+(2.0D0*A1/3.0D0
     1+1.5D0*B1)*DL+(2.0D0*A0/3.0D0+1.5D0*B0)
      RETURN
      END

*-MB----------------------------------------------------------------------
*-MB  \~{A}^{S,(2)}_{gq,H} eq. (B.5) in BMSN

      REAL*8 FUNCTION A2GQ(Z,FS2,HM2,OR)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CHOICE/NOP,NL,MODE
      L=IDINT(OR)
      DLM=DLOG(1.0D0-Z)
      DLM2=DLM*DLM
      DL=DLOG(HM2/FS2)
      DL2=DL*DL
C DOUBLE LOGARITHMIC TERMS LN^2(M^2/MU^2)
C C_F.T_f - PART
      A2=16.0D0/Z/3.0D0-16.0D0/3.0D0+8.0D0*Z/3.0D0
      A1=0.0D0
      A0=0.0D0
      IF(L.LE.1) GO TO 1
C SINGLE LOGARITHMIC TERMS LN(M^2/MU^2)
C C_F.T_f - PART
      A1=160.0D0/Z/9.0D0-160.0D0/9.0D0+128.0D0*Z/9.0D0
     1+(32.0D0/Z/3.0D0-32.0D0/3.0D0+16.0D0*Z/3.0D0)*DLM
      A0=0.0D0
      IF(L.LE.2) GO TO 1
C NON LOGARITHMIC TERMS
C C_F.T_f - PART
      A0=4.0D0*(2.0D0/Z-2.0D0+Z)*DLM2/3.0D0+8.0D0*(10.0D0/Z
     1-10.0D0+8.0D0*Z)*DLM/9.0D0+(448.0D0/Z-448.0D0+344.0D0*Z)
     2/27.0D0
    1 A2GQ=2.0D0*(A2*DL2+A1*DL+A0)/3.0D0
      RETURN
      END

*-MB----------------------------------------------------------------------
*-MB  \~{A}^{S,(2)}_{gq,H} eq. (B.5) in BMSN

      REAL*8 FUNCTION AQGAQQ(Z,FS2,HM2)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 WGPLG
      ZETA2=1.6449340668D0
      CF=4.0D0/3.0D0
      DL1=DLOG(HM2/FS2)
*-MB  Convert to dble to avoid compiler warnings (01-01-2012)      
      SPZ=dble(WGPLG(1,1,1.0D0-Z))
      DL=DLOG(Z)
      DLM=DLOG(1.0D0-Z)
      A1=4.0D0*(1.0D0-2.0D0*Z+2.0D0*Z*Z)*(DLM-DL)+2.0D0*(1.0D0
     1-2.0D0*Z)*DL-1.0D0+4.0D0*Z
      A2=4.0D0*(1.0D0-2.0D0*Z+4.0D0*Z*Z)*(SPZ+DL*DLM)-4.0D0
     1*(1.0D0-2.0D0*Z+2.0D0*Z*Z)*DLM*DLM+4.0D0*(1.0D0-3.0D0*Z
     2+Z*Z)*DLM+2.0D0*(1.0D0+2.0D0*Z-2.0D0*Z*Z)*DL
     3+2.0D0*(2.0D0-5.0D0*Z+5.0D0*Z*Z)
      AQGAQQ=-2.0D0*CF*DL1*(-A1*DL1+A2)
      RETURN
      END

