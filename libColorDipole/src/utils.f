
C--   This is the file utils.f containing MBUTIL utility routines
C--
C--   subroutine smb_vfill(aa,n,val)
C--   subroutine smb_ifill(ia,n,ival)
C--   subroutine smb_rsort(a,n)
C--   subroutine smb_asort(a,n,m)
C--   real function rmb_urand(iy)
C--   double precision function dmb_dilog(x)
C--   double precision function dmb_gamma(x)
C--   subroutine smb_deriv(f,x,delta,dfdx,rerr)
C--   double precision function dmb_gauss(f,a,b,eps)
C--   subroutine smb_dminv(n,a,idim,ir,ifail)
C--   subroutine smb_dfact(n,a,idim,ir,ifail,det,jfail)
C--   subroutine smb_dfinv(n,a,idim,ir)
C--   subroutine smb_dsinv(n,a,idim,ifail)
C--   subroutine smb_dseqn(n,a,idim,ifail,k,b)

C     ==============================
      subroutine smb_vfill(aa,n,val)
C     ==============================

C--   Set n elements of AA to val
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      implicit double precision (a-h,o-z)

      dimension aa(n)

      do i = 1,n
        aa(i) = val
      enddo

      return
      end

C     ===============================
      subroutine smb_ifill(ia,n,ival)
C     ===============================

C--   Set n elements of IA to ival
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      implicit double precision (a-h,o-z)

      dimension ia(n)

      do i = 1,n
        ia(i) = ival
      enddo

      return
      end

C     =========================
      subroutine smb_rsort(a,n)
C     =========================
 
C--   Cernlib routine FLPSOR M103
C--
C--   Taken from /afs/cern.ch/asis/share/cern/2001/src/...
C--           .../packlib/kernlib/kerngen/tcgen/flpsor.F

C
C   SORT THE ONE-DIMENSIONAL FLOATING POINT ARRAY A(1),...,A(N) BY
C   INCREASING VALUES
C
C-    PROGRAM  M103  TAKEN FROM CERN PROGRAM LIBRARY,  29-APR-78
C
      DIMENSION A(N)
*mb   COMMON /SLATE/ LT(20),RT(20)
      DIMENSION LT(20),RT(20)
      INTEGER R,RT
C
      LEVEL=1
      LT(1)=1
      RT(1)=N
   10 L=LT(LEVEL)
      R=RT(LEVEL)
      LEVEL=LEVEL-1
   20 IF(R.GT.L) GO TO 200
      IF(LEVEL) 50,50,10
C
C   SUBDIVIDE THE INTERVAL L,R
C     L : LOWER LIMIT OF THE INTERVAL (INPUT)
C     R : UPPER LIMIT OF THE INTERVAL (INPUT)
C     J : UPPER LIMIT OF LOWER SUB-INTERVAL (OUTPUT)
C     I : LOWER LIMIT OF UPPER SUB-INTERVAL (OUTPUT)
C
  200 I=L
      J=R
      M=(L+R)/2
      X=A(M)
  220 IF(A(I).GE.X) GO TO 230
      I=I+1
      GO TO 220
  230 IF(A(J).LE.X) GO TO 231
      J=J-1
      GO TO 230
C
  231 IF(I.GT.J) GO TO 232
      W=A(I)
      A(I)=A(J)
      A(J)=W
      I=I+1
      J=J-1
      IF(I.LE.J) GO TO 220
C
  232 LEVEL=LEVEL+1
      IF((R-I).GE.(J-L)) GO TO 30
      LT(LEVEL)=L
      RT(LEVEL)=J
      L=I
      GO TO 20
   30 LT(LEVEL)=I
      RT(LEVEL)=R
      R=J
      GO TO 20
   50 RETURN
      END

C     ===========================
      subroutine smb_asort(a,n,m)
C     ===========================
 
C--   Extension of Cernlib routine FLPSOR M103
C--
C--   Sort real array a into itself but weed out equal elements.
C--   On exit A(1) < A(2) < .... < A(m), with m <= n.
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      dimension a(n)

      call smb_rsort(a,n)

      m    = 1
      do i = 2,n
        if(a(i).ne.a(m)) then
          m    = m+1
          a(m) = a(i)
        endif
      enddo

      return
      end

C     ===========================
      real function rmb_urand(iy)
C     ===========================

C--   Uniform random number generator taken from netlib.

      integer  iy
c
c  urand is a uniform random number generator based  on  theory  and
c  suggestions  given  in  d.e. knuth (1969),  vol  2.   the integer  iy
c  should be initialized to an arbitrary integer prior to the first call
c  to urand.  the calling program should  not  alter  the  value  of  iy
c  between  subsequent calls to urand.  values of urand will be returned
c  in the interval (0,1).
c
      integer  ia,ic,itwo,m2,m,mic
      double precision  halfm
      real  s
      double precision  datan,dsqrt
      data m2/0/,itwo/2/
*mb   just to have no compiler complaints about uninitialized variables
      data ia/0/,ic/0/,mic/0/,s/0./
*mb 
      if (m2 .ne. 0) go to 20
c
c  if first entry, compute machine integer word length
c
      m = 1
   10 m2 = m
      m = itwo*m2
      if (m .gt. m2) go to 10
      halfm = m2
c
c  compute multiplier and increment for linear congruential method
c
      ia = 8*idint(halfm*datan(1.d0)/8.d0) + 5
      ic = 2*idint(halfm*(0.5d0-dsqrt(3.d0)/6.d0)) + 1
      mic = (m2 - ic) + m2
c
c  s is the scale factor for converting to floating point
c
      s = real(0.5D0/halfm)
c
c  compute next random number
c
   20 iy = iy*ia
c
c  the following statement is for computers which do not allow
c  integer overflow on addition
c
      if (iy .gt. mic) iy = (iy - m2) - m2
c
      iy = iy + ic
c
c  the following statement is for computers where the
c  word length for addition is greater than for multiplication
c
      if (iy/2 .gt. m2) iy = (iy - m2) - m2
c
c  the following statement is for computers where integer
c  overflow affects the sign bit
c
      if (iy .lt. 0) iy = (iy + m2) + m2
      rmb_urand = float(iy)*s

      return
      end
        
C     ======================================
      double precision function dmb_dilog(x)
C     ======================================

C--   Cernlib routine DDILOG C332
C--
C--   Taken from /afs/cern.ch/asis/share/cern/2001/src/...
C--           .../mathlib/gen/c/dilog64.F

      implicit double precision (a-h,o-z)

      DIMENSION C(0:19)

      PARAMETER (Z1 = 1, HF = Z1/2)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (PI3 = PI**2/3, PI6 = PI**2/6, PI12 = PI**2/12)

      DATA C( 0) / 0.42996 69356 08136 97D0/
      DATA C( 1) / 0.40975 98753 30771 05D0/
      DATA C( 2) /-0.01858 84366 50145 92D0/
      DATA C( 3) / 0.00145 75108 40622 68D0/
      DATA C( 4) /-0.00014 30418 44423 40D0/
      DATA C( 5) / 0.00001 58841 55418 80D0/
      DATA C( 6) /-0.00000 19078 49593 87D0/
      DATA C( 7) / 0.00000 02419 51808 54D0/
      DATA C( 8) /-0.00000 00319 33412 74D0/
      DATA C( 9) / 0.00000 00043 45450 63D0/
      DATA C(10) /-0.00000 00006 05784 80D0/
      DATA C(11) / 0.00000 00000 86120 98D0/
      DATA C(12) /-0.00000 00000 12443 32D0/
      DATA C(13) / 0.00000 00000 01822 56D0/
      DATA C(14) /-0.00000 00000 00270 07D0/
      DATA C(15) / 0.00000 00000 00040 42D0/
      DATA C(16) /-0.00000 00000 00006 10D0/
      DATA C(17) / 0.00000 00000 00000 93D0/
      DATA C(18) /-0.00000 00000 00000 14D0/
      DATA C(19) /+0.00000 00000 00000 02D0/

      IF(X .EQ. 1) THEN
       H=PI6
      ELSEIF(X .EQ. -1) THEN
       H=-PI12
      ELSE
       T=-X
       IF(T .LE. -2) THEN
        Y=-1/(1+T)
        S=1
        A=-PI3+HF*(LOG(-T)**2-LOG(1+1/T)**2)
       ELSEIF(T .LT. -1) THEN
        Y=-1-T
        S=-1
        A=LOG(-T)
        A=-PI6+A*(A+LOG(1+1/T))
       ELSE IF(T .LE. -HF) THEN
        Y=-(1+T)/T
        S=1
        A=LOG(-T)
        A=-PI6+A*(-HF*A+LOG(1+T))
       ELSE IF(T .LT. 0) THEN
        Y=-T/(1+T)
        S=-1
        A=HF*LOG(1+T)**2
       ELSE IF(T .LE. 1) THEN
        Y=T
        S=1
        A=0
       ELSE
        Y=1/T
        S=-1
        A=PI6+HF*LOG(T)**2
       ENDIF
       H=Y+Y-1
       ALFA=H+H
       B1=0
       B2=0
       DO 1 I = 19,0,-1
       B0=C(I)+ALFA*B1-B2
       B2=B1
    1  B1=B0
       H=-(S*(B0-H*B2)+A)
      ENDIF
      DMB_DILOG=H
      RETURN
      END

C     ======================================
      double precision function dmb_gamma(x)
C     ======================================
 
C--   Cernlib routine DGAMMA C302
C--
C--   Taken from /afs/cern.ch/asis/share/cern/2001/src/...
C--           .../mathlib/gen/c/gamm64.F

      implicit double precision (a-h,o-z)
C
      DIMENSION C(0:15)

      DATA C( 0) /3.65738 77250 83382 44D0/
      DATA C( 1) /1.95754 34566 61268 27D0/
      DATA C( 2) /0.33829 71138 26160 39D0/
      DATA C( 3) /0.04208 95127 65575 49D0/
      DATA C( 4) /0.00428 76504 82129 09D0/
      DATA C( 5) /0.00036 52121 69294 62D0/
      DATA C( 6) /0.00002 74006 42226 42D0/
      DATA C( 7) /0.00000 18124 02333 65D0/
      DATA C( 8) /0.00000 01096 57758 66D0/
      DATA C( 9) /0.00000 00059 87184 05D0/
      DATA C(10) /0.00000 00003 07690 81D0/
      DATA C(11) /0.00000 00000 14317 93D0/
      DATA C(12) /0.00000 00000 00651 09D0/
      DATA C(13) /0.00000 00000 00025 96D0/
      DATA C(14) /0.00000 00000 00001 11D0/
      DATA C(15) /0.00000 00000 00000 04D0/

      U=X
      IF(U .LE. 0) THEN
       WRITE(6,'(/'' DMB_GAMMA: negative argument ='',E15.5,'//
     +         '  '' ---> STOP'')') U
       STOP
      ENDIF
      F=1
      IF(U .LT. 3) THEN
       DO 1 I = 1,INT(4-U)
       F=F/U
    1  U=U+1
      ELSE
       DO 2 I = 1,INT(U-3)
       U=U-1
    2  F=F*U
      END IF
      H=U+U-7
      ALFA=H+H
      B1=0
      B2=0
      DO 3 I = 15,0,-1
      B0=C(I)+ALFA*B1-B2
      B2=B1
    3 B1=B0
      DMB_GAMMA=F*(B0-H*B2)
      RETURN
      END

C     =========================================
      subroutine smb_deriv(f,x,delta,dfdx,rerr)
C     =========================================
 
C--   Cernlib routine DDERIV D401
C--
C--   Taken from /afs/cern.ch/asis/share/cern/2001/src/...
C--           .../mathlib/gen/d/deriv64.F and gausscod.inc

      implicit double precision (a-h,o-z)

      external f

C     Computes the derivative f'(x) of f(x) at x = X. Based on
C     H. Rutishauser, Ausdehnung des Rombergschen Prinzips
C     (Extension of Romberg's Principle), Numer. Math. 5 (1963) 48-54

      DIMENSION DX(0:9),W(0:9,3),T(0:9,0:9),A(0:9)
      LOGICAL LEV(0:9),LMT

      PARAMETER (EPS = 5D-14)
      PARAMETER (Z1 = 1, S = Z1/10)

      DATA DX /0.0256D0, 0.0192D0, 0.0128D0, 0.0096D0, 0.0064D0,
     1         0.0048D0, 0.0032D0, 0.0024D0, 0.0016D0, 0.0012D0/

      DATA (LEV(K),K=0,8,2) /5*.TRUE./
      DATA (LEV(K),K=1,9,2) /5*.FALSE./

      DATA W(1,1) /1.33333 33333 333333D+00/
      DATA W(3,1) /1.06666 66666 666667D+00/
      DATA W(5,1) /1.01587 30158 730159D+00/
      DATA W(7,1) /1.00392 15686 274510D+00/

      DATA W(2,1) /3.33333 33333 333333D-01/
      DATA W(4,1) /6.66666 66666 666667D-02/
      DATA W(6,1) /1.58730 15873 015873D-02/
      DATA W(8,1) /3.92156 86274 509804D-03/

      DATA W(0,2) /2.28571 42857 142857D+00/
      DATA W(2,2) /1.16363 63636 363636D+00/
      DATA W(4,2) /1.03643 72469 635628D+00/
      DATA W(6,2) /1.00886 69950 738916D+00/
      DATA W(8,2) /1.00220 21042 329337D+00/

      DATA W(1,2) /1.28571 42857 142857D+00/
      DATA W(3,2) /1.63636 36363 636364D-01/
      DATA W(5,2) /3.64372 46963 562753D-02/
      DATA W(7,2) /8.86699 50738 916256D-03/
      DATA W(9,2) /2.20210 42329 336922D-03/

      DATA W(0,3) /1.80000 00000 000000D+00/
      DATA W(2,3) /1.12500 00000 000000D+00/
      DATA W(4,3) /1.02857 14285 714286D+00/
      DATA W(6,3) /1.00699 30069 930070D+00/
      DATA W(8,3) /1.00173 91304 347826D+00/

      DATA W(1,3) /8.00000 00000 000000D-01/
      DATA W(3,3) /1.25000 00000 000000D-01/
      DATA W(5,3) /2.85714 28571 428571D-02/
      DATA W(7,3) /6.99300 69930 069930D-03/
      DATA W(9,3) /1.73913 04347 826087D-03/

      DEL=10*ABS(DELTA)
      IS=10

    4 IS=IS-1
      DEL=S*DEL
      IF(IS .EQ. 0 .OR. X+DEL*DX(9) .EQ. X) THEN
       DELTA=DEL
       DFDX=0
       RERR=1
*mb    WRITE(ERRTXT,101) X
*mb    CALL MTLPRT(NAME,'D401.1',ERRTXT)
       WRITE(6,'(/'' SMB_DDERIV: failure for X = '',D15.8,'//
     +         '  '' ---> STOP'')') X
       STOP
      ENDIF
      DO 1 K = 0,9
      H=DEL*DX(K)
      T(K,0)=(F(X+H)-F(X-H))/(H+H)
    1 A(K)=T(K,0)

      IF(A(0) .GE. A(9)) THEN
       DO 5 K = 0,9
    5  A(K)=-A(K)
      ENDIF

      LMT=.TRUE.
      DO 3 K = 1,9
      H=A(K-1)-A(K)
    3 LMT=LMT .AND. (H .LE. 0 .OR. ABS(H) .LE. EPS*ABS(A(K)))
      IF(.NOT.LMT) GO TO 4

      DO 2 M = 1,9
      DO 2 K = 0,9-M
      IF(LEV(M)) THEN
       T(K,M)=W(M-1,1)*T(K+1,M-1)-W(M,1)*T(K,M-1)
      ELSEIF(LEV(K)) THEN
       T(K,M)=W(M-1,2)*T(K+1,M-1)-W(M,2)*T(K,M-1)
      ELSE
       T(K,M)=W(M-1,3)*T(K+1,M-1)-W(M,3)*T(K,M-1)
      ENDIF
    2 CONTINUE
      DFDX=T(0,9)
      RERR=0
      IF(DFDX .NE. 0) RERR=(DFDX-T(0,8))/DFDX
      DELTA=DEL
      RETURN
      END

C     ==============================================
      double precision function dmb_gauss(f,a,b,eps)
C     ==============================================
 
C--   Cernlib routine DGAUSS D103
C--
C--   Taken from /afs/cern.ch/asis/share/cern/2001/src/...
C--           .../mathlib/gen/d/gauss64.F and gausscod.inc

      implicit double precision (a-h,o-z)

      external f

      DIMENSION W(12),X(12)

      PARAMETER (Z1 = 1, HF = Z1/2, CST = 5*Z1/1000)

      DATA X( 1) /9.6028985649753623D-1/, W( 1) /1.0122853629037626D-1/
      DATA X( 2) /7.9666647741362674D-1/, W( 2) /2.2238103445337447D-1/
      DATA X( 3) /5.2553240991632899D-1/, W( 3) /3.1370664587788729D-1/
      DATA X( 4) /1.8343464249564980D-1/, W( 4) /3.6268378337836198D-1/
      DATA X( 5) /9.8940093499164993D-1/, W( 5) /2.7152459411754095D-2/
      DATA X( 6) /9.4457502307323258D-1/, W( 6) /6.2253523938647893D-2/
      DATA X( 7) /8.6563120238783174D-1/, W( 7) /9.5158511682492785D-2/
      DATA X( 8) /7.5540440835500303D-1/, W( 8) /1.2462897125553387D-1/
      DATA X( 9) /6.1787624440264375D-1/, W( 9) /1.4959598881657673D-1/
      DATA X(10) /4.5801677765722739D-1/, W(10) /1.6915651939500254D-1/
      DATA X(11) /2.8160355077925891D-1/, W(11) /1.8260341504492359D-1/
      DATA X(12) /9.5012509837637440D-2/, W(12) /1.8945061045506850D-1/

      H=0
      IF(B .EQ. A) GO TO 99
      CONST=CST/ABS(B-A)
      BB=A
    1 AA=BB
      BB=B
    2 C1=HF*(BB+AA)
      C2=HF*(BB-AA)
      S8=0
      DO 3 I = 1,4
      U=C2*X(I)
    3 S8=S8+W(I)*(F(C1+U)+F(C1-U))
      S16=0
      DO 4 I = 5,12
      U=C2*X(I)
    4 S16=S16+W(I)*(F(C1+U)+F(C1-U))
      S16=C2*S16
      IF(ABS(S16-C2*S8) .LE. EPS*(1+ABS(S16))) THEN
       H=H+S16
       IF(BB .NE. B) GO TO 1
      ELSE
       BB=C1
       IF(1+CONST*ABS(C2) .NE. 1) GO TO 2
       H=0
*mb    CALL MTLPRT(NAME,'D103.1','TOO HIGH ACCURACY REQUIRED')
*mb    GO TO 99
       WRITE(6,'(/'' DMB_GAUSS: too high accuracy required'','//
     +         '  '' ---> STOP'')')
       STOP
      END IF

   99 DMB_GAUSS=H
      RETURN
      END

C     =======================================
      subroutine smb_dminv(n,a,idim,ir,ifail)
C     =======================================
 
C--   Cernlib routine DINV F010
C--
C--   Taken from /afs/cern.ch/asis/share/cern/2001/src/...
C--           .../packlib/kernlib/kernnum/f010fort/dinv.F
C--
C--   Calls F011 DFACT and DFINV (also included in mbclib.f)

*mb   REAL R(N),T1,T2,T3
      INTEGER IR(N)
      REAL T1,T2,T3
      DOUBLE PRECISION A(IDIM,N),DET,TEMP,S,
     $                 C11,C12,C13,C21,C22,C23,C31,C32,C33
C
C     ******************************************************************
C
C     REPLACES A BY ITS INVERSE.
C
C     (PARAMETERS AS FOR DEQINV.)
C
C     CALLS ... DFACT, DFINV, F010PR, ABEND.
C
C     ******************************************************************
C
C  TEST FOR PARAMETER ERRORS.
C
      IF((N.LT.1).OR.(N.GT.IDIM)) GO TO 7
C
C  TEST FOR N.LE.3.
C
      IF(N.GT.3) GO TO 6
      IFAIL=0
      IF(N.LT.3) GO TO 4
C
C  N=3 CASE.
C
C     COMPUTE COFACTORS.
      C11=A(2,2)*A(3,3)-A(2,3)*A(3,2)
      C12=A(2,3)*A(3,1)-A(2,1)*A(3,3)
      C13=A(2,1)*A(3,2)-A(2,2)*A(3,1)
      C21=A(3,2)*A(1,3)-A(3,3)*A(1,2)
      C22=A(3,3)*A(1,1)-A(3,1)*A(1,3)
      C23=A(3,1)*A(1,2)-A(3,2)*A(1,1)
      C31=A(1,2)*A(2,3)-A(1,3)*A(2,2)
      C32=A(1,3)*A(2,1)-A(1,1)*A(2,3)
      C33=A(1,1)*A(2,2)-A(1,2)*A(2,1)
      T1=ABS(SNGL(A(1,1)))
      T2=ABS(SNGL(A(2,1)))
      T3=ABS(SNGL(A(3,1)))
C
C     (SET TEMP=PIVOT AND DET=PIVOT*DET.)
      IF(T1.GE.T2) GO TO 1
         IF(T3.GE.T2) GO TO 2
C        (PIVOT IS A21)
            TEMP=A(2,1)
            DET=C13*C32-C12*C33
            GO TO 3
    1 IF(T3.GE.T1) GO TO 2
C     (PIVOT IS A11)
         TEMP=A(1,1)
         DET=C22*C33-C23*C32
         GO TO 3
C     (PIVOT IS A31)
    2    TEMP=A(3,1)
         DET=C23*C12-C22*C13
C
C     SET ELEMENTS OF INVERSE IN A.
    3 IF(DET.EQ.0D0) GO TO 8
      S=TEMP/DET
      A(1,1)=S*C11
      A(1,2)=S*C21
      A(1,3)=S*C31
      A(2,1)=S*C12
      A(2,2)=S*C22
      A(2,3)=S*C32
      A(3,1)=S*C13
      A(3,2)=S*C23
      A(3,3)=S*C33
      RETURN
C
    4 IF(N.LT.2) GO TO 5
C
C  N=2 CASE BY CRAMERS RULE.
C
      DET=A(1,1)*A(2,2)-A(1,2)*A(2,1)
      IF(DET.EQ.0D0) GO TO 8
      S=1D0/DET
      C11   =S*A(2,2)
      A(1,2)=-S*A(1,2)
      A(2,1)=-S*A(2,1)
      A(2,2)=S*A(1,1)
      A(1,1)=C11
      RETURN
C
C  N=1 CASE.
C
    5 IF(A(1,1).EQ.0D0) GO TO 8
      A(1,1)=1D0/A(1,1)
      RETURN
C
C  N.GT.3 CASES.  FACTORIZE MATRIX AND INVERT.
C
*mb 6 CALL SMB_DFACT(N,A,IDIM,R,IFAIL,DET,JFAIL)
    6 CALL SMB_DFACT(N,A,IDIM,IR,IFAIL,DET,JFAIL)
      IF(IFAIL.NE.0) RETURN
*mb   CALL SMB_DFINV(N,A,IDIM,R)
      CALL SMB_DFINV(N,A,IDIM,IR)
      RETURN
C
C  ERROR EXITS.
C
    7 IFAIL=+1
*mb   CALL F010PR(NAME,N,IDIM,K,KPRNT)
      RETURN
C
    8 IFAIL=-1
      RETURN
C
      END

C     =================================================
      subroutine smb_dfact(n,a,idim,ir,ifail,det,jfail)
C     =================================================

C--   Cernlib routine DFACT F011
C--
C--   Taken from /afs/cern.ch/asis/share/cern/2001/src/...
C--           .../packlib/kernlib/kernnum/f011fort/dfact.F and fact.inc
C--
C--   Used in SMB_DINV

      INTEGER             IR(*),    IPAIRF
      DOUBLE PRECISION    A(IDIM,*),DET,      ZERO,     ONE,X,Y,TF
      REAL                G1,       G2
      REAL                PIVOTF,   P,        Q,        SIZEF,  T
      DOUBLE PRECISION    S11, S12, DOTF
      IPAIRF(J,K)  =  J*2**12 + K
      PIVOTF(X)    =  ABS(SNGL(X))
      SIZEF(X)     =  ABS(SNGL(X))
      DOTF(X,Y,S11)  =  X * Y + S11

      DATA      G1, G2              /  1.E-19,  1.E19  /

      DATA      ZERO, ONE           /  0.D0, 1.D0  /
      DATA      NORMAL, IMPOSS      /  0, -1  /
      DATA      JRANGE, JOVER, JUNDER  /  0, +1, -1  /
*
* fact.inc
*
      IF(IDIM .GE. N  .AND.  N .GT. 0)  GOTO 110
         WRITE(6,'('' SMB_DFACT n ='',I10,'' not in range [ 1 -'',
     +          I10,'' ] ---> STOP'')') N, IDIM
         RETURN
 110  IFAIL  =  NORMAL
      JFAIL  =  JRANGE
      NXCH   =  0
      DET    =  ONE
      DO 144    J  =  1, N
         K  =  J
         P  =  PIVOTF(A(J,J))
         IF(J .EQ. N)  GOTO 122
         JP1  =  J+1
         DO 121    I  =  JP1, N
            Q  =  PIVOTF(A(I,J))
            IF(Q .LE. P)  GOTO 121
               K  =  I
               P  =  Q
 121        CONTINUE
         IF(K .NE. J)  GOTO 123
 122     IF(P .GT. 0.)  GOTO 130
            DET    =  ZERO
            IFAIL  =  IMPOSS
            JFAIL  =  JRANGE
            RETURN
 123     DO 124    L  =  1, N
            TF      =  A(J,L)
            A(J,L)  =  A(K,L)
            A(K,L)  =  TF
 124        CONTINUE
         NXCH      =  NXCH + 1
         IR(NXCH)  =  IPAIRF(J,K)
 130     DET     =  DET * A(J,J)
         A(J,J)  =  ONE / A(J,J)
         T  =  SIZEF(DET)
         IF(T .LT. G1)  THEN
            DET    =  ZERO
            IF(JFAIL .EQ. JRANGE)  JFAIL  =  JUNDER
         ELSEIF(T .GT. G2)  THEN
            DET    =  ONE
            IF(JFAIL .EQ. JRANGE)  JFAIL  =  JOVER
         ENDIF
         IF(J .EQ. N)  GOTO 144
         JM1  =  J-1
         JP1  =  J+1
         DO 143   K  =  JP1, N
            S11  =  -A(J,K)
            S12  =  -A(K,J+1)
            IF(J .EQ. 1)  GOTO 142
            DO 141  I  =  1, JM1
               S11  =  DOTF(A(I,K),A(J,I),S11)
               S12  =  DOTF(A(I,J+1),A(K,I),S12)
 141           CONTINUE
 142        A(J,K)    =  -S11 * A(J,J)
            A(K,J+1)  =  -DOTF(A(J,J+1),A(K,J),S12)
 143        CONTINUE
 144     CONTINUE
      IF(MOD(NXCH,2) .NE. 0)  DET  =  -DET
      IF(JFAIL .NE. JRANGE)   DET  =  ZERO
      IR(N)  =  NXCH

      RETURN
      END

C     =================================
      subroutine smb_dfinv(n,a,idim,ir)
C     =================================

C--   Cernlib routine DFINV F011
C--
C--   Taken from /afs/cern.ch/asis/share/cern/2001/src/...
C--           .../packlib/kernlib/kernnum/f011fort/dfinv.F and finv.inc
C--
C--   Used in SMB_DINV

      INTEGER             IR(*)
      DOUBLE PRECISION    A(IDIM,*),ZERO,     X, Y, TI
      DOUBLE PRECISION    S31, S32, S33, S34, DOTF
      DOTF(X,Y,S31)  =  X*Y + S31
      DATA      ZERO      /  0.D0  /

      IF(IDIM .GE. N  .AND.  N .GT. 0)  GOTO 310
         WRITE(6,'('' SMB_DFINV n ='',I10,'' not in range [ 1 -'',
     +          I10,'' ] ---> STOP'')') N, IDIM
         RETURN
 310  IF(N .EQ. 1)  RETURN
      A(2,1)  =  -A(2,2) * DOTF(A(1,1),A(2,1),ZERO)
      A(1,2)  =  -A(1,2)
      IF(N .EQ. 2)  GOTO 330
      DO 314    I  =  3, N
         IM2  =  I-2
         DO 312 J  =  1, IM2
            S31  =  ZERO
            S32  =  A(J,I)
            DO 311  K  =  J, IM2
               S31  =  DOTF(A(K,J),A(I,K),S31)
               S32  =  DOTF(A(J,K+1),A(K+1,I),S32)
 311           CONTINUE
            A(I,J)  =  -A(I,I) * DOTF(A(I-1,J),A(I,I-1),S31)
            A(J,I)  =  -S32
 312        CONTINUE
         A(I,I-1)  =  -A(I,I) * DOTF(A(I-1,I-1),A(I,I-1),ZERO)
         A(I-1,I)  =  -A(I-1,I)
 314     CONTINUE
 330  NM1  =  N-1
      DO 335   I  =  1, NM1
         NMI  =  N-I
         DO 332   J  =  1, I
            S33  =  A(I,J)
            DO 331   K  =  1, NMI
               S33  =  DOTF(A(I+K,J),A(I,I+K),S33)
 331           CONTINUE
            A(I,J)  =  S33
 332        CONTINUE
         DO 334   J  =  1, NMI
            S34  =  ZERO
            DO 333   K  =  J, NMI
               S34  =  DOTF(A(I+K,I+J),A(I,I+K),S34)
 333           CONTINUE
            A(I,I+J)  =  S34
 334        CONTINUE
 335     CONTINUE
      NXCH  =  IR(N)
      IF(NXCH .EQ. 0)  RETURN
        DO 342 M  =  1, NXCH
         K   =  NXCH - M+1
         IJ  =  IR(K)
         I   =  IJ / 4096
         J   =  MOD(IJ,4096)
         DO 341  K  =  1, N
            TI      =  A(K,I)
            A(K,I)  =  A(K,J)
            A(K,J)  =  TI
 341        CONTINUE
 342     CONTINUE

      RETURN
      END

C     ====================================
      subroutine smb_dsinv(n,a,idim,ifail)
C     ====================================

C--   Cernlib routine DSINV F012
C--
C--   Taken from /afs/cern.ch/asis/share/cern/2001/src/...
C--           .../packlib/kernlib/kernnum/f012fort/dsinv.F
C--                                                sfact.inc 
C--                                                sfinv.inc

      DOUBLE PRECISION    A(IDIM,*),  ZERO,  ONE,  X, Y
      DOUBLE PRECISION    S1, S31, S32, S33,  DOTF
      DOTF(X,Y,S1)  =  X * Y + S1
      DATA      ZERO, ONE           /  0.D0, 1.D0 /

      IF(IDIM .LT. N  .OR.  N .LE. 0)  GOTO 900

*
* sfact.inc
*
      IFAIL  =  0
      JP1    =  0
      DO 144    J  =  1, N
         IF((A(J,J)) .LE. ZERO)  GOTO 150
         A(J,J)  =  ONE / A(J,J)
         IF(J .EQ. N)  GOTO 199
         JP1  =  J+1
         DO 143   L  =  JP1, N
            A(J,L)  =  A(J,J)*A(L,J)
            S1      =  -A(L,J+1)
            DO 141  I  =  1, J
               S1  =  DOTF(A(L,I),A(I,J+1),S1)
 141           CONTINUE
            A(L,J+1)  =  -S1
 143        CONTINUE
 144     CONTINUE
 150  IFAIL  =  -1
      RETURN
 199  CONTINUE

*
* sfinv.inc
*
      IF(N .EQ. 1)  GOTO 399
      A(1,2)  =  -A(1,2)
      A(2,1)  =   A(1,2)*A(2,2)
      IF(N .EQ. 2)  GOTO 320
      DO 314    J  =  3, N
         JM2  =  J - 2
         DO 312 K  =  1, JM2
            S31  =  A(K,J)
            DO 311  I  =  K, JM2
               S31  =  DOTF(A(K,I+1),A(I+1,J),S31)
 311           CONTINUE
            A(K,J)  =  -S31
            A(J,K)  =  -S31*A(J,J)
 312        CONTINUE
         A(J-1,J)  =  -A(J-1,J)
         A(J,J-1)  =   A(J-1,J)*A(J,J)
 314     CONTINUE
 320  J  =  1
 323     S33  =  A(J,J)
         IF(J .EQ. N)  GOTO 325
         JP1  =  J + 1
         DO 324 I  =  JP1, N
            S33  =  DOTF(A(J,I),A(I,J),S33)
 324        CONTINUE
 325     A(J,J)  =  S33
      JM1  =  J
      J    =  JP1
         DO 328 K  =  1, JM1
            S32  =  ZERO
            DO 327  I  =  J, N
               S32  =  DOTF(A(K,I),A(I,J),S32)
 327           CONTINUE
            A(K,J)  =  S32
            A(J,K)  =  S32
 328        CONTINUE
      IF(J .LT. N)  GOTO 323
 399  CONTINUE
      RETURN

 900  CONTINUE
         WRITE(6,'('' SMB_DSINV n ='',I10,'' not in range [ 1 -'',
     +           I10,'' ] ---> STOP'')') N, IDIM

      STOP  

      END

C     ========================================
      subroutine smb_dseqn(n,a,idim,ifail,k,b)
C     ========================================

C--   Cernlib routine DSEQN F012
C--
C--   Taken from /afs/cern.ch/asis/share/cern/2001/src/...
C--           .../packlib/kernlib/kernnum/f012fort/dseqn.F
C--                                                sfact.inc 
C--                                                sfeqn.inc

          DOUBLE PRECISION    A(IDIM,*), B(IDIM,*),  ONE,  X, Y
          DOUBLE PRECISION    S1, S21, S22,       DOTF
          DOTF(X,Y,S1)  =  X * Y + S1
          DATA      ZERO, ONE           /  0.D0, 1.D0 /
          IF(IDIM .LT. N  .OR.  N .LE. 0  .OR.  K .LT. 0)  GOTO 900
*
* sfact.inc
*
          IFAIL  =  0
          DO 144    J  =  1, N
             IF((A(J,J)) .LE. ZERO)  GOTO 150
             A(J,J)  =  ONE / A(J,J)
             IF(J .EQ. N)  GOTO 199
             JP1  =  J+1
             DO 143   L  =  JP1, N
                A(J,L)  =  A(J,J)*A(L,J)
                S1      =  -A(L,J+1)
                DO 141  I  =  1, J
                   S1  =  DOTF(A(L,I),A(I,J+1),S1)
 141               CONTINUE
                A(L,J+1)  =  -S1
 143            CONTINUE
 144         CONTINUE
 150      IFAIL  =  -1
          RETURN
 199      CONTINUE
*
* sfeqn.inc
*
          IF(K .LE. 0)  GOTO 299
          DO 220    L  =  1, K
             B(1,L)  =  A(1,1)*B(1,L)
 220         CONTINUE
          IF(N .EQ. 1)  GOTO 299
          DO 243    L  =  1, K
             DO 232   I  =  2, N
                IM1  =  I-1
                S21  =  - B(I,L)
                DO 231   J  =  1, IM1
                   S21  =  DOTF(A(I,J),B(J,L),S21)
 231               CONTINUE
                B(I,L)  =  - A(I,I)*S21
 232            CONTINUE
             NM1  =  N-1
             DO 242   I  =  1, NM1
                NMI  =  N-I
                S22  =  - B(NMI,L)
                DO 241   J  =  1, I
                   NMJP1  =  N - J+1
                   S22    =  DOTF(A(NMI,NMJP1),B(NMJP1,L),S22)
 241               CONTINUE
                B(NMI,L)  =  - S22
 242            CONTINUE
 243         CONTINUE
 299      CONTINUE

          RETURN

 900     CONTINUE
         WRITE(6,'('' SMB_DSEQN inconsistent input n, idim, k ='',
     +            3I10,'' ---> STOP'')') N,IDIM,K

      STOP  

      RETURN
      END

C     ================================================
      double precision function dmb_vnorm(vec,n,inorm)
C     ================================================

C--   Calculate ||v||_p = [Sum_i |v_i|^p]^(1/p) for p >= 1
C--                        Max_i |v_i|          otherwise
C--
C--   vec    (in)  input vector of at least dimension n
C--   n      (in)  number of elements in vec
C--   inorm  (in)  norm definition 
C--
C--   remark  n=1  corresponds to city-block normalisation, that is,
C--                the distance to walk from A to B in New York.
C--           n=2  corresponds to the familiar Euclidian norm.
C--           n<=0 corresponds to ||v||_{infinity} = Max |v_i|.

C--   Author: Michiel Botje h24@nikhef.nl   05-10-10
              
      implicit double precision (a-h,o-z)
      
      dimension vec(*)
      
      if(n.le.0) stop ' DMB_VNORM input vector dimension .le. zero' 
      
      vn = 0.D0
      
      if    (inorm .eq. 1) then        !Cityblock
        do i = 1,n
          vn = vn + abs(vec(i))
        enddo
      elseif(inorm .eq. 2) then        !Euclidian 
        do i = 1,n
          vn = vn + vec(i)*vec(i)
        enddo
        vn = sqrt(vn)
      elseif(inorm .ge. 3) then        !P-norm
        do i = 1,n
          vn = vn + abs(vec(i))**inorm
        enddo
        vn = vn**(1.D0/inorm)
      else                             !Major
        do i = 1,n
          vn = max(vn,abs(vec(i)))
        enddo
      endif
      
      dmb_vnorm = vn
      
      return
      end        

      
