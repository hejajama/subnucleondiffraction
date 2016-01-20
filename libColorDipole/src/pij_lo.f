
C     ========================================
      DOUBLE PRECISION FUNCTION dqcP0GGS(X,NF)
C     ========================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      include 'qconst.inc'

!       IDUM     = NF        !avoid compiler warning      
      dqcP0GGS = 1.D0 / ( 1.D0 - X )
 
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION dqcP0GGR(X,NF)
C     ========================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      include 'qconst.inc'

!       IDUM     = NF        !avoid compiler warning       
      dqcP0GGR = 6.D0 * X 
 
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION dqcP0GGA(X,NF)
C     ========================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      include 'qconst.inc'

!       IDUM     = NF        !avoid compiler warning       
      ONEMX    = 1.D0 - X
      dqcP0GGA = 6.D0 * ( ONEMX/X + X*ONEMX )
 
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION dqcP0GFA(X,NF)
C     ========================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      include 'qconst.inc'
 
!       IDUM     = NF        !avoid compiler warning 
      dqcP0GFA = 0!4.D0 * ( 1.D0 + (1.D0-X)*(1.D0-X) ) / ( 3.D0*X )
 
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION dqcP0FGA(X,NF)
C     ========================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      include 'qconst.inc'
 
!       IDUM     = NF        !avoid compiler warning 
      dqcP0FGA = 0!0.5D0 * ( X*X + (1.D0-X)*(1.D0-X) )
 
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION dqcP0FFS(X,NF)
C     ========================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      include 'qconst.inc'
 
!       IDUM     = NF        !avoid compiler warning 
      dqcP0FFS = 0!1.D0 / (1.D0-X)
 
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION dqcP0FFR(X,NF)
C     ========================================

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      include 'qconst.inc'

!       IDUM     = NF        !avoid compiler warning 
      dqcP0FFR = 0!4.D0 * ( 1.D0 + X*X ) / 3.D0
 
      RETURN
      END
