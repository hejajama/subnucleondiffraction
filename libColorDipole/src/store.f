
C--   This is the file store.f containing the MBUTIL linear store routines
C--
C--   subroutine smb_bkmat(imin,imax,karr,n,it1,it2)
C--   integer function imb_index(iarr,karr,n)
   
C     ==============================================
      subroutine smb_bkmat(imin,imax,karr,n,it1,it2)
C     ==============================================

C--   Book an n-dimensional matrix M in a linear store S.
C--   In Fortran syntax M is defined as
C--
C--           dimension M(i1:j1,i2:j2,...,in:jn)
C--
C--   An element of the matrix M is stored in S through  
C--
C--   M(m1,m2,...,mn) = S(k0 + k1*m1 + ... + kn*mn). 
C--
C--   This routine calculates the factors k0, ..., kn, given the lower
C--   and upper limits of each index, and the address of first element of 
C--   the matrix M.
C--
C--   Input:  imin = array of n lower index values.
C--           imax = array of n upper index values.
C--           n    = dimension of matrix M.
C--           it1  = address of first element i.e. address of M(i1,...,in).
C--
C--   Output  karr = weight factors; should be dimensioned karrr(0:n).
C--           it2  = address of last element i.e. address of M(j1,...,jn).
C--
C--   Author: Michiel Botje h24@nikhef.nl   22-07-09

      implicit double precision (a-h,o-z)

      dimension imin(n), imax(n), karr(0:n)

C--   Check limits of index range (i .le. j):
      do i = 1,n
        if(imin(i).gt.imax(i)) then
          write(6,'(/'' SMB_BKMAT: lower .gt. upper index ---> STOP'')')
          stop
        endif
      enddo

      karr(0) = it1-imin(1)
      karr(1) = 1
      do i = 2,n
        karr(i) = karr(i-1)*(imax(i-1)-imin(i-1)+1)
        karr(0) = karr(0)-imin(i)*karr(i)
      enddo

      it2 = it1 + karr(n)*(imax(n)-imin(n)+1)-1
      
C--   Weed dummy indices, that is, those for which imin = imax
C--   For these indices karr(i)*ival is added to karr(0) after
C--   which karr(i) is set to zero. This allows for fast addressing since
C--   you can ignore the null coefficients and thereby save multiplications.
C--   Provided of course that you know beforehand which indices are dummy.      

      do i = 1,n
        if(imin(i).eq.imax(i)) then
          karr(0) = karr(0) + imin(i)*karr(i)
          karr(i) = 0
        endif
      enddo
          
      return 
      end

C     =======================================
      integer function imb_index(iarr,karr,n)
C     =======================================

C--   Calculate the address of M(m1,m2,...,mn). 
C--   The matrix M is defined in s/r smb_bkmat.
C--
C--   Input:  iarr = array of n index values.
C--           karr = weight factors calculated in s/r smb_bkmat.
C--           n    = dimension of matrix M.
C--
C--   Output: imb_index = address of M(i1,i2,...,in).
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      implicit double precision (a-h,o-z)

      dimension iarr(n),karr(0:n)

      ia = karr(0)
      do i = 1,n
        ia = ia + karr(i)*iarr(i)
      enddo
      imb_index = ia

      return 
      end

