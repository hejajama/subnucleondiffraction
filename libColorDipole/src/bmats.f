
C--   This is the file bmats.f containing the MBUTIL band/triangular matrix routines
C--
C--   subroutine smb_LMeqs(a,na,m,x,b,n,ierr)
C--   subroutine smb_UMeqs(a,na,m,x,b,n,ierr)
C--   integer function imb_LLadr(i,j,m,n)
C--   subroutine smb_LLeqs(a,m,x,b,n,ierr)
C--   integer function imb_ULadr(i,j,m,n)
C--   subroutine smb_ULeqs(a,m,x,b,n,ierr)
C--   integer function imb_LTadr(i,j,m,n)
C--   subroutine smb_LTeqs(a,m,x,b,n,ierr)
C--   integer function imb_UTadr(i,j,m,n)
C--   subroutine smb_UTeqs(a,m,x,b,n,ierr)
C--   integer function imb_LBadr(i,j,m,n)
C--   subroutine smb_LBeqs(a,m,x,b,n,ierr)
C--   integer function imb_UBadr(i,j,m,n)
C--   subroutine smb_UBeqs(a,m,x,b,n,ierr)

C     =======================================
      subroutine smb_LMeqs(a,na,m,x,b,n,ierr)
C     ======================================= 

C--   Solves the n x n matrix equation Ax = b where A is a lower diagonal
C--   band matrix. 
C--
C--                                  | a11                 |
C--                                  | a21 a22             |
C--   For 5 x 5 and bandwidth 3, A = | a31 a32 a33         |
C--                                  |     a42 a43 a44     |
C--                                  |         a53 a54 a55 |
C--
C--   a    (in)  input 2-dim matrix with n*n submatrix filled as example above
C--   na   (in)  dimension a(na,na) declared in the calling routine (na >= n)
C--   m    (in)  bandwidth (3 in the example above, m <= n)
C--   x(n) (out) output solution vector
C--   b(n) (in)  input right-hand side vector
C--   n    (in)  dimension of a, x and b (5 in the example above)
C--   ierr (out) 0 = OK, 1 = singular matrix
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      implicit double precision (a-h,o-z)

      dimension a(na,na), x(*), b(*)

      if(a(1,1).eq.0.D0) then
        ierr = 1
        return
      endif
      ierr = 0

      x(1) = b(1)/a(1,1)

      do i = 2,n
        sum = 0.D0
        j1  = max(i+1-m,1)
        do j = j1,i-1
          sum = sum + x(j)*a(i,j)
        enddo
        if(a(i,i).eq.0.D0) then
          ierr = 1
          return
        endif
        x(i) = (b(i)-sum)/a(i,i)
      enddo

      return
      end

C     =======================================
      subroutine smb_UMeqs(a,na,m,x,b,n,ierr)
C     ======================================= 

C--   Solves the n x n matrix equation Ax = b where A is a upper diagonal
C--   band matrix.
C--
C--                                  | a11 a12 a13         |
C--                                  |     a22 a23 a24     |
C--   For 5 x 5 and bandwidth 3, A = |         a33 a34 a35 |
C--                                  |             a44 a45 |
C--                                  |                 a55 |
C--
C--   a    (in)  input 2-dim matrix with n*n submatrix filled as example above
C--   na   (in)  dimension a(na,na) declared in the calling routine (na >= n)
C--   m    (in)  bandwidth (3 in the example above, m <= n)
C--   x(n) (out) output solution vector
C--   b(n) (in)  input right-hand side vector
C--   n    (in)  dimension of a, x and b (5 in the example above)
C--   ierr (out) 0 = OK, 1 = singular matrix
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      implicit double precision (a-h,o-z)

      dimension a(na,na), x(*), b(*)
      
      if(a(n,n).eq.0.D0) then
        ierr = 1
        return
      endif
      ierr = 0
      x(n) = b(n)/a(n,n)

      do i = n-1,1,-1
        sum = 0.D0
        j2  = min(i-1+m,n)
        do j = i+1,j2
          sum = sum + x(j)*a(i,j)
        enddo
        if(a(i,i).eq.0.D0) then
          ierr = 1
          return
        endif
        x(i) = (b(i)-sum)/a(i,i)
      enddo

      return
      end

C     ===================================
      integer function imb_LLadr(i,j,m,n)
C     ===================================

      implicit double precision (a-h,o-z)

C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09      

      if(i.ge.1.and.i.le.n.and.j.ge.1.and.j.le.n .and.
     &   (i-j).ge.0.and.(i-j).le.min(m-1,n-1)) then
        imb_LLadr = i+(j-1)*n
      else
        imb_LLadr = 0
      endif

      return
      end

C     ====================================
      subroutine smb_LLeqs(a,m,x,b,n,ierr)
C     ==================================== 

C--   Solves the n x n matrix equation Ax = b where A is a lower diagonal
C--   band matrix. The matrix A is stored in a linear store
C--
C--                                  | 1            |
C--                                  | 2 7          |
C--   For 5 x 5 and bandwidth 3, A = | 3 8 13       |
C--                                  |   9 14 19    |
C--                                  |     15 20 25 |
C--
C--   a    (in)  input 1-dim linear store dim(n*n) or larger
C--   m    (in)  bandwidth (3 in the example above, m <= n)
C--   x(n) (out) output solution vector
C--   b(n) (in)  input right-hand side vector
C--   n    (in)  dimension of A, x and b (5 in the example above)
C--   ierr (out) 0 = OK, 1 = singular matrix
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      implicit double precision (a-h,o-z)

      dimension a(*), x(*), b(*)

      if(a(1).eq.0.D0) then
        ierr = 1
        return
      endif
      ierr = 0

      x(1)  = b(1)/a(1)
      jstep = n

      do i = 2,n
        sum = 0.D0
        j1  = max(i+1-m,1)
        ia  = i + (j1-1)*n
        do j = j1,i-1
          sum = sum + x(j)*a(ia)
          ia  = ia + jstep
        enddo
        div = a(i+(i-1)*n)
        if(div.eq.0.D0) then
          ierr = 1
          return
        endif
        x(i) = (b(i)-sum)/div
      enddo

      return
      end

C     ===================================
      integer function imb_ULadr(i,j,m,n)
C     ===================================

      implicit double precision (a-h,o-z)
      
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09      

      if(i.ge.1.and.i.le.n.and.j.ge.1.and.j.le.n .and.
     &   (j-i).ge.0.and.(j-i).le.min(m-1,n-1)) then
        imb_ULadr = i+(j-1)*n
      else
        imb_ULadr = 0
      endif

      return
      end

C     ====================================
      subroutine smb_ULeqs(a,m,x,b,n,ierr)
C     ==================================== 

C--   Solves the n x n matrix equation Ax = b where A is a upper diagonal
C--   band matrix. The matrix A is stored in a linear store.
C--
C--                                  | 1 6 11       |
C--                                  |   7 12 17    |
C--   For 5 x 5 and bandwidth 3, A = |     13 18 23 |
C--                                  |        19 24 |
C--                                  |           25 |
C--
C--   a    (in)  input 1-dim linear store dim(n*n) or larger
C--   m    (in)  bandwidth (3 in the example above, m <= n)
C--   x(n) (out) output solution vector
C--   b(n) (in)  input right-hand side vector
C--   n    (in)  dimension of A, x and b (5 in the example above)
C--   ierr (out) 0 = OK, 1 = singular matrix
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      implicit double precision (a-h,o-z)

      dimension a(*), x(*), b(*)

      if(a(n*n).eq.0.D0) then
        ierr = 1
        return
      endif
      ierr = 0

      x(n)  = b(n)/a(n*n)
      jstep = n

      do i = n-1,1,-1
        sum = 0.D0
        j2  = min(i-1+m,n)
        ia  = i*(n+1) 
        do j = i+1,j2
          sum = sum + x(j)*a(ia)
          ia  = ia + jstep
        enddo
        div = a(i+(i-1)*n)
        if(div.eq.0.D0) then
          ierr = 1
          return
        endif
        x(i) = (b(i)-sum)/div
      enddo

      return
      end

C     ===================================
      integer function imb_LTadr(i,j,m,n)
C     ===================================

      implicit double precision (a-h,o-z)
      
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09      

      if(i.ge.1.and.i.le.n.and.j.ge.1.and.j.le.n .and.
     &   (i-j).ge.0.and.(i-j).le.min(m-1,n-1)) then
        imb_LTadr = i*(i-1)/2 + j
      else
        imb_LTadr = 0
      endif

      return
      end

C     ====================================
      subroutine smb_LTeqs(a,m,x,b,n,ierr)
C     ==================================== 

C--   Solves the n x n matrix equation Ax = b where A is a lower diagonal
C--   band matrix. The matrix A is stored lower triangle packed
C--
C--                                  | 1             |
C--                                  | 2  3          |
C--   For 5 x 5 and bandwidth 3, A = | 4  5  6       |
C--                                  |    8  9 10    |
C--                                  |      13 14 15 |
C--
C--   a    (in)  input 1-dim linear store dim(n*n) or larger
C--   m    (in)  bandwidth (3 in the example above, m <= n)
C--   x(n) (out) output solution vector
C--   b(n) (in)  input right-hand side vector
C--   n    (in)  dimension of A, x and b (5 in the example above)
C--   ierr (out) 0 = OK, 1 = singular matrix
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      implicit double precision (a-h,o-z)

      dimension a(*), x(*), b(*)

      if(a(1).eq.0.D0) then
        ierr = 1
        return
      endif
      ierr = 0

      x(1)  = b(1)/a(1)
      jstep = 1

      do i = 2,n
        sum = 0.D0
        j1  = max(i+1-m,1)
        ia  = i*(i-1)/2 + j1
        do j = j1,i-1
          sum = sum + x(j)*a(ia)
          ia  = ia + jstep
        enddo
        div = a(i*(i+1)/2)
        if(div.eq.0.D0) then
          ierr = 1
          return
        endif
        x(i) = (b(i)-sum)/div
      enddo

      return
      end

C     ===================================
      integer function imb_UTadr(i,j,m,n)
C     ===================================

      implicit double precision (a-h,o-z)
      
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09      

      if(i.ge.1.and.i.le.n.and.j.ge.1.and.j.le.n .and.
     &   (j-i).ge.0.and.(j-i).le.min(m-1,n-1)) then
        imb_UTadr = (n+1-i)*(n-i)/2 + n + 1 - j
      else
        imb_UTadr = 0
      endif

      return
      end

C     ====================================
      subroutine smb_UTeqs(a,m,x,b,n,ierr)
C     ==================================== 

C--   Solves the n x n matrix equation Ax = b where A is a upper diagonal
C--   band matrix. The matrix A is stored upper triangle packed.
C--
C--                                  | 15 14 13       |
C--                                  |    10  9  8    |
C--   For 5 x 5 and bandwidth 3, A = |        6  5  4 |
C--                                  |           3  2 |
C--                                  |              1 |
C--
C--   a    (in)  input 1-dim linear store dim(n*n) or larger
C--   m    (in)  bandwidth (3 in the example above, m <= n)
C--   x(n) (out) output solution vector
C--   b(n) (in)  input right-hand side vector
C--   n    (in)  dimension of A, x and b (5 in the example above)
C--   ierr (out) 0 = OK, 1 = singular matrix
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      implicit double precision (a-h,o-z)

      dimension a(*), x(*), b(*)

      if(a(1).eq.0.D0) then
        ierr = 1
        return
      endif
      ierr = 0

      x(n)  = b(n)/a(1)
      jstep = -1

      do i = n-1,1,-1
        sum = 0.D0
        j2  = min(i-1+m,n)
        ia  = (n-i)*(n+3-i)/2 
        do j = i+1,j2
          sum = sum + x(j)*a(ia)
          ia  = ia + jstep
        enddo
        div = a((n+1-i)*(n+2-i)/2)
        if(div.eq.0.D0) then
          ierr = 1
          return
        endif
        x(i) = (b(i)-sum)/div
      enddo

      return
      end

C     ===================================
      integer function imb_LBadr(i,j,m,n)
C     ===================================

      implicit double precision (a-h,o-z)
      
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09      

      if(i.ge.1.and.i.le.n.and.j.ge.1.and.j.le.n .and.
     &   (i-j).ge.0.and.(i-j).le.min(m-1,n-1)) then
        imb_LBadr = (i-j)*n + i
      else
        imb_LBadr = 0
      endif

      return
      end

C     ====================================
      subroutine smb_LBeqs(a,m,x,b,n,ierr)
C     ==================================== 

C--   Solves the n x n matrix equation Ax = b where A is a lower diagonal
C--   band matrix. The matrix A is stored lower band
C--
C--                                  |  1             |
C--                                  |  7  2          |
C--   For 5 x 5 and bandwidth 3, A = | 13  8  3       |
C--                                  |    14  9  4    |
C--                                  |       15 10  5 |
C--
C--   a    (in)  input 1-dim linear store dim(n*n) or larger
C--   m    (in)  bandwidth (3 in the example above, m <= n)
C--   x(n) (out) output solution vector
C--   b(n) (in)  input right-hand side vector
C--   n    (in)  dimension of A, x and b (5 in the example above)
C--   ierr (out) 0 = OK, 1 = singular matrix
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      implicit double precision (a-h,o-z)

      dimension a(*), x(*), b(*)

      if(a(1).eq.0.D0) then
        ierr = 1
        return
      endif
      ierr = 0

      x(1)  = b(1)/a(1)
      jstep = -n

      do i = 2,n
        sum = 0.D0
        j1  = max(i+1-m,1)
        ia  = (i-j1)*n+i
        do j = j1,i-1
          sum = sum + x(j)*a(ia)
          ia  = ia + jstep
        enddo
        div = a(i)
        if(div.eq.0.D0) then
          ierr = 1
          return
        endif
        x(i) = (b(i)-sum)/div
      enddo

      return
      end

C     ===================================
      integer function imb_UBadr(i,j,m,n)
C     ===================================

      implicit double precision (a-h,o-z)
      
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09      

      if(i.ge.1.and.i.le.n.and.j.ge.1.and.j.le.n .and.
     &   (j-i).ge.0.and.(j-i).le.min(m-1,n-1)) then
        imb_UBadr = (j-i)*n + j
      else
        imb_UBadr = 0
      endif

      return
      end

C     ====================================
      subroutine smb_UBeqs(a,m,x,b,n,ierr)
C     ==================================== 

C--   Solves the n x n matrix equation Ax = b where A is a upper diagonal
C--   band matrix. The matrix A is stored upper triangle packed.
C--
C--                                  |  1  7 13       |
C--                                  |     2  8 14    |
C--   For 5 x 5 and bandwidth 3, A = |        3  9 15 |
C--                                  |           4 10 |
C--                                  |              5 |
C--
C--   a    (in)  input 1-dim linear store dim(n*n) or larger
C--   m    (in)  bandwidth (3 in the example above, m <= n)
C--   x(n) (out) output solution vector
C--   b(n) (in)  input right-hand side vector
C--   n    (in)  dimension of A, x and b (5 in the example above)
C--   ierr (out) 0 = OK, 1 = singular matrix
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      implicit double precision (a-h,o-z)

      dimension a(*), x(*), b(*)

      if(a(n).eq.0.D0) then
        ierr = 1
        return
      endif
      ierr = 0

      x(n)  = b(n)/a(n)
      jstep = n+1

      do i = n-1,1,-1
        sum = 0.D0
        j2  = min(i-1+m,n)
        ia  = n+i+1 
        do j = i+1,j2
          sum = sum + x(j)*a(ia)
          ia  = ia + jstep
        enddo
        div = a(i)
        if(div.eq.0.D0) then
          ierr = 1
          return
        endif
        x(i) = (b(i)-sum)/div
      enddo

      return
      end
