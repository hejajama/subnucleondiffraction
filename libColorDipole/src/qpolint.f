
C--   This is the file qpolint.f containing interpolation routines 

C--   subroutine sqcWeedIt(yin,tin,nin,yout,tout,infout,nout)
C--   subroutine sqcIntLst(jset,id,iy1,iy2,iz1,iz2,y,t,f,n)
C--   subroutine sqcIntPol(jset,id,iy1,iy2,iz1,iz2,y,t,f)
C--   subroutine sqcZmesh(y,t,margin,iymi,iyma,izmi,izma)
C--   subroutine sqcMarkit(mark,ylst,tlst,margin,iy1,iy2,iz1,iz2,n)

C     =======================================================
      subroutine sqcWeedIt(yin,tin,nin,yout,tout,infout,nout)
C     =======================================================

C--   Weed points outside x-mu2 grid, including x = 1
C--
C--   yin(i),tin(i)    (in)  : list of y-t points
C--   nin              (in)  : number of points in yin,tin
C--   yout(j),tout(j)  (out) : list of y-t points inside grid
C--   infout(j)        (out) : input index i from output index j
C--   nout             (out) : number of points in yout,tout

      implicit double precision (a-h,o-z)
      
      logical lqcInside
      
      dimension yin(*),tin(*),yout(*),tout(*),infout(*)
      
      nout  = 0
      do i = 1,nin
        xx = exp(-yin(i))
        qq = exp(tin(i))        
        if(lqcInside(xx,qq,ifail)) then
          nout         = nout+1   
          yout(nout)   = yin(i)
          tout(nout)   = tin(i)
          infout(nout) = i
        endif
      enddo

      return
      end

C     =====================================================
      subroutine sqcIntLst(jset,id,iy1,iy2,iz1,iz2,y,t,f,n)
C     =====================================================

C--   List of polynomial interpolations on a table 
C--
C--   jset       (in)  :  pdf set identifier
C--   id         (in)  :  identifier of table to be interpolated
C--   iy1(i)     (in)  :  lower limit of interpolation mesh in y
C--   iy2(i)     (in)  :  upper limit of interpolation mesh in y
C--   iz1(i)     (in)  :  lower limit of interpolation mesh in z
C--   iz2(i)     (in)  :  upper limit of interpolation mesh in z
C--   y(i),t(i)  (in)  :  interpolation point
C--   f(i)       (out) :  value of the interpolation polynomial P(y,t)
C--   n          (in)  :  number of items in the list

      implicit double precision (a-h,o-z)
      
      dimension iy1(*),iy2(*),iz1(*),iz2(*),y(*),t(*),f(*)
      
      do i = 1,n
        call sqcIntPol(jset,id,iy1(i),iy2(i),iz1(i),iz2(i),
     +                 y(i),t(i),f(i))
      enddo
      
      return
      end

C     ===================================================
      subroutine sqcIntPol(jset,id,iy1,iy2,iz1,iz2,y,t,f)
C     ===================================================

C--   Polynomial interpolation on a table 
C--
C--   jset (in)  :  pdf set identifier
C--   id   (in)  :  identifier of table to be interpolated
C--   iy1  (in)  :  lower limit of interpolation mesh in y
C--   iy2  (in)  :  upper limit of interpolation mesh in y
C--   iz1  (in)  :  lower limit of interpolation mesh in z
C--   iz2  (in)  :  upper limit of interpolation mesh in z
C--   y,t  (in)  :  interpolation point
C--   f    (out) :  value of the interpolation polynomial P(y,t)
C--
C--   NB: y,t should be within the interpolation mesh, otherwise
C--       we have extrapolation. The routine sqcZmesh figures out
C--       what the interpolation mesh boundaries are for given y,t

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension work(iom0)

      ny = iy2-iy1+1
      nz = iz2-iz1+1
C--   Do nz interpolations in y
      j  = 0
      do iz = iz1,iz2
        j    = j+1
        iadr = iqcPdfIjkl(iy1,iz,id,jset)
*mb     call sqcPolint(ygrid2(iy1),stor7(iadr),ny,y,work(j),eps)
        work(j) = dqcPolint(y,ygrid2(iy1),stor7(iadr),ny)
      enddo
C--   Do one interpolation in t
*mb   call sqcPolint(zgrid2(iz1),work,nz,t,f,eps)
      f = dqcPolint(t,zgrid2(iz1),work,nz)
      
      return
      end
      
C     ===================================================
      subroutine sqcZmesh(y,t,margin,iymi,iyma,izmi,izma)
C     ===================================================

C--   Find the boundaries of a iosp*3 interpolation mesh in (iy,iz) 
C--   around y and t. This is not entirely trivial because the mesh
C--   should not  cross the grid boundaries or thresholds. It may
C--   happen that the mesh is too wide to fit between two thresholds or
C--   between a threshold and a boundary in which case we must adjust
C--   the width of the mesh. The interpolation point (y,t) must be
C--   inside the boundaries of the qcdnum y-t grid.
C--   The margin parameter defines the distance between the upper edge
C--   of the mesh and a threshold. Margin = 0 the upper edge can
C--   touch the threshold, margin = 1 the edge stays one point away
C--   from the threshold. Since the thresholda are at least two points
C--   apart, there should be no problem with margin = 1.

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qmaps8.inc'

      if(margin.ne.0 .and. margin.ne.1)
     +   stop  'sqcZmesh: invalid margin'

C--   Get bin in y
      iy = iqcFindIy(y)
      if(iy.eq.-1) stop 'sqcZmesh: y out of range ---> STOP'
C--   Mesh width is equal to the spline order ioy2
      iyma = min(iy+ioy2-1,nyy2(0))
      iymi = max(iyma-ioy2+1,0)

C--   Get bin in t
      it = iqcItfrmt(t)
      if(it.eq.0)  stop 'sqcZmesh: t out of range ---> STOP'
      iz  = izfit2(it)
      nf  = nffiz2(iz)
      iz1 = izminf8(nf)
      iz2 = izmaxf8(nf)
*mb      iz1 = 0
*mb      iz2 = 0
*mbC--   Find upper and lower limit of the nf range
*mb      do i = 1,nlist8
*mb        if(nflist8(i).eq.nf) then
*mb          iz1 = izmin8(i)
*mb          iz2 = izmax8(i)
*mb        endif
*mb      enddo
*mb      if(iz1.eq.0 .or. iz2.eq.0) 
*mb     &   stop 'sqcZmesh: zero nf range  ---> STOP'
C--   Mesh width is always 3 for interpolation in t
      izma = min(iz+2,iz2-margin)
      izmi = max(izma-2,iz1)
      if(izma.le.izmi)
     +  stop 'sqcZmesh: zero or negative mesh width in t ---> STOP'
      return
      end
      
C     =============================================================
      subroutine sqcMarkit(mark,ylst,tlst,margin,iy1,iy2,iz1,iz2,n)
C     =============================================================

C--   Process list of interpolation points and mark gridpoints
C--
C--   mark    (out)  : logical table containing the marks in y and z
C--   ylst(n) (in)   : list of y points
C--   tlst(n) (in)   : list of t points
C--   margin  (in)   : distance away from thresholds [0,1]
C--   iy1(n)  (out)  : list of lower mesh limits in y
C--   iy2(n)  (out)  : list of upper mesh limits in y
C--   iz1(n)  (out)  : list of lower mesh limits in z
C--   iz2(n)  (out)  : list of upper mesh limits in z
C--   n       (in)   : number of points in the list

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension ylst(*),tlst(*),iy1(*),iy2(*),iz1(*),iz2(*)
      logical   mark(0:mxx0,0:mqq0+7)

C--   Initialize
      do j = 0,mqq0+7
        do i = 0,mxx0
          mark(i,j) = .false.
        enddo
      enddo
C--   Loop over interpolation points
      do i = 1,n
        call sqcZmesh(ylst(i),tlst(i),margin,iya,iyb,iza,izb)
        iy1(i) = iya
        iy2(i) = iyb
        iz1(i) = iza
        iz2(i) = izb
        do iz = iza,izb
          do iy = iya,iyb
            mark(iy,iz) = .true.
          enddo
        enddo
      enddo

      return 
      end
