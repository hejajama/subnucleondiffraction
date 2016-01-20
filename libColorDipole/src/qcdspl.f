
C--   This is the file qcdspl.f containing the qcdnum spline utilities
C--   ----------------------------------------------------------------
C--   Set of routines to generate and access Bsplines for a given 
C--   order (kk), grid definition (yy) and grid multiplicity (mm) 
C--
C--   y-grid   -------------------------------------------------------
C--   subroutine sqcSpyIni(kk,yy,mm,ny,ns,nc)
C--   double precision function dqcBsplyy(idk,isp,y)
C--   double precision function dqcBsplyi(idk,isp,iy,y)
C--   double precision function dqcDsplyi(idk,isp,iy,y)
C--   double precision function dqcBspliy(idk,isp,iy)
C--   subroutine sqcByjLim(idk,iy,jmin,jmax)
C--   q-grid   -------------------------------------------------------
C--   subroutine sqcSpqIni(kk,qq,mm,nq,ns,nc)
C--   double precision function dqcBsplqq(isp,q,nq)
C--   double precision function dqcDsplqq(isp,q,nq)
C--   double precision function dqcBsplqi(isp,iq,q)
C--   double precision function dqcDsplqi(isp,iq,q)
C--   double precision function dqcBspliq(isp,iq)
C--   subroutine sqcBqjLim(iq,jmin,jmax)
C--   y and q grid   --------------------------------------------
C--   subroutine sqcGetTau(k,x,m,jtfx,nx,t,ixft,mt,nt,ierr)
C--   subroutine sqcSrange(k,ixft,nt,imi,ima,nx,ierr)
C--   subroutine sqcSplCat(k,t,ic,nt,nc,ierr)
C--   subroutine sqcBsplin(iord,x,tau,ntau,bsder,iom,ntm,imi,ima,ierr)
C--   subroutine sqcFilCat(k,xx,jtfx,nx,tt,ic,nt,der,cat,mk,mt,nc,ierr)
C--   double precision function 
C--                    dqcBsplxx(k,is,x,ix,jtfx,mi,ma,nx,tt,ic,cat,mk,mt)
C--   double precision function 
C--                    dqcDsplxx(k,is,x,ix,jtfx,mi,ma,nx,tt,ic,cat,mk,mt)
C--   integer function iqcBGetIx(x,xx,nx)

C--   ====================================================================
C--   1. Routines in the variable y 
C--   ====================================================================

C     =======================================
      subroutine sqcSpyIni(kk,yy,mm,ny,ns,nc)
C     =======================================

C--   Calculate Bsplines in y for later interpolation
C--
C--   kk  (in)  Spline order
C--   yy  (in)  Y-grid
C--   mm  (in)  Y-grid multiplicities
C--   ny  (in)  Number of y-grid points
C--   ns  (out) Number of generated Bsplines
C--   nc  (out) Number of splines in the catalog
C--
C--   In the following the user grid is called 'y': it is characterized
C--   set of y-points and by a multiplicity (m) for each point. This
C--   generates internally a working grid 'tau' (or 't') where each y
C--   is repeated m times. The pointer lists are shown in the example below
C--
C--   iy              1  2  3  4  5  6  7    ... y-grid index
C--   yy(1--ny)      y1 y2 y3 y4 y5 y6 y7    ... y-grid points
C--   mm(1--ny)       1  1  0  1  3  1  1    ... multiplicity
C--   jtfy(1--ny)     1  2  2  3  6  7  8    ... largest jt containig yi
C--
C--   jt              1  2  3  4  5  6  7  8 ... t-grid index
C--   tauy(1--nt)    y1 y2 y4 y5 y5 y5 y6 y7 ... t-grid points
C--   iyft(1--nt)     1  2  4  5  5  5  6  7 ... pointer to y-grid
C--
C--   Note from the above that it is possible to set the multiplicity 
C--   m_i = 0 in which case the corresponding grid-point y_i is ignored.
 
      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qcdspl.inc'

      dimension yy(*),mm(*),der(iom0,mxx0)

C--   Storage index
      if(kk.ne.2 .and. kk.ne.3) stop 
     +            'sqcSpyIni: spline order not 2 or 3 ---> STOP' 
      idk = kk-1
C--   Check dimensions
      if(kk.gt.iom0)        stop 
     +            'sqcSpyIni: spline order too large ---> STOP'
      if(ny.gt.mxx0-2*iom0) stop 
     +            'sqcSpyIni: too many y-points ---> STOP'
C--   Copy to common block
      nny(idk) = ny
      kky(idk) = kk
      do i = 1,ny
        grdy(i,idk) = yy(i)
      enddo
C--   Setup node point grid (tau) and a few pointer arrays (see above)
      call sqcGetTau(kk,yy,mm,jtfy(1,idk),ny,tauy(1,idk),iyft(1,idk),
     +               mxx0,nty(idk),ierr)
C--   Minby(iy), maxby(iy) are the first, last non-zero Bspline in bin iy
      call sqcSrange(kk,iyft(1,idk),nty(idk),minby(1,idk),maxby(1,idk),
     +               ny,ierr)
C--   Make index of catalog of Bsplines with different shapes:
C--     icaty(jt) is the address of Bspline jt in the catalog; 
C--     ncy is the total number of entries in the catalog
      call sqcSplCat(kk,tauy(1,idk),icaty(1,idk),nty(idk),ncy(idk),ierr)
C--   Now fill the catalog  
      call sqcFilCat(
     +     kk,yy,jtfy(1,idk),ny,tauy(1,idk),icaty(1,idk),nty(idk),der,
     +     caty(1,1,1,idk),iom0,mxx0,ncy(idk),ierr)
C--   Number of generated Bsplines
      ns = nty(idk)-kk
      nc = ncy(idk)

      return
      end

C     ==============================================
      double precision function dqcBsplyy(idk,isp,y)
C     ==============================================

C--   Calculate Bspline isp at y (used in the weight table routines) 

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qcdspl.inc'

      iy        = iqcBGetIx(y,grdy(1,idk),nny(idk))
      dqcBsplyy = dqcBsplxx(kky(idk),isp,y,iy,jtfy(1,idk),minby(1,idk),
     +                      maxby(1,idk),nny(idk),tauy(1,idk),
     +                      icaty(1,idk),caty(1,1,1,idk),iom0,mxx0)

      return
      end


C     =================================================
      double precision function dqcBsplyi(idk,isp,iy,y)
C     =================================================

C--   Calculate Bspline isp at y in bin iy
C--   Thus it is assumed that the index iy has been determined elsewhere

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qcdspl.inc'
      include 'qgrid2.inc'

      dqcBsplyi = dqcBsplxx(kky(idk),isp,y,iy,jtfy(1,idk),minby(1,idk),
     +                      maxby(1,idk),nny(idk),
     +                      tauy(1,idk),
     +                      icaty(1,idk),caty(1,1,1,idk),iom0,mxx0)

      return
      end

C     =================================================
      double precision function dqcDsplyi(idk,isp,iy,y)
C     =================================================

C--   Calculate dB_isp/dy at y in bin iy
C--   Thus it is assumed that the index iy has been determined elsewhere

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qcdspl.inc'

       dqcDsplyi = dqcDsplxx(kky(idk),isp,y,iy,jtfy(1,idk),minby(1,idk),
     +                      maxby(1,idk),nny(idk),tauy(1,idk),
     +                      icaty(1,idk),caty(1,1,1,idk),iom0,mxx0)

      return
      end


C     ===============================================
      double precision function dqcBspliy(idk,isp,iy)
C     ===============================================

C--   Calculate Bspline isp at iy

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qcdspl.inc'

      dqcBspliy = dqcBsplix(kky(idk),isp,iy,jtfy(1,idk),minby(1,idk),
     +                      maxby(1,idk),nny(idk),tauy(1,idk),
     +                      icaty(1,idk),caty(1,1,1,idk),iom0,mxx0)

      return
      end

C     ======================================
      subroutine sqcByjLim(idk,iy,jmin,jmax)
C     ======================================

C--   Index of first and last bspline with support in bin iy

C--   iy    (in)   ybin index
C--   jmin  (out)  first spline with support in iy
C--   jmax  (out)  last  spline with support in iy

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qcdspl.inc'

      jmin = minby(iy,idk)
      jmax = maxby(iy,idk)

      return
      end

C--   ====================================================================
C--   2. Routines in the variable q 
C--   ====================================================================

C     =======================================
      subroutine sqcSpqIni(kk,qq,mm,nq,ns,nc)
C     =======================================

C--   Calculate Bsplines in Q2 for later interpolation
C--
C--   kk  (in)  Spline order
C--   qq  (in)  Q-grid
C--   mm  (in)  Q-grid multiplicities
C--   nq  (in)  Number of q-grid points
C--   ns  (out) Number of generated Bsplines
C--   nc  (out) Number of splines in the catalog
C--
C--   The pointer lists are shown in the example below
C--
C--   iq              1  2  3  4  5  6  7    ... q-grid index
C--   qq(nq)         q1 q2 q3 q4 q5 q6 q7    ... q-grid points
C--   mm(nq)          1  1  0  1  3  1  1    ... multiplicity
C--   jtfq(nq)        1  2  2  3  6  7  8    ... largest jt containig qi
C--
C--   jt              1  2  3  4  5  6  7  8 ... t-grid index
C--   tauq(nt)       q1 q2 q4 q5 q5 q5 q6 q7 ... t-grid points
C--   iqft(nt)        1  2  4  5  5  5  6  7 ... pointer to q-grid
C--
C--   Note from the above that it is possible to set the multiplicity 
C--   m_i = 0 in which case the corresponding grid-point q_i is ignored.

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qcdspl.inc'

      dimension qq(*),mm(*),der(iom0,mqq0)

C--   Check dimensions
      if(kk.gt.iom0)        stop 
     +              'sqcSpqIni: spline order too large ---> STOP'
      if(nq.gt.mqq0-2*iom0) stop 
     +              'sqcSpqIni: too many q-points ---> STOP'
C--   Copy to common block
      nnq = nq
      kkq = kk
      do i = 1,nq
        grdq(i) = qq(i)
      enddo
C--   Setup node point grid (tau) and a few pointer arrays (see above)
      call sqcGetTau(kk,qq,mm,jtfq,nq,tauq,iqft,mqq0,ntq,ierr)
C--   Minbq(iq), maxbq(iq) are the first, last non-zero Bspline in bin iq
      call sqcSrange(kk,iqft,ntq,minbq,maxbq,nq,ierr)
C--   Make index of catalog of Bsplines with different shapes:
C--     icatq(jt) is the address of Bspline jt in the catalog; 
C--     ncq is the total number of entries in the catalog
      call sqcSplCat(kk,tauq,icatq,ntq,ncq,ierr)
C--   Now fill the catalog  
      call sqcFilCat(
     +     kk,qq,jtfq,nq,tauq,icatq,ntq,der,catq,iom0,mqq0,ncq,ierr) 
C--   Number of generated Bsplines
      ns = ntq-kk
      nc = ncq

      return
      end

C     =============================================
      double precision function dqcBsplqq(isp,q,nq)
C     =============================================

C--   Calculate Bspline isp at q (used to setup transformation matrix)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qcdspl.inc'

      iq         = iqcBGetIx(q,grdq,nq)
      dqcBsplqq  = dqcBsplxx(kkq,isp,q,iq,jtfq,minbq,maxbq,nnq,tauq,
     +                       icatq,catq,iom0,mqq0)

      return
      end

C     =============================================
      double precision function dqcDsplqq(isp,q,nq)
C     =============================================

C--   Calculate dB_isp/dq at q (used to setup transformation matrix)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qcdspl.inc'

      iq         = iqcBGetIx(q,grdq,nq)
      dqcDsplqq  = dqcDsplxx(kkq,isp,q,iq,jtfq,minbq,maxbq,nnq,tauq,
     +                       icatq,catq,iom0,mqq0)

      return
      end

C     =============================================
      double precision function dqcBsplqi(isp,iq,q)
C     =============================================

C--   Calculate Bspline isp at q in bin iq

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qcdspl.inc'

      dqcBsplqi  = dqcBsplxx(kkq,isp,q,iq,jtfq,minbq,maxbq,nnq,tauq,
     +                       icatq,catq,iom0,mqq0)

      return
      end

C     =============================================
      double precision function dqcDsplqi(isp,iq,q)
C     =============================================

C--   Calculate dB_isp/dq at q

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qcdspl.inc'

      dqcDsplqi  = dqcDsplxx(kkq,isp,q,iq,jtfq,minbq,maxbq,nnq,tauq,
     +                       icatq,catq,iom0,mqq0)

      return
      end


C     ===========================================
      double precision function dqcBspliq(isp,iq)
C     ===========================================

C--   Calculate Bspline isp at iq

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qcdspl.inc'

      dqcBspliq  = dqcBsplix(kkq,isp,iq,jtfq,minbq,maxbq,nnq,tauq,
     +                       icatq,catq,iom0,mqq0)

      return
      end

C     ==================================
      subroutine sqcBqjLim(iq,jmin,jmax)
C     ==================================

C--   Index of first and last bspline with support in bin iq

C--   iq    (in)   qbin index
C--   jmin  (out)  first spline with support in iq
C--   jmax  (out)  last  spline with support in iq

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qcdspl.inc'

      jmin = minbq(iq)
      jmax = maxbq(iq)

      return
      end

C--   ====================================================================
C--   3. General workhorse routines which can be used for both y and q
C--      but have clumsy argument lists because all the bookkeeping arrays 
C--      must be passed....... 
C--   ====================================================================

C     =====================================================
      subroutine sqcGetTau(k,x,m,jtfx,nx,t,ixft,mt,nt,ierr)
C     =====================================================

C--   Create node point list t(j) where each x(i) is listed m(i) times
C--   The pointer lists are shown in the example below
C--
C--   ix          1  2  3  4  5  6  7    ... x-grid index
C--   x          x1 x2 x3 x4 x5 x6 x7    ... x-grid points
C--   m           1  1  0  1  3  1  1    ... multiplicity
C--   jtfx        1  2  2  3  6  7  8    ... jt-bin of xi
C--
C--   jt          1  2  3  4  5  6  7  8 ... t-grid index
C--   t          x1 x2 x4 x5 x5 x5 x6 x7 ... t-grid points
C--   ixft        1  2  4  5  5  5  6  7 ... pointer to x-grid
C--
C--   k     (in)    spline order
C--   x     (in)    list of x-grid points
C--   m     (in)    list of multiplicities
C--   jtfx  (out)   jt(i) is the largest index j with t(j) = x(i)
C--   nx    (in)    number of elements of x, m and jtfx
C--   t     (out)   expanded list (tau) of node points
C--   ixft  (out)   pointer to x-grid
C--   mt    (in)    dimension of t and ixft in the calling routine
C--   nt    (out)   number of tau gridpoints = sum_i m_i
C--   ierr  (out)   0 = all OK
C--                 1 = less than two x-points
C--                 2 = xi not in ascending order
C--                 3 = one or both end-point multiplicities zero
C--                 4 = m(i) > k
C--                 5 = run out of space

      implicit double precision (a-h,o-z)
    
      dimension x(*),m(*),jtfx(*),t(*),ixft(*)

      ierr = 0
C--   Enough x-points?
      if(nx.lt.2) then
        ierr = 1
        stop 'sqcGetTau: nx .lt. 2 ---> STOP'
      endif
C--   x in strictly ascending order?
      do ix = 2,nx
        if(x(ix).le.x(ix-1)) then
          ierr = 2
          stop 'sqcGetTau: x not in ascending order ---> STOP'
        endif
      enddo
C--   No multiplicty 0 at the end points, thank you
      if((m(1).le.0) .or. (m(nx).le.0)) then
        ierr = 3
        stop 'sqcGetTau: m(1) or m(nx) are zero ---> STOP'
      endif 
      nt   = 0      
      do i = 1,nx
        jtfx(i) = nt
C--     Multiplicity cannot be larger than spline order
        if(m(i).gt.k) then
          ierr = 4
          stop 'sqcGetTau: mult larger than spline order ---> STOP'
        endif
        do j = 1,m(i)
          nt = nt+1
C--       Check that we do not run out of space
          if(nt.gt.mt) then 
            ierr = 5
            stop 'sqcGetTau: too many points in t-grid ---> STOP'
          endif
          t(nt)    = x(i)
          ixft(nt) = i
          jtfx(i)  = nt
        enddo
      enddo

      return
      end 

C     ===============================================
      subroutine sqcSrange(k,ixft,nt,imi,ima,nx,ierr)
C     ===============================================

C--   Find for each bin in x the lowest and highest index of the splines
C--   which have their support (i.e. are non-zero) in this bin.
C--   Thus: the range of splines covering x-bin (i) is imi(i) to ima(i).
C--
C--   k     (in)    spline order
C--   ixft  (in)    map from jt to ix [from sqcGetTau]
C--   nt    (in)    number of elements of tau (and ixft)
C--   imi   (out)   list of lower spline indices for xbin i
C--   ima   (out)   list of upper spline indices for xbin i
C--   nx    (in)    number of x-grid points
C--   ierr  (out)   0 = all OK
C--                 1 = x-index out of range

      implicit double precision (a-h,o-z)
    
      dimension ixft(*),imi(*),ima(*)

      ierr = 0
      
      do i = 1,nx
        imi(i) = nt+1
        ima(i) = 0
      enddo

      do i = 1,nt-k
C--     what is the start of the support?
        ix1 = ixft(i)
        if(ix1.gt.nx) then
          ierr = 1
          stop 'sqcSrange: ix1 out of range ---> STOP'
        endif
C--     what is the end of the support?
        ix2 = ixft(i+k)
        if(ix2.gt.nx) then
          ierr = 1
          stop 'sqcSrange: ix2 out of range ---> STOP'
        endif
C--     so, spline i covers bins ix1 to ix2-1
        do j = ix1,ix2-1 
          imi(j) = min(imi(j),i)
          ima(j) = max(ima(j),i)
        enddo
      enddo

      return
      end

C     =======================================
      subroutine sqcSplCat(k,t,ic,nt,nc,ierr)
C     =======================================

C--   Set-up a Catalog by weeding-out identical Bsplines.
C--   Bsplines i and j are identical if shifted [tau(i),...,tau(i+k)]-tau(i)
C--   equals [tau(j),...,tau(j+k)]-tau(j). Different splines make a new entry
C--   in the catalog. For each Bspline (i), ic(i) gives the entry in the
C--   catalog.
C--
C--   k     (in)    spline order
C--   t     (in)    tau-grid [from sqcGetTau]
C--   ic    (out)   ic(i) is catalog index of bspline i  
C--   nt    (in)    number of points in the tau-grid (and of ic)
C--   nc    (out)   number of entries in the catalog (number of
C--                 Bsplines with a different shape)
C--   ierr  (out)   always 0 (this routine cannot fail)

      implicit double precision (a-h,o-z)

      logical lqcAcomp
      dimension t(*),ic(*)

      data epsi /1.D-10/

      ierr = 0

      nc    = 1
      ic(1) = nc

      do i = 2,nt-k
        icount = 0
        do j = 1,k+1
C--       current shifted tau
          tau1 = t(i+j-1)-t(i)
C--       previous shifted tau
          tau2 = t(i+j-2)-t(i-1)
C--       same?
          if(lqcAcomp(tau1,tau2,epsi)) icount = icount+1
        enddo
        if(icount.eq.k+1) then
C--       no new entry in the catalog
          ic(i) = nc
        else
C--       new entry in the catalog
          nc    = nc+1
          ic(i) = nc
        endif
      enddo

      return
      end

C     ================================================================
      subroutine sqcBsplin(iord,x,tau,ntau,bsder,iom,ntm,imi,ima,ierr)
C     ================================================================

C--   Calculate all non-zero normalized bsplines at x and also the
C--   derivatives (up to order iord-1) and store results in bsder.
C--   Called for each x-grid point in s/r sqcFilCat.
C--
C--   Iterative algorithm taken from L.L. Shumaker,
C--   'Spline functions: basic theory', Krieger Publishing Company
C--   (1993), p192, algorithm 5.5.
C--
C--   iord  (in)  order of the spline (2=linear, 3=quadratic etc)   
C--   x     (in)  value of x
C--   tau   (in)  array of tau-nodes
C--   ntau  (in)  number of tau-nodes 
C--   bsder (out) bsder(j,i) j=1,...,iord is the derivative index: 
C--                          1=value, 2=1st derivative etc;
C--                          i=1,...,ntau-iord is the bspline index
C--   iom   (in)  first dim of bsder in the calling routine (.ge.iord)
C--   ntm   (in)  second dimension of bsder (.ge.ntau)
C--   imi   (out) index of the first nonzero spline at x
C--   ima   (out) index of the last  nonzero spline at x
C--   ierr  (out) 0 all OK
C--               1 iom (first dim of bsder) too small
C--               2 ntm (2nd   dim of bsder) too small
C--               3 x out of range (all bsplines set to zero)             

      implicit double precision (a-h,o-z)

      dimension tau(*),bsder(iom,ntm)

C--   Check dimensions
      imi = 0
      ima = 0
      if(iord.gt.iom) then
        stop 'sqcBsplin: first dim of bsder too small ---> STOP'
      endif
      if(ntau.gt.ntm) then
        stop 'sqcBsplin: sedond dim of bsder too small ---> STOP'
      endif

C--   Initialization
      ierr = 0
      do j = 1,iom
        do i = 1,ntm
          bsder(j,i) = 0.D0
        enddo
      enddo

C--   Find bin
      ix = 0
      do i = ntau-1,1,-1
        if(x.ge.tau(i).and.x.lt.tau(i+1)) then
          ix = i
          goto 10
        endif
      enddo
 10   continue
      if(ix.eq.0) then
        ierr = 3
        return
      endif

C--   First order normalized bspline
      bsder(1,ix) = 1.D0
      if(iord.eq.1) return

C--   Get iord-1 unnormalized bsplines by iteration
      bsder(1,ix) = 1.D0/(tau(ix+1)-tau(ix))
      do k = 2,iord-1
        ix1 = max(ix-k+1,1)
        do i = ix1,ix
          d  = tau(i+k)-tau(i)
          t  = (x-tau(i))/d
C--       The derivatives.....
          do j = k,2,-1
            bsder(j,i) = (k-1)*
     +                       (bsder(j-1,i)-bsder(j-1,i+1))/d
          enddo
C--       The bsplines.....
          bsder(1,i) = t*bsder(1,i) + (1.D0-t)*bsder(1,i+1)
        enddo
      enddo

C--   Final iteration: get iord normalized bsplines
      imi  = max(ix-iord+1,1)
      ima  = min(ix,ntau-iord)
      do i = imi,ima
C--     The derivatives.....
        do j = iord,2,-1
          bsder(j,i) = (iord-1)*
     +                     (bsder(j-1,i)-bsder(j-1,i+1))
        enddo
C--     The bsplines.....
        bsder(1,i) = 
     +  (x-tau(i))*bsder(1,i) + (tau(i+iord)-x)*bsder(1,i+1)
      enddo

      return
      end

C     =================================================================
      subroutine sqcFilCat(k,xx,jtfx,nx,tt,ic,nt,der,cat,mk,mt,nc,ierr)
C     =================================================================

C--   Fill spline catalog (3-dim array cat) defined as follows: 
C--   Let x fall in region i of the support of Bspline j, 
C--   for example for order k = 3:
C--
C--             |<-------- support of B_j  ------->|
C--             |                                  |
C--             t_j ---- t_{j+1} ---- t_{j+2} ---- t_{j+3}
C--                  ^             ^           ^
C--                  |             |           |
C--                 i=1           i=2         i=3
C--
C--   Let n be the address of B_j in the catalog
C--   Then, for x in the region i:
C--
C--         B_j(x) = Sum_{m=1}^k cat(m,i,n) (x-t_{j+i-1})^{m-1}
C--
C--   k    (in)  spline order
C--   xx   (in)  x grid
C--   jtfx (in)  array which maps ix onto jtau (from s/r sqcGetTau)
C--   nx   (in)  number of x-grid points
C--   tt   (in)  tau-grid (as obtained by s/r sqcGetTau)
C--   ic   (in)  index array (as obtained by s/r sqcSplcat)
C--   nt   (in)  number of points in tt and ic
C--   der  (out) work array dimensioned der(mk,mt) in the calling routine
C--   cat  (out) spline catalog dimensioned cat(mk,mk,mt)
C--   mk   (in)  first dim of der and first two dims of cat
C--   mt   (in)  last dim of der and cat
C--   nc   (out) number of entries in the catalog
C--   ierr (out) 0 all OK
C--              1 indexing error --> should never occur....

      implicit double precision (a-h,o-z)

      dimension xx(*),jtfx(*),tt(*),ic(*),der(mk,mt),cat(mk,mk,mt)

*mb   xx is apparently a dummy variable --> get ridofit at some point
      dum = xx(1)   !get rid of compiler warning in poormans solution 

C--   Initialization
      ierr = 0
      do kk = 1,mt
        do jj = 1,mk
          do ii = 1,mk
            cat(ii,jj,kk) = 0.D0
          enddo
        enddo
      enddo

C--   Get values and derivatives of all non-zero bsplines at each xi
      nc = 0
      do ix = 1,nx-1
        jt = jtfx(ix)
        x  = tt(jt)
        call sqcBsplin(k,x,tt,nt,der,mk,mt,imi,ima,ierr)
        if(ierr.ne.0) then
          stop 'sqcFilCat: invalid call to sqcBsplin ---> STOP'
        endif
        do is = imi,ima
C--       Find region within support of the spline
          js  = jt-is+1
          if(js.le.0.or.js.gt.k) then
            ierr = 1
            stop 'sqcFilCat: indexing error ---> STOP'
          endif
          id = ic(is)
          nc = max(nc,id)
C--       Value
          cat(1,js,id) = der(1,is)
C--       Derivatives divided by n!
          fac = 1.D0
          do j = 2,k
            cat(j,js,id) = der(j,is)/fac
            fac          = fac*j
          enddo
        enddo
      enddo

      return
      end

C     ==================================================
      double precision function 
     +  dqcBsplix(k,is,ix,jtfx,mi,ma,nx,tt,ic,cat,mk,mt)
C     ==================================================

C--   Calculate Bspline 'is' at ix
C--
C--   k     (in) spline order
C--   is    (in) Bspline index (can be anything)
C--   ix    (in) x-bin index should be set by iqcBGetIx before
C--   jtfx  (in) array which maps ix onto jtau       (from s/r sqcGetTau)
C--   mi    (in) array with lower spline indices     (from s/r sqcSrange)
C--   ma    (in) array with upper spline indices     (from s/r sqcSrange)
C--   nx    (in) number of points in the x-grid
C--   tt    (in) t-grid
C--   ic    (in) pointer array to the catalog        (from s/r sqcSplCat)
C--   cat   (in) catalog cat(mk,mk,mt)               (from s/r sqcFilCat)
C--   mk    (in) first two dims of cat as declared in the calling routine 
C--   mt    (in) last      dim  of cat as declared in the calling routine 

      implicit double precision (a-h,o-z)

      dimension jtfx(*),mi(*),ma(*),tt(*),ic(*),cat(mk,mk,mt)

*mb   nx and tt are dummies --> getrid of these at some point
      idum = nx       !get rid of compiler warning
      dum  = tt(1)    !get rid of compiler warning

      dqcBsplix = 0.D0
C--   ix in range?
      if(ix.eq.0) return
C--   Check if spline index is in range
      if((is.lt.mi(ix)).or.(is.gt.ma(ix))) return
C--   Find entry in the spline catalog
      ii = ic(is)
C--   Find region within support of the spline
      jt  = jtfx(ix)
      js  = jt-is+1
C--   Check if in range (debug code if not ...)
      if(js.le.0.or.js.gt.k) stop 'Index error in dqcBsplix ---> STOP'
C--   Pick constant term of the Taylor expansion
      dqcBsplix = cat(1,js,ii)

      return
      end

C     ====================================================
      double precision function 
     +  dqcBsplxx(k,is,x,ix,jtfx,mi,ma,nx,tt,ic,cat,mk,mt)
C     ====================================================

C--   Calculate Bspline 'is' at x using piecewise polynomial expansion
C--
C--   k     (in) spline order
C--   is    (in) Bspline index (can be anything)
C--   x     (in) x value       
C--   ix    (in) x-bin index should be set by iqcBGetIx before
C--   jtfx  (in) array which maps ix onto jtau       (from s/r sqcGetTau)
C--   mi    (in) array with lower spline indices     (from s/r sqcSrange)
C--   ma    (in) array with upper spline indices     (from s/r sqcSrange)
C--   nx    (in) number of points in the x-grid
C--   tt    (in) t-grid
C--   ic    (in) pointer array to the catalog        (from s/r sqcSplCat)
C--   cat   (in) catalog cat(mk,mk,mt)               (from s/r sqcFilCat)
C--   mk    (in) first two dims of cat as declared in the calling routine 
C--   mt    (in) last      dim  of cat as declared in the calling routine 

      implicit double precision (a-h,o-z)

      dimension jtfx(*),mi(*),ma(*),tt(*),ic(*),cat(mk,mk,mt)

*mb   Dummy argument nx
      idum = nx

      dqcBsplxx = 0.D0
C--   x in range?
      if(ix.eq.0) return
C--   Check if spline index is in range
      if((is.lt.mi(ix)).or.(is.gt.ma(ix))) return
C--   Find entry in the spline catalog
      ii = ic(is)
C--   Find region within support of the spline
      jt  = jtfx(ix)
      js  = jt-is+1
C--   Check if in range (debug code if not ...)
      if(js.le.0.or.js.gt.k) stop 'Index error in dqcBsplxx ---> STOP'
C--   Interpolation variable
      r = x-tt(jt)

C--   Calculate Taylor expansion using Horners scheme
      spl = cat(k,js,ii)
      do m = k-1,1,-1
        spl = spl*r + cat(m,js,ii)
      enddo
      dqcBsplxx = spl 

      return
      end

C     ====================================================
      double precision function 
     +  dqcDsplxx(k,is,x,ix,jtfx,mi,ma,nx,tt,ic,cat,mk,mt)
C     ====================================================

C--   Calculate dBspline_is/dx at x using piecewise polynomial expansion
C--
C--   k     (in) spline order
C--   is    (in) Bspline index (can be anything)
C--   x     (in) x value       
C--   ix    (in) x-bin index should be set by iqcBGetIx before
C--   jtfx  (in) array which maps ix onto jtau       (from s/r sqcGetTau)
C--   mi    (in) array with lower spline indices     (from s/r sqcSrange)
C--   ma    (in) array with upper spline indices     (from s/r sqcSrange)
C--   nx    (in) number of points in the x-grid
C--   tt    (in) t-grid
C--   ic    (in) pointer array to the catalog        (from s/r sqcSplCat)
C--   cat   (in) catalog cat(mk,mk,mt)               (from s/r sqcFilCat)
C--   mk    (in) first two dims of cat as declared in the calling routine 
C--   mt    (in) last      dim  of cat as declared in the calling routine 

      implicit double precision (a-h,o-z)

      dimension jtfx(*),mi(*),ma(*),tt(*),ic(*),cat(mk,mk,mt)

*mb   Dummy argument nx
      idum = nx

      dqcDsplxx = 0.D0
C--   x in range?
      if(ix.eq.0) return
C--   Check if spline index is in range
      if((is.lt.mi(ix)).or.(is.gt.ma(ix))) return
C--   Find entry in the spline catalog
      ii = ic(is)
C--   Find region within support of the spline
      jt  = jtfx(ix)
      js  = jt-is+1
C--   Check if in range (debug code if not ...)
      if(js.le.0.or.js.gt.k) stop 'Index error in dqcBsplxx ---> STOP'
C--   Interpolation variable
      r = x-tt(jt)

C--   Calculate Taylor expansion using Horners scheme
      spl = (k-1)*cat(k,js,ii)
      do m = k-2,1,-1
        spl = spl*r + m*cat(m+1,js,ii)
      enddo
      dqcDsplxx = spl 

      return
      end

C     ===================================
      integer function iqcBGetIx(x,xx,nx)
C     ===================================

C--   Search for index ix with xx(ix) .le. x
C--   If, within tolerance epsi, x = xx(nx) then ix = nx-1
C--   If x is out of range then ix is set to zero
C--
C--   x     (in) x value (can be anything)
C--   xx    (in) array with x-grid points
C--   nx    (in) number of points in the x-grid

      implicit double precision (a-h,o-z)
      logical lqcAcomp

      dimension xx(*)

      data epsi /1.D-10/

C--   Search for bin ix
      iqcBGetIx = 0
      do i = 1,nx-1
        if((x.ge.xx(i)).and.(x.lt.xx(i+1))) then
          iqcBGetIx = i
          goto 10
        endif
      enddo
 10   continue
C--   x not in range but check if x = xx(nx)
      if(iqcBGetIx.eq.0) then 
        if(lqcAcomp(x,xx(nx),epsi)) then
          iqcBGetIx = nx-1
        else
          return
        endif
      endif

      return
      end
