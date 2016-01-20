
C--   This is the file usrgrd.f containing the qcdnum grid routines

C--   subroutine gxmake(xmi,iwt,n,nxin,nxout,iosp)
C--   integer function ixfrmx(x)
C--   logical function xxatix(x,ix)
C--   double precision function xfrmix(ix)
C--   subroutine gxcopy(array,n,nx)
C--   subroutine gqmake(qq,ww,n,nqin,nqout)
C--   integer function iqfrmq(q)
C--   logical function qqatiq(q,iq)
C--   double precision function qfrmiq(iq)
C--   subroutine gqcopy(array,n,nq)
C--   subroutine grpars(nx,xmi,xma,nq,qmi,qma,iosp)
C--
C--   subroutine setcut(xmi,qmi,qma,roots)
C--   subroutine getcut(xmi,qmi,qma,roots)
C--   logical function lpassc(x,qmu2,ifail,jchk)
      
C==   ===============================================================
C==   x-Grid routines ===============================================
C==   ===============================================================

C     ============================================
      subroutine gxmake(xmi,iwx,n,nxin,nxout,iosp)
C     ============================================

C--   Define logarithmic x-grid
C--
C--   xmi(n)   = (in)  list of lowest gridpoints for each subgrid
C--   iwx(n)   = (in)  list of point density weights for each subgrid
C--   n        = (in)  number of subgrids
C--   nxin     = (in)  requested number of grid points
C--   nxout    = (out) generated number of grid points
C--   iosp     = (in)  spline interpolation 2=lin, 3=quad

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qpdfs7.inc'
      include 'qmaps8.inc'
 
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension xmi(*)   ,iwx(*)
      dimension yma(mxg0),iwy(mxg0)

      character*80 subnam
      data subnam /'GXMAKE ( XMI, IWT, NGR, NXIN, NXOUT, IORD )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check user input
C--   1-check if spline order OK
      call sqcIlele(subnam,'IORD',2,iosp,3,
     + 'Only linear (2) or quadratic (3) interpolation is allowed')
C--   2-Check if number of subgrids OK
      call sqcIlele(subnam,'NGR',1,n,mxg0,
     + 'Remark: you can increase mxg0 in qcdnum.inc and recompile')
C--   3-Check if number of gridpoints OK 
      call sqcIlele(subnam,'NXIN', max(iosp,n), nxin, mxx0-11,
     + 'Remark: you can increase mxx0 in qcdnum.inc and recompile')
C--   4-Check if all xmi(i) are in range
      do i = 1,n
        call sqcDltlt(subnam,'XMI(i)',0.D0,xmi(i),1.D0,
     +  'At least one of the XMI(i) outside allowed range')
      enddo
C--   5-Check if all xmi(i) are in ascending order
      if(n.ge.2) then
        do i = 2,n
          if(xmi(i).le.xmi(i-1)) call sqcErrMsg(subnam,
     +    'XMI(i) not in ascending order')
        enddo
      endif
C--   6-Check that all weights are ascending integer multiples
      if(iwx(1).le.0) call sqcErrMsg(subnam, 
     +       'Zero or negative weight encountered')
      do i = 2,n
        if(iwx(i).le.0) call sqcErrMsg(subnam, 
     +     'Zero or negative weight encountered')
        irat = iwx(i)/iwx(i-1)
        if(iwx(i).ne.irat*iwx(i-1)) call sqcErrMsg(subnam, 
     +     'Weights are not ascending integer multiples')
      enddo

C--   Transform x-grid to y-grid
      do i = 1,n
        j = n+1-i
        yma(j) = -log(xmi(i))
        iwy(j) =  iwx(i)
      enddo
C--   Do the work
      call sqcGryDef(yma,iwy,n,nxin,nxout,iosp)

C--   Invalidate weight tables
      Lwtini8 = .false.
C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
      
C     ==========================
      integer function ixfrmx(x)
C     ==========================

C--   Gives binnumber ix, given x
C--   ix = 0 if x < xmin, x >= 1 or if xgrid not defined.

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
 
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      character*80 subnam
      data subnam /'IXFRMX ( X )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

      ixfrmx = 0
C--   Check status bits
      call sqcChekit(1,ichk,jbit)
      if(jbit.ne.0)                return
C--   Check user input
      if(x.le.0.D0 .or. x.ge.1.D0) return
C--   Go...
      y         = -log(x)
      iy        = iqcFindIy(y)
      if(iqcYhitIy(y,iy).eq.1) then
        ixfrmx = nyy2(0) + 1 - iy
      else
        ixfrmx = nyy2(0) - iy
      endif
 
      return
      end
      
C     =============================
      logical function xxatix(x,ix)
C     =============================

C--   True if x is at gridpoint ix.
C--   False of not at gridpoint, x or ix out of range or no xgrid defined

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
 
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      character*80 subnam
      data subnam /'XXATIX ( X, IX )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

      xxatix = .false.
C--   Check status bits
      call sqcChekit(1,ichk,jbit)
      if(jbit.ne.0)                  return
C--   Check user input
      ymax = ygrid2(nyy2(0))
      xmin = exp(-ymax)
      if(x.lt.xmin .or. x.ge.1.D0)   return
      if(ix.lt.1 .or. ix.gt.nyy2(0)) return
C--   Go...
      ihit  = iqcYhitIy(-log(x),nyy2(0)+1-ix)
      if(ihit.eq.1) then
        xxatix = .true.
      else
        xxatix = .false.     
      endif

      return
      end
      
C     ====================================
      double precision function xfrmix(ix)
C     ====================================

C--   Get value of x-grid point ix.
C--   Returns 0.D0 if ix out of range [1,nyy2] or xgrid not defined

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
 
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      character*80 subnam
      data subnam /'XFRMIX ( IX )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

      xfrmix = 0.D0
C--   Check status bits
      call sqcChekit(1,ichk,jbit)
      if(jbit.ne.0)               return
C--   Check user input
      if(ix.lt.1 .or. ix.gt.nyy2(0)) return
C--   Go...     
      iy     = nyy2(0) + 1 - ix
      yy     = ygrid2(iy)
      xfrmix = exp(-yy)
      
      return
      end

C     =============================
      subroutine gxcopy(array,n,nx)
C     =============================

C--   Copy x grid to local array
C--
C--   Input      n = dimension of array as declared in the calling routine
C--   Output array = target array declared in the calling routine
C--             nx = number of x grid point copied to array

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension array(*)

      character*80 subnam
      data subnam /'GXCOPY ( XARRAY, N, NX )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check user input
      call sqcIlele(subnam,'N',nyy2(0),n,100000,
     &              'XARRAY not large enough to contain x-grid')

      nx = nyy2(0)     
      do ix = 1,nx
        iy         = nyy2(0) + 1 - ix
        yy         = ygrid2(iy)
        array(ix)  = exp(-yy)
      enddo
      
      return
      end

C==   ===============================================================
C==   Q2-Grid routines ==============================================
C==   ===============================================================

C     =====================================
      subroutine gqmake(qq,ww,n,nqin,nqout)
C     =====================================

C--   Define logarithmic Q2-grid
C--
C--   qq(n)   (in)  List of Q2 values. qq(1) and qq(n) are the grid limits
C--   ww(n)   (in)  Weights: generated point density beween qq(i) and qq(i+1)
C--                 will be proportional to ww(i)   
C--   n       (in)  Number of points in qq and ww (>=2)
C--   nqin    (in)  Requested number of grid points. If <=n then qq will be
C--                 copied to the internal grid. If < 0 loglog spacing
C--   nqout   (out) Number of generated grid points

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qpdfs7.inc'
      include 'qmaps8.inc'
 
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension qq(*),ww(*),qlog(mqq0)
      
      logical loglog

      character*80 subnam
      data subnam /'GQMAKE ( QARR, WGT, N, NQIN, NQOUT )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check user input
      call sqcIlele(subnam,'N',2,n,mqq0-5,
     + 'Remark: You can increase mqq0 in qcdnum.inc and recompile')
      call sqcIlele(subnam,'NQIN',n,nqin,mqq0-10,
     + 'Remark: You can increase mqq0 in qcdnum.inc and recompile')
C--   Q2 should be larger than 0.1
      call sqcDltlt(subnam,'QARR(1)',qlimd6,qq(1),qlimu6,
     + 'Remark: the allowed range can be changed by a call to SETVAL')
C--   Spacing should be larger than 0.01 GeV2 --> force ascending order
      do i = 2,n
        if( (qq(i-1)+1.D-2) .ge. qq(i) ) then 
          call sqcErrMsg(subnam, 
     +   'QARR(i) not ascending or spaced by less than 0.01 GeV2')
        endif
        call sqcDltlt(subnam,'QARR(i)',qlimd6,qq(i),qlimu6,
     + 'Remark: these Q2 limits can be changed by a call to SETVAL')
      enddo
C--   Weights should be between 0.1 and 10
      do i = 1,n-1
        call sqcDlele(subnam,'WGT(i)',0.1D0,ww(i),10.D0,
     +        'Weights should be in a reasonable range [0.1,10]')
      enddo
C--   Safety margin for max number of points
      call sqcIlele(subnam,'NQIN',5-mqq0,nqin,mqq0-5,
     + 'Remark: You can increase mqq0 in qcdnum.inc and recompile')

C--   Do the work
      do i = 1,n
        qlog(i) = log(qq(i))
      enddo
C--   Log or loglog, thats the question
C--   Loglog is disabled in the check above since nqin < 0 is not allowed
C--   It is disabled because the loglog option is not fully checked if OK
      if(nqin.gt.0) then       
        nq2    =  nqin
        loglog = .false.
      else
        nq2    = -nqin  
        loglog = .true.
      endif  
      call sqcGrTdef(qlog,ww,n,nq2,loglog,jerr)
      if(jerr.ne.0) then 
         write(lunerr1,*) 'sqcGrTdef jerr = ',jerr,' ---> STOP'
         stop
      endif
      nqout = nq2
C--   Safety margin for max number of points
      call sqcIlele(subnam,'NQOUT',2,nqout,mqq0-5,
     & 'Remark: You can increase mqq0 in qcdnum.inc and recompile')

C--   Update status bits
      call sqcSetflg(iset,idel,0)
C--   Invalidate weight tables
      Lwtini8 = .false.      
C--   Invalidate astable, flavor map and tsplines
      Lastab8 = .false.
      Lnfmap8 = .false.
    
      return
      end
      
C     ==========================
      integer function iqfrmq(q)
C     ==========================

C--   Get index iq of first gridpoint below Q2.
C--   iq = 0 id q < qmin, q > qmax or no Q2 grid defined

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
 
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      logical lqcAxlty, lqcAxgty

      character*80 subnam
      data subnam /'IQFRMQ ( Q2 )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

      iqfrmq = 0
C--   Check status bits
      call sqcChekit(1,ichk,jbit)
      if(jbit.ne.0)                         return
C--   Check user input
      if(q.le.0.D0)                         return
      t = log(q)
      if(lqcAxlty(t,tgrid2(1)   ,aepsi6))   return
      if(lqcAxgty(t,tgrid2(ntt2),aepsi6))   return
C--   Get binnumber
      iqfrmq = iqcItfrmt(t)  

      return
      end
      
C     =============================
      logical function qqatiq(q,iq)
C     =============================

C--   True if q at grid point iq
C--   False if q not at iq, q or iq not in range, or if grid not defined

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
 
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      logical lqcAxlty, lqcAxgty

      character*80 subnam
      data subnam /'QQATIQ ( Q2, IQ )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

      qqatiq = .false.
C--   Check status bits
      call sqcChekit(1,ichk,jbit)
      if(jbit.ne.0)                         return
C--   Check user input
      if(q.le.0.D0)                         return
      t = log(q)
      if(lqcAxlty(t,tgrid2(1)   ,aepsi6))   return
      if(lqcAxgty(t,tgrid2(ntt2),aepsi6))   return
      if(iq.lt.1 .or. iq.gt.ntt2)           return
C--   Go...
      ihit = iqcThitit(t,iq)
      if(ihit.eq.1) then
        qqatiq = .true.
      else
        qqatiq = .false.     
      endif

      return
      end
      
C     ====================================
      double precision function qfrmiq(iq)
C     ====================================

C--   Get value of Q2-grid point iq.
C--   Returns 0.D0 if iq not in range [1,ntt2] or if Q2 grid not defined

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
 
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      character*80 subnam
      data subnam /'QFRMIQ ( IQ )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

      qfrmiq = 0.D0
C--   Check status bits
      call sqcChekit(1,ichk,jbit)
      if(jbit.ne.0)                  return
C--   Check user input
      if(iq.lt.1 .or. iq.gt.ntt2)    return
C--   Go...     
      qfrmiq = exp(tgrid2(iq))
      
      return
      end
      
C     =============================
      subroutine gqcopy(array,n,nq)
C     =============================

C--   Copy q2 grid to local array
C--
C--   Input      n = dimension of array as declared in the calling routine
C--   Output array = target array declared in the calling routine
C--             nq = number of q2 grid point copied to array

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
 
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension array(*)

      character*80 subnam
      data subnam /'GQCOPY ( QARRAY, N, NQ )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check user input
      call sqcIlele(subnam,'N',ntt2,n,100000,
     &              'QARRAY not large enough to contain Q2-grid')

      nq = ntt2
      do iq = 1,nq     
        array(iq)= exp(tgrid2(iq))
      enddo
      
      return
      stop

      end

C==   ===============================================================
C==   General access to the x-Q2 grid ===============================
C==   ===============================================================

C     =============================================
      subroutine grpars(nx,xmi,xma,nq,qmi,qma,iosp)
C     =============================================

C--   Returns the number of points and the limits of the x-Q2 grid.

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
 
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      character*80 subnam
      data subnam /'GRPARS ( NX, XMI, XMA, NQ, QMI, QMA, IORD )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

      nx   = nyy2(0)
      xmi  = exp(-ymax2(0))
      xma  = 1.D0
      nq   = ntt2
      qmi  = exp(tgrid2(1))
      qma  = exp(tgrid2(ntt2))
      iosp = ioy2
      
      return
      end
      
C==   ===============================================================
C==   Evolution cuts ================================================
C==   ===============================================================

C     ======================================
      subroutine setcut(xmi,q2mi,q2ma,roots)
C     ======================================

C--   Set kinematic cuts

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      character*80 subnam
      data subnam /'SETCUT ( XMI, Q2MI, Q2MA, DUMMY )'/
      
C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
      
C--   Do the work
*mb      call sqcSetCut(xmi,q2mi,q2ma,roots)
      dummy = roots  !avoid compiler warning
      call sqcSetCut(xmi,q2mi,q2ma,0.D0)        
      
C--   Update status bits
      call sqcSetflg(iset,idel,0)      

      return
      end

C     ======================================
      subroutine getcut(xmi,q2mi,q2ma,roots)
C     ======================================

C--   Set kinematic cuts

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      character*80 subnam
      data subnam /'GETCUT ( XMI, Q2MI, Q2MA, DUMMY )'/
      
C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
      
C--   Do the work
      xmi   = xminc2
      q2mi  = qminc2
      q2ma  = qmaxc2
      roots = 0.D0
*mb      roots = sqrt(smaxc2)     

      return
      end
      
C     ==========================================      
      logical function lpassc(x,qmu2,ifail,jchk)
C     ========================================== 

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qsnam3.inc'
      
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      character*80 subnam
      data subnam /'LPASSC ( X, QMU2, IFAIL, ICHK )'/
      
      logical lqcInside
      
C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
      
C--   Do the work
      lpassc = lqcInside(x,qmu2,ifail)
      
      if(lpassc .or. jchk.ne.1) return
      
      len = imb_lenoc(usrnam3)
      if(len.ne.0) then
        call sqcCutMsg(usrnam3,x,qmu2,ifail,1)
      else  
        call sqcCutMsg(subnam,x,qmu2,ifail,1)
      endif  
            
      return
      end

      
     
