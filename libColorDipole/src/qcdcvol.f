
C--   This is the file qcdcvol.f containing convolution and 
C--   interpolation routines used by the convolution engine
C--
C--   double precision function dqcFcrossC(w,id,idi,idt,jy,it)
C--   double precision function dqcFcrossK(w,idw,jset,idf,iy,it)
C--   double precision function dqcFcrossF(w,idw,jset,ida,idb,iy,it)
C--   subroutine sqcStInterp(stfun,xlist,qlist,flist,n,ierr)
C--   subroutine sqcStIntMpt(stfun,xlist,qlist,flist,n,ierr)


C==   ==============================================================
C==   Convolution ==================================================
C==   ==============================================================

C     ========================================================
      double precision function dqcFcrossC(w,id,idi,idt,jy,it)
C     ========================================================

C--   Dummy routine

      implicit double precision (a-h,o-z)

      dimension w(*)

      dum        = w(1)
      idum       = id
      idum       = idi
      idum       = idt
      idum       = jy
      idum       = it
      dqcFcrossC = 0.D0

      return
      end
      
C     ==========================================================
      double precision function dqcFcrossK(w,idw,jset,idf,iy,it)
C     ========================================================== 

C--   Calculate convolution F cross K

C--   w     (in)  :  store previously filled with weight tables
C--   idw   (in)  :  table with weights K
C--   jset  (in)  :  pdf set identifier
C--   idf   (in)  :  identifier of F
C--   iy,it (in)  :  grid point           
      
      implicit double precision(a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      
      dimension w(*),coef(mxx0)
      
      iz = izfit2(it)
      nf = nffiz2(iz)

C--   Get spline coefficients
      call sqcGetSplA(jset,idf,iy,iz,ig,iyg,coef)
      
C--   Weight table base address
      iwk = iqcWaddr(w,1,it,nf,ig,idw)-1

C--   Convolute               
      fxk = 0.D0
      do j = 1,iyg
        fxk = fxk + coef(j)*w(iwk+iyg-j+1)
      enddo
                      
      dqcFcrossK = fxk
      
      return
      end

C     ==============================================================
      double precision function dqcFcrossF(w,idw,jset,ida,idb,iy,it)
C     ============================================================== 

C--   Calculate convolution f_a cross f_b

C--   w     (in)  :  store previously filled with weight tables
C--   idw   (in)  :  table with fxf weights
C--   jset  (in)  :  pdf set identifier
C--   ida   (in)  :  identifier of pdf f_a
C--   idb   (in)  :  identifier of pdf f_b 
C--   iy,it (in)  :  grid point           
      
      implicit double precision(a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      
      dimension w(*),coefa(mxx0),coefb(mxx0)
      
      iz = izfit2(it)
      nf = nffiz2(iz)

C--   Get spline coefficients
      call sqcGetSplA(jset,ida,iy,iz,ig,iyg,coefa)
      call sqcGetSplA(jset,idb,iy,iz,ig,iyg,coefb)
      
C--   Weight table base address
      iwx = iqcWaddr(w,1,it,nf,ig,idw)-1

C--   Convolute               
      fxf = 0.D0
      do j = 1,iyg
        Aj = coefa(j)
        do k = 1,iyg-j+1
          Bk  = coefb(k) 
          fxf = fxf + Aj*Bk*w(iwx+iyg-j-k+2)
        enddo
      enddo
                      
      dqcFcrossF = fxf
      
      return
      end

C==   ==============================================================
C==   Interpolation ================================================
C==   ==============================================================

C     ======================================================
      subroutine sqcStInterp(stfun,xlist,qlist,flist,n,ierr)
C     ======================================================

C--   Interpolation of stfun(ix,iq) for arbitrary large n
C--   This routine does not suffer from the mpt0 limit (by buffering)
C--
C--   stfun      (in)  :     Function of (ix,iq), declared external
C--   xlist(n)   (in)  :     List of x-values
C--   qlist(n)   (in)  :     List of qmu2 values
C--   flist(n)   (out) :     List of interpolated structure functions
C--   n          (in)  :     Number of items in the list
C--   ierr       (out) :     -1   Count mpt0 exceeded (never occurs)
C--                           0   OK
C--                           1   One or more x,q outside grid 
C--

      implicit double precision(a-h,o-z)

      include 'qcdnum.inc'
      
      external stfun
      
      dimension xx(mpt0),qq(mpt0)
      dimension xlist(*), qlist(*), flist(*)
      
C--   Fill output array flist in batches of mpt0 words
      ipt = 0
      jj  = 0      
      do i = 1,n
        ipt     = ipt+1
        xx(ipt) = xlist(i)
        qq(ipt) = qlist(i)
        if(ipt.eq.mpt0) then
          call sqcStIntMpt(stfun,xx,qq,flist(jj*mpt0+1),mpt0,ierr)
          ipt = 0
          jj  = jj+1
        endif
      enddo
C--   Flush remaining ipt points
      if(ipt.ne.0) then
        call sqcStIntMpt(stfun,xx,qq,flist(jj*mpt0+1),ipt,ierr)
      endif

      return
      end

C     ======================================================
      subroutine sqcStIntMpt(stfun,xlist,qlist,flist,n,ierr)
C     ======================================================

C--   Interpolation of stfun(ix,iq)
C--
C--   stfun      (in)  :     Function of (ix,iq), declared external
C--   xlist(n)   (in)  :     List of x-values
C--   qlist(n)   (in)  :     List of qmu2 values
C--   flist(n)   (out) :     List of interpolated structure functions
C--   n          (in)  :     Number of items in the list
C--   ierr       (out) :     -1   Count mpt0 exceeded
C--                           0   OK
C--                           1   One or more x,q outside grid  

      implicit double precision(a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      external stfun

      dimension xlist(*), qlist(*), flist(*)
      dimension yin(mpt0),tin(mpt0),yout(mpt0),tout(mpt0),infout(mpt0)
      dimension iy1(mpt0),iy2(mpt0),iz1(mpt0),iz2(mpt0),fout(mpt0)
      logical   mark(0:mxx0,0:mqq0+7)

C--   Check
      if(n.gt.mpt0) then
         ierr = -1
         return
      endif
C--   Convert to y-t
      do i = 1,n
        yin(i) = -log(xlist(i))
        tin(i) =  log(qlist(i))
      enddo
C--   Weed points outside grid or cuts
      call sqcWeedIt(yin,tin,n,yout,tout,infout,nout)
C--   Not all points inside grid
      ierr = 0
      if(nout.ne.n) ierr = 1                       
C--   Mark grid points
      margin = 1
      call sqcMarkit(mark,yout,tout,margin,iy1,iy2,iz1,iz2,nout)
C--   Find pdf set for temp buffer
      jset = 0
      do j = 1,9
        if(mxord7(j).ne.0) jset = j 
      enddo 
      if(jset.eq.0) stop 'sqcStInterp: no temp buffer available'      
      idt = -2
C--   Store structure function in scratch buffer
      do it = 1,ntt2
        iz = izfit2(it)
        ia = iqcPdfIjkl(1,iz,idt,jset)-1
        do iy = 1,nyy2(0)
          ia = ia+1
          if(mark(iy,iz)) then
            ix        = nyy2(0)+1-iy
            stor7(ia) = stfun(ix,it)
          endif
        enddo    
      enddo
C--   Interpolate
      call sqcIntLst(jset,idt,iy1,iy2,iz1,iz2,yout,tout,fout,nout)
C--   Preset
      do i = 1,n
        flist(i) = qnull6
      enddo
C--   Return interpolated result
      do i = 1,nout
        flist(infout(i)) = fout(i)
      enddo                

      return
      end
