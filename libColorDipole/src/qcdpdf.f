
C--   This is the file qcdpdf.f containing handling of pdfs

C--   subroutine sqcPdfMat
C--   subroutine sqcPdIdef(tmatpq,ierr)
C--   subroutine sqcAllInp(itype,func)
C--   subroutine sqcEfrmP(pval,eval)
C--   double precision function dqcEifrmP(i,pval)
C--   subroutine sqcPdfInp(subr,jset,del,epsi,nwlast,ierr)
C--   double precision function dqcXSplne(jset,id,y,it)
C--   double precision function dqcEpmyt(jset,id,yy,tt)
C--   double precision function dqcEpmij(jset,id,iy,it)
C--   subroutine sqcQpQmyt(jset,id,y,t,qplus,qminu)
C--   subroutine sqcQpQmij(jset,id,iy,it,qplus,qminu)
C--   double precision function dqcSumQQByt(jset,def,y,t)
C--   double precision function dqcSumQQBij(jset,def,iy,it)
C--   double precision function dqcSumQQGyt(jset,def,y,t)
C--   subroutine sqcAllQpQmyt(jset,y,t,val)
C--   subroutine sqcAllQpQmij(jset,iy,it,val)
C--   double precision function dqcSplChk(jset,id,it)
C--   subroutine sqcEfromQQ(qvec,evec,nf)
C--   subroutine sqcEweedQQ(qvec,wt,id,n,nf)
C--   subroutine sqcPdfTab(jset,def,xx,nxx,qq,nqq,fff)
C--   subroutine sqcPdfLin(jset,def,iymi,iyma,izmi,izma,idout)

C===================================================================
C==   Pdf linear combinations ======================================
C===================================================================

C     ====================
      subroutine sqcPdfMat
C     ====================

C--   Setup the matrices which relate the flavour basis |q> to the
C--   to the evolution basis |e>
C-- 
C--   |ei> = tmateq(i,j) |qj>
C--   |qj> = tmatqe(i,j) |ei>    NB: |q> = |q +- qbar>
C--
C--            1  2  3  4  5  6  7  8  9 10 11 12
C--   |ei>     s  2+ 3+ 4+ 5+ 6+ v  2- 3- 4- 5- 6-
C--   |qi>     d+ u+ s+ c+ b+ t+ d- u- s- c- b- t-
C--
C--   Called by sqc_qcinit

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension imateq(6,6), work(6,6)

      data imateq /
C--      d  u  s  c  b  t
C--      1  2  3  4  5  6 
     +   1, 1, 1, 1, 1, 1, !si   =  1
     +  -1, 1, 0, 0, 0, 0, !ns1  =  2
     +   1, 1,-2, 0, 0, 0, !ns2  =  3
     +   1, 1, 1,-3, 0, 0, !ns3  =  4
     +   1, 1, 1, 1,-4, 0, !ns4  =  5
     +   1, 1, 1, 1, 1,-5/ !ns5  =  6

C--   Transform imateq to math notation (swap indices) -> umateq 
      do i = 1,6
        do j = 1,6
          umateq(i,j) = imateq(j,i)
        enddo
      enddo
C--   Loop over 3,..,6 flavors
      do nf = 3,6
C--     Invert the nf*nf submatrix of umateq
        call sqcOrtInv(umateq,work,6,nf)
C--     Store the result in common block 
        do i = 1,6
          do j = 1,6
            umatqe(i,j,nf) = work(i,j)
          enddo
        enddo 
      enddo
C--   Umatqe depends on nf but for tmatqe we can take nf = 6 because
C--   all 12 pdfs are always correctly filled by evolff even when nf < 6.
C--   Fill matrix tmateq and tmatqe
      do i = 1,6
        do j = 1,6
          tmateq(i,j)     = umateq(i,j)
          tmateq(i,j+6)   = 0.D0
          tmateq(i+6,j)   = 0.D0
          tmateq(i+6,j+6) = umateq(i,j)
          tmatqe(i,j)     = umatqe(i,j,6)
          tmatqe(i,j+6)   = 0.D0
          tmatqe(i+6,j)   = 0.D0
          tmatqe(i+6,j+6) = umatqe(i,j,6)
        enddo
      enddo

      return
      end

C     =================================
      subroutine sqcPdIdef(tmatpq,ierr)
C     =================================
  
C--   The set of input distns defined (by the user) on the q+-
C--   basis is stored in the matrix tmatpq(i,j), where i = 1,...,12 is
C--   the user's input pdf index and j = 1,...12 is the q+- index: 
C--            1  2  3  4  5  6  7  8  9 10 11 12
C--           d+ u+ s+ c+ b+ t+ d- u- s- c- b- t-
C--
C--   Here are the various transformations performed in QCDNUM:
C--
C--   |pi> = tmatpq(i,j) |qj>  pdfs written on the q+- basis (input)
C--   |qi> = tmatqp(i,j) |pj>  inverse of the above
C--   |qi> = tmatqe(i,j) |ej>  q->e transformation matrix (from sqcPdfMat)
C--   |ei> = tmateq(i,j) |qj>  inverse of above (from sqcPdfMat)
C--   |pi> = tmatpe(i,j) |ej>  pdfs written on the si/ns basis
C--   |ei> = tmatep(i,j) |pj>  si/ns basis as lin comb of input pdfs 
C--
C--   Given the matrices tmatpq, tmatqe and tmateq the routine calculates
C--   the matrices tmatpe and tmatep

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension tmatpq(12,12),tmatqp(12,12)
      dimension iwork(12)

C--   Invert tmatpq
      do i = 1,12
        do j = 1,12
          tmatqp(i,j)   = tmatpq(i,j)
        enddo
      enddo
      call smb_dminv(12,tmatqp,12,iwork,ierr)
C--   Oh, lala .... that doesnt look good      
      if(ierr.ne.0) return

C--   tmatpe = tmatpq*tmatqe
      do i = 1,12
        do j = 1,12
          sum = 0.D0
          do k = 1,12
            sum = sum + tmatpq(i,k)*tmatqe(k,j)
          enddo
          tmatpe(i,j) = sum
        enddo
      enddo
C--   tmatep = tmateq*tmatqp
      do i = 1,12
        do j = 1,12
          sum = 0.D0
          do k = 1,12
            sum = sum + tmateq(i,k)*tmatqp(k,j)
          enddo
          tmatep(i,j) = sum
        enddo
      enddo

      return
      end

C===================================================================
C==   Pdf input routines ===========================================
C===================================================================

C     ================================
      subroutine sqcAllInp(itype,func)
C     ================================

C--   Calculate 2nf+1 pdfs provided by the user in func(j,x) and 
C--   store these in      |p> (j=0,...,12)  
C--   Then transform to   |e> (i=0,...,12)
C--   Put the transformed pdfs in stor7(iy,it=0,i), i = 0,...,12
C--   NB: pdfs stored in bin it = 0 serve as evolution start values

      implicit double precision (a-h,o-z)

      external func

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension pval(0:12), eval(0:12)

C--   Initialize
      nfmax = max(nfix6,3)
      do i = 0,12
        pval(i) = 0.D0
        eval(i) = 0.D0
      enddo

C--   Loop over y gridpoints
      do iy = 1,nyy2(0)
        y  = ygrid2(iy)
        x  = exp(-y)
        do j = 0,2*nfmax
          pval(j) = func(j,x)
        enddo
C--     Transform
        call sqcEfrmP(pval,eval)
C--     Store
        do id = 0,12
          iadr = iqcPdfIjkl(iy,0,id,itype)
          stor7(iadr) = eval(id)
        enddo
      enddo

      return
      end

C     ==============================
      subroutine sqcEfrmP(pval,eval)
C     ==============================

C--   Transform 13 parton values |p> (j=0,...,12)
C--   to the evolution basis     |e> (i=0,...,12)
C--
C--   Input:  pval(0:12)  input values  |p>
C--   Output: eval(0:12)  output values |e>

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension pval(0:12), eval(0:12)

      do i = 0,12
        eval(i) = dqcEifrmP(i,pval)
      enddo

      return
      end

C     ===========================================
      double precision function dqcEifrmP(i,pval)
C     ===========================================

C--   Transform 13 parton values            |pj> (j=0,...,12)
C--   to one element of the evolution basis |ei> (i=0,...,12)
C--
C--   Input:  i          = index of |ei> basis element (i = 0 = gluon)
C--           pval(0:12) = input vector |p>  (j = 0 = gluon)
C--
C--   Output: dqcEifrmP  = value of evolution basis element |ei>

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension pval(0:12)

      if(i.eq.0) then
C--     Gluon
        dqcEifrmP = pval(0)
      else
C--     |ei>  = tmatep(i,j) |pj>
        sum = 0.D0
        do j = 1,12
          sum = sum + tmatep(i,j)*pval(j)
        enddo
        dqcEifrmP = sum
      endif

      return
      end
      
C     ====================================================
      subroutine sqcPdfInp(subr,jset,del,epsi,nwlast,ierr)
C     ====================================================

C--   Fill set of pdf tables. If jset does not exist, it is created 

C--   subr    (external) user supplied subroutine
C--   jset    (in)       pdf set [5-9]
C--   del     (in)       threshold offset mu2h*(1+-del)
C--   epsi    (out)      estimate of spline accuracy
C--   nwlast  (out)      last word occupied in the store (<0 no space)
C--   ierr    (out)      0=OK
C--
C--   subr(x,mu2,qqbar) should return qqbar(-6:6)   ... ub, db, gl, d, u, ...

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'qmaps8.inc'
      
      dimension qqbar(-6:6),qpm(12),epm(0:12)
      
      external subr

C--   Initialize
      ierr   = 0
      epsi    = 0.D0
      nwlast = nwlast7

C--   Yes/no initialize, thats the question
      if(.not.Lwtini8) call sqcIniWt
      if(.not.Lnfmap8) call sqcNfTab(0)
      if(.not.Lastab8) call sqcAlfTab(iord6)
      if(.not.Levcut8) call sqcEvCuts(0) 

      if(mxord7(jset).eq.0) then
C--     Book new set of tables      
        call sqcPdfIni(ipoint7,nwlast,idmin,idmax,lpdf)
C--     nwlast = last word used < 0 no space
        if(nwlast.le.0) return     
C--     Store pdf table offset
        knul7(jset) = kpdf7(0)
C--     Flag jset exists
        mxord7(jset) = 1
C--     Store last word used
        nwlast7 = nwlast
C--     Pointer to first word of next batch 
        ipoint7 = nwlast+1
      endif         

C--   Address increment
      inc = iqcPdfIjkl(1,1,2,jset)-iqcPdfIjkl(1,1,1,jset)
C--   Loop over iz and it      
      do iz = izmic2,izmac2
        it = itfiz2(iz)
        qi = exp(tgrid2(it))
C--     Detect threshold crossing: isign =  ... 0, 0, 0, -1, +1, 0, 0 ... 
        if(iz.ne.1 .and. iz.ne.nzz2) then
          isign = 2*nffiz2(iz)-nffiz2(iz-1)-nffiz2(iz+1)
        else
          isign = 0
        endif  
        qi = qi * (1.D0+isign*del)        
        do iy = 1,iymac2
          xi = exp(-ygrid2(iy))
C--       Get |qqbar>
          call subr(xi,qi,qqbar)
C--       Convert to |q+->
          do i = 1,6
            qpm(i)   = qqbar(i) + qqbar(-i)
            qpm(6+i) = qqbar(i) - qqbar(-i)
          enddo
C--       Gluon
          epm(0) = qqbar(0)
C--       Convert to |e+->          
          do i = 1,12
            epm(i) = 0.D0
            do j = 1,12
              epm(i) = epm(i) + tmateq(i,j)*qpm(j)
            enddo
          enddo
C--       Fill tables                                               
          ia  = iqcPdfIjkl(iy,iz,0,jset)-inc
          do i = 0,12
            ia        = ia+inc
            stor7(ia) = epm(i)
          enddo  
C--     End of loop over y
        enddo  
C--   End of loop over t
      enddo
      
C--   Check for spline oscillations
      epsi = 0.D0 
      do id = 0,12
        iq1 = itfiz2(izmic2)
        iq2 = itfiz2(izmac2)
        do iq = iq1,iq2     
          eps  = dqcSplChk(jset,id,iq)
          epsi = max(epsi,eps)
        enddo  
      enddo

      return
      end

C===================================================================
C==   Pdf output routines ==========================================
C===================================================================

C     =================================================      
      double precision function dqcXSplne(jset,id,y,it)
C     =================================================

C--   Spline interpolation in x

C--   jset   (in)  : Pdf set index [1-9]
C--   id     (in)  : Si/ns pdf index [0-12]
C--   y      (in)  : value of y
C--   it     (in)  : t grid point

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      
      logical lqcAxeqy
      
      dimension acoef(mxx0)
      
C--   Catch y = 0  (x = 1)      
      if(lqcAxeqy(y,0.D0,aepsi6)) then
        dqcXSplne = 0.D0
        return
      endif  
C--   Spline storage index
      idk = ioy2-1      
C--   Find y grid index in main grid G0
      iy = iqcFindIy(y)
C--   z grid index
      iz = izfit2(it) 
C--   Convert to A coefficients
      call sqcGetSplA(jset,id,iy,iz,ig,iyg,acoef)
C--   We did hit the end-point: iyg --> iyg-1
      iyg = min(iyg,nyy2(ig)-1)
C--   Spline index range in y
      call sqcByjLim(idk,iyg+1,minby,maxby)  !iy --> iy+1
C--   Loop over spline y
      val = 0.D0
      do jy = minby,maxby
        yjm1  = (jy-1)*dely2(ig)
        spy = dqcBsplyy(idk,1,(y-yjm1)/dely2(ig))
        val = val + acoef(jy)*spy
      enddo

      dqcXSplne = val

      return
      end
      
C     =================================================
      double precision function dqcEpmyt(jset,id,yy,tt)
C     =================================================

C--   Value of |gluon> or |e+-> at y and t
C--
C--   id  =  0   1   2   3   4   5   6   7   8   9  10  11  12 
C--          g  e1+ e2+ e3+ e4+ e5+ e6+ e1- e2- e3- e4- e5- e6-
C--
      
      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'qpars6.inc'
      
      logical lqcAxeqy

C--   Catch y = 0  (x = 1)      
      if(lqcAxeqy(yy,0.D0,aepsi6)) then
        dqcEpmyt = 0.D0
      else
        call sqcZmesh(yy,tt,0,iy1,iy2,iz1,iz2)
        call sqcIntPol(jset,id,iy1,iy2,iz1,iz2,yy,tt,dqcEpmyt)
      endif

      return
      end
      
C     =================================================
      double precision function dqcEpmij(jset,id,iy,it)
C     =================================================

C--   Value of |gluon> or |e+-> at iy and it
C--
C--   id  =  0   1   2   3   4   5   6   7   8   9  10  11  12 
C--          g  e1+ e2+ e3+ e4+ e5+ e6+ e1- e2- e3- e4- e5- e6-
C--
      
      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      iz = izfit2(it)
      ia = iqcPdfIjkl(iy,iz,id,jset)
      dqcEpmij = stor7(ia)

      return
      end      

C     =============================================
      subroutine sqcQpQmyt(jset,id,y,t,qplus,qminu)
C     =============================================

C--   Return gluon or qplus and qminus
C--   id     =  0  1  2  3  4  5  6
C--   val    =  g  d  u  s  c  b  t

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'qpars6.inc'
      
      call sqcZmesh(y,t,0,iy1,iy2,iz1,iz2)
      
C--   Output gluon
      if(id.eq.0) then
        call sqcIntPol(jset,id,iy1,iy2,iz1,iz2,y,t,qplus)
        qminu = qplus
      else
C--     Figure out number of flavors
        it = iqcItfrmt(t)
        if(it.eq.0)  stop 'sqcQpmyt: t out of range ---> STOP'
        iz  = izfit2(it)
        nf  = nffiz2(iz)
C--     Output |qi> = umatqe(i,j,nf) |ej>
C--     i  =  0  1  2  3  4  5  6  7  8  9 10 11 12 
C--     qi =  g  d+ u+ s+ c+ b+ t+ d- u- s- c- b- t-
C--     q+qbar and q-qbar
        qplus = 0.D0
        qminu = 0.D0
        do j = 1,nf
          call sqcIntPol(jset,j  ,iy1,iy2,iz1,iz2,y,t,eplus)
          call sqcIntPol(jset,j+6,iy1,iy2,iz1,iz2,y,t,eminu)
          qplus = qplus + umatqe(id,j,nf)*eplus
          qminu = qminu + umatqe(id,j,nf)*eminu
        enddo
      endif

      return
      end
      
C     ===============================================
      subroutine sqcQpQmij(jset,id,iy,it,qplus,qminu)
C     ===============================================

C--   Return gluon or qplus and qminus
C--   id     =  0  1  2  3  4  5  6
C--   val    =  g  d  u  s  c  b  t

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      
      iz = izfit2(it)
C--   Output gluon
      if(id.eq.0) then
        ia    = iqcPdfIjkl(iy,iz,id,jset)
        qplus = stor7(ia)
        qminu = qplus
      else
C--     Figure out number of flavors
        nf  = nffiz2(iz)
C--     Output |qi> = umatqe(i,j,nf) |ej>
C--     i  =  0  1  2  3  4  5  6  7  8  9 10 11 12 
C--     qi =  g  d+ u+ s+ c+ b+ t+ d- u- s- c- b- t-
C--     q+qbar and q-qbar
        qplus = 0.D0
        qminu = 0.D0
        do j = 1,nf
          ia    = iqcPdfIjkl(iy,iz,j  ,jset)
          eplus = stor7(ia)
          ia    = iqcPdfIjkl(iy,iz,j+6,jset)
          eminu = stor7(ia)
          qplus = qplus + umatqe(id,j,nf)*eplus
          qminu = qminu + umatqe(id,j,nf)*eminu
        enddo
      endif

      return
      end

C     ===================================================
      double precision function dqcSumQQByt(jset,def,y,t)
C     ===================================================

C--   Return linear combination of q and qbar
C--   id     =  -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
C--   def    =  tb bb cb sb ub db  g  d  u  s  c  b  t

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'
      
      dimension def(-6:6)
      
      sum = 0.D0
      do id = 1,6
        if(def(id).ne.0.D0 .or. def(-id).ne.0.D0) then
          call sqcQpQmyt(jset,id,y,t,qplus,qminu)
          wpl = 0.5D0 * (def(id)+def(-id))
          wmi = 0.5D0 * (def(id)-def(-id))
          sum = sum + wpl*qplus + wmi*qminu
        endif
      enddo

      dqcSumQQByt = sum

      return
      end

C     =====================================================
      double precision function dqcSumQQBij(jset,def,iy,it)
C     =====================================================

C--   Return linear combination of q and qbar
C--   id     =  -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
C--   def    =  tb bb cb sb ub db  g  d  u  s  c  b  t

      implicit double precision (a-h,o-z)

      dimension def(-6:6)
 
      sum = 0.D0
      do id = 1,6
        if(def(id).ne.0.D0 .or. def(-id).ne.0.D0) then
          call sqcQpQmij(jset,id,iy,it,qplus,qminu)
          wpl = 0.5D0 * (def(id)+def(-id))
          wmi = 0.5D0 * (def(id)-def(-id))
          sum = sum + wpl*qplus + wmi*qminu
        endif
      enddo

      dqcSumQQBij = sum

      return
      end
      
C     ===================================================
      double precision function dqcSumQQGyt(jset,def,y,t)
C     ===================================================

C--   Return linear combination of q, qbar and gluon
C--   id     =  -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
C--   def    =  tb bb cb sb ub db  g  d  u  s  c  b  t

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qpars6.inc'

      dimension def(-6:6)
      
      sum = 0.D0
      do id = 0,6
        if(def(id).ne.0.D0 .or. def(-id).ne.0.D0) then
          call sqcQpQmyt(jset,id,y,t,qplus,qminu)
          wpl = 0.5D0 * (def(id)+def(-id))
          wmi = 0.5D0 * (def(id)-def(-id))
          sum = sum + wpl*qplus + wmi*qminu
        endif
      enddo

      dqcSumQQGyt = sum

      return
      end

C     =====================================
      subroutine sqcAllQpQmyt(jset,y,t,val)
C     =====================================

C--   i      =  0  1  2  3  4  5  6  7  8  9 10 11 12 
C--   val(i) =  g  d+ u+ s+ c+ b+ t+ d- u- s- c- b- t-

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'qpars6.inc'

      dimension val(0:12)
      
      call sqcZmesh(y,t,0,iy1,iy2,iz1,iz2)
      
C--   Output gluon
      call sqcIntPol(jset,0,iy1,iy2,iz1,iz2,y,t,val(0))
C--   Output |qi> = tmatqe(i,j) |ej>
      do i = 1,12
        val(i) = 0.D0
      enddo
      do j = 1,12
        call sqcIntPol(jset,j,iy1,iy2,iz1,iz2,y,t,funj)
        do i = 1,12
          val(i) = val(i) + tmatqe(i,j)*funj
        enddo
      enddo

      return
      end
      
C     =======================================
      subroutine sqcAllQpQmij(jset,iy,it,val)
C     =======================================

C--   i      =  0  1  2  3  4  5  6  7  8  9 10 11 12 
C--   val(i) =  g  d+ u+ s+ c+ b+ t+ d- u- s- c- b- t-

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension val(0:12)

      iz = izfit2(it)
C--   Output gluon
      ia = iqcPdfIjkl(iy,iz,0,jset)
      val(0) = stor7(ia)
C--   Output |qi> = tmatqe(i,j) |ej>
      do i = 1,12
        val(i) = 0.D0
      enddo
      do j = 1,12
        ia   = iqcPdfIjkl(iy,iz,j,jset)
        funj = stor7(ia)
        do i = 1,12
          val(i) = val(i) + tmatqe(i,j)*funj
        enddo
      enddo

      return
      end
      
C     ===============================================
      double precision function dqcSplChk(jset,id,it)
C     ===============================================

C--   Returns epsi = || quad-lin || interpolation at midpoints
C--
C--   jset     (in) pdf set [1-9]
C--   id       (in) pdf identifier 0=gl, 1-12=epm singlet/non-singlet
C--   it       (in) mu2 grid point

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qgrid2.inc'
      
      dimension acoef(mxx0), epsi(mxx0)
 
C--   In  itialize 
      dqcSplChk = 0.D0     
C--   Linear interpolation
      if(ioy2.ne.3) return
C--   Find z-bin 
      iz = izfit2(it)       
C--   Loop over subgrids
      do jg = 1,nyg2
        iy = iyma2(jg)
        call sqcGetSplA(jset,id,iy,iz,ig,ny,acoef)
C--     Debug checks
        if(ig.ne.jg)        stop 'dqcSplChk: ig not jg'
        if(ny.ne.nyy2(jg))  stop 'dqcSplChk: ny not nyy2(jg)'
C--     Now get vector of deviations 
        nyma = iqcIyMaxG(iymac2,ig)       
        call sqcDHalf(ioy2,acoef,epsi,nyma)
C--     Max deviation
        do iy = 1,nyma
          dqcSplChk = max(dqcSplChk,abs(epsi(iy)))
        enddo  
      enddo
      
      return
      end
    
C     ===================================
      subroutine sqcEfromQQ(qvec,evec,nf)
C     ===================================

C--   Transform coefficients from flavor basis to si/ns basis
C--
C--         -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
C--   qvec  tb bb cb sb ub db  g  d  u  s  c  b  t
C--
C--   |qpm>  0  1  2  3  4  5  6  7  8  9 10 11 12
C--          g  d+ u+ s+ c+ b+ t+ d- u- s- c- b- t-
C--
C--   |epm>  0  1  2  3  4  5  6  7  8  9 10 11 12
C--          g  si 2+ 3+ 4+ 5+ 6+ v  2- 3- 4- 5- 6-
C--
C--   qvec(-6:6)   (in)  Coefficients in flavour space
C--                        qvec(0) is ignored --> quarks only!
C--   evec(12)     (out) Coefficients in si/ns space
C--   nf           (in)  Number of active flavours [3-6]

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpdfs7.inc'
      
      dimension qvec(-6:6),qpm(12),evec(12)

C--   Get qpm coefficient from q,qbar coefficients (eq 2.20)
      do i = 1,nf
        qpm(i  ) = 0.5*(qvec(i) + qvec(-i))
        qpm(6+i) = 0.5*(qvec(i) - qvec(-i))
      enddo
C--   Transform (eq 2.29)
      do i = 1,12
        evec(i) = 0.D0
      enddo
      do i = 1,nf
        diplu = 0.D0
        dimin = 0.D0
        do j = 1,nf
          diplu = diplu + qpm(j  )*umatqe(j,i,nf)
          dimin = dimin + qpm(6+j)*umatqe(j,i,nf)
        enddo
        evec(i  ) = diplu
        evec(i+6) = dimin
      enddo

      return
      end
      
C     ======================================
      subroutine sqcEweedQQ(qvec,wt,id,n,nf)
C     ======================================

C--   Transform coefficients from flavor basis to si/ns basis
C--   and weed out zero-value coefficients
C--
C--         -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
C--   qvec  tb bb cb sb ub db  g  d  u  s  c  b  t
C--
C--   qvec(-6:6)   (in)  Coefficients in flavour space
C--                        qvec(0) = nonzero --> pick gluon
C--                        qvec(0) = zero    --> pick quarks
C--   wt(i)        (out) Weight of i-th si/ns pdf
C--   id(i)        (out) Identifier of i-th si/ns pdf
C--   n            (out) Number of terms in weighted sum
C--   nf           (in)  Number of active flavours [3-6]

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'
      
      logical lqcAcomp
      
      dimension qvec(-6:6),wt(12),id(12)
      
      if(.not.lqcAcomp(qvec(0),0.D0,aepsi6)) then
C--     Gluon      
        wt(1) = qvec(0)
        id(1) = 0
        n     = 1
      else
C--     Quarks 
        call sqcEfromQQ(qvec,wt,nf)
        n = 0
        do i = 1,12 
          if(.not.lqcAcomp(wt(i),0.D0,aepsi6)) then
            n     = n + 1
            wt(n) = wt(i)
            id(n) = i
          endif              
        enddo
      endif          

      return
      end
      
C     ================================================      
      subroutine sqcPdfTab(jset,def,xx,nxx,qq,nqq,fff)
C     ================================================

C--   Return pdfs interpolated on a x-mu2 grid 
C--
C--   jset         (in)   pdf set [1-9]
C--   def(-6:6)    (in)   array of pdf coefficients
C--                       def(0) = 0.D0   select quarks
C--                       def(0) # 0.D0   select gluons 
C--   xx(i)        (in)   table of x-values
C--   nxx          (in)   number of x-values
C--   qq(i)        (in)   table of mu2 values
C--   nqq          (in)   number of mu2 values
C--   fff(nxx*nqq) (out)  interpolated pdfs (linear store)
C--
C--   Remark: xx(i) and qq(i) are all supposed to be within the grid
C--   Remark: result fff(nx*nq) is stored linearly 

      implicit double precision (a-h,o-z)
      
      dimension def(-6:6)
      dimension xx(*), qq(*), fff(*)
      
C--   Scratch buffer
      idout = -1      
      
C--   Find subgrid range in x and q2      
      ymi    = -log(xx(nxx))
      yma    = -log(xx(1))
      tmi    =  log(qq(1))
      tma    =  log(qq(nqq))
      margin = 0
      call sqcZmesh(ymi,tmi,margin,iymi,iydum,izmi,izdum)
      call sqcZmesh(yma,tma,margin,iydum,iyma,izdum,izma)
C--   Calculate lincomb of pdfs inside subgrid and store in idout
      call sqcPdfLin(jset,def,iymi,iyma,izmi,izma,idout)
C--   Loop over x and q and interpolate idout, store result in fff
      ia = 0
      do iq = 1,nqq
        t = log(qq(iq))
        do ix = 1,nxx
          ia = ia + 1
          y  = -log(xx(ix))
          call sqcZmesh(y,t,margin,iy1,iy2,iz1,iz2)
          call sqcIntPol(jset,idout,iy1,iy2,iz1,iz2,y,t,fff(ia))
        enddo         
      enddo
      
      return
      
      end

C     ========================================================      
      subroutine sqcPdfLin(jset,def,iymi,iyma,izmi,izma,idout)
C     ======================================================== 

C--   Store linear combination of pdfs in idout
C--
C--   jset         (in)   pdf set [1-9]
C--   def(-6:6)    (in)   array of pdf coefficients
C--                       def(0) = 0.D0   select quarks
C--                       def(0) # 0.D0   select gluons 
C--   iymi,iyma    (in)   range in iy
C--   izmi,izma    (in)   range in iz
C--   idout        (in)   output table identifier

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      
      dimension def(-6:6), wt(12,3:6), id(12,3:6), nn(3:6)
      
C--   Get coefficients in si/ns space      
      do nf = 3,6
        call sqcEweedQQ(def,wt(1,nf),id(1,nf),nn(nf),nf)
      enddo 

C--   Initialize output table
      call sqcPreSet(jset,idout,0.D0)      
C--   Sum linear combination
      do iz = izmi,izma 
C--     Output base address      
        jo = iqcPdfIjkl(iymi,iz,idout,jset)-1
C--     Figure out number of flavors
        nf  = nffiz2(iz)
C--     Loop over input identifiers        
        do i = 1,nn(nf)
C--       Output base address
          io = jo        
C--       Input base address        
          ia = iqcPdfIjkl(iymi,iz,id(i,nf),jset)-1          
          do iy = iymi,iyma
            ia        = ia + 1
            io        = io + 1
            stor7(io) = stor7(io) + wt(i,nf)*stor7(ia)
          enddo
        enddo
      enddo
          
      return      
      end
