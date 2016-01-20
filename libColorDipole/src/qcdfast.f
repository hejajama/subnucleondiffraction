
C--   file qcdfast.f containing fast convolution routines

C--   subroutine sqcFastBook(nwords)
C--   subroutine sqcSetMark(xlist,qlist,n,margin,ierr)
C--   subroutine sqcFastPdf(jset,coef,id,idense)
C--   subroutine sqcFastAdd(jset1,id,wt,n,jset2,id2,
C--                         nzlist,izlist,nylist,iylist)
C--   subroutine sqcFastFxK(w,idwt,idi,ido,idense,ierr)
C--   subroutine sqcFccAtIt(w,idwt,idi,ido,list,nl,iz)
C--   subroutine sqcFastWgt(w,idwt,iz,nf,ig,wmat)
C--   subroutine sqcFastFxF(w,idx,ida,idb,idout,idense)
C--   subroutine sqcFcfAtIt(w,idx,ida,idb,ido,list,nl,iz)
C--   subroutine sqcFastKin(id,fun)
C--   subroutine sqcFastCpy(id1,id2,iadd,idense)
C--   subroutine sqcFastFxq(jset,id,stf,n)

C     ==============================================================
C     Fast structure function calculation ==========================
C     ==============================================================

C     ==============================
      subroutine sqcFastBook(nwords)
C     ==============================

C--   Book set of scratch tables
C--  
C--   nwords (out) : last word used in the store (<0 no space)

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

C--   Check tables already booked      
      if(mxord7(0).ne.0) then
        nwords = nwlast7
        return
      endif 
      
C--   Book tables
      call sqcPdfBuf(nscratch6,ipoint7,nwords,knul)
C--   Not enough space
      if(nwords.lt.0) return
C--   Update pointers
      nwlast7   = nwords
      ipoint7   = nwords+1
      mxord7(0) = nscratch6
      knul7(0)  = knul
      
      return
      end          

C     ================================================
      subroutine sqcSetMark(xlist,qlist,n,margin,ierr)
C     ================================================

C--   Mark grid points
C--
C--   Input:  xlist      =  list of x   points
C--           qlist      =  list of mu2 points
C--           n          =  number of items in xlist and qlist
C--           margin     =  points to stay away from threshold [0,1]
C--
C--   Output: ierr       =  0    all OK
C--                         1    at least one x,q outside cuts
C--   Output in qfast9.inc
C--           mark9(0:mxx0,0:mqq0+7) = logical array with marked grid points
C--           xlst9,qlst9,nxq9       = copy of xlist, qlist and n 
C--           ylst9,tlst9,nyt9       = contains only points inside grid
C--           ixqfyt9                = index from ylst,tlst --> xlst,qlst 
C--           iy19,iy29,iz19,iz29    = interpolation mesh for each entry
C--                                    in ylst,tlst          

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qfast9.inc'

      dimension xlist(*), qlist(*)
      logical   lqcInside

C--   Check
      if(n.gt.mpt0) stop 'sqcSetMark: too many points n ---> STOP'
      
C--   Copy to common block
      nxq9 = n
      nyt9 = 0
      ierr = 0
      do i = 1,n
        xlst9(i) = xlist(i)
        qlst9(i) = qlist(i)
C--     Process only points inside grid        
        if(lqcInside(xlist(i),qlist(i),ifail)) then
          nyt9         =  nyt9+1   
          ylst9(nyt9)  = -log(xlist(i))
          tlst9(nyt9)  =  log(qlist(i))
C--       Remember position in original list
          ixqfyt9(nyt9) = i
        else
          ierr = 1            
        endif
      enddo

C--   Clear markers
      do j = 0,mqq0+7
        do i = 0,mxx0
          mark9(i,j) = .false.
        enddo
      enddo

C--   Put markers in the grid
      call sqcMarkit(
     +     mark9,ylst9,tlst9,margin,iy19,iy29,iz19,iz29,nyt9)
      
C--   Set up sparse and dense lists
C--   Loop over z points
      nz    = 0
      iymax = 0
      do iz = 1,nzz2
        ny = 0
        do iy = 1,nyy2(0)
          if(mark9(iy,iz)) then
            ny    = ny+1
            iymax = iy
          endif  
        enddo
        if(ny.ne.0) then
          nz           = nz+1
          izlist9(nz)  = iz
          nyslist9(nz) = ny
          nydlist9(nz) = iymax
        endif
        ny = 0
        do iy = 1,iymax
          iydlist9(iy,nz) = iy
          if(mark9(iy,iz)) then
            ny              = ny+1
            iyslist9(ny,nz) = iy
          endif                    
        enddo
      enddo
      nzlist9 = nz
        
      return
      end

C     ==========================================
      subroutine sqcFastPdf(jset,coef,id,idense)
C     ==========================================

C--   Copy linear combination of pdfs to id
C--
C--   jset             (in) : pdf set indentifier
C--   coef(0:12,3:6)   (in) : coefficients of the linear combination
C--   id               (in) : scratch buffer identifier
C--   idense           (in) : 0/1 no/yes dense table

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'

      logical lqcAcomp
      dimension coef(0:12,3:6),cvec(12,3:6),ivec(12)
      
C--   Initialize
      call sqcPreset(0,id,0.D0)
C--   Weedout zero coefficients
      nvec = 0
      do i = 0,12
        istore = 0
        do  nf = 3,6
          if(.not.lqcAcomp(coef(i,nf),0.D0,aepsi6)) istore = 1
        enddo
        if(istore.eq.1) then
          nvec       = nvec+1
          if(nvec.gt.12) stop 'sqcFastPdf: nvec larger than 12'
          ivec(nvec) = i
          do nf = 3,6
            cvec(nvec,nf) = coef(i,nf)
          enddo
        endif
      enddo 
C--   Summitup 
      if(nvec.ne.0) then
        if(idense.eq.1) then
          call sqcFastAdd(jset,ivec,cvec,nvec,0,id,
     +                    nzlist9,izlist9,nydlist9,iydlist9)
        else
          call sqcFastAdd(jset,ivec,cvec,nvec,0,id,
     +                    nzlist9,izlist9,nyslist9,iyslist9)     
        endif
      endif
      
      return
      end

C     ==================================================
      subroutine sqcFastAdd(jset1,id,wt,n,jset2,id2,
     +                      nzlist,izlist,nylist,iylist)
C     ==================================================

C--   Weigted sum of tables with weights dependent on nf.

C--   jset1         pdf set identifier of input table 
C--   id(i)         list of n input identifiers declared id(12) in calling routine
C--   wt(i,nf)      list of weights declared wt(12,3:6) in calling routine 
C--   n             index i in id(i) and wt(i,nf) runs from 1,...,n
C--   jset2         pdf set identifier of the output table
C--   id2           identifier where weighted sum will be stored
C--   nzlist        number of z grid points to loop over
C--   izlist(j)     grid point jz
C--   nylist(j)     number of y grid points to loop over at jz
C--   iylist((i,j)  grid point iy at jz

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      dimension id(12),wt(12,3:6)
      dimension izlist(*),nylist(*),iylist(mxx0,*)

C--   Protect
      if(n.lt.1.or.n.gt.12) stop 
     +      'sqcFastAdd: n not in range [1,12] ---> STOP'
      do i = 1,n
        if(jset1.eq.jset2 .and. id(i).eq.id2) stop 
     +      'sqcFastAdd: attempt to overwrite input id ---> STOP'
      enddo
C--   Initialize output table
      call sqcPreset(jset2,id2,0.D0)
C--   Loop over table id's and fill target id
      do k = 1,n
        do j = 1,nzlist
          iz  = izlist(j)
          ia1 = iqcPdfIjkl(1,iz,id(k),jset1)-1
          ia2 = iqcPdfIjkl(1,iz,id2  ,jset2)-1
          nf  = nffiz2(iz)
          wgt = wt(k,nf)
          do i = 1,nylist(j)
            iy = iylist(i,j)
            stor7(ia2+iy) = stor7(ia2+iy) + wgt*stor7(ia1+iy)
          enddo
        enddo
      enddo
      
      return
      end

C     =================================================
      subroutine sqcFastFxK(w,idwt,idi,ido,idense,ierr)
C     =================================================

C--   Calculate FcrossK for a list of y,t values
C--
C--   w        (in)  : store with weight tables
C--   idwt(4)  (in)  : id(1)-(3) table ids LO, NLO, NNLO (0=no table) 
C--                    id(4)     leading power of alfas  (0 or 1)
C--   idi      (in)  : identifier of table with input pdf
C--   ido      (in)  : identifier of table with output convolution
C--   idense   (in)  : 0/1 = no/yes dense output
C--   ierr     (out) : 0 = OK
C--                    1 = at least one grid point below alphaslim

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qfast9.inc'

      dimension w(*),idwt(*)

C--   Loop over iz-bins
      ierr = 0
      do j = 1,nzlist9
        iz = izlist9(j)
        it = itfiz2(iz)
C--     No alphas at such low t-bin
        if(it.lt.itmin6) ierr = 1
C--     Calculate stf for marked points and put result in ido
        if(idense.eq.0) then
          call sqcFccAtIt(w,idwt,idi,ido,iyslist9(1,j),nyslist9(j),iz)
        else
          call sqcFccAtIt(w,idwt,idi,ido,iydlist9(1,j),nydlist9(j),iz) 
        endif
      enddo  

      return
      end

C     ================================================
      subroutine sqcFccAtIt(w,idwt,idi,ido,list,nl,iz)
C     ================================================

C--   Calculate FcrossK for a list of grid points iy at fixed iz

C--   Input: w        store containing weight tables
C--          idwt(4)  weight table ids
C--          idi      identifier of input  pdf table
C--          ido      identifier of output fxk table
C--          list     list of grid points in y
C--          nl       number of grid points in the list
C--          iz       t-grid index

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qmaps8.inc'

      dimension w(*),coef(mxx0)
      dimension list(mxx0)
      dimension idwt(4)
      dimension wmat(mxx0)

C--   Value of t
      it  = itfiz2(iz)
      tt  = tgrid2(it)
      nf  = nffiz2(iz) 

C--   Base address
      iao = iqcPdfIjkl(1,iz,ido,0)-1

C--   Check list in ascending order
      if(list(nl).lt.list(1)) stop 'sqcFccAtIt: wrong y-loop'      
C--   Loop over all requested grid points
      iglast = 0
C--   Must loop from high iy to low iy!      
      do ii = nl,1,-1
        jy  = list(ii)
        if(jy.eq.0) then
C--       Lower edge of y-grid: set convolution to zero
          fyj = 0.D0
        else
C--       Gridpoint above lower edge of y-grid
          yj  = ygrid2(jy)
          ig  = iqcFindIg(yj)
C--       New subgrid encountered
          if(ig.ne.iglast) then
C--         Fill weight table
            call sqcFastWgt(w,idwt,iz,nf,ig,wmat)
C--         Convert F values to A values
            call sqcGetSplA(0,idi,jy,iz,jg,jyg,coef)
            iglast = ig
          endif
C--       Done with setting-up subgrid; get y-index in subgrid 
          ky  = iqcIyfrmY(yj,dely2(ig),nyy2(ig))
C--       Convolution loop in subgrid
          fyj = 0.D0
          do i = 1,ky
            fyj = fyj + wmat(ky+1-i)*coef(i)
          enddo
        endif
C--     Store FxK in output table
        if(it.ge.itlow8) then
          stor7(iao+jy) = fyj      !alphas available for this t-bin
        else
          stor7(iao+jy) = qnull6   !no alphas for this t-bin
        endif
C--     End of loop over requested grid points
      enddo

      return
      end

C     ===========================================
      subroutine sqcFastWgt(w,idwt,iz,nf,ig,wmat)
C     ===========================================

C--   Return weight table or perturbative expansion of tables
C--
C--   Input:  w        store containing weight tables
C--           idwt(4)  id(1)-(3) table ids LO, NLO, NNLO (0=no table)
C--                    id(4)     leading power of alfas  (0 or 1)
C--           iz       z-grid index
C--           nf       number of flavors
C--           ig       subgrid index
C--   Output: wmat     weight table

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qmaps8.inc'

      dimension w(*),idwt(*),wmat(*),mi(5),ma(5)

      it = itfiz2(iz)
      
      do iy = 1,nyy2(ig)
        wmat(iy) = 0.D0
      enddo
      
C--   Make sure astable exists
      if(.not.Lastab8) call sqcAlfTab(iord6)

      do io = 1,iord6
        id = idwt(io)
        if(id.ne.0) then
          ityp = id/100
          call sqcTabLims(w,ityp,mi,ma,jerr)
          if(jerr.ne.0) stop 'sqcFastWt: no weight table of this type'
          jt = max(it,mi(2))
          jt = min(jt,ma(2))
          kf = max(nf,mi(3))
          kf = min(kf,ma(3))
          ia = iqcWaddr(w,1,jt,kf,ig,id)-1
C--       Fill weight table
          if(idwt(4).eq.0) then
C--         Multiply LO,NLO,NNLO by 1,as,as2            
            if(io.eq.1) then
              do iy = 1,nyy2(ig)
                ia = ia + 1
                wmat(iy) = wmat(iy) + w(ia)
              enddo  
            else
              do iy = 1,nyy2(ig)
                ia = ia + 1
                wmat(iy) = wmat(iy) + antab8(iz,-(io-1))*w(ia)
              enddo
            endif
          else
C--         Multiply LO,NLO,NNLO by as,as2,as3
            do iy = 1,nyy2(ig)
              ia = ia + 1 
              wmat(iy) = wmat(iy) + antab8(iz,io)*w(ia)
            enddo                 
          endif
        endif
      enddo

      return
      end
      
C     =================================================
      subroutine sqcFastFxF(w,idx,ida,idb,idout,idense)
C     =================================================

C--   Calculate FcrossF for a list of y,t values
C--
C--   w        (in)  : store with weight tables
C--   idx      (in)  : weight table previously filled by MakeWtX
C--   ida,idb  (in)  : identifiers of tables with input pdfs
C--   idout    (in)  : identifier of table with output convolution
C--   idense   (in)  : 0/1 = no/yes dense output
C--   ierr     (out) : 0 = OK
C--                    1 = at least one grid point below alphaslim

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qfast9.inc'
      
      dimension w(*)

C--   Loop over iz-bins
      do j = 1,nzlist9
        iz = izlist9(j)
        it = itfiz2(iz)
C--     Calculate stf for marked points and put result in idout
        if(idense.eq.0) then
          call sqcFcfAtIt(w,idx,ida,idb,idout,
     +                    iyslist9(1,j),nyslist9(j),iz)
        else
          call sqcFcfAtIt(w,idx,ida,idb,idout,
     +                    iydlist9(1,j),nydlist9(j),iz) 
        endif
      enddo
      
      return
      end

C     ===================================================      
      subroutine sqcFcfAtIt(w,idx,ida,idb,ido,list,nl,iz)
C     ===================================================

C--   Calculate FcrossF for a list of grid points iy at fixed iz

C--   Input: w        store containing weight tables
C--          idx      weight table id
C--          ida,idb  identifiers of input pdf tables
C--          ido      identifier of output fxf table
C--          list     list of grid points in y
C--          nl       number of grid points in the list
C--          iz       t-grid index

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      
      dimension w(*),coefa(mxx0),coefb(mxx0)
      dimension list(mxx0)
      
C--   Value of t
      it  = itfiz2(iz)
      tt  = tgrid2(it)
      nf  = nffiz2(iz)

C--   Base address
      iao = iqcPdfIjkl(1,iz,ido,0)-1
      
C--   Check list in ascending order
      if(list(nl).lt.list(1)) stop 'sqcFcfAtIt: wrong y-loop'
C--   Loop over all requested grid points
      iglast = 0
C--   Must loop from high iy to low iy!      
      do ii = nl,1,-1
        jy  = list(ii)
        if(jy.eq.0) then
C--       Lower edge of y-grid: set convolution to zero
          fxf = 0.D0
        else
C--       Gridpoint above lower edge of y-grid
          yj  = ygrid2(jy)
          ig  = iqcFindIg(yj)
          iwx = 0                   !avoid compiler warning
C--       New subgrid encountered
          if(ig.ne.iglast) then
C--         Convert F values to A values
            call sqcGetSplA(0,ida,jy,iz,jg,jyg,coefa)
            call sqcGetSplA(0,idb,jy,iz,jg,jyg,coefb)
            iglast = ig
C--         Weight table base address
            iwx = iqcWaddr(w,1,it,nf,ig,idx)-1
          endif
C--       Done with setting-up subgrid; get y-index in subgrid 
          ky  = iqcIyfrmY(yj,dely2(ig),nyy2(ig))
C--       Convolution loop in subgrid
          fxf = 0.D0
          do j = 1,ky
            Aj = coefa(j)
            do k = 1,ky-j+1
              Bk  = coefb(k) 
              fxf = fxf + Aj*Bk*w(iwx+ky-j-k+2)
            enddo
          enddo
        endif
C--     Store FxF in output table
        stor7(iao+jy) = fxf
C--     End of loop over requested grid points
      enddo                         
      
      return
      end            
      
C     =============================
      subroutine sqcFastKin(id,fun)
C     =============================

C--   Multiply stf  with x-q dependent factor
C--
C--   Input:  id   =  Stf scratch table identifier
C--           fun  =  function of ix, iq, nf and ithres

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'

      external fun

C--   Loop over grid points
      if(isparse9(id).eq.1) then                      !sparse output
        do j     = 1,nzlist9
          iz     = izlist9(j)
          it     = itfiz2(iz)
          nf     = nffiz2(iz)
          ithres = 0
          if(iz.ne.1 .and. iz.ne.nzz2) then
            if(nffiz2(iz+1).eq.nffiz2(iz)+1) ithres = -1
            if(nffiz2(iz-1).eq.nffiz2(iz)-1) ithres =  1
          endif  
          ia = iqcPdfijkl(1,iz,id,0)-1
          do i = 1,nyslist9(j)
            iy = iyslist9(i,j)
            ix = nyy2(0)-iy+1
            stor7(ia+iy) = stor7(ia+iy)*fun(ix,it,nf,ithres)
          enddo
        enddo
      else                                           !dense output
        do j = 1,nzlist9
          iz     = izlist9(j)
          it     = itfiz2(iz)
          nf     = nffiz2(iz)
          ithres = 0
          if(iz.ne.1 .and. iz.ne.nzz2) then
            if(nffiz2(iz+1).eq.nffiz2(iz)+1) ithres = -1
            if(nffiz2(iz-1).eq.nffiz2(iz)-1) ithres =  1
          endif 
          ia = iqcPdfijkl(1,iz,id,0)-1
          do i = 1,nydlist9(j)
            iy = iydlist9(i,j)
            ix = nyy2(0)-iy+1
            stor7(ia+iy) = stor7(ia+iy)*fun(ix,it,nf,ithres)
          enddo
        enddo  
      endif  

      return
      end

C     ==========================================
      subroutine sqcFastCpy(id1,id2,iadd,idense)
C     ==========================================

C--   Add content of stf table id1 to that of id2

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'
      
      logical first
      save    first
      data    first /.true./
      
C--   Initialize at first call
      if(first) then
        call sqcPreset(0,id2,0.D0)
        first = .false.
      endif

      if(idense.eq.0) then                          !sparse output
      
        if(iadd.eq.-1) then
          do j = 1,nzlist9
            iz    = izlist9(j) 
            ia1   = iqcPdfijkl(1,iz,id1,0)-1
            ia2   = iqcPdfijkl(1,iz,id2,0)-1
            do i = 1,nyslist9(j)
              iy = iyslist9(i,j) 
              stor7(ia2+iy) = stor7(ia2+iy)-stor7(ia1+iy)
            enddo
          enddo
        elseif(iadd.eq.0) then
          do j = 1,nzlist9
            iz    = izlist9(j) 
            ia1   = iqcPdfijkl(1,iz,id1,0)-1
            ia2   = iqcPdfijkl(1,iz,id2,0)-1
            do i = 1,nyslist9(j)
              iy = iyslist9(i,j)
              stor7(ia2+iy) = stor7(ia1+iy)
            enddo
          enddo  
        elseif(iadd.eq.1) then
          do j = 1,nzlist9
            iz    = izlist9(j) 
            ia1   = iqcPdfijkl(1,iz,id1,0)-1
            ia2   = iqcPdfijkl(1,iz,id2,0)-1
            do i = 1,nyslist9(j)
              iy = iyslist9(i,j)
              stor7(ia2+iy) = stor7(ia2+iy)+stor7(ia1+iy)
            enddo
          enddo  
        else
          stop 'sqcFastCpy: invalid iadd'   
        endif
      
      else                                          !dense output
      
        if(iadd.eq.-1) then
          do j = 1,nzlist9
            iz    = izlist9(j) 
            ia1   = iqcPdfijkl(1,iz,id1,0)-1
            ia2   = iqcPdfijkl(1,iz,id2,0)-1
            do i = 1,nydlist9(j)
              iy = iydlist9(i,j) 
              stor7(ia2+iy) = stor7(ia2+iy)-stor7(ia1+iy)
            enddo
          enddo
        elseif(iadd.eq.0) then
          do j = 1,nzlist9
            iz    = izlist9(j) 
            ia1   = iqcPdfijkl(1,iz,id1,0)-1
            ia2   = iqcPdfijkl(1,iz,id2,0)-1
            do i = 1,nydlist9(j)
              iy = iydlist9(i,j)
              stor7(ia2+iy) = stor7(ia1+iy)
            enddo
          enddo  
        elseif(iadd.eq.1) then
          do j = 1,nzlist9
            iz    = izlist9(j) 
            ia1   = iqcPdfijkl(1,iz,id1,0)-1
            ia2   = iqcPdfijkl(1,iz,id2,0)-1
            do i = 1,nydlist9(j)
              iy = iydlist9(i,j)
              stor7(ia2+iy) = stor7(ia2+iy)+stor7(ia1+iy)
            enddo
          enddo  
        else
          stop 'sqcFastCpy: invalid iadd'   
        endif              
        
      endif
      
      return
      end
      
C     ====================================
      subroutine sqcFastFxq(jset,id,stf,n)
C     ====================================

C--   Interpolation

C--   jset    (in)  : pdf set identifier [0-9]
C--   id      (in)  : identifier of (scratch) table
C--   stf(i)  (out) : list of interpolated results
C--   n       (in)  : dimension of stf in the calling routine
C--
C--   Input via qfast9.inc : this common block is filled by sqcSetMark
C--   sqcSetMark is called by FastIni
C--   sqcSetMark is called by sqcStInterp is called by StfunXq
C--
C--   ylst9(i)  (in) :  list of y points
C--   tlst9(i)  (in) :  list of t points
C--   iy19(i)   (in) :  list of lower mesh limits in y
C--   iy29(i)   (in) :  list of upper mesh limits in y
C--   iz19(i)   (in) :  list of lower mesh limits in z
C--   iz29(i)   (in) :  list of upper mesh limits in z
C--   ixqfyt(i) (in) :  pointer to index in stf(n)
C--   nyt9      (in) :  number of items in the lists above

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'qfast9.inc'

      dimension stf(*)
      
C--   Preset for points outside grid
      do i = 1,n
        stf(i) = qnull6
      enddo   
C--   Now interpolate for points inside grid
      do i = 1,min(n,nyt9)
        call sqcIntPol(jset,
     +    id,iy19(i),iy29(i),iz19(i),iz29(i),ylst9(i),tlst9(i),val)
        stf(ixqfyt9(i)) = val
      enddo
      
      return
      end
