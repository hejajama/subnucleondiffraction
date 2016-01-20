
C--   file usrfast.f containing fast convolution

C--   subroutine sqcUFBook(subnam)
C--   subroutine sqcUFIni(subnam,xlist,qlist,n,jchk)
C--   subroutine FastIni(xlist,qlist,n,jchk)
C--   subroutine FastClr(id)
C--   subroutine FastEpm(jset,id1,id2)
C--   subroutine FastSns(jset,qvec,isel,id)
C--   subroutine FastSum(iset,coef,id)
C--   subroutine FastFxK(w,idwt,id1,id2)
C--   subroutine FastFxF(w,idx,ida,idb,jdo)
C--   subroutine FastKin(id,fun)
C--   subroutine FastCpy(id1,id2,iadd)
C--   subroutine FastFxq(id,stf,n)

C==   ===============================================================
C==   Fast structure functions ======================================
C==   ===============================================================

C     ============================      
      subroutine sqcUFBook(subnam)
C     ============================

C--   Book scratch buffers
C--   Stop with error message if not enough space
C--   Called by FASTINI, PDFLST and PDFGRD

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      
      character*80 subnam 
      character*60 emsg
      character*10 etxt
      
C--   Book scratch buffers (or flag as empty if already booked)
      call sqcFastBook(nwords)
C--   Make sure that there is enough space to hold first word of next set
      nwfirst = abs(nwords)+1
      if(nwfirst.gt.nwf0) then
        call smb_itoch(nwfirst,etxt,ltxt)
        write(emsg,'(''Need at least '',A,
     &    '' words --> increase NWF0 '',
     &    ''in qcdnum.inc'')') etxt(1:ltxt)
        call sqcErrMsg(subnam,emsg)
      endif  
      
      return
      
      end
      
C     ==============================================      
      subroutine sqcUFIni(subnam,xlist,qlist,n,jchk)
C     ==============================================

C--   Setup list of interpolation points and clear buffers
C--   If jchk.ne.0 stop with error message when outside grid
C--   Called by FASTINI, PDFLST and PDFGRD

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'
      
      character*80 subnam
      
      dimension xlist(*),qlist(*)
      
      margin = 0
      call sqcSetMark(xlist,qlist,n,margin,ierr)
C--   At least one x,qmu2 outside grid
      if(jchk.eq.1 .and. ierr.eq.1) then
        call sqcErrMsg(subnam,'At least one x, mu2 outside cuts')
      endif
C--   Clear buffers
      do i = 1,mxord7(0)
         call sqcPreset(0,i,0.D0)
         isparse9(i) = 0          
      enddo

      return
      
      end
      
C     ======================================
      subroutine FastIni(xlist,qlist,n,jchk)
C     ======================================

C--   Pass list of x, mu2 values
C--   Also book scratch tables, if not already done
C--
C--   xlist = list of x values
C--   qlist = list of mu2 values
C--   n     = number of items in the list
C--   jchk  = 0 do not check xlist,qlist
C--         = 1 insist that all x,mu2 are inside grid   

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension xlist(*),qlist(*)
      
      character*80 subnam
      data subnam /'FASTINI ( X, QMU2, N, ICHK )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check user input 
      call sqcIlele(subnam,'N',1,n,mpt0,
     +    'You can increase mpt0 in qcdnum.inc (not recommended)')
C--   Book scratch buffers (or flag as empty if already booked)     
      call sqcUFBook(subnam)
C--   Setup list of interpolation points and clear buffers
      call sqcUFIni(subnam,xlist,qlist,n,jchk)

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
      
C     ======================
      subroutine FastClr(id)
C     ======================

C--   Clear scratch buffer. If id = 0, clear all buffers  

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      character*80 subnam
      data subnam /'FASTCLR ( ID )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check identifier
      call sqcIlele(subnam,'ID',0,id,mxord7(0),' ')
      if(id.eq.0) then
        idmin = 1
        idmax = mxord7(0)
      else
        idmin = id
        idmax = id
      endif
C--   Clear buffers
      do i = idmin,idmax
         call sqcPreset(0,i,0.D0)
         isparse9(i) = 0       
      enddo
C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
      
C     ==================================
      subroutine FastEpm(jset,idf,jdout)
C     ==================================

C--   Copy gluon or quark basis pdf to scratch table jdout
C--
C--   jset  (in) : pdf set identifier [1-9]
C--   idf   (in) : pdf indentifier [0,12]
C--   jdout (in) : id of output scratch table [1-mxord7(jset)]
C--                > 0 dense table
C--                < 0 sparse table

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension coef(0:12,3:6)

      character*80 subnam
      data subnam /'FASTEPM ( ISET, IDF, IDOUT )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Output identifier
      idout = abs(jdout)      
C--   Check jset      
      call sqcIlele(subnam,'ISET',1,jset,mset0,' ') 
C--   Check status bits
      call sqcChkflg(jset,ichk,subnam)
C--   Check identifier
      call sqcIlele(subnam,'IDF',0,idf,12,' ')
      call sqcIlele(subnam,'IDOUT',1,idout,mxord7(0),' ')
C--   Initialize
      isparse9(idout) = 0                             !empty table        
C--   Do the work
      do j = 3,6
        do i = 0,12
          coef(i,j) = 0.D0
        enddo
        coef(idf,j) = 1.D0
      enddo 
      if(jdout.gt.0) then
        isparse9(idout) = 2                           !dense table
        call sqcFastPdf(jset,coef,idout,1)
      else
        isparse9(idout) = 1                           !sparse table
        call sqcFastPdf(jset,coef,idout,0)
      endif   
      
C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
      
C     =====================================
      subroutine FastSns(jset,qvec,isel,jd)
C     =====================================

C--   Copy gluon or selected si/ns component to scratch table jd
C--
C--   jset        (in) : pdf set identifier [1-9]
C--   qvec(-6:6)  (in) : lin combination of q, qbar
C--   isel        (in) : selection flag [0-7]
C--   jd          (in) : id of output pdf table [1-mxord7(jset)]
C--                      > 0 dense table
C--                      < 0 sparse table

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension qvec(-6:6),coef(0:12,3:6)
      dimension mask(0:12,0:7)
C--              g s           v
C--              0 1 2 3 4 5 6 7 8 9 0 1 2      
      data mask /1,0,0,0,0,0,0,0,0,0,0,0,0,     !0=gluon
     +           0,1,0,0,0,0,0,0,0,0,0,0,0,     !1=singlet
     +           0,0,1,1,1,1,1,0,0,0,0,0,0,     !2=ns+
     +           0,0,0,0,0,0,0,1,0,0,0,0,0,     !3=valence
     +           0,0,0,0,0,0,0,0,1,1,1,1,1,     !4=ns-
     +           0,0,0,0,0,0,0,1,1,1,1,1,1,     !5=v and ns-
     +           0,0,1,1,1,1,1,1,1,1,1,1,1,     !6=all ns
     +           0,1,1,1,1,1,1,1,1,1,1,1,1   /  !7=all quarks

      character*80 subnam
      data subnam /'FASTSNS( ISET, DEF, ISEL, ID )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Output identifier
      id = abs(jd)      
C--   Check jset      
      call sqcIlele(subnam,'ISET',1,jset,mset0,' ') 
C--   Check status bits
      call sqcChkflg(jset,ichk,subnam)
C--   Check selection flag
      call sqcIlele(subnam,'ISEL',0,isel,7,' ')
C--   Check identifier
      call sqcIlele(subnam,'ID',1,id,mxord7(0),' ')
C--   Initialize
      isparse9(id) = 0                                !empty table       
C--   Do the work
      do nf = 3,6
        call sqcEfromQQ(qvec, coef(1,nf), nf)         !quark coefficients
        coef(0,nf) = coef(1,nf)                       !gluon = singlet
        do i = 0,12
          coef(i,nf) = coef(i,nf)*mask(i,isel)        !apply selection mask
        enddo
      enddo 
      if(jd.gt.0) then
        isparse9(id) = 2                              !dense table
        call sqcFastPdf(jset,coef,id,1)
      else
        isparse9(id) = 1                              !sparse table
        call sqcFastPdf(jset,coef,id,0)
      endif    

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
      
C     ================================
      subroutine FastSum(jset,coef,jd)
C     ================================

C--   Copy gluon or linear combination of quarks to scratch table id
C--
C--   jset            (in) : pdf set identifier [1-9]
C--   coef(0:12,3:6)  (in) : table of coefficients
C--   jd              (in) : id of output pdf table [1-mxord7(jset)]
C--                          > 0 dense table
C--                          < 0 sparse table

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension coef(0:12,3:6)

      character*80 subnam
      data subnam /'FASTSUM ( ISET, COEF, ID )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Output identifier
      id = abs(jd)      
C--   Check jset      
      call sqcIlele(subnam,'ISET',1,jset,mset0,' ') 
C--   Check status bits
      call sqcChkflg(jset,ichk,subnam)
C--   Check identifier
      call sqcIlele(subnam,'ID',1,id,mxord7(0),' ')
C--   Initialize
      isparse9(id) = 0                                !empty table       
C--   Do the work
      if(jd.gt.0) then
        isparse9(id) = 2                              !dense table
        call sqcFastPdf(jset,coef,id,1)
      else
        isparse9(id) = 1                              !sparse table
        call sqcFastPdf(jset,coef,id,0)
      endif 

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
      
C     ==================================
      subroutine FastFxK(w,idwt,id1,jd2)
C     ==================================

C--   Convolution F cross K
C--
C--   w       (in) : store filled with weight tables
C--   idwt(4) (in) : id(1)-(3) table ids LO, NLO, NNLO (0=no table) 
C--                  id(4)     leading power of alfas  (0 or 1)
C--   id1     (in) : scratch table w/pdfs previously filled by fastinp
C--   jd2     (in) : output scratch table
C--                  > 0 sparse table
C--                  < 0 sparse table

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension w(*),idwt(4)

      character*80 subnam
      data subnam /'FASTFXK ( W, IDW, ID1, ID2 )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Output identifier
      id2 = abs(jd2)
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check weight table identifiers
      do i = 1,3
        if(idwt(i).ne.0) then 
          ierr = iqcCheckId(w,idwt(i))
          if(ierr.ne.0) call sqcIdEmsg(subnam,'IDW(i)',idwt(i),ierr)   
        endif
      enddo
C--   Check power of alphas
      call sqcIlele(subnam,'IDW(4)',0,idwt(4),1,' ')
C--   No overwrite, thank you
      if(id1.eq.id2) then
        call sqcErrMsg(subnam,'ID1 cannot be equal to ID2')
      endif
C--   Check identifiers           
      call sqcIlele(subnam,'ID1',1,id1,mxord7(0),' ')
      call sqcIlele(subnam,'ID2',1,id2,mxord7(0),' ')
C--   Check id1 not empty or sparse
      if(isparse9(id1).eq.0) then
        call sqcErrMsg(subnam,'ID1 empty buffer')
      endif  
C--   Check id1 not sparse
      if(isparse9(id1).eq.1) then
        call sqcErrMsg(subnam,'ID1 sparse buffer')
      endif        
C--   Initialize
      isparse9(id2) = 0                                !empty table 
C--   Do the work
      if(jd2.gt.0) then
        isparse9(id2) = 1                              !sparse table
        call sqcFastFxK(w,idwt,id1,id2,0,ierr)
      else
        isparse9(id2) = 2                              !dense table
        call sqcFastFxK(w,idwt,id1,id2,1,ierr)
      endif 

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
      
C     =====================================
      subroutine FastFxF(w,idx,ida,idb,jdo)
C     =====================================

C--   Convolution F cross F
C--
C--   w       (in) : store filled with weight tables
C--   idx     (in) : weight table filled by makewtx
C--   ida,idb (in) : scratch tables with pdfs
C--   jdo     (in) : output scratch table
C--                  > 0 sparse table
C--                  < 0 sparse table

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension w(*)

      character*80 subnam
      data subnam /'FASTFXF ( W, IDX, IDA, IDB, IDOUT )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Output identifier
      ido = abs(jdo)
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check weight table identifier
      ierr = iqcCheckId(w,idx)
      if(ierr.ne.0) call sqcIdEmsg(subnam,'IDX',idx,ierr)   
C--   No overwrite, thank you
      if(ido.eq.ida .or. ido.eq.idb) then
        call sqcErrMsg(subnam,'IDOUT cannot be equal to IDA or IDB')
      endif
C--   Check identifiers           
      call sqcIlele(subnam,'IDA'  ,1,ida,mxord7(0),' ')
      call sqcIlele(subnam,'IDB'  ,1,idb,mxord7(0),' ')
      call sqcIlele(subnam,'IDOUT',1,ido,mxord7(0),' ')
C--   Check ida,b not empty
      if(isparse9(ida).eq.0) then
        call sqcErrMsg(subnam,'IDA empty buffer')
      endif
      if(isparse9(idb).eq.0) then
        call sqcErrMsg(subnam,'IDB empty buffer')
      endif   
C--   Check ida,b not sparse
      if(isparse9(ida).eq.1) then
        call sqcErrMsg(subnam,'IDA sparse buffer')
      endif
      if(isparse9(idb).eq.1) then
        call sqcErrMsg(subnam,'IDB sparse buffer')
      endif        
C--   Initialize
      isparse9(ido) = 0                                !empty table 
C--   Do the work
      if(jdo.gt.0) then
        isparse9(ido) = 1                              !sparse table
        call sqcFastFxF(w,idx,ida,idb,ido,0)
      else
        isparse9(ido) = 2                              !dense table
        call sqcFastFxF(w,idx,ida,idb,ido,1)
      endif 

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C     ==========================
      subroutine FastKin(id,fun)
C     ==========================

C--   Multiply contents of pdf table id by fun(ix,iq)
C--
C--   Input:   id   = scratch table filled with convolutions
C--            fun  = user defined function fun(ix,iq,nf,ithresh)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      external fun

      character*80 subnam
      data subnam /'FASTKIN ( ID, FUN )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check identifier
      call sqcIlele(subnam,'ID',1,id,mxord7(0),' ')
C--   Check id not empty
      if(isparse9(id).eq.0) then
        call sqcErrMsg(subnam,'ID empty buffer')
      endif       
      
C--   Do the work
      call sqcFastKin(id,fun)

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
      
C     ================================
      subroutine FastCpy(id1,id2,iadd)
C     ================================

C--   Add contents of id1 to id2

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      character*80 subnam
      data subnam /'FASTCPY ( ID1, ID2, IADD )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check identifiers
      if(id1.eq.id2) then
        call sqcErrMsg(subnam,'ID1 cannot be equal to ID2')
      endif
      call sqcIlele(subnam,'ID1',1,id1,mxord7(0),' ')
      call sqcIlele(subnam,'ID2',1,id2,mxord7(0),' ')
C--   Check id1 not empty
      if(isparse9(id1).eq.0) then
        call sqcErrMsg(subnam,'ID1 empty buffer')
      endif       
C--   Check iadd
      call sqcIlele(subnam,'IADD',-1,iadd,1,' ')
C--   Sparse or dense output thats the question
      if(isparse9(id2).eq.0 .or. iadd.eq.0) then
        isparse9(id2) = isparse9(id1)
      else  
        isparse9(id2) = min(isparse9(id1),isparse9(id2))
      endif             
C--   Do the work
      call sqcFastCpy(id1,id2,iadd,isparse9(id2)-1)

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C     ============================
      subroutine FastFxq(id,stf,n)
C     ============================

C--   Interpolation
C--
C--   Input:   id    = scratch table filled with structure function
C--            xlst9 = list of x values in /qfast9/
C--            qlst9 = list of qmu2 values in /qfast9/ 
C--            n     = number of interpolations requested = min(n,nlst9)
C--
C--   Output:  stf   = array of interpolated stf values

C--   NB: the list of x,qmu2 values is entered via s/r fastini

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'
      include 'qfast9.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension stf(*)

      character*80 subnam
      data subnam /'FASTFXQ ( ID, STF, N )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check identifier
      call sqcIlele(subnam,'ID',1,id,mxord7(0),' ')
C--   Check id not empty
      if(isparse9(id).eq.0) then
        call sqcErrMsg(subnam,'ID empty buffer')
      endif       
      call sqcIlele(subnam,'N',1,n,mpt0,
     +        'You can increase MPT0 in QCDNUM.INC and recompile')
C--   Do the work
      call sqcFastFxq(0,id,stf,n)

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
