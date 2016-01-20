
C--   file usrcvol.f containing convolution engine user interface

C--   subroutine setUmsg(name)
C--   subroutine clrUmsg

C--   subroutine BookTab(w,nw,jtypes,nwds)

C--   subroutine SetWpar(w,par,n)
C--   subroutine GetWPar(w,par,n)

C--   subroutine TabDump(w,lun,file,key)

C--   subroutine TabRead(w,nw,lun,file,key,nwords,ierr)

C--   subroutine MakeWtA(w,id,afun,achi)
C--   subroutine MakeWtB(w,id,bfun,achi,ndel)
C--   subroutine MakeWRS(w,id,rfun,sfun,achi,ndel)
C--   subroutine MakeWtD(w,id,dfun,achi)
C--   subroutine MakeWtX(w,id)

C--   subroutine ScaleWt(w,c,id)
C--   subroutine CopyWgt(w,jd1,id2,iadd)
C--   subroutine WcrossW(w,jda,jdb,idc,iadd)
C--   subroutine WtimesF(w,fun,jd1,id2,iadd)

C--   integer function idSpfun(string,iord,jset)
C--
C--   double precision function FcrossK(w,idw,jset,idf,ix,iq)
C--   double precision function FcrossF(w,idw,jset,ida,idb,ix,iq)
C--   integer function Nflavor(iq)
C--   double precision function GetAlfN(iq,n,ierr)
C--   subroutine EfromQQ(qvec,evec,nf)
C--   subroutine QQfromE(evec,qvec,nf)
C--
C--   subroutine StfunXq(stfun,x,qmu2,stf,n,jchk)

C==   ===============================================================
C==   Set and clear package subroutine names ========================
C==   ===============================================================
            
C     ========================
      subroutine setUmsg(name)
C     ========================

C--   Set name of user subroutine

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qsnam3.inc'
      
      logical first
      save    first
      data    first /.true./

      character*(*) name
      
      character*80 subnam
      data subnam /'SETUMSG ( CHNAME )'/

C--   Check if QCDNUM is initialized
      if(first) then
        call sqcChkIni(subnam)
        first = .false.
      endif
    
      call smb_cfill(' ',usrnam3)
      len = min(imb_lenoc(name),80)
      usrnam3(1:len) = name(1:len)
      
      return
      end
            
C     ==================
      subroutine clrUmsg
C     ==================

C--   Set name of user subroutine

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qsnam3.inc'
      
      logical first
      save    first
      data    first /.true./
      
      character*80 subnam
      data subnam /'CLRUMSG'/

C--   Check if QCDNUM is initialized
      if(first) then
        call sqcChkIni(subnam)
        first = .false.
      endif

      call smb_cfill(' ',usrnam3)
      
      return
      end

C==   ===============================================================
C==   Store partition, weight calculation, dump and read ============
C==   ===============================================================
            
C     ======================================
      subroutine BookTab(w,nw,jtypes,nwords)
C     ======================================

C--   Partition the store into tables
C--
C--   Input:  w        array dimensioned nw in the calling routine
C--           nw       number of words in w
C--           jtypes   jtypes(i) = requested # of type-i tables  i = 1,...,4
C--   Output: nwords   number of words needed. < 0 if not enough space

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension w(*), jtypes(*), itypes(4)
      character*50 emsg

      character*80 subnam
      data subnam /'BOOKTAB ( W, NW, ITYPES, NWDS )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
      
C--   Copy input array to avoid changing input
      do i = 1,4
        itypes(i) = jtypes(i)
      enddo  

C--   Check user input
      call sqcIlele(subnam,'ITYPES(1)',0,itypes(1),99,' ')
      call sqcIlele(subnam,'ITYPES(2)',0,itypes(2),97,' ')
      call sqcIlele(subnam,'ITYPES(3)',0,itypes(3),99,' ')
      call sqcIlele(subnam,'ITYPES(4)',0,itypes(4),99,' ')
      if(itypes(1).eq.0 .and. itypes(2).eq.0 .and.
     +   itypes(3).eq.0 .and. itypes(4).eq.0      ) then
         call sqcErrMsg(subnam,'No tables to book')
      endif
C--   Add two type-2 scratch tables
      itypes(2) = itypes(2)+2                       
C--   Do the work
      call sqcBookTab(w,nw,itypes,nwords,ierr)
C--   No valid type encountered
      if(ierr.eq.1) then
        call sqcErrMsg(subnam,'No tables to book')
      endif
C--   Not enough storage
      if(nwords.le.0) then
        write(emsg,'(''Increase NW to at least '',I10)') -nwords 
        call sqcErrMsg(subnam,emsg)
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
      
C     ===========================
      subroutine SetWpar(w,par,n)
C     ===========================

C--   Write a set of parameters in the store
C--
C--   w        (in) : store dimensioned in the calling routine
C--   par(n)   (in) : input set of parameters
C--   n        (in) : number of parameters to be stored

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension w(*),par(*)

      character*80 subnam
      data subnam /'SETWPAR ( W, PAR, N )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check store partitioned
      if(anint(w(1)).ne.123456) then
        call sqcErrMsg(subnam,'Store not partitioned into tables')
      endif
C--   Check enough words
      niw = int(anint(w(3)))
      call sqcIlele(subnam,'N',1,n,niw,
     +              'You can increase MIW0 in qcdnum.inc')
C--   Do the work
      do i = 1,n
        w(3+i) = par(i)
      enddo  

C--   Update status bits
      call sqcSetflg(iset,idel,0)     

      return
      end
      
C     ===========================
      subroutine GetWpar(w,par,n)
C     ===========================

C--   Retrieve a set of parameters from the store
C--
C--   w     (in)  :  store dimensioned in the calling routine
C--   n     (in)  :  number of parameters to be retrieved
C--   par   (out) :  output set of parameters 

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension w(*),par(*)

      character*80 subnam
      data subnam /'GETWPAR ( W, PAR, N )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check store partitioned
      if(anint(w(1)).ne.123456) then
        call sqcErrMsg(subnam,'Store not partitioned into tables')
      endif
C--   Check enough words
      niw = int(anint(w(3)))
      call sqcIlele(subnam,'N',1,n,niw,
     +              'You can increase MIW0 in qcdnum.inc')
C--   Do the work
      do i = 1,n
        par(i) = w(3+i)
      enddo  

C--   Update status bits
      call sqcSetflg(iset,idel,0)     

      return
      end
      
C     ==================================
      subroutine TabDump(w,lun,file,key)
C     ==================================

C--   Dump store to disk
C--
C--   Input:  w        store dimensioned in the calling routine
C--           lun      logical unit number
C--           file     output file name
C--           key      key string to be written on the file

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension w(*)
      character*(*) file, key 

      character*80 subnam
      data subnam /'TABDUMP ( W, LUN, FILE, KEY )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check store partitioned
      if(anint(w(1)).ne.123456) then
        call sqcErrMsg(subnam,'Store not partitioned into tables')
      endif
C--   Open output file
      open(unit=lun,file=file,form='unformatted',status='unknown',
     +     err=500)
C--   Do the work
      call sqcDumpTab(w,lun,key,ierr)
      close(lun)
      if(ierr.ne.0) goto 501

      write(lunerr1 ,'(/'' TABDUMP: tables written to '',A/)')
     +      file 

C--   Update status bits
      call sqcSetflg(iset,idel,0)     

      return

  500 continue
C--   Open error
      call sqcErrMsg(subnam,'Cannot open output file')
      return
  501 continue
C--   Write error
      call sqcErrMsg(subnam,'Write error on output file')
      return

      end
      
     

C     =================================================
      subroutine TabRead(w,nw,lun,file,key,nwords,ierr)
C     =================================================

C--   Read store from disk
C--
C--   w       (in)  : store dimensioned in the calling routine
C--   nw      (in)  : dimension of w
C--   lun     (in)  : input logical unit number
C--   file    (in)  : input file name
C--   key     (in)  : key to match that on the file
C--   nwords  (out) : number of words read in
C--   ierr    (out) :   0  : all OK
C--                     1  : read error
C--                     2  : problem with QCDNUM version
C--                     3  : problem with file stamp
C--                     4  : x-mu2 grid not the same 

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension w(*)
      character*(*) file, key

      character*50 emsg

      character*80 subnam
      data subnam /'TABREAD ( W, NW, LUN, FILE, KEY, NWORDS, IERR )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Open input file
      open(unit=lun,file=file,form='unformatted',status='old',
     +     err=500)
C--   Do the work
      call sqcReadTab(w,nw,lun,key,nwords,ierr)
      close(lun)
C--   Not enough space
      if(ierr.eq.5) then
        write(emsg,'(''Increase NW to at least '',I10)') nwords 
        call sqcErrMsg(subnam,emsg)
      endif
      if(ierr.eq.0) then
        write(lunerr1 ,'(/'' TABREAD: tables read in from '',A/)')
     +        file
      endif
      
C--   Update status bits
      call sqcSetflg(iset,idel,0)      

      return

  500 continue
C--   Read error
      ierr = 1
      return

      end

C     ==================================
      subroutine MakeWtA(w,id,afun,achi)
C     ==================================

C--   Make weight table for regular contribution A
C--
C--   w     (in) :   store dimensioned in the calling routine
C--   id    (in) :   table identifier
C--   afun  (in) :   function declared external in the calling routine 
C--   achi  (in) :   function declared external in the calling routine 

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'

      external afun, achi

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension w(*)

      character*80 subnam
      data subnam /'MAKEWTA ( W, ID, AFUN, ACHI )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
      
C--   Check id
      ierr = iqcCheckId(w,id)
      if(ierr.gt.0) call sqcIdEmsg(subnam,'ID',id,ierr)      

C--   Do the work
      call sqcUweitA(w,id,ioy2,afun,achi,ierr)
      if(ierr.eq.1) then 
        call sqcErrMsg(subnam,
     +                'Function achi(qmu2) < 1 encountered')
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C     =======================================
      subroutine MakeWtB(w,id,bfun,achi,ndel)
C     =======================================

C--   Make weight table for singular contribution B
C--
C--   w      (in) :  store dimensioned in the calling routine
C--   id     (in) :  table identifier
C--   bfun   (in) :  function declared external in the calling routine
C--   achi   (in) :  function declared external in the calling routine
C--   ndel   (in) :  1 = no delta(1-x) contribution; otherwise yes

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'

      external bfun, achi

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension w(*)

      character*80 subnam
      data subnam  /'MAKEWTB ( W, ID, BFUN, ACHI, NODELTA )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
      
C--   Check id
      ierr = iqcCheckId(w,id)
      if(ierr.gt.0) call sqcIdEmsg(subnam,'ID',id,ierr)      

C--   Do the work (note flip ndel 0 <--> 1)
      call sqcUweitB(w,id,ioy2,bfun,achi,1-ndel,ierr)
      if(ierr.eq.1) then 
        call sqcErrMsg(subnam,
     +                'Function achi(qmu2) < 1 encountered')
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C     ============================================
      subroutine MakeWRS(w,id,rfun,sfun,achi,ndel)
C     ============================================

C--   Make weight table for contribution RS
C--
C--   w      (in) :  store dimensioned in the calling routine
C--   id     (in) :  table identifier
C--   rfun   (in) :  function declared external in the calling routine
C--   sfun   (in) :  function declared external in the calling routine 
C--   achi   (in) :  function declared external in the calling routine
C--   ndel   (in) :  1 = no delta(1-x) contribution; otherwise yes

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'

      external rfun, sfun, achi

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension w(*)

      character*80 subnam
      data subnam /'MAKEWRS ( W, ID, RFUN, SFUN, ACHI , NODELTA )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
      
C--   Check id
      ierr = iqcCheckId(w,id)
      if(ierr.gt.0) call sqcIdEmsg(subnam,'ID',id,ierr)      

C--   Do the work (note flip ndel 0 <--> 1)
      call sqcUwgtRS(w,id,ioy2,rfun,sfun,achi,1-ndel,ierr)
      if(ierr.eq.1) then 
        call sqcErrMsg(subnam,
     +                'Function achi(qmu2) < 1 encountered')
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C     ==================================
      subroutine MakeWtD(w,id,dfun,achi)
C     ==================================

C--   Make weight table for factor*delta(1-x) contribution
C--
C--   w      (in) :  store dimensioned in the calling routine
C--   id     (in) :  table identifier
C--   dfun   (in) :  function declared external in the calling routine
C--   achi   (in) :  function declared external in the calling routine


      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'

      external dfun,achi

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension w(*)

      character*80 subnam
      data subnam  /'MAKEWTD ( W, ID, DFUN, ACHI )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check id
      ierr = iqcCheckId(w,id)
      if(ierr.gt.0) call sqcIdEmsg(subnam,'ID',id,ierr)            

C--   Do the work
      call sqcUweitD(w,id,ioy2,dfun,achi,ierr)
      if(ierr.eq.1) then 
        call sqcErrMsg(subnam,
     +                'Function achi(qmu2) < 1 encountered')
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
      
C     ========================
      subroutine MakeWtX(w,id)
C     ========================

C--   Make weight table for FxF convolution
C--
C--   w    (in) :    store dimensioned in the calling routine
C--   id   (in) :    table identifier

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension w(*)

      character*80 subnam
      data subnam /'MAKEWTX ( W, ID )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
      
C--   Check id
      ierr = iqcCheckId(w,id)
      if(ierr.gt.0) call sqcIdEmsg(subnam,'ID',id,ierr)      

C--   Do the work
      call sqcUweitX(w,id,ioy2,ierr)
      if(ierr.eq.1) then 
        call sqcErrMsg(subnam,'Error condition encountered')
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
      
C     ==========================
      subroutine ScaleWt(w,c,id)
C     ==========================

C--   Multiply weight table by a constant
C--
C--   w   (in)     store dimensioned in the calling routine
C--   c   (in)     constant
C--   id  (in)     table identifier

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension w(*)

      character*80 subnam
      data subnam /'SCALEWT ( W, C, ID )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check id      
      ierr = iqcCheckId(w,id)
      if(ierr.ne.0) call sqcIdEmsg(subnam,'ID',id,ierr)

C--   Do the work
      call sqcScaleWt(w,c,id)

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C     ==================================
      subroutine CopyWgt(w,jd1,id2,iadd)
C     ==================================

C--   Add content of id1 to id2
C--
C--   w    (in) : store declared and partitioned in the calling routine
C--   jd1  (in) : input table id, < 0 address splitting function table
C--   id2  (in) : output table id in the store w
C--   iadd (in) : -1,0,1 = subtract/copy/add id1 to id2

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension w(*)

C--   iotyp(ityp_in,ityp_out) = 0 --> combination not allowed      
      dimension iotyp(4,4)
C--   ityp_in      1  2  3  4              
      data iotyp / 1, 0, 0, 0,    !ityp_out 1
     +             1, 1, 0, 0,    !ityp_out 2
     +             1, 0, 1, 0,    !ityp_out 3
     +             1, 1, 1, 1 /   !ityp_out 4

      character*10 number
      character*80 emsg
      character*80 subnam
      data subnam /'COPYWGT ( W, ID1, ID2, IADD )'/

C--   Avoid compiler warning
      ifirst = 0
C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
      
C--   Overwrite, no thank you...
      if(id2.eq.jd1) then
        call sqcErrMsg(subnam,'ID2 cannot be equal to ID1')
      endif      

C--   Check iadd
      call sqcIlele(subnam,'IADD',-1,iadd,1,' ')

C--   Check id1            
      if(jd1.lt.0) then
C--     Splitting function table      
        call sqcChekPij(jd1,id1,jset,ifirst,ierr)
        if(ierr.ne.0) call sqcIdEmsg(subnam,'ID1',jd1,ierr)
        ityp1 = id1/100
      else
C--     Table in store w      
        id1 = jd1
C--     Check id1
        ierr = iqcCheckId(w,id1)
        if(ierr.ne.0) call sqcIdEmsg(subnam,'ID1',id1,ierr)
        ityp1 = id1/100
      endif

C--   Check id2
      ierr = iqcCheckId(w,id2)
      if(ierr.gt.0) call sqcIdEmsg(subnam,'ID2',id2,ierr)
      ityp2 = id2/100
      
C--   Check input/output table types      
      if(iotyp(ityp1,ityp2) .eq. 0) then
        call smb_itoch(id2,number,lnum)
        write(emsg,'(
     +   ''ID2 = '',A,'' : table type incompatible with ID1'')'),
     +        number(1:lnum)
        call sqcErrMsg(subnam,emsg)
      endif
            
C--   Do the work
      if(jd1.lt.0) then
        call sqcCopyWt(stor7(ifirst),id1,w,id2,iadd)
      else
        call sqcCopyWt(w,id1,w,id2,iadd)
      endif    

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C     ======================================
      subroutine WcrossW(w,jda,jdb,idc,iadd)
C     ======================================

C--   Make weight table for convolution C = A cross B

C--   w    (in) : store declared and partitioned in the calling routine
C--   jda  (in) : input table id, < 0 address splitting function table
C--   jdb  (in) : input table id, < 0 address splitting function table
C--   idc  (in) : output table id in the store w
C--   iadd (in) : -1,0,1 = subtract/copy/add id1 to id2

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'
      
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension w(*)
      
C--   iotyp(ityp_in,ityp_out) = 0 --> combination not allowed      
      dimension iotyp(4,4)
C--   ityp_in      1  2  3  4              
      data iotyp / 1, 0, 0, 0,    !ityp_out 1
     +             1, 1, 0, 0,    !ityp_out 2
     +             1, 0, 1, 0,    !ityp_out 3
     +             1, 1, 1, 1 /   !ityp_out 4

      character*10 number
      character*80 emsg
      character*80 subnam
      data subnam /'WCROSSW ( W, IDA, IDB, IDC, IADD )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Overwrite, no thank you...
      if(idc.eq.jda .or. idc.eq.jdb) then
        call sqcErrMsg(subnam,'IDC cannot be equal to IDA or IDB')
      endif
      
C--   Check iadd
      call sqcIlele(subnam,'IADD',-1,iadd,1,' ')
      
C--   Check ida            
      if(jda.lt.0) then
C--     Splitting function table      
        call sqcChekPij(jda,ida,jseta,ifirsta,ierr)
        if(ierr.ne.0) call sqcIdEmsg(subnam,'IDA',jda,ierr)
        itypa = ida/100
      else
C--     Table in store w      
        ida = jda
C--     Check ida
        ierr = iqcCheckId(w,ida)
        if(ierr.ne.0) call sqcIdEmsg(subnam,'IDA',ida,ierr)
        itypa = ida/100
      endif                          

C--   Check idb            
      if(jdb.lt.0) then
C--     Splitting function table      
        call sqcChekPij(jdb,idb,jsetb,ifirstb,ierr)
        if(ierr.ne.0) call sqcIdEmsg(subnam,'IDB',jdb,ierr)
        itypb = idb/100
      else
C--     Table in store w      
        idb = jdb
C--     Check idb
        ierr = iqcCheckId(w,idb)
        if(ierr.ne.0) call sqcIdEmsg(subnam,'IDB',idb,ierr)
        itypb = idb/100
      endif                          

C--   Check idc
      ierr = iqcCheckId(w,idc)
      if(ierr.gt.0) call sqcIdEmsg(subnam,'IDC',idc,ierr)
      itypc = idc/100
      
C--   Check input/output table types      
      if(iotyp(itypa,itypc).eq.0   .or.
     +   iotyp(itypb,itypc).eq.0 ) then
        call smb_itoch(idc,number,lnum)
        write(emsg,'(
     + ''IDC = '',A,'' : table type incompatible with IDA or IDB'')'),
     +        number(1:lnum)
        call sqcErrMsg(subnam,emsg)
      endif 
            
C--   Do the work, idt = -1,-2 are scratch tables
      idt1 = -1
      idt2 = -2
      if(jda.lt.0 .and. jdb.lt.0) then
        call sqcWcrossW(stor7(ifirsta),ida,
     +                  stor7(ifirstb),idb,w,idc,idt1,idt2,iadd)
      elseif(jda.lt.0 .and. jdb.gt.0) then
        call sqcWcrossW(stor7(ifirsta),ida,
     +                  w             ,idb,w,idc,idt1,idt2,iadd)
      elseif(jda.gt.0 .and. jdb.lt.0) then
        call sqcWcrossW(w             ,ida,
     +                  stor7(ifirstb),idb,w,idc,idt1,idt2,iadd)
      else 
        call sqcWcrossW(w,ida,w,idb,w,idc,idt1,idt2,iadd)
      endif  

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
      
C     ======================================
      subroutine WtimesF(w,fun,jd1,id2,iadd)
C     ======================================

C--   Multiply id1 by fun(x,muf,mur,nf,alfas/2pi) and store result in id2
C--
C--   w    (in) :    store dimensioned in the calling routine
C--   fun  (in) :    user defined function
C--   jd1  (in) :    input table (< 0 means splitting function table) 
C--   id2  (in) :    output table identifier
C--   iadd (in) :    -1,0,1  subtract,store,add result to id2  

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'
      
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension w(*)
      
C--   iotyp(ityp_in,ityp_out) = 0 --> combination not allowed      
      dimension iotyp(4,4)
C--   ityp_in      1  2  3  4              
      data iotyp / 1, 0, 0, 0,    !ityp_out 1
     +             1, 1, 0, 0,    !ityp_out 2
     +             1, 0, 1, 0,    !ityp_out 3
     +             1, 1, 1, 1 /   !ityp_out 4

      external fun

      character*10 number
      character*80 emsg
      character*80 subnam
      data subnam /'WTIMESF ( W, FUN, ID1, ID2, IADD )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check iadd
      call sqcIlele(subnam,'IADD',-1,iadd,1,' ')
      
C--   Check id1            
      if(jd1.lt.0) then
C--     Splitting function table      
        call sqcChekPij(jd1,id1,jset,ifirst,ierr)
        if(ierr.ne.0) call sqcIdEmsg(subnam,'ID1',jd1,ierr)
        ityp1 = id1/100
      else
C--     Table in store w      
        id1 = jd1
C--     Check id1
        ierr = iqcCheckId(w,id1)
        if(ierr.ne.0) call sqcIdEmsg(subnam,'ID1',id1,ierr)
        ityp1 = id1/100
      endif
      
C--   Check id2
      ierr = iqcCheckId(w,id2)
      if(ierr.gt.0) call sqcIdEmsg(subnam,'ID2',id2,ierr)
      ityp2 = id2/100      
      
      if(iotyp(ityp1,ityp2) .eq. 0) then
        call smb_itoch(idc,number,lnum)
        write(emsg,'(
     + ''ID2 = '',A,'' : table type incompatible with ID1'')'),
     +        number(1:lnum)
        call sqcErrMsg(subnam,emsg)
      endif      
      
C--   Do the work
      if(jd1.lt.0) then
        call sqcWtimesF(fun,stor7(ifirst),id1,w,id2,iadd)
      else
        call sqcWtimesF(fun,w,id1,w,id2,iadd)
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
      
C     ==========================================
      integer function idSpfun(string,iord,jset)
C     ==========================================

C--   Returns splitting function table index encoded as -(1000*jset+idspl)
C--   Returns -1 if splitting function table not available

      implicit double precision (a-h,o-z)      
      
      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'
      
      character*(*) string
      character*3   input
      character*3   ptab(7)
C--                 1     2     3     4     5     6     7
      data ptab / 'PQQ','PQG','PGQ','PGG','PMI','PPL','PVA'/
      
      idspfun = -1
      if(jset.lt.1. .or. jset.gt.4) return
      if(mxord7(jset).eq.0)         return

      input        = '   '
      len          = min(imb_lenoc(string),3)
      input(1:len) = string(1:len)
      call smb_cltou(input)
      id = 0
      do i = 1,7
        if(input.eq.ptab(i))   id = idPij7(i,iord,jset)
      enddo
      
      if(id.eq.0) then
        idspfun = -1
      else
        idspfun = -(id+1000*jset)  
      endif

      return
      end
      
C==   ===============================================================
C==   Convolution engine ============================================
C==   ===============================================================

C     =======================================================
      double precision function FcrossK(w,idw,jset,idf,ix,iq)
C     =======================================================

C--   Convolution F cross K
C--
C--   w     (in) :   store dimensioned in the calling routine
C--   idw   (in) :   weight table identifier
C--   jset  (in) :   pdf set identifier
C--   idf   (in) :   identifier F
C--   ix    (in) :   x-grid index
C--   iq    (in) :   mu2 grid index

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension w(*)

      character*80 subnam
      data subnam /'FCROSSK ( W, IDW, ISET, IDF, IX, IQ )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check iset (should be 1-9 when evol in new format)
      call sqcIlele(subnam,'ISET',1,jset,mset0,' ')

C--   Check status bits
      call sqcChkflg(jset,ichk,subnam)
      
C--   Check idw
      ierr = iqcCheckId(w,idw)
      if(ierr.ne.0) call sqcIdEmsg(subnam,'IDW',id,ierr)      

C--   Check cuts
      ixmi = nyy2(0)+1-iymac2
      iqmi = itfiz2(izmic2)
      iqma = itfiz2(izmac2)
      call sqcIlele(subnam,'IDF',0,idf,12,' ')
      call sqcIlele(subnam,'IX',ixmi,ix,nyy2(0)+1,' ')
      call sqcIlele(subnam,'IQ',iqmi,iq,iqma,' ')

C--   Do the work
      iy = nyy2(0) + 1 - ix
      FcrossK = dqcFcrossK(w,idw,jset,idf,iy,iq)

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end    
      
C     ===========================================================
      double precision function FcrossF(w,idw,jset,ida,idb,ix,iq)
C     ===========================================================

C--   Convolution Fa cross Fb
C--
C--   w     (in) :   store dimensioned in the calling routine
C--   idw   (in) :   weight table identifier
C--   jset  (in) :   pdf set identifier
C--   ida   (in) :   identifier Fa
C--   idb   (in) :   identifier Fb
C--   ix    (in) :   x-grid index
C--   iq    (in) :   mu2 grid index

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension w(*)

      character*80 subnam
      data subnam /'FCROSSF ( W, IDW, ISET, IDA, IDB, IX, IQ )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check iset (should be 1-9 when evol in new format)
      call sqcIlele(subnam,'ISET',1,jset,mset0,' ')

C--   Check status bits
      call sqcChkflg(jset,ichk,subnam)

C--   Check idw
      ierr = iqcCheckId(w,idw)
      if(ierr.ne.0) call sqcIdEmsg(subnam,'IDW',id,ierr)
      
C--   Check ranges
      call sqcIlele(subnam,'IDA',0,ida,12,' ')
      call sqcIlele(subnam,'IDB',0,idb,12,' ')

C--   Check cuts
      ixmi = nyy2(0)+1-iymac2
      iqmi = itfiz2(izmic2)
      iqma = itfiz2(izmac2)      
      call sqcIlele(subnam,'IX',ixmi,ix,nyy2(0)+1,' ')
      call sqcIlele(subnam,'IQ',iqmi,iq,iqma,' ')

C--   Do the work
      iy = nyy2(0) + 1 - ix
      FcrossF = dqcFcrossF(w,idw,jset,ida,idb,iy,iq)

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end              
    
C     ============================
      integer function Nflavor(iq)
C     ============================

C--   Returns number of flavors

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qmaps8.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      character*80 subnam
      data subnam  /'NFLAVOR ( IQ )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check range
      call sqcIlele(subnam,'IQ',1,iq,ntt2,' ')
C--   Make sure Nfmap exists 
      if(.not.Lnfmap8) call sqcNfTab(0)     
C--   Do the work
      iz       = izfit2(iq)
      nflavor  = nffiz2(iz) 
C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
    
C     ============================================
      double precision function GetAlfN(jq,n,ierr)
C     ============================================

C--   Returns (as/2pi)^n
C--
C--   ierr = 1  iq close or below Lambda^2
C--          2  iq outside grid

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qmaps8.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      character*80 subnam
      data subnam /'GETALFN ( IQ, N, IERR )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Make sure Nfmap exists 
      if(.not.Lnfmap8) call sqcNfTab(0)      
C--   Make sure astable exists
      if(.not.Lastab8) call sqcAlfTab(iord6)    
C--   Check range
      call sqcIlele(subnam,'N',-2,n,20,' ')
C--   Check iq inside grid 
      iq = abs(jq)
      if(iq.lt.1 .or. iq.gt.ntt2) then
        getalfn = qnull6
        ierr    = 2
        return
      endif
C--   Not below Lambda, thank you      
      if(iq.lt.itlow8) then
        getalfn = qnull6
        ierr    = 1
        return
      endif  
C--   Find out which iz index to take
      iz = izfit2(iq)
C--   Take iz-1 at threshold if jq preceeded by minus sign
      if(jq.lt.0 .and. iz.ne.1) then
        if(nffiz2(iz-1).eq.nffiz2(iz)-1) iz=iz-1
      endif
C--   Do the work
      ierr = 0
      if(n.eq.0)                      then
        getalfn = 1.D0
      elseif(n.lt.0)                  then
        getalfn = antab8(iz,n)  
      elseif(n.gt.0 .and. n.le.iord6) then
        getalfn = antab8(iz,n)
      else  
        getalfn = (antab8(iz,0))**n
      endif
 
C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
    
C     ================================
      subroutine EfromQQ(qvec,evec,nf)
C     ================================

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

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpdfs7.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension qvec(-6:6), evec(12)

      character*80 subnam
      data subnam /'EFROMQQ ( QVEC, EVEC, NF )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Do the work      
      call sqcEfromQQ(qvec,evec,nf)
C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C     ================================
      subroutine QQfromE(evec,qvec,nf)
C     ================================

C--   Transform si/ns basis to flavor basis
C--
C--   |epm>  0  1  2  3  4  5  6  7  8  9 10 11 12
C--          g  si 2+ 3+ 4+ 5+ 6+ v  2- 3- 4- 5- 6-
C--
C--   |qpm>  0  1  2  3  4  5  6  7  8  9 10 11 12
C--          g  d+ u+ s+ c+ b+ t+ d- u- s- c- b- t-
C--
C--         -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
C--   qvec  tb bb cb sb ub db  g  d  u  s  c  b  t

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpdfs7.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension evec(12), qpm(12), qvec(-6:6)

      character*80 subnam
      data subnam /'QQFROME ( EVEC, QVEC, NF )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Transform  (eq 2.29)
      do i = 1,nf
        biplu = 0.D0
        bimin = 0.D0
        do j = 1,nf
          biplu = biplu + evec(j  )*umateq(j,i)
          bimin = bimin + evec(6+j)*umateq(j,i)
        enddo
        qpm(i  ) = biplu
        qpm(i+6) = bimin
      enddo
C--   Get q,qbar coeff from qpm coeff  (eq 2.20)
      do i = 0,12
        qvec(i) = 0.D0
      enddo
      do i = 1,nf
        qvec( i) = qpm(i) + qpm(6+i)
        qvec(-i) = qpm(i) - qpm(6+i)
      enddo
 
C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end

C==   ===============================================================
C==   Interpolation =================================================
C==   ===============================================================

C     ===========================================
      subroutine StfunXq(stfun,x,qmu2,stf,n,jchk)
C     ===========================================

C--   Interpolate stfun(ix,iq)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

      external stfun
      dimension x(*),qmu2(*),stf(*)

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      character*80 subnam
      data subnam /'STFUNXQ ( STFUN, X, QMU2, STF, N, ICHK )'/

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
C--   Do the work
      call sqcStInterp(stfun,x,qmu2,stf,n,ierr)
C--   Again, check on number of interpolations
      if(ierr.eq.-1) then
        call sqcErrMsg(subnam,
     + 'Too many interpolations: increase mpt0 in qcdnum.inc')
      endif
C--   Outside grid
      if(jchk.ne.0 .and. ierr.eq.1) then
        call sqcErrMsg(subnam,'At least one x or mu2 outside grid')
      endif
C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
