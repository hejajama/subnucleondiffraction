
C--   This is the file usrini.f containing initialization routines

C--   subroutine qcinit(lun,fname)
C--   subroutine setlun(lun,fname)
C--   subroutine setval(chopt,dval)
C--   subroutine getval(chopt,dval)
C--   subroutine setint(chopt,ival)
C--   subroutine getint(chopt,ival)
C--   subroutine setord(iord)
C--   subroutine getord(iord)
C--   subroutine setalf(as,22)
C--   subroutine getalf(as,r2)
C--   subroutine setcbt(nfix,iqc,iqb,iqt) 
C--   subroutine getcbt(nfix,iqc,iqb,iqt) 
C--   subroutine setabr(ar,br)
C--   subroutine getabr(ar,br)

C==   ===============================================================
C==   Initialization ================================================
C==   ===============================================================

C     ============================
      subroutine qcinit(lun,fname)
C     ============================

C--   Initialize qcdnum

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qvers1.inc'
      include 'qsflg4.inc'
      include 'qibit4.inc'

      character*(*) fname

C--   Error; -6 is allowed for std output without banner
      if(lun.le.0.and.lun.ne.-6) goto 500

C--   Version   12345678 
      ivers1 =  170005
      cvers1 = '17-00/05'
      cdate1 = '10-04-12'
      
C--   Initialize qcdnum status
      do j = 1,mset0
        do i = 1,mbp0
          istat4(i,j) = 0
        enddo
      enddo  

C--   Set initialization flag
      iniflg4 = 123456
      
C--   Assign status bits
      call sqcBitIni

C--   Initialize constants
      call sqcIniCns

C--   Setup transformation from udscbt basis to si/ns basis
      call sqcPdfMat
      
C--   Initialize weight tables and pdf sets      
      call sqcIniWt

C--   Set a few status bits
      do i = 1,mset0
        call sqcSetbit(ibinit4,istat4(1,i),mbp0)   !initialization done
      enddo  

C--   Set output stream and print banner
      call sqcSetLun(abs(lun),fname)
C--   No banner if lun = -6
      if(lun.ne.-6) call sqcBanner(lunerr1)
      call sqcReftoo(lunerr1)

      return

 500  continue

C--   Error message
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in QCINIT ( LUN, FNAME ) ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')
      write(lunerr1,*) 'LUN = ',lun,' should be positive'

      stop
      end
      
C     ============================
      subroutine setlun(lun,fname)
C     ============================

C--   (Re)Set logical unit number and output file name.

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'

      logical first
      save    first
      data    first /.true./

      character*(*) fname

      character*80 subnam
      data subnam /'SETLUN ( LUN, FNAME )'/

C--   Check if QCDNUM is initialized
      if(first) then
        call sqcChkIni(subnam)
        first = .false.
      endif
C--   Check if input is in allowed range
      call sqcIlele(subnam,'LUN',1,lun,99,
     +             'LUN should be between 1 and 99')

C--   Do the work
      call sqcSetLun(lun,fname)

      return
      end

C==   ===============================================================
C==   Routines to set and get parameters ============================
C==   ===============================================================

C     =============================
      subroutine setval(chopt,dval)
C     =============================

C--   Set double precision value in /qpars6/

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
      include 'qmaps8.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      character*(*) chopt
      character*4   opt
      character*20  message
      data message /'    : Unknown option'/

      character*80 subnam
      data subnam /'SETVAL ( CHOPT, DVAL )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first  = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Use mbutil for character string manipulations
      len        = min(imb_lenoc(chopt),4)
C--   Avoid changing the argument so make a copy 
      opt(1:len) = chopt(1:len)
C--   OK now convert to upper case
      call smb_cltou(opt)
C--   Do the work
      if    (opt(1:len).eq.'EPSI') then
        call sqcDlele(subnam,'EPSI',1.D-10,dval,1.D-4,' ')
        aepsi6 = dval
        repsi6 = dval
      elseif(opt(1:len).eq.'EPSG') then
        call sqcDlele(subnam,'EPSG',1.D-9,dval,1.D-1,' ')
        gepsi6 = dval
      elseif(opt(1:len).eq.'ELIM') then
        call sqcDlele(subnam,'ELIM',-1.D10,dval,1.D10,' ')
        dflim6 = dval  
      elseif(opt(1:len).eq.'ALIM') then
        call sqcDlele(subnam,'ALIM',1.D-10,dval,1.D10,' ')
        aslim6 = dval
C--     Invalidate alphas table
        Lastab8 = .false.
      elseif(opt(1:len).eq.'QMIN') then
        call sqcDlele(subnam,'QMIN',1.D-1,dval,qlimu6,' ')
        qlimd6 = dval
      elseif(opt(1:len).eq.'QMAX') then
        call sqcDlele(subnam,'QMAX',qlimd6,dval,1.D11,' ')
        qlimu6 = dval
      elseif(opt(1:len).eq.'NULL') then
        qnull6 = dval
      else
        message(1:len) = opt(1:len)
        call sqcErrMsg(subnam,message)
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
      
C     =============================
      subroutine getval(chopt,dval)
C     =============================

C--   Get double precision value from /qpars6/

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      logical first
      save    first
      data    first /.true./

      character*(*) chopt
      character*4   opt
      character*20  message
      data message /'    : Unknown option'/

      character*80 subnam
      data subnam /'GETVAL ( CHOPT, DVAL )'/

C--   Check if QCDNUM is initialized
      if(first) then
        call sqcChkIni(subnam)
        first = .false.
      endif

C--   Use mbutil for character string manipulations
      len        = min(imb_lenoc(chopt),4)
C--   Avoid changing the argument so make a copy 
      opt(1:len) = chopt(1:len)
C--   OK now convert to upper case
      call smb_cltou(opt)
C--   Do the work
      if    (opt(1:len).eq.'EPSI') then
        dval = repsi6
      elseif(opt(1:len).eq.'EPSG') then
        dval = gepsi6
      elseif(opt(1:len).eq.'ELIM') then
        dval = dflim6  
      elseif(opt(1:len).eq.'ALIM') then
        dval = aslim6
      elseif(opt(1:len).eq.'QMIN') then
        dval = qlimd6
      elseif(opt(1:len).eq.'QMAX') then
        dval = qlimu6
      elseif(opt(1:len).eq.'NULL') then
        dval = qnull6
      else
        message(1:len) = opt(1:len)
        call sqcErrMsg(subnam,message)
      endif

      return
      end

C     =============================
      subroutine setint(chopt,ival)
C     =============================

C--   Set integer value in /qpars6/

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qsflg4.inc'
      include 'qibit4.inc'
      include 'qpars6.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      character*(*) chopt
      character*4   opt
      character*20  message
      data message /'    : Unknown option'/

      character*80 subnam
      data subnam /'SETINT ( CHOPT, IVAL )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first  = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Use mbutil for character string manipulations
      len        = min(imb_lenoc(chopt),4)
C--   Avoid changing the argument so make a copy 
      opt(1:len) = chopt(1:len)
C--   OK now convert to upper case
      call smb_cltou(opt)
C--   Do the work
      if    (opt(1:len).eq.'ITER') then
        call sqcIlele(subnam,'ITER',-9999,ival,5,' ')
        niter6 = ival
      elseif(opt(1:len).eq.'NTAB') then
        call sqcIlele(subnam,'NTAB',1,ival,mbf0,
     +   'You can increase mbf0 in qcdnum.inc and recompile')
        nscratch6 = ival
      else
        message(1:len) = opt(1:len)
        call sqcErrMsg(subnam,message)
      endif

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
      
C     =============================
      subroutine getint(chopt,ival)
C     =============================

C--   Get integer value from /qpars6/, /qluns1/ or /qcdnum/

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qvers1.inc'
      include 'qluns1.inc'
      include 'qsflg4.inc'
      include 'qpars6.inc'

      logical first
      save    first
      data    first /.true./

      character*(*) chopt
      character*4   opt
      character*20  message
      data message /'    : Unknown option'/

      character*80 subnam
      data subnam /'GETINT ( CHOPT, IVAL )'/
      
C--   Use mbutil for character string manipulations
      len        = min(imb_lenoc(chopt),4)
C--   Avoid changing the argument so make a copy 
      opt(1:len) = chopt(1:len)
C--   OK now convert to upper case
      call smb_cltou(opt)

C--   Set version number to zero if qcdnum is not initialised
      if(opt(1:len).eq.'VERS' .and. iniflg4.ne.123456) then
        ival = 0
        return
      endif  

C--   Check if QCDNUM is initialized
      if(first) then
        call sqcChkIni(subnam)
        first = .false.
      endif

C--   Do the work
      if    (opt(1:len).eq.'ITER') then
        ival = niter6
      elseif(opt(1:len).eq.'LUNQ') then
        ival = lunerr1
      elseif(opt(1:len).eq.'NTAB') then
        ival = nscratch6
      elseif(opt(1:len).eq.'MXG0') then  
        ival = mxg0
      elseif(opt(1:len).eq.'MXX0') then   
        ival = mxx0
      elseif(opt(1:len).eq.'MQQ0') then   
        ival = mqq0
      elseif(opt(1:len).eq.'MPT0') then   
        ival = mpt0
      elseif(opt(1:len).eq.'MIW0') then   
        ival = miw0
      elseif(opt(1:len).eq.'MBF0') then   
        ival = mbf0
      elseif(opt(1:len).eq.'NWF0') then   
        ival = nwf0
      elseif(opt(1:len).eq.'VERS') then   
        ival = ivers1
      else
        message(1:len) = opt(1:len)
        call sqcErrMsg(subnam,message)
      endif

      return
      end

C     =======================
      subroutine setord(iord)
C     =======================

C--   Set order 1,2,3 = LO, NLO, NNLO.

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
      include 'qpdfs7.inc'
      include 'qmaps8.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      character*80 subnam
      data subnam /'SETORD ( IORD )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first  = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check user input
      call sqcIlele(subnam,'IORD',1,iord,3,' ')
C--   Do the work
      iord6   = iord

C--   Invalidate alphas table
      Lastab8 = .false.

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
      
C     =======================
      subroutine getord(iord)
C     =======================

C--   Get current value of iord = 1,2,3 for LO, NLO, NNLO.

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      logical first
      save    first
      data    first /.true./

      character*80 subnam
      data subnam /'GETORD ( IORD )'/

C--   Check if QCDNUM is initialized
      if(first) then
        call sqcChkIni(subnam)
        first = .false.
      endif

C--   Do the work
      iord = iord6

      return
      end

C     ========================
      subroutine setalf(as,r2)
C     ========================

C--   Set input value of alpha_s(r2)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
      include 'qpdfs7.inc'
      include 'qmaps8.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      character*80 subnam
      data subnam /'SETALF ( AS, R2 )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first  = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

C--   Check user input
      call sqcDlele(subnam,'AS',1.D-5,as,aslim6,
     + 'Remark: the upper AS limit can be changed by a call to SETVAL')
      call sqcDlele(subnam,'R2',qlimd6,abs(r2),qlimu6,
     + 'Remark: these R2 limits can be changed by a call to SETVAL')

C--   Do the work
      alfq06  = as
      q0alf6  = r2

C--   Invalidate alphas table
      Lastab8 = .false.

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
      
C     ========================
      subroutine getalf(as,r2)
C     ========================

C--   Get input value of alpha_s(R2)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      logical first
      save    first
      data    first /.true./

      character*80 subnam
      data subnam /'GETALF ( AS, R2 )'/

C--   Check if QCDNUM is initialized
      if(first) then
        call sqcChkIni(subnam)
        first = .false.
      endif

C--   Do the work
      as = alfq06
      r2 = q0alf6

      return
      end

C     ===================================
      subroutine setcbt(nfix,iqc,iqb,iqt)
C     ===================================

C--   If 3 .le. nfix .le. 6 set FFNS, otherwise VFNS with
C--   c, b and t thresholds defined on the factorization scale.

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

      character*80 subnam
      data subnam /'SETCBT ( NFIX, IQC, IQB, IQT )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first   = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

      if(nfix.ge.3 .and. nfix.le.6) then
C--     FFNS         
        call sqcThrFFNS(nfix)
      else
C--     VFNS
        if((iqc.le.0 .or. iqc.gt.ntt2) .and.
     +     (iqb.le.0 .or. iqb.gt.ntt2) .and.
     +     (iqt.le.0 .or. iqt.gt.ntt2))        then
C--       OK to have iqc, iqb and iqt outside grid
          icc = 0
          ibb = 0
          itt = 0
        else
C--       Check charm
          icc = iqc
          call sqcIlele(subnam,'IQC',1,icc,ntt2,
     +    'Charm threshold should be inside range of Q2 grid')
          call sqcIlele(subnam,'IQC',2,icc,ntt2,
     +    'Charm threshold cannot be at lowest Q2 grid point')
          if((iqb.le.0 .or. iqb.gt.ntt2) .and.
     +       (iqt.le.0 .or. iqt.gt.ntt2))        then
C--         OK to have both iqb and iqt outside grid
            ibb = 0
            itt = 0
          else
C--         Check bottom
            ibb = iqb
            call sqcIlele(subnam,'IQB',iqc+2,ibb,ntt2,
     +      'IQB must be at least IQC+2')
            if((iqt.le.0 .or. iqt.gt.ntt2))      then
C--           OK to have iqt outside grid
              itt = 0
            else
C--           Check top
              itt = iqt
              call sqcIlele(subnam,'IQT',iqb+2,itt,ntt2,
     +        'IQT must be at least IQB+2')
            endif
          endif
        endif
C--     Seems OK...
        call sqcThrVFNS(icc,ibb,itt)
      endif

C--   Invalidate flavour map and Tsplines
      Lnfmap8 = .false.
C--   Invalidate alphas map
      Lastab8 = .false.

C--   Update status bits
      call sqcSetflg(iset,idel,0)

      return
      end
      
C     ===================================      
      subroutine mixfns(nfix,r2c,r2b,r2t)
C     =================================== 

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qpars6.inc'
      include 'qmaps8.inc'
      
            logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      character*80 subnam
      data subnam /'MIXFNS ( NFIX, R2C, R2B, R2T )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first   = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
      
C--   Check nfix
      call sqcIlele(subnam,'NFIX',3,nfix,6,' ')       
      

      if    (r2c.gt.0.D0) then
C--     r2c is lowest in range --> check r2c,b,t in ascending order      
        if(r2b.lt.1.01*r2c) 
     +   call sqcErrMsg(subnam,'r2b smaller than 1.01*r2c')
        if(r2t.lt.1.01*r2b) 
     +   call sqcErrMsg(subnam,'r2t smaller than 1.01*r2b')
      elseif(r2b.gt.0.D0) then
C--     r2b is lowest in range --> check r2b,t in ascending order      
         if(r2t.lt.1.01*r2b) 
     +   call sqcErrMsg(subnam,'r2t smaller than 1.01*r2b')
      endif        
     
C--   Seems OK...
      call sqcThrMFNS(nfix,r2c,r2b,r2t)
     
C--   Invalidate flavour map and Tsplines
      Lnfmap8 = .false.
C--   Invalidate alphas map
      Lastab8 = .false.

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end                                            

C     ===================================
      subroutine getcbt(nfix,q2c,q2b,q2t)
C     ===================================

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      character*80 subnam
      data subnam /'GETCBT ( NFIX, Q2C, Q2B, Q2T )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first   = .false.
      endif
      
C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
      
C--   Do the work
      if(Lmfns6) then
C--     MFNS            
        q2c  = rthrs6(4)
        q2b  = rthrs6(5)
        q2t  = rthrs6(6)
        nfix = -nfix6
      else
C--     FFNS or VFNS            
        q2c  = qthrs6(4)
        q2b  = qthrs6(5)
        q2t  = qthrs6(6)
        nfix = nfix6
      endif  
        
      return
      end

C     ========================
      subroutine setabr(ar,br)
C     ========================

C--   Set mu2r = ar*mu2f + br

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
      include 'qpdfs7.inc'
      include 'qmaps8.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      character*80 subnam
      data subnam /'SETABR ( AR, BR )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first   = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
C--   Check user input
      call sqcDlele(subnam,'AR',1.D-2,ar,1.D2,' ')
      call sqcDlele(subnam,'BR',-1.D2,br,1.D2,' ')
C--   Do the work
      aar6 = ar
      bbr6 = br
C--   Invalidate alphas table
      Lastab8 = .false.

C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
      
C     ========================
      subroutine getabr(ar,br)
C     ========================

C--   Get current values of Q2charm, Q2bottom and Q2top.

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      logical first
      save    first
      data    first /.true./

      character*80 subnam
      data subnam /'GETABR ( AR, BR )'/

C--   Check if QCDNUM is initialized
      if(first) then
        call sqcChkIni(subnam)
        first = .false.
      endif

C--   Do the work
      ar = aar6
      br = bbr6

      return
      end
