
C--   This is the file usrwgt.f containing the qcdnum weight routines

C--   subroutine fillwt(iselect,idmin,idmax,nwlast)
C--   subroutine fillwc(usub,idmin,idmax,nwlast)
C--   subroutine dmpwgt(jset,lun,file)
C--   subroutine readwt(lun,file,idmin,idmax,nwlast,ierr)
C--   subroutine nwused(nwtot,nwuse,nwtab)

C==   ===============================================================
C==   Weight calculation, dump and read =============================
C==   ===============================================================
            
C     =============================================
      subroutine fillwt(iselect,idmin,idmax,nwlast)
C     =============================================

C--   Fill standard weight tables
C--
C--   iselect (in)     1 = Unpolarised 
C--                    2 = Polarised
C--                    3 = Fragmentation functions (timelike)
C--                 else = Unpolarised
C--   idmin   (out)        first pdf identifier (always 0 = gluon)
C--   idmax   (out)        last  pdf identifier
C--   nwlast  (out)        last word used in the store < 0 no space

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'qmaps8.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)
      
      character*60 emsg
      character*10 etxt
      
      external sqcFilWU, sqcFilWP, sqcFilWF

      character*80 subnam
      data subnam /'FILLWT ( ITYPE, IDMIN, IDMAX, NWDS )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
      
C--   Yes/no re-initialize, thats the question
      if(.not.Lwtini8) then
        call sqcIniWt
      endif 
      
      if(iselect.eq.2) then
      
C--     Calculate polarised weight tables
        write(lunerr1,'(/
     +        '' FILLWT: start polarised weight calculations'')')
        write(lunerr1,'( '' Subgrids'',I5, 
     +        '' Subgrid points'',100I5)'),nyg2,(nyy2(i),i=1,nyg2)
        call sqcFilWt(sqcFilWP,2,nwlast,ierr)
        jselect = 2
        
      elseif(iselect.eq.3) then
      
C--     Calculate timelike weight tables (fragmentation functions)
        write(lunerr1,'(/
     +        '' FILLWT: start fragmentation weight calculations'')')
        write(lunerr1,'( '' Subgrids'',I5, 
     +        '' Subgrid points'',100I5)'),nyg2,(nyy2(i),i=1,nyg2)
        call sqcFilWt(sqcFilWF,3,nwlast,ierr)
        jselect = 3
      
      else
            
C--     Calculate unpolarised weight tables (default)
        write(lunerr1,'(/
     +        '' FILLWT: start unpolarised weight calculations'')')
        write(lunerr1,'( '' Subgrids'',I5, 
     +        '' Subgrid points'',100I5)'),nyg2,(nyy2(i),i=1,nyg2)
        call sqcFilWt(sqcFilWU,1,nwlast,ierr)
        jselect = 1
      
      endif
      
C--   Weight tables already existed so inform that nothing has been done
      if(ierr.eq.-1) then
        write(lunerr1,'( 
     +           '' Tables already exist --> nothing to be done'')')
      endif       
        
C--   Make sure that there is enough space to hold first word of next set
      nwfirst = abs(nwlast)+1
      if(nwfirst.gt.nwf0) then
        call smb_itoch(nwfirst,etxt,ltxt)
        write(emsg,'(''Need at least '',A,
     &    '' words --> increase NWF0 '',
     &    ''in qcdnum.inc'')') etxt(1:ltxt)
        call sqcErrMsg(subnam,emsg)
      endif
      write(lunerr1,'('' FILLWT: weight calculations completed''/)')
      
C--   Idmin and Idmax
      idmin  = 0
      idmax  = idmx7

C--   NFmap, splines etc
      if(.not.Lnfmap8) call sqcNFtab(0)

C--   Update status bits
      call sqcSetflg(iset,idel,jselect)
      
      return
      end
      
C     ==========================================
      subroutine fillwc(usub,idmin,idmax,nwlast)
C     ==========================================

C--   Fill standard weight tables
C--
C--   usub    (in)         user defined subroutine     
C--   idmin   (out)        first pdf identifier (always 0 = gluon)
C--   idmax   (out)        last  pdf identifier
C--   nwlast  (out)        last word used in the store < 0 no space

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'qmaps8.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)
      
      character*60 emsg
      character*10 etxt
      
      external usub

      character*80 subnam
      data subnam /'FILLWC ( MYSUB, IDMIN, IDMAX, NWDS )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
      
C--   Yes/no re-initialize, thats the question
      if(.not.Lwtini8) then
        call sqcIniWt
      endif 
            
C--   Calculate custom weight tables
      write(lunerr1,'(/
     +      '' FILLWC: start custom weight calculations'')')
      write(lunerr1,'( '' Subgrids'',I5, 
     +      '' Subgrid points'',100I5)'),nyg2,(nyy2(i),i=1,nyg2)
      call sqcFilWt(usub,4,nwlast,ierr)
      
C--   Custom weight tables already exist so thats an error
      if(ierr.eq.-1) then
        call sqcErrMsg(subnam,'Custom tables already exist')
      endif
C--   Mxord not in range 1-3
      if(ierr.eq.-2) then
        call sqcErrMsg(subnam,'Maxord not in range [1-3]')
      endif             
        
C--   Make sure that there is enough space to hold first word of next set
      nwfirst = abs(nwlast)+1
      if(nwfirst.gt.nwf0) then
        call smb_itoch(nwfirst,etxt,ltxt)
        write(emsg,'(''Need at least '',A,
     &    '' words --> increase NWF0 '',
     &    ''in qcdnum.inc'')') etxt(1:ltxt)
        call sqcErrMsg(subnam,emsg)
      endif
      write(lunerr1,'('' FILLWC: weight calculations completed''/)')
      
C--   Idmin and Idmax
      idmin  = 0
      idmax  = idmx7

C--   NFmap, splines etc
      if(.not.Lnfmap8) call sqcNFtab(0)

C--   Update status bits
      call sqcSetflg(iset,idel,4)
      
      return
      end      
      
C     ================================
      subroutine dmpwgt(jset,lun,file)
C     ================================

C--   Dump weights to a disk file (unformatted write)
C--
C--   jset   (in)   0     = dump all pdf sets except custom
C--                 other = dump jset

      implicit double precision (a-h,o-z)

      character*(*) file

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)
      
      character*5 etxt
      character*13 txt(4)
C--               1234567890123      
      data txt / 'unpolarised  ',
     +           'polarised    ',
     +           'fragmentation',
     +           'custom       '  /
      
      character*80 subnam
      data subnam /'DMPWGT ( ITYPE, LUN, FILE )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif
C--   Check qcdnum initialized
      call sqcChkIni(subnam)
C--   Check logical unit number
      if(lun.le.0 .or. lun.eq.6) then
        call smb_itoch(lun,etxt,ltxt)
        call sqcErrMsg(subnam,
     +    'Invalid logical unit number lun = '//etxt(1:ltxt))      
      endif
C--   Check jset
      call sqcIlele(subnam,'ISET',0,jset,4,' ')      
      if(jset.eq.0) then
        write(lunerr1 ,'(/
     +               '' DUMPWT: dump all standard weight tables'')')
      else
C--     Check status bits
        call sqcChkflg(jset,ichk,subnam)      
        leng = imb_lenoc(txt(jset))
        write(lunerr1 ,'(/'' DUMPWT: dump '',A,'' weight tables'')')
     +  txt(jset)(1:leng)
      endif
C--   Now dumpit 
      open(unit=lun,file=file,form='unformatted',
     +     status='unknown',err=500)
      call sqcDumpWt(lun,jset,' ',ierr)
      close(lun)
      if(ierr.eq.1) goto 501
      if(ierr.eq.2) goto 502

      write(lunerr1 ,'(''         weights written to '',A/)') file
     
C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return

  500 continue
C--   Open error
      call sqcErrMsg(subnam,'Cannot open output weight file')
      return
  501 continue
C--   Write error
      call sqcErrMsg(subnam,'Write error on output weight file')
      return
  502 continue
C--   No weight tables available
      call sqcErrMsg(subnam,'No weight tables available')
      return      
  
      end

C     ===================================================
      subroutine readwt(lun,file,idmin,idmax,nwlast,ierr)
C     ===================================================

C--   Read weights from a disk file (unformatted read) and partition
C--   the store (that is, partiton weight tables and also pdf tables)
C--
C--   lun     (in)   nput logical unit number
C--   file    (in)   input file name
C--   idmin   (out)  index of first pdf table (always 0 = gluon)
C--   idmax   (out)  index of last  pdf table
C--   nwlast  (out)  last word occupied in store < 0 not enough space 
C--   ierr    (out)  0 = all OK
C--                  1 = read error
C--                  2 = problem with QCDNUM version
C--                  3 = key mismatch
C--                  4 = x-mu2 grid not the same

      implicit double precision (a-h,o-z)

      character*(*) file

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'qmaps8.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)
      
      character*60 emsg
      character*10 etxt
      character*13 txt(4)
C--               1234567890123      
      data txt / 'unpolarised  ',
     +           'polarised    ',
     +           'fragmentation',
     +           'custom       '  /
      
      dimension iread(mset0)

      character*80 subnam
      data subnam /'READWT ( LUN, FNAM, IDMI, IDMA, NWDS, IERR )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)
      
C--   Yes/no re-initialize, thats the question
      if(.not.Lwtini8) then
        call sqcIniWt
      endif      
      
C--   Read weight tables
      write(lunerr1 ,'(/'' READWT: open file '',A)')
     +      file
      open(unit=lun,file=file,form='unformatted',status='old',
     +     err=500)
      call sqcReadWt(lun,' ',nwlast,iread,ierr)
      close(lun)
      
C--   First word of next set (next call to fillwt or readwt)
      nwfirst = abs(nwlast)+1    
      
C--   Not enough space
      if(nwfirst.gt.nwf0) then
        call smb_itoch(nwfirst,etxt,ltxt)
        write(emsg,'(''Need at least '',A,
     &  '' words --> increase NWF0 '',
     &  ''in qcdnum.inc'')') etxt(1:ltxt) 
        call sqcErrMsg(subnam,emsg)
      endif 

C--   Idmin and Idmax
      idmin  = 0
      idmax  = idmx7     

C--   Update status bits
      do i = 1,mset0
        if(iread(i).eq.1) then
          call sqcSetflg(iset,idel,i)
          leng = imb_lenoc(txt(i)) 
          write(lunerr1,'(''         read '', A, 
     +     '' weight tables'')') txt(i)(1:leng)
        elseif(iread(i).eq.-1) then
          leng = imb_lenoc(txt(i)) 
          write(lunerr1,'(9X,A, '' tables already exist'', 
     +     '' --> nothing done'')') txt(i)(1:leng)
        endif  
      enddo
      write(lunerr1,'(/)')
      
      return

  500 continue
C--   Read error
      ierr = 1
      return

      end
      
C     ====================================      
      subroutine nwused(nwtot,nwuse,nwtab)      
C     ====================================

C--   Get store parameters

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      character*80 subnam
      data subnam /'NWUSED ( NWTOT, NWUSE, NWTAB )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first  = .false.
      endif

C--   Check status bits
      call sqcChkflg(1,ichk,subnam)

      nwtot = nwf0
      nwuse = nwlast7      
      nwtab = leng7
      
C--   Update status bits
      call sqcSetflg(iset,idel,0)
      
      return
      end
      
