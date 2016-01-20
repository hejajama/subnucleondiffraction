
C--   This is the file qcdwgtq.f with Pij and Aij filling routines

C--   integer function iqcPijklm(iy,it,nf,ig,id)
C--   subroutine sqcIniWt
C--   subroutine sqcFilWt(filit,iset,nwlast,ierr)
C--   subroutine sqcFilWU(w,nw,nwords,idpij,mxord,ierr)
C--   subroutine sqcFilWP(w,nw,nwords,idpij,mxord,ierr)
C--   subroutine sqcFilWF(w,nw,nwords,idpij,mxord,ierr)
C--   subroutine sqcReadPij(w,nw,nwords,idpij,mxord,ierr)
C--   subroutine sqcDumpWt(lun,iset,key,ierr)
C--   subroutine sqcDumpW(lun,w,nwords,idpij,mxord,ierr)
C--   subroutine sqcReadWt(lun,key,nwlast,iread,ierr)

C     ==========================================
      integer function iqcPijklm(iy,it,nf,ig,id)
C     ==========================================

C--   Pij address function

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qstor7.inc'
    
      iqcPijklm =
     +    iqcWaddr(stor7(ifirst7),iy,it,nf,ig,id) + ifirst7 - 1

      return
      end
      
C     ===================      
      subroutine sqcIniWt
C     =================== 

C--   Initializes weight tables

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qstor7.inc'
      include 'qmaps8.inc'
      
C--   Default position of current set of tables      
      ifirst7 = 1
C--   Default pointer
      ipoint7 = 1            
C--   Table bookkeeping arrays      
      do k = 1,mset0
        do j = 1,3
          do i = 1,mpp0
            idPij7(i,j,k) = 0
          enddo
        enddo
      enddo
      do i = 1,maa0
        idAij7(i) = 0
      enddo
      do j = 2,3
        do i = 1,mset0
          ifst7(i,j) = 0
          ilst7(i,j) = 0
        enddo
      enddo    
      do i = 0,mset0
        knul7(i)  = 0
        mxord7(i) = 0
      enddo      
      
C--   Flag weight tables initialised      
      Lwtini8 = .true.
      
      return
      end                                        
          
C     ===========================================
      subroutine sqcFilWt(filit,iset,nwlast,ierr)
C     ===========================================

C--   Fill weight tables and book associated pdf tables

C--   filit    (in)   subroutine declared external in the calling routine
C--   iset     (in)   1=unpol, 2=pol, 3=frag, 4=custom 
C--   nwlast   (out)  last  word used < 0 not enough space
C--   ierr     (out)  -1 iset already exisits --> do nothing
C--                   -2 mxord not in range [1-3]

C--   The subroutine filit should have the following syntax
C--
C--        subroutine filit(w,nw,nwords,idpij,mxord,ierr)
C--
C--        w           (in)   store, declared w(*) in filit
C--        nw          (in)   number of words available in the store
C--        nwords      (out)  number of words used < 0 not enough space
C--        idpij(7,3)  (out)  Pij(i,iord) table ids, i = QQ,QG,GQ,GG,N+,N-,NV
C--        mxord       (out)  maximum perturbative order LO,NLO,NNLO
C--        ierr        (out)  >0 error condition
C--                           -1 weight tables already exist, do nothing
       

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'qmaps8.inc'
      
      dimension idpij(7,3)
      
      external filit

C--   Nonzero mxord7 indicates that the job is already done      
      if(mxord7(iset).ne.0) then
        ierr   = -1
        nwlast = nwlast7
        return
      endif  
      
C--   Initialize table indices
      do j = 1,3
        do i = 1,7
            idpij(i,j) = 0
        enddo
      enddo      

C--   Loop over lin, quad interpolation 
      iorem = ioy2     
      do ioy2 = 2,iorem
C--     First word of batch of tables      
        ifst7(iset,ioy2) = ipoint7
C--     Words left in the store        
        nwsize = nwf0-ipoint7+1
C--     Fill tables (filit is a generic s/r name passed as an argument)
C--        stor7(ipoint7)  (in)   store
C--        nwsize      (in)   number of words currently available
C--        nwords      (out)  number of words used < 0 not enough space
C--        idpij       (out)  list of Pij table identifiers
C--        mxord       (out)  maximum perturbative order LO,NLO,NNLO
C--        ierr        (out)  0=OK, otherwise not enough space
        call filit(stor7(ipoint7),nwsize,nwords,idpij,mxord,ierr)
C--     Last word occupied in the store        
        nwlast = ipoint7 + abs(nwords)-1
        ilst7(iset,ioy2) = nwlast
C--     Not enough space
        if(nwords.lt.0) goto 10         
C--     First word of next batch of tables
        ipoint7 = nwlast+1        
C--     End of loop over ioy         
      enddo 
      
  10  continue
C--   Restore interpolation order
      ioy2 = iorem        
C--   First word of pdf tables
      ipoint7 = nwlast+1
C--   nwlast = last word used < 0 no space
      call sqcPdfIni(ipoint7,nwlast,idmin,idmax,lpdf)
      
C--   Store table indices        
      do j = 1,3
        do i = 1,7
            idpij7(i,j,iset) = idpij(i,j)
        enddo
      enddo      
C--   Store pdf table offset
      knul7(iset) = kpdf7(0)
C--   Store perturbative order
      mxord7(iset) = mxord
C--   Store last word used
      nwlast7 = nwlast
C--   Pointer to first word of next batch 
      ipoint7 = nwlast+1
C--   Check mxord in range 1-3     
      if(mxord.lt.1 .or. mxord.gt.3) ierr = -2

      return
      end
      
C     =================================================
      subroutine sqcFilWU(w,nw,nwords,idpij,mxord,ierr)
C     =================================================

C--   Fill Pij tables unpolarised

C--   w           (in)   store
C--   nw          (in)   number of words available
C--   nwords      (out)  number of words used < 0 not enough space
C--   idpij       (out)  list of Pij table identifiers
C--   mxord       (out)  maximum perturbative order LO,NLO,NNLO
C--   ierr        (out)  0=OK, otherwise not enough space

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      
      dimension w(*)
      dimension idpij(7,3)
      dimension itypes(4)
      data itypes /0,0,0,0/
      
      external dqcAchi
      external dqcPQQ0R, dqcPQQ0S, dqcPQQ0D            ! (1,1) = PQQ0
      external dqcPQG0A                                ! (2,1) = PQG0
      external dqcPGQ0A                                ! (3,1) = PGQ0
      external dqcPGG0A, dqcPGG0R, dqcPGG0S, dqcPGG0D  ! (4,1) = PGG0
      
      external dqcPPL1A, dqcPPL1B                      ! (5,2) = PPL1
      external dqcPMI1B                                ! (6,2) = PMI1
      external dqcPQQ1A, dqcPQQ1B                      ! (1,2) = PQQ1
      external dqcPQG1A                                ! (2,2) = PQG1
      external dqcPGQ1A                                ! (3,2) = PGQ1
      external dqcPGG1A, dqcPGG1B                      ! (4,2) = PGG1 
      
      external dqcPPL2A, dqcPPL2B, dqcPPL2D            ! (5,3) = PPL2
      external dqcPMI2A, dqcPMI2B, dqcPMI2D            ! (6,3) = PMI2
      external dqcPVA2A                                ! (7,3) = PVA2
      external dqcPQQ2A                                ! (1,3) = PQQ2
      external dqcPQG2A                                ! (2,3) = PQG2
      external dqcPGQ2A                                ! (3,3) = PGQ2
      external dqcPGG2A, dqcPGG2B, dqcPGG2D            ! (4,3) = PGG2
      
      external dqcAGQ2A                                !  AGQ2
      external dqcAGG2A, dqcAGG2B, dqcAGG2D            !  AGG2
      external dqcAQQ2A, dqcAQQ2B, dqcAQQ2D            !  AQQ2
      external dqcAHQ2A                                !  AHQ2
      external dqcAHG2A, dqcAHG2D                      !  AHG2

C--   Max perturbative order
      mxord = 3 
           
C--   Partition
      itypes(1) = 5
      itypes(2) = 17 
      call sqcBookTab(w,nw,itypes,nwords,ierr) 
      if(ierr.ne.0) return
      
C--   LO
      write(lunerr1,'('' Pij LO    for ospline = '',I1)') ioy2
      idPij(1,1)     = 201 !PQQ LO
      call sqcUwgtRS(w,201,ioy2,dqcPQQ0R,dqcPQQ0S,dqcAchi,1,ie)   !1=delta
      call sqcUweitD(w,201,ioy2,dqcPQQ0D,dqcAchi,ie)
      idPij(2,1)     = 202 !PQG LO
      call sqcUweitA(w,202,ioy2,dqcPQG0A,dqcAchi,ie)
      idPij(3,1)     = 203 !PGQ LO
      call sqcUweitA(w,203,ioy2,dqcPGQ0A,dqcAchi,ie)
      idPij(4,1)     = 204 !PGG LO
      call sqcUweitA(w,204,ioy2,dqcPGG0A,dqcAchi,ie)
      call sqcUwgtRS(w,204,ioy2,dqcPGG0R,dqcPGG0S,dqcAchi,1,ie)   !1=delta
      call sqcUweitD(w,204,ioy2,dqcPGG0D,dqcAchi,ie)
      idPij(5,1)     = 201 !PPL LO
      idPij(6,1)     = 201 !PMI LO
      idPij(7,1)     = 201 !PVA LO
C--   NLO
      write(lunerr1,'('' Pij NLO   for ospline = '',I1)') ioy2
      idPij(5,2)     = 205 !PPL NLO
      call sqcUweitA(w,205,ioy2,dqcPPL1A,dqcAchi,ie)
      call sqcUweitB(w,205,ioy2,dqcPPL1B,dqcAchi,1,ie)            !1=delta
      idPij(6,2)     = 206 !PMI NLO
      idPij(7,2)     = 206 !PVA NLO
      call sqcUweitB(w,206,ioy2,dqcPMI1B,dqcAchi,1,ie)            !1=delta
      idPij(1,2)     = 207 !PQQ NLO
      call sqcUweitA(w,207,ioy2,dqcPQQ1A,dqcAchi,ie)
      call sqcUweitB(w,207,ioy2,dqcPQQ1B,dqcAchi,1,ie)            !1=delta
      idPij(2,2)     = 208 !PQG NLO
      call sqcUweitA(w,208,ioy2,dqcPQG1A,dqcAchi,ie)
      idPij(3,2)     = 209 !PGQ NLO
      call sqcUweitA(w,209,ioy2,dqcPGQ1A,dqcAchi,ie)
      idPij(4,2)     = 210 !PGG NLO
      call sqcUweitA(w,210,ioy2,dqcPGG1A,dqcAchi,ie)
      call sqcUweitB(w,210,ioy2,dqcPGG1B,dqcAchi,1,ie)            !1=delta
C--   NNLO
      write(lunerr1,'('' Pij NNLO  for ospline = '',I1)') ioy2
      idPij(5,3)     = 211 !PPL NNLO
      call sqcUweitA(w,211,ioy2,dqcPPL2A,dqcAchi,ie)
      call sqcUweitB(w,211,ioy2,dqcPPL2B,dqcAchi,0,ie)            !0=nodelta
      call sqcUweitD(w,211,ioy2,dqcPPL2D,dqcAchi,ie)
      idPij(6,3)     = 212 !PMI NNLO
      call sqcUweitA(w,212,ioy2,dqcPMI2A,dqcAchi,ie)
      call sqcUweitB(w,212,ioy2,dqcPMI2B,dqcAchi,0,ie)            !0=nodelta
      call sqcUweitD(w,212,ioy2,dqcPMI2D,dqcAchi,ie)
      idPij(7,3)     = 213 !PVA NNLO
      call sqcCopyWt(w,212,w,213,0)
      call sqcUweitA(w,213,ioy2,dqcPVA2A,dqcAchi,ie)
      idPij(1,3)     = 214 !PQQ NNLO
      call sqcCopyWt(w,211,w,214,0)
      call sqcUweitA(w,214,ioy2,dqcPQQ2A,dqcAchi,ie)
      idPij(2,3)     = 215 !PQG NNLO
      call sqcUweitA(w,215,ioy2,dqcPQG2A,dqcAchi,ie)
      idPij(3,3)     = 216 !PGQ NNLO            
      call sqcUweitA(w,216,ioy2,dqcPGQ2A,dqcAchi,ie)
      idPij(4,3)     = 217 !PGG NNLO
      call sqcUweitA(w,217,ioy2,dqcPGG2A,dqcAchi,ie)
      call sqcUweitB(w,217,ioy2,dqcPGG2B,dqcAchi,0,ie)            !0=nodelta
      call sqcUweitD(w,217,ioy2,dqcPGG2D,dqcAchi,ie)
C--   NNLO      
      write(lunerr1,'('' Aij NNLO  for ospline = '',I1)') ioy2
      idAij7(1)      = 101 !AGQ NNLO
      call sqcUweitA(w,101,ioy2,dqcAGQ2A,dqcAchi,ie)
      idAij7(2)      = 102 !AGG NNLO
      call sqcUweitA(w,102,ioy2,dqcAGG2A,dqcAchi,ie)
      call sqcUweitB(w,102,ioy2,dqcAGG2B,dqcAchi,0,ie)            !0=nodelta
      call sqcUweitD(w,102,ioy2,dqcAGG2D,dqcAchi,ie)
      idAij7(3)      = 103 !AQQ NNLO
      call sqcUweitA(w,103,ioy2,dqcAQQ2A,dqcAchi,ie)
      call sqcUweitB(w,103,ioy2,dqcAQQ2B,dqcAchi,0,ie)            !0=nodelta
      call sqcUweitD(w,103,ioy2,dqcAQQ2D,dqcAchi,ie)
      idAij7(4)      = 104 !AHQ NNLO
      call sqcUweitA(w,104,ioy2,dqcAHQ2A,dqcAchi,ie)
      idAij7(5)      = 105 !AHG NNLO
      call sqcUweitA(w,105,ioy2,dqcAHG2A,dqcAchi,ie)
      call sqcUweitD(w,105,ioy2,dqcAHG2D,dqcAchi,ie)

      return
      end
      
C     =================================================
      subroutine sqcFilWP(w,nw,nwords,idpij,mxord,ierr)
C     =================================================

C--   Fill Pij tables polarised

C--   w           (in)   store
C--   nw          (in)   number of words available
C--   nwords      (out)  number of words used < 0 not enough space
C--   idpij       (out)  list of Pij table identifiers
C--   mxord       (out)  maximum perturbative order LO,NLO,NNLO
C--   ierr        (out)  0=OK, otherwise not enough space

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      
      dimension w(*)
      dimension idpij(7,3)
      dimension itypes(4)
      data itypes /0,0,0,0/
      
      external dqcAchi
      external dqcDPQQ0A, dqcDPQQ0B, dqcDPQQ0D            ! (1,1) = DPQQ0
      external dqcDPQG0A                                  ! (2,1) = DPQG0
      external dqcDPGQ0A                                  ! (3,1) = DPGQ0
      external dqcDPGG0A, dqcDPGG0B, dqcDPGG0D            ! (4,1) = DPGG0
      
      external dqcDPPL1B                                  ! (5,2) = DPPL1
      external dqcDPMI1A, dqcDPMI1B                       ! (6,2) = DPMI1
      external dqcDPQS1A                                  ! (1,2) = DPQS1
      external dqcDPQG1A                                  ! (2,2) = DPQG1
      external dqcDPGQ1A                                  ! (3,2) = DPGQ1
      external dqcDPGG1A, dqcDPGG1R, dqcDPGG1S, dqcDPGG1D ! (4,2) = DPGG1
      
C--   Max perturbative order
      mxord = 2 

C--   Partition      
      itypes(1) = 2
      itypes(2) = 8
      call sqcBookTab(w,nw,itypes,nwords,ierr) 
      if(ierr.ne.0) return

C--   LO
      write(lunerr1,'('' Pij LO    for ospline = '',I1)') ioy2
      idPij(1,1)     = 101 !DPQQ LO
      call sqcUweitA(w,101,ioy2,dqcDPQQ0A,dqcAchi,ie)
      call sqcUweitB(w,101,ioy2,dqcDPQQ0B,dqcAchi,1,ie)           !1=delta
      call sqcUweitD(w,101,ioy2,dqcDPQQ0D,dqcAchi,ie)
      idPij(2,1)     = 201 !DPQG LO
      call sqcUweitA(w,201,ioy2,dqcDPQG0A,dqcAchi,ie)
      idPij(3,1)     = 102 !DPGQ LO
      call sqcUweitA(w,102,ioy2,dqcDPGQ0A,dqcAchi,ie)
      idPij(4,1)     = 202 !DPGG LO
      call sqcUweitA(w,202,ioy2,dqcDPGG0A,dqcAchi,ie)
      call sqcUweitB(w,202,ioy2,dqcDPGG0B,dqcAchi,1,ie)           !1=delta
      call sqcUweitD(w,202,ioy2,dqcDPGG0D,dqcAchi,ie)
      idPij(5,1)     = 101 !DPPL LO
      idPij(6,1)     = 101 !DPMI LO
      idPij(7,1)     = 101 !DPVA LO
      
C--   NLO
      write(lunerr1,'('' Pij NLO   for ospline = '',I1)') ioy2
      idPij(5,2)     = 203 !DPPL NLO
      call sqcUweitB(w,203,ioy2,dqcDPPL1B,dqcAchi,1,ie)           !1=delta
      idPij(6,2)     = 204 !DPMI NLO
      call sqcUweitA(w,204,ioy2,dqcDPMI1A,dqcAchi,ie)
      call sqcUweitB(w,204,ioy2,dqcDPMI1B,dqcAchi,1,ie)           !1=delta
      idPij(7,2)     = 204 !DPVA NLO
      idPij(1,2)     = 205 !PQQ NLO
      call sqcCopyWt(w,203,w,205,0)
      call sqcUweitA(w,205,ioy2,dqcDPQS1A,dqcAchi,ie)
      idPij(2,2)     = 206 !PQG NLO
      call sqcUweitA(w,206,ioy2,dqcDPQG1A,dqcAchi,ie)
      idPij(3,2)     = 207 !PGQ NLO
      call sqcUweitA(w,207,ioy2,dqcDPGQ1A,dqcAchi,ie)
      idPij(4,2)     = 208 !PGG NLO
      call sqcUweitA(w,208,ioy2,dqcDPGG1A,dqcAchi,ie)
      call sqcUwgtRS(w,208,ioy2,dqcDPGG1R,dqcDPGG1S,dqcAchi,1,ie) !1=delta
      call sqcUweitD(w,208,ioy2,dqcDPGG1D,dqcAchi,ie)

      return
      end
      
C     =================================================
      subroutine sqcFilWF(w,nw,nwords,idpij,mxord,ierr)
C     =================================================

C--   Fill Pij tables timelike  (fragmentation functions)

C--   w           (in)   store
C--   nw          (in)   number of words available
C--   nwords      (out)  number of words used < 0 not enough space
C--   idpij       (out)  list of Pij table identifiers
C--   mxord       (out)  maximum perturbative order LO,NLO,NNLO
C--   ierr        (out)  0=OK, otherwise not enough space

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      
      dimension w(*)
      dimension idpij(7,3)
      dimension itypes(4)
      data itypes /0,0,0,0/
      
      external dqcAchi
      external dqcPQQ0R, dqcPQQ0S, dqcPQQ0D            ! (1,1) = PQQ0
      external dqcTQG0A                                ! (2,1) = PQG0
      external dqcTGQ0A                                ! (3,1) = PGQ0
      external dqcPGG0A, dqcPGG0R, dqcPGG0S, dqcPGG0D  ! (4,1) = PGG0
      
      external dqcTPL1A, dqcTPL1B                      ! (5,2) = PPL1
      external dqcTMI1B                                ! (6,2) = PMI1
      external dqcTQQ1A, dqcTQQ1B                      ! (1,2) = PQQ1
      external dqcTQG1A                                ! (2,2) = PQG1
      external dqcTGQ1A                                ! (3,2) = PGQ1
      external dqcTGG1A, dqcTGG1B                      ! (4,2) = PGG1 
      
C--   Max perturbative order
      mxord = 2 
           
C--   Partition
      itypes(2) = 10 
      call sqcBookTab(w,nw,itypes,nwords,ierr) 
      if(ierr.ne.0) return
      
C--   LO (QQ and GG same as spacelike unpolarised)
      write(lunerr1,'('' Pij LO    for ospline = '',I1)') ioy2
      idPij(1,1)     = 201 !PQQ LO
      call sqcUwgtRS(w,201,ioy2,dqcPQQ0R,dqcPQQ0S,dqcAchi,1,ie)   !1=delta
      call sqcUweitD(w,201,ioy2,dqcPQQ0D,dqcAchi,ie)
      idPij(2,1)     = 202 !PQG LO
      call sqcUweitA(w,202,ioy2,dqcTQG0A,dqcAchi,ie)
      idPij(3,1)     = 203 !PGQ LO
      call sqcUweitA(w,203,ioy2,dqcTGQ0A,dqcAchi,ie)
      idPij(4,1)     = 204 !PGG LO
      call sqcUweitA(w,204,ioy2,dqcPGG0A,dqcAchi,ie)
      call sqcUwgtRS(w,204,ioy2,dqcPGG0R,dqcPGG0S,dqcAchi,1,ie)   !1=delta
      call sqcUweitD(w,204,ioy2,dqcPGG0D,dqcAchi,ie)
      idPij(5,1)     = 201 !PPL LO
      idPij(6,1)     = 201 !PMI LO
      idPij(7,1)     = 201 !PVA LO
C--   NLO
      write(lunerr1,'('' Pij NLO   for ospline = '',I1)') ioy2
      idPij(5,2)     = 205 !PPL NLO
      call sqcUweitA(w,205,ioy2,dqcTPL1A,dqcAchi,ie)
      call sqcUweitB(w,205,ioy2,dqcTPL1B,dqcAchi,1,ie)            !1=delta
      idPij(6,2)     = 206 !PMI NLO
      idPij(7,2)     = 206 !PVA NLO
      call sqcUweitB(w,206,ioy2,dqcTMI1B,dqcAchi,1,ie)            !1=delta
      idPij(1,2)     = 207 !PQQ NLO
      call sqcUweitA(w,207,ioy2,dqcTQQ1A,dqcAchi,ie)
      call sqcUweitB(w,207,ioy2,dqcTQQ1B,dqcAchi,1,ie)            !1=delta
      idPij(2,2)     = 208 !PQG NLO
      call sqcUweitA(w,208,ioy2,dqcTQG1A,dqcAchi,ie)
      idPij(3,2)     = 209 !PGQ NLO
      call sqcUweitA(w,209,ioy2,dqcTGQ1A,dqcAchi,ie)
      idPij(4,2)     = 210 !PGG NLO
      call sqcUweitA(w,210,ioy2,dqcTGG1A,dqcAchi,ie)
      call sqcUweitB(w,210,ioy2,dqcTGG1B,dqcAchi,1,ie)            !1=delta      

      return
      end
      
C     ===================================================
      subroutine sqcReadPij(w,nw,nwords,idpij,mxord,ierr)
C     ===================================================

C--   Fill Pij tables by reading them from disk

C--   w           (in)   store
C--   nw          (in)   number of words available
C--   nwords      (out)  number of words used < 0 not enough space
C--   idpij       (out)  list of Pij table identifiers
C--   mxord       (out)  maximum perturbative order LO,NLO,NNLO
C--   ierr        (out)  0 = all OK
C--                      1 = read error
C--                      5 = not enough space
C--
C--   The logical unit number is passed via /qluns1/lunwgt1


      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      
      dimension w(*)
      dimension idpij(7,3)

      ierr = 0
      
      read(lunwgt1,err=500,end=500) nwords
      if(nwords.le.0) goto 500
      read(lunwgt1,err=500,end=500) idpij, mxord

C--   Not enough space
      if(nwords.gt.nw) then
        nwords = -nwords 
        ierr   = 5
        return
      endif
      
C--   Now read the store      
      read(lunwgt1,err=500,end=500) (w(i),i=1,nwords)

      return

  500 continue
C--   Read error
      ierr = 1

      return
      end
      
C     ===============================
      subroutine sqcRdumPij(lun,ierr)
C     ===============================

C--   Dummy read of pij tables, used to position the file

C--   ierr        (out)  0 = all OK
C--                      1 = read error

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qluns1.inc'
      
      dimension idpij(7,3)

      ierr = 0
      
      read(lun,err=500,end=500) nwords
      if(nwords.le.0) goto 500
      read(lun,err=500,end=500) idpij, mxord      
      read(lun,err=500,end=500) (dummy,i=1,nwords)

      return

  500 continue
C--   Read error
      ierr = 1

      return
      end
      
C     =======================================      
      subroutine sqcDumpWt(lun,jset,key,ierr)
C     =======================================

C--   Write all sets of Pij tables to disk 
C--
C--   lun    (in)  logical unit number
C--   jset   (in)  = 0 write all exising tables except custom tables
C--                # 0 write pdf set jset 
C--   key    (in)  key character string
C--   ierr   (out) 0 = OK
C--                1 = write error
C--                2 = no sets to be written

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qvers1.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      
      dimension     idpij(7,3)
      character*(*) key
      character*50  keyout
      
C--   Initialize
      call sqcSetKey(key,keyout)      
C--   Dump QCDNUM version
      write(lun,err=500) cvers1, cdate1
C--   Dump key
      write(lun,err=500) keyout      
C--   Some array sizes
      write(lun,err=500) mxg0, mxx0, mqq0, miw0
C--   Relevant grid parameters
      write(lun,err=500) nyy2, nyg2, ioy2, dely2
      write(lun,err=500) ntt2
      write(lun,err=500) (tgrid2(i),i=1,ntt2)
C--   Now flush out the store
      ierr = 2
C--   Loop over all requested sets
      if(jset.eq.0)     then
        iset1 = 1
        iset2 = 3
      else
        iset1 = jset
        iset2 = jset
      endif                  
      do iset = iset1,iset2
        if(mxord7(iset).ge.1 .and. mxord7(iset).le.3) then
C--       Set is avaialble, copy idPij
          ierr = 0
          do j = 1,3
            do i = 1,7
              idpij(i,j) = idPij7(i,j,iset)
            enddo
          enddo
          mxord = mxord7(iset)
C--       Write iset
          write(lun,err=500) iset
C--       Write idAij
          write(lun,err=500) idAij7          
C--       Loop over interpolation orders
          do ioy = 2,ioy2
            nw1    = ifst7(iset,ioy)
            nwords = ilst7(iset,ioy)-ifst7(iset,ioy)+1
            call sqcDumpW(lun,stor7(nw1),nwords,idpij,mxord,jerr)
            if(jerr.ne.0) goto 500
          enddo          
        endif
      enddo                   
C--   End of loop over sets, write an 'end-of-file'              
      iset = 0
      write(lun,err=500) iset
      
      return
      
  500 continue
C--   Write error
      ierr = 1

      return
      end  

C     ==================================================
      subroutine sqcDumpW(lun,w,nwords,idpij,mxord,ierr)
C     ==================================================

C--   Write one set of Pij tables to disk

C--   lun         (in)   logical unit number
C--   w           (in)   store
C--   nwords      (in)   number of words to be written 
C--   idpij       (in)   list of Pij table identifiers
C--   mxord       (in)   maximum perturbative order LO,NLO,NNLO
C--   ierr        (out)  0 = OK
C--                      1 = write error

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qstor7.inc'
      
      dimension w(*)
      dimension idpij(7,3)

      ierr = 0
      
      write(lun,err=500) nwords
      write(lun,err=500) idpij, mxord
      write(lun,err=500) (w(i),i=1,nwords)

      return

  500 continue
C--   Write error
      ierr = 1

      return
      end                 
      
C     ===============================================
      subroutine sqcReadWt(lun,key,nwlast,iread,ierr)
C     ===============================================

C--   Read weights from a disk file (unformatted read) and partition
C--   the store (that is, partiton weight tables and also pdf tables)
C--  
C--   lun          (in)   input logical unit number
C--   key          (in)   key character string
C--   nwlast       (out)  last word used in the store < 0 no space
C--   iread(mset0) (out)  list of pdf sets read in -1 set already exists
C--                                                 0 not on disk
C--                                                 1 set read in from disk
C--   ierr         (out)  0 = all OK
C--                       1 = read error
C--                       2 = problem with QCDNUM version
C--                       3 = key mismatch
C--                       4 = x-mu2 grid not the same
C--                       5 = not enough space

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qvers1.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qmaps8.inc'
      
      character*10  cversr
      character*8   cdater
      character*50  keyred
      character*(*) key
      logical lqcSjekey
      dimension nyyr(0:mxg0),delyr(0:mxg0)
      dimension tgridr(mqq0),iread(mset0)
      
      external sqcReadPij
      
C--   Pass lun via common block
      lunwgt1 = lun
C--   Initialize      
      nwlast  = 0
      ierr    = 0
      do i = 1,mset0
        iread(i) = 0
      enddo  
      
C--   QCDNUM version
      read(lun,err=500,end=500) cversr, cdater
      if(cversr.ne.cvers1 .or. cdater.ne.cdate1) then
        ierr = 2
        return
      endif
C--   Key
      read(lun,err=500,end=500) keyred
      if(.not.lqcSjekey(key,keyred)) then
        ierr = 3
        return
      endif    
C--   Some array sizes
      read(lun,err=500,end=500) mxgr, mxxr, mqqr, miwr
      if(mxgr.ne.mxg0 .or. mxxr.ne.mxx0 .or. 
     +   mqqr.ne.mqq0 .or. miwr.ne.miw0) then
        ierr = 2
        return
      endif
C--   Relevant x-grid parameters
      read(lun,err=500,end=500) nyyr, nygr, ioyr, delyr
      if(nygr.ne.nyg2 .or. ioyr.ne.ioy2) then
        ierr = 4
        return
      endif
      do i = 0,mxg0
        if(nyyr(i).ne.nyy2(i) .or. delyr(i).ne.dely2(i)) then
          ierr = 4
          return
        endif
      enddo
C--   Mu2 grid
      read(lun,err=500,end=500) nttr
      if(nttr .ne. ntt2) then
        ierr = 4
        return
      endif
      read(lun,err=500,end=500) (tgridr(i),i=1,ntt2)
      do i = 1,ntt2
        if(tgrid2(i).ne.tgridr(i)) then
          ierr = 4
          return
        endif
      enddo
      
C--   Now read Pij sets
  10  continue  
      read(lun,err=500) iset
C--   No more sets to read      
      if(iset.le.0 .or. iset.ge. 5) return
      iread(iset) = 1
C--   Read idAij      
      read(lun,err=500) idAij7
C--   Read Pij      
      call sqcFilWt(sqcReadPij,iset,nwlast,ierr)
      if(ierr.eq.-1) then
C--     No-op because Pij tables did already exist --> do dummy read
        do ioy = 2,ioy2
          call sqcRdumPij(lun,jerr)
          if(jerr.ne.0) goto 500
        enddo  
        iread(iset) = -1
        ierr = 0
      endif
C--   Not enough space
      if(nwlast.le.0) then
        ierr = 5
        return
      endif  
      goto 10

C--   NFmap, splines etc
      if(.not.Lnfmap8) call sqcNFtab(0)
      
      return
      
C--   Read error      
 500  continue     
      ierr = 1
      
      return
      end
      