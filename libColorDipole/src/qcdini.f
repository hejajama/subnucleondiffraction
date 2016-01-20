
C--   This is the file qcdini.f containing the qcdnum initialization routines

C--   subroutine sqcIniCns
C--   subroutine sqcSetLun(lun,fname)
C--   subroutine sqcBanner(lun)

C===================================================================
C==   Initialization routines ======================================
C===================================================================

C     ====================
      subroutine sqcIniCns
C     ====================

C---  Initialize constants.
C---  Called by sqc_qcinit.
 
      implicit double precision (a-h,o-z)

      include 'qcdnum.inc' 
      include 'qconst.inc'
      include 'qvers1.inc'
      include 'qgrid2.inc'
      include 'qsnam3.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'qmaps8.inc'
      include 'qfast9.inc'

C--   Fixed parameters
C--   ----------------

      pi     = 3.14159265359
      proton = 0.9382796
      eutron = 0.9395731
      ucleon = (proton + eutron) / 2.

C--   These constants are inhereted from Ouarou and Virchaux (original QCDNUM) 
      c1s3   = 1./3.
      c2s3   = 2./3.
      c4s3   = 4./3.
      c5s3   = 5./3.
      c8s3   = 8./3.
      c14s3  = 14./3.
      c16s3  = 16./3.
      c20s3  = 20./3.
      c28s3  = 28./3.
      c38s3  = 38./3.
      c40s3  = 40./3.
      c44s3  = 44./3.
      c52s3  = 52./3.
      c136s3 = 136./3.
      c11s6  = 11./6.
      c2s9   = 2./9.
      c4s9   = 4./9.
      c10s9  = 10./9.
      c14s9  = 14./9.
      c16s9  = 16./9.
      c40s9  = 40./9.
      c44s9  = 44./9.
      c62s9  = 62./9.
      c112s9 = 112./9.
      c182s9 = 182./9.
      c11s12 = 11./12.
      c35s18 = 35./18.
      c11s3  = 11./3.
      c22s3  = 22./3.
      c61s12 = 61./12.
      c215s1 = 215./12.
      c29s12 = 29./12.
      cpi2s3 = pi**2/3.
      cpia   = 67./18. - cpi2s3/2.
      cpib   = 4.*cpi2s3
      cpic   = 17./18. + 3.5*cpi2s3
      cpid   = 367./36. - cpi2s3
      cpie   = 5. - cpi2s3
      cpif   = cpi2s3 - 218./9.

      cca    = 3.
      ccf    = (cca*cca-1.)/(2.*cca)
      ctf    = 0.5
      catf   = cca*ctf
      cftf   = ccf*ctf

C--   Initialize storage
C--   ------------------
      do i = 1,nwf0
        stor7(i) = 0.D0
      enddo
      
C--   No grid available
      Lygrid2 = .false.
      Ltgrid2 = .false.
C--   Weight tables not initialized      
      Lwtini8 = .false.
C--   No flavor map available
      Lnfmap8 = .false.
C--   No cut map available
      Levcut8 = .false.      
C--   No alphas table available
      Lastab8 = .false.
C--   No fast structure functions
      nmax9 = 0
      nnff9 = 0

C--   Blank addon package subroutine name (smb_cfill comes from mbutil)
C--   -----------------------------------------------------------------
      call smb_cfill(' ',usrnam3)

C--   Default parameters which can be changed by the user
C--   ---------------------------------------------------

C--   Precisions and limits
      aepsi6 = 1.D-9          !absolute precision of comparison
      repsi6 = 1.D-9          !relative precision of comparison
      gepsi6 = 1.D-7          !requested accuracy Gauss integration
      dflim6 = 0.5D0          !max abs deviation of spline from function
      aslim6 = 10.D0          !largest allowed value of alphas
      qnull6 = 1.D11          !qcdnum null value
      qlimd6 = 0.1D0          !lowest allowed mu2 value
      qlimu6 = 1.D11          !largest allowed mu2 value
      niter6 = 1              !number of dnward evolution iterations

C--   Default FFNS with nf = 3
      nfix6     = 3
      call sqcThrFFNS(nfix6)

      q0alf6 = 8315.25D0      !MZ^2
      alfq06 = 0.118D0        !alphas(MZ^2)

C--   Default for renormalization scale r2 = aar6*f2 + bbr6
      aar6   = 1.D0
      bbr6   = 0.D0
      itmin6 = 1
      itlow8 = 1

C--   Default is NLO
      iord6  = 2
      
C--   Default number of scratch tables for convolution engine
      nscratch6 = 5      

      return
      end

C     ===============================
      subroutine sqcSetLun(lun,fname)
C     ===============================

C--   Redirect QCDNUM output

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*(*) fname

      lunerr1 = lun
      if(lun.ne.6) then
        open(unit=lun,file=fname,status='unknown')
      endif

      return
      end

C     =========================
      subroutine sqcBanner(lun)
C     =========================

C--   QCDNUM banner printout

      implicit double precision (a-h,o-z)
      include 'qvers1.inc'
	return
      
      write(lun,'('' '')')
      write(lun,'(''                  ///                 '',
     &            ''                 .().                 '')')
      write(lun,'(''                 (..)                 '',
     &            ''                 (--)                 '')')
      write(lun,'(''  +----------ooO--()--Ooo-------------'',
     &            ''-------------ooO------Ooo---------+   '')')
      write(lun,'(''  |                                   '',
     &            ''                                  |   '')')
      write(lun,'(''  |    #####      ######    ######    '',
     &            '' ##    ##   ##    ##   ##     ##  |   '')')
      write(lun,'(''  |   ##   ##    ##    ##   ##   ##   '',
     &            '' ###   ##   ##    ##   ###   ###  |   '')')
      write(lun,'(''  |  ##     ##   ##    ##   ##    ##  '',
     &            '' ####  ##   ##    ##   #### ####  |   '')')
      write(lun,'(''  |  ##     ##   ##         ##    ##  '',
     &            '' ## ## ##   ##    ##   ## ### ##  |   '')')
      write(lun,'(''  |  ##     ##   ##         ##    ##  '',
     &            '' ##  ####   ##    ##   ##  #  ##  |   '')')
      write(lun,'(''  |   ##   ##    ##    ##   ##   ##   '',
     &            '' ##   ###   ##    ##   ##     ##  |   '')')
      write(lun,'(''  |    #####      ######    ######    '',
     &            '' ##    ##    ######    ##     ##  |   '')')
      write(lun,'(''  |        ##                         '',
     &            ''                                  |   '')')
      write(lun,'(''  |                                   '',
     &            ''                                  |   '')')
      write(lun,'(''  |    Version '',A10,''  '',A8,''   '',
     &            ''      Author m.botje@nikhef.nl    |   '')') 
     &            cvers1,cdate1
      write(lun,'(''  |                                   '',
     &            ''                                  |   '')')
      write(lun,'(''  +-----------------------------------'',
     &            ''----------------------------------+   ''//)')
     
      return
      end
      
C     =========================
      subroutine sqcReftoo(lun)
C     =========================

C--   QCDNUM banner printout
       return
     
      write(lun,'(''  +-----------------------------------'',
     &            ''----------------------------------+   '')')
      write(lun,'(''  |                                   '',
     &            ''                                  |   '')')
      write(lun,'(''  |    If you use QCDNUM, please refer'',
     &            '' to:                              |   '')')
      write(lun,'(''  |                                   '',
     &            ''                                  |   '')')
      write(lun,'(''  |    M. Botje, Comput. Phys. Commun.'',
     &            '' 182(2011)490, arXiV:1005.1481    |   '')')
      write(lun,'(''  |                                   '',
     &            ''                                  |   '')')
      write(lun,'(''  +-----------------------------------'',
     &            ''----------------------------------+   '')')
    

      return
      end
