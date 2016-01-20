
C--   This is the file usrevol.f containing the evolution routines

C--   double precision function rfromf(f2)
C--   double precision function ffromr(r2)
C--   function asfunc(r2,nfout,ierr)
C--   double precision function EvolAs(iord,as0,r20,r2,iqcdnum,nfout,ierr)
C--   subroutine EvolFG(itype,func,def,iq0,epsi)


C==   ===============================================================
C==   Scale transformations =========================================
C==   ===============================================================
      
C     ====================================
      double precision function rfromf(f2)
C     ====================================

C--   Convert factorization scale f2 to renormalization scale r2.

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      logical first
      save    first
      data    first /.true./

      character*80 subnam
      data subnam /'RFROMF ( F2 )'/

C--   Check if QCDNUM is initialized
      if(first) then
        call sqcChkIni(subnam)
        first = .false.
      endif

C--   Transform
      rfromf = aar6*f2 + bbr6

      return
      end
      
C     ====================================
      double precision function ffromr(r2)
C     ====================================

C--   Convert renormalization scale r2 to factorization scale f2.

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      logical first
      save    first
      data    first /.true./

      character*80 subnam
      data subnam /'FFROMR ( R2 )'/

C--   Check if QCDNUM is initialized
      if(first) then
        call sqcChkIni(subnam)
        first = .false.
      endif

C--   Transform
      ffromr = (r2-bbr6)/aar6

      return
      end
      
C==   ===============================================================
C==   Alpha-s =======================================================
C==   ===============================================================
            
C     ===============================================
      double precision function asfunc(r2,nfout,ierr)
C     ===============================================

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'

      dimension rmsq(4:6)
 
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      character*80 subnam
      data subnam /'ASFUNC ( R2, NFOUT, IERR )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

      jord = iord6
      alf0 = alfq06
      ralf = q0alf6
      do i = 4,6
        rmsq(i) = aar6*qthrs6(i) + bbr6
      enddo

      asfunc = 
     & dqcAsEvol(r2,ralf,alf0,rmsq,jord,nfout,ierr)

      return
      end
            
C     ===========================================================
      double precision function 
     +                 EvolAs(iord,as0,r20,r2,iqcdnum,nfout,ierr)
C     ===========================================================

C--   Evolve alphas(r2) at order iord.
C--   This function is not documented in the write-up but kept for
C--   backward compatibility
C--
C--   Input   iord     1 = LO, 2 = NLO, 3 = NNLO, other = QCDNUM setting
C--           as0      Starting alphas, QCDNUM setting if .le. 0  or .gt. 1
C--           r20      Starting R2,     QCDNUM setting if .le. 0.1
C--           r2       End scale should be .gt. 0.1
C--           iqcdnum  0 take thesholds as set by setthr (or by default)
C--                    1 take internal thresholds (rounded down to grid points)
C--   
C--   Output  nfout    number of flavors at r2
C--           ierr     0 all OK
C--                    1 r2 at or below lambda2
C--                    2 no internal thresholds available

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'

      dimension rmsq(4:6)
 
      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      character*80 subnam
      data subnam 
     +    /'EVOLAS ( IORD, AS0, R20, R2, INTERN, NFOUT, IERR )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Order taken from user input or from current QCDNUM setting
      if(iord.ge.1 .and. iord.le.3) then
        jord = iord
      else
        jord = iord6
      endif
C--   Alpha0 taken from user or from current QCDNUM setting
      if(as0.gt.0.D0 .and. as0.lt.1.D0) then
        alf0 = as0
      else
        alf0 = alfq06
      endif
C--   R20 taken from user or from current QCDNUM setting
      if(abs(r20).gt.0.1) then
        ralf = r20
      else
        ralf = q0alf6
      endif 
C--   Check which thresholds to take
      if(iqcdnum.ne.1) then
C--     Use the thresholds defined by setthrs or by default
C--     Convert these from factorization to renormalization scale
        do i = 4,6
          rmsq(i) = aar6*qthrs6(i) + bbr6
        enddo
      else
C--     Use the internal QCDNUM values
C--     First, check that grid is available (no error msg)
        call sqcChekit(1,ichk,jbit)
C--     error
        if(jbit.ne.0) then
          evolas = qnull6
          ierr   = 2
          return
        endif
C--     OK, now get the pole mass squared from QCDNUM
        do i = 4,6
          rmsq(i) = aar6*qthrs6(i) + bbr6
        enddo
      endif

      evolas = 
     & dqcAsEvol(r2,ralf,alf0,rmsq,jord,nfout,ierr)
      
      return
      end
          
C==   ===============================================================
C==   PDF evolution =================================================
C==   ===============================================================
            
      
C     ==========================================
      subroutine EvolFG(itype,func,def,iq0,epsi)
C     ==========================================

C--   Evolve all pdfs.
C-- 
C--   In:   itype       :  Type 1=unpol, 2=pol, 3=timelike, 4=custom 
C--         func(ipf,x) :  Function that returns a set of input pdfs
C--         def(i,ipdf) :  Flavor decomposition of each pdf returned by func()
C--         iq0         :  Starting scale
C--   Out:  epsi        :  Largest deviation quad - lin at midpoint
C--   
C--   Flavor indices :  tb bb cb sb ub db  g  d  u  s  c  b  t
C--                     -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
C--   q+- indices    :   g  d+ u+ s+ c+ b+ t+ d- u- s- c- b- t-
C--                      0  1  2  3  4  5  6  7  8  9 10 11 12  
C--   Si/ns  indices :   g  si 2+ 3+ 4+ 5+ 6+ va 2- 3- 4- 5- 6-
C--                      0  1  2  3  4  5  6  7  8  9 10 11 12    

      implicit double precision (a-h,o-z)

      external func

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'qmaps8.inc'

      logical first
      save    first
      data    first /.true./

      save      ichk,       iset,       idel
      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

      dimension def(-6:6,12), pdef(12,12)
      
*mb
*mb      call sqcDebug('EVOLFG')
*mb 

      character*80 subnam
      data subnam /'EVOLFG ( ITYPE, FUNC, DEF, IQ0, EPSI )'/

C--   Initialize flagbits
      if(first) then
        call sqcMakeFl(subnam,ichk,iset,idel)
        first = .false.
      endif

C--   Check evolution type
      call sqcIlele(subnam,'ITYPE',1,itype,4,' ')      

C--   Check status bits
      call sqcChkflg(itype,ichk,subnam)

C--   Setup flavor map
      if(.not.Lnfmap8) call sqcNfTab(0)
      
C--   Check perturbative order
      call sqcIlele(subnam,'IORD',1,iord6,mxord7(itype),
     &              'Too large order, please call SETORD')       

C--   Check user input
      if(nfix6.eq.0) then
C--     VFNS
        call sqcIlele(subnam,'IQ0',1,iq0,itchm2-1,
     &                'IQ0 should be below charm threshold')
      else
C--     FFNS
        call sqcIlele(subnam,'IQ0',1,iq0,ntt2,' ')
      endif

C--   Initialize pdef 
      do i = 1,12
        do j = 1,12
          pdef(i,j) = 0.D0
        enddo
      enddo
C--   Predefine heavy quark qplus  in case nfmax < 4
      pdef( 7, 4) = 1.D0      !cplus
      pdef( 8,10) = 1.D0      !cmin
      pdef( 9, 5) = 1.D0      !bplus
      pdef(10,11) = 1.D0      !bmin
      pdef(11, 6) = 1.D0      !tplus
      pdef(12,12) = 1.D0      !tmin
C--   Build q+- matrix (pdef) from user input (def)
C--   Note index swap def(iflavor,ipdf) --> pdef(ipdf,iflavor)
      nfmax = max(nfix6,3)
      do ipdf = 1,2*nfmax
        do j = 1,6
          pdef(ipdf,j)   = 0.5*(def(j,ipdf)+def(-j,ipdf))
          pdef(ipdf,j+6) = 0.5*(def(j,ipdf)-def(-j,ipdf))
        enddo
      enddo 

C--   Now tell qcdnum about these definitions
      call sqcPdIdef(pdef,ierr)
C--   Well that did not go OK...
      if(ierr.ne.0) call sqcErrMsg(subnam,
     +             'Input quark densities not linearly independent')

C--   Enter input distributions
      call sqcAllInp(itype,func)
C--   Cut range           
      iqminc = itfiz2(izmic2)
      iqmaxc = itfiz2(izmac2)
C--   Extend cut range to include iq0      
      iq1    = min(iq0,iqminc)
      iq2    = max(iq0,iqmaxc)
C--   Off we go...
*      call sqcEvolve(itype,iord6,iq0,1,ntt2,epsi,ierr)
      call sqcEvolve(itype,iord6,iq0,iq1,iq2,epsi,ierr)
      if(ierr.eq.1) call sqcErrMsg(subnam,
     +                   'Attempt to evolve with too large alpha-s')
C--   Check max deviation
      if(dflim6.gt.0.D0 .and. epsi.gt.dflim6) call sqcErrMsg(subnam, 
     +          'Possible spline oscillation detected')
C--   Update status bits
      call sqcSetflg(iset,idel,itype)
      
      return
      end
      

      
