
C--   This is the file qcdevol.f containing the qcdnum evolution routines

C--   subroutine sqcEvolve(itype,iord,it0,it1,it2,epsm,ierr)

C--   subroutine sqcNStart(itype,idout,idin,ids,idg,iyg,iord,dlam,it0)
C--   subroutine sqcGridns(itype,
C--              ipdf,ids,idg,ityp,iyg,iord,it0,it1,it2,eps,ierr)
C--   subroutine sqcNSevnf(itype,ipdf,ityp,iyg,iord,nf,iw1,iw2)
C--   subroutine sqcNSder(itype,ipdf,ityp,iyg,iord,nf,iwa,iwader)
C--   subroutine sqcNSjup(itype,idq,ids,idg,iyg,iord,dlam,ny,iwin,iwout)
C--   subroutine sqcNSjdn(itype,idq,ids,idg,iyg,iord,dlam,ny,iwin,iwout)
C--   subroutine sqcNSStoreStart(itype,ipdf,iy1,iy2,it0)
C--   subroutine sqcNSNewStart(itype,ipdf,iy1,iy2,it0,epsi)
C--   subroutine sqcNSRestoreStart(itype,ipdf,iy1,iy2,it0)
C--
C--   subroutine sqcGridsg(itype,
C--                        idf,idg,iyg,iord,it0,it1,it2,eps,ierr)
C--   subroutine sqcSGevnf(itype,idf,idg,iyg,iord,nf,iw1,iw2)
C--   subroutine sqcSGder(itype,idf,idg,iyg,iord,nf,iwa,iwader)
C--   subroutine sqcSGjup(itype,ids,idg,iyg,iord,ny,iwin,iwout)
C--   subroutine sqcSGjdn(itype,ids,idg,iyg,iord,ny,iwin,iwout)
C--   subroutine sqcSGStoreStart(itype,ids,idg,iy1,iy2,it0)
C--   subroutine sqcSGNewStart(itype,ids,idg,iy1,iy2,it0,epsi)
C--   subroutine sqcSGRestoreStart(itype,ids,idg,iy1,iy2,it0)
C--
C--   double precision function dqcGetEps(itype,id,ny,it)
C--
C--   subroutine sqcNfTab(it0)
C--   subroutine sqcEvCuts(iprint)
C--   subroutine sqcTimeFac(facL,facQ)
C--   integer function iqcIyMaxG(iymax0,ig)
C--   subroutine sqcEvLims(it0,it1,it2,
C--              iwu1,iwu2,nflu,nup,iwd1,iwd2,nfld,ndn,ibl,ibu)


C     ================================================================
C     Steering routine to evolve all pdfs in the FFNS and VFNS
C     ================================================================

C     ======================================================
      subroutine sqcEvolve(itype,iord,it0,it1,it2,epsm,ierr)
C     ======================================================

C--   Steering routine to evolve all pdfs
C--   This routine does not set the starting values
C--   epsm  (out) max deviation quad - lin at midpoints 
C--   ierr  (out) 1 = no alphas available at it1

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'qmaps8.inc'

*mbdebug
      common /debug/lundbg, igdbg, iqdbg
*mbdebug

      dimension dlam(4:6)
      dimension lambda(4:6,12) !heavy quark weights in NNLO jumps
C--          nf =   4  5  6
      data lambda / 1, 1, 1,   ! 1 = e1+ singlet
     +              0, 0, 0,   ! 2 = e2+ updown ns+
     +              0, 0, 0,   ! 3 = e3+ strange ns+
     +             -3, 0, 0,   ! 4 = e4+ charm ns+
     +              0,-4, 0,   ! 5 = e5+ bottom ns+
     +              0, 0,-5,   ! 6 = e6+ top ns+
     +              0, 0, 0,   ! 7 = e1- valence
     +              0, 0, 0,   ! 8 = e2- updown ns-
     +              0, 0, 0,   ! 9 = e3- strange ns- 
     +              0, 0, 0,   !10 = e4- charm ns- 
     +              0, 0, 0,   !11 = e5- bottom ns- 
     +              0, 0, 0 /  !12 = e6- top ns- 

C--   Setup a flavor map
      if(.not.Lnfmap8) call sqcNfTab(0)
C--   Setup kinematic plane
      if(.not.Levcut8) call sqcEvCuts(0)      

C--   VFNS checks
      if(nfix6.eq.0) then
C--     To calculate ns jumps in NNLO we need the singlet at nf-1. This 
C--     implies that for charm a 3-flavor singlet is used. A 3-flavor
C--     singlet is only guaranteed to be available when the charm threshold
C--     is at the second grid point or larger.
        if(itchm2.lt.2 .and. iord.eq.3) stop 
     +     'sqcEvolve: itchm2 .lt. 2 not allowed in NNLO ---> STOP'
C--     Furthermore it0 must be below the charm threshold
        if(it0.ge.itchm2) stop 
     +     'sqcEvolve: it0 at or above itchm2 ---> STOP'
      endif

C--   Pdf indices in main storage 
      idg    =  0
      idf    =  1
      iudpl  =  2
      isspl  =  3
      ichpl  =  4
      ibopl  =  5
      itopl  =  6
      idv    =  7
      iudmi  =  8
      issmi  =  9
      ichmi  = 10
      ibomi  = 11
      itomi  = 12
      
      epsm = 0.D0
      
C--   iz range
      iz1 = izfit2(it1)
      iz2 = izfit2(it2)

C--   Loop over subgrids
      do ig = 1,nyg2
      
C--     Upper y index in subgrid      
        nyg = iqcIyMaxG(iymac2,ig)

C--  A. Singlet/gluon evolution
C--     -----------------------
        iftmp  = -4  !singlet subgrid identifier
        igtmp  = -3  !gluon   subgrid identifier
C--     Copy starting values at it=0 from G0 to Gi 
        call sqcG0toGi(itype,idf,iftmp,ig,nyg,0)
        call sqcG0toGi(itype,idg,igtmp,ig,nyg,0)
C--     Convert starting values at it=0 from F to A
        call sqcGiFtoA(itype,iftmp,iftmp,nyg,0,0)
        call sqcGiFtoA(itype,igtmp,igtmp,nyg,0,0)
C--     Evolve
        call sqcGridsg(itype,
     +                 iftmp,igtmp,ig,iord,it0,it1,it2,eps,ierr)
        if(ierr.ne.0) return
        epsm = max(epsm,eps)

C--      Convert A in Gi to F and copy to  G0
         call sqcAitoF0(itype,iftmp,ig,nyg,iz1,iz2,idf)
         call sqcAitoF0(itype,igtmp,ig,nyg,iz1,iz2,idg)

C--  B. NS+ evolution
C--     -------------
        ityp   = 1      !NS+

C--  B1 ud+ evolution
C--     -------------
        ipdf    = iudpl  !ud+
        itemp   = -2     !subgrid identifier (temporary storage)
        dlam(4) = lambda(4,2)
        dlam(5) = lambda(5,2)
        dlam(6) = lambda(6,2)
C--     Copy startvalues at it=0 from G0 to Gi
        call sqcG0toGi(itype,ipdf,itemp,ig,nyg,0)
C--     Convert starting values at it=0 from F to A
        call sqcGiFtoA(itype,itemp,itemp,nyg,0,0)
C--     Evolve
        call sqcGridns(itype,
     +       itemp,iftmp,igtmp,ityp,ig,iord,dlam,it0,it1,it2,eps,ierr)
        if(ierr.ne.0) return
        epsm = max(epsm,eps)
C--     Convert A in Gi to F and copy to  G0
        call sqcAitoF0(itype,itemp,ig,nyg,iz1,iz2,ipdf)        

C--  B2 s+ evolution
C--     ------------
        ipdf    = isspl  !s+
        itemp   = -2     !subgrid identifier (temporary storage)
        dlam(4) = lambda(4,3)
        dlam(5) = lambda(5,3)
        dlam(6) = lambda(6,3)
C--     Copy startvalues at it=0 from G0 to Gi
        call sqcG0toGi(itype,ipdf,itemp,ig,nyg,0)
C--     Convert starting values at it=0 from F to A
        call sqcGiFtoA(itype,itemp,itemp,nyg,0,0)
C--     Evolve
        call sqcGridns(itype,
     +       itemp,iftmp,igtmp,ityp,ig,iord,dlam,it0,it1,it2,eps,ierr)
        if(ierr.ne.0) return
        epsm = max(epsm,eps)
C--     Convert A in Gi to F and copy to  G0
        call sqcAitoF0(itype,itemp,ig,nyg,iz1,iz2,ipdf)         

C--  B3 c+ evolution
C--     ------------
        itc     = max(it1,itchm2)
        ipdf    = ichpl  !c+
        itemp   = -2     !subgrid identifier (temporary storage)
        dlam(4) = lambda(4,4)
        dlam(5) = lambda(5,4)
        dlam(6) = lambda(6,4)
        if(nfix6.eq.0) then
C--       VFNS: Copy singlet to Gi
          call sqcNStart(itype,
     +                   itemp,iftmp,iftmp,igtmp,ig,iord,dlam(4),itc)
C--       VFNS: Evolve
          if(itc.lt.it2) then
            call sqcGridns(itype,
     +       itemp,iftmp,igtmp,ityp,ig,iord,dlam,itc,itc,it2,eps,ierr)
            if(ierr.ne.0) return
            epsm = max(epsm,eps)
          endif
        elseif(nfix6.ge.4) then
C--       FFNS: Copy startvalues at it=0 from G0 to Gi
          call sqcG0toGi(itype,ipdf,itemp,ig,nyg,0)
C--       FFNS: Convert starting values at it=0 from F to A
          call sqcGiFtoA(itype,itemp,itemp,nyg,0,0)
C--       FFNS: Evolve
          call sqcGridns(itype,
     +       itemp,iftmp,igtmp,ityp,ig,iord,dlam,it0,it1,it2,eps,ierr)
          if(ierr.ne.0) return
          epsm = max(epsm,eps)
        else
C--       FFNS: copy singlet to Gi
          call sqcPdfCop(itype,iftmp,itemp)
*mb          call sqcT1toT2(itype,iftmp,itemp,1,nyg,iz1,iz2)          
        endif
C--     Convert A in Gi to F and copy to  G0
        call sqcAitoF0(itype,itemp,ig,nyg,iz1,iz2,ipdf)         

C--  B4 b+ evolution
C--     ------------
        itb     = max(it1,itbot2)
        ipdf    = ibopl  !b+
        itemp   = -2     !subgrid identifier (temporary storage)
        dlam(4) = lambda(4,5)
        dlam(5) = lambda(5,5)
        dlam(6) = lambda(6,5)
        if(nfix6.eq.0) then
C--       VFNS: Copy singlet to Gi
          call sqcNStart(itype,
     +                   itemp,iftmp,iftmp,igtmp,ig,iord,dlam(5),itb)
C--       VFNS: Evolve
          if(itb.lt.it2) then
            call sqcGridns(itype,
     +       itemp,iftmp,igtmp,ityp,ig,iord,dlam,itb,itb,it2,eps,ierr)
            if(ierr.ne.0) return
            epsm = max(epsm,eps)
          endif
        elseif(nfix6.ge.5) then
C--       FFNS: Copy startvalues at it=0 from G0 to Gi
          call sqcG0toGi(itype,ipdf,itemp,ig,nyg,0)
C--       FFNS: Convert starting values at it=0 from F to A
          call sqcGiFtoA(itype,itemp,itemp,nyg,0,0)
C--       FFNS: Evolve
          call sqcGridns(itype,
     +       itemp,iftmp,igtmp,ityp,ig,iord,dlam,it0,it1,it2,eps,ierr)
          if(ierr.ne.0) return
          epsm = max(epsm,eps)
        else
C--       FFNS: copy singlet to Gi
          call sqcPdfCop(itype,iftmp,itemp)
*mb          call sqcT1toT2(itype,iftmp,itemp,1,nyg,iz1,iz2)
        endif
C--     Convert A in Gi to F and copy to  G0
        call sqcAitoF0(itype,itemp,ig,nyg,iz1,iz2,ipdf)         

C--  B5 t+ evolution
C--     ------------
        itt     = max(it1,ittop2)
        ipdf    = itopl  !t+
        itemp   = -2     !subgrid identifier (temporary storage)
        dlam(4) = lambda(4,6)
        dlam(5) = lambda(5,6)
        dlam(6) = lambda(6,6)
        if(nfix6.eq.0) then
C--       VFNS: Copy singlet to Gi
          call sqcNStart(itype,
     +                   itemp,iftmp,iftmp,igtmp,ig,iord,dlam(6),itt)
C--       VFNS: Evolve
          if(itt.lt.it2) then 
            call sqcGridns(itype,
     +       itemp,iftmp,igtmp,ityp,ig,iord,dlam,itt,itt,it2,eps,ierr)
            if(ierr.ne.0) return
            epsm = max(epsm,eps)
          endif
        elseif(nfix6.ge.6) then
C--       FFNS: Copy startvalues at it=0 from G0 to Gi
          call sqcG0toGi(itype,ipdf,itemp,ig,nyg,0)
C--       FFNS: Convert starting values at it=0 from F to A
          call sqcGiFtoA(itype,itemp,itemp,nyg,0,0)
C--       FFNS: Evolve
          call sqcGridns(itype,
     +       itemp,iftmp,igtmp,ityp,ig,iord,dlam,it0,it1,it2,eps,ierr)
          if(ierr.ne.0) return
          epsm = max(epsm,eps)
        else
C--       FFNS: copy singlet to Gi
          call sqcPdfCop(itype,iftmp,itemp)
*mb          call sqcT1toT2(itype,iftmp,itemp,1,nyg,iz1,iz2)
        endif
C--     Convert A in Gi to F and copy to  G0
        call sqcAitoF0(itype,itemp,ig,nyg,iz1,iz2,ipdf)         

C--  C. NSV evolution
C--     -------------
        ityp    = 3      !NSV
        ipdf    = idv    !valence
        ivtmp   = -3     !subgrid identifier (temporary storage)
        dlam(4) = lambda(4,7)
        dlam(5) = lambda(5,7)
        dlam(6) = lambda(6,7)
C--     Copy startvalues at it=0 from G0 to Gi
        call sqcG0toGi(itype,ipdf,ivtmp,ig,nyg,0)
C--     Convert starting values at it=0 from F to A
        call sqcGiFtoA(itype,ivtmp,ivtmp,nyg,0,0)
C--     Evolve
        call sqcGridns(itype,
     +        ivtmp,0,0,ityp,ig,iord,dlam,it0,it1,it2,eps,ierr)
        if(ierr.ne.0) return
        epsm = max(epsm,eps)
C--     Convert A in Gi to F and copy to  G0
        call sqcAitoF0(itype,ivtmp,ig,nyg,iz1,iz2,ipdf)         

C--  D. NS- evolution
C--     -------------
        ityp   = 2      !NS-

C--  D1 ud- evolution
C--     -------------
        ipdf    = iudmi  !ud-
        itemp   = -2     !subgrid identifier (temporary storage)
        dlam(4) = lambda(4,8)
        dlam(5) = lambda(5,8)
        dlam(6) = lambda(6,8)
C--     Copy startvalues at it=0 from G0 to Gi
        call sqcG0toGi(itype,ipdf,itemp,ig,nyg,0)
C--     Convert starting values at it=0 from F to A
        call sqcGiFtoA(itype,itemp,itemp,nyg,0,0)
C--     Evolve
        call sqcGridns(itype,
     +        itemp,0,0,ityp,ig,iord,dlam,it0,it1,it2,eps,ierr)
        if(ierr.ne.0) return
        epsm = max(epsm,eps)
C--     Convert A in Gi to F and copy to  G0
        call sqcAitoF0(itype,itemp,ig,nyg,iz1,iz2,ipdf) 

C--  D2 s- evolution
C--     ------------
        ipdf    = issmi  !s-
        itemp   = -2     !subgrid identifier (temporary storage)
        dlam(4) = lambda(4,9)
        dlam(5) = lambda(5,9)
        dlam(6) = lambda(6,9)
C--     Copy startvalues from G0 to Gi
        call sqcG0toGi(itype,ipdf,itemp,ig,nyg,0)
C--     Convert starting values at it=0 from F to A
        call sqcGiFtoA(itype,itemp,itemp,nyg,0,0)
C--     Evolve
        call sqcGridns(itype,
     +        itemp,0,0,ityp,ig,iord,dlam,it0,it1,it2,eps,ierr)
        if(ierr.ne.0) return
        epsm = max(epsm,eps)
C--     Convert A in Gi to F and copy to  G0
        call sqcAitoF0(itype,itemp,ig,nyg,iz1,iz2,ipdf)         

C--  D3 c- evolution (if VFNS then only in NLLO)
C--     ----------------------------------------
        itc     = max(it1,itchm2)
        ipdf    = ichmi  !c-
        itemp   = -2     !subgrid identifier (temporary storage)
        dlam(4) = lambda(4,10)
        dlam(5) = lambda(5,10)
        dlam(6) = lambda(6,10)
        if(nfix6.eq.0) then
C--       VFNS: Copy valence to Gi
          call sqcNStart(itype,itemp,ivtmp,0,0,ig,iord,dlam(4),itc)
C--       VFNS: Evolve
          if(itc.lt.it2 .and. iord.eq.3) then
            call sqcGridns(itype,
     +          itemp,0,0,ityp,ig,iord,dlam,itc,itc,it2,eps,ierr)
            if(ierr.ne.0) return
            epsm = max(epsm,eps)
          endif
        elseif(nfix6.ge.4) then
C--       FFNS: Copy startvalues at it=0 from G0 to Gi
          call sqcG0toGi(itype,ipdf,itemp,ig,nyg,0)
C--       FFNS: Convert starting values at it=0 from F to A
          call sqcGiFtoA(itype,itemp,itemp,nyg,0,0)
C--       FFNS: Evolve
          call sqcGridns(itype,
     +          itemp,0,0,ityp,ig,iord,dlam,it0,it1,it2,eps,ierr)
          if(ierr.ne.0) return
          epsm = max(epsm,eps)
        else
C--       FFNS: copy valence to Gi
          call sqcPdfCop(itype,ivtmp,itemp)
*mb          call sqcT1toT2(itype,ivtmp,itemp,1,nyg,iz1,iz2)
        endif
C--     Convert A in Gi to F and copy to  G0
        call sqcAitoF0(itype,itemp,ig,nyg,iz1,iz2,ipdf)         

C--  D4 b- evolution (if VFNS then only in NLLO)
C--     ----------------------------------------
        itb     = max(it1,itbot2)
        ipdf    = ibomi  !b-
        itemp   = -2     !subgrid identifier (temporary storage)
        dlam(4) = lambda(4,11)
        dlam(5) = lambda(5,11)
        dlam(6) = lambda(6,11)
        if(nfix6.eq.0) then
C--       VFNS: Copy valence to Gi
          call sqcNStart(itype,itemp,ivtmp,0,0,ig,iord,dlam(5),itb)
C--       VFNS: Evolve
          if(itb.lt.it2 .and. iord.eq.3) then 
            call sqcGridns(itype,
     +          itemp,0,0,ityp,ig,iord,dlam,itb,itb,it2,eps,ierr)
            if(ierr.ne.0) return
            epsm = max(epsm,eps)
          endif
        elseif(nfix6.ge.5) then
C--       FFNS: Copy startvalues at it=0 from G0 to Gi
          call sqcG0toGi(itype,ipdf,itemp,ig,nyg,0)
C--       FFNS: Convert starting values at it=0 from F to A
          call sqcGiFtoA(itype,itemp,itemp,nyg,0,0)
C--       FFNS: Evolve
          call sqcGridns(itype,
     +          itemp,0,0,ityp,ig,iord,dlam,it0,it1,it2,eps,ierr)
          if(ierr.ne.0) return
          epsm = max(epsm,eps)
        else
C--       FFNS: copy valence to Gi
          call sqcPdfCop(itype,ivtmp,itemp)
*mb          call sqcT1toT2(itype,ivtmp,itemp,1,nyg,iz1,iz2)
        endif
C--     Convert A in Gi to F and copy to  G0
        call sqcAitoF0(itype,itemp,ig,nyg,iz1,iz2,ipdf)         
        
C--  D5 t- evolution (if VFNS then only in NLLO)
C--     ----------------------------------------
        itt     = max(it1,ittop2)
        ipdf    = itomi  !t-
        itemp   = -2     !subgrid identifier (temporary storage)
        dlam(4) = lambda(4,12)
        dlam(5) = lambda(5,12)
        dlam(6) = lambda(6,12)
        if(nfix6.eq.0) then
C--       VFNS: Copy valence to Gi
          call sqcNStart(itype,itemp,ivtmp,0,0,ig,iord,dlam(6),itt)
C--       VFNS: Evolve
          if(itt.lt.it2 .and. iord.eq.3) then
            call sqcGridns(itype,
     +          itemp,0,0,ityp,ig,iord,dlam,itt,itt,it2,eps,ierr)
            if(ierr.ne.0) return
            epsm = max(epsm,eps)
          endif
        elseif(nfix6.ge.6) then
C--       FFNS: Copy startvalues at it=0 from G0 to Gi
          call sqcG0toGi(itype,ipdf,itemp,ig,nyg,0)
C--       FFNS: Convert starting values at it=0 from F to A
          call sqcGiFtoA(itype,itemp,itemp,nyg,0,0)
C--       FFNS: Evolve
          call sqcGridns(itype,
     +          itemp,0,0,ityp,ig,iord,dlam,it0,it1,it2,eps,ierr)
          if(ierr.ne.0) return
          epsm = max(epsm,eps)
        else
C--       FFNS: copy valence to Gi
          call sqcPdfCop(itype,ivtmp,itemp)
*mb          call sqcT1toT2(itype,ivtmp,itemp,1,nyg,iz1,iz2)
        endif
C--     Convert A in Gi to F and copy to  G0
        call sqcAitoF0(itype,itemp,ig,nyg,iz1,iz2,ipdf)         
        
C--   End of loop over subgrids
      enddo

      return
      end

C     ================================================================
C     Nonsinglet evolution routines
C     ================================================================

C     ================================================================
      subroutine sqcNStart(itype,idout,idin,ids,idg,iyg,iord,dlam,it0)
C     ================================================================

C--   Copy singlet or valence and set proper startvalue for VFNS heavy quarks

C--   idout    (in)    output nonsinglet table
C--   idin     (in)    input singlet or valence table
C--   ids      (in)    singlet table
C--   idg      (in)    gluon table
C--   iyg      (in)    subgrid index
C--   iord     (in)    1 = LO, 2 = NLO, 3 = NNLO
C--   dlam     (in)    heavy flavor weight
C--   it0      (in)    starting point

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qmaps8.inc'
      
C--   Copy input to output table (input is singlet or valence)
      call sqcPdfCop(itype,idin,idout)
*mb      call sqcT1toT2(itype,idin,idout,1,nyy2(iyg),1,nzz2)
C--   Find out in which region it0 falls
      ibin0 = 0
      do i = 1,nlist8
        if(it0.eq.itmin8(i)) ibin0 = i
      enddo
C--   Bin not found, dont set startvalue
      if(ibin0.eq.0) return
C--   Calculate iz0; izmin-itmin gives the offset between iz and it
      iz0 = it0 + izmin8(ibin0)-itmin8(ibin0)
C--   LO, NLO: copy starting value to bin 0 and then exit
      if(iord.le.2) then
        call sqcPCopjj(itype,idout,iz0,idout,0)
        return
      endif
C--   NNLO: copy starting value to bin 0 but also apply jumps
C--   Start scale must be at least in the second region
      if(ibin0.eq.1) stop 'sqcNStart: NNLO ibin0 .eq. 1 ---> STOP' 
C--   Where to get the startvalue from
      iz1 = izmax8(ibin0-1) !Thats why ibin0 must be larger than 1
C--   Add discontinuity and store in iz0
      call sqcNSjup
     +        (itype,idout,ids,idg,iyg,iord,dlam,nyy2(iyg),iz1,iz0)
C--   Copy starting value to bin 0 (thats where sqcGridns looks for it)
      call sqcPCopjj(itype,idout,iz0,idout,0)

      return
      end

C     ===============================================================
      subroutine sqcGridns(itype,
     +          ipdf,ids,idg,ityp,iyg,iord,dlam,it0,it1,it2,eps,ierr)
C     ===============================================================

C--   Steering routine for nonsinglet evolution on an equidistant y-subgrid
C--   This routine handles crossing of the flavor thresholds
C--   This routine does not set the starting values
C--
C--   1. it1 must be at lower boundary or at one of the thresholds
C--   2. it2 must be at upper boundary
C--   3. it0 must be >= it1 and <= it2
C--
C--   ipdf      (in) nonsinglet pdf index
C--   ids       (in) singlet index
C--   idg       (in) gluon index 
C--   ityp      (in) 1=NS+, 2=NS-, 3=NSV
C--   iyg       (in) subgrid index 1,...,nyg2
C--   iord      (in) 1=LO , 2=NLO, 3=NNLO
C--   dlam(4:6) (in) heavy quark weights in NNLO nonsinglet jumps
C--   it0       (in) starting point          
C--   it1       (in) lower limit of evolution
C--   it2       (in) upper limit of evolution
C--   eps      (out) max deviation quad - lin interpolation
C--   ierr     (out) 1 = no alphas available at it1
C--
C--   Remark: ids, idg and dlam are used to calculate the NNLO jumps
C--           If, on input, ids = idg then  dont calculate the heavy
C--           flavor contribution to the jumps

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'qmaps8.inc'

      dimension dlam(4:6)
      dimension izu1(4),izu2(4),nflu(4),izd1(4),izd2(4),nfld(4)
      
C--   Setup a flavor map for up and downward evolution
      if(.not.Lnfmap8) call sqcNfTab(0)
C--   Fill tables with (alpha_s/2pi)^iord if not done already
      if(.not.Lastab8) call sqcAlfTab(iord)
C--   Check if alphas available at lowest t-bin
      if(it1.lt.itlow8) then
        ierr = 1
        return
      endif
      ierr = 0
      eps  = 0.D0

C--   Get evolution limits (kinda joblist of how to proceed)
      call sqcEvLims(it0,it1,it2,
     +               izu1,izu2,nflu,nup,izd1,izd2,nfld,ndn,idl,idu)
     
C--   Get upper y-index
      iymax = iqcIyMaxG(iymac2,iyg)     

C--   Upward evolutions
      do i = 1,nup
        if(i.eq.1) then
C--       Copy starting value from bin 0 to bin iz
          call sqcPCopjj(itype,ipdf,0,ipdf,izu1(i))
        else
C--       Copy starting value from previous evolution and add discontinuity
          dlm = dlam(nflu(i))
          call sqcNSjup(itype,
     +         ipdf,ids,idg,iyg,iord,dlm,iymax,izu2(i-1),izu1(i))
*mb     +         ipdf,ids,idg,iyg,iord,dlm,nyy2(iyg),izu2(i-1),izu1(i))
        endif
C--     Evolve upward with fixed nf
        call sqcNSevnf(itype,
     +                 ipdf,ityp,iyg,iord,nflu(i),izu1(i),izu2(i))
      enddo

C--   Downward evolutions (always w/linear interpolation)
      do i = 1,ndn
      
        if(i.eq.1) then
C--       Copy starting value from bin 0 to bin iz
          call sqcPCopjj(itype,ipdf,0,ipdf,izd1(i))
        else
C--       Copy starting value from previous evolution and add discontinuity
          dlm = dlam(nfld(i))
          call sqcNSjdn(itype,
     +         ipdf,ids,idg,iyg,iord,dlm,iymax,izd2(i-1),izd1(i))
*mb     +         ipdf,ids,idg,iyg,iord,dlm,nyy2(iyg),izd2(i-1),izd1(i))
        endif
        
        if(ioy2.eq.2) then
        
C--       Evolve downward with current ioy2 = lin
          call sqcNSevnf(itype,
     +               ipdf,ityp,iyg,iord,nfld(i),izd1(i),izd2(i))
     
        elseif(ioy2.eq.3 .and. niter6.lt.0) then
        
C--       Evolve downward with current ioy2 = quad
          call sqcNSevnf(itype,
     +               ipdf,ityp,iyg,iord,nfld(i),izd1(i),izd2(i))
          
        elseif(ioy2.eq.3 .and. niter6.eq.0) then
        
C--       Switch to linear interpolation
          call sqcGiQtoL(itype,ipdf,ipdf,iymax,izd1(i),izd1(i))        
C--       Evolve downward with linear interpolation
          call sqcNSevnf(itype,
     +               ipdf,ityp,iyg,iord,nfld(i),izd1(i),izd2(i))
C--       Switch back to quad interpolation
          call sqcGiLtoQ(itype,ipdf,ipdf,iymax,izd1(i),izd1(i))        

        elseif(ioy2.eq.3 .and. niter6.gt.0) then
        
C--       Switch to linear interpolation
          call sqcGiQtoL(itype,ipdf,ipdf,iymax,izd1(i),izd1(i))
C--       Evolve downward with lin interpolation
          call sqcNSevnf(itype,
     +               ipdf,ityp,iyg,iord,nfld(i),izd1(i),izd2(i))
C--       Remember startvalue
          call sqcNSStoreStart(itype,ipdf,1,iymax,izd1(i))
*mb          call sqcNSStoreStart(itype,ipdf,1,nyy2(iyg),izd1(i))
C--       Switch to quad interpolation
          call sqcGiLtoQ(itype,ipdf,ipdf,iymax,izd2(i),izd2(i))
C--       Evolve upward with quad interpolation
          call sqcNSevnf(itype,
     +               ipdf,ityp,iyg,iord,nfld(i),izd2(i),izd1(i))
C--       Iteration loop
          do iter = 0,niter6
C--         Switch to linear interpolation
            call sqcGiQtoL(itype,ipdf,ipdf,iymax,izd1(i),izd1(i))
C--         New starting value
            call sqcNSNewStart(itype,ipdf,1,iymax,izd1(i),eps)
*mb            call sqcNSNewStart(itype,ipdf,1,nyy2(iyg),izd1(i),eps)
C--         Finished?
            if(iter.eq.niter6) then
C--           Restore startvalue
              call sqcNSRestoreStart(itype,ipdf,1,iymax,izd1(i))
*mb              call sqcNSRestoreStart(itype,ipdf,1,nyy2(iyg),izd1(i))
C--           Switch to quad interpolation
              call sqcGiLtoQ(itype,ipdf,ipdf,iymax,izd1(i),izd1(i))
            else
C--           Evolve downward with lin interpolation
              call sqcNSevnf(itype,
     +                   ipdf,ityp,iyg,iord,nfld(i),izd1(i),izd2(i))
C--           Switch to quad interpolation
              call sqcGiLtoQ(itype,ipdf,ipdf,iymax,izd2(i),izd2(i))
C--           Evolve upward with quad interpolation
              call sqcNSevnf(itype,
     +                   ipdf,ityp,iyg,iord,nfld(i),izd2(i),izd1(i))
            endif 
          enddo
        endif
        
C--   End of loop over downward evolutions        
      enddo
      
      eps0 = dqcGetEps(itype,ipdf,iymax,it0)
      eps1 = dqcGetEps(itype,ipdf,iymax,it1)
      eps2 = dqcGetEps(itype,ipdf,iymax,it2)
      eps  = max(eps0,eps1,eps2)

      return
      end

C     =========================================================
      subroutine sqcNSevnf(itype,ipdf,ityp,iyg,iord,nf,iz1,iz2)
C     =========================================================

C--   Non-singlet evolution from iz1 to iz2 at fixed nf
C--
C--   Input:  ipdf = index of pdf to be evolved
C--           ityp = 1=NS+, 2=NS-, 3=NSV
C--           iyg  = subgrid index 1,...,nyg2
C--           iord = 1=LO,  2=NLO, 3=NNLO
C--           nf   = number of flavours
C--           iz1  = start of evolution
C--           iz2  = end of evolution
C--
C--   Output:  Spline coefficients in /qstor7/stor7 (linear store)
C--
C--   NB: the routines sqcNSmult and sqcNSeqs can be found in qcdutil.f

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'qmaps8.inc'

      dimension idwt(3)
C--               NS+ NS- NSV
      data idwt /  5,  6,  7  /

      dimension sbar(mxx0),ssum(mxx0),vmat(mxx0),bvec(mxx0),hvec(mxx0)
      
*mb
*mb      logical first,lprint
*mb      save first,lprint,nflast
*mb      data first /.true./
*mb      
*mb      if(first) then
*mb        first  = .false.
*mb        lprint = .true.
*mb        nflast = nf 
*mb      endif
*mb      
*mb      if(nf.ne.nflast .and. lprint) then
*mb        write(6,*) 'EVOL: change number of flavours'
*mb        lprint = .false.
*mb      endif
*mb                                
      
C--   Initialization
C--   --------------
C--   Interpolation index
      idk     = ioy2-1
C--   Weight table offset     
      ifirst7 = ifst7(itype,ioy2)
C--   Direction of evolution (isign) and first point after it1 (next)
      isign = 1
      next  = iz1+1
      if(iz2.lt.iz1) then 
        isign = -1
        next  = iz1-1
      endif
C--   Setup sbar
      do i = 1,nyy2(iyg)
        sbar(i) = 0.D0
        ssum(i) = 0.D0
      enddo
C--   Index limits of region iyg on subgrid iyg
      iy1 = 1
      iy2 = iqcIyMaxG(iymac2,iyg)
*      iy2 = nyy2(iyg)

C--   Calculate vector b at input scale iz1
C--   -------------------------------------
C--   Grid spacing delta = z(next)-z (divided by 2)
      delt = 0.5*abs(zgrid2(next)-zgrid2(iz1))
C--   Transformation matrix divided by delta (Sbar)        
      do i = 1,nmaty2
        sbar(i) = smaty2(i)/delt
      enddo
C--   Find t-index
      it1 = itfiz2(iz1)      
C--   Weight matrix and  V = Sbar+W at t1
      do iy = 1,iy2
        vmat(iy) = sbar(iy)
      enddo
      do k = 1,iord
        id = idPij7(idwt(ityp),k,itype)
        as = antab8(iz1,k)
        ia = iqcPijklm(1,it1,nf,iyg,id)-1
        do iy = 1,iy2
          ia = ia+1
          vmat(iy) = vmat(iy) + isign*as*stor7(ia)
        enddo  
      enddo      
C--   Address of pdf(iy=1,iz1,ipdf) in the store
      iadr = iqcPdfIjkl(1,iz1,ipdf,itype)
C--   Calculate V a = b
      call sqcNSmult(vmat,iy2,stor7(iadr),bvec,iy2)

C--   Evolution loop over t
C--   ---------------------
      do iz = next,iz2,isign
C--     Find t-index
        it = itfiz2(iz)
*mbC--     Yes/no apply roots cut, thats the question
*mb        if(isign.eq.1) then
*mb          iy2 = iqcIyMaxG(ismac2(iz),iyg)  !apply roots cut 4upward
*mb        else
*mb          iy2 = iqcIyMaxG(iymac2,iyg)  !not apply roots cut 4dnward
*mb        endif                 
C--     Weight matrix and V matrix at t
        do iy = 1,iy2
          vmat(iy) = sbar(iy)
        enddo
        do k = 1,iord
          id = idPij7(idwt(ityp),k,itype)
          as = antab8(iz,k)
          ia = iqcPijklm(1,it,nf,iyg,id)-1
          do iy = 1,iy2
            ia = ia+1
            vmat(iy) = vmat(iy) - isign*as*stor7(ia)
          enddo  
        enddo
C--     Address of pdf(1,iz,ipdf) in the store
        iadr = iqcPdfIjkl(1,iz,ipdf,itype)
C--     Solve V a = b
        call sqcNSeqs(vmat,iy2,stor7(iadr),bvec,iy2)
C--     Update b for the next iteration: not at last iteration of the loop
        if(iz.ne.iz2) then
C--       Grid spacing delta = z(next)-z (divided by 2)
          delt = 0.5*abs(zgrid2(iz+isign)-zgrid2(iz))
C--       Sum of current and next sbar; store next sbar       
          do i = 1,nmaty2
            sbnext  = smaty2(i)/delt
            ssum(i) = sbar(i)+sbnext
            sbar(i) = sbnext
          enddo
C--       Calculate Ssum a = h
          call sqcNSmult(ssum,nmaty2,stor7(iadr),hvec,iy2)
C--       Update b
          do iy = 1,iy2
            bvec(iy) = hvec(iy)-bvec(iy)
          enddo
        endif
      enddo

      return
      end

C     ===============================================================
      subroutine 
     +        sqcNSjup(itype,idq,ids,idg,iyg,iord,dlam,ny,izin,izout)
C     ===============================================================

C--   Calculate NNLO jump in nonsinglet for upward evolution
C--
C--   Input:  idq    = index of nonsinglet pdf
C--           ids    = index of singlet pdf
C--           idg    = index of gluon pdf
C--           iyg    = subgrid index 1,...,nyg2         
C--           iord   = 1=LO,  2=NLO, 3=NNLO
C--           dlam   = weight of heavy quark discontinuity
C--           ny     = upper index limit yloop
C--           izin   = z-grid index where A(nf) is to be found 
C--           izout  = z-grid index where A(nf+1) will be stored
C--
C--   Output: Spline coefficients in /qstor7/stor7 (linear store)
C--
C--   Remark: (1) Acts as a do-nothing when iord < 3 (except copy)
C--           (2) When ids = idg the heavy flavor contribution is set to zero

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
*      include 'qwtab7.inc'
      include 'qpdfs7.inc'
      include 'qmaps8.inc'

      data idqq, idhq, idhg / 3, 4, 5 /

      dimension wmatqq(mxx0), wmathq(mxx0), wmathg(mxx0)
      dimension qjump(mxx0) , sjump(mxx0) , gjump(mxx0)
      
C--   Do nothing except copy izin to izout
      if(iord.lt.3) then
        call sqcPCopjj(itype,idq,izin,idq,izout)
        return
      endif

C--   Weight table offset      
      ifirst7 = ifst7(itype,ioy2)

C--   Alphas to use
*      assq = antab8(izout,2)
      assq = antab8(izout,0)*antab8(izout,0)

C--   Find t-index   
      itin = itfiz2(izin)

C--   Calculate quark weight Wqq
      nf = 3  !dummy variable since the A-weights do not depend on nf
      ia = iqcPijklm(1,itin,nf,iyg,idAij7(idqq))-1
      do iy = 1,ny
        wmatqq(iy) = assq*stor7(ia+iy)
      enddo
C--   Add transformation matrix to Wqq
      do i = 1,nmaty2
        wmatqq(i) = wmatqq(i)+smaty2(i)
      enddo
C--   Calculate qjump = W*a and store in buffer
      iaq   = iqcPdfIjkl(1,izin,idq,itype) !address of input nonsinglet
      call sqcNSmult(wmatqq,ny,stor7(iaq),qjump,ny)
C--   Add heavy flavor but only if ids .ne. idg
      if(ids.ne.idg) then
C--     Weight matrix
        iaq = iqcPijklm(1,itin,nf,iyg,idAij7(idhq))-1
        iag = iqcPijklm(1,itin,nf,iyg,idAij7(idhg))-1
        do iy = 1,ny
          wmathq(iy) = assq*stor7(iaq+iy)
          wmathg(iy) = assq*stor7(iag+iy)
        enddo
        ias = iqcPdfIjkl(1,izin,ids,itype) !address of input singlet
        iag = iqcPdfIjkl(1,izin,idg,itype) !address of input gluon
C--     Calculate jump = W*a and store in buffers sjump and gjump
        call sqcNSmult(wmathq,ny,stor7(ias),sjump,ny)
        call sqcNSmult(wmathg,ny,stor7(iag),gjump,ny)
C--     Add weighted heavy quark jump
        do i = 1,ny
          qjump(i) = qjump(i) + dlam * (sjump(i)+gjump(i))
        enddo
      endif
C--   Calculate anew by solving S*anew = qnew
      iaq  = iqcPdfIjkl(1,izout,idq,itype) !address of output nonsinglet
      call sqcNSeqs(smaty2,nmaty2,stor7(iaq),qjump,ny)

      return
      end

C     ===============================================================
      subroutine 
     +        sqcNSjdn(itype,idq,ids,idg,iyg,iord,dlam,ny,izin,izout)
C     ===============================================================

C--   Calculate NNLO jump in nonsinglet for downward evolution
C--
C--   Input:  idq    = index of nonsinglet pdf
C--           ids    = index of singlet pdf
C--           idg    = index of gluon pdf
C--           iyg    = subgrid index 1,...,nyg2         
C--           iord   = 1=LO,  2=NLO, 3=NNLO
C--           dlam   = weight of heavy quark discontinuity
C--           ny     = upper index limit yloop
C--           izin   = z-grid index where A(nf+1) is to be found 
C--           izout  = z-grid index where A(nf) will be stored
C--
C--   Output: Spline coefficients in /qstor7/stor7 (linear store)
C--
C--   Remark: (1) Acts as a do-nothing when iord < 3
C--           (2) When ids = idg the heavy flavor contribution is set to zero

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
*      include 'qwtab7.inc'
      include 'qpdfs7.inc'
      include 'qmaps8.inc'

      data idqq, idhq, idhg / 3, 4, 5 /

      dimension wmatqq(mxx0), wmathq(mxx0), wmathg(mxx0)
      dimension qstor(mxx0) , sstor(mxx0) , gstor(mxx0) 

C--   Do nothing except copy izin to izout
      if(iord.lt.3) then
        call sqcPCopjj(itype,idq,izin,idq,izout)
        return
      endif

C--   Weight table offset      
      ifirst7 = ifst7(itype,ioy2)

C--   Alphas to use
*      assq = antab8(izin,2)
      assq = antab8(izin,0)*antab8(izin,0)
      
C--   Find t-index
      itin = itfiz2(izin)      

C--   Calculate quark weight Wqq
      nf = 3  !dummy variable since the A-weights do not depend on nf
      ia = iqcPijklm(1,itin,nf,iyg,idAij7(idqq))-1
      do iy = 1,ny
        wmatqq(iy) = assq*stor7(ia+iy)
      enddo
C--   Add transformation matrix to Wqq
      do i = 1,nmaty2
        wmatqq(i) = wmatqq(i)+smaty2(i)
      enddo
C--   Calculate qold = S*aold and store in buffer 
      iaq   = iqcPdfIjkl(1,izin,idq,itype) !address of input nonsinglet
      call sqcNSmult(smaty2,nmaty2,stor7(iaq),qstor,ny)
C--   Subtract heavy flavor but only if ids .ne. idg
      if(ids.ne.idg) then
C--     Weight matrix
        iaq = iqcPijklm(iy,itin,nf,iyg,idAij7(idhq))-1
        iag = iqcPijklm(iy,itin,nf,iyg,idAij7(idhg))-1
        do iy = 1,ny
          wmathq(iy) = assq*stor7(iaq+iy)
          wmathg(iy) = assq*stor7(iag+iy)
        enddo
        ias = iqcPdfIjkl(1,izin,ids,itype) !address of input singlet
        iag = iqcPdfIjkl(1,izin,idg,itype) !address of input gluon
C--     Calculate jump = W*a and store in buffers sstor and gstor
        call sqcNSmult(wmathq,ny,stor7(ias),sstor,ny)
        call sqcNSmult(wmathg,ny,stor7(iag),gstor,ny)
C--     Subtract weighted heavy quark jump
        do i = 1,ny
          qstor(i) = qstor(i) - dlam * (sstor(i)+gstor(i))
        enddo
      endif
C--   Calculate anew by solving Wqq*anew = qstor
      iaq  = iqcPdfIjkl(1,izout,idq,itype) !address of output nonsinglet
      call sqcNSeqs(wmatqq,ny,stor7(iaq),qstor,ny)

      return
      end

C     ==================================================
      subroutine sqcNSStoreStart(itype,ipdf,iy1,iy2,iz0)
C     ==================================================

C--   Store  startvalue into a buffer
C--
C--   ipdf   Nonsinglet table
C--   iy1    First yloop index
C--   iy2    Last yloop index
C--   iz0    t-bin containing startvalue  

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'

      common /stbuf/ qtarg(mxx0), gtarg(mxx0),
     +               qlast(mxx0), glast(mxx0)

C--   Calculate base address
      ia  = iqcPdfIjkl(iy1,iz0,ipdf,itype)-1
C--   Store startvalue
      do j = iy1,iy2
        ia = ia+1
        qtarg(j) = stor7(ia)
        qlast(j) = stor7(ia)
      enddo

      return
      end

C     ====================================================
      subroutine sqcNSNewStart(itype,ipdf,iy1,iy2,iz0,eps)
C     ====================================================

C--   Set new startvalue
C--
C--   ipdf   Nonsinglet table
C--   iy1    First yloop index
C--   iy2    Last yloop index
C--   iz0    z-bin containing startvalue
C--   eps    (out) Max |new-original|  

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'

      common /stbuf/ qtarg(mxx0), gtarg(mxx0),
     +               qlast(mxx0), glast(mxx0)

C--   Calculate base address
      ia = iqcPdfIjkl(iy1,iz0,ipdf,itype)-1
C--   New startvalue
      eps = -999.D0
      do j = iy1,iy2
        ia = ia+1
        dif = stor7(ia)-qtarg(j)
        eps = max(eps, abs(dif))
        stor7(ia) = qlast(j)-dif
        qlast(j)  = stor7(ia)
      enddo

      return
      end

C     ====================================================
      subroutine sqcNSRestoreStart(itype,ipdf,iy1,iy2,iz0)
C     ====================================================

C--   Restore  startvalue
C--
C--   ipdf   Nonsinglet table
C--   iy1    First yloop index
C--   iy2    Last yloop index
C--   iz0    z-bin containing startvalue  

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'

      common /stbuf/ qtarg(mxx0), gtarg(mxx0),
     +               qlast(mxx0), glast(mxx0)

C--   Calculate base address
      ia = iqcPdfIjkl(iy1,iz0,ipdf,itype)-1
C--   Restore startvalue
      do j = iy1,iy2
        ia = ia+1
        stor7(ia) = qtarg(j)
      enddo

      return
      end

C     ================================================================
C     Singlet/Gluon evolution routines
C     ================================================================

C     ===========================================================
      subroutine sqcGridsg(itype,
     +                     idf,idg,iyg,iord,it0,it1,it2,eps,ierr)
C     ===========================================================

C--   Steering routine for singlet/gluon evolution on an equidistant subgrid
C--   This routine handles crossing of the flavor thresholds
C--   This routine does not set the starting values
C--
C--   1. it1 must be at lower boundary or at one of the thresholds
C--   2. it2 must be at upper boundary
C--   3. it0 must be >= it1 and <= it2
C--
C--   idf    (in) pdf index for singlet
C--   idg    (in) pdf index for gluon
C--   iyg    (in) subgrid index 1,...,nyg2
C--   iord   (in) 1-LO , 2=NLO, 3=NNLO
C--   it0    (in) starting point           ) 
C--   it1    (in) lower limit of evolution ) it1 <= it0 <= it2
C--   it2    (in) upper limit of evolution )
C--   eps   (out) max deviation quad - lin interpolation
C--   ierr  (out) 1 = no alphas available at it1

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'qmaps8.inc'

      dimension izu1(4),izu2(4),nflu(4),izd1(4),izd2(4),nfld(4)
      
C--   Setup a flavor map for up and downward evolution
      if(.not.Lnfmap8) call sqcNfTab(0)
C--   Fill tables with (alpha_s/2pi)^iord if not done already
      if(.not.Lastab8) call sqcAlfTab(iord)
C--   Check if alphas available at lowest t-bin
      if(it1.lt.itlow8) then
        ierr = 1
        return
      endif
      ierr = 0
      eps  = 0.D0

C--   Get evolution limits (kinda joblist of how to proceed)
      call sqcEvLims(it0,it1,it2,
     +               izu1,izu2,nflu,nup,izd1,izd2,nfld,ndn,ibl,ibu)
     
C--   Get upper y-index
      iymax = iqcIyMaxG(iymac2,iyg)

C--   Upward evolutions
      do i = 1,nup
        if(i.eq.1) then
C--       Copy starting value from bin 0 to bin iz
          call sqcPCopjj(itype,idf,0,idf,izu1(i))
          call sqcPCopjj(itype,idg,0,idg,izu1(i))
        else
C--       Copy starting value from previous evolution and add discontinuity
          call sqcSGjup
*mb     +           (itype,idf,idg,iyg,iord,nyy2(iyg),izu2(i-1),izu1(i))
     +           (itype,idf,idg,iyg,iord,iymax,izu2(i-1),izu1(i))
        endif
C--     Evolve upward with fixed nf
        call sqcSGevnf(itype,idf,idg,iyg,iord,nflu(i),izu1(i),izu2(i))
      enddo

C--   Downward evolutions (always w/linear interpolation) 
      do i = 1,ndn
      
        if(i.eq.1) then
C--       Copy starting value from bin 0 to bin iz
          call sqcPCopjj(itype,idf,0,idf,izd1(i))
          call sqcPCopjj(itype,idg,0,idg,izd1(i))
        else
C--       Copy starting value from previous evolution and add discontinuity
          call sqcSGjdn
*mb     +           (itype,idf,idg,iyg,iord,nyy2(iyg),izd2(i-1),izd1(i))
     +           (itype,idf,idg,iyg,iord,iymax,izd2(i-1),izd1(i))
        endif
        
        if(ioy2.eq.2) then
        
C--       Evolve downward with current ioy2 = lin
          call sqcSGevnf(itype,
     +               idf,idg,iyg,iord,nfld(i),izd1(i),izd2(i)) 
     
        elseif(ioy2.eq.3 .and. niter6.lt.0) then   
                 
C--       Evolve downward with current ioy2 = quad 
          call sqcSGevnf(itype,
     +               idf,idg,iyg,iord,nfld(i),izd1(i),izd2(i))
     
        elseif(ioy2.eq.3 .and. niter6.eq.0) then
        
C--       Switch to linear interpolation                         
          call sqcGiQtoL(itype,idf,idf,iymax,izd1(i),izd1(i))
          call sqcGiQtoL(itype,idg,idg,iymax,izd1(i),izd1(i))
C--       Evolve downward with linear interpolation 
          call sqcSGevnf(itype,
     +               idf,idg,iyg,iord,nfld(i),izd1(i),izd2(i))
C--       Switch back to quad interpolation              
          call sqcGiLtoQ(itype,idf,idf,iymax,izd1(i),izd1(i))
          call sqcGiLtoQ(itype,idg,idg,iymax,izd1(i),izd1(i))
          
       elseif(ioy2.eq.3 .and. niter6.gt.0) then
       
C--       Switch to linear interpolation       
          call sqcGiQtoL(itype,idf,idf,iymax,izd1(i),izd1(i))
          call sqcGiQtoL(itype,idg,idg,iymax,izd1(i),izd1(i))
C--       Evolve downward with lin interpolation          
          call sqcSGevnf(itype,
     +               idf,idg,iyg,iord,nfld(i),izd1(i),izd2(i))
C--       Remember startvalue     
*mb          call sqcSGStoreStart(itype,idf,idg,1,nyy2(iyg),izd1(i))
          call sqcSGStoreStart(itype,idf,idg,1,iymax,izd1(i))
C--       Switch to quad interpolation
          call sqcGiLtoQ(itype,idf,idf,iymax,izd2(i),izd2(i))
          call sqcGiLtoQ(itype,idg,idg,iymax,izd2(i),izd2(i))
C--       Evolve upward with quad interpolation     
          call sqcSGevnf(itype,
     +               idf,idg,iyg,iord,nfld(i),izd2(i),izd1(i))
C--       Iteration loop
          do iter = 0,niter6     
C--         Switch to linear interpolation     
            call sqcGiQtoL(itype,idf,idf,iymax,izd1(i),izd1(i))
            call sqcGiQtoL(itype,idg,idg,iymax,izd1(i),izd1(i))
C--         New starting value          
            call sqcSGNewStart
*mb     +              (itype,idf,idg,1,nyy2(iyg),izd1(i),eps)
     +              (itype,idf,idg,1,iymax,izd1(i),eps)
C--         Finished?
            if(iter.eq.niter6) then
C--           Restore startvalue          
              call sqcSGRestoreStart
*mb     +                      (itype,idf,idg,1,nyy2(iyg),izd1(i))
     +                      (itype,idf,idg,1,iymax,izd1(i))
C--           Switch to quad interpolation              
              call sqcGiLtoQ(itype,idf,idf,iymax,izd1(i),izd1(i))
              call sqcGiLtoQ(itype,idg,idg,iymax,izd1(i),izd1(i))
            else
C--           Evolve downward with lin interpolation   
              call sqcSGevnf(itype,
     +                   idf,idg,iyg,iord,nfld(i),izd1(i),izd2(i))
C--           Switch to quad interpolation     
              call sqcGiLtoQ(itype,idf,idf,iymax,izd2(i),izd2(i))
              call sqcGiLtoQ(itype,idg,idg,iymax,izd2(i),izd2(i))
C--           Evolve upward with quad interpolation            
              call sqcSGevnf(itype,
     +                   idf,idg,iyg,iord,nfld(i),izd2(i),izd1(i))
            endif 
          enddo
        endif

C--   End of loop over downward evolutions        
      enddo
      
      epf0 = dqcGetEps(itype,idf,iymax,it0)
      epf1 = dqcGetEps(itype,idf,iymax,it1)
      epf2 = dqcGetEps(itype,idf,iymax,it2)
      epg0 = dqcGetEps(itype,idg,iymax,it0)
      epg1 = dqcGetEps(itype,idg,iymax,it1)
      epg2 = dqcGetEps(itype,idg,iymax,it2)
      eps  = max(epf0,epf1,epf2,epg0,epg1,epg2)

      return
      end

C     =======================================================
      subroutine sqcSGevnf(itype,idf,idg,iyg,iord,nf,iz1,iz2)
C     =======================================================

C--   Singlet-gluon evolution from iz1 to iz2 at fixed nf
C--
C--   Input:  idf  = index of singlet pdf to be evolved
C--           idg  = index of gluon pdf
C--           iyg  = subgrid index 1,...,nyg2         
C--           iord = 1=LO,  2=NLO, 3=NNLO
C--           nf   = number of flavours
C--           iz1  = start of evolution
C--           iz2  = end of evolution
C--
C--   Output:  Spline coefficients in /qstor7/stor7 (linear store)
C--
C--   NB: the routines sqcSGmult and sqcSGeqs can be found in qcdutil.f

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
*      include 'qwtab7.inc'
      include 'qpdfs7.inc'
      include 'qmaps8.inc'

      dimension idwt(4)
C--                 QQ  QG  GQ  GG
      data idwt /    1,  2,  3,  4  /

      dimension sbar(mxx0,4),ssum(mxx0,4),vmat(mxx0,4)
      dimension bf(mxx0),bg(mxx0),hf(mxx0),hg(mxx0)

C--   Initialization
C--   --------------
C--   Interpolation index
      idk     = ioy2-1
C--   Weight table offset      
      ifirst7 = ifst7(itype,ioy2)
      
C--   Direction of evolution (isign) and first point after it1 (next)
      isign = 1
      next  = iz1+1
      if(iz2.lt.iz1) then 
        isign = -1
        next  = iz1-1
      endif
C--   Setup sbar
      do j = 1,4
        do i = 1,nyy2(iyg)
          sbar(i,j) = 0.D0
          ssum(i,j) = 0.D0
        enddo
      enddo
      
C--   Set here yloop index range      
      iy1 = 1
*      iy2 = nyy2(iyg)
      iy2 = iqcIyMaxG(iymac2,iyg)

C--   Calculate vector b at input scale iz1
C--   -------------------------------------
C--   Grid spacing delta = z(next)-z (divided by 2)
      delt = 0.5*abs(zgrid2(next)-zgrid2(iz1))
C--   Transformation matrix divided by delta (Sbar)        
      do i = 1,nmaty2
        sbar(i,1) = smaty2(i)/delt
        sbar(i,4) = smaty2(i)/delt
      enddo
C--   Weight matrix and  V = Sbar+W at t1
      do ityp = 1,4
        do iy = 1,iy2
          vmat(iy,ityp) = sbar(iy,ityp)
        enddo
C--     Find t-index
        it1 = itfiz2(iz1)        
        do k = 1,iord
          id = idPij7(idwt(ityp),k,itype)
          as = antab8(iz1,k)
          ia = iqcPijklm(1,it1,nf,iyg,id)-1
          do iy = 1,iy2
            ia = ia+1
            vmat(iy,ityp) = vmat(iy,ityp) + isign*as*stor7(ia)
          enddo
        enddo  
      enddo
C--   Address of pdf(iy=1,iz1,ipdf) in the store
      iaf = iqcPdfIjkl(1,iz1,idf,itype)
      iag = iqcPdfIjkl(1,iz1,idg,itype)
C--   Calculate V a = b
      call sqcSGmult(vmat(1,1),vmat(1,2),vmat(1,3),vmat(1,4),iy2,
     +               stor7(iaf),stor7(iag),bf,bg,iy2)

C--   Evolution loop over t
C--   ---------------------
      do iz = next,iz2,isign      
C--     Find t-index
        it = itfiz2(iz)
*mbC--     Yes/no apply roots cut, thats the question
*mb        if(isign.eq.1) then
*mb          iy2 = iqcIyMaxG(ismac2(iz),iyg)  !apply roots cut 4upward
*mb        else
*mb          iy2 = iqcIyMaxG(iymac2,iyg)  !not apply roots cut 4dnward
*mb        endif         
C--     Weight matrix and V matrix at t
        do ityp = 1,4
          do iy = 1,iy2
            vmat(iy,ityp) = sbar(iy,ityp)
          enddo
          do k = 1,iord
            id = idPij7(idwt(ityp),k,itype)
            as = antab8(iz,k)
            ia = iqcPijklm(1,it,nf,iyg,id)-1
            do iy = 1,iy2
              ia = ia+1
              vmat(iy,ityp) = vmat(iy,ityp) - isign*as*stor7(ia)
            enddo  
          enddo  
        enddo
C--     Address of pdf(1,iz,ipdf) in the store
        iaf = iqcPdfIjkl(1,iz,idf,itype)
        iag = iqcPdfIjkl(1,iz,idg,itype)
C--     Solve V a = b
        call sqcSGeqs(vmat(1,1),vmat(1,2),vmat(1,3),vmat(1,4),
     +                stor7(iaf),stor7(iag),bf,bg,iy2)
C--     Update b for the next iteration: not at last iteration of the loop
        if(iz.ne.iz2) then
C--       Grid spacing delta = z(next)-z (divided by 2)
          delt = 0.5*abs(zgrid2(iz+isign)-zgrid2(iz))
C--       Sum of current and next sbar; store next sbar       
          do i = 1,nmaty2
            sbnext    = smaty2(i)/delt
            ssum(i,1) = sbar(i,1)+sbnext
            ssum(i,4) = sbar(i,4)+sbnext
            sbar(i,1) = sbnext
            sbar(i,4) = sbnext
          enddo
C--       Calculate Ssum a = h
          call sqcSGmult(ssum(1,1),ssum(1,2),ssum(1,3),ssum(1,4),
     +                   nmaty2,stor7(iaf),stor7(iag),hf,hg,iy2)
C--       Update b
          do iy = 1,iy2
            bf(iy) = hf(iy)-bf(iy)
            bg(iy) = hg(iy)-bg(iy)
          enddo
        endif
      enddo

C--   Thats it...

      return
      end

C     =========================================================
      subroutine sqcSGjup(itype,ids,idg,iyg,iord,ny,izin,izout)
C     =========================================================

C--   Calculate NNLO singlet and gluon jump for upward evolution 
C--
C--   Input:  ids    = index of singlet pdf
C--           idg    = index of gluon pdf
C--           iyg    = subgrid index 1,...,nyg2         
C--           iord   = 1=LO,  2=NLO, 3=NNLO
C--           ny     = upper index yloop
C--           izin   = z-grid index where A(nf) is to be found 
C--           izout  = z-grid index where A(nf+1) will be stored
C--
C--   Output: Spline coefficients in /qstor7/stor7 (linear store)
C--
C--   Remark: Acts as a do-nothing when iord < 3 or when ids = idg

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
*      include 'qwtab7.inc'
      include 'qpdfs7.inc'
      include 'qmaps8.inc'

      data idgq, idgg, idqq, idhq, idhg / 1, 2, 3, 4, 5 /

      dimension wmatqq(mxx0), wmatqg(mxx0), wmatgq(mxx0), wmatgg(mxx0)
      dimension fjump(mxx0) , gjump(mxx0)
      
C--   Do nothing except copy izin to izout
      if(iord.lt.3 .or. ids.eq.idg) then
        call sqcPCopjj(itype,ids,izin,ids,izout)
        call sqcPCopjj(itype,idg,izin,idg,izout)
        return
      endif

C--   Weight table offset
      ifirst7 = ifst7(itype,ioy2)

C--   Alphas to use
*      assq = antab8(izout,2)
      assq = antab8(izout,0)*antab8(izout,0)

C--   Find t-index
      itin = itfiz2(izin)      

C--   Weight matrix
      nf = 3  !dummy variable since the A-weights do not depend on nf
      iaqq = iqcPijklm(1,itin,nf,iyg,idAij7(idqq))-1
      iahq = iqcPijklm(1,itin,nf,iyg,idAij7(idhq))-1
      iahg = iqcPijklm(1,itin,nf,iyg,idAij7(idhg))-1
      iagq = iqcPijklm(1,itin,nf,iyg,idAij7(idgq))-1
      iagg = iqcPijklm(1,itin,nf,iyg,idAij7(idgg))-1
      do iy = 1,ny
        wmatqq(iy) = assq * (stor7(iaqq+iy) + stor7(iahq+iy))
        wmatqg(iy) = assq *  stor7(iahg+iy)
        wmatgq(iy) = assq *  stor7(iagq+iy)
        wmatgg(iy) = assq *  stor7(iagg+iy)
      enddo
C--   Add transformation matrix to Wqq and Wgg
      do i = 1,nmaty2
        wmatqq(i) = wmatqq(i)+smaty2(i)
        wmatgg(i) = wmatgg(i)+smaty2(i)
      enddo
C--   Address of a(iy=1,izin,ipdf) in the store
      ias   = iqcPdfIjkl(1,izin,ids,itype) !address of input quark
      iag   = iqcPdfIjkl(1,izin,idg,itype) !address of input gluon
C--   Calculate fnew = W*a and store in buffers
      call sqcSGmult(wmatqq,wmatqg,wmatgq,wmatgg,ny,
     &               stor7(ias),stor7(iag),fjump,gjump,ny)
C--   Calculate anew by solving S*anew = fnew
      ias   = iqcPdfIjkl(1,izout,ids,itype) !address of output singlet
      iag   = iqcPdfIjkl(1,izout,idg,itype) !address of output gluon
      call sqcNSeqs(smaty2,nmaty2,stor7(ias),fjump,ny)
      call sqcNSeqs(smaty2,nmaty2,stor7(iag),gjump,ny)

      return
      end

C     =========================================================
      subroutine sqcSGjdn(itype,ids,idg,iyg,iord,ny,izin,izout)
C     =========================================================

C--   Calculate NNLO singlet and gluon jump for downward evolution 
C--
C--   Input:  ids    = index of singlet pdf
C--           idg    = index of gluon pdf
C--           iyg    = subgrid index 1,...,nyg2         
C--           iord   = 1=LO,  2=NLO, 3=NNLO
C--           ny     = upper index yloop
C--           izin   = z-grid index where A(nf+1) is to be found 
C--           izout  = z-grid index where A(nf) will be stored
C--
C--   Output: Spline coefficients in /qstor7/stor7 (linear store)
C--
C--   Remark: Acts as a do-nothing when iord < 3 or when ids = idg

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
*      include 'qwtab7.inc'
      include 'qpdfs7.inc'
      include 'qmaps8.inc'

      data idgq, idgg, idqq, idhq, idhg / 1, 2, 3, 4, 5 /

      dimension wmatqq(mxx0), wmatqg(mxx0), wmatgq(mxx0), wmatgg(mxx0)
      dimension fstor(mxx0) , gstor(mxx0) 

C--   Do nothing except copy izin to izout
      if(iord.lt.3 .or. ids.eq.idg) then
        call sqcPCopjj(itype,ids,izin,ids,izout)
        call sqcPCopjj(itype,idg,izin,idg,izout)
        return
      endif

C--   Weight table offset      
      ifirst7 = ifst7(itype,ioy2)

C--   Alphas to use
*      assq = antab8(izin,2)
      assq = antab8(izin,0)*antab8(izin,0)
      
C--   Find t-index
      itin = itfiz2(izin)      

C--   Weight matrix
      nf = 3  !dummy variable since the A-weights do not depend on nf
      iaqq = iqcPijklm(1,itin,nf,iyg,idAij7(idqq))-1
      iahq = iqcPijklm(1,itin,nf,iyg,idAij7(idhq))-1
      iahg = iqcPijklm(1,itin,nf,iyg,idAij7(idhg))-1
      iagq = iqcPijklm(1,itin,nf,iyg,idAij7(idgq))-1
      iagg = iqcPijklm(1,itin,nf,iyg,idAij7(idgg))-1
      do iy = 1,ny
        wmatqq(iy) = assq * (stor7(iaqq+iy) + stor7(iahq+iy))
        wmatqg(iy) = assq *  stor7(iahg+iy)
        wmatgq(iy) = assq *  stor7(iagq+iy)
        wmatgg(iy) = assq *  stor7(iagg+iy)
      enddo
C--   Add transformation matrix to Wqq and Wgg
      do i = 1,nmaty2
        wmatqq(i) = wmatqq(i)+smaty2(i)
        wmatgg(i) = wmatgg(i)+smaty2(i)
      enddo
C--   Address of a(iy=1,izin,ipdf) in the store
      ias   = iqcPdfIjkl(1,izin,ids,itype) !address of input quark
      iag   = iqcPdfIjkl(1,izin,idg,itype) !address of input gluon
C--   Calculate fold = S*a and store in buffers
      call sqcNSmult(smaty2,nmaty2,stor7(ias),fstor,ny)
      call sqcNSmult(smaty2,nmaty2,stor7(iag),gstor,ny)
C--   Calculate anew by solving W*anew = fold 
      ias   = iqcPdfIjkl(1,izout,ids,itype) !address of output singlet
      iag   = iqcPdfIjkl(1,izout,idg,itype) !address of output gluon
      call sqcSGeqs(wmatqq,wmatqg,wmatgq,wmatgg,
     &              stor7(ias),stor7(iag),fstor,gstor,ny)

      return
      end

C     =====================================================
      subroutine sqcSGStoreStart(itype,ids,idg,iy1,iy2,iz0)
C     =====================================================

C--   Store  startvalues into a buffer
C--
C--   ids    Singlet table
C--   idg    Gluon table
C--   iy1    First yloop index
C--   iy2    Last  yloop index
C--   iz0    z-bin containing startvalue  

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'

      common /stbuf/ qtarg(mxx0), gtarg(mxx0),
     +               qlast(mxx0), glast(mxx0)

C--   Calculate base addresses
      ias  = iqcPdfIjkl(iy1,iz0,ids,itype)-1
      iag  = iqcPdfIjkl(iy1,iz0,idg,itype)-1
C--   Store startvalues
      do j = iy1,iy2
        ias = ias+1
        iag = iag+1
        qtarg(j) = stor7(ias)
        gtarg(j) = stor7(iag)
        qlast(j) = stor7(ias)
        glast(j) = stor7(iag)
      enddo

      return
      end

C     =======================================================
      subroutine sqcSGNewStart(itype,ids,idg,iy1,iy2,iz0,eps)
C     =======================================================

C--   Set new startvalues
C--
C--   ids    Singlet table
C--   idg    Gluon table
C--   iy1    First yloop index
C--   iy2    Last yloop index
C--   iz0    z-bin containing startvalue
C--   eps    (out) Max |new-original|  

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'

      common /stbuf/ qtarg(mxx0), gtarg(mxx0),
     +               qlast(mxx0), glast(mxx0)

C--   Calculate base addresses
      ias = iqcPdfIjkl(iy1,iz0,ids,itype)-1
      iag = iqcPdfIjkl(iy1,iz0,idg,itype)-1
C--   New startvalues
      eps = -999.D0
      do j = iy1,iy2
        ias  = ias+1
        iag  = iag+1
        difs = stor7(ias)-qtarg(j)
        difg = stor7(iag)-gtarg(j)
        eps  = max(eps, abs(difs))
        eps  = max(eps, abs(difg))
        stor7(ias) = qlast(j)-difs
        stor7(iag) = glast(j)-difg
        qlast(j)   = stor7(ias)
        glast(j)   = stor7(iag)
      enddo

      return
      end

C     =======================================================
      subroutine sqcSGRestoreStart(itype,ids,idg,iy1,iy2,iz0)
C     =======================================================

C--   Restore  startvalues
C--
C--   ids    Singlet table
C--   idg    Gluon table
C--   iy1    First yloop index
C--   iy2    Last yloop index
C--   iz0    z-bin containing startvalue  

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'

      common /stbuf/ qtarg(mxx0), gtarg(mxx0),
     +               qlast(mxx0), glast(mxx0)

C--   Calculate base addresses
      ias = iqcPdfIjkl(iy1,iz0,ids,itype)-1
      iag = iqcPdfIjkl(iy1,iz0,idg,itype)-1
C--   Restore startvalues
      do j = iy1,iy2
        ias = ias+1
        iag = iag+1
        stor7(ias) = qtarg(j)
        stor7(iag) = gtarg(j)
      enddo

      return
      end

C     ===================================================
      double precision function dqcGetEps(itype,id,ny,it)
C     ===================================================

C--   Get max deviation quad-lin at midpoints
C--   This works on a subgrid with A coefficients
C--
C--   itype   (in)  : Pdf type identifier
C--   id      (in)  : Pdf identifier
C--   ny      (in)  : upper loop index in y
C--   it      (in)  : t-grid point

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      
      dimension epsi(mxx0)
      
      dqcGetEps = 0.D0
      if(ioy2.ne.3) return                

C--   Base address
      iz = izfit2(it)
      ia = iqcPdfIjkl(1,iz,id,itype)
C--   Now get vector of deviations        
      call sqcDHalf(ioy2,stor7(ia),epsi,ny)
C--   Max deviation
      do iy = 1,ny
        dqcGetEps = max(dqcGetEps,abs(epsi(iy)))
      enddo
      
      return
      end  

C     ================================================================
C     Setup pointer definitions
C     ================================================================

C     ========================
      subroutine sqcNfTab(it0)
C     ========================

C--   it0 = 0: no debug printout
C--   it0 # 0: generate debug printout. it0 is the starting point of a fake
C--            evolution in this printout
C-- 
C--   The lists created by this routine are explained by the following
C--   example of 3 flavors on an 8-point t-grid: nf = 3 from t1-t3,
C--   nf = 4 from t3-t6 and nf = 5 from t6-t8 (note the one-point overlap)
C--   So we have three flavors (3,4,5) which are put in a list (nflist8)
C--   The ranges of these flavors are also put in lists (itmin8, itmax8)
C--      
C--   it --->    1  2  3  4  5  6  7  8      nlist8 = 3 --> 1  2  3  4
C--    t        t1 t2 t3 t4 t5 t6 t7 t8     nflist8         3  4  5  0
C--   nfmap8(3)  1  1  1  0  0  0  0  0      itmin8         1  3  6  0
C--   nfmap8(4)  0  0  1  1  1  1  0  0      itmax8         3  6  8  0
C--   nfmap8(5)  0  0  0  0  0  1  1  1      
C--   nfmap8(6)  0  0  0  0  0  0  0  0      
C--
C--   Zgrid is a working grid with each threshold doubled
C--
C--   iz --->    1  2  3  4  5  6  7  8  9 10       i --->  1  2  3  4
C--   nf         3  3  3  4  4  4  4  5  5  5      nflist8  3  4  5  0
C--   zgrid     t1 t2 t3 t3 t4 t5 t6 t6 t7 t8       izmin8  1  4  8  0
C--                                                 izmax8  3  7 10  0
C--
C--                                                nf --->  3  4  5  6
C--                                                izminf8  1  4  8  0
C--                                                izmaxf8  3  7 10  0
C--   
C--   The following pointer arrays are used to map the indices:
C--
C--   i --->     1  2  3  4  5  6  7  8  9 10
C--   itfiz2     1  2  3  3  4  5  6  6  7  8             (it from iz)
C--   izfit2     1  2  4  5  6  8  9 10                   (iz from it)
C--   nffiz2     3  3  3  4  4  4  4  5  5  5             (nf from iz) 

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qmaps8.inc'

      dimension mm(mqq0),itmin(3:6)

C--   Initialize
      do i = 1,4
        nf = i+2
        do it = 0,mqq0+1
          nfmap8(it,nf) = 0
        enddo
        nflist8(i) = 0
        itmin8(i)  = 0
        itmax8(i)  = 0
        izmin8(i)  = 0
        izmax8(i)  = 0
      enddo
      do i = 3,6
        izminf8(i) = 0
        izmaxf8(i) = 0
      enddo  

C--   Set grid thresholds at the first gridpoint below the real thresholds
      tmax = tgrid2(ntt2)
      tchm = tthrs6(4)
      ichm = iqcItfrmT(tchm)
C--   iqcItfrmT returns zero if tchm above endpoint: we dont want this
C--   we also dont want a flavor change at the end point
      if(tchm.ge.tmax) ichm = ntt2+1
      tbot = tthrs6(5)
      ibot = iqcItfrmT(tbot)
      if(tbot.ge.tmax) ibot = ntt2+1
      ttop = tthrs6(6)
      itop = iqcItfrmT(ttop)
      if(ttop.ge.tmax) itop = ntt2+1
C--   Store for later use
      itchm2 = ichm
      itbot2 = ibot
      ittop2 = itop
    
C--   Fill nfmap
      do it = 1,ntt2
        if(it.le.ichm)                then
          nfmap8(it,3) = 1
        endif
        if(it.ge.ichm.and.it.le.ibot) then
          nfmap8(it,4) = 1
        endif
        if(it.ge.ibot.and.it.le.itop) then
          nfmap8(it,5) = 1
        endif
        if(it.ge.itop)                then
          nfmap8(it,6) = 1
        endif
      enddo

C--   Multiplicity = 3 at boundaries and flavor thesholds, m = 1 otherwise
      do i = 1,ntt2
        mm(i) = 1
      enddo
      mm(1) = 3
      if(ichm.ge.1.and.ichm.le.ntt2) mm(ichm) = 3
      if(ibot.ge.1.and.ibot.le.ntt2) mm(ibot) = 3
      if(itop.ge.1.and.itop.le.ntt2) mm(itop) = 3
      mm(ntt2) = 3

C--   Fill pointers in a fake evolution flavor by flavor upwards from it0 = 1
      iw     = 0
      nlist8 = 0
      do nf = 3,6
        itmin(nf) = 1
      enddo
      do it = 1,ntt2
        do nf = 3,6
          if(nfmap8(it-1,nf).eq.0 .and. nfmap8(it,nf).eq.1) then
C--         We hit the starting point of the evolution with nf flavors
            itmin(nf) = it
          endif
          if(nfmap8(it,nf).eq.1 .and. nfmap8(it+1,nf).eq.0) then
C--         We hit the end point of the evolution with nf flavors but
C--         evolve only when the end point is not equal to the start point
            if(it.gt.itmin(nf)) then
              nlist8          = nlist8+1
              nflist8(nlist8) = nf
              itmin8(nlist8)  = itmin(nf)
              itmax8(nlist8)  = it
C--           1st point of evolution
              it1             = itmin(nf)
C--           End of fake evolution loop
            endif
          endif 
        enddo
      enddo
      
C--   Fill pointers z grid
      iz     = 0
      nl     = 0
      do nf = 3,6
        do it = 1,ntt2
          if(nfmap8(it,nf).eq.1) then
            iz = iz+1
            if(nfmap8(it-1,nf).eq.0 .and. it.ne.ntt2) then
C--           First point with nf flavors
              nl          = nl+1 
              izmin8(nl)  = iz
              izminf8(nf) = iz
            endif
            if(nfmap8(it+1,nf).eq.0) then
C--           Last point with nf flavors
              izmax8(nl)  = iz
              izmaxf8(nf) = iz
            endif
            zgrid2(iz) = tgrid2(it)
            itfiz2(iz) = it
            izfit2(it) = iz
            nffiz2(iz) = nf               
          endif
        enddo
      enddo
C--   Store total number of z grid points
      nzz2 = iz

C--   Calculate B-spline basis for t-interpolation
      kk = 3    !k = 3 --> quadratic interpolation in t
      call sqcSpqIni(kk,tgrid2,mm,ntt2,npars8,nc)
      
C--   Now put cuts on the kinematic plane      
      call sqcEvCuts(it0)

C--   Done...
      Lnfmap8 = .true.
      if(it0.eq.0) return

C--   Debug printout
      write(6,*) 'tchm,tbot,ttop',tchm,tbot,ttop
      write(6,*) 'ichm,ibot,itop',ichm,ibot,itop
      write(6,*) 'nflavors,npars',nlist8,npars8
C--   Print nfmap
      write(6,'(/'' it             tt  3  4  5  6  m'')')
      do i = 1,ntt2
        write(6,'(I3,E15.5,5I3)') i,tgrid2(i),(nfmap8(i,j),j=3,6),
     &         mm(i)
      enddo
C--   Print ranges
      write(6,'(/''nflist8 '',4I4)') nflist8
      write(6,'( ''itmin8  '',4I4)') itmin8
      write(6,'( ''itmax8  '',4I4)') itmax8
      write(6,'( ''izmin8  '',4I4)') izmin8
      write(6,'( ''izmax8  '',4I4)') izmax8
      write(6,'(/''nf--->     3   4   5   6'')')
      write(6,'( ''izminf8 '',4I4)') izminf8
      write(6,'( ''izmaxf8 '',4I4)') izmaxf8  
C--   Pointer arrays
      write(6,'(/''i ---> '',23I3)') (i,i=1,23)
      write(6,'( ''izfit  '',23I3)') (izfit2(i),i=1,23)
      write(6,'( ''itfiz  '',23I3)') (itfiz2(i),i=1,23)
      write(6,'( ''nffiz  '',23I3)') (nffiz2(i),i=1,23)
C--   Fake evolution loop for testing
C--   Figure out which flavor range we are, take the highest flavor
      write(6,*) ' '
      write(6,*) 'it0 =',it0
      ibin0 = 0
      do i = 1,nlist8
        if(itmin8(i).le.it0 .and. it0.le.itmax8(i)) ibin0 = i
      enddo
      if(ibin0.eq.0) stop 'sqcNfTab: it0 out of range ---> STOP'
C--   Upward from it0, evolve flavor by flavor
      do i = ibin0,nlist8
        it1 = max(it0,itmin8(i))
        it2 = itmax8(i)
        nf  = nflist8(i)
        if(it1.lt.it2) write(6,*) 'evolup nf, it1, it2',nf,it1,it2
      enddo
C--   Downward from it0, evolve flavor by flavor
      do i = ibin0,1,-1
        it1 = min(it0,itmax8(i))
        it2 = itmin8(i)
        nf  = nflist8(i)
        if(it2.lt.it1) write(6,*) 'evoldn nf, it1, it2',nf,it1,it2
      enddo

      return
      end
      
      
C     ============================      
      subroutine sqcEvCuts(iprint)
C     ============================

C--   Set cuts on the kinematic plane
C--
C--   iprint  (in)   debug printout if .ne. 0
C--   output in common qgrid2.inc
C--      iymac2     = largest  y-index within cuts
C--      izmic2     = smallest z-index within cuts
C--      izmac2     = largest  x-index within cuts
C--      ismac2(iz) = -1   --> outside z-range [izmic2,izmac2]
C--                 =  iy  --> max y-index below roots cut  

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qmaps8.inc'
      
      character*80 string

*mb   Roots cut disabled      
*mb      dimension ylist(mqq0),tlist(mqq0)
*mb      dimension iy1(mqq0),iy2(mqq0),iz1(mqq0),iz2(mqq0)
*mb      logical   mark(0:mxx0,0:mqq0+7)
      
      margin = 0
      
C--   Roots cut set out-of range by default
      do iz = 1,nzz2
        ismac2(iz) = -1
      enddo
      
C--   Get ymax and zmin,zmax limits
      call sqcZmesh(ymaxc2,tminc2,margin,iymi,iymac2,izmic2,izma)
      call sqcZmesh(ymaxc2,tmaxc2,margin,iymi,iymac2,izmi,izmac2)
      
*mbC--   Roots cut 
*mb      basically works but disable since not worth it
*mb      nlist = 0
*mb      do iz = 1,nzz2
*mb        tt = zgrid2(iz)
*mb        qq           = exp(tt)
*mb        xx           = qq/smaxc2
*mb        xx           = max(xx,xminc2)
*mb        if(xx .le. xmaxc2) then
*mb          nlist        =  nlist+1
*mb          ylist(nlist) = -log(xx)
*mb          tlist(nlist) =  tt
*mb        endif  
*mb      enddo
*mb      if(nlist.eq.0) stop 'sqcEvCuts: empty y-z region after cuts'
*mb      
*mb      call sqcMarkit(mark,ylist,tlist,margin,iy1,iy2,iz1,iz2,nlist)
*mb      
*mb      izma = -1
*mb      do iz = 1,nzz2
*mb        iyma = -1
*mb        do iy = 1,iymac2
*mb          if(mark(iy,iz)) iyma = iy
*mb        enddo
*mb        ismac2(iz) = iyma 
*mb        if(iyma.ne.-1) izma = iz 
*mb      enddo
*mb      if(izma.eq.-1) stop 'sqcEvCuts: empty z region after cuts'
*mb      izmac2 = izma
      
C--   Done...      
      Levcut8 = .true.
      
      if(iprint.eq.0) return
      
C--   Debug printout
      write(6,'(/'' vertical iy, horizontal iz''/)')
      do iy = nyy2(0),1,-1
        do iz = 1,nzz2
          string(iz:iz) = '.'
          if(ismac2(iz).eq.iy) string(iz:iz) = '*'
        enddo
        write(6,'(I4,1X,A)') iy,string(1:nzz2)
      enddo
      write(6,'(/)')    

      return
      end
      
C     ================================      
      subroutine sqcTimeFac(facL,facQ)
C     ================================

C--   facL   (out)  linear time gain factor
C--   facQ   (out)  quadratic time gain factor

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qgrid2.inc'
      
      totL = nyy2(0)*nzz2
      totQ = nyy2(0)*totL
      
      cutL = iymac2*(izmac2-izmic2+1)
      cutQ = iymac2*cutL
      
*mb      cutL = 0.D0
*mb      cutQ = 0.D0
      
*mb      do iz = izmic2,izmac2
*mb        if(ismac2(iz).eq.-1) stop 'sqcTimeFac: iz out of range'
*mb        cutL = cutL + ismac2(iz)
*mb        cutQ = cutQ + ismac2(iz)*ismac2(iz)
*mb      enddo
      
      facL = cutL/totL
      facQ = cutQ/totQ
      
      return
      end              
      
C     =====================================      
      integer function iqcIyMaxG(iymax0,ig)
C     =====================================

C--   Translate y-max index in G0 to an y-max index in Gi
C--
C--   iymax0     (in)    [1,nyy2(0)]    y-max index in G0 
C--   ig         (in)    [1,nyg2]       subgrid index 
C--   iqcIyMaxG  (out)   [1,nyy2(ig)]   y-max index in subgrid

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qgrid2.inc'
      
      ymax = ygrid2(iymax0)
      iy   = iqcIyfrmY(ymax,dely2(ig),nyy2(ig)) !find index below ymax  
      if(iy.eq.-1) then 
        iy = nyy2(ig)               !index out of subgrid range
      else
        iy = min(iy+1,nyy2(ig))     !put index above ymax
      endif
      
      iqcIyMaxG = iy
      
      return
      end        
      
C     =============================================================
      subroutine sqcEvLims(it0,it1,it2,
     +               izu1,izu2,nflu,nup,izd1,izd2,nfld,ndn,ibl,ibu)
C     =============================================================

C--   Setup ranges for up and downward evolution 
C--
C--   it0         (in)   tgrid start point
C--   it1         (in)   lower tgrid limit of the evolution
C--   it2         (in)   upper tgrid limit of the evolution
C--   izu1(4)     (out)  lower zgrid limits of upward evolutions
C--   izu2(4)     (out)  upper zgrid limits of upward evolutions
C--   nflu(4)     (out)  number of flavors of each upward evolution
C--   nup         (out)  number of upward evolutions
C--   izd1(4)     (out)  upper zgrid limits of downward evolutions
C--   izd2(4)     (out)  lower zgrid limits of downward evolutions
C--   nfld(4)     (out)  number of flavors of each downward evolution
C--   ndn         (out)  number of downward evolutions
C--   ibl         (out)  lowest nf bin accessed
C--   ibu         (out)  highest nf bin accessed

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qmaps8.inc'

      dimension izu1(4),izu2(4),nflu(4),izd1(4),izd2(4),nfld(4)

C--   Setup a flavor map for up and downward evolution
      if(.not.Lnfmap8) call sqcNfTab(0)
      if(.not.Levcut8) call sqcEvCuts(0)

C--   Initialize
      nup = 0
      ndn = 0
      idl = 0
      idu = 0
      do i = 1,4
        izu1(i) = 0
        izu2(i) = 0
        nflu(i) = 0
        izd1(i) = 0
        izd2(i) = 0
        nfld(i) = 0
      enddo
C--   Adjust upper and lower evolution limits to grid limits
      jt1    = max(it1,1)
      jt2    = min(it2,ntt2)
C--   Check it0 is in between grid limits
      if(it0.lt.jt1 .or. it0.gt.jt2) return !it0 out of range 

C--   Find out in which region it0, jt1 and jt2 falls
      ibin0 = 0
      ibin1 = 0
      ibin2 = 0
      do i = 1,nlist8
        if(itmin8(i).le.it0 .and. it0.le.itmax8(i)) ibin0 = i
        if(itmin8(i).le.jt1 .and. jt1.le.itmax8(i)) ibin1 = i
        if(itmin8(i).le.jt2 .and. jt2.le.itmax8(i)) ibin2 = i
      enddo
      if(ibin0.eq.0) return   !it0 out of range
      if(ibin1.eq.0) return   !it1 out of range
      if(ibin2.eq.0) return   !it2 out of range
C--   Calculate iz0; izmin-itmin gives the offset between iz and it
      iz0 = it0 + izmin8(ibin0)-itmin8(ibin0)
      izd = it1 + izmin8(ibin1)-itmin8(ibin1)
      izu = it2 + izmin8(ibin2)-itmin8(ibin2)
C--   Upward from iz0, evolve flavor by flavor
      do i = ibin0,ibin2
        iz1 = max(iz0,izmin8(i))
        iz2 = min(izmax8(i),izu)
        nf  = nflist8(i)
        if(iz1.lt.iz2) then
          nup       = nup+1
          izu1(nup) = iz1
          izu2(nup) = iz2
          nflu(nup) = nf
        endif
      enddo
C--   Downward from iz0, evolve flavor by flavor
      do i = ibin0,ibin1,-1
        iz1 = min(iz0,izmax8(i))
        iz2 = max(izmin8(i),izd)
        nf  = nflist8(i)
        if(iz2.lt.iz1) then
          ndn       = ndn+1
          izd1(ndn) = iz1
          izd2(ndn) = iz2
          nfld(ndn) = nf
        endif
      enddo
C--   Range of regions accessed
      ibl = ibin1
      ibu = ibin2

      return
      end
      
      
      
      