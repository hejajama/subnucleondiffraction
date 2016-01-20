
C--   This is the file qcdalf.f containing the qcdnum alpha_s routines

C--   subroutine                sqcAlfTab(iord)
C--   double precision function dqcAsEvol(rs1,rs0,as0,rmsq,iord,nff,ierr)
C--   double precision function dqcA0ToA1(rs1,rs0,as0,rmsq,iord,nff,ierr)
C--   double precision function dqcAlfNew(alfnf,r2,f2,iord,nfjump,asjump)
C--   double precision function dqcGetAm(asmin)
C--   double precision function dqcAlfar(r2,r20,alp0,nf,io,ierr)
C--   double precision function dqcAjump(alfnf,r2,f2,iord)
C--   integer function          iqcGetNf(scalesq,thrs,ihit)
C--   subroutine                sqcGetLim(r1,r2,thrs,n,rr1,rr2,nfl,iupd)

C===================================================================
C==   Alpha_s evolution routines ===================================
C===================================================================

C     ==========================
      subroutine sqcAlfTab(iord)
C     ==========================

C--   Make lookup table of as/2pi, (as/2pi)^2 and (as/2pi)^3.
C--
C--   Input:           iord        = order of evolution.
C--           /qpars6/ q0alf6      = r20
C--                    alfq06      = value of alphas at r20
C--                    aar6, bbr6  = transform f2 to r2
C--                    qthrs6(4:6) = thresholds on f2
C--                    rthrs6(4:6) = thresholds on r2
C--                    nfix6       = 0 when in vfns
C--           /qgrid2/ tgrid2      = table of ln(Q2) values.
C--
C--   Output: /qmaps8/ astab8      = table of alphas values.
C--                    itlow8      = lowest gridpoint with valid alphas
C--                    itmin6      = idem
C--
C--                    ierr    # 0 --> two or more thresholds map onto the
C--                              same gridpoint --> fatal error (no table)
C--
C--   sqcNfTab should have been called before

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qmaps8.inc'

      data pi /3.14159265358979D0/

C--   Make sure sqcNfTab is called
      if(.not.Lnfmap8) call sqcNfTab(0) 

C--   Update the thresholds on the renormalization scale (not for MVNS)
      if(.not.Lmfns6) call sqcRmass2(qthrs6,rthrs6)

C--   Cross check
      if(iord.ne.iord6) stop 
     +            'sqcAlfTab: inconsistent QCD order ---> STOP'

      itlow8 = 1
      itmin6 = 1

C--   Table in new format

      do iz = 1,nzz2

        it    = itfiz2(iz)
        qf2   = exp(tgrid2(it))
        qr2   = aar6*qf2+bbr6
        jerr  = 1
        as    = qnull6
        nfiz  = nffiz2(iz)
        as    = dqcAsEvol(qr2,q0alf6,alfq06,rthrs6,iord,nf,jerr)        

        if(jerr.ne.0 .or. as.gt.aslim6) then  
C--       Scale is too low or alphas too large
          itlow8        = it+1
          itmin6        = it+1
          antab8(iz, 0) = qnull6
          antab8(iz, 1) = qnull6
          antab8(iz, 2) = qnull6
          antab8(iz, 3) = qnull6
          antab8(iz,-1) = qnull6
          antab8(iz,-2) = qnull6
          antab8(iz,-3) = qnull6
        else
C--       If nf(iz+1) = nf(iz)+1 we are at a threshold with (3,4,5) convention
          if(iz.ne.nzz2 .and. (nffiz2(iz+1).eq.nfiz+1)) then
            as    = dqcAsEvol(-qr2,q0alf6,alfq06,rthrs6,iord,nf,jerr)
          endif
C--       Check
          if(.not.Lmfns6 .and. nf.ne.nfiz) 
     +                    stop 'sqcAlfTab: problem with nf'          
C--       Qcdnum uses alphas/2pi
          as  = as/(2.D0*pi)
          as2 = as*as
          as3 = as2*as
C--       Calculate beta functions for nf flavors
          b0 = 11.D0/2.D0-nf/3.D0
          b1 = 51.D0/2.D0-nf*19.D0/6.D0
          b2 = 2857.D0/16.D0-nf*5033.D0/144.D0+nf*nf*325.D0/432.D0
C--       Now Taylor expand and truncate to appropriate order
C--       See Section 2.3 in the QCDNUM manual
          fac = log(qf2/qr2)
          if(iord.eq.1) then
            antab8(iz, 0) = as
            antab8(iz, 1) = as
            antab8(iz, 2) = 0.D0
            antab8(iz, 3) = 0.D0
            antab8(iz,-1) = 0.D0
            antab8(iz,-2) = 0.D0
            antab8(iz,-3) = 0.D0
          elseif(iord.eq.2) then
            antab8(iz, 0) = as
            antab8(iz, 1) = as-b0*fac*as2
            antab8(iz, 2) = as2
            antab8(iz, 3) = 0.D0
            antab8(iz,-1) = as
            antab8(iz,-2) = 0.D0
            antab8(iz,-3) = 0.D0
          elseif(iord.eq.3) then
            antab8(iz, 0) = as
            antab8(iz, 1) = as-b0*fac*as2-(b1*fac-b0*b0*fac*fac)*as3
            antab8(iz, 2) = as2-2.D0*b0*fac*as3
            antab8(iz, 3) = as3
            antab8(iz,-1) = as-b0*fac*as2
            antab8(iz,-2) = as2
            antab8(iz,-3) = 0.D0
          else
            stop 'sqcAlfTab: unknown order (iord)'
          endif
        endif

      enddo

C--   Set status flag that alphas table is available
      Lastab8 = .true.

      return
      end

C     ===============================================================
      double precision function
     +                      dqcAsEvol(rs1,rs0,as0,rmsq,iord,nff,ierr)
C     ===============================================================

C--   Interface routine to dqcA0ToA1 allowing for different conventions
C--   of flavor settings at the thresholds. The conventions are:
C--       QCDNUM   nf = (4,5,6) at thresholds (c,b,t)
C--       PEGASUS  nf = (3,4,5) at thresholds (c,b,t)
C--   If the input (output) scale is prepended by a minus sign we
C--   take the PEGASUS convention otherwise the QCDNUM convention.
C--
C--   Input:   rs1         =   ren scale^2 of output alpha_s.
C--            rs0         =   ren scale^2 of input  alpha_s.
C--            as0         =   input alpha_s.
C--            rmsq(4:6)   =   thresholds defined on renormalization scale.
C--            iord        =   1,2,3 = LO,NLO,NNLO.
C--
C--   Output:  nff         =   # flavors at rs1.
C--            ierr        =   0, all OK.
C--                            1, rs1 too close to lambda2.
C--            dqcAsEvol   =   alphas at rs1 if ierr = 0.
C--                            some large value if ierr # 0.

      implicit double precision (a-h,o-z)
  
      include 'qcdnum.inc'
      include 'qpars6.inc'

      dimension rmsq(4:6)

C--   Catch too low scales
      dqcAsEvol = qnull6
      ierr      = 1
      if(abs(rs1).lt.0.1D0) return
      if(abs(rs0).lt.0.1D0) return      

C--   Initialize
      ierr = 0
      asin = as0

C--   Do we have (3,4,5) convention on input?
      if(rs0.lt.0.D0) then
C--     Yes, now check if rs0 is at a threshold
        nf0 = iqcGetNf(abs(rs0),rmsq,ihit)
        if(ihit.ne.0) then
C--       Hit a threshold; apply discontinuity nf --> nf+1
          nfjump = 1
          fs0    = (abs(rs0)-bbr6)/aar6
          asin   = dqcAlfNew(as0,abs(rs0),fs0,iord,nfjump,asjump)
        endif
      endif

C--   At this point we have the (4,5,6) convention 
      asout = dqcA0ToA1(abs(rs1),abs(rs0),asin,rmsq,iord,nff,ierr)
        
C--   Do we have (3,4,5) convention on output?
      if(rs1.lt.0.D0) then
C--     Yes, now check if rs1 is at a threshold
        nf1 = iqcGetNf(abs(rs1),rmsq,ihit)
        if(ihit.ne.0) then
C--       Hit a threshold; calculate discontinuity nf --> nf-1
          nfjump = -1
          fs1    = (abs(rs1)-bbr6)/aar6
          asout  = dqcAlfNew(asout,abs(rs1),fs1,iord,nfjump,asjump)
          nff    = nff-1
        endif
      endif

      dqcAsEvol = asout

      return
      end

C     ===============================================================
      double precision function
     +                      dqcA0ToA1(rs1,rs0,as0,rmsq,iord,nff,ierr)
C     ===============================================================

C--   Steering routine to carry alpha_s evolution over flavor thresholds.
C--
C--   Input:   rs1         =   ren scale^2 of output alpha_s.
C--            rs0         =   ren scale^2 of input  alpha_s.
C--            as0         =   input alpha_s.
C--            rmsq(4:6)   =   thresholds defined on renormalization scale.
C--            iord        =   1,2,3 = LO,NLO,NNLO.
C--
C--   Output:  nff         =   # flavors at rs1.
C--            ierr        =   0, all OK.
C--                            1, rs1 too close to lambda2.
C--            dqcA0ToA1   =   alphas at rs1 if ierr = 0.
C--                            some large value if ierr # 0.
C--
C--   Remark:  This routine assumes nf = (4,5,6) at the thresholds (c,b,t) 

      implicit double precision (a-h,o-z)
  
      include 'qcdnum.inc'
      include 'qpars6.inc'

      dimension rmsq(4:6),rr1(4),rr2(4),nfl(4)

C--   Get evolution limits
      call sqcGetLim(rs0,rs1,rmsq,nn,rr1,rr2,nfl,iupd)

C--   Nothing to evolve
      if(iupd.eq.0) then
        dqcA0ToA1 = as0
        nff       = iqcGetNf(rs0,rmsq,ihit)
        ierr      = 0
        return
      endif

C--   Upward evolution 
      if(iupd.eq.1) then
        asend = as0
        do i = 1,nn
          rend  = rr2(i)
          asend = dqcAlfaR(rend,rr1(i),asend,nfl(i),iord,jerr)
C--       Error
          if(jerr.ne.0) goto 200
C--       Is the endpoint at a flavor threshold?
          nfdum = iqcGetNf(rend,rmsq,ihit)
C--       Add discontinuity if the endpoint is at a flavor threshold
          if(ihit.ne.0) then
            nfjump = 1
            fend   = (rend-bbr6)/aar6
            asend  = dqcAlfNew(asend,rend,fend,iord,nfjump,asjump)
          endif
        enddo

C--   Downward evolution
      else
        asend = as0
        do i = 1,nn
C--       Is the starting point at a flavor threshold?
          rstart = rr1(i)
          nfdum  = iqcGetNf(rstart,rmsq,ihit)
C--       Subtract discontinuity if the start point is at a flavor threshold
          if(ihit.ne.0) then 
            nfjump = -1
            fstart = (rstart-bbr6)/aar6
            asend  = dqcAlfNew(asend,rstart,fstart,iord,nfjump,asjump)
          endif
C--       evolve
          asend = dqcAlfaR(rr2(i),rstart,asend,nfl(i),iord,jerr)
C--       error
          if(jerr.ne.0) goto 200
        enddo
      endif

C--   Evolution done w/o error
      dqcA0ToA1 = asend
      nff       = iqcGetNf(rs1,rmsq,ihit)
      ierr      = jerr
      return

C--   Attempt to evolve below Lambda
 200  continue
      dqcA0ToA1 = qnull6
      nff       = 0
      ierr      = jerr

      return 
      end

C     ============================================================
      double precision function 
     +                   dqcAlfNew(alfnf,r2,f2,iord,nfjump,asjump)
C     ============================================================

C--   Calculate at a flavor threshold as(nf+1) or as(nf-1) from as(nf)
C--
C--   Input:    alfnf      alphas for nf flavors
C--             r2         renormalization scale where alphas is given
C--             f2         factorization   scale where alphas is given
C--             iord       1 = LO, 2 = NLO, 3 = NNLO
C--             nfjump    -1 calculate as(nf-1) for downward evolution
C--                       +1 calculate as(nf+1) for upward   evolution
C--
C--   Output:   dqcAlfNew  new value of alphas
C--             asjump     as(nf+1) - as(nf)   or   as(nf) - as(nf-1)
  
      implicit double precision (a-h,o-z)
  
      external dqcGetAm

      common /apass/ asplus,r2scale,f2scale,iorder

C--   No jump at LO
      if(iord.eq.1) then
        dqcAlfNew = alfnf
        asjump    = 0.D0
        return
      endif
C--   Jump for upward evolution from nf to nf+1
      if(nfjump.eq. 1) then
        asnew  = alfnf + dqcAjump(alfnf,r2,f2,iord)
        asjump = asnew - alfnf
C--   Jump for downward evolution from nf to nf-1
      elseif(nfjump.eq.-1) then
C--     calculate as(nf-1) by solving the equation
C--     alf(nf) = alf(nf-1) + ajump[alf(nf-1)] 
C--     In present notation this reads
C--     alfnf   = asnew     + dqcAjump(asnew)
C--     Thus we have to solve 
C--     dqcGetAm(asnew) = alfnf - asnew - dqcAjump(asnew) = 0
C--     As a first step, pass values needed by dqcGetAm by /apass/
        asplus  = alfnf
        r2scale = r2
        f2scale = f2
        iorder  = iord
C--     Now make an initial guess of interval containing asnew
        ami = 0.95 * asplus
        ama = 1.05 * asplus
C--     Find the bracket by extending the initial guess, if needed
        call sqcBrackit(dqcGetAm,ami,ama,ierr)
        if(ierr.ne.0) stop 
     +    'dqcA0ToA1: cant bracket alfas downward evolution ---> STOP'
C--     Bracket found, now find asmin within eps = 10^-6
        asnew = dqcBiSect(dqcGetAm,ami,ama,1.D-6,ierr)
        if(ierr.ne.0) stop 
     +    'dqcA0ToA1: cant find as(nf-1) within tolerance ---> STOP'
        asjump = alfnf - asnew 
      else
        stop 'dqcAlfNew: invalid nfjump'
      endif

      dqcAlfNew = asnew

      return
      end

C     =========================================
      double precision function dqcGetAm(asmin)
C     =========================================

C--   Solve the alphas jump: asmin = asplus - Delta(asmin) or equivalently:
C--   dqcGetAm(asmin) = asplus-asmin-dqcAjump(asmin,r2,f2,iord) = 0.
C--
C--   input               asmin    alphas just before the nf threshold 
C--   input via /apass/   asplus   alphas just after  the nf threshold
C--                       r2scale  renormalization scale
C--                       f2scale  factorization   scale
C--                       iord     1 = LO, 2 = NLO, 3 = NNLO

      implicit double precision (a-h,o-z)

      common /apass/ asplus,r2scale,f2scale,iorder

      dqcGetAm = asplus-asmin-dqcAjump(asmin,r2scale,f2scale,iorder)

      return
      end

C     ==========================================================
      double precision function dqcAlfar(r2,r20,alp0,nf,io,ierr)
C     ==========================================================

C--   Code taken from Andreas Vogt 27-05-2001
C--
C--   Calculate alphas from Runge Kutta integration of RGE.
C--
C--   r2     (in)     renormalization scale (GeV2)
C--   r20    (in)     input scale
C--   alp0   (in)     alphas at input scale
C--   nf     (in)     [3-6] number of flavors
C--   io     (in)     [1-4] order (1=LO, 2=NLO, 3= NNLO, 4= N3LO)
C--   ierr  (out)     0 = all OK
C--                   1 = r2 too close to lambda2 (da/dlnr2 too large)
C--
C--   MB 12-08-02: Correct nasty bug: at the first call nf was redefined
C--                in the loop calculating beta ---> rename nf to mf in
C--                this loop.
C--   MB 02-09-10: Raise cut on da/dlnr2 to 10^4 (was 0.3)

      implicit double precision (a-h,o-z)

      save      beta
      dimension beta (3:6, 0:3)

      data nsteps  /10/
      data pi      /3.14159265358979D0/
      data dadrlim /1.0D4/

      dimension ic(0:3,1:4)
      data ic / 1, 0, 0, 0,
     +          1, 1, 0, 0,
     +          1, 1, 1, 0,
     +          1, 1, 1, 1  /

      logical first
      save    first
      data    first /.true./
      
*mb
*mb      logical lprint
*mb      save lprint,nflast
*mb      


      fbeta(a) = - a**2 * (beta(nf,0)*ic(0,io) + 
     +                a * (beta(nf,1)*ic(1,io) +
     +                a * (beta(nf,2)*ic(2,io) +
     +                a *  beta(nf,3)*ic(3,io)     )))
     
*mb      
*mb      if(first) then
*mb        lprint = .true.
*mb        nflast = nf
*mb      endif
*mb      
*mb      if(nf.ne.nflast .and. lprint) then
*mb        write(6,*) 'ALFAS: change number of flavours'
*mb        lprint = .false.
*mb      endif
*mb        

      if(first) then
        do k1 = 3,6
          mf  = k1
          mf2 = mf * mf
          mf3 = mf * mf2
          beta(mf,0) = 11.0000 - .666667* mf
          beta(mf,1) = 102.000 - 12.6667* mf
          beta(mf,2) = 1428.50 - 279.611* mf + 6.01852* mf2
          beta(mf,3) = 29243.0 - 6946.30* mf + 405.089* mf2 + 
     +                                         1.49931* mf3
        enddo
        first = .false.
      endif

      dlq   = dlog(r2/r20)/nsteps
      as    = alp0 / (4.* pi)
      ierr  = 0
      do k1 = 1, nsteps
        bet = fbeta(as)
        if(bet.lt.-dadrlim) then
          dqcAlfar = 0.D0
          ierr   = 1
          return
        endif
        xk0 = dlq * bet
        xk1 = dlq * fbeta(as+0.5*xk0)
        xk2 = dlq * fbeta(as+0.5*xk1)
        xk3 = dlq * fbeta(as+xk2)
        as  = as + 1./6.d0 * (xk0 + 2.* xk1 + 2.* xk2 + xk3) 
      enddo

      if(abs(fbeta(as)).gt.dadrlim) then
        dqcAlfar = 0.D0
        ierr   = 1
        return
      endif

C--   This routine calculated alphas/4pi but we return alphas
      dqcAlfar = as*4.*pi

      return
      end

C     ====================================================
      double precision function dqcAjump(alfnf,r2,f2,iord)
C     ====================================================

C--   Calculate the discontinuity in the alpha_s evolution at the
C--   renormalization scale r2.
C--   See eq. (2.41) in the Pegasus writeup hep-ph/0408244.
C--
C--       alf(r2,nf+1) = alf(r2,nf) + dqcAjump[alf(nf),r2,f2,iord]
C--
C--   Input:  alfnf  =  value of alphas for nf flavors
C--           r2     =  renormalization scale where alfnf is given
C--           f2     =  factorization scale corresponding to r2
C--         iord     =  1,2,3 for LO, NLO, NNLO
C--
C--   Output: dqcAjump = alf(nf+1)-alf(nf) at the scale r2   

      implicit double precision (a-h,o-z)
      logical first

      dimension c(2,0:2)
      save c

      data first /.true./
      data pi /3.14159265358979D0/

      if(first) then
        c(1,0) = 0.D0
        c(1,1) = 2.D0/3.D0
        c(2,0) = 14.D0/3.D0
        c(2,1) = 38.D0/3.D0
        c(2,2) = 4.D0/9.D0
        first  = .false.
      endif

C--   No jump at LO
      if(iord.le.1) then
        dqcAjump = 0.D0
        return
      endif

C--   The formulas apply to alphas/4pi
      alfin = alfnf/(4*pi)
      ratlg = log(r2/f2)
      ajump = 0.D0
      aterm = alfin
      do n = 1,iord-1
        cterm = 0.D0
        rterm = 1.D0
        do j = 0,n
          cterm = cterm + c(n,j)*rterm
          rterm = rterm*ratlg
        enddo
        aterm = aterm*alfin
        ajump = ajump + aterm*cterm
      enddo
C--   Multiply by 4pi to get the jump in alphas
      dqcAjump = ajump*4*pi

      return
      end

C     ============================================
      integer function iqcGetNf(scalesq,thrs,ihit)
C     ============================================

C--   Find the number of flavors at scalesq
C--
C--   Input:  scalesq           input scale (squared)
C--           thrs(4:6)         thresholds  (defined on that scale)
C--   Output: iqcGetNf          number of flavors at scalesq
C--           ihit              1/0 if scale is yes/no at the threshold
 
      implicit double precision (a-h,o-z)
      logical lqcRcomp

      include 'qcdnum.inc'
      include 'qpars6.inc'

      dimension thrs(4:6)

      nf   = 3
      ihit = 0
      do i = 4,6
        if(scalesq.ge.thrs(i)) nf = i
C--     do we hit a threshold?
        if(lqcRcomp(scalesq,thrs(i),repsi6)) then
          nf   = i
          ihit = 1
        endif
      enddo

      iqcGetNf = nf

      return
      end

C     ===================================================
      subroutine sqcGetLim(r1,r2,thrs,n,rr1,rr2,nfl,iupd)
C     ===================================================

C--   Setup evolution limits
C--
C--   Input:  r1                start scale
C--           r2                end   scale
C--           thrs(4:6)         thresholds  (defined on that scale)
C--   Output: n                 number of regions
C--           rr1(4)            rr1(i) i = 1,..n start point for region i
C--           rr2(4)            rr2(i) i = 1,..n end   point for region i
C--           nfl(4)            nfl(i) i = 1,..n # flavors   for region i
C--           iupd              -1 downward evolution
C--                              0 no evolution (r1 = r2)
C--                             +1 upward evolution
 
      implicit double precision (a-h,o-z)
      logical lqcRxeqy, lqcRxlty

      include 'qcdnum.inc'
      include 'qpars6.inc'

      dimension thrs(4:6),rr1(4),rr2(4),nfl(4)

C--   Initialize
      do i = 1,4
        rr1(i) = 0.D0
        rr2(i) = 0.D0
        nfl(i) = 0
      enddo

C--   No evolution, thank you
      if(lqcRxeqy(r1,r2,repsi6)) then
        iupd = 0
        return
      endif

C--   Evolution direction
      if(lqcRxlty(r1,r2,repsi6)) then
        iupd = 1
      else
        iupd = -1
      endif

C--   Number of flavors at the startpoint and the endpoint
      nf1 = iqcGetNf(r1,thrs,ihit1)
      nf2 = iqcGetNf(r2,thrs,ihit2)


      if(iupd.eq.1) then

C--     Upward evolution
C--     If the end point is at a flavor threshold nf2 --> nf2-1
        if(ihit2.ne.0) nf2 = nf2-1
        n = 0
        do i = nf1,nf2
          n = n+1
          if(i.eq.3) then
            rr1(n) = r1
            rr2(n) = min(r2,thrs(4))
            nfl(n) = 3
          elseif(i.ge.4 .and. i.le.5) then
            rr1(n) = max(r1,thrs(i))
            rr1(n) = min(rr1(n),thrs(i+1))
            rr2(n) = max(r2,thrs(i))
            rr2(n) = min(rr2(n),thrs(i+1))
            nfl(n) = i
          elseif(i.eq.6) then
            rr1(n) = max(r1,thrs(6))
            rr2(n) = r2
            nfl(n) = 6
          endif
        enddo

      else

C--     Downward evolution
C--     If the start point is at a flavor threshold nf1 --> nf1-1
        if(ihit1.ne.0) nf1 = nf1-1
        n = 0
        do i = nf1,nf2,-1
          n = n+1
          if(i.eq.6) then
            rr1(n) = r1
            rr2(n) = max(r2,thrs(6))
            nfl(n) = 6
          elseif(i.ge.4 .and. i.le.5) then
            rr1(n) = max(r1,thrs(i))
            rr1(n) = min(rr1(n),thrs(i+1))
            rr2(n) = max(r2,thrs(i))
            rr2(n) = min(rr2(n),thrs(i+1))
            nfl(n) = i
          elseif(i.eq.3) then
            rr1(n) = min(r1,thrs(4))
            rr2(n) = r2
            nfl(n) = 3
          endif
        enddo
      
      endif  

      return
      end
