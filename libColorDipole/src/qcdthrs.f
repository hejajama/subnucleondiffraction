
C--  Threshold routines on the file qcdthrs.f
C--  ----------------------------------------

C--   subroutine sqcThrFFNS(nf)
C--   subroutine sqcThrMFNS(nf,rc,rb,rt)
C--   subroutine sqcThrVFNS(iqc,iqb,iqt)
C--   subroutine sqcRmass2(fthr,rthr)

C     =========================
      subroutine sqcThrFFNS(nf)
C     =========================

C--   Define thresholds for the fixed flavor number scheme
C--
C--   Input   nf           number of flavors [3-6]
C--   Output  qthrs6(4:6)  thresholds on the mu2 grid
C--           tthrs6(4:6)  thresholds on t = log(mu2) grid
C--

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      if    (nf.eq.3) then
        qthrs6(4) = 4000.0*qlimu6
        qthrs6(5) = 5000.0*qlimu6
        qthrs6(6) = 6000.0*qlimu6
      elseif(nf.eq.4) then
        qthrs6(4) = 0.0004*qlimd6
        qthrs6(5) = 5000.0*qlimu6
        qthrs6(6) = 6000.0*qlimu6
      elseif(nf.eq.5) then
        qthrs6(4) = 0.0004*qlimd6
        qthrs6(5) = 0.0005*qlimd6
        qthrs6(6) = 6000.0*qlimu6
      elseif(nf.eq.6) then
        qthrs6(4) = 0.0004*qlimd6
        qthrs6(5) = 0.0005*qlimd6
        qthrs6(6) = 0.0006*qlimd6
      endif

      do i = 4,6
        rthrs6(i) = qthrs6(i)
        tthrs6(i) = log(qthrs6(i))
      enddo

      nfix6  = nf
      Lmfns6 = .false.
  
      return
      end
      
C     ==================================      
      subroutine sqcThrMFNS(nf,rc,rb,rt)
C     ================================== 

C--   Mixed flavour number scheme thresholds

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'
      
C--   Set nfix6 and thresholds on the fact scale      
      call sqcThrFFNS(nf)
C--   Set thresholds on the renormalisation scale
      rthrs6(4) = rc
      rthrs6(5) = rb
      rthrs6(6) = rt
      
C--   Flag MVNS
      Lmfns6 = .true.     
      
      return
      end

C     ==================================
      subroutine sqcThrVFNS(iqc,iqb,iqt)
C     ==================================

C--   Set thresholds in the vfns
C--
C--   Input   iqc,iqb,iqt  
C--   Output  qthrs6(4:6)  thresholds on the mu2 grid
C--
C--   Can ony be called when t-grid is defined

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

C--   Initialize all thresholds above grid
      qthrs6(4) = 4000.0*qlimu6
      qthrs6(5) = 5000.0*qlimu6
      qthrs6(6) = 6000.0*qlimu6

C--   Check iqc, iqb, iqt in succession and drop out if not in range 
      if(iqc.ge.1 .and. iqc.le.ntt2) then
C--     iqc in range
        qthrs6(4) = exp(tgrid2(iqc))
        if(iqb.gt.iqc .and. iqb.le.ntt2) then
C--       iqb in range
          qthrs6(5) = exp(tgrid2(iqb))
          if(iqt.gt.iqb .and. iqt.le.ntt2) then
C--         iqt in range
            qthrs6(6) = exp(tgrid2(iqt))
          endif
        endif
      endif

      do i = 4,6
        tthrs6(i) = log(qthrs6(i))
      enddo     

      nfix6  = 0
      Lmfns6 = .false.

      return
      end

C     ===============================
      subroutine sqcRmass2(fthr,rthr)
C     ===============================

C--   Convert threshold on the factorization scale to the renormalization
C--   scale, using the scale factors stored in qpars6 
C--
C--   Input   fthr(4:6)  thresholds fcbt2 defined on the fact  scale
C--   Output  rmas(4:6)  thresholds rcbt2 defined on the renor scale

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qpars6.inc'

      dimension fthr(4:6),rthr(4:6)

      do i = 4,6
        rthr(i) = aar6*fthr(i) + bbr6
      enddo

      return
      end
