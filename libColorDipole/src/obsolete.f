
C--   This is the file obsolete.f with obsolete routines

C--   function alfunc(iord,as0,r20,r2,ierr)                 [obsolete]
C--   double precision function AsEvol
C--                  (iord,as0,r20,r2,iqcdnum,nfout,ierr)   [obsolete]
C--   subroutine Evolve(func,def,iq0,epsi)                  [obsolete]
C--   subroutine EvolFF(func,def,iq0,epsi)                  [obsolete}
C--   subroutine NSevol(ityp,func,idf,iq0,epsi)             [obsolete]
C--   subroutine EvolNS(ityp,func,idf,iq0,epsi)             [obsolete]
C--   subroutine SGEvol(funf,idf,fung,idg,iq0,epsi)         [obsolete]
C--   subroutine EvolSG(funf,idf,fung,idg,iq0,epsi)         [obsolete]
C--   subroutine getids(idmin,idmax,nwords)                 [obsolete]
C--   subroutine dumpwt(lun,file)                           [obsolete]
C--   subroutine FastFac(id,funxq)                          [obsolete]
C--   subroutine FtimesW(w,fun,jd1,id2,iadd)                [obsolete]
C--   double precision function FcrossC(w,id,idf,ix,iq)     [obsolete]
C--   subroutine SetPdf(itype)                              [obsolete]
C--   subroutine GetPdf(itype)                              [obsolete]
C--   subroutine FastPdf(id,coef)                           [obsolete]
C--   subroutine FastFcC(w,id,id1,id2)                      [obsolete]
C--   subroutine FastAdd(id1,id2)                           [obsolete]
C--   subroutine setabq(aq,bq)                              [obsolete]
C--   subroutine getabq(aq,bq)                              [obsolete]
C--   subroutine setthr(nfix,q2c,q2b,q2t)                   [obsolete]
C--   subroutine getthr(nfix,q2c,q2b,q2t,iqcdnum)           [obsolete]
C--   subroutine pdfsum(id,wt,n,idout)                      [obsolete]
C--   double precision function pdfval(idf,xx,qq,mode)      [obsolete]
C--   double precision function pgluon(x,qmu2,jchk)         [obsolete]
C--   double precision function onepdf(x,qmu2,def,jchk)     [obsolete]
C--   subroutine allpdf(x,qmu2,pdf,imode)                   [obsolete]
C--   subroutine pdfsxq(xx,qq,pdf,jchk)                     [obsolete]
C--   subroutine pdfsij(ix,iq,pdf,jchk)                     [obsolete]
C--   double precision function sgnsxq(id,xx,qq,jchk)       [obsolete]
C--   double precision function sgnsij(id,ix,iq,jchk)       [obsolete]
C--   double precision function pfunxq(id,xx,qq,jchk)       [obsolete]
C--   double precision function pfunij(id,ix,iq,jchk)       [obsolete]
C--   double precision function psumxq(def,xx,qq,jchk)      [obsolete]
C--   double precision function psumij(def,ix,iq,jchk)      [obsolete]
C--   subroutine stfval(istf,ityp,id,x,q,f,n,mode)          [obsolete]
C--   subroutine StrFun(istf,user,x,q,f,n,mode)             [obsolete]
C--   subroutine StampIt('string')                          [obsolete]
C--   subroutine DumpTab(w,lun,file)                        [obsolete]
C--   subroutine ReadTab(w,nw,lun,file,nwords,ierr)         [obsolete]
C--   integer function idSplij(string)                      [obsolete]

C     ======================================================
      double precision function alfunc(iord,as0,r20,r2,ierr)
C     ======================================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'ALFUNC ( IORD, AS0, R20, R2, IERR )'/

      alfunc = 0.D0
      idum   = iord
      dum    = as0
      dum    = r20
      dum    = r2
      ierr   = 999

      call sqcErrMsg(subnam,
     +              'Alfunc obsolete, pls use EvolAs instead')

      return
      end
      
C     ==========================================================
      double precision function AsEvol
     +                      (iord,as0,r20,r2,iqcdnum,nfout,ierr)
C     ==========================================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam 
     +  /'ASEVOL ( IORD, AS0, R20, R2, INTERN, NFOUT, IERR )'/

      asevol = 0.D0
      idum   = iord
      dum    = as0
      dum    = r20
      dum    = r2
      idum   = iqcdnum
      nfout  = 0
      ierr   = 999
           
      call sqcErrMsg(subnam,
     +              'AsEvol obsolete, pls use EvolAs instead')

      return
      end      
      
C     ====================================
      subroutine Evolve(func,def,iq0,epsi)
C     ====================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'EVOLVE ( FUNC, DEF, IQ0, EPSI )'/

      dum  = func
      dum  = def
      idum = iq0
      dum  = epsi

      call sqcErrMsg(subnam,
     +              'Evolve obsolete, pls use EvolFG instead')

      return
      end

C     ====================================
      subroutine EvolFF(func,def,iq0,epsi)
C     ====================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'EVOLFF ( FUNC, DEF, IQ0, EPSI )'/

      dum  = func
      dum  = def
      idum = iq0
      dum  = epsi

      call sqcErrMsg(subnam,
     +              'EVOLFF obsolete, pls use EVOLFG instead')

      return
      end   
      
C     =========================================
      subroutine NSevol(ityp,func,idf,iq0,epsi)
C     =========================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'NSEVOL ( ITYP, FUNC, IDF, IQ0, EPSI )'/

      call sqcErrMsg(subnam,'NSEVOL obsolete, use EVOLFG instead')

      idum = ityp
      dum  = func
      idum = idf
      idum = iq0
      epsi = 1.D11

      return
      end

C     =========================================
      subroutine EvolNS(ityp,func,idf,iq0,epsi)
C     =========================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'EVOLNS ( ITYP, FUNC, IDF, IQ0, EPSI )'/

      call sqcErrMsg(subnam,'EVOLNS obsolete, use EVOLFG instead')

      idum = ityp
      dum  = func
      idum = idf
      idum = iq0
      epsi = 1.D11

      return
      end
            
C     =============================================
      subroutine SGEvol(funf,idf,fung,idg,iq0,epsi)
C     =============================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'SGEVOL ( SFUN, IDS, GFUN, IDG, IQ0, EPSI )'/

      call sqcErrMsg(subnam,'SGEVOL obsolete, use EVOLFG instead')

      dum  = funf
      idum = idf
      dum  = fung
      idum = idg
      idum = iq0
      epsi = 1.D11

      return
      end

C     =============================================
      subroutine EvolSG(funf,idf,fung,idg,iq0,epsi)
C     =============================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'EVOLSG ( SFUN, IDS, GFUN, IDG, IQ0, EPSI )'/

      call sqcErrMsg(subnam,'EVOLSG obsolete, use EVOLFG instead')

      dum  = funf
      idum = idf
      dum  = fung
      idum = idg
      idum = iq0
      epsi = 1.D11

      return
      end                  

C     =====================================
      subroutine getids(idmin,idmax,nwords)
C     =====================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      
      character*80 subnam
      data subnam /'GETIDS ( IDMIN, IDMAX, NWORDS )'/

      call sqcErrMsg(subnam,
     +   'GETIDS obsolete, please use NWUSED( nwtot, nwuse, nwtab )')

      idum  = idmin               !avoid compiler warning
      idum  = idmax               !avoid compiler warning
      idum  = nwords              !avoid compiler warning

      return
      end

C     ===========================
      subroutine dumpwt(lun,file)
C     ===========================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      
      character*(*) file

      character*80 subnam
      data subnam /'DUMPWT ( LUN, FILE )'/

      call sqcErrMsg(subnam,
     +   'DUMPWT obsolete, please use DMPWGT(itype,lun,file)')

      idum  = lun               !avoid compiler warning
      leng  = imb_lenoc(file)   !avoid compiler warning

      return
      end

C     ============================
      subroutine FastFac(id,funxq)
C     ============================

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'FASTFAC ( ID, FUNXQ )'/

      call sqcErrMsg(subnam,
     + 'FASTFAC obsolete, please use FASTKIN')

      idum  = id              !avoid compiler warning
      ddum  = funxq           !avoid compiler warning

      return
      end

C     ======================================
      subroutine FtimesW(w,fun,jd1,id2,iadd)
C     ======================================

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'FTIMESW ( W, FUN, ID1, ID2, IADD )'/

      call sqcErrMsg(subnam,
     + 'FTIMESW obsolete, please use WTIMESF')

      ddum  = w                !avoid compiler warning
      ddum  = fun              !avoid compiler warning
      idum  = jd1              !avoid compiler warning
      idum  = id2              !avoid compiler warning
      idum  = iadd             !avoid compiler warning

      return
      end
      
C     =================================================
      double precision function FcrossC(w,id,idf,ix,iq)
C     =================================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      
      character*80 subnam
      data subnam /'FCROSSC ( W, IDW, IDF, IX, IQ )'/

      
      dum     = w
      idum    = id
      idum    = idf
      idum    = ix
      idum    = iq
      FcrossC = 0.D0

      call sqcErrMsg(subnam,'FCROSSC obsolete please use FCROSSK')

      return
      end      

C     ========================
      subroutine SetPdf(itype)
C     ========================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'SETPDF ( ITYPE )'/

      call sqcErrMsg(subnam,'SETPDF obsolete')

      idum  = itype            !avoid compiler warning

      return
      end
      
C     ========================
      subroutine GetPdf(itype)
C     ========================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'GETPDF ( ITYPE )'/

      call sqcErrMsg(subnam,'GETPDF obsolete')

      idum  = itype            !avoid compiler warning

      return
      end      

C     ===========================
      subroutine FastPdf(id,coef)
C     ===========================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'FASTPDF ( ID, COEF )'/

      call sqcErrMsg(subnam,
     +   'FASTPDF obsolete, please use FASTSUM( iset, coef, id )')

      idum  = id               !avoid compiler warning
      ddum  = coef             !avoid compiler warning

      return
      end

C     ==================================
      subroutine FastFcC(w,idwt,id1,id2)
C     ==================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'FASTFCC ( W, ID, ID1, ID2 )'/

      call sqcErrMsg(subnam,
     +   'FASTFCC obsolete, please use FASTFXK')

      ddum  = w          !avoid compiler warning
      idum  = idwt       !avoid compiler warning
      idum  = id1       !avoid compiler warning
      idum  = id2       !avoid compiler warning

      return
      end

C     ===========================
      subroutine FastAdd(id1,id2)
C     ===========================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'FASTADD ( ID1, ID2 )'/

      call sqcErrMsg(subnam,
     +   'FASTADD obsolete, please use FASTCPY ( id1, id2, iadd )')

      idum  = id1       !avoid compiler warning
      idum  = id2       !avoid compiler warning

      return
      end

C     ========================
      subroutine setabq(aq,bq)
C     ========================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'SETABQ ( AQ, BQ )'/

      call sqcErrMsg(subnam,'SETABQ obsolete')

      dum = aq
      dum = bq
      
      return
      end
      
C     ========================
      subroutine getabq(aq,bq)
C     ========================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'GETABQ ( AQ, BQ )'/

      call sqcErrMsg(subnam,'GETABQ obsolete')

      dum = aq
      dum = bq
      
      return
      end

C     ===================================
      subroutine setthr(nfix,q2c,q2b,q2t)
C     ===================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'SETTHR ( NFIX, Q2C, Q2B, Q2T )'/

      call sqcErrMsg(subnam,
     +              'SETTHR obsolete, pls use SETCBT instead')

      idum = nfix
      dum  = q2c
      dum  = q2t
      dum  = q2b
      
      return
      end
      
C     ==========================================
      subroutine getthr(nfix,q2c,q2b,q2t,intern)
C     ==========================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'GETTHR ( NFIX, Q2C, Q2B, Q2T, INTERN )'/

      call sqcErrMsg(subnam,
     +              'GETTHR obsolete, pls use GETCBT instead')

      idum = nfix
      dum  = q2c
      dum  = q2b
      dum  = q2t
      idum = intern

      return
      end

C     ================================
      subroutine pdfsum(id,wt,n,idout)
C     ================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'PDFSUM ( IDIN, WGT, N, IDOUT )'/

      call sqcErrMsg(subnam,'PDFSUM obsolete')

      idum  = id
      dum   = wt
      idum  = n
      idout = 999 

      return
      end
      
C     ================================================
      double precision function pdfval(idf,xx,qq,mode)
C     ================================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'PDFVAL ( ID, X, Q2, MODE )'/

      pdfval = 0.D0
      idum   = idf
      dum    = xx
      dum    = qq
      idum   = mode

      call sqcErrMsg(subnam,'PDFVAL obsolete')

      return
      end
            
C     ============================================
      double precision function pgluon(xx,qq,jchk)
C     ============================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'PGLUON ( X, QMU2, ICHK )'/

      pgluon = 0.D0
      dum    = xx
      dum    = qq
      idfum  = jchk
   
      call sqcErrMsg(subnam,'PGLUON obsolete')

      return
      end

C     ================================================
      double precision function onepdf(xx,qq,def,jchk)
C     ================================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'ONEPDF ( X, QMU2, DEF, ICHK )'/

      onepdf = 0.D0
      dum    = xx
      dum    = qq
      dum    = def
      jchk   = 999

      call sqcErrMsg(subnam,'ONEPDF obsolete')

      return
      end

C     ==================================
      subroutine allpdf(xx,qq,pdf,imode)
C     ==================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'ALLPDF ( X, QMU2, VAL, IMODE )'/

      dum   = xx
      dum   = qq
      dum   =  pdf
      idum  = imode

      call sqcErrMsg(subnam,'ALLPDF obsolete please call FPDFXQ')

      return
      end
      
C     ==================================
      subroutine pdfsxq(xx,qq,pdf,jchk)
C     ==================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'PDFSXQ ( X, QMU2, PDF, ICHK )'/

      dum   = xx
      dum   = qq
      dum   = pdf
      idum  = jchk

      call sqcErrMsg(subnam,'PDFSXQ obsolete please call FPDFXQ')

      return
      end
      
C     ==================================
      subroutine pdfsij(ix,iq,pdf,jchk)
C     ==================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'PDFSIJ ( IX, IQ, PDF, ICHK )'/

      idum  = ix
      idum  = iq
      dum   = pdf
      idum  = jchk

      call sqcErrMsg(subnam,'PDFSIJ obsolete please call FPDFIJ')

      return
      end                   

C     ===============================================
      double precision function sgnsxq(id,xx,qq,jchk)
C     ===============================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'SGNSXQ ( ID, X, QMU2, ICHK )'/

      sgnsxq = 0.D0
      idum   = id
      dum    = xx
      dum    = qq
      idum   = jchk

      call sqcErrMsg(subnam,'SGNSXQ obsolete please call FSNSXQ')

      return
      end
       
C     ===============================================
      double precision function sgnsij(id,ix,iq,jchk)
C     ===============================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'SGNSIJ ( ID, IX, IQ, ICHK )'/

      sgnsij = 0.D0
      idum   = id
      idum   = ix
      idum   = iq
      idum   = jchk

      call sqcErrMsg(subnam,'SGNSIJ obsolete please call FSNSIJ')

      return
      end                    

C     ===============================================
      double precision function pfunxq(id,xx,qq,jchk)
C     ===============================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'PFUNXQ ( ID, X, QMU2, ICHK )'/

      pfunxq = 0.D0
      idum   = id
      dum    = xx
      dum    = qq
      idum   = jchk

      call sqcErrMsg(subnam,'PFUNXQ obsolete please call FVALXQ')

      return
      end
      
C     ===============================================
      double precision function pfunij(id,ix,iq,jchk)
C     ===============================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'PFUNIJ ( ID, IX, IQ, ICHK )'/

      pfunij = 0.D0
      idum   = id
      idum   = ix
      idum   = iq
      idum   = jchk

      call sqcErrMsg(subnam,'PFUNIJ obsolete please call FVALIJ')

      return
      end         

C     ================================================
      double precision function psumxq(def,xx,qq,jchk)
C     ================================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'PSUMXQ ( DEF, X, QMU2, ICHK )'/

      psumxq = 0.D0
      dum    = def
      dum    = xx
      dum    = qq
      idum   = jchk

      call sqcErrMsg(subnam,'PSUMXQ obsolete please call FSUMXQ')

      return
      end
                        
C     ================================================
      double precision function psumij(def,ix,iq,jchk)
C     ================================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'PSUMIJ ( DEF, IX, IQ, ICHK )'/

      psumij = 0.D0
      dum    = def
      idum   = ix
      idum   = iq
      idum   = jchk

      call sqcErrMsg(subnam,'PSUMIJ obsolete please call FSUMIJ')

      return
      end        

C     ============================================
      subroutine stfval(istf,ityp,id,x,q,f,n,mode)
C     ============================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      character*80 subnam
      data subnam /'STFVAL ( ISTF, ITYP, ID, X, Q2, F, n, MODE )'/

      idum = istf
      idum = ityp
      idum = id
      dum  = x
      dum  = q
      dum  = f
      idum = n
      idum = mode

      call sqcErrMsg(subnam,'STFVal obsolete, pls use ZMSTF package')

      return
      end

C     =========================================    
      subroutine StrFun(istf,user,x,q,f,n,mode)
C     =========================================

C--   Obsolete

      implicit double precision (a-h,o-z)
      
      include 'qluns1.inc'

      character*80 subnam
      data subnam /'STRFUN ( ISTF, DEF, X, Q2, F, N, MODE )'/

      idum = istf
      dum  = user
      dum  = x
      dum  = q
      dum  = f
      idum = n
      idum = mode
      
      call sqcErrMsg(subnam,'STRFUN obsolete, pls use ZMSTF package')
      
      return
      end

C     ==========================
      subroutine StampIt(string)
C     ==========================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      
      character*(*) string

      character*80 subnam
      data subnam /'STAMPIT ( STRING )'/
      len = imb_lenoc(string)   !avoid compiler warning
      call sqcErrMsg(subnam,'STAMPIT obsolete')

      return
      end

C     ==============================
      subroutine DumpTab(w,lun,file)
C     ==============================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      
      dimension w(*)
      character*(*) file

      character*80 subnam
      data subnam /'DUMPTAB ( W, LUN, FILE )'/
      dum = w(1)            !avoid compiler warning
      kun = lun             !avoid compiler warning
      len = imb_lenoc(file) !avoid compiler warning
      call sqcErrMsg(subnam,'DUMPTAB obsolete, use TABDUMP instead')

      return
      end      

C     =============================================
      subroutine ReadTab(w,nw,lun,file,nwords,ierr)
C     =============================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      
      dimension w(*)
      character*(*) file

      character*80 subnam
      data subnam /'READTAB ( W, NW, LUN, FILE, NWORDS, IERR )'/
      
      dum = w(1)            !avoid compiler warning
      iii = nw              !avoid compiler warning
      kun = lun             !avoid compiler warning
      len = imb_lenoc(file) !avoid compiler warning
      jjj = nwords          !avoid compiler warning
      kkk = ierr            !avoid compiler warning
      call sqcErrMsg(subnam,'READTAB obsolete, use TABREAD instead')

      return
      end
      
C     ================================
      integer function idSplij(string)
C     ================================

C--   Obsolete

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      
      character*80 subnam
      data subnam /'IDSPLIJ ( STRING )'/
      
      character*(*) string

      len     = imb_lenoc(string)   !avoid compiler warning
      idSplij = 0

      call sqcErrMsg(subnam,'IDSPLIJ obsolete please use IDSPFUN')

      return
      end
      


      
      
