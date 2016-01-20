
C--   This is the file usrerr.f containing the error bits and messages

C--   subroutine sqcBitini
C--   subroutine sqcMakefl(subnam,ichk,iset,idel,iall)
C--   subroutine sqcChkIni(subnam)
C--   subroutine sqcChkflg(jset,ichk,subnam)
C--   subroutine sqcChekit(jset,ichk,jbit)
C--   subroutine sqcSetflg(iset,idel,iall)
C--   subroutine sqcSetbit(ibit,iword,n)
C--   subroutine sqcDelbit(ibit,iword,n)
C--   integer function iqcGetbit(ibit,iword,n)
C--   subroutine sqcIlele(subnam,parnam,imin,ival,imax,comment)
C--   subroutine sqcIlelt(subnam,parnam,imin,ival,imax,comment)
C--   subroutine sqcIltle(subnam,parnam,imin,ival,imax,comment)
C--   subroutine sqcIltlt(subnam,parnam,imin,ival,imax,comment)
C--   subroutine sqcDlele(subnam,parnam,dmin,dval,dmax,comment)
C--   subroutine sqcDlelt(subnam,parnam,dmin,dval,dmax,comment)
C--   subroutine sqcDltle(subnam,parnam,dmin,dval,dmax,comment)
C--   subroutine sqcDltlt(subnam,parnam,dmin,dval,dmax,comment)
C--   subroutine sqcErrMsg(subnam,message)
C--   subroutine sqcSubMsg(subnam,message)
C--   subroutine sqcIdEmsg(subnam,parnam,idval,ierr)
C--   subroutine sqcCutMsg(subnam,x,qmu2,ifail,noextra)

C==   ===============================================================
C==   Manage status bits and error messages =========================
C==   ===============================================================
      
C     ====================
      subroutine sqcBitini
C     ====================

C--   Assign bits and define error messages
C--   Called by qcinit

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qibit4.inc'

C--   !!!!Warning!!!! increase dimension of errmsg3 if # bits > 10
      character*35    errmsg3,commsg3
      common /qemsg3/ errmsg3(10),commsg3(10)

C--                      '-----------------------------------'
      ibInit4          = 1
      errmsg3(ibInit4) = 'QCDNUM not initialized             '
      commsg3(ibInit4) = 'Please call QCINIT                 ' 

      ibXgri4          = 2
      errmsg3(ibXgri4) = 'No x-grid available                '
      commsg3(ibXgri4) = 'Please call GXMAKE                 '

      ibQgri4          = 3
      errmsg3(ibQgri4) = 'No Q2-grid available               '
      commsg3(ibQgri4) = 'Please call GQMAKE                 '

      ibWeit4          = 4
      errmsg3(ibWeit4) = 'No weight tables available         '
      commsg3(ibWeit4) = 'Please call FILLWT or READWT       '

      ibPdfs4          = 5
      errmsg3(ibPdfs4) = 'No pdfs available                  '
      commsg3(ibPdfs4) = 'Please call EVOLFF or PDFINP       '
      
      ibUtab4          = 6
      errmsg3(ibUtab4) = 'No user tables initialized         '
      commsg3(ibUtab4) = 'Please call BOOKTAB                '

      ibUwgt4          = 7
      errmsg3(ibUwgt4) = 'No user weights available          '
      commsg3(ibUwgt4) = 'Please call MAKEWTx or READTAB     '
      
      ibFbuf4          = 8
      errmsg3(ibFbuf4) = 'No scratch buffers avaialble       '
      commsg3(ibFbuf4) = 'Please call FASTINI                '
      
C--   NB: if new bits are defined like ibXXX4 then do not forget
C--       to declare them in qibit4.inc !!!!       

      return
      end

C     ===========================================
      subroutine sqcMakefl(subnam,ichk,iset,idel)
C     ===========================================

C--   Set bitpatterns for user subroutines (large switchyard)
C--   Called at the first invocation of each user routine
C--
C--   Input:  subnam subroutine name
C--
C--   Output: ichk(mbp0)  pattern of check-bits (check memory bit) 
C--           iset(mbp0)  pattern of set-bits   (set   memory bit)
C--           idel(mbp0)  pattern of del-bits   (unset memory bit)
C--
C--   Common: The bit assignments are defined in sqcBitini and stored
C--           in the common block /qibit4/    

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qsflg4.inc'
      include 'qibit4.inc'
      
      character*(*)   subnam

      dimension ichk(mbp0), iset(mbp0), idel(mbp0)

C--   Check if qcdnum is initialized
      call sqcChkIni(subnam)

C--   Initialize
      do i = 1,mbp0
        ichk(i) = 0
        iset(i) = 0
        idel(i) = 0
      enddo

      if    (subnam(1:6) .eq. 'SETVAL') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:6) .eq. 'SETINT') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:6) .eq. 'SETORD') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

        call sqcSetbit(ibPdfs4,idel,mbp0)   !invalidate pdfs 
        
      elseif(subnam(1:6) .eq. 'SETALF') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

        call sqcSetbit(ibPdfs4,idel,mbp0)   !invalidate pdfs 
        
      elseif(subnam(1:6) .eq. 'SETCBT') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

        call sqcSetbit(ibPdfs4,idel,mbp0)   !invalidate pdfs 
        
      elseif(subnam(1:6) .eq. 'MIXFNS') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

        call sqcSetbit(ibPdfs4,idel,mbp0)   !invalidate pdfs 
        
      elseif(subnam(1:6) .eq. 'GETCBT') then
C--   ----------------------------------
        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:6) .eq. 'SETABR') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

        call sqcSetbit(ibPdfs4,idel,mbp0)   !invalidate pdfs 
        
      elseif(subnam(1:6) .eq. 'GXMAKE') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

        call sqcSetbit(ibXgri4,iset,mbp0)   !validate x-grid

        call sqcSetbit(ibWeit4,idel,mbp0)   !invalidate weight tables
        call sqcSetbit(ibPdfs4,idel,mbp0)   !invalidate pdfs 
        call sqcSetbit(ibUtab4,idel,mbp0)   !invalidate user tables
        call sqcSetbit(ibUwgt4,idel,mbp0)   !invalidate user weights
        
      elseif(subnam(1:6) .eq. 'IXFRMX') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists

      elseif(subnam(1:6) .eq. 'XFRMIX') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists

      elseif(subnam(1:6) .eq. 'XXATIX') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists

      elseif(subnam(1:6) .eq. 'GQMAKE') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

        call sqcSetbit(ibQgri4,iset,mbp0)   !validate q2-grid 
        
        call sqcSetbit(ibWeit4,idel,mbp0)   !invalidate weight tables
        call sqcSetbit(ibPdfs4,idel,mbp0)   !invalidate pdfs
        call sqcSetbit(ibUtab4,idel,mbp0)   !invalidate user tables 
        call sqcSetbit(ibUwgt4,idel,mbp0)   !invalidate user weights
        
      elseif(subnam(1:6) .eq. 'IQFRMQ') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:6) .eq. 'QFRMIQ') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:6) .eq. 'QQATIQ') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:6) .eq. 'GRPARS') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:6) .eq. 'GXCOPY') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists

      elseif(subnam(1:6) .eq. 'GQCOPY') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:6) .eq. 'FILLWT') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

        call sqcSetbit(ibWeit4,iset,mbp0)   !validate weights

        call sqcSetbit(ibPdfs4,idel,mbp0)   !invalidate pdfs
        
      elseif(subnam(1:6) .eq. 'FILLWC') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

        call sqcSetbit(ibWeit4,iset,mbp0)   !validate weights 

        call sqcSetbit(ibPdfs4,idel,mbp0)   !invalidate pdfs

      elseif(subnam(1:6) .eq. 'DMPWGT') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibWeit4,ichk,mbp0)   !check weights exist

      elseif(subnam(1:6) .eq. 'READWT') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

        call sqcSetbit(ibWeit4,iset,mbp0)   !validate weights 
        
        call sqcSetbit(ibPdfs4,idel,mbp0)   !invalidate pdfs

      elseif(subnam(1:6) .eq. 'NWUSED') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        
      elseif(subnam(1:6) .eq. 'SETCUT') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        
        call sqcSetbit(ibPdfs4,idel,mbp0)   !invalidate pdfs

      elseif(subnam(1:6) .eq. 'GETCUT') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        
      elseif(subnam(1:6) .eq. 'LPASSC') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists        
        
      elseif(subnam(1:6) .eq. 'ASFUNC') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:6) .eq. 'EVOLAS') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:6) .eq. 'EVOLFG') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        call sqcSetbit(ibWeit4,ichk,mbp0)   !check weight tables exist

        call sqcSetbit(ibPdfs4,iset,mbp0)   !validate pdfs
        
      elseif(subnam(1:6) .eq. 'PDFINP') then
C--   ----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

        call sqcSetbit(ibPdfs4,iset,mbp0)   !validate pdfs
        
      elseif(subnam(1:6) .eq. 'FPDFXQ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist

      elseif(subnam(1:6) .eq. 'FPDFIJ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist

      elseif(subnam(1:6) .eq. 'FSNSXQ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist

      elseif(subnam(1:6) .eq. 'FSNSIJ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist

      elseif(subnam(1:6) .eq. 'FVALXQ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist

      elseif(subnam(1:6) .eq. 'FVALIJ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist

      elseif(subnam(1:6) .eq. 'FSUMXQ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist

      elseif(subnam(1:6) .eq. 'FSUMIJ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist

      elseif(subnam(1:6) .eq. 'FSPLNE') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist        
        
      elseif(subnam(1:6) .eq. 'SPLCHK') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist
        
      elseif(subnam(1:6) .eq. 'PDFLST') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist        
        
      elseif(subnam(1:6) .eq. 'PDFTAB') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist

      elseif(subnam(1:6) .eq. 'STRFUN') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist

      elseif(subnam(1:7) .eq. 'BOOKTAB') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        
        call sqcSetbit(ibUtab4,iset,mbp0)   !validate user tables 
        
        call sqcSetbit(ibUwgt4,idel,mbp0)   !invalidate user weights
        
      elseif(subnam(1:7) .eq. 'SETWPAR') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibUtab4,ichk,mbp0)   !check user tables exist

      elseif(subnam(1:7) .eq. 'GETWPAR') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibUtab4,ichk,mbp0)   !check user tables exist

      elseif(subnam(1:7) .eq. 'TABDUMP') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        call sqcSetbit(ibUwgt4,ichk,mbp0)   !check user weights exist

      elseif(subnam(1:7) .eq. 'TABREAD') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

        call sqcSetbit(ibUtab4,iset,mbp0)   !validate user tables 
        call sqcSetbit(ibUwgt4,iset,mbp0)   !validate user weights 

      elseif(subnam(1:7) .eq. 'MAKEWTA') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibUtab4,ichk,mbp0)   !check user tables exist 
        
        call sqcSetbit(ibUwgt4,iset,mbp0)   !validate user weights
        
      elseif(subnam(1:7) .eq. 'MAKEWTB') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibUtab4,ichk,mbp0)   !check user tables exist
        
        call sqcSetbit(ibUwgt4,iset,mbp0)   !validate user weights
        
      elseif(subnam(1:7) .eq. 'MAKEWRS') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibUtab4,ichk,mbp0)   !check user tables exist

        call sqcSetbit(ibUwgt4,iset,mbp0)   !validate user weights 
        
      elseif(subnam(1:7) .eq. 'MAKEWTD') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibUtab4,ichk,mbp0)   !check user tables exist 

        call sqcSetbit(ibUwgt4,iset,mbp0)   !validate user weights
        
      elseif(subnam(1:7) .eq. 'MAKEWTX') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibUtab4,ichk,mbp0)   !check user tables exist 

        call sqcSetbit(ibUwgt4,iset,mbp0)   !validate user weights 
        
      elseif(subnam(1:7) .eq. 'SCALEWT') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibUtab4,ichk,mbp0)   !check user tables exist

      elseif(subnam(1:7) .eq. 'COPYWGT') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibUtab4,ichk,mbp0)   !check user tables initialized
        call sqcSetbit(ibUtab4,ichk,mbp0)   !check user tables exist
        
      elseif(subnam(1:7) .eq. 'WCROSSW') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibUtab4,ichk,mbp0)   !check user tables initialized
        call sqcSetbit(ibUtab4,ichk,mbp0)   !check user tables exist
        
      elseif(subnam(1:7) .eq. 'WTIMESF') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibUtab4,ichk,mbp0)   !check user tables initialized
        call sqcSetbit(ibUtab4,ichk,mbp0)   !check user tables exist
                 
      elseif(subnam(1:7) .eq. 'FCROSSK') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist
        call sqcSetbit(ibUtab4,ichk,mbp0)   !check user tables exist
        
      elseif(subnam(1:7) .eq. 'FCROSSF') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist
        call sqcSetbit(ibUtab4,ichk,mbp0)   !check user tables exist

      elseif(subnam(1:7) .eq. 'NFLAVOR') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:7) .eq. 'GETALFN') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists

      elseif(subnam(1:7) .eq. 'EFROMQQ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:7) .eq. 'QQFROME') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized

      elseif(subnam(1:7) .eq. 'STFUNXQ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist

      elseif(subnam(1:7) .eq. 'FASTINI') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibXgri4,ichk,mbp0)   !check xgrid exists
        call sqcSetbit(ibQgri4,ichk,mbp0)   !check qgrid exists
        
        call sqcSetbit(ibFbuf4,iset,mbp0)   !validate fast buffers
        
      elseif(subnam(1:7) .eq. 'FASTCLR') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibFbuf4,ichk,mbp0)   !check fast buffers exist
        
      elseif(subnam(1:7) .eq. 'FASTEPM') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibFbuf4,ichk,mbp0)   !check fast buffers exist
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist 
        
      elseif(subnam(1:7) .eq. 'FASTSNS') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibFbuf4,ichk,mbp0)   !check fast buffers exist
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist

      elseif(subnam(1:7) .eq. 'FASTSUM') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibFbuf4,ichk,mbp0)   !check fast buffers exist
        call sqcSetbit(ibPdfs4,ichk,mbp0)   !check pdfs exist
        
      elseif(subnam(1:7) .eq. 'FASTFXK') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibFbuf4,ichk,mbp0)   !check fast buffers exist
        
      elseif(subnam(1:7) .eq. 'FASTFXF') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibFbuf4,ichk,mbp0)   !check fast buffers exist

      elseif(subnam(1:7) .eq. 'FASTKIN') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibFbuf4,ichk,mbp0)   !check fast buffers exist
        
      elseif(subnam(1:7) .eq. 'FASTCPY') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibFbuf4,ichk,mbp0)   !check fast buffers exist

      elseif(subnam(1:7) .eq. 'FASTFXQ') then
C--   -----------------------------------

        call sqcSetbit(ibInit4,ichk,mbp0)   !check qcdnum initialized
        call sqcSetbit(ibFbuf4,ichk,mbp0)   !check fast buffers exist

      else
C--   ----

        goto 510

      endif
C--   -----

      return

C--   Error messages

 510  continue
      write(lunerr1,'(/'' sqcMAKEFL: unknown subroutine '',A10,'//
     +        '  '' ---> STOP'')') subnam(1:7)
      stop

      end

C     ============================
      subroutine sqcChkIni(subnam)
C     ============================

C--   Check if QCDNUM is initialized

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qsflg4.inc'
      
      character*(*)   subnam

      if(iniflg4.ne.123456) then
C--     Error message
        leng = imb_lenoc(subnam)
        write(lundef1,'(/1X,70(''-''))')
        write(lundef1,'('' Error in '',A,'' ---> STOP'')')
     +      subnam(1:leng)
        write(lundef1,'( 1X,70(''-''))')
        write(lundef1,'(
     +    '' QCDNUM not initialized (no call to QCINIT)'')')
        stop
      endif

      return
      end

C     ======================================
      subroutine sqcChkflg(jset,ichk,subnam)
C     ======================================

C--   Check status bits; abort with error message if not OK.
C--   Called by each user subroutine.
C--
C--   Input:  ichk(mbp0)            pattern of status bits to be checked
C--
C--   Common: istat4(mbp0,mset0)    QCDNUM status passed through /qcsflg/
C--           subnam(1:7)           name of calling routine

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qsnam3.inc'
      include 'qibit4.inc'
      include 'qsflg4.inc'
      include 'qstor7.inc'

      character*(*)   subnam
      character*35    errmsg3,commsg3
      common /qemsg3/ errmsg3(10),commsg3(10)
      
      character*35 etxt1(9)
      data etxt1 /'No unpolarised weights available   ',
     +            'No polarised weights available     ',
     +            'No fragmentation weights available ',
     +            'No custom weight tables available  ',
     +            'Cant evolve external pdf set 5     ', 
     +            'Cant evolve external pdf set 6     ', 
     +            'Cant evolve external pdf set 7     ', 
     +            'Cant evolve external pdf set 8     ', 
     +            'Cant evolve external pdf set 9     '/  
     
      character*35 etxt2(9)
      data etxt2 /'No unpolarised pdfs available      ',
     +            'No polarised pdfs available        ',
     +            'No fragmentation funcs available   ',
     +            'No custom pdfs available           ',
     +            'No external pdf set 5 available    ',
     +            'No external pdf set 6 available    ',
     +            'No external pdf set 7 available    ',
     +            'No external pdf set 8 available    ',
     +            'No external pdf set 9 available    '/  

      dimension ichk(mbp0)

C--   Loop over status words and check status bits
      do i = 1,mbp0
        ierr = imb_test1(ichk(i),istat4(i,jset))
        if(ierr.ne.0) then
          jword = i
          goto  500
        endif
      enddo

      return

C--   Error messages

 500  continue

C--   First status bit which failed
      jbit = 0
      do i = 1,32
        if(imb_gbitn(ichk(jword),i)      .eq.1 .and.
     +     imb_gbitn(istat4(jword,jset),i).eq.0   ) then
           jbit = i
           goto 510
        endif
      enddo

 510  continue
C--   Error message
      leng = imb_lenoc(subnam)
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,'('' Error in '',A,'' ---> STOP'')') 
     +     subnam(1:leng)
      write(lunerr1,'( 1X,70(''-''))')
      if(jbit.eq.0) then
        write(lunerr1,'('' No error message found'')')
      elseif(jbit.eq.ibWeit4) then
        write(lunerr1,'(1X,A35)') etxt1(jset)
        write(lunerr1,'(1X,A35)') commsg3((jword-1)*32+jbit)
      elseif(jbit.eq.ibPdfs4) then
        write(lunerr1,'(1X,A35)') etxt2(jset)
        write(lunerr1,'(1X,A35)') commsg3((jword-1)*32+jbit)
      elseif(jbit.eq.ibInit4) then
        write(lundef1,'(1X,A35)') errmsg3((jword-1)*32+jbit)
        write(lundef1,'(1X,A35)') commsg3((jword-1)*32+jbit)
      else
        write(lunerr1,'(1X,A35)') errmsg3((jword-1)*32+jbit)
        write(lunerr1,'(1X,A35)') commsg3((jword-1)*32+jbit)
      endif
      leng = imb_lenoc(usrnam3)
      if(leng.gt.0) then
        write(lunerr1,'(/1X,A,'' was called by '',A)')
     +       subnam(1:7),usrnam3(1:leng)
      endif

      stop

      end

C     ====================================
      subroutine sqcChekit(jset,ichk,jbit)
C     ====================================

C--   Check status bits as in sqcChkflg but no abort with error message.
C--   Called by each user function which should run quiet.
C--
C--   Input:  ichk(mbp0)            pattern of status bits to be checked
C--   Output: jbit                  status bit which failed: 0 = OK
C--
C--   Common: istat4(mbp0,mset0)    QCDNUM status passed through /qcsflg/

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qsflg4.inc'
      include 'qstor7.inc'

      dimension ichk(mbp0)

C--   Initialize
      jbit = 0

C--   Loop over status words and check status bits
      do i = 1,mbp0
        jerr = imb_test1(ichk(i),istat4(i,jset))
        if(jerr.ne.0) then
          jword = i
          goto  500
        endif
      enddo

      return

 500  continue
C--   Search for first status bit which failed
      jbit = 0
      do i = 1,32
        if(imb_gbitn(ichk(jword),i)      .eq.1 .and.
     +     imb_gbitn(istat4(jword,jset),i).eq.0   ) then
           jbit = i
           goto 510
        endif
      enddo

 510  continue
      return
      end
      
C     ====================================
      subroutine sqcSetflg(iset,idel,ipdf)
C     ====================================

C--   Update QCDNUM status words.
C--   Called by each user subroutine.
C--
C--   Input:  iset(mbp0)           pattern of status bits to be set to 1
C--           idel(mbp0)           pattern of status bits to be set to 0
C--           ipdf = 0             update status word of all pdf sets
C--                # 0             update status word of ipdf set
C--
C--   Common: istat4(mbp0,mset0)   QCDNUM status passed through /qcsflg/

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qsflg4.inc'
      include 'qstor7.inc'

      dimension iset(mbp0), idel(mbp0)
      
      if(ipdf.eq.0) then
        j1 = 1
        j2 = mset0
      else
        j1 = ipdf
        j2 = ipdf
      endif    

C--   Loop over status words
      do j = j1,j2
        do i = 1,mbp0
C--       Set bits to one
          istat4(i,j) = ior(iset(i),istat4(i,j))
C--       Set bits to zero
          istat4(i,j) = iand(not(idel(i)),istat4(i,j))
        enddo
      enddo

      return
      end

C==   ===============================================================
C==   Set and get status bits =======================================
C==   ===============================================================

C     ==================================
      subroutine sqcSetbit(ibit,iword,n)
C     ==================================

C--   Set bit in sequence of n integer words to 1.
C--
C--   Input:   ibit     integer in the range [1,n*32]
C--            iword(n) array of n 32bit integers

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      dimension iword(n)

C--   Which word?
      nwd = (ibit-1)/32 + 1
C--   Check
      if(nwd.lt.1 .or. nwd.gt.n) goto 500

C--   Which bit?
      ibt = mod(ibit-1,32) + 1
C--   Check
      if(ibt.lt.1 .or. ibt.gt.32) goto 510

C--   Set bit
      call smb_sbit1(iword(nwd),ibt)

      return

C--   Error messages

 500  continue
      write(lunerr1,
     +    '(/'' sqcSETBIT: iwd .gt. maxwd '',2I15,'//
     +    '  '' ---> STOP'')') nwd, n
      write(lunerr1,*) ' Input ibit = ', ibit
      write(lunerr1,*) ' Input n    = ', n
      stop

 510  continue
      write(lunerr1,
     +    '(/'' sqcSETBIT: ibt not in range [1,32] '',I5,'//
     +    '  '' ---> STOP'')') ibt
      write(lunerr1,*) ' Input  ibit = ', ibit
      write(lunerr1,*) ' Input  n    = ', n
      write(lunerr1,*) ' Output ibt  = ', ibt
      stop

      end

C     ==================================
      subroutine sqcDelbit(ibit,iword,n)
C     ==================================

C--   Set bit in sequence of n integer words to 0.
C--
C--   Input:   ibit     integer in the range [1,n*32]
C--            iword(n) array of n 32bit integers

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      dimension iword(n)

C--   Which word?
      nwd = (ibit-1)/32 + 1
C--   Check
      if(nwd.lt.1 .or. nwd.gt.n) goto 500

C--   Which bit?
      ibt = mod(ibit-1,32) + 1
C--   Check
      if(ibt.lt.1 .or. ibt.gt.32) goto 510

C--   Set bit
      call smb_sbit0(iword(nwd),ibt)

      return

C--   Error messages

 500  continue
      write(lunerr1,
     +    '(/'' sqcDELBIT: iwd .gt. maxwd '',2I5,'//
     +    '  '' ---> STOP'')') nwd, n
      stop

 510  continue
      write(lunerr1,
     +    '(/'' sqcDELBIT: ibt not in range [1,32] '',I5,'//
     +    '  '' ---> STOP'')') ibt
      stop

      end
 
C     ========================================
      integer function iqcGetbit(ibit,iword,n)
C     ========================================

C--   Get bit in sequence of n integer words.
C--
C--   Input:   ibit     integer in the range [1,n*32]
C--            iword(n) array of n 32bit integers

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'

      dimension iword(n)

C--   Which word?
      nwd = (ibit-1)/32 + 1
C--   Check
      if(nwd.lt.1 .or. nwd.gt.n) goto 500

C--   Which bit?
      ibt = mod(ibit-1,32) + 1
C--   Check
      if(ibt.lt.1 .or. ibt.gt.32) goto 510

C--   Value of bit
      iqcGetbit = imb_gbitn(iword(nwd),ibt)

      return

C--   Error messages

 500  continue
      write(lunerr1,
     +    '(/'' iqcGETBIT: iwd .gt. maxwd '',2I5,'//
     +    '  '' ---> STOP'')') nwd, n
      stop

 510  continue
      write(lunerr1,
     +    '(/'' iqcGETBIT: ibt not in range [1,32] '',I5,'//
     +    '  '' ---> STOP'')') ibt
      stop

      end

C==   ===============================================================
C==   Range checks and associated messages  =========================
C==   ===============================================================

C     =========================================================
      subroutine sqcIlele(subnam,parnam,imin,ival,imax,comment)
C     =========================================================

C--   OK if imin .le. ival .le. imax

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      include 'qsnam3.inc'

      character*10  cmin, cval, cmax
      character*(*) parnam, comment
      character*(*) subnam

C--   First for the check
      if( (imin.le.ival) .and. (ival.le.imax) ) return
C--   Now for the error message
      leng = imb_lenoc(subnam)
      call smb_ItoCh(imin,cmin,lmin)
      call smb_ItoCh(ival,cval,lval)
      call smb_ItoCh(imax,cmax,lmax)
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in ', subnam(1:leng), ' ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')
      write(lunerr1,*) parnam,' = ',cval(1:lval),' not in range [ ',
     +                 cmin(1:lmin),' , ',cmax(1:lmax),' ]'
      write(lunerr1,*) comment
      leng = imb_lenoc(usrnam3)
      if(leng.gt.0) then
        write(lunerr1,*) ' '
        write(lunerr1,*) subnam(1:7),
     +  ' was called by ',usrnam3(1:leng)
      endif

      stop
      end
      
C     =========================================================
      subroutine sqcIlelt(subnam,parnam,imin,ival,imax,comment)
C     =========================================================

C--   OK if imin .le. ival .lt. imax

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      include 'qsnam3.inc'

      character*10  cmin, cval, cmax
      character*(*) parnam, comment
      character*(*) subnam

C--   First for the check
      if( (imin.le.ival) .and. (ival.lt.imax) ) return
C--   Now for the error message
      leng = imb_lenoc(subnam)
      call smb_ItoCh(imin,cmin,lmin)
      call smb_ItoCh(ival,cval,lval)
      call smb_ItoCh(imax,cmax,lmax)
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in ', subnam(1:leng), ' ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')
      write(lunerr1,*) parnam,' = ',cval(1:lval),' not in range [ ',
     +                 cmin(1:lmin),' , ',cmax(1:lmax),' )'     
      write(lunerr1,*) comment
      leng = imb_lenoc(usrnam3)
      if(leng.gt.0) then
        write(lunerr1,*) ' '
        write(lunerr1,*) subnam(1:7),
     +  ' was called by ',usrnam3(1:leng)
      endif

      stop
      end
      
C     =========================================================
      subroutine sqcIltle(subnam,parnam,imin,ival,imax,comment)
C     =========================================================

C--   OK if imin .lt. ival .le. imax

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      include 'qsnam3.inc'

      character*10  cmin, cval, cmax
      character*(*) parnam, comment
      character*(*) subnam

C--   First for the check
      if( (imin.lt.ival) .and. (ival.le.imax) ) return
C--   Now for the error message
      leng = imb_lenoc(subnam)
      call smb_ItoCh(imin,cmin,lmin)
      call smb_ItoCh(ival,cval,lval)
      call smb_ItoCh(imax,cmax,lmax)
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in ', subnam(1:leng), ' ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')
      write(lunerr1,*) parnam,' = ',cval(1:lval),' not in range ( ',
     +                 cmin(1:lmin),' , ',cmax(1:lmax),' ]'     
      write(lunerr1,*) comment
      leng = imb_lenoc(usrnam3)
      if(leng.gt.0) then
        write(lunerr1,*) ' '
        write(lunerr1,*) subnam(1:7),
     +  ' was called by ',usrnam3(1:leng)
      endif

      stop
      end
      
C     =========================================================
      subroutine sqcIltlt(subnam,parnam,imin,ival,imax,comment)
C     =========================================================

C--   OK if imin .lt. ival .lt. imax

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      include 'qsnam3.inc'

      character*10  cmin, cval, cmax
      character*(*) parnam, comment
      character*(*) subnam

C--   First for the check
      if( (imin.lt.ival) .and. (ival.lt.imax) ) return
C--   Now for the error message
      leng = imb_lenoc(subnam)
      call smb_ItoCh(imin,cmin,lmin)
      call smb_ItoCh(ival,cval,lval)
      call smb_ItoCh(imax,cmax,lmax)      
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in ', subnam(1:leng), ' ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')
      write(lunerr1,*) parnam,' = ',cval(1:lval),' not in range ( ',
     +                 cmin(1:lmin),' , ',cmax(1:lmax),' )'     
      write(lunerr1,*) comment
      leng = imb_lenoc(usrnam3)
      if(leng.gt.0) then
        write(lunerr1,*) ' '
        write(lunerr1,*) subnam(1:7),
     +  ' was called by ',usrnam3(1:leng)
      endif

      stop
      end
      
C     =========================================================
      subroutine sqcDlele(subnam,parnam,dmin,dval,dmax,comment)
C     =========================================================

C--   OK if dmin .le. dval .le. dmax
C--   We use lqcRcomp to check equality within relative tolerance repsi6

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      include 'qsnam3.inc'
      include 'qpars6.inc'

      logical       lqcRcomp
      character*(*) parnam, comment
      character*(*) subnam

C--   First for the check
      if(lqcRcomp(dval,dmin,repsi6))            return
      if(lqcRcomp(dval,dmax,repsi6))            return
      if( (dmin.le.dval) .and. (dval.le.dmax) ) return
C--   Now for the error message
      leng = imb_lenoc(subnam)
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in ', subnam(1:leng), ' ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')
      write(lunerr1,'( 1X,A,'' = '',G11.4,'' not in range [ '',G11.4,
     +                  '' , '',G11.4,'' ]'')') parnam,dval,dmin,dmax
      write(lunerr1,*) comment
      leng = imb_lenoc(usrnam3)
      if(leng.gt.0) then
        write(lunerr1,*) ' '
        write(lunerr1,*) subnam(1:7),
     +  ' was called by ',usrnam3(1:leng)
      endif

      stop
      end
      
C     =========================================================
      subroutine sqcDlelt(subnam,parnam,dmin,dval,dmax,comment)
C     =========================================================

C--   OK if dmin .le. dval .lt. dmax
C--   We use lqcRcomp to check equality within relative tolerance repsi6

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      include 'qsnam3.inc'
      include 'qpars6.inc'

      logical       lqcRcomp
      character*(*) parnam, comment
      character*(*) subnam

C--   First for the check
      if(lqcRcomp(dval,dmin,repsi6))            return
      if(lqcRcomp(dval,dmax,repsi6))            goto 500
      if( (dmin.le.dval) .and. (dval.lt.dmax) ) return
C--   We dont pass the check
 500  continue
C--   Now for the error message
      leng = imb_lenoc(subnam)
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in ', subnam(1:leng), ' ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')
      write(lunerr1,'( 1X,A,'' = '',G11.4,'' not in range [ '',G11.4,
     +                  '' , '',G11.4,'' )'')') parnam,dval,dmin,dmax
      write(lunerr1,*) comment
      leng = imb_lenoc(usrnam3)
      if(leng.gt.0) then
        write(lunerr1,*) ' '
        write(lunerr1,*) subnam(1:7),
     +  ' was called by ',usrnam3(1:leng)
      endif

      stop
      end
      
C     =========================================================
      subroutine sqcDltle(subnam,parnam,dmin,dval,dmax,comment)
C     =========================================================

C--   OK if dmin .lt. dval .le. dmax
C--   We use lqcRcomp to check equality within relative tolerance repsi6

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      include 'qsnam3.inc'
      include 'qpars6.inc'

      logical       lqcRcomp
      character*(*) parnam, comment
      character*(*) subnam

C--   First for the check
      if(lqcRcomp(dval,dmin,repsi6))            goto 500
      if(lqcRcomp(dval,dmax,repsi6))            return
      if( (dmin.lt.dval) .and. (dval.le.dmax) ) return
C--   We dont pass the check
 500  continue
C--   Now for the error message
      leng = imb_lenoc(subnam)
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in ', subnam(1:leng), ' ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')
      write(lunerr1,'( 1X,A,'' = '',G11.4,'' not in range ( '',G11.4,
     +                  '' , '',G11.4,'' ]'')') parnam,dval,dmin,dmax
      write(lunerr1,*) comment
      leng = imb_lenoc(usrnam3)
      if(leng.gt.0) then
        write(lunerr1,*) ' '
        write(lunerr1,*) subnam(1:7),
     +  ' was called by ',usrnam3(1:leng)
      endif

      stop
      end

C     =========================================================
      subroutine sqcDltlt(subnam,parnam,dmin,dval,dmax,comment)
C     =========================================================

C--   OK if dmin .lt. dval .lt. dmax
C--   We use lqcRcomp to check equality within relative tolerance repsi6

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      include 'qsnam3.inc'
      include 'qpars6.inc'

      logical       lqcRcomp
      character*(*) parnam, comment
      character*(*) subnam

C--   First for the check
      if(lqcRcomp(dval,dmin,repsi6))            goto 500
      if(lqcRcomp(dval,dmax,repsi6))            goto 500
      if( (dmin.lt.dval) .and. (dval.lt.dmax) ) return
C--   We dont pass the check
 500  continue
C--   Now for the error message
      leng = imb_lenoc(subnam)
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in ', subnam(1:leng), ' ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')
      write(lunerr1,'( 1X,A,'' = '',G12.4,'' not in range ( '',G12.4,
     +                  '' , '',G12.4,'' )'')') parnam,dval,dmin,dmax
      write(lunerr1,*) comment
      leng = imb_lenoc(usrnam3)
      if(leng.gt.0) then
        write(lunerr1,*) ' '
        write(lunerr1,*) subnam(1:7),
     +  ' was called by ',usrnam3(1:leng)
      endif

      stop
      end

C     ====================================
      subroutine sqcErrMsg(subnam,message)
C     ====================================

C--   Print a message

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      include 'qsnam3.inc'

      character*(*) message
      character*(*) subnam

C--   Now for the error message
      leng = imb_lenoc(subnam)
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in ', subnam(1:leng), ' ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')
      write(lunerr1,*) message
C--   Extra message      
      leng = imb_lenoc(usrnam3)
      if(leng.gt.0) then
        write(lunerr1,*) ' '
        write(lunerr1,*) subnam(1:7),
     +  ' was called by ',usrnam3(1:leng)
      endif

      stop
      end
      
C     ==============================================
      subroutine sqcIdEmsg(subnam,parnam,idval,ierr)
C     ==============================================

C--   Print messages corresponding to iqcCheckId error code
C--
C--   subnam   (in)  : subroputine name
C--   parnam   (in)  : identifier name
C--   idval    (in)  : identifier value
C--   ierr     (in)  : iqcCheckId error code

      implicit double precision (a-h,o-z)

      include 'qluns1.inc'
      include 'qsnam3.inc'

      character*(*) subnam
      character*(*) parnam
      
      character*10 number
      character*41 emsg(-1:4)
C--              '12345678901234567890123456789012345678901'
      data emsg /'Empty weight table',                            ! -1
     +           '  ',                                            !  0
     +           'No tables booked please call BOOKTAB',          !  1
     +           'Table type not 1, 2, 3 or 4',                   !  2
     +           'No tables booked of this type',                 !  3
     +           'Refers to non-existent weight table' /          !  4

      if(ierr.eq. 0) return
      if(ierr.lt.-1) return
      if(ierr.gt. 6) return
      
C--   Now for the error message
      leng = imb_lenoc(subnam)
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in ', subnam(1:leng), ' ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')
      lpar = imb_lenoc(parnam)
      call smb_ItoCh(idval,number,lnum)
      leng = imb_lenoc(emsg(ierr))
      write(lunerr1,'(A,'' = '',A,2X,A)') parnam(1:lpar),
     +               number(1:lnum),emsg(ierr)(1:leng)
C--   Extra message      
      leng = imb_lenoc(usrnam3)            
      if(leng.gt.0) then
        write(lunerr1,*) ' '
        write(lunerr1,*) subnam(1:7),
     +  ' was called by ',usrnam3(1:leng)
      endif

      stop
      end
      
C     ==============================================      
      subroutine sqcCutMsg(subnam,x,q,ifail,noextra)
C     ==============================================

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qgrid2.inc'
      include 'qsnam3.inc'
      
      character*(*) subnam
      
C--   Now for the error message      
      leng = imb_lenoc(subnam)
      write(lunerr1,'(/1X,70(''-''))')
      write(lunerr1,*) 'Error in ', subnam(1:leng), ' ---> STOP'
      write(lunerr1,'( 1X,70(''-''))')
      
      if(ifail.eq.0 .or. ifail.eq.1) then
        write(lunerr1,
     +   '('' X = '',1PE11.3,'' fails range ['',1PE11.3,
     +     '', 1.0 ]'')')  x,xminc2
      elseif(ifail.eq.2) then
        write(lunerr1,
     +   '('' Mu2 = '',1PE11.3,'' fails qmin cut'',1PE11.3)') q,qminc2
      elseif(ifail.eq.3) then
        write(lunerr1,
     +   '('' Mu2 = '',1PE11.3,'' fails qmax cut'',1PE11.3)') q,qmaxc2 
      else
        write(lunerr1,'(''Unknown ifail'')')
      endif
      
C--   Extra message      
      leng = imb_lenoc(usrnam3)            
      if(leng.gt.0 .and. noextra.ne.1) then
        write(lunerr1,*) ' '
        write(lunerr1,*) subnam(1:7),
     +  ' was called by ',usrnam3(1:leng)
      endif

      stop
      end
                                

