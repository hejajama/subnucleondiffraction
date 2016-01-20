
C     ===========================
      subroutine sqcDebug(srname)
C     ===========================

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      
      character*(*) srname
      
      leng = imb_lenoc(srname)
      
      write(6,'(/A,'' ifirst7, iset7  ='',3I6)') srname,
     +          ifirst7,iset7,kpdf7(0)
      write(6,'(6X,'' ifst7(iset,ioy) ='',8I6)') ifst7
      write(6,'(6X,'' ilst7(iset,ioy) ='',8I6)') ilst7
      write(6,'(6X,'' mord7(iset)     ='',4I6)') mxord7
      write(6,'(6X,'' knul7(iset)     ='',4I6)') knul7
      
      return
      end
  