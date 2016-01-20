
C--   This is the file qcdstore.f containing handling of the pdf store

C--   subroutine sqcPdfIni(iwmin,iwmax,idmin,idmax,lpdf)
C--   subroutine sqcPdfBuf(ntab,iwmin,iwmax,kanul)
C--   integer function iqcPdfijkl(i,j,k,iset)
C--
C--   subroutine sqcPreSet(jset,id1,val)
C--   subroutine sqcPSetjj(jset,id,it,val)
C--   subroutine sqcPCopjj(jset,id1,j1,id2,j2)
C--   subroutine sqcPdfCop(jset,id1,id2)
C--   subroutine sqcT1toT2(jset,idin,idout,iy1,iy2,iz1,iz2)
C--
C--   subroutine sqcG0toGi(jset,idg0,idgi,ig,nyg,iz)
C--   subroutine sqcGitoG0(jset,idgi,ig,idg0)
C--   subroutine sqcGiFtoA(jset,id1,id2,nyg,iw1,iw2)
C--   subroutine sqcGiAtoF(jset,id1,id2,nyg,iw1,iw2)
C--
C--   subroutine sqcGiLtoQ(jset,id1,id2,nyg,iz1,iz2)
C--   subroutine sqcGiQtoL(jset,id1,id2,nyg,iz1,iz2)
C--
C--   subroutine sqcAllAtoF(jset,id1,id2)
C--   subroutine sqcGetSplA(jset,id,iy,iz,igout,iyout,aout)
C--   subroutine sqcAitoF0(jset,idgi,ig,idg0)


C===================================================================
C==   Partition and address the pdf store ==========================
C===================================================================

C     ==================================================
      subroutine sqcPdfIni(iwmin,iwmax,idmin,idmax,lpdf)
C     ==================================================

C--   Partition stor7(nwf0) into  tables pdf(0:nyy2,0:ntt2+7,-4:idmax)
C--   id = -4, -3, -2, -1 are reserved as working arrays
C--   id =  0,...,idmax are for the user
C--
C--   Input:  iwmin = address first word of first pdf 
C--   Output: iwmax = address last  word of last  pdf < 0 no space
C--           idmin = first user identifier (always idmin = 0)
C--           idmax = last user identifier  
C--           lpdf  = number of words used by each pdf (nyy2+1)*(ntt2+8)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

C--   Size of pdf array and # of tables
      lpdf  = (nstory2+1)*(ntt2+8)   !number of words per table
      leng7 = lpdf                   !number of words per table
      npdf7 = 13 + nwrk0             !number of tables
      idmin = 0                      !gluon identifier
      idmax = npdf7-nwrk0-1          !largest identifier
      idmx7 = idmax
C--   Partition the array stor7(1:nwf0) --> 
C--                       pdf( iy=0:nstory2, it=0:ntt2+7, id=-nwrk0:idmax )
C--   Warning: if table dimensions change then make also the change in
C--            subroutine sqcBookTab in file qcdustf.f
      ipdf7(1) = 0
      jpdf7(1) = nstory2
      ipdf7(2) = 0
      jpdf7(2) = ntt2+7
      ipdf7(3) = -nwrk0
      jpdf7(3) = idmax
C--   smb_bkmat comes from the mbclib package
      call smb_bkmat(ipdf7,jpdf7,kpdf7,3,iwmin,iwmax)
C--   iwmax is the max address --> check again that it fits in nwf0
      if(iwmax.gt.nwf0) iwmax = -iwmax

      return
      end
      
C     ============================================
      subroutine sqcPdfBuf(ntab,iwmin,iwmax,kanul)
C     ============================================

C--   Book ntab pdf tables in the store
C--
C--   ntab   (in)  : number of tables to be booked
C--   iwmin  (in)  : address of first word of table  
C--   iwmax  (out) : address of last  word of table (<0 no space)
C--   kanul  (out) : k(0) address parameter


      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpdfs7.inc'
      
C--   Size of pdf array and # of tables
      lpdf  = (nstory2+1)*(ntt2+8)   !number of words per table
      idmin = 1                      !first identifier
      idmax = ntab                   !last  identifier
C--   Partition the array stor7(1:nwf0) --> 
C--                       pdf( iy=0:nstory2, it=0:ntt2+7, id=-nwrk0:idmax )
C--   Warning: if table dimensions change then make also the change in
C--            subroutine sqcBookTab in file qcdustf.f
      ipdf7(1) = 0
      jpdf7(1) = nstory2
      ipdf7(2) = 0
      jpdf7(2) = ntt2+7
      ipdf7(3) = idmin
      jpdf7(3) = idmax
C--   smb_bkmat comes from the mbclib package
      call smb_bkmat(ipdf7,jpdf7,kpdf7,3,iwmin,iwmax)
      kanul  = kpdf7(0)
C--   iwmax is the max address --> check that it fits in nwf0
      if(iwmax.gt.nwf0) iwmax = -iwmax

      return
      end

C     =======================================
      integer function iqcPdfijkl(i,j,k,iset)
C     =======================================

C--   Find the address (ia) in stor7(1:nwf0).
C--   The address arithmetic is set-up in s/r sqcPdfIni.
C--   The constants kpdf7 are stored in /qpdfs7/
C--
C--   i [0,nstory2] = first  index of pdf array (adresses iy)
C--   j [0,ntt2+7]  = second index (addresses it)
C--   k [-4,idmx7]  = third  index (addresses id)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      iqcPdfijkl = knul7(iset)+kpdf7(1)*i+kpdf7(2)*j+kpdf7(3)*k

      return
      end      

C===================================================================
C==   Basic operations on the store ================================
C===================================================================

C     ==================================
      subroutine sqcPreSet(jset,id1,val)
C     ==================================

C--   Set all entries of id1 to val

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

C--   Calculate base address
      iy1 = ipdf7(1)
      iy2 = jpdf7(1)
      it1 = ipdf7(2)
      it2 = jpdf7(2)
      ia1 = iqcPdfIjkl(iy1,it1,id1,jset)-1
C--   Check leng7
      nwd = (iy2-iy1+1)*(it2-it1+1)
      if(nwd.ne.leng7) stop 'sqcPreset: nwd not equal to leng7'      
C--   Set leng7 words to val; leng7 is the length of each table (in words)
      do i = 1,leng7
        ia1        = ia1+1
        stor7(ia1) = val
      enddo

      return
      end
      
C     ====================================
      subroutine sqcPSetjj(jset,id,it,val)
C     ====================================

C--   Set all entries of id at bin it to val

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

C--   Calculate base address
      iy1 = ipdf7(1)
      iy2 = jpdf7(1)
      ia1 = iqcPdfIjkl(iy1,it,id,jset)-1
C--   Set value for all iy
      do i = iy1,iy2
        ia1        = ia1+1
        stor7(ia1) = val
      enddo

      return
      end      

C     ========================================
      subroutine sqcPCopjj(jset,id1,j1,id2,j2)
C     ========================================

C--   Copy column j1 of table id1 to column j2 of table id2

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

C--   Calculate base addresses
      iy1 = ipdf7(1)
      iy2 = jpdf7(1)
      it1 = ipdf7(2)
      it2 = jpdf7(2)
      ia1 = iqcPdfIjkl(iy1,j1,id1,jset)-1
      ia2 = iqcPdfIjkl(iy1,j2,id2,jset)-1
C--   Copy column
      do i = iy1,iy2
        ia1        = ia1+1
        ia2        = ia2+1
        stor7(ia2) = stor7(ia1)
      enddo

      return
      end

C     ==================================
      subroutine sqcPdfCop(jset,id1,id2)
C     ==================================

C--   Copy table id1 to table id2

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      if(id1.eq.id2) return
C--   Calculate base addresses
      iy1 = ipdf7(1)
      iy2 = jpdf7(1)
      it1 = ipdf7(2)
      it2 = jpdf7(2)
      ia1 = iqcPdfIjkl(iy1,it1,id1,jset)-1
      ia2 = iqcPdfIjkl(iy1,it1,id2,jset)-1
C--   Check table size
      nwd = (iy2-iy1+1)*(it2-it1+1)
      if(nwd.ne.leng7) stop 'sqcPdfCop: nwd not equal to leng7'      
C--   Copy leng7 words; leng7 is the length of each table (in words)
      do i = 1,leng7
        ia1        = ia1+1
        ia2        = ia2+1
        stor7(ia2) = stor7(ia1)
      enddo

      return
      end
      
C     =====================================================
      subroutine sqcT1toT2(jset,idin,idout,iy1,iy2,iz1,iz2)
C     =====================================================

C--   Copy given range of table idin to table idout
C--
C--   jset     (in)   Pdf set identifier
C--   idin     (in)   Table 1 identifier
C--   idout    (in)   Table 2 identifier
C--   iy1,iy2  (in)   iy range to copy
C--   iz1,iz2  (in)   iz range to copy

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      if(idin.eq.idout) return
C--   Calculate increments
      ia0  = iqcPdfIjkl(1,1,idin,jset)
      incy = iqcPdfIjkl(2,1,idin,jset) - ia0 
      incz = iqcPdfIjkl(1,2,idin,jset) - ia0 
      iain = iqcPdfIjkl(iy1,iz1,idin ,jset)
      idif = iqcPdfIjkl(iy1,iz1,idout,jset)-iain
      iaz  = iain-incz
      do iz = iz1,iz2
        iaz = iaz + incz
        iay = iaz - incy
        do iy = iy1,iy2 
          iay        = iay + incy
          ia2        = iay + idif
          stor7(ia2) = stor7(iay)
        enddo  
      enddo

      return
      end
      

C===================================================================
C==   Transformations between interpolation and evolution grids ====
C===================================================================

C     ==============================================
      subroutine sqcG0toGi(jset,idg0,idgi,ig,nyg,iz)
C     ==============================================

C--   Copy t-bin iz from interpolation grid G0 to subgrid Gi.
C--
C--   (in)  idg0  identifier of the nonequidistant interpolation grid G0
C--   (in)  idgi  identifier of the equidistant evolution subgrid
C--   (in)  ig    index (i) = 1,...,nyg2  of the subgrid
C--   (in)  nyg   Upper yloop index in subgrid
C--   (in)  iz    z-bin in G0   to be copied 
C--
C--   For the y-grid index parameters see sqcGyMake in qcdgrd.f

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

      iai = iqcPdfIjkl(1,iz,idgi,jset)-1
      ia0 = iqcPdfIjkl(1,iz,idg0,jset)-1
C--   Loop over ybins in Grid i
      do iyi = 1,nyg
        iy0            = iy0fiyg2(iyi,ig)
        stor7(iai+iyi) = stor7(ia0+iy0)
      enddo
     
      return
      end

C     =======================================
      subroutine sqcGitoG0(jset,idgi,ig,idg0)
C     =======================================

C--   Copy F coefficients from Gi to interpolation grid G0
C--   The bin it=0 is not copied
C--
C--   (in)  idgi  identifier of the equidistant evolution subgrid
C--   (in)  ig    index (i) = 1,...,nyg2  of the subgrid
C--   (in)  idg0  identifier of the nonequidistant interpolation grid G0
C--
C--   For the y-grid index parameters see sqcGyMake in qcdgrd.f

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

C--   Loop over z-grid
      do iz = 1,nzz2
C--     Base address in grid G0 and Gi
        ia0 = iqcPdfIjkl(1,iz,idg0,jset)-1
        iai = iqcPdfIjkl(1,iz,idgi,jset)-1
C--     Loop over subgrid points to copy
        do iyi = jymi2(ig),nyy2(ig)
          iy0            = iy0fiyg2(iyi,ig)
          stor7(ia0+iy0) = stor7(iai+iyi)
        enddo                          
      enddo

      return
      end

C     ==============================================
      subroutine sqcGiFtoA(jset,id1,id2,nyg,iz1,iz2)
C     ==============================================

C--   Convert t-bins iz1 to iz2 from F to A and store result in id2
C--
C--   id1   Input  subgrid table
C--   id2   Output subgrid table
C--   nyg   Upper limit yloop
C--   iz1   First t-bin to convert 
C--   iz2   Last  t-bin to convert
C--   
C--   It is allowed to have id1 = id2 (in-place conversion)
C--   To convert all t-bins, set iz1 (iz2) very large negative (positive)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

C--   Index ranges
      iy1 = 1
      iy2 = nyg
      it1 = max(ipdf7(2),iz1)
      it2 = min(jpdf7(2),iz2)
C--   Calculate base addresses
      inc = iqcPdfIjkl(iy1,it1+1,id1,jset)-
     +      iqcPdfIjkl(iy1,it1  ,id1,jset)
      ia1 = iqcPdfIjkl(iy1,it1,id1,jset)-inc
      ia2 = iqcPdfIjkl(iy1,it1,id2,jset)-inc
C--   Convert columns f --> a by solving f = Sa
      do j = it1,it2
        ia1 = ia1+inc
        ia2 = ia2+inc
        call sqcNSeqs(smaty2,nmaty2,stor7(ia2),stor7(ia1),iy2)
      enddo

      return
      end

C     =============================================
      subroutine sqcGiAtoF(jset,id1,id2,nyg,iz1,iz2)
C     =============================================

C--   Convert t-bins iz1 to iz2 from A to F and store result in id2
C--
C--   id1   Input  subgrid table
C--   id2   Output subgrid table
C--   nyg   Upper limit yloop
C--   iz1   First t-bin to convert 
C--   iz2   Last  t-bin to convert
C--   
C--   It is allowed to have id1 = id2 (in-place conversion via buffer)
C--   To convert all t-bins, set iz1 (iz2) very large negative (positive)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

C--   Index ranges
      iy1 = 1
      iy2 = nyg
      it1 = max(ipdf7(2),iz1)
      it2 = min(jpdf7(2),iz2)
C--   Calculate base addresses
      inc = iqcPdfIjkl(iy1,it1+1,id1,jset)-
     +      iqcPdfIjkl(iy1,it1  ,id1,jset)
      ia1 = iqcPdfIjkl(iy1,it1,id1,jset)-inc
      ia2 = iqcPdfIjkl(iy1,it1,id2,jset)-inc
C--   Convert columns a --> f by multiplying f = Sa
      do j = it1,it2
        ia1 = ia1+inc
        ia2 = ia2+inc
        call sqcNSmult(smaty2,nmaty2,stor7(ia1),bufy7,iy2)
        do i = 1,iy2
          stor7(ia2-1+i) = bufy7(i)
        enddo
      enddo

      return
      end

C     ==============================================
      subroutine sqcGiLtoQ(jset,id1,id2,nyg,iz1,iz2)
C     ==============================================

C--   Switch to quadratic interpolation and convert A coefficients
C--   of bins iz1 to iz2 from lin to quad 
C--
C--   id1   Input  subgrid table
C--   id2   Output subgrid table
C--   nyg   Upper limit yloop
C--   iz1   First t-bin to convert 
C--   iz2   Last  t-bin to convert
C--   
C--   It is allowed to have id1 = id2 (in-place conversion)
C--   To convert all t-bins, set iz1 (iz2) very large negative (positive)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'

C--   A coefficients are already quadratic
      if(ioy2.eq.3) then
*mb        call sqcPdfCop(jset,id1,id2)
        call sqcT1toT2(jset,id1,id2,1,nyg,iz1,iz2)
        return
      endif  
C--   Convert lin A to F
      call sqcGiAtoF(jset,id1,id2,nyg,iz1,iz2)
C--   Set quad interpolation
      ioy2 = 3
      call sqcGryMat(ioy2)
C--   Convert F to quad A
      call sqcGiFtoA(jset,id2,id2,nyg,iz1,iz2)
      
      return
      end

C     ==============================================
      subroutine sqcGiQtoL(jset,id1,id2,nyg,iz1,iz2)
C     ==============================================

C--   Switch to linear interpolation and convert A coefficients
C--   of bins iz1 to iz2 from quad to lin 
C--
C--   id1   Input  subgrid table
C--   id2   Output subgrid table
C--   nyg   Upper limit yloop
C--   iz1   First t-bin to convert 
C--   iz2   Last  t-bin to convert
C--   
C--   It is allowed to have id1 = id2 (in-place conversion)
C--   To convert all t-bins, set iz1 (iz2) very large negative (positive)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'

C--   A coefficients are already linear
      if(ioy2.eq.2) then
        call sqcPdfCop(jset,id1,id2)
        return
      endif  
C--   Convert lin A to F
      call sqcGiAtoF(jset,id1,id2,nyg,iz1,iz2)
C--   Set lin interpolation
      ioy2 = 2
      call sqcGryMat(ioy2)
C--   Convert F to quad A
      call sqcGiFtoA(jset,id1,id2,nyg,iz1,iz2)
      
      return
      end

C===================================================================
C==   Spline coefficient transformations ===========================
C===================================================================

C     ===================================
      subroutine sqcAllAtoF(jset,id1,id2)
C     ===================================

C--   f = Sa:  a in id1 --> f in id2
C--   It is allowed to have id1 = id2 (inplace conversion via buffer)
C--   This routine works on an equidistant evolution grid
      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      include 'qgrid2.inc'

C--   Calculate base addresses
      iy1 = 1             !start at iy=1 not iy=0
      iy2 = jpdf7(1)
      it1 = ipdf7(2)
      it2 = jpdf7(2)
      inc = iqcPdfIjkl(iy1,it1+1,id1,jset) -
     +      iqcPdfIjkl(iy1,it1,id1,jset)
      ia1 = iqcPdfIjkl(iy1,it1,id1,jset)-inc
      ia2 = iqcPdfIjkl(iy1,it1,id2,jset)-inc
C--   Convert columns a --> f = Sa
      do j = it1,it2
        ia1 = ia1+inc
        ia2 = ia2+inc
        call sqcNSmult(smaty2,nmaty2,stor7(ia1),bufy7,nyy2(0))
        do i = 1,nyy2(0)
          stor7(ia2-1+i) = bufy7(i)
        enddo
      enddo

      return
      end

C     ================================================      
      subroutine sqcGetSplA(jset,id,iy,iz,ig,iyg,aout)
C     ================================================

C--   id     (in)   pdf table index
C--   iy     (in)   y-grid index in main grid G0
C--   iz     (in)   z-grid index in main grid G0
C--   ig     (out)  subgrid index
C--   iyg    (out)  y-grid index in subgrid Gi
C--   aout   (out)  array aout(1,...,iyg) of spline coefficients
C--                 defined on the subgrid Gi 

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'
      
      dimension aout(*),buf(mxx0)
      
C--   Find subgrid
      ig = 1
      do i = 2,nyg2
        if(iy.gt.iyma2(i-1)) ig = i
      enddo
C--   Find y-grid index in subgrid (basically iy = y/del)
      iyg = iqcIyfrmY(ygrid2(iy),dely2(ig),nyy2(ig))
C--   Cross-check
      if(iy.ne.iy0fiyg2(iyg,ig)) then
       stop 'sqcGetSplA: problem y index in subgrid'
      endif
C--   Base address
      ia0 = iqcPdfIjkl(1,iz,id,jset)-1
C--   Copy pdf values to buffer
      do jy = 1,iyg
        iy0     = iy0fiyg2(jy,ig)
        buf(jy) = stor7(ia0+iy0)
      enddo
C--   Now convert pdf values to spline coefficients
      call sqcNSeqs(smaty2,nmaty2,aout,buf,iyg)
      
      return
      end

C     ===================================================
      subroutine sqcAitoF0(jset,idgi,ig,nyg,iz1,iz2,idg0)
C     ===================================================

C--   Transform A coefficients in Gi to F values in G0
C--   The bin it=0 is not copied
C--
C--   (in)  jset  pdf set identifier
C--   (in)  idgi  identifier of the equidistant evolution subgrid
C--   (in)  ig    index of the subgrid
C--   (in)  nyg   upper loop index in subgrid
C--   (in)  iz1   lower z index
C--   (in)  iz2   upper z index
C--   (in)  idg0  identifier of the interpolation grid G0
C--
C--   For the y-grid index parameters see sqcGyMake in qcdgrd.f

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
      include 'qpdfs7.inc'

*C--   Loop over z-grid
*      do iz = 1,nzz2
*C--     Base address in grid G0 and Gi
*        ia0 = iqcPdfIjkl(1,iz,idg0,jset)-1
*        iai = iqcPdfIjkl(1,iz,idgi,jset)
*        call sqcNSmult(smaty2,nmaty2,stor7(iai),bufy7,nyy2(ig))
*        do iyi = jymi2(ig),nyy2(ig)
*          iy0            = iy0fiyg2(iyi,ig)
*          stor7(ia0+iy0) = bufy7(iyi)
*        enddo  
*      enddo
      
C--   Loop over z-grid
      do iz = iz1,iz2
C--     Base address in grid G0 and Gi
        ia0 = iqcPdfIjkl(1,iz,idg0,jset)-1
        iai = iqcPdfIjkl(1,iz,idgi,jset)
        call sqcNSmult(smaty2,nmaty2,stor7(iai),bufy7,nyg)
        do iyi = jymi2(ig),nyg
          iy0            = iy0fiyg2(iyi,ig)
          stor7(ia0+iy0) = bufy7(iyi)
        enddo  
      enddo      
      
      return
      end

