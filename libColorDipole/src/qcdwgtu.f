
C--   file qcdwgtu.f containing (user) weight routines
C--
C--   subroutine sqcBookTab(w,nw,itypes,nwords,ierr)
C--   subroutine sqcSetKey(keyin,keyout)
C--   logical function lqcSjekey(key1,key2)
C--   subroutine sqcDumpTab(w,lun,key,ierr)
C--   subroutine sqcReadTab(w,nw,lun,key,nwords,ierr)

C--   integer function iqcWaddr(w,i,j,k,l,m)
C--   integer function iqcWCadr(w,i,j,k,l,m)
C--   integer function iqcW1ijk(w,i,j,k)
C--   integer function iqcW2ijkl(w,i,j,k,l)
C--   integer function iqcW3ijkl(w,i,j,k,l)
C--   integer function iqcW4ijklm(w,i,j,k,l,m)

C--   subroutine sqcTabLims(w,ityp,imin,imax,ierr)
C--   subroutine sqcSetNoWt(w,ityp)
C--   integer function iqcCheckId(w,id)
C--   subroutine sqcChekPij(idin,idout,jset,ierr)

C--   subroutine sqcUweitA(w,id,ioy,afun,achi,ierr)
C--   subroutine sqcUweitB(w,id,ioy,bfun,achi,idel,ierr)
C--   subroutine sqcUwgtRS(w,id,ioy,rfun,sfun,achi,idel,ierr)
C--   subroutine sqcUweitD(w,id,ioy,dfun,achi,ierr)
C--   subroutine sqcUweitX(w,id,ioy,ierr)

C--   subroutine sqcScaleWt(w,c,id)
C--   subroutine sqcCopyWt(w1,id1,w2,id2,iadd)
C--   subroutine sqcWcrossW(wa,ida,wb,idb,wc,idc,idt1,idt2,iadd)
C--   subroutine sqcWtimesF(w,fun,id1,id2,iadd)

C--   double precision function dqcUIgauss(pfun,ti,nf,a,b)
C--   double precision function 
C--  +   dqcUAgauss(idk,afun,yi,ti,nf,a,b,del)
C--   double precision function 
C--  +   dqcUBgauss(idk,bfun,yi,ti,nf,a,b,del)
C--   double precision function 
C--  +   dqcURSgaus(idk,rfun,sfun,yi,ti,nf,a,b,del)

C==   ==============================================================
C==   Store partition, dump and read ===============================
C==   ==============================================================

C     ==============================================
      subroutine sqcBookTab(w,nw,itypes,nwords,ierr)
C     ==============================================

C--   Partition store in tables of different types
C--
C--   Input:  w         array dimensioned nw in the calling routine
C--           nw        number of words in w
C--           itypes(4) number of tables of type (i)
C--   Output: nwords    number of words needed < 0 if not enough space
C--           ierr      0 all OK
C--                     1 no valid itypes encountered
C--                     2 not enough space               
C--
C--   Layout of the store:
C--
C--   first word           :   partition code 123456
C--   second word          :   number of words used
C--   third word           :   number of info words (miw0)
C--   next  miw0 words     :   info words
C--   next  4 words        :   table addresses of each set (0 = absent)
C--   next  16 words       :   imin,imax,karr of first 5-dim set 
C--   next  words          :   tables of the first set
C--   next  16 words       :   imin,imax,karr of second 5-dim set
C--   next  words          :   tables of the second set
C--   etc.
C--
C--   Indexing of tables: always 5 indices with one or more dummies
C--                  1        2        3       4        5  
C--   type 1 : w( ix[1,nx] iq[1,1]  nf[3,3] ig[1,ng] id[101,199] )
C--   type 2 : w( ix[1,nx] iq[1,1]  nf[3,6] ig[1,ng] id[201,299] )
C--   type 3 : w( ix[1,nx] iq[1,nq] nf[3,3] ig[1,ng] id[301,399] )
C--   type 4 : w( ix[1,nx] iq[1,nq] nf[3,6] ig[1,ng] id[401,499] )

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'

      dimension w(*), itypes(*)
      dimension imin(5),imax(5),karr(0:5)

C--   Initialize store
      do i = 1,nw
        w(i) = 0.D0
      enddo
      do i = 1,5
        imin(i) = 0
        imax(i) = 0
        karr(i) = 0
      enddo
      karr(0) = 0
      ifirstw = 7 + miw0 + 1
      nwheadr = 16
      ierr    = 0
      jwrite  = 0
C--   Loop over table types
      do ityp = 1,4
        iwrite = 0
        if(ityp.eq.1 .and. 
     +    itypes(1).ge.1 .and. itypes(1).le.99) then
          imin(1) =  1
          imax(1) =  nyy2(0)
          imin(2) =  1
          imax(2) =  1
          imin(3) =  3
          imax(3) =  3
          imin(4) =  1
          imax(4) =  nyg2
          imin(5) =  101
          imax(5) =  100+itypes(1)
          iwrite  =  1
          jwrite  =  1
        elseif(ityp.eq.2 .and.
     +    itypes(2).ge.1 .and. itypes(2).le.99) then
          imin(1) =  1
          imax(1) =  nyy2(0)
          imin(2) =  1
          imax(2) =  1
          imin(3) =  3
          imax(3) =  6
          imin(4) =  1
          imax(4) =  nyg2
          imin(5) =  201
          imax(5) =  200+itypes(2)
          iwrite  =  1
          jwrite  =  1
        elseif(ityp.eq.3 .and.
     +    itypes(3).ge.1 .and. itypes(3).le.99) then
          imin(1) =  1
          imax(1) =  nyy2(0)
          imin(2) =  1
          imax(2) =  ntt2
          imin(3) =  3
          imax(3) =  3
          imin(4) =  1
          imax(4) =  nyg2
          imin(5) =  301
          imax(5) =  300+itypes(3)
          iwrite  =  1
          jwrite  =  1
        elseif(ityp.eq.4 .and.
     +    itypes(4).ge.1 .and. itypes(4).le.99) then
          imin(1) =  1
          imax(1) =  nyy2(0)
          imin(2) =  1
          imax(2) =  ntt2
          imin(3) =  3
          imax(3) =  6
          imin(4) =  1
          imax(4) =  nyg2
          imin(5) =  401
          imax(5) =  400+itypes(4)
          iwrite  =  1
          jwrite  =  1
        endif
        if(iwrite.eq.1) then
C--       Table partition
          call smb_bkmat(imin,imax,karr,5,ifirstw+nwheadr,ilastw)
          if(ilastw.le.nw) then
C--         Store partition definition (if enough space)
            iw = ifirstw-1
            do i = 1,5
              iw = iw+1
              w(iw) = dble(imin(i))
              iw = iw+1
              w(iw) = dble(imax(i))
            enddo
            do i = 0,5
              iw = iw+1
              w(iw) = dble(karr(i))
            enddo
C--         Store base address
            w(ityp+3+miw0) = dble(ifirstw)
          endif
C--       Start of next set of tables
          ifirstw = ilastw+1
        endif
C--   End of loop over table types
      enddo

C--   No table type found
      if(jwrite.eq.0) then
         ierr   =  1
         nwords = -1
C--   Enough words? 
      elseif(ilastw.le.nw) then
        nwords = ilastw
        w(1)   = dble(123456)
        w(2)   = dble(nwords)
        w(3)   = miw0
C--     Flag tables as unfilled
        call sqcSetNoWt(w,1)
        call sqcSetNoWt(w,2)
        call sqcSetNoWt(w,3)
        call sqcSetNoWt(w,4)
      else
        nwords = -ilastw
        ierr   = 2
      endif

      return
      end
      
C     ==================================
      subroutine sqcSetKey(keyin,keyout)
C     ==================================

C--   Left adjust, truncate to 50 chars and convert to upper case
C--
C--   ierr   (out)  0 = OK
C--                 1 = key truncated

      implicit double precision (a-h,o-z)
      
      character*(*) keyin
      character*50  keyout
      
      call smb_cfill(' ',keyout)
      i1   = imb_frstc(keyin)
      i2   = imb_lenoc(keyin)
      if(i1.eq.i2) return           !empty string
      leng = min(i2-i1+1,50)
      keyout = keyin(i1:i1+leng-1)
      call smb_cltou(keyout)
      
      return
      end
      
C     =====================================      
      logical function lqcSjekey(key1,key2)
C     =====================================

C--   True if two keys match

      implicit double precision(a-h,o-z)
      
      character*(*) key1,key2
      character*50  kkk1,kkk2              

      call sqcSetKey(key1,kkk1)
      call sqcsetKey(key2,kkk2)
      
      lqcSjekey = .false.
      if(kkk1.eq.kkk2) lqcSjekey = .true.
      
      return
      end
      
C     =====================================
      subroutine sqcDumpTab(w,lun,key,ierr)
C     =====================================

C--   Dump store on logical unit number lun
C--
C--   ierr = 0  : all OK
C--          1  : write error

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qvers1.inc'
      include 'qgrid2.inc'

      dimension w(*)
      character*(*) key
      character*50  keyout

C--   Initialize
      ierr = 0
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
      nwords = int(anint(w(2)))
      write(lun,err=500) nwords
      write(lun,err=500) (w(i),i=1,nwords)

      return

 500  continue
C--   Write error
      ierr = 1
      return

      end

C     ===============================================
      subroutine sqcReadTab(w,nw,lun,key,nwords,ierr)
C     ===============================================

C--   Read store from logical unit number lun.
C--
C--   Ierr = 0  : all OK
C--          1  : read error
C--          2  : problem with QCDNUM version
C--          3  : key mismatch
C--          4  : x-mu2 grid not the same
C--          5  : store too small

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qvers1.inc'
      include 'qgrid2.inc'

      character*10  cversr
      character*8   cdater
      character*50  keyred
      character*(*) key
      logical lqcSjekey
      dimension nyyr(0:mxg0),delyr(0:mxg0)
      dimension tgridr(mqq0)

      dimension w(*)

C--   Initialize
      nwords = 0
      ierr   = 0
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
*mb          write(6,*) 'i,tgrid2,tgridr = ',i,tgrid2(i),tgridr(i)
          return
        endif
      enddo
      
C--   Now read the store
      read(lun,err=500,end=500) nwords
      if(nwords.gt.nw) then
C--     Store too small
        ierr = 5
        return
      endif
      read(lun,err=500,end=500) (w(i),i=1,nwords)

      return

  500 continue
C--   Read error
      ierr = 1
      return

      end

C==   ==============================================================
C==   Indexing =====================================================
C==   ==============================================================

C     ======================================
      integer function iqcWaddr(w,i,j,k,l,m)
C     ======================================

C--   Find the address (ia) in store w: acts as a do-nothing if store
C--   not partitioned or if no tables of type ityp exist
C--
C--   Meaning of indices:  i   j   k   l   m
C--                       ix  iq  nf  ig  id

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension w(*)

      iqcWaddr = 0
      if(anint(w(1)).ne.123456)    return
      ityp = m/100
      if(ityp.le.0 .or. ityp.ge.5) return
      iw = int(anint(w(ityp+3+miw0)))
      if(iw.eq.0)                  return
      iqcWaddr = int(anint(w(iw+10)) + 
     +           anint(w(iw+11))*i + anint(w(iw+12))*j +
     +           anint(w(iw+13))*k + anint(w(iw+14))*l +
     +           anint(w(iw+15))*m)
     
      return
      end 

C     ======================================
      integer function iqcWCadr(w,i,j,k,l,m)
C     ======================================

C--   As iqcWaddr but with array boundary check 

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension w(*)

      iqcWaddr = 0
      if(anint(w(1)).ne.123456)    stop
     +  'iqcWCadr: store not partitioned'
      ityp = m/100
      if(ityp.le.0 .or. ityp.ge.5) stop
     +  'iqcWCadr: impossible table type'
      iw = int(anint(w(ityp+3+miw0)))
      if(iw.eq.0)                  stop
     +  'iqcWCadr: table type not in store'
      if(i.lt.anint(w(iw   )) .or. i.gt.anint(w(iw+1 ))) stop
     +  'iqcWCadr: index 1 (i) out of range'
      if(j.lt.anint(w(iw+2 )) .or. j.gt.anint(w(iw+3 ))) stop
     +  'iqcWCadr: index 2 (j) out of range'
      if(k.lt.anint(w(iw+4 )) .or. k.gt.anint(w(iw+5 ))) stop
     +  'iqcWCadr: index 3 (k) out of range'
      if(l.lt.anint(w(iw+6 )) .or. l.gt.anint(w(iw+7 ))) stop
     +  'iqcWCadr: index 4 (l) out of range'
      if(m.lt.anint(w(iw+8 )) .or. m.gt.anint(w(iw+9 ))) stop
     +  'iqcWCadr: index 5 (m) out of range'
      iqcWCadr = int(anint(w(iw+10)) + 
     +           anint(w(iw+11))*i + anint(w(iw+12))*j +
     +           anint(w(iw+13))*k + anint(w(iw+14))*l +
     +           anint(w(iw+15))*m)
     
      return
      end 

C     ==================================
      integer function iqcW1ijk(w,i,j,k)
C     ==================================

C--   Fast type 1 indexing: acts as a do-nothing if store
C--   not partitioned or if no tables of type 1 exist

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension w(*)

C--                         ix iq nf ig id
C--   iqcW1ijk = iqcWaddr(w, i, x, x, j, k)
C--                         11 12 13 14 15

      iqcW1ijk = 0
      if(anint(w(1)).ne.123456)    return
      ityp = k/100
      if(ityp.le.0 .or. ityp.ge.5) return
      iw = int(anint(w(ityp+3+miw0)))
      if(iw.eq.0)                  return
      iqcW1ijk = int(anint(w(iw+10)) + 
     +           anint(w(iw+11))*i + 
     +           anint(w(iw+14))*j +
     +           anint(w(iw+15))*k)
     
      return
      end

C     =====================================
      integer function iqcW2ijkl(w,i,j,k,l)
C     =====================================

C--   Convenient type 2 indexing: acts as a do-nothing if store
C--   not partitioned or if no tables of type 2 exist

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension w(*)

C--                          ix iq nf ig id
C--   iqcW2ijkl = iqcWaddr(w, i, x, j, k, l)
C--                          11 12 13 14 15
      
      iqcW2ijkl = 0
      if(anint(w(1)).ne.123456)    return
      ityp = l/100
      if(ityp.le.0 .or. ityp.ge.5) return
      iw = int(anint(w(ityp+3+miw0)))
      if(iw.eq.0)                  return
      iqcW2ijkl = int(anint(w(iw+10)) + 
     +            anint(w(iw+11))*i + 
     +            anint(w(iw+13))*j +
     +            anint(w(iw+14))*k +
     +            anint(w(iw+15))*l)
     
      return
      end 
      
C     =====================================
      integer function iqcW3ijkl(w,i,j,k,l)
C     =====================================

C--   Convenient type 3 indexing: acts as a do-nothing if store
C--   not partitioned or if no tables of type 3 exist

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension w(*)

C--                          ix iq nf ig id
C--   iqcW3ijkl = iqcWaddr(w, i, j, x, k, l)
C--                          11 12 13 14 15 
      
      iqcW3ijkl = 0
      if(anint(w(1)).ne.123456)    return
      ityp = l/100
      if(ityp.le.0 .or. ityp.ge.5) return
      iw = int(anint(w(ityp+3+miw0)))
      if(iw.eq.0)                  return
      iqcW3ijkl = int(anint(w(iw+10)) + 
     +            anint(w(iw+11))*i + 
     +            anint(w(iw+12))*j +
     +            anint(w(iw+14))*k +
     +            anint(w(iw+15))*l)
     
      return
      end       

C     ========================================
      integer function iqcW4ijklm(w,i,j,k,l,m)
C     ========================================

C--   Convenient type 4 indexing: acts as a do-nothing if store
C--   not partitioned or if no tables of type 4 exist

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension w(*)

C--                           ix iq nf ig id
      iqcW4ijklm = iqcWaddr(w, i, j, k, l, m)
     
      return
      end 

C==   ==============================================================
C==   Utilities and checks =========================================
C==   ==============================================================        

C     ============================================
      subroutine sqcTabLims(w,ityp,imin,imax,ierr)
C     ============================================

C--   Get table limits from the store
C--
C--   Input  :  w        store dimensioned in the calling routine
C--             ityp     table type
C--
C--   Output :  imin(i)  lower index limits
C--             imax(i)  upper index limits
C--             ierr     1 = no partition or tables ityp do not exist

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension w(*),imin(*),imax(*)

      ierr = 1
      if(anint(w(1)).ne.123456)    return
      if(ityp.le.0 .or. ityp.ge.5) return
      ia = int(anint(w(ityp+3+miw0)))
      if(ia.eq.0)                  return
      ierr = 0
      ia   = ia-1
      do i = 1,5
        ia = ia+1
        imin(i) = int(anint(w(ia)))
        ia = ia+1
        imax(i) = int(anint(w(ia)))
      enddo

      return
      end

C     =============================
      subroutine sqcSetNoWt(w,ityp)
C     =============================

C--   Write 123456 in the first word of all tables of
C--   type ityp. This flags the table as not filled

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'

      dimension w(*),mi(5),ma(5)

      call sqcTabLims(w,ityp,mi,ma,ierr)

C--   No tables
      if(ierr.eq.1) return
C--   Loop over tables
      do i = mi(5),ma(5)
          ia    = iqcWaddr(w,mi(1),mi(2),mi(3),mi(4),i)
          w(ia) = 123456
      enddo

      return
      end
    
C     =================================      
      integer function iqcCheckId(w,id)
C     =================================

C--   Checks if table id exists in store w
C--
C--   -1   : Table exists but is empty
C--    0   : OK
C--    1   : No tables booked
C--    2   : Table type does not exist
C--    3   : No tables booked of requested type
C--    4   : Id does not exist

      implicit double precision (a-h,o-z)
      
      dimension w(*)
      dimension mi(5),ma(5)
      
C--   Check store partitioned
      if(anint(w(1)).ne.123456) then
        iqcCheckId = 1
        return
      endif
C--   Check global range
      if(id.lt.101 .or. id.gt.499) then
        iqcCheckId = 2
        return
      endif      
C--   Target table type
      ityp = id/100
C--   Get index limits      
      call sqcTabLims(w,ityp,mi,ma,ierr)
C--   No type-x tables in the store
      if(ierr.ne.0) then
        iqcCheckId = 3
        return
      endif
C--   Check id is in range
      if(id.lt.mi(5) .or. id.gt.ma(5)) then
        iqcCheckId = 4
        return
      endif
C--   Check table is empty
      ia = iqcWaddr(w,mi(1),mi(2),mi(3),mi(4),id)
      if(anint(w(ia)) .eq. 123456) then
        iqcCheckId = -1 
        return
      endif
C--   OK ...
      iqcCheckId = 0      
      
      return
      end
      
C     ===================================================      
      subroutine sqcChekPij(idin,idout,jset,ioffset,ierr)
C     ===================================================

C--   To recognise an internal splitting function table
C--   of jset = 1(unpol), 2(pol), 3(timelike), 4(custom),
C--   the identifier is coded as -(1000*jset+id). 
C--   This routine decodes an input index idin into jset and
C--   idout, and checks if jset is OK. Then it calls iqcCheckid
C--   to check if idout is OK
C--
C--   idin    (in)  : coded identifier -(1000*jset+idout)
C--   idout   (out) : identifier of splitting function table 
C--   jset    (out) : evolution type unpol/pol/timelike/custom
C--   ioffset (out) : first word of jset in stor7
C--   ierr    (out) : -1 to 4 as given by iqcCheckId

      implicit double precision (a-h,o-z)
      
      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qstor7.inc'
                                 
      jset =abs(idin)/1000
      if(jset.lt.1 .or. jset.gt.4) then
        ierr = 4
        return
      endif
      if(mxord7(jset).eq.0) then
        ierr = 4
        return
      endif
      ioffset = ifst7(jset,ioy2)
      idout   = abs(idin)-1000*jset
      ierr    = iqcCheckId(stor7(ioffset),idout)

      return
      end

C==   ==============================================================
C==   Weight calculation ===========================================
C==   ==============================================================

C     =============================================
      subroutine sqcUweitA(w,id,ioy,afun,achi,ierr)
C     =============================================

C--   Calculate weights for regular piece

C--   ierr = 0  : OK
C--          1  : achi(qmu2) < 1

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

      dimension w(*)
      dimension mi(5),ma(5)

      external afun, achi

      logical lqcAcomp

      ityp = id/100
     
C--   Get index limits      
      call sqcTabLims(w,ityp,mi,ma,ierr)

C--   No table
      if(ierr.ne.0) stop 'sqcUweitA: no weight table'

C--   Set first word of table to zero
      ia = iqcWaddr(w,mi(1),mi(2),mi(3),mi(4),id)
      if(anint(w(ia)) .eq. 123456) then
        w(ia) = 0.D0
*mb        write(6,*) 'sqcUweitA: reset table id = ',id
      endif

C--   Prepare fast indexing
      inc1 = iqcWaddr(w,1,0,0,0,id)-iqcWaddr(w,0,0,0,0,id) 
      inc2 = iqcWaddr(w,0,1,0,0,id)-iqcWaddr(w,0,0,0,0,id)
      inc3 = iqcWaddr(w,0,0,1,0,id)-iqcWaddr(w,0,0,0,0,id)
      inc4 = iqcWaddr(w,0,0,0,1,id)-iqcWaddr(w,0,0,0,0,id)
      ia4  = iqcWaddr(w,mi(1),mi(2),mi(3),mi(4),id)-inc4
C--   Nested loop: depending on the table type, several of these are 
C--   one-trip loops. The index ranges are set by sqcTabLims above
C--   Loop over subgrids
      do ig = mi(4),ma(4)
        ia4 = ia4+inc4
        ia3 = ia4-inc3
        del = dely2(ig)
C--     Loop over number of flavors
        do nf = mi(3),ma(3)
          ia3 = ia3+inc3
          ia2 = ia3-inc2
C--       Loop over t-grid
          do it = mi(2),ma(2)
            ia2 = ia2+inc2
            ia1 = ia2-inc1
            ti  = tgrid2(it)
            aa  = achi(exp(ti))
C--         Check
            if(lqcAcomp(aa,1.D0,aepsi6)) then
              aa = 1.D0
            elseif(aa.lt.1.D0)           then
              ierr = 1
              return
            endif
            bb  = log(aa)
C--         Loop over y-subgrid
            do iy = 1,nyy2(ig)
              yi   = iy*del
              yb   = yi-bb
              weit = 0.D0
              if(yb.gt.0.D0) then
                a    = 0.D0             !lower integration limit
                b    = min(yb,ioy*del)  !upper integration limit
                weit = dqcUAgauss(ioy-1,afun,yb,ti,nf,a,b,del)
                weit = weit/aa
              endif
C--           Add weight
              ia1    = ia1+inc1
              w(ia1) = w(ia1)+weit
C--         End of loop over y-subgrid
            enddo
C--       End of loop over t-grid
          enddo
C--     End of loop over flavors
        enddo
C--   End of loop over subgrids
      enddo

      return
      end

C     ==================================================
      subroutine sqcUweitB(w,id,ioy,bfun,achi,idel,ierr)
C     ==================================================

C--   Calculate weights for singular piece

C--   idel = 0  : dont include delta(1-x) piece
C--             : otherwise include this piece
C--   ierr = 0  : OK
C--          1  : achi(qmu2) < 1

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

      dimension w(*)
      dimension mi(5),ma(5)

      external bfun, achi

      logical lqcAcomp

      ityp = id/100

C--   Get index limits      
      call sqcTabLims(w,ityp,mi,ma,ierr)

C--   No table
      if(ierr.ne.0) stop 'sqcUweitb: no weight table'

C--   Set first word of table to zero
      ia = iqcWaddr(w,mi(1),mi(2),mi(3),mi(4),id)
      if(anint(w(ia)) .eq. 123456) then
        w(ia) = 0.D0
*mb        write(6,*) 'sqcUweitB: reset table id = ',id
      endif

C--   Prepare fast indexing
      inc1 = iqcWaddr(w,1,0,0,0,id)-iqcWaddr(w,0,0,0,0,id) 
      inc2 = iqcWaddr(w,0,1,0,0,id)-iqcWaddr(w,0,0,0,0,id)
      inc3 = iqcWaddr(w,0,0,1,0,id)-iqcWaddr(w,0,0,0,0,id)
      inc4 = iqcWaddr(w,0,0,0,1,id)-iqcWaddr(w,0,0,0,0,id)
      ia4  = iqcWaddr(w,mi(1),mi(2),mi(3),mi(4),id)-inc4
C--   Nested loop: depending on the table type, several of these are 
C--   one-trip loops. The index ranges are set by sqcTabLims above
C--   Loop over subgrids
      do ig = mi(4),ma(4)
        ia4 = ia4+inc4
        ia3 = ia4-inc3
        del = dely2(ig)
C--     Loop over number of flavors
        do nf = mi(3),ma(3)
          ia3 = ia3+inc3
          ia2 = ia3-inc2
C--       Loop over t-grid
          do it = mi(2),ma(2)
            ia2 = ia2+inc2
            ia1 = ia2-inc1
            ti = tgrid2(it)
            aa  = achi(exp(ti))
C--         Check
            if(lqcAcomp(aa,1.D0,aepsi6)) then
              aa = 1.D0
            elseif(aa.lt.1.D0)           then
              ierr = 1
              return
            endif       
            bb  = log(aa)
C--         Loop over y-subgrid
            do iy = 1,nyy2(ig)
              yi     = iy*del
              yb     = yi-bb
              weit   = 0.D0
              if(yb.gt.0.D0) then
                xi     = exp(-yi)
                ax     = aa*xi
                a      = 0.D0             !lower integration limit
                b      = min(yb,ioy*del)  !upper integration limit
                wgt1   = dqcUBgauss(ioy-1,bfun,yb,ti,nf,a,b,del)
                if(idel.ne.0) then
                  wgt2   = dqcBsplyy(ioy-1,1,yb/del) * 
     +                     dqcUIgauss(bfun,ti,nf,0.D0,ax)
                else
                  wgt2 = 0.D0
                endif
                weit   = (wgt1-wgt2)/aa
              endif
C--           Add weight
              ia1    = ia1+inc1
              w(ia1) = w(ia1)+weit
C--         End of loop over y-subgrid
            enddo
C--       End of loop over t-grid
          enddo
C--     End of loop over flavors
        enddo
C--   End of loop over subgrids
      enddo

      return
      end

C     =======================================================
      subroutine sqcUwgtRS(w,id,ioy,rfun,sfun,achi,idel,ierr)
C     =======================================================

C--   Calculate weights for RS piece

C--   idel = 0  : dont include delta(1-x) piece
C--             : otherwise include this piece
C--   ierr = 0  : OK
C--          1  : achi(qmu2) < 1

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

      dimension w(*)
      dimension mi(5),ma(5)

      external rfun, sfun, achi

      logical lqcAcomp

      ityp = id/100

C--   Get index limits      
      call sqcTabLims(w,ityp,mi,ma,ierr)

C--   No table
      if(ierr.ne.0) stop 'sqcUwgtRS: no weight table'

C--   Set first word of table to zero
      ia = iqcWaddr(w,mi(1),mi(2),mi(3),mi(4),id)
      if(anint(w(ia)) .eq. 123456) then
        w(ia) = 0.D0
*mb        write(6,*) 'sqcUwgtRS: reset table id = ',id
      endif

C--   Prepare fast indexing
      inc1 = iqcWaddr(w,1,0,0,0,id)-iqcWaddr(w,0,0,0,0,id) 
      inc2 = iqcWaddr(w,0,1,0,0,id)-iqcWaddr(w,0,0,0,0,id)
      inc3 = iqcWaddr(w,0,0,1,0,id)-iqcWaddr(w,0,0,0,0,id)
      inc4 = iqcWaddr(w,0,0,0,1,id)-iqcWaddr(w,0,0,0,0,id)
      ia4  = iqcWaddr(w,mi(1),mi(2),mi(3),mi(4),id)-inc4
C--   Nested loop: depending on the table type, several of these are 
C--   one-trip loops. The index ranges are set by sqcTabLims above
C--   Loop over subgrids
      do ig = mi(4),ma(4)
        ia4  = ia4+inc4
        ia3  = ia4-inc3
        del  = dely2(ig)
C--     Loop over number of flavors
        do nf = mi(3),ma(3)
          ia3 = ia3+inc3
          ia2 = ia3-inc2
C--       Loop over t-grid
          do it = mi(2),ma(2)
            ia2 = ia2+inc2
            ia1 = ia2-inc1
            ti = tgrid2(it)
            aa  = achi(exp(ti))
C--         Check
            if(lqcAcomp(aa,1.D0,aepsi6)) then
              aa = 1.D0
            elseif(aa.lt.1.D0)           then
              ierr = 1
              return
            endif      
            bb  = log(aa)
C--         Loop over y-subgrid
            do iy = 1,nyy2(ig)
              yi   = iy*del
              yb   = yi-bb
              weit = 0.D0
              if(yb.gt.0.D0) then
                xi   = exp(-yi)
                ax   = aa*xi
                a    = 0.D0             !lower integration limit
                b    = min(yb,ioy*del)  !upper integration limit
                wgt1 = dqcURSgaus(
     +                 ioy-1,rfun,sfun,yb,ti,nf,a,b,del)
                if(idel.ne.0) then
                  wgt2 = rfun(1.D0,ti,nf)          *
     +                   dqcBsplyy(ioy-1,1,yb/del) *
     +                   dqcUIgauss(sfun,ti,nf,0.D0,ax)
                else
                  wgt2 = 0.D0
                endif
                weit = (wgt1-wgt2)/aa
              endif
C--           Add weight
              ia1    = ia1+inc1
              w(ia1) = w(ia1)+weit
C--         End of loop over y-subgrid
            enddo
C--       End of loop over t-grid
          enddo
C--     End of loop over flavors
        enddo
C--   End of loop over subgrids
      enddo

      return
      end

C     =============================================
      subroutine sqcUweitD(w,id,ioy,dfun,achi,ierr)
C     =============================================

C--   Calculate weights for D(x)*delta(1-x)
C--
C--   ierr = 0  : OK
C--          1  : achi(qmu2) < 1

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

      dimension w(*)
      dimension mi(5),ma(5)

      external dfun, achi

      logical lqcAcomp

      ityp = id/100

C--   Get index limits      
      call sqcTabLims(w,ityp,mi,ma,ierr)

C--   No table
      if(ierr.ne.0) stop 'sqcUweitD: no weight table'

C--   Set first word of table to zero
      ia = iqcWaddr(w,mi(1),mi(2),mi(3),mi(4),id)
      if(anint(w(ia)) .eq. 123456) then
        w(ia) = 0.D0
*mb        write(6,*) 'sqcUweitD: reset table id = ',id
      endif

C--   Prepare fast indexing
      inc1 = iqcWaddr(w,1,0,0,0,id)-iqcWaddr(w,0,0,0,0,id) 
      inc2 = iqcWaddr(w,0,1,0,0,id)-iqcWaddr(w,0,0,0,0,id)
      inc3 = iqcWaddr(w,0,0,1,0,id)-iqcWaddr(w,0,0,0,0,id)
      inc4 = iqcWaddr(w,0,0,0,1,id)-iqcWaddr(w,0,0,0,0,id)
      ia4  = iqcWaddr(w,mi(1),mi(2),mi(3),mi(4),id)-inc4
C--   Nested loop: depending on the table type, several of these are 
C--   one-trip loops. The index ranges are set by sqcTabLims above
C--   Loop over subgrids
      do ig = mi(4),ma(4)
        ia4  = ia4+inc4
        ia3  = ia4-inc3
        del  = dely2(ig)
C--     Loop over number of flavors
        do nf = mi(3),ma(3)
          ia3 = ia3+inc3
          ia2 = ia3-inc2
C--       Loop over t-grid
          do it = mi(2),ma(2)
            ia2 = ia2+inc2
            ia1 = ia2-inc1
            ti  = tgrid2(it)
            aa  = achi(exp(ti))
C--         Check
            if(lqcAcomp(aa,1.D0,aepsi6)) then
              aa = 1.D0
            elseif(aa.lt.1.D0)           then
              ierr = 1
              return
            endif
            bb  = log(aa)
C--         Loop over y-subgrid
            do iy = 1,nyy2(ig)
              yi   = iy*del
              yb   = yi-bb
              wgt  = 0.D0
              if(yb.gt.0.D0) then
                wgt  = dfun(exp(-yb),exp(ti),nf) * 
     +                 dqcBsplyy(ioy-1,1,yb/del)
                wgt  = wgt/aa
              endif
C--           Add weight
              ia1    = ia1+inc1
              w(ia1) = w(ia1)+wgt
C--         End of loop over y-subgrid
            enddo
C--       End of loop over t-grid
          enddo
C--     End of loop over flavors
        enddo
C--   End of loop over subgrids
      enddo

      return
      end
      
C     ===================================
      subroutine sqcUweitX(w,id,ioy,ierr)
C     ===================================

C--   Calculate weights for convolution F cross F

C--   ierr = 0  : OK
C--          1  : error (should never occur)

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

      dimension w(*)
      dimension mi(5),ma(5)

      ityp = id/100
     
C--   Get index limits      
      call sqcTabLims(w,ityp,mi,ma,ierr)

C--   No table
      if(ierr.ne.0) stop 'sqcUweitX: no weight table'

C--   Set first word of table to zero
      ia = iqcWaddr(w,mi(1),mi(2),mi(3),mi(4),id)
      if(anint(w(ia)) .eq. 123456) then
        w(ia) = 0.D0
*mb        write(6,*) 'sqcUweitX: reset table id = ',id
      endif

C--   Prepare fast indexing
      inc1 = iqcWaddr(w,1,0,0,0,id)-iqcWaddr(w,0,0,0,0,id) 
      inc2 = iqcWaddr(w,0,1,0,0,id)-iqcWaddr(w,0,0,0,0,id)
      inc3 = iqcWaddr(w,0,0,1,0,id)-iqcWaddr(w,0,0,0,0,id)
      inc4 = iqcWaddr(w,0,0,0,1,id)-iqcWaddr(w,0,0,0,0,id)
      ia4  = iqcWaddr(w,mi(1),mi(2),mi(3),mi(4),id)-inc4
C--   Nested loop: depending on the table type, several of these are 
C--   one-trip loops. The index ranges are set by sqcTabLims above
C--   Loop over subgrids
      do ig = mi(4),ma(4)
        ia4 = ia4+inc4
        ia3 = ia4-inc3
        del = dely2(ig)
C--     Loop over number of flavors
        do nf = mi(3),ma(3)
          ia3 = ia3+inc3
          ia2 = ia3-inc2
C--       Loop over t-grid
          do it = mi(2),ma(2)
            ia2 = ia2+inc2
            ia1 = ia2-inc1
            ti  = tgrid2(it)
C--         Loop over y-subgrid
            do iy = 1,nyy2(ig)
              yi   = iy*del
              weit = 0.D0
              a    = 0.D0             !lower integration limit
              b    = min(yi,ioy*del)  !upper integration limit
              weit = dqcUXgauss(ioy-1,yi,a,b,del)
C--           Store weight
              ia1    = ia1+inc1
              w(ia1) = weit
C--         End of loop over y-subgrid
            enddo
C--       End of loop over t-grid
          enddo
C--     End of loop over flavors
        enddo
C--   End of loop over subgrids
      enddo

      return
      end
      
C==   ==============================================================
C==   Operations on weight tables ==================================
C==   ==============================================================      

C     =============================
      subroutine sqcScaleWt(w,c,id)
C     =============================

C--   Multiply weight table by constant factor c

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

      dimension w(*)
      dimension mi(5),ma(5)

      ityp = id/100

C--   Get index limits      
      call sqcTabLims(w,ityp,mi,ma,ierr)

C--   No table
      if(ierr.ne.0) stop 'sqcScaleWt: no weight table'

C--   First and last word of table
      ia1 = iqcWaddr(w,mi(1),mi(2),mi(3),mi(4),id)
      ia2 = iqcWaddr(w,ma(1),ma(2),ma(3),ma(4),id)

      do ia = ia1,ia2
        w(ia) = c*w(ia)
      enddo

      return
      end
      
C     ========================================
      subroutine sqcCopyWt(w1,id1,w2,id2,iadd)
C     ========================================

C--   Copy contents of table id1 in w1 to table id2 in w2
C--   iadd = -1,0,1    subtract/copy/add  id1 to id2
C--   Type of id2 should be .ge. type id1 but this is not checked

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'

      dimension w1(*),w2(*)
      dimension mi1(5),ma1(5),mi2(5),ma2(5)

C--   Get index limits      
      ityp1 = id1/100
      call sqcTabLims(w1,ityp1,mi1,ma1,ierr)
      ityp2 = id2/100
      call sqcTabLims(w2,ityp2,mi2,ma2,ierr)

C--   First word of output table
      ia2 = iqcWaddr(w2,mi2(1),mi2(2),mi2(3),mi2(4),id2)
C--   Set first word of table to zero
      if(anint(w2(ia2)) .eq. 123456) then
        w2(ia2) = 0.D0
*mb        write(6,*) 'sqcCopyWt: reset table id = ',id2
      endif
      
C--   Loop over output table
      do ig2 = mi2(4),ma2(4)              !loop over grids
C--     Bracket ig1 limits
        ig1 = max(mi1(4),ig2)
        ig1 = min(ma1(4),ig1)
        do nf2 = mi2(3),ma2(3)            !loop over flavors
C--       Bracket nf1  limits 
          nf1 = max(mi1(3),nf2)
          nf1 = min(ma1(3),nf1)
          do iq2 = mi2(2),ma2(2)          !loop over mu2
C--         Bracket iq1 limits
            iq1 = max(mi1(2),iq2)
            iq1 = min(ma1(2),iq1)
C--         Base address of first x-bin
            ia1 = iqcWaddr(w1,mi1(1),iq1,nf1,ig1,id1)-1 
            ia2 = iqcWaddr(w2,mi2(1),iq2,nf2,ig2,id2)-1
            if(iadd.eq.-1) then
C--           id2 = id2 - id1
              do i = mi2(1),ma2(1)
                ia1     = ia1+1
                ia2     = ia2+1
                w2(ia2) = w2(ia2) - w1(ia1)
              enddo              
            elseif(iadd.eq.0) then
C--           id2 = id1
              do i = mi2(1),ma2(1)
                ia1    = ia1+1
                ia2    = ia2+1
                w2(ia2) = w1(ia1)
              enddo               
            elseif(iadd.eq.+1) then
C--           id2 = id2 + id1
              do i = mi2(1),ma2(1)
                ia1 = ia1+1
                ia2 = ia2+1
                w2(ia2) = w2(ia2) + w1(ia1)
              enddo   
            else
              stop 'sqcCopyWt: invalid iadd'
            endif  
          enddo
        enddo
      enddo

      return
      end
      
C     ==========================================================
      subroutine sqcWcrossW(wa,ida,wb,idb,wc,idc,idt1,idt2,iadd)
C     ==========================================================

C--   Store in idc the weight of ida * idb :   Wc = Wa S^{-1} Wb
C--   wa, wb and wc are the stores of ida, idb and idc
C--   idt1 and idt2 are scratch buffers

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      dimension wa(*), wb(*), wc(*)
      dimension mia(5),maa(5),mib(5),mab(5),mic(5),mac(5)
      
C--   Find pdf set for temp buffers
      jset = 0
      do j = 1,9
        if(mxord7(j).ne.0) jset = j 
      enddo 
      if(jset.eq.0) stop 'sqcWcrossW: no pdf set available'         
      
C--   Get index limits
      itypa = ida/100
      call sqcTabLims(wa,itypa,mia,maa,ierr) 
      itypb = idb/100
      call sqcTabLims(wb,itypb,mib,mab,ierr) 
      itypc = idc/100
      call sqcTabLims(wc,itypc,mic,mac,ierr)
C--   Set first word of table to zero
      ia = iqcWaddr(wc,mic(1),mic(2),mic(3),mic(4),idc)
      if(anint(wc(ia)) .eq. 123456) then
        wc(ia) = 0.D0
      endif
      
C--   Loop over output table
      do igc = mic(4),mac(4)              !loop over grids
C--     Bracket iga and igb limits
        iga = max(mia(4),igc)
        iga = min(maa(4),iga)
        igb = max(mib(4),igc)
        igb = min(mab(4),igb)
        do nfc = mic(3),mac(3)            !loop over flavors
C--       Bracket nfa and nfb limits 
          nfa = max(mia(3),nfc)
          nfa = min(maa(3),nfa)
          nfb = max(mib(3),nfc)
          nfb = min(mab(3),nfb)       
          do iqc = mic(2),mac(2)          !loop over mu2
C--         Bracket iqa and iqb limits
            iqa = max(mia(2),iqc)
            iqa = min(maa(2),iqa)
            iqb = max(mib(2),iqc)
            iqb = min(mab(2),iqb)
C--         Address of first x-bin
            ia1 = iqcWaddr(wa,mia(1),iqa,nfa,iga,ida) 
            ib1 = iqcWaddr(wb,mib(1),iqb,nfb,igb,idb) 
            ic1 = iqcWaddr(wc,mic(1),iqc,nfc,igc,idc)
            it1 = iqcPdfIjkl(1,1,idt1,jset)
            it2 = iqcPdfIjkl(1,1,idt2,jset)            
C--         First do Z = S^{-1}B   (store in it1)
            call sqcABmult(sinvy2(1),wb(ib1),stor7(it1),nyy2(0))
C--         Now do C = AZ (store in it2)
            call sqcABmult(wa(ia1),stor7(it1),stor7(it2),nyy2(0))
            ic1 = ic1-1
            it2 = it2-1
            if(iadd.eq.-1) then
C--           C = C - convol
              do i = mic(1),mac(1)
                ic1     = ic1+1
                it2     = it2+1
                wc(ic1) = wc(ic1) - stor7(it2)
              enddo              
            elseif(iadd.eq.0) then
C--           C = convol
              do i = mic(1),mac(1)
                ic1     = ic1+1
                it2     = it2+1
                wc(ic1) = stor7(it2)
              enddo               
            elseif(iadd.eq.+1) then
C--           C = C + convol
              do i = mic(1),mac(1)
                ic1     = ic1+1
                it2     = it2+1
                wc(ic1) = wc(ic1) + stor7(it2)
              enddo   
            else
              stop 'sqcWcrossW: invalid iadd'
            endif  
          enddo
        enddo
      enddo

      return
      end

C     =============================================
      subroutine sqcWtimesF(fun,w1,id1,w2,id2,iadd) 
C     =============================================

C--   Multiply id1 by fun(iq,nf) and store result in id2
C--   iadd = -1,0,1   subtract,store,add  result to id2

      implicit double precision (a-h,o-z)

      include 'qcdnum.inc'
      include 'qgrid2.inc'
      include 'qpars6.inc'
      include 'qstor7.inc'

      dimension w1(*),w2(*)
      dimension mi1(5),ma1(5),mi2(5),ma2(5)
      
      external fun
      
C--   Get index limits
      ityp1 = id1/100
      call sqcTabLims(w1,ityp1,mi1,ma1,ierr)
      ityp2 = id2/100
      call sqcTabLims(w2,ityp2,mi2,ma2,ierr) 
C--   Set first word of output table to zero
      ia = iqcWaddr(w2,mi2(1),mi2(2),mi2(3),mi2(4),id2)
      if(anint(w2(ia)) .eq. 123456) then
        w2(ia) = 0.D0
      endif
      
C--   Loop over output table
      do ig2 = mi2(4),ma2(4)              !loop over grids
C--     Bracket ig1 limits
        ig1 = max(mi1(4),ig2)
        ig1 = min(ma1(4),ig1)
        do nf2 = mi2(3),ma2(3)            !loop over flavors
C--       Bracket nf1  limits 
          nf1 = max(mi1(3),nf2)
          nf1 = min(ma1(3),nf1)
          do iq2 = mi2(2),ma2(2)          !loop over mu2
C--         Bracket iq1 limits
            iq1 = max(mi1(2),iq2)
            iq1 = min(ma1(2),iq1)
C--         Base address of first x-bin  (with boundary check for the moment)
            ia1 = iqcWaddr(w1,mi1(1),iq1,nf1,ig1,id1)-1 
            ia2 = iqcWaddr(w2,mi2(1),iq2,nf2,ig2,id2)-1
            fac = fun(iq2,nf2) 
            if(iadd.eq.-1) then
C--           id2 = id2 - f*id1
              do i = mi2(1),ma2(1)
                ia1     = ia1+1
                ia2     = ia2+1
                w2(ia2) = w2(ia2) - fac*w1(ia1)
              enddo              
            elseif(iadd.eq.0) then
C--           id2 = f*id1
              do i = mi2(1),ma2(1)
                ia1     = ia1+1
                ia2     = ia2+1
                w2(ia2) = fac*w1(ia1)
              enddo               
            elseif(iadd.eq.+1) then
C--           id2 = id2 + f*id1
              do i = mi2(1),ma2(1)
                ia1     = ia1+1
                ia2     = ia2+1
                w2(ia2) = w2(ia2) + fac*w1(ia1)
              enddo   
            else
              stop 'sqcWtimesF: invalid iadd'
            endif  
          enddo
        enddo
      enddo

      return
      end

C==   ==============================================================
C==   Gauss integration ============================================
C==   ==============================================================

C     ====================================================
      double precision function dqcUIgauss(pfun,ti,nf,a,b)
C     ====================================================
 
C--   Integrate coeffient function pfun(x,q,nf) from a to b.
C--   NB: integral over x and not over y = -ln x!
C--
C--   Modified Cernlib routine DGAUSS D103.
C--   Taken from /afs/cern.ch/asis/share/cern/2001/src/...
C--           .../mathlib/gen/d/gauss64.F and gausscod.inc

      implicit double precision (a-h,o-z)

C--   MB extension
      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
      external pfun
C--   End of MB extension

      DIMENSION W(12),X(12)

      PARAMETER (Z1 = 1, HF = Z1/2, CST = 5*Z1/1000)

      DATA X( 1) /9.6028985649753623D-1/, W( 1) /1.0122853629037626D-1/
      DATA X( 2) /7.9666647741362674D-1/, W( 2) /2.2238103445337447D-1/
      DATA X( 3) /5.2553240991632899D-1/, W( 3) /3.1370664587788729D-1/
      DATA X( 4) /1.8343464249564980D-1/, W( 4) /3.6268378337836198D-1/
      DATA X( 5) /9.8940093499164993D-1/, W( 5) /2.7152459411754095D-2/
      DATA X( 6) /9.4457502307323258D-1/, W( 6) /6.2253523938647893D-2/
      DATA X( 7) /8.6563120238783174D-1/, W( 7) /9.5158511682492785D-2/
      DATA X( 8) /7.5540440835500303D-1/, W( 8) /1.2462897125553387D-1/
      DATA X( 9) /6.1787624440264375D-1/, W( 9) /1.4959598881657673D-1/
      DATA X(10) /4.5801677765722739D-1/, W(10) /1.6915651939500254D-1/
      DATA X(11) /2.8160355077925891D-1/, W(11) /1.8260341504492359D-1/
      DATA X(12) /9.5012509837637440D-2/, W(12) /1.8945061045506850D-1/

C--   MB extension
      f(u) = pfun(u,exp(ti),nf)
      
      eps  = gepsi6
C--   End of MB extension
 
      H=0
      IF(B .EQ. A) GO TO 99
      CONST=CST/ABS(B-A)
      BB=A
    1 AA=BB
      BB=B
    2 C1=HF*(BB+AA)
      C2=HF*(BB-AA)
      S8=0
      DO 3 I = 1,4
      U=C2*X(I)
    3 S8=S8+W(I)*(F(C1+U)+F(C1-U))
      S16=0
      DO 4 I = 5,12
      U=C2*X(I)
    4 S16=S16+W(I)*(F(C1+U)+F(C1-U))
      S16=C2*S16
      IF(ABS(S16-C2*S8) .LE. EPS*(1+ABS(S16))) THEN
       H=H+S16
       IF(BB .NE. B) GO TO 1
      ELSE
       BB=C1
       IF(1+CONST*ABS(C2) .NE. 1) GO TO 2
*mb    H=0
C--    MB extension
       WRITE(lunerr1,'(/'' dqcUIgauss: too high accuracy required'','//
     +         '  '' ---> STOP'')')
       STOP
*mb    ierr = 1
C--    End of MB extension
       goto 99
      END IF

   99 dqcUIgauss=H

      RETURN
      END

C     ========================================
      double precision function 
     +   dqcUAgauss(idk,afun,yi,ti,nf,a,b,del)
C     ========================================

C--   Integrate regular piece (A) of coeffient function.
C--
C--   I(yi) = Int_a^b dz Abar(yi-z) B1(z) with Abar(z) = exp(-z)A(exp(-z)).
C--
C--   idk        = spline order iord-1: 1=linear and 2=quadratic
C--   afun       = coefficient function versus x = exp(-y)
C--   yi         = upper limit of convolution (exp(-yi) is argument of afun)
C--   ti         = t-value (exp(t) is argument of afun)      
C--   nf         = number of flavors (argument of afun) 
C--   a          = lower limit of integral = 0.D0
C--   b          = upper limit of integral = min[yi,iord*del]
C--   del        = grid spacing
C--
C--   Modified Cernlib routine DGAUSS D103.
C--   Taken from /afs/cern.ch/asis/share/cern/2001/src/...
C--           .../mathlib/gen/d/gauss64.F and gausscod.inc

      implicit double precision (a-h,o-z)

C--   MB extension
      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
      external afun
C--   End of MB extension

      DIMENSION W(12),X(12)

      PARAMETER (Z1 = 1, HF = Z1/2, CST = 5*Z1/1000)

      DATA X( 1) /9.6028985649753623D-1/, W( 1) /1.0122853629037626D-1/
      DATA X( 2) /7.9666647741362674D-1/, W( 2) /2.2238103445337447D-1/
      DATA X( 3) /5.2553240991632899D-1/, W( 3) /3.1370664587788729D-1/
      DATA X( 4) /1.8343464249564980D-1/, W( 4) /3.6268378337836198D-1/
      DATA X( 5) /9.8940093499164993D-1/, W( 5) /2.7152459411754095D-2/
      DATA X( 6) /9.4457502307323258D-1/, W( 6) /6.2253523938647893D-2/
      DATA X( 7) /8.6563120238783174D-1/, W( 7) /9.5158511682492785D-2/
      DATA X( 8) /7.5540440835500303D-1/, W( 8) /1.2462897125553387D-1/
      DATA X( 9) /6.1787624440264375D-1/, W( 9) /1.4959598881657673D-1/
      DATA X(10) /4.5801677765722739D-1/, W(10) /1.6915651939500254D-1/
      DATA X(11) /2.8160355077925891D-1/, W(11) /1.8260341504492359D-1/
      DATA X(12) /9.5012509837637440D-2/, W(12) /1.8945061045506850D-1/

C--   MB extension
C--   Define function to integrate
      f(u) = dqcBsplyy(idk,1,u/del) * exp(-(yi-u)) * 
     +       afun(exp(-(yi-u)),exp(ti),nf)
C--   Check limits
      if(b.le.a) then
        dqcUAGauss = 0.D0
        return
      endif
      eps  = gepsi6
C--   End of MB extension

      H=0
      IF(B .EQ. A) GO TO 99
      CONST=CST/ABS(B-A)
      BB=A
    1 AA=BB
      BB=B
    2 C1=HF*(BB+AA)
      C2=HF*(BB-AA)
      S8=0
      DO 3 I = 1,4
      U=C2*X(I)
    3 S8=S8+W(I)*(F(C1+U)+F(C1-U))
      S16=0
      DO 4 I = 5,12
      U=C2*X(I)
    4 S16=S16+W(I)*(F(C1+U)+F(C1-U))
      S16=C2*S16
      IF(ABS(S16-C2*S8) .LE. EPS*(1+ABS(S16))) THEN
       H=H+S16
       IF(BB .NE. B) GO TO 1
      ELSE
       BB=C1
       IF(1+CONST*ABS(C2) .NE. 1) GO TO 2
*mb    H=0
C--    MB extension
       WRITE(lunerr1,'(/'' dqcUAgauss: too high accuracy required'','//
     +               '  '' ---> STOP'')')
       STOP
*mb    ierr = 1
C--    End of MB extension
       goto 99
      END IF

   99 dqcUAgauss=H

      RETURN
      END

C     ========================================
      double precision function 
     +   dqcUBgauss(idk,bfun,yi,ti,nf,a,b,del)
C     ========================================
 
C--   Integrate singular piece (B) of splitting or coeffient function.
C--
C--   I(yi) = Int_a^b dz Bbar(yi-z)[B1(z)-B1(yi)]
C--   with Bbar(z) = exp(-z)B(exp(-z)).
C--
C--   idk        = spline order iord-1: 1=linear and 2=quadratic
C--   bfun       = coefficient function versus x = exp(-y)
C--   yi         = upper limit of convolution (exp(-yi) is argument of bfun)
C--   ti         = t-value (exp(t) is argument of bfun)      
C--   nf         = number of flavors (argument of bfun) 
C--   a          = lower limit of integral = 0.D0
C--   b          = upper limit of integral = min[yi,iord*del]
C--   del        = grid spacing
C--
C--   Modified Cernlib routine DGAUSS D103.
C--   Taken from /afs/cern.ch/asis/share/cern/2001/src/...
C--           .../mathlib/gen/d/gauss64.F and gausscod.inc

      implicit double precision (a-h,o-z)

C--   MB extension
      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
      external bfun
C--   End of MB extension

      DIMENSION W(12),X(12)

      PARAMETER (Z1 = 1, HF = Z1/2, CST = 5*Z1/1000)

      DATA X( 1) /9.6028985649753623D-1/, W( 1) /1.0122853629037626D-1/
      DATA X( 2) /7.9666647741362674D-1/, W( 2) /2.2238103445337447D-1/
      DATA X( 3) /5.2553240991632899D-1/, W( 3) /3.1370664587788729D-1/
      DATA X( 4) /1.8343464249564980D-1/, W( 4) /3.6268378337836198D-1/
      DATA X( 5) /9.8940093499164993D-1/, W( 5) /2.7152459411754095D-2/
      DATA X( 6) /9.4457502307323258D-1/, W( 6) /6.2253523938647893D-2/
      DATA X( 7) /8.6563120238783174D-1/, W( 7) /9.5158511682492785D-2/
      DATA X( 8) /7.5540440835500303D-1/, W( 8) /1.2462897125553387D-1/
      DATA X( 9) /6.1787624440264375D-1/, W( 9) /1.4959598881657673D-1/
      DATA X(10) /4.5801677765722739D-1/, W(10) /1.6915651939500254D-1/
      DATA X(11) /2.8160355077925891D-1/, W(11) /1.8260341504492359D-1/
      DATA X(12) /9.5012509837637440D-2/, W(12) /1.8945061045506850D-1/

C--   MB extension
C--   Define function to integrate
      f(u) = (dqcBsplyy(idk,1,u/del) - dqcBsplyy(idk,1,yi/del)) * 
     +  exp(-(yi-u)) * bfun(exp(-(yi-u)),exp(ti),nf)
C--   Check limits
      if(b.le.a) then
        dqcUBGauss = 0.D0
        return
      endif
      eps  = gepsi6
C--   End of MB extension
 
      H=0
      IF(B .EQ. A) GO TO 99
      CONST=CST/ABS(B-A)
      BB=A
    1 AA=BB
      BB=B
    2 C1=HF*(BB+AA)
      C2=HF*(BB-AA)
      S8=0
      DO 3 I = 1,4
      U=C2*X(I)
    3 S8=S8+W(I)*(F(C1+U)+F(C1-U))
      S16=0
      DO 4 I = 5,12
      U=C2*X(I)
    4 S16=S16+W(I)*(F(C1+U)+F(C1-U))
      S16=C2*S16
      IF(ABS(S16-C2*S8) .LE. EPS*(1+ABS(S16))) THEN
       H=H+S16
       IF(BB .NE. B) GO TO 1
      ELSE
       BB=C1
       IF(1+CONST*ABS(C2) .NE. 1) GO TO 2
*mb    H=0
C--    MB extension
       WRITE(lunerr1,'(/'' dqcUBgauss: too high accuracy required'','//
     +         '  '' ---> STOP'')')
       STOP
*mb    ierr = 1
C--    End of MB extension
       goto 99
      END IF

   99 dqcUBgauss=H

      RETURN
      END

C     =============================================
      double precision function 
     +   dqcURSgaus(idk,rfun,sfun,yi,ti,nf,a,b,del)
C     =============================================
 
C--   Integrate singular*regular piece (R*S) of splitting function.
C--
C--   I(yi) = Int_a^b dz Sbar(yi-z)[Rbar(yi-z)B1(z)-Rbar(0)B1(y)]
C--   with Sbar(z) = exp(-z)S(exp(-z)) and Rbar(z) = R(exp(-z)).
C--
C--   idk  = spline order iord-1: 1=linear and 2=quadratic
C--   rfun = coefficient function versus x = exp(-y)
C--   sfun = coefficient function versus x = exp(-y)
C--   yi   = upper limit of convolution (exp(-yi) is argument of rfun)
C--   ti   = t-value (exp(t) is argument of rfun)      
C--   nf   = number of flavors (argument of rfun) 
C--   a    = lower limit of integral = 0.D0
C--   b    = upper limit of integral = min[yi,iord*del]
C--   del  = grid spacing
C--
C--   Modified Cernlib routine DGAUSS D103.
C--   Taken from /afs/cern.ch/asis/share/cern/2001/src/...
C--           .../mathlib/gen/d/gauss64.F and gausscod.inc

      implicit double precision (a-h,o-z)

C--   MB extension
      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
      external rfun,sfun
C--   End of MB extension

      DIMENSION W(12),X(12)

      PARAMETER (Z1 = 1, HF = Z1/2, CST = 5*Z1/1000)

      DATA X( 1) /9.6028985649753623D-1/, W( 1) /1.0122853629037626D-1/
      DATA X( 2) /7.9666647741362674D-1/, W( 2) /2.2238103445337447D-1/
      DATA X( 3) /5.2553240991632899D-1/, W( 3) /3.1370664587788729D-1/
      DATA X( 4) /1.8343464249564980D-1/, W( 4) /3.6268378337836198D-1/
      DATA X( 5) /9.8940093499164993D-1/, W( 5) /2.7152459411754095D-2/
      DATA X( 6) /9.4457502307323258D-1/, W( 6) /6.2253523938647893D-2/
      DATA X( 7) /8.6563120238783174D-1/, W( 7) /9.5158511682492785D-2/
      DATA X( 8) /7.5540440835500303D-1/, W( 8) /1.2462897125553387D-1/
      DATA X( 9) /6.1787624440264375D-1/, W( 9) /1.4959598881657673D-1/
      DATA X(10) /4.5801677765722739D-1/, W(10) /1.6915651939500254D-1/
      DATA X(11) /2.8160355077925891D-1/, W(11) /1.8260341504492359D-1/
      DATA X(12) /9.5012509837637440D-2/, W(12) /1.8945061045506850D-1/

C--   MB extension
C--   Define function to integrate
      f(u) = 
     +  (rfun(exp(-(yi-u)),exp(ti),nf)*dqcBsplyy(idk,1,u/del) -
     +   rfun(1.D0,exp(ti),nf)*dqcBsplyy(idk,1,yi/del)) * 
     +   exp(-(yi-u)) * sfun(exp(-(yi-u)),exp(ti),nf)
C--   Check limits
      if(b.le.a) then
        dqcURSgaus = 0.D0
        return
      endif
      eps  = gepsi6
C--   End of MB extension
 
      H=0
      IF(B .EQ. A) GO TO 99
      CONST=CST/ABS(B-A)
      BB=A
    1 AA=BB
      BB=B
    2 C1=HF*(BB+AA)
      C2=HF*(BB-AA)
      S8=0
      DO 3 I = 1,4
      U=C2*X(I)
    3 S8=S8+W(I)*(F(C1+U)+F(C1-U))
      S16=0
      DO 4 I = 5,12
      U=C2*X(I)
    4 S16=S16+W(I)*(F(C1+U)+F(C1-U))
      S16=C2*S16
      IF(ABS(S16-C2*S8) .LE. EPS*(1+ABS(S16))) THEN
       H=H+S16
       IF(BB .NE. B) GO TO 1
      ELSE
       BB=C1
       IF(1+CONST*ABS(C2) .NE. 1) GO TO 2
*mb    H=0
C--    MB extension
       WRITE(lunerr1,'(/'' dqcURSgaus: too high accuracy required'','//
     +         '  '' ---> STOP'')')
       STOP
*mb    ierr = 1
C--    End of MB extension
       goto 99
      END IF

   99 dqcURSgaus=H

      RETURN
      END
      
C     ====================================================
      double precision function dqcUXgauss(idk,yi,a,b,del)
C     ====================================================

C--   Integrate product of B-spline functions for convolution FcrossF.
C--
C--   I(yi) = Int_a^b dz B1(z) B1(yi-z).
C--
C--   idk        = spline order iord-1: 1=linear and 2=quadratic
C--   yi         = upper limit of convolution       
C--   a          = lower limit of integral = 0.D0
C--   b          = upper limit of integral = min[yi,iord*del]
C--   del        = grid spacing
C--
C--   Modified Cernlib routine DGAUSS D103.
C--   Taken from /afs/cern.ch/asis/share/cern/2001/src/...
C--           .../mathlib/gen/d/gauss64.F and gausscod.inc

      implicit double precision (a-h,o-z)

C--   MB extension
      include 'qcdnum.inc'
      include 'qluns1.inc'
      include 'qpars6.inc'
C--   End of MB extension

      DIMENSION W(12),X(12)

      PARAMETER (Z1 = 1, HF = Z1/2, CST = 5*Z1/1000)

      DATA X( 1) /9.6028985649753623D-1/, W( 1) /1.0122853629037626D-1/
      DATA X( 2) /7.9666647741362674D-1/, W( 2) /2.2238103445337447D-1/
      DATA X( 3) /5.2553240991632899D-1/, W( 3) /3.1370664587788729D-1/
      DATA X( 4) /1.8343464249564980D-1/, W( 4) /3.6268378337836198D-1/
      DATA X( 5) /9.8940093499164993D-1/, W( 5) /2.7152459411754095D-2/
      DATA X( 6) /9.4457502307323258D-1/, W( 6) /6.2253523938647893D-2/
      DATA X( 7) /8.6563120238783174D-1/, W( 7) /9.5158511682492785D-2/
      DATA X( 8) /7.5540440835500303D-1/, W( 8) /1.2462897125553387D-1/
      DATA X( 9) /6.1787624440264375D-1/, W( 9) /1.4959598881657673D-1/
      DATA X(10) /4.5801677765722739D-1/, W(10) /1.6915651939500254D-1/
      DATA X(11) /2.8160355077925891D-1/, W(11) /1.8260341504492359D-1/
      DATA X(12) /9.5012509837637440D-2/, W(12) /1.8945061045506850D-1/

C--   MB extension
C--   Define function to integrate
      f(u) = dqcBsplyy(idk,1,u/del) * dqcBsplyy(idk,1,(yi-u)/del)
C--   Check limits
      if(b.le.a) then
        dqcUXGauss = 0.D0
        return
      endif
      eps  = gepsi6
C--   End of MB extension

      H=0
      IF(B .EQ. A) GO TO 99
      CONST=CST/ABS(B-A)
      BB=A
    1 AA=BB
      BB=B
    2 C1=HF*(BB+AA)
      C2=HF*(BB-AA)
      S8=0
      DO 3 I = 1,4
      U=C2*X(I)
    3 S8=S8+W(I)*(F(C1+U)+F(C1-U))
      S16=0
      DO 4 I = 5,12
      U=C2*X(I)
    4 S16=S16+W(I)*(F(C1+U)+F(C1-U))
      S16=C2*S16
      IF(ABS(S16-C2*S8) .LE. EPS*(1+ABS(S16))) THEN
       H=H+S16
       IF(BB .NE. B) GO TO 1
      ELSE
       BB=C1
       IF(1+CONST*ABS(C2) .NE. 1) GO TO 2
*mb    H=0
C--    MB extension
       WRITE(lunerr1,'(/'' dqcUXgauss: too high accuracy required'','//
     +               '  '' ---> STOP'')')
       STOP
*mb    ierr = 1
C--    End of MB extension
       goto 99
      END IF

   99 dqcUXgauss=H

      RETURN
      END
       

