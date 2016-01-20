
C--   This is the file mchar.f containing the MBUTIL character manipulation routines 
C--
C--   subroutine smb_cfill(char,cstring)
C--   subroutine smb_cleft(cstring)
C--   subroutine smb_crght(cstring)
C--   subroutine smb_cltou(cstring)
C--   subroutine smb_cutol(cstring)
C--   integer function imb_lenoc(cstring)
C--   integer function imb_frstc(cstring)
C--   logical function lmb_compc(string1,string2,n1,n2)
C--   logical function lmb_match(string,substr,wdcard)
C--   subroutine smb_itoch(in,chout,leng)
      
C     ==================================
      subroutine smb_cfill(char,cstring)
C     ==================================

C--   Input:  character string char.
C--           character string cstring.
C--           On exit all characters in cstring will be set to char(1:1).
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      character*(*) char, cstring

      l = len(cstring)

      do i = 1,l
        cstring(i:i) = char(1:1)
      enddo

      return
      end

C     =============================
      subroutine smb_cleft(cstring)
C     =============================

C--   Left adjust character string cstring.
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      character*(*) cstring

      max         = len(cstring)
      if(max.le.0)  return
      i1          = imb_frstc(cstring)
      i2          = imb_lenoc(cstring)
      k           = 0
      do i = i1,i2
        k = k+1
        cstring(k:k) = cstring(i:i)
      enddo
      do i = k+1,max
        cstring(i:i) = ' '
      enddo

      return
      end

C     =============================
      subroutine smb_crght(cstring)
C     =============================

C--   Right adjust character string cstring.
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      character*(*) cstring

      max         = len(cstring)
      if(max.le.0)  return
      i1          = imb_frstc(cstring)
      i2          = imb_lenoc(cstring)
      k           = max+1
      do i = i2,i1,-1
        k = k-1
        cstring(k:k) = cstring(i:i)
      enddo
      do i = k-1,1,-1
        cstring(i:i) = ' '
      enddo

      return
      end

C     =============================
      subroutine smb_cltou(cstring)
C     =============================

C--   Input: character string cstring.
C--          On exit all characters in cstring will be set
C--          to upper case.
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      character*(*) cstring
      character*26  charl, charu
      data charl   /'abcdefghijklmnopqrstuvwxyz'/
      data charu   /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      j          = len(cstring)

      do i = 1,j
        do k = 1,26
          if(cstring(i:i).eq.charl(k:k)) then
            cstring(i:i) = charu(k:k)
          endif
        enddo
      enddo

      return
      end

C     =============================
      subroutine smb_cutol(cstring)
C     =============================

C--   Input: character string cstring.
C--          On exit all characters in cstring will be set
C--          to lower case.
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      character*(*) cstring
      character*26  charl, charu
      data charl   /'abcdefghijklmnopqrstuvwxyz'/
      data charu   /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      j          = len(cstring)

      do i = 1,j
        do k = 1,26
          if(cstring(i:i).eq.charu(k:k)) then
            cstring(i:i) = charl(k:k)
          endif
        enddo
      enddo

      return
      end

C     ===================================
      integer function imb_lenoc(cstring)
C     ===================================

C--   Input:  cstring       Character string.
C--   Output: imb_lenoc     Position of the last non-blank
C--                         character in cstring.
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      character*(*) cstring

      j         = len(cstring)
      imb_lenoc = 0

      do i = j,1,-1
        if(cstring(i:i).ne.' ') then
          imb_lenoc = i
          return
        endif
      enddo

      return
      end

C     ===================================
      integer function imb_frstc(cstring)
C     ===================================

C--   Input:  cstring       Character string.
C--   Output: imb_frstc     Position of the first non-blank
C--                         character in cstring.
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      character*(*) cstring

      j         = len(cstring)
      imb_frstc = 0

      do i = 1,j
        if(cstring(i:i).ne.' ') then
          imb_frstc = i
          return
        endif
      enddo

      return
      end

C     =================================================
      logical function lmb_compc(string1,string2,n1,n2)
C     =================================================

C--   Do a case independent comparision of the characters
C--   string1(n1:n2) and string2(n1:n2).
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      character*(*) string1, string2
      character*1   ch1    , ch2

      lmb_compc = .false.
      if((n1.le.0).or.(n2.le.0).or.(n2.lt.n1)) return
      len1 = imb_lenoc(string1)
      if(len1.lt.n2)                           return
      len2 = imb_lenoc(string2)
      if(len2.lt.n2)                           return

      lmb_compc = .true.
      do i = n1,n2
        ch1 = string1(i:i)
        ch2 = string2(i:i)
        call smb_cltou(ch1)
        call smb_cltou(ch2)
        if(ch1.ne.ch2) goto 900
      enddo

      return

 900  continue

      lmb_compc = .false.

      return
      end

C     ================================================
      logical function lmb_match(string,substr,wdcard)
C     ================================================

C--   Translates string, substr to upper case and returns
C--   .true. if substr is contained in string, .false. otherwise. 
C--   If string and/or substr are null strings, lmb_match = .FALSE.
C--   Leading blanks in substr are stripped off.
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      character*(*) string,substr,wdcard
      character*80  str,sub,sbb
      character*1   wcd

      lmb_match = .false.
      ipos      = 0

      len1 = imb_lenoc(string)
      if(len1.eq.0.or.len1.gt.80)   return
      len2 = imb_lenoc(substr)
      if(len2.eq.0.or.len2.gt.80)   return

C--   Avoid modifying the input argument(s)
      call smb_cfill(' ',str)
      call smb_cfill(' ',sub)
      str(1:len1) = string(1:len1)
      sub(1:len2) = substr(1:len2)
      wcd         = wdcard(1:1)
      call smb_cltou(str)
      call smb_cltou(sub)
      call smb_cltou(wcd)
      call smb_cleft(sub)
      len2 = imb_lenoc(sub)
      if(len2.gt.len1) return

      i2 = len1-len2+1
      do i = 1,i2
        iend = i-1+len2
        do k = i,iend
          l = k-i+1
          sbb(l:l) = sub(l:l)
          if(wcd.ne.' '.and.sub(l:l).eq.wcd) sbb(l:l) = str(k:k)
        enddo
        if(str(i:iend).eq.sbb(1:len2)) lmb_match = .true.
      enddo

      return
      end
      
C     ===================================
      subroutine smb_itoch(in,chout,leng)
C     ===================================

C--   Format integer as character string with length leng
C--
C--   in      (in)     input integer
C--   chout   (out)    character string
C--   leng    (out)    length of character string
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      implicit double precision(a-h,o-z)
      
      character chout*(*)
      character f1*2, f2*10, f3*1
      parameter(f1 = '(I')
      parameter(f3 = ')')
      
      call smb_cfill(' ',chout)
      call smb_cfill(' ',f2)
      lmax = len(chout)
      din  = in 

      if(in.eq.0) then
        leng = 1
      elseif(in.gt.0) then
        leng = int(log10(din)+1)
      elseif(in.lt.0) then
        leng = int(log10(abs(din))+2)
      endif
      
      if(leng.gt.lmax) then
         call smb_cfill('*',chout)
         leng = lmax
         return
      else
         write(f2,'(I10)') leng
         lf2 = imb_lenoc(f2)
      endif

      write(unit=chout,fmt=f1//f2(1:lf2)//f3) in
      
      return
      end                        
      
