
C--   This is the file mbits.f containing the MBUTIL bit manipulation routines
C--
C--   subroutine smb_sbit1(i,n)
C--   subroutine smb_sbit0(i,n)
C--   integer function imb_gbitn(i,n)
C--   integer function imb_sbits(cpatt)
C--   subroutine smb_gbits(i,cpatt)
C--   integer function imb_test0(mask,i)
C--   integer function imb_test1(mask,i)

C     =========================
      subroutine smb_sbit1(i,n)
C     =========================

C--   Set bit n of 32 bit integer i to 1.
C--
C--   n =  1: set LSB i.e. rightmost bit.
C--   n = 32: set MSB i.e. leftmost  bit.
C--
C--   Assumes that the 32 bit integer representation is such that
C--   that i = 1 has all bits set to zero except LSB = 1, i.e.:
C--
C--           33322222222221111111111000000000
C--           21098765432109876543210987654321
C--   i = 1 = 00000000000000000000000000000001
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      data i1 /1/

      i = ior(i,ishft(i1,n-1))

      return
      end

C     =========================
      subroutine smb_sbit0(i,n)
C     =========================

C--   Set bit n of 32 bit integer i to 0.
C--
C--   n =  1: set LSB i.e. rightmost bit.
C--   n = 32: set MSB i.e. leftmost  bit.
C--
C--   Assumes that the 32 bit integer representation is such that
C--   that i = 1 has all bits set to zero except LSB = 1, i.e.:
C--
C--           33322222222221111111111000000000
C--           21098765432109876543210987654321
C--   i = 1 = 00000000000000000000000000000001
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      data i1 /1/

      i = iand(i,not(ishft(i1,n-1)))

      return
      end

C     ===============================
      integer function imb_gbitn(i,n)
C     ===============================

C--   Get bit n of 32 bit integer i.
C--
C--   n =  1: get LSB i.e. rightmost bit.
C--   n = 32: get MSB i.e. leftmost  bit.
C--
C--   Assumes that the 32 bit integer representation is such that
C--   that i = 1 has all bits set to zero except LSB = 1, i.e.:
C--
C--           33322222222221111111111000000000
C--           21098765432109876543210987654321
C--   i = 1 = 00000000000000000000000000000001
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      data i1 /1/

      imb_gbitn = ishft(iand(i,ishft(i1,n-1)),1-n)

      return
      end

C     =================================
      integer function imb_sbits(cpatt)
C     =================================

C--   Set bitpattern of 32 bit integer i
C--
C--   Input:  character*32 cpatt containing the bitpattern.
C--   Output: integer i (32 bits) with bits set according to cpatt.
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      character*32 cpatt

      do j = 1,32
        k = 33-j
        if(cpatt(k:k).eq.'0') then
          call smb_sbit0(imb_sbits,j)
        else
          call smb_sbit1(imb_sbits,j)
        endif
      enddo

      return
      end

C     =============================
      subroutine smb_gbits(i,cpatt)
C     =============================

C--   Get bitpattern of 32 bit integer i
C--
C--   Input:  integer i (32 bits).
C--   Output: character*32 cpatt containing the bitpattern.
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      character*32 cpatt

      do j = 1,32
        k = 33-j
        l = imb_gbitn(i,j)
        if(l.eq.0) then
          cpatt(k:k) = '0'
        else
          cpatt(k:k) = '1'
        endif
      enddo

      return
      end

C     ==================================
      integer function imb_test0(mask,i)
C     ==================================

C--   Test pattern of 'zero' bits in 32 bit integer i
C--
C--   Input:  integer mask with bit n = 0(1) -> ignore(test) bit n of i.
C--           integer i containing the bitpattern to be tested.
C--   Output: integer imb_test0 = 0 if all tested bits of i are 0;
C--                             # 0 if not all tested bits of i are 0.
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      imb_test0 = iand(mask,i)

      return
      end

C     ==================================
      integer function imb_test1(mask,i)
C     ==================================

C--   Test pattern of 'one' bits in 32 bit integer i
C--
C--   Input:  integer mask with bit n = 0(1) -> ignore(test) bit n of i.
C--           integer i containing the bitpattern to be tested.
C--   Output: integer imb_test1 = 0 if all tested bits of i are 1;
C--                             # 0 if not all tested bits of i are 1.
C--
C--   Author: Michiel Botje h24@nikhef.nl   13-04-09

      imb_test1 = iand(mask,not(i))

      return
      end
