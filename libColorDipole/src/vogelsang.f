
C--   This is the file vogelsang.f with the polarised splitting functions
C--
C--   double precision function dqcDPQQ0A(x,qmu2,nf)  !Vogelsamg Eq. (28)
C--   double precision function dqcDPQQ0B(x,qmu2,nf)  !Vogelsang Eq. (28)
C--   double precision function dqcDPQQ0D(x,qmu2,nf)  !Vogelsang Eq. (28)
C--   double precision function dqcDPQG0A(x,qmu2,nf)  !Vogelsang Eq. (29)
C--   double precision function dqcDPGQ0A(x,qmu2,nf)  !Vogelsang Eq. (30)
C--   double precision function dqcDPGG0A(x,qmu2,nf)  !Vogelsang Eq. (31)
C--   double precision function dqcDPGG0B(x,qmu2,nf)  !Vogelsang Eq. (31)
C--   double precision function dqcDPGG0D(x,qmu2,nf)  !Vogelsang Eq. (31)
C--
C--   double precision function dqcDPPL1B(x,qmu2,nf)  !Vogelsang Eq. (44)
C--   double precision function dqcDPMI1A(x,qmu2,nf)  !Vogelsang Eq. (44)
C--   double precision function dqcDPMI1B(x,qmu2,nf)  !Vogelsang Eq. (44)
C--   double precision function dqcDPQS1A(x,qmu2,nf)  !Vogelsang Eq. (45)
C--   double precision function delPQG(x)             !Vogelsang Eq. (43)
C--   double precision function EStwo(x)              !Vogelsang Eq. (52)
C--   double precision function dqcDPQG1A(x,qmu2,nf)  !Vogelsang Eq. (46)
C--   double precision function delPGQ(x)             !Vogelsang Eq. (43)
C--   double precision function dqcDPGQ1A(x,qmu2,nf)  !Vogelsang Eq. (47)
C--   double precision function delPGGA(x)            !Vogelsang Eq. (43)
C--   double precision function delPGGB(x)            !Vogelsang Eq. (43)
C--   double precision function dqcDPGG1A(x,qmu2,nf)  !Vogelsang Eq. (48)
C--   double precision function dqcDPGG1R(x,qmu2,nf)  !Vogelsang Eq. (48)
C--   double precision function dqcDPGG1S(x,qmu2,nf)  !Vogelsang Eq. (48)
C--   double precision function dqcDPGG1D(x,qmu2,nf)  !Vogelsang Eq. (48)
C--
C--   LO and NLO polarised splitting functions typed in from
C--   W. Vogelsang, Nucl. Phys. B475 (1996) 47, hep-ph/9603366
C--
C--   All splitting functions are given x, qmu2, and nf as arguments,
C--   some of which are dummy (w/dummy assignment in body of function to
C--   avoid compiler warnings)
C-- 
C--   Naming convention:  dqcDPxxxA()   regular piece
C--                       dqcDPxxxB()   singular piece
C--                       dqcDPxxxR()   regular piece of RS product
C--                       dqcDPxxxS()   singular piece of RS product
C--                       dqcDPxxxD()   multiplies delta(1-x)

C------------------------------------------------------------------------

C--   LO polarised splitting functions

C     ==============================================
      double precision function dqcDPQQ0A(x,qmu2,nf)  !Vogelsamg Eq. (28)
C     ==============================================

      implicit double precision (a-h,o-z)
      
      parameter (CF = 4.D0/3.D0)
      
      ddum   = qmu2
      idum   = nf
      dqcDPQQ0A = CF * (-1.D0-x)
      
      return
      end
      
C     ==============================================
      double precision function dqcDPQQ0B(x,qmu2,nf)  !Vogelsang Eq. (28) 
C     ==============================================

      implicit double precision (a-h,o-z)
      
      parameter (CF = 4.D0/3.D0)

      ddum   = qmu2
      idum   = nf
      dqcDPQQ0B = CF * 2.D0 / (1.D0-x)
      
      return
      end  
      
C     ==============================================
      double precision function dqcDPQQ0D(x,qmu2,nf)  !Vogelsang Eq. (28)
C     ==============================================

      implicit double precision (a-h,o-z)
      
      parameter (CF = 4.D0/3.D0)

      ddum   = x
      ddum   = qmu2
      idum   = nf
      dqcDPQQ0D = CF * 3.D0 / 2.D0
      
      return
      end
      
C     ==============================================
      double precision function dqcDPQG0A(x,qmu2,nf)  !Vogelsang Eq. (29)
C     ==============================================

      implicit double precision (a-h,o-z)
            
      parameter (TR = 1.D0/2.D0)
      
      ddum   = qmu2
      TF     = TR*nf
      dqcDPQG0A = 2.D0 * TF * (2.D0*x-1.D0)
      
      return
      end
      
C     ==============================================
      double precision function dqcDPGQ0A(x,qmu2,nf)  !Vogelsang Eq. (30)
C     ==============================================

      implicit double precision (a-h,o-z)
      
      parameter (CF = 4.D0/3.D0)

      ddum   = qmu2
      idum   = nf
      dqcDPGQ0A = CF * (2.D0-x)
      
      return
      end
      
C     ==============================================
      double precision function dqcDPGG0A(x,qmu2,nf)  !Vogelsang Eq. (31)
C     ==============================================

      implicit double precision (a-h,o-z)
      
      parameter (NC = 3)

      ddum   = qmu2
      idum   = nf
      dqcDPGG0A = 2.D0 * NC * (-2.D0*x+1.D0)
      
      return
      end
      
C     ==============================================
      double precision function dqcDPGG0B(x,qmu2,nf)  !Vogelsang Eq. (31)
C     ==============================================

      implicit double precision (a-h,o-z)
      
      parameter (NC = 3)

      ddum   = qmu2
      idum   = nf
      dqcDPGG0B = 2.D0 * NC / (1.D0-x)
      
      return
      end  
      
C     ==============================================
      double precision function dqcDPGG0D(x,qmu2,nf)  !Vogelsang Eq. (31)
C     ==============================================

      implicit double precision (a-h,o-z)
      
      parameter (NC = 3)
      parameter (TR = 1.D0/2.D0)
      
      ddum = x
      ddum = qmu2
      TF   = TR*nf
      B0   = NC*11.D0/3.D0 - TF*4.D0/3.D0

      dqcDPGG0D = B0/2.D0 
      
      return
      end

C------------------------------------------------------------------------

C--   NLO polarised splitting functions

C     ==============================================
      double precision function dqcDPPL1B(x,qmu2,nf)  !Vogelsang Eq. (44)
C     ==============================================

      implicit double precision (a-h,o-z)
      
      dqcDPPL1B = dqcPMI1B(x,qmu2,nf)
      
      return
      end
      

C     ==============================================
      double precision function dqcDPMI1A(x,qmu2,nf)  !Vogelsang Eq. (44)
C     ==============================================

      implicit double precision (a-h,o-z)
      
      dqcDPMI1A = dqcPPL1A(x,qmu2,nf)
      
      return
      end     
      
C     ==============================================
      double precision function dqcDPMI1B(x,qmu2,nf)  !Vogelsang Eq. (44)
C     ==============================================

      implicit double precision (a-h,o-z)
      
      dqcDPMI1B = dqcPPL1B(x,qmu2,nf)
      
      return
      end                                     

C     ==============================================
      double precision function dqcDPQS1A(x,qmu2,nf)  !Vogelsang Eq. (45)
C     ==============================================

      implicit double precision (a-h,o-z)
      
      parameter (CF = 4.D0/3.D0)
      parameter (TR = 1.D0/2.D0)
      
      ddum = qmu2
      TF   = TR*nf
      dlx  = log(x) 
      
      dqcDPQS1A = 2.D0 * CF * TF * 
     &          ( 1.D0-x - (1.D0-3.D0*x)*dlx - (1.D0+x)*dlx*dlx )
     
      return
      end
      
C     ===================================      
      double precision function delPQG(x)  !Vogelsang Eq. (43)
C     ===================================

      implicit double precision (a-h,o-z)
      
      delPQG = 2.D0*x - 1.D0
      
      return
      end
      
C     ==================================   
      double precision function EStwo(x)  !Vogelsang Eq. (52)
C     ==================================

      implicit double precision (a-h,o-z)
      
      parameter (PI = 3.14159265)
      
      dlx    = log(x)
      dlopx  = log(1.D0+x)
      EStwo  = -2.D0*dmb_dilog(-x) - 2.D0*dlx*dlopx +
     &          dlx*dlx/2.D0 - PI*PI/6.D0 
      
      return
      end                        
      
C     ==============================================      
      double precision function dqcDPQG1A(x,qmu2,nf)  !Vogelsang Eq. (46)
C     ==============================================

      implicit double precision (a-h,o-z)
      
      parameter (CF = 4.D0/3.D0)
      parameter (TR = 1.D0/2.D0)
      parameter (NC = 3)
      parameter (PI = 3.14159265)
      
      ddum  = qmu2
      TF    = TR*nf
      dlx   = log(x)
      omx   = 1.D0-x
      dlomx = log(omx) 
      
      dqcDPQG1A = 
     &     CF * TF * 
     &    ( -22.D0 + 27.D0*x - 9.D0*dlx + 8.D0*omx*dlomx +
     &      delPQG(x) * ( 2.D0*dlomx*dlomx - 4.D0*dlomx*dlx +
     &                    dlx*dlx - 2.D0*PI*PI/3.D0 ) ) 
     &   + NC * TF *
     &    ( 2.D0*(12.D0-11.D0*x) - 8.D0*omx*dlomx + 
     &      2.D0*(1.D0+8.D0*x)*dlx -
     &      2.D0*(dlomx*dlomx - PI*PI/6.D0)*delPQG(x) -
     &     (2.D0*EStwo(x) - 3.D0*dlx*dlx)*delPQG(-x) )  
      
      return
      end
      
C     ===================================      
      double precision function delPGQ(x)  !Vogelsang Eq. (43)
C     ===================================

      implicit double precision (a-h,o-z)
      
      delPGQ = 2.D0 - x
      
      return
      end
      
C     ==============================================      
      double precision function dqcDPGQ1A(x,qmu2,nf)  !Vogelsang Eq. (47)
C     ==============================================

      implicit double precision (a-h,o-z)
      
      parameter (CF = 4.D0/3.D0)
      parameter (TR = 1.D0/2.D0)
      parameter (NC = 3)
      parameter (PI = 3.14159265)
      
      ddum  = qmu2
      TF    = TR*nf
      dlx   = log(x)
      omx   = 1.D0-x
      dlomx = log(omx)
      
      dqcDPGQ1A = 
     &    CF * TF *
     &   ( -4.D0*(x+4.D0)/9.D0 - 4.D0*delPGQ(x)*dlomx/3.D0 )
     &  + CF * CF *
     &   ( -1.D0/2.D0 - (4.D0-x)*dlx/2.D0 - delPGQ(-x)*dlomx +
     &     (-4.D0 - dlomx*dlomx + dlx*dlx/2.D0)*delPGQ(x) )
     &  + CF * NC *
     &   ( (4.D0-13.D0*x)*dlx + (10.D0+x)*dlomx/3.D0 +
     &     (41.D0+35.D0*x)/9.D0 + 
     &     (-2.D0*EStwo(x) + 3.D0*dlx*dlx)*delPGQ(-x)/2.D0 +
     &     (dlomx*dlomx - 2.D0*dlomx*dlx - PI*PI/6.D0)*delPGQ(x) )
     
      return
      end
     
C     ====================================      
      double precision function delPGGA(x)  !Vogelsang Eq. (43)
C     ====================================

      implicit double precision (a-h,o-z)
      
      delPGGA = -2.D0*x + 1.D0
      
      return
      end
      
     
C     ====================================      
      double precision function delPGGB(x)  !Vogelsang Eq. (43)
C     ====================================

      implicit double precision (a-h,o-z)
      
      delPGGB = 1.D0/(1.D0-x)
      
      return
      end 
      
C     ==============================================      
      double precision function dqcDPGG1A(x,qmu2,nf)  !Vogelsang Eq. (48)
C     ==============================================

      implicit double precision (a-h,o-z)
      
      parameter (CF = 4.D0/3.D0)
      parameter (TR = 1.D0/2.D0)
      parameter (NC = 3)
      parameter (PI = 3.14159265)
      
      ddum  = qmu2
      TF    = TR*nf
      dlx   = log(x)
      omx   = 1.D0-x
      opx   = 1.D0+x
      dlomx = log(omx)
      
      dqcDPGG1A = 
     &    -NC * TF *
     &   ( 4.D0*omx + 4.D0*opx*dlx/3.D0 + 
     &     20.D0*delPGGA(x)/9.D0 )
     &  -  CF * TF *
     &   ( 10.D0*omx + 2.D0*(5.D0-x)*dlx + 2.D0*opx*dlx*dlx )
     &  +  NC * NC *
     &   ( (29.D0-67.D0*x)*dlx/3.D0 - 19.D0*omx/2.D0 +
     &     4.D0*opx*dlx*dlx - 
     &     2.D0*EStwo(x)*(delPGGA(-x)+delPGGB(-x)) +
     &     (67.D0/9.D0 - 4.D0*dlomx*dlx + dlx*dlx - 
     &      PI*PI/3.D0)*delPGGA(x) )
     
      return
      end

C     ==============================================      
      double precision function dqcDPGG1R(x,qmu2,nf)  !Vogelsang Eq. (48)
C     ==============================================

      implicit double precision (a-h,o-z)
      
      parameter (CF = 4.D0/3.D0)
      parameter (TR = 1.D0/2.D0)
      parameter (NC = 3)
      parameter (PI = 3.14159265)
      
      ddum  = qmu2
      TF    = TR*nf
      dlx   = log(x)
      omx   = 1.D0-x
      dlomx = log(omx)
C--   Set ln(x)ln(1-x) = 0 for x = 1       
      if(x.eq.1.D0) then
        product = 0.D0
      else
        product = dlomx*dlx
      endif    
      
      dqcDPGG1R = -NC * TF * 20.D0 / 9.D0 +
     &             NC * NC * ( 67.D0/9.D0 - 4.D0*product + 
     &             dlx*dlx - PI*PI/3.D0 )
     
      return
      end
      
C     ==============================================      
      double precision function dqcDPGG1S(x,qmu2,nf)  !Vogelsang Eq. (48)
C     ==============================================

      implicit double precision (a-h,o-z)
      
      ddum      = qmu2
      idum      = nf
      dqcDPGG1S = delPGGB(x)
      
      return
      end

C     ==============================================      
      double precision function dqcDPGG1D(x,qmu2,nf)  !Vogelsang Eq. (48)
C     ==============================================

      implicit double precision (a-h,o-z)
      
      parameter (CF = 4.D0/3.D0)
      parameter (TR = 1.D0/2.D0)
      parameter (NC = 3)
      parameter (PI = 3.14159265)
      parameter (Z3 = 1.202057)
      
      ddum  = x
      ddum  = qmu2
      TF    = TR*nf
      
      dqcDPGG1D = -NC*TF*4.D0/3.D0 - CF*TF + 
     &             NC*NC*(3.D0*Z3 + 8.D0/3.D0)
     
      return
      end      
