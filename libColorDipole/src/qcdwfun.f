
C--   This is qcdwfun.f containing the splitting function interfaces
C--   The splitting functions themselves are in the pij subdirectory
C--
C--   double precision function dqcAchi(qmu2)
C--
C--   LO spacelike      
C--   --------------------------------------------- 
C--   double precision function dqcPQQ0R(x,qmu2,nf)
C--   double precision function dqcPQQ0S(x,qmu2,nf)
C--   double precision function dqcPQQ0D(x,qmu2,nf)
C--   double precision function dqcPQG0A(x,qmu2,nf)
C--   double precision function dqcPGQ0A(x,qmu2,nf)
C--   double precision function dqcPGG0A(x,qmu2,nf)
C--   double precision function dqcPGG0R(x,qmu2,nf)
C--   double precision function dqcPGG0S(x,qmu2,nf)
C--   double precision function dqcPGG0D(x,qmu2,nf)
C--
C--   LO timelike  (swap QG and GQ)
C--   --------------------------------------------- 
C--   double precision function dqcTQG0A(x,qmu2,nf)
C--   double precision function dqcTGQ0A(x,qmu2,nf)
C--
C--   NLO spacelike
C--   ---------------------------------------------  
C--   double precision function dqcPPL1A(x,qmu2,nf)
C--   double precision function dqcPPL1B(x,qmu2,nf)
C--   double precision function dqcPMI1B(x,qmu2,nf)
C--   double precision function dqcPQQ1A(x,qmu2,nf)
C--   double precision function dqcPQQ1B(x,qmu2,nf)
C--   double precision function dqcPQG1A(x,qmu2,nf)
C--   double precision function dqcPGQ1A(x,qmu2,nf)
C--   double precision function dqcPGG1A(x,qmu2,nf)
C--   double precision function dqcPGG1B(x,qmu2,nf)
C--
C--   NLO timelike (fragmetation function evolution)     
C--   ---------------------------------------------- 
C--   double precision function dqcTPL1A(x,qmu2,nf)
C--   double precision function dqcTPL1B(x,qmu2,nf)
C--   double precision function dqcTMI1B(x,qmu2,nf)
C--   double precision function dqcTQQ1A(x,qmu2,nf)
C--   double precision function dqcTQQ1B(x,qmu2,nf)
C--   double precision function dqcTQG1A(x,qmu2,nf)
C--   double precision function dqcTGQ1A(x,qmu2,nf)
C--   double precision function dqcTGG1A(x,qmu2,nf)
C--   double precision function dqcTGG1B(x,qmu2,nf)
C--
C--   NNLO spacelike (no timelike available)     
C--   --------------------------------------------- 
C--   double precision function dqcPPL2A(x,qmu2,nf)
C--   double precision function dqcPPL2B(x,qmu2,nf)
C--   double precision function dqcPPL2D(x,qmu2,nf)
C--   double precision function dqcPMI2A(x,qmu2,nf)
C--   double precision function dqcPMI2B(x,qmu2,nf)
C--   double precision function dqcPMI2D(x,qmu2,nf)
C--   double precision function dqcPVA2A(x,qmu2,nf)
C--   double precision function dqcPQQ2A(x,qmu2,nf)
C--   double precision function dqcPQG2A(x,qmu2,nf) 
C--   double precision function dqcPGQ2A(x,qmu2,nf)
C--   double precision function dqcPGG2A(x,qmu2,nf)
C--   double precision function dqcPGG2B(x,qmu2,nf)
C--   double precision function dqcPGG2D(x,qmu2,nf)
C--
C--   NNLO A coefficients for discontinuities     
C--   --------------------------------------------- 
C--   double precision function dqcAGQ2A(x,qmu2,nf)
C--   double precision function dqcAGG2A(x,qmu2,nf)
C--   double precision function dqcAGG2B(x,qmu2,nf)
C--   double precision function dqcAGG2D(x,qmu2,nf)
C--   double precision function dqcAQQ2A(x,qmu2,nf)
C--   double precision function dqcAQQ2B(x,qmu2,nf)
C--   double precision function dqcAQQ2D(x,qmu2,nf)
C--   double precision function dqcAHQ2A(x,qmu2,nf)
C--   double precision function dqcAHG2A(x,qmu2,nf)
C--   double precision function dqcAHG2D(x,qmu2,nf)

C     =======================================
      double precision function dqcAchi(qmu2)
C     =======================================

      implicit double precision (a-h,o-z)
      
      qq       = qmu2    !avoid compiler warning
      dqcAchi  = 1.D0
      
      return
      end
      
C--   ---------------------------------------------      
C--   LO spacelike      
C--   ---------------------------------------------      

C     =============================================
      double precision function dqcPQQ0R(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)
      
      qq       = qmu2    !avoid compiler warning
      dqcPQQ0R = dqcP0FFR(x,nf)
      
      return
      end
      
C     =============================================
      double precision function dqcPQQ0S(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)
      
      qq       = qmu2    !avoid compiler warning
      dqcPQQ0S = dqcP0FFS(x,nf)
      
      return
      end
      
C     =============================================
      double precision function dqcPQQ0D(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)
      
      xx       = x       !avoid compiler warning
      qq       = qmu2    !avoid compiler warning
      jf       = nf      !avoid compiler warning
      dqcPQQ0D = 2.D0
      
      return
      end      
      
C     =============================================
      double precision function dqcPQG0A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)
      
      qq       = qmu2    !avoid compiler warning
      dqcPQG0A = 2.D0 * nf * dqcP0FGA(x,nf)
      
      return
      end
      
C     =============================================
      double precision function dqcPGQ0A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)
      
      qq       = qmu2    !avoid compiler warning
      dqcPGQ0A = dqcP0GFA(x,nf)
      
      return
      end
      
C     =============================================
      double precision function dqcPGG0A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)
      
      qq       = qmu2    !avoid compiler warning
      dqcPGG0A = dqcP0GGA(x,nf)
      
      return
      end 
      
C     =============================================
      double precision function dqcPGG0R(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)
      
      qq       = qmu2    !avoid compiler warning
      dqcPGG0R = dqcP0GGR(x,nf)
      
      return
      end
      
C     =============================================
      double precision function dqcPGG0S(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)
      
      qq       = qmu2    !avoid compiler warning
      dqcPGG0S = dqcP0GGS(x,nf)
      
      return
      end
      
C     =============================================
      double precision function dqcPGG0D(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)
      
      xx       = x       !avoid compiler warning
      qq       = qmu2    !avoid compiler warning
      dqcPGG0D = 6.D0*(11.D0/12.D0 - nf/18.D0)
      
      return
      end
      
C--   ---------------------------------------------      
C--   LO timelike  (swap QG and GQ)
C--   ---------------------------------------------       
      
C     =============================================
      double precision function dqcTQG0A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)
      
      qq       = qmu2    !avoid compiler warning
      dqcTQG0A = 2.D0 * nf * dqcP0GFA(x,nf)
      
      return
      end
      
C     =============================================
      double precision function dqcTGQ0A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)
      
      qq       = qmu2    !avoid compiler warning
      dqcTGQ0A = dqcP0FGA(x,nf)
      
      return
      end
      

C--   ---------------------------------------------      
C--   NLO spacelike
C--   ---------------------------------------------      
      
C     =============================================
      double precision function dqcPPL1A(x,qmu2,nf)
C     =============================================
      
      implicit double precision (a-h,o-z)
      
      qq       = qmu2    !avoid compiler warning
      dqcPPL1A = pp1sfunc(x,nf)-pm1sfunc(x,nf)
      
      return
      end
      
C     =============================================
      double precision function dqcPPL1B(x,qmu2,nf)
C     =============================================
      
      implicit double precision (a-h,o-z)
      
      qq       = qmu2    !avoid compiler warning
      dqcPPL1B = pm1sfunc(x,nf)
      
      return
      end
      
C     =============================================
      double precision function dqcPMI1B(x,qmu2,nf)
C     =============================================
      
      implicit double precision (a-h,o-z)
      
      qq       = qmu2    !avoid compiler warning
      dqcPMI1B = pm1sfunc(x,nf)
      
      return
      end

C     =============================================
      double precision function dqcPQQ1A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcPQQ1A = ff1sfunc(x,nf)-xf1tfunc(x,nf)
      
      return
      end
      
C     =============================================
      double precision function dqcPQQ1B(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcPQQ1B = xf1tfunc(x,nf)
      
      return
      end 
      
C     =============================================
      double precision function dqcPQG1A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
C--   Here we swap notation QG (QCDNUM) = GF (Furmanski&Petronzio)      
      dqcPQG1A = gf1sfunc(x,nf)
      
      return
      end
      
C     =============================================
      double precision function dqcPGQ1A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
C--   Here we swap notation GQ (QCDNUM) = FG (Furmanski&Petronzio)      
      dqcPGQ1A = fg1sfunc(x,nf)
      
      return
      end

C     =============================================
      double precision function dqcPGG1A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcPGG1A = gg1sfunc(x,nf)-xg1tfunc(x,nf)
      
      return
      end
      
C     =============================================
      double precision function dqcPGG1B(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcPGG1B = xg1tfunc(x,nf)
      
      return
      end

C--   ----------------------------------------------      
C--   NLO timelike (fragmetation function evolution)     
C--   ----------------------------------------------      
      
C     =============================================
      double precision function dqcTPL1A(x,qmu2,nf)
C     =============================================
      
      implicit double precision (a-h,o-z)
      
      qq       = qmu2    !avoid compiler warning
      dqcTPL1A = pp1tfunc(x,nf)-pm1tfunc(x,nf)
      
      return
      end
      
C     =============================================
      double precision function dqcTPL1B(x,qmu2,nf)
C     =============================================
      
      implicit double precision (a-h,o-z)
      
      qq       = qmu2    !avoid compiler warning
      dqcTPL1B = pm1tfunc(x,nf)
      
      return
      end
      
C     =============================================
      double precision function dqcTMI1B(x,qmu2,nf)
C     =============================================
      
      implicit double precision (a-h,o-z)
      
      qq       = qmu2    !avoid compiler warning
      dqcTMI1B = pm1tfunc(x,nf)
      
      return
      end

C     =============================================
      double precision function dqcTQQ1A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcTQQ1A = ff1tfunc(x,nf)-xf1tfunc(x,nf)
      
      return
      end
      
C     =============================================
      double precision function dqcTQQ1B(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcTQQ1B = xf1tfunc(x,nf)
      
      return
      end 
      
C     =============================================
      double precision function dqcTQG1A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
C--   Here we swap notation QG (QCDNUM) = GF (Furmanski&Petronzio)      
      dqcTQG1A = gf1tfunc(x,nf)
      
      return
      end
      
C     =============================================
      double precision function dqcTGQ1A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
C--   Here we swap notation GQ (QCDNUM) = FG (Furmanski&Petronzio)      
      dqcTGQ1A = fg1tfunc(x,nf)
      
      return
      end

C     =============================================
      double precision function dqcTGG1A(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcTGG1A = gg1tfunc(x,nf)-xg1tfunc(x,nf)
      
      return
      end
      
C     =============================================
      double precision function dqcTGG1B(x,qmu2,nf)
C     =============================================

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcTGG1B = xg1tfunc(x,nf)
      
      return
      end
      
C--   ---------------------------------------------      
C--   NNLO spacelike (no timelike available)     
C--   ---------------------------------------------       
      
C     =============================================
      double precision function dqcPPL2A(x,qmu2,nf)
C     =============================================

C--   Division by 8 because the splitting functions assume
C--   an expansion in alphas/4pi while QCDNUM uses alphas/2pi.

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcPPL2A = p2nspa(x,nf)/8.D0
      
      return
      end 
      
C     =============================================
      double precision function dqcPPL2B(x,qmu2,nf)
C     =============================================

C--   Division by 8 because the splitting functions assume
C--   an expansion in alphas/4pi while QCDNUM uses alphas/2pi.

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcPPL2B = p2nsb(x,nf)/8.D0
      
      return
      end 
      
C     =============================================
      double precision function dqcPPL2D(x,qmu2,nf)
C     =============================================

C--   Division by 8 because the splitting functions assume
C--   an expansion in alphas/4pi while QCDNUM uses alphas/2pi.

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcPPL2D = p2nspc(x,nf)/8.D0
      
      return
      end 
      
C     =============================================
      double precision function dqcPMI2A(x,qmu2,nf)
C     =============================================

C--   Division by 8 because the splitting functions assume
C--   an expansion in alphas/4pi while QCDNUM uses alphas/2pi.

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcPMI2A = p2nsma(x,nf)/8.D0
      
      return
      end 
      
C     =============================================
      double precision function dqcPMI2B(x,qmu2,nf)
C     =============================================

C--   Division by 8 because the splitting functions assume
C--   an expansion in alphas/4pi while QCDNUM uses alphas/2pi.

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcPMI2B = p2nsb(x,nf)/8.D0
      
      return
      end 
      
C     =============================================
      double precision function dqcPMI2D(x,qmu2,nf)
C     =============================================

C--   Division by 8 because the splitting functions assume
C--   an expansion in alphas/4pi while QCDNUM uses alphas/2pi.

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcPMI2D = p2nsmc(x,nf)/8.D0
      
      return
      end 
      
C     =============================================
      double precision function dqcPVA2A(x,qmu2,nf)
C     =============================================

C--   Division by 8 because the splitting functions assume
C--   an expansion in alphas/4pi while QCDNUM uses alphas/2pi.

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcPVA2A = p2nssa(x,nf)/8.D0
      
      return
      end   
      
C     =============================================
      double precision function dqcPQQ2A(x,qmu2,nf)
C     =============================================

C--   Division by 8 because the splitting functions assume
C--   an expansion in alphas/4pi while QCDNUM uses alphas/2pi.

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcPQQ2A = p2psa(x,nf)/8.D0
      
      return
      end  
      
C     =============================================
      double precision function dqcPQG2A(x,qmu2,nf)
C     =============================================

C--   Division by 8 because the splitting functions assume
C--   an expansion in alphas/4pi while QCDNUM uses alphas/2pi.

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcPQG2A = p2qga(x,nf)/8.D0
      
      return
      end  
      
C     =============================================
      double precision function dqcPGQ2A(x,qmu2,nf)
C     =============================================

C--   Division by 8 because the splitting functions assume
C--   an expansion in alphas/4pi while QCDNUM uses alphas/2pi.

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcPGQ2A = p2gqa(x,nf)/8.D0
      
      return
      end
      
C     =============================================
      double precision function dqcPGG2A(x,qmu2,nf)
C     =============================================

C--   Division by 8 because the splitting functions assume
C--   an expansion in alphas/4pi while QCDNUM uses alphas/2pi.

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcPGG2A = p2gga(x,nf)/8.D0
      
      return
      end  
      
C     =============================================
      double precision function dqcPGG2B(x,qmu2,nf)
C     =============================================

C--   Division by 8 because the splitting functions assume
C--   an expansion in alphas/4pi while QCDNUM uses alphas/2pi.

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcPGG2B = p2ggb(x,nf)/8.D0
      
      return
      end  
      
C     =============================================
      double precision function dqcPGG2D(x,qmu2,nf)
C     =============================================

C--   Division by 8 because the splitting functions assume
C--   an expansion in alphas/4pi while QCDNUM uses alphas/2pi.

      implicit double precision (a-h,o-z)

      qq       = qmu2    !avoid compiler warning
      dqcPGG2D = p2ggc(x,nf)/8.D0
      
      return
      end
      
C--   ---------------------------------------------      
C--   NNLO A coefficients for discontinuities     
C--   ---------------------------------------------       
      
C     =============================================
      double precision function dqcAGQ2A(x,qmu2,nf)
C     =============================================

C--   Division by 4 because the splitting functions assume
C--   an expansion in alphas/4pi while QCDNUM uses alphas/2pi.

      implicit double precision (a-h,o-z)

      qq        = qmu2    !avoid compiler warning      
      idum      = nf      !avoid compiler warning
      dqcAGQ2A  = a2gq(x,10.D0,10.D0,10.D0)/4.D0
      
      return
      end 
      
C     =============================================
      double precision function dqcAGG2A(x,qmu2,nf)
C     =============================================

C--   Division by 4 because the splitting functions assume
C--   an expansion in alphas/4pi while QCDNUM uses alphas/2pi.

      implicit double precision (a-h,o-z)

      qq        = qmu2    !avoid compiler warning      
      idum      = nf      !avoid compiler warning
      dqcAGG2A  = a2gg(x,10.D0,10.D0,10.D0)/4.D0
      
      return
      end  
      
C     =============================================
      double precision function dqcAGG2B(x,qmu2,nf)
C     =============================================

C--   Division by 4 because the splitting functions assume
C--   an expansion in alphas/4pi while QCDNUM uses alphas/2pi.

      implicit double precision (a-h,o-z)

      qq        = qmu2    !avoid compiler warning      
      idum      = nf      !avoid compiler warning
      dqcAGG2B  = softg(x,10.D0,10.D0,10.D0)/4.D0
      
      return
      end 
      
C     =============================================
      double precision function dqcAGG2D(x,qmu2,nf)
C     =============================================

C--   Division by 4 because the splitting functions assume
C--   an expansion in alphas/4pi while QCDNUM uses alphas/2pi.

      implicit double precision (a-h,o-z)

      qq        = qmu2    !avoid compiler warning      
      idum      = nf      !avoid compiler warning
      dqcAGG2D  = corg2(x,10.D0,10.D0,10.D0)/4.D0
      
      return
      end
      
C     =============================================
      double precision function dqcAQQ2A(x,qmu2,nf)
C     =============================================

C--   Division by 4 because the splitting functions assume
C--   an expansion in alphas/4pi while QCDNUM uses alphas/2pi.

      implicit double precision (a-h,o-z)

      qq        = qmu2    !avoid compiler warning      
      idum      = nf      !avoid compiler warning
      dqcAQQ2A  = a2qqns(x,10.D0,10.D0,10.D0)/4.D0
      
      return
      end  
      
C     =============================================
      double precision function dqcAQQ2B(x,qmu2,nf)
C     =============================================

C--   Division by 4 because the splitting functions assume
C--   an expansion in alphas/4pi while QCDNUM uses alphas/2pi.

      implicit double precision (a-h,o-z)

      qq        = qmu2    !avoid compiler warning      
      idum      = nf      !avoid compiler warning
      dqcAQQ2B  = softq2(x,10.D0,10.D0,10.D0)/4.D0
      
      return
      end 
      
C     =============================================
      double precision function dqcAQQ2D(x,qmu2,nf)
C     =============================================

C--   Division by 4 because the splitting functions assume
C--   an expansion in alphas/4pi while QCDNUM uses alphas/2pi.

      implicit double precision (a-h,o-z)

      qq        = qmu2    !avoid compiler warning      
      idum      = nf      !avoid compiler warning
      dqcAQQ2D  = corq2(x,10.D0,10.D0,10.D0)/4.D0
      
      return
      end  
      
C     =============================================
      double precision function dqcAHQ2A(x,qmu2,nf)
C     =============================================

C--   Division by 4 because the splitting functions assume
C--   an expansion in alphas/4pi while QCDNUM uses alphas/2pi.

      implicit double precision (a-h,o-z)

      qq        = qmu2    !avoid compiler warning      
      idum      = nf      !avoid compiler warning
      dqcAHQ2A  = a2qq(x,10.D0,10.D0,10.D0)/4.D0
      
      return
      end 
      
C     =============================================
      double precision function dqcAHG2A(x,qmu2,nf)
C     =============================================

C--   Division by 4 because the splitting functions assume
C--   an expansion in alphas/4pi while QCDNUM uses alphas/2pi.

      implicit double precision (a-h,o-z)

      qq        = qmu2    !avoid compiler warning      
      idum      = nf      !avoid compiler warning
      dqcAHG2A  = a2hga(x,10.D0,10.D0,10.D0)/4.D0
      
      return
      end 
      
C     =============================================
      double precision function dqcAHG2D(x,qmu2,nf)
C     =============================================

C--   Division by 4 because the splitting functions assume
C--   an expansion in alphas/4pi while QCDNUM uses alphas/2pi.

      implicit double precision (a-h,o-z)

      qq        = qmu2    !avoid compiler warning      
      idum      = nf      !avoid compiler warning
      dqcAHG2D  = a2hgc(x,10.D0,10.D0,10.D0)/4.D0
      
      return
      end 
      
      