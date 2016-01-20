!---------------------------------------------------------------------------  
! DESCRIPTION: 
!> Elastic dipole amplitude (in the forward limit total dipole cross-section)
!> @brief
!> Elastic dipole amplitude (in the forward limit total dipole cross-section)
!> 
!> @param[in] xBj Bjorken  <i>x</i>
!> @param[in] r transverse dipole size
!> @param[in] b impact parameter
!> @param[in] param which parametrization to choose (1: m_c=1.27, 2: m_c=1.4 parameter set; default is the first one)
!> @return dipole amplitude
!--------------------------------------------------------------------------- 

	DOUBLE PRECISION FUNCTION dipole_amplitude(xBj, r, b, param)
	IMPLICIT NONE
	INTEGER NPAR, PARAM, SCHEME
	PARAMETER (NPAR=3)
	DOUBLE PRECISION PARS(NPAR)
	DOUBLE PRECISION xBj, r, b
	DOUBLE PRECISION Pi, mu0Squared
	DOUBLE PRECISION B_eff, B_factor
	DOUBLE PRECISION xgValue, mu2, alpha_S_VALUE, FVALXQ
        DOUBLE PRECISION m_N, m_q, m_c
	Data Pi/3.1415926/
	DATA B_eff/4/
	SAVE mu0Squared
        COMMON /MASSES/ m_N, m_q, m_c
	COMMON /CURRENT_SCHEME/SCHEME

	IF(SCHEME.ne.param)THEN
	 IF(PARAM.ne.1.and.PARAM.ne.2)PARAM=1
	 SCHEME=PARAM
         IF(SCHEME.eq.1)THEN
	 PARS(1)=2.308059
	 PARS(2)=-0.057584
	 PARS(3)=1.513333
	 m_c=1.27
	 ELSE
	 PARS(1)=2.373
	 PARS(2)=-0.052
	 PARS(3)=1.428
	 m_c=1.4
	 ENDIF
	call init_gluons(npar, pars) ! This is sloooow
	mu0Squared=PARS(3)
	ENDIF


	B_factor=exp(-b**2/(2*B_eff))/(2*Pi*B_eff)


	if(r.le.0)THEN ! This is unphysical/nonperturbative regime, model is not valid there
	 dipole_amplitude=0
	 RETURN
	ENDIF

	mu2=mu0Squared+4./r**2
	if(mu2.gt.1E7)THEN ! Too large value of scale should be consistent with upper scale in gluonEvolve; Q2<mu2 make no sense
	 mu2=1E7
	ENDIF
	xgValue=FVALXQ( 1, 0, xBj, mu2, 0 )
	IF(xgValue.gt.0.99E11)THEN !! We are outside evaluated region for PDFs
	 dipole_amplitude=0
	 RETURN
	ENDIF

!       LO Nf=4 scheme,  Lambda_QCD = 0.156447
	alpha_S_VALUE=1.50796/Log(mu2/0.0244757)
	dipole_amplitude=2*(1-exp(-Pi**2*r**2/6*alpha_S_VALUE
     *         *xgValue*B_factor))
	RETURN
	END

