! This file illustrates how to use dipole cross-sections
      program example
      implicit none
      DOUBLE PRECISION xBJ, r, b, res, fm
      DOUBLE PRECISION dipole_amplitude
      INTEGER I, parameterSet

	fm=1/0.19733
	r=1.*fm
	b=1.*fm
	parameterSet=1 ! 1=parameter set with m_c=1.27; 2: with m_c=1.4

	WRITE(*,2)
	WRITE(*,*)"       x                d^2 sigma/db^2    "
	WRITE(*,2)
	DO I=-50, -20
	xBj=10.**(I/10.)
	res=dipole_amplitude(xBj, r, b, parameterSet)
	WRITE(*,*)xBj,res
	ENDDO

	xBj=1e-4
	WRITE(*,2)
	WRITE(*,*)"       r, fm            d^2 sigma/db^2    "
	WRITE(*,2)
	DO I=0,30
	r=I/5.
	res=dipole_amplitude(xBj, r*fm, b, parameterSet)
	WRITE(*,*)r,res
	ENDDO

	r=1*fm
	WRITE(*,2)
	WRITE(*,*)"       b, fm            d^2 sigma/db^2    "
	WRITE(*,2)
	xBj=1e-4
	DO I=0,25
	b=I/20.
	res=dipole_amplitude(xBj, r, b*fm, parameterSet)
	WRITE(*,*)b,res
	ENDDO


1	FORMAT(2G8.4)
2	FORMAT("_______________________________________________")
      end
