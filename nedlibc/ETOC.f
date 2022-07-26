	SUBROUTINE ETOC(E,C)
	REAL*4 E(3),C(3),TEMP
C
C	E	INPUT	vector in ecliptic coordinates (J2000.0)
C	C	OUTPUT	vector in celestial coordinates (J2000.0)
C
	PARAMETER (COBL=.917482058,SOBL=.39777716)
	C(1)=E(1)
	TEMP=COBL*E(3)+SOBL*E(2)
	C(2)=COBL*E(2)-SOBL*E(3)
	C(3)=TEMP
	RETURN
	END
