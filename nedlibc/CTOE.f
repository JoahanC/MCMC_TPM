	SUBROUTINE CTOE(C,E)
	REAL*4 C(3),E(3),TEMP
C
C	C	INPUT	vector in celestial coordinates
C	E	OUTPUT	vector in ecliptic coordinates
C
	PARAMETER (COBL=.917482058,SOBL=.39777716)
	E(1)=C(1)
	TEMP=COBL*C(3)-SOBL*C(2)
	E(2)=COBL*C(2)+SOBL*C(3)
	E(3)=TEMP
	RETURN
	END
