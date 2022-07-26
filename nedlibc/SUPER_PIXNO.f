	SUBROUTINE SUPER_PIXNO(C,NSUPER)
	REAL*4 C(3)	!Unit vector input
	INTEGER*4 NSUPER	!Output pixel number 0,...,6*2**28-1
	PARAMETER (IT14=16384,IT28=2**28)
C
C	PIXELIZES SPHERE TO 20" ACCURACY, EQUAL-AREA PIXELS USING
C	QUADRILATERIZED SPHERICAL CUBE
C
	INTEGER*4 IX,IY,K,IP
	SAVE /COBE_SPRBIT/
	COMMON /COBE_SPRBIT/ IX(128),IY(128)
	DATA IX(128) /0/
	IF(IX(128).EQ.0) THEN			!Set up the lookup tables
	   DO I=1,128				!for converting x,y into
	      J=I-1				!pixel numbers
	      K=0
	      IP=1
10	      IF(J.EQ.0) THEN
		IX(I)=K
		IY(I)=2*K
	      ELSE
		ID=MOD(J,2)
		J=J/2
		K=IP*ID+K
		IP=IP*4
		GO TO 10
	      ENDIF
	   ENDDO
	ENDIF
	CALL AXISXY(C,NFACE,X,Y)
	I=IT14*X
	J=IT14*Y
	I=MIN(IT14-1,I)			!Roundoff protection for cases very
	J=MIN(IT14-1,J)			!near an edge of the cube
	IL=MOD(I,128)
	IH=I/128
	JL=MOD(J,128)
	JH=J/128
	NSUPER=NFACE*IT28+IX(IL+1)+IY(JL+1)+IT14*(IX(IH+1)+IY(JH+1))
	RETURN
	END
