	SUBROUTINE SPRBIT_SET
C
C	sets up bit table
C
	INTEGER*4 IX,IY,K,IP
	SAVE /SPRBIT/
	COMMON /SPRBIT/ IX(128),IY(128)
	DATA IX(128) /0/
	WRITE (*,*) 'SPRBIT_SET called'
	IF (IX(128).EQ.0) THEN			!Set up the lookup tables
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
	RETURN
	END
