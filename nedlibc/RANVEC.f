	SUBROUTINE RANVEC(U)
C	make a random unit vector
	REAL*4 U(3)
	EXTERNAL RAN
	REAL RAN
	INTEGER I, IR,JR
	REAL RR
C
	RR = 4
	DO WHILE (RR.GT.1.0 .OR. RR.EQ.0.0)
	  RR = 0
	  DO I=1,3
	    U(I) = 2*RAN(IR,JR)-1
	    RR = RR+U(I)**2
	  ENDDO
	ENDDO
	RR = SQRT(RR)
	DO I=1,3
	  U(I) = U(I)/RR
	ENDDO
	RETURN
	END
	
