	FUNCTION RAN2(JDUM)
C
C	from the Numerical Recipes column in Computers in Physics
C	Sep/Oct 1992, vol 6, p 522.  W. H. Press & S. A. Teukolsky
C
C	if argument is negative, sequence is reset based on the value
C	and the argument is changed to a positive value
C	otherwise the argument JDUM has no effect and is not affected
C
	IMPLICIT NONE
	INTEGER JDUM,IDUM,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,
	1	NTAB,NDIV
	REAL RAN2,AM,EPS,RNMX
	PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.0/IM1,IMM1=IM1-1,
	1	IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
	2	IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2E-7,
	3	RNMX=1.0-EPS)
	INTEGER IDUM2,J,K,IV(NTAB),IY
	SAVE IV,IY,IDUM,IDUM2
	DATA IDUM2 /123456789/, IV/NTAB*0/, IY/0/
	IF (JDUM.LE.0) THEN
	  IDUM = MAX(-JDUM,1)
	  IDUM2 = IDUM
	  DO J=NTAB+8,1,-1
	    K = IDUM/IQ1
	    IDUM = IA1*(IDUM-K*IQ1)-K*IR1
	    IF (IDUM.LT.0) IDUM = IDUM+IM1
	    IF (J.LE.NTAB) IV(J) = IDUM
	  ENDDO
	  JDUM = IDUM
	  IY = IV(1)
	ENDIF
	K = IDUM/IQ1
	IDUM = IA1*(IDUM-K*IQ1)-K*IR1
	IF (IDUM.LT.0) IDUM = IDUM+IM1
	K = IDUM2/IQ2
	IDUM2 = IA2*(IDUM2-K*IQ2)-K*IR2
	IF (IDUM2.LT.0) IDUM2 = IDUM2+IM2
	J = 1+IY/NDIV
	IY = IV(J)-IDUM2
	IV(J) = IDUM
	IF (IY.LT.1) IY = IY+IMM1
	RAN2 = MIN(AM*IY,RNMX)
	RETURN
	END
