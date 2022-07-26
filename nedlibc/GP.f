	FUNCTION GP(XX)
C
C	for XX in units of standard deviations, returns the probability
C	that a Gaussian is less than XX
C
	PARAMETER (GLN = 0.57236494292159)
        PARAMETER (ITMAX=150,EPS=3.E-7)
	PARAMETER (A = 0.5)
C
	IF (XX.LT.-5.6) THEN
	  GP = 0.
	  RETURN
	ENDIF
	IF (XX.GT.5.6) THEN
	  GP = 1.0
	  RETURN
	ENDIF
	IF (XX.GT.0.) THEN
	  IS = 1
	ELSE
	  IS = -1
	ENDIF
	X = XX*XX/2
	IF (X.LT.1.5) THEN
	  AP=A
	  SUM=1./A
	  DEL=SUM
	  DO N=1,ITMAX
	    AP=AP+1.
	    DEL=DEL*X/AP
	    SUM=SUM+DEL
	    IF (ABS(DEL).LT.ABS(SUM)*EPS) GO TO 1
	  ENDDO
          STOP 'gser'
1	  V = SUM*EXP(-X+A*LOG(X)-GLN)
	ELSE
	  GOLD=0.
	  A0=1.
	  A1=X
	  B0=0.
	  B1=1.
	  FAC=1.
	  DO N=1,ITMAX
	    AN=FLOAT(N)
	    ANA=AN-A
	    A0=(A1+A0*ANA)*FAC
	    B0=(B1+B0*ANA)*FAC
	    ANF=AN*FAC
	    A1=X*A0+ANF*A1
	    B1=X*B0+ANF*B1
	    IF (A1.NE.0.) THEN
	      FAC=1./A1
	      G=B1*FAC
	      IF (ABS((G-GOLD)/G).LT.EPS) GO TO 2
	      GOLD=G
	    ENDIF
	  ENDDO
	  STOP 'gcf'
2	  V = EXP(-X+A*ALOG(X)-GLN)*G
	  V = 1-V
	ENDIF
	GP = 0.5+IS*V/2
	RETURN
	END
