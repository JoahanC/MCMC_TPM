	FUNCTION GAUSRN(IR,JR)
	EXTERNAL RAN
C       ISW not equal to 0 if we have a precomputed value
	DATA ISW /0/
	SAVE ISW,Y,F
	IF (ISW.EQ.0) GO TO 10
	GAUSRN = Y*F
	ISW = 0
	RETURN
C       compute two Gaussian random variables
C       first find X,Y uniform in the unit circle
10	X = 2.*RAN(IR,JR)-1.
	Y = 2.*RAN(IR,JR)-1.
	RR = X*X+Y*Y
	IF (RR.GE.1.) GO TO 10
C       prob of r < Rg for 2-D Gaussian is 1-exp(-0.5*Rg^2)
C       prob of r < Ru for uniform on unit circle is Ru^2
C       set probabilities equal to transform from unit circle to 2-D gaussian
C       giving Rg = SQRT(-2ln(1-Ru^2)), so F = Rg/Ru as given below
	IF (RR.LT.0.01) THEN
	  F = SQRT(2.*(1+RR*(0.5+RR/3.)))
	ELSE
	  F = SQRT(-2.*ALOG(1.-RR)/RR)
	ENDIF
C       first Gaussian random variable is F*X, with F*Y saved for next call
	GAUSRN = F*X
	ISW = 1
	RETURN
	END
