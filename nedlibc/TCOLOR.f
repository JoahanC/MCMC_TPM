	REAL FUNCTION TCOLOR(NU1,NU2,RATIO)
	REAL NU1,NU2,RATIO
C
C	NU1 is frequency 1 in 1/cm
C	NU2 is frequency 2 in 1/cm
C	RATIO is B_nu2(Tc)/B_nu1(Tc)
C
	REAL HCK
	PARAMETER (HCK=1.4387752)
	REAL N1,N2
C	test for bad input
	IF (NU1.LE.0..OR.NU2.LE.0..OR.NU1.EQ.NU2.OR.RATIO.LE.0.) THEN
	  WRITE (*,*) 'Bad input to TCOLOR:',NU1,NU2,RATIO
	  STOP
	ENDIF
C	copy input to local variables
	N1 = NU1
	N2 = NU2
	RAT = RATIO
C	put in order to NU1 < NU2
	IF (N1.GT.N2) THEN
	  N2 = NU1
	  N1 = NU2
	  RAT = 1/RATIO
	ENDIF
C	use B = 1/T as the variable to solve for
	BL = 0
	YL = LOG((N2/N1)**2/RAT)
C	test for steeper than Rayleigh-Jeans, return "infinity" if so
	TCOLOR = 1.E10
	IF (YL.LE.0.) RETURN
C	pick the low T (high B) end to avoid overflowing the EXP function
	BH = MIN(30/((N2-N1)*HCK),140/((N1+N2)*HCK))
	YH = LOG((N2/N1)**3*(EXP(HCK*N1*BH)-1)/(EXP(HCK*N2*BH)-1)/RAT)
C
	NIT = 0
100	SLOPE = (YH-YL)/(BH-BL)
	DB = YL/SLOPE
	BT = BL-DB
	YT = LOG((N2/N1)**3*(EXP(HCK*N1*BT)-1)/(EXP(HCK*N2*BT)-1)/RAT)
	IF (ABS(YT).GT.1.E-4) THEN
	  IF (ABS(YH).GT.ABS(YL)) THEN
	    BH = BT
	    YH = YT
	  ELSE
	    BL = BT
	    YL = YT
	  ENDIF
	  NIT = NIT+1
	  IF (NIT.LT.100) GO TO 100
	ENDIF
	TCOLOR = 1/BT
	RETURN
	END
