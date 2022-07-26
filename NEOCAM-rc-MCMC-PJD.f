C	NEOCAM-rc-MCMC-PJD.f
C	gfortran -o NEOCAM-rc-MCMC-PJD NEOCAM-rc-MCMC-PJD.f ~/Ned.a
C
C	modified 5-JUL-2016 to use 6 parameters: 
C	RA,Dec[rad],ln(pV),ln(period[hr]),ln(J in J/[m^2*K*rt(sec)]),ln(D[km])
C
C	modified 12-AUG-2014 to use fluxes instead of magnitudes in the
C	chi^2 calculation.  This is done by robust averaging the "area"
C	which is D^2 with D in km.
C
C	input format remains the same, but use a negative sigma(mag) for
C	band with a negative central value for the flux
C	so input m = mzp-2.5*log_10(|F|) and 
C	sigma_m = sgn(F)*2.5*log_10(1+sigma(F)/|F|)
C
	PARAMETER (NP=6,NDMX=9,NB=2)
	REAL*4 PHT(NDMX),DET(NDMX),RAT(NDMX),DECT(NDMX),DST(NDMX)
	REAL*4 SIGT(NB,NDMX),MAGT(NB,NDMX)
	REAL*8 T0(NDMX)
	COMMON /NRMLPTS/ T0,H,DH,PER,MAGT,SIGT,DST,DET,PHT,RAT,DECT,NDM
	COMMON /PREDICTED/ WMAGT(NB,NDMX)
	REAL*8 MJD,JD(NDMX)
	REAL*8 T1JAN2010,T,JD1JAN2010
	REAL*8 TSTART,TWALL
	REAL*8 JULIAN_DAY,GMT2T68
	REAL*4 EA(3),A(3),E(3)
	REAL*4 WMAG(NB),SIG(NB)
	REAL*4 P(NP)
	REAL*4 P1(3),P2(3),P3(3)
	CHARACTER*11 GMT
	CHARACTER*10 NAME
C
	COMMON /WATERMARK/ BESTEVER,PB(NP)
	PARAMETER (NW=100)
	REAL*4 PW(NP,NW),CH2(NW)
	INTEGER NTRY(NW)
	INTEGER PPMIL
	COMMON /WALKERS/ PW,CH2,NTRY
C
	REAL*4 SIGMIN(NB) /NB*0.02/
	EXTERNAL RAN
C
	CALL RANDOMIZE
	BESTEVER = 1.E30
	NDM=0
	WRITE (*,FMT='(A,$)') 'Enter H,sig(H),Period,Name:'
	READ (*,FMT='(3F7.3,A)') H,DH,PER,NAME
	WRITE (*,FMT='(3F7.3,1X,A)') H,DH,PER,NAME
	CALL TIME_TAG(2010,1,1,0,0,0.,T1JAN2010)
	JD1JAN2010 = JULIAN_DAY(T1JAN2010)
	DTR = 45/ATAN(1.)
	WRITE (*,FMT='(A)') 'Enter MJD,RA,Dec,Delta,NC12 with sigmas:'
	WRITE (*,FMT='(A)') '   MJD      RA      DEC     DE      DS'//
	1	'    PHASE NC12 & sigmas'
10	READ (*,900,END=100) MJD,RA,DEC,DELTA,(WMAG(I),SIG(I),I=1,NB)
900	FORMAT(F11.5,2F11.7,F8.4,4(F7.3,F6.3))
	NDM=NDM+1
	DO I=1,NB
	  MAGT(I,NDM) = WMAG(I)
	  SIGT(I,NDM) = SQRT(SIG(I)**2+SIGMIN(I)**2)
	  IF (SIG(I).LT.0.) SIGT(I,NDM) = -SIGT(I,NDM)
	ENDDO
	RAT(NDM) = RA
	DECT(NDM) = DEC
	DET(NDM) = DELTA
	T0(NDM) = 86400*(MJD+2400000.5D0-JD1JAN2010)+T1JAN2010
	CALL EPHEM(-1,T0(NDM),E)
	CALL RADEC_TO_C(RAT(NDM)/DTR,DECT(NDM)/DTR,EA)
	DO J=1,3
	  E(J)=-E(J)
	  EA(J) = EA(J)*DET(NDM)
	  A(J) = EA(J)+E(J)
	ENDDO
	DST(NDM) = SQRT(DOT(A,A))
	PHT(NDM) = DTR*ACOS(DOT(EA,A)/(DST(NDM)*DET(NDM)))
	WRITE (*,901) MJD,RAT(NDM),DECT(NDM),DET(NDM),DST(NDM),
	1	PHT(NDM),(WMAG(I),SIG(I),I=1,NB)
901	FORMAT(F9.3,2F8.3,2F8.5,F7.2,4(F7.3,F6.3))
	GO TO 10
100	CLOSE (UNIT=1)
	WRITE (*,*) NDM
C	Monte Carlo Markov chain
C
	WRITE (21,*) ' Monte Carlo Markov Chain'
	CLOSE (UNIT=21)
	TSTART = GMT2T68('NOW')
	IMC = 0
	TWALL = GMT2T68('NOW')-TSTART
	WRITE (*,FMT='(I4,F10.0,2H P)') IMC,TWALL
	WRITE (*,*)
C	setup the initial simplex of walkers
	DO I=1,NW
	  RADIUS = 4
	  DO WHILE (RADIUS.GT.1.5)
	    DO J=1,3
	      PW(J,I) = 3*RAN(IR,JR)-1.5
	    ENDDO
	    RADIUS = SQRT(DOT(PW(1,I),PW(1,I)))
	  ENDDO
C	to use RA,DEC,ln(p_V),ln(Period),ln(J),ln(D) as parameters
	  CALL C_TO_RADEC(PW(1,I),RA1,DEC1)
	  PW(1,I) = RA1
	  PW(2,I) = DEC1
C	albedo uniform in log 0.01 to 1.01
150	  PW(3,I) = ALOG(RAN(IR,JR)+0.01)
C	periods 0.0001 hours to 1000 hours
	  IF (PER.LE.0.) THEN
	    PW(4,I) = -4*ALOG(10.)+7*ALOG(10.)*RAN(IR,JR)
	  ELSE
	    PW(4,I) = ALOG(PER)
	  ENDIF
C	J = sqrt(C\rho\kappa) from 5 to 2500
	  PW(5,I) = ALOG(5.)+RAN(IR,IR)*ALOG(500.)
C	diameter from 1 meter to 1000 km
	  PW(6,I) = ALOG(.001)+ALOG(1.E6)*RAN(IR,JR)
	  DO J=1,6
	    P(J) = PW(J,I)
	  ENDDO
	  EP = PJD_PV_PENALTY(P)
	  U = EXP(-0.5*EP)
C	  IF (RAN(IR,JR).GT.U) GO TO 150
	  D1 = EXP(P(6))
	  CH2(I) = ERROR(P,1.E30)
	  E1 = CH2(I)
C	reset diameter for better start
	  SWDD = 0
	  SW = 0
	  IF (DH.LT.8.) THEN
	    SW = SW+1
C	SDD is the ln[(D from H & PV)/D]
	    SWDD=SWDD+ALOG(1329./10.**(H/5)/SQRT(EXP(PW(3,I))))-PW(6,I)
	  ENDIF
	  DO K=1,NDM
	    DO J=1,NB
	      IF (SIGT(J,K).LT.8.) THEN
	        SW=SW+1
C	add 0.4605*(pred-obs) mag so if D is too big, get a negative
	        SWDD = SWDD+0.4605*(WMAGT(J,K)-MAGT(J,K))
	      ENDIF
	    ENDDO
	  ENDDO
	  PW(6,I) = (SWDD/SW)+PW(6,I)
	  P(6) = (SWDD/SW)+P(6)
	  D2 = EXP(P(6))
	  CH2(I) = ERROR(P,1.E30)
	  E2 = CH2(I)
C	write (*,*) d1,e1,d2,e2
C
	  NTRY(I) = 0
	  WRITE (*,FMT='(I3,6F8.3,F15.2)') I,(PW(J,I),J=1,6),CH2(I)
	ENDDO
	WRITE (*,*)
	DO IMC=1,10000
	  DO KW=1,NW
	    JW = KW
	    DO WHILE (JW.EQ.KW)
	      XW = 1.+NW*RAN(IR,JR)
	      JW = XW
	      JW = MIN(NW,MAX(1,JW))
	    ENDDO
	    U = RAN(IR,JR)
C	dz/sqrt(z) \propto d[sqrt(z)]
	    AA = 2
	    Z = SQRT(1/AA)+(SQRT(AA)-SQRT(1/AA))*U
C	median = 1.5*sqrt(0.5) = 1.06066 above, 1.125 below
	    Z = Z*Z
C	lowest value of Z is 1/AA, highest is AA
C	1st %-tile is 0.51005, 99th %-tile is 1.98005
C	z=0, stay at JWth walker, Z=1 move to KWth walker
	    DO I=3,6
	      P(I) = PW(I,JW)+Z*(PW(I,KW)-PW(I,JW))
	    ENDDO
C	stretching in angle code
	    CALL RADEC_TO_C(PW(1,KW),PW(2,KW),P2)
	    CALL RADEC_TO_C(PW(1,JW),PW(2,JW),P1)
	    TH = 2*ASIN(SQRT(DVS(P2,P1)/4))
C	    WRITE (*,FMT='(6F9.5,F11.5)') P1,P2,TH
	    CALL ORTH2(P1,P2,P3)
	    ZTH = Z*TH
	    CZTH = COS(ZTH)
	    SZTH = SIN(ZTH)
C	z=0, stay at P1 which is the JWth walker, z=1 move to KWth walker
	    DO I=1,3
	      P3(I) = CZTH*P1(I)+SZTH*P2(I)
	    ENDDO
	    TH = 2*ASIN(SQRT(DVS(P3,P1)/4))
C	    IF (ABS(TH-ZTH).GT.0.001) WRITE (*,*) '-----> ',ZTH,TH
C	    WRITE (*,FMT='(F8.5,3F9.5)') Z,P3
	    CALL C_TO_RADEC(P3,P(1),P(2))
	    TWALL = GMT2T68('NOW')-TSTART
	    IF (KW.NE.1) OPEN (UNIT=23,ACCESS='APPEND')
	    WRITE (23,FMT='(I5,2I4,7F8.4,F9.0)') IMC,KW,JW,Z,P,TWALL
	    CHSQ = ERROR(P,CH2(KW)+100)
	    Q = 0
	    XP = -0.5*(CHSQ-CH2(KW))
	    IF (XP.GT.50.) Q = 1
	    IF (ABS(XP).LT.50.) 
	1	Q = MIN(1.0,ABS(SZTH/SIN(TH))*Z**4*EXP(XP))
	    U = RAN(IR,JR)
	    WRITE (23,FMT='(4F13.5)') CHSQ,CH2(KW),Q,U
	    CLOSE (UNIT=23)
	    IF (U.LT.Q) THEN
C	accept
	      OPEN (UNIT=21,ACCESS='APPEND')
C	write out the old position
	      WRITE (21,FMT='(F8.5,5F9.5,F17.5,I5)') 
	1	(PW(I,KW),I=1,NP),CH2(KW),NTRY(KW)
	      CLOSE (UNIT=21)
	      WRITE (*,FMT='(F7.4,5F8.4,F17.5,I4,I3)') 
	1	(PW(I,KW),I=1,NP),CH2(KW),NTRY(KW),(KW-1)
C	move to the new position
	      DO I=1,NP
	        PW(I,KW) = P(I)
	      ENDDO
	      CH2(KW) = CHSQ
	      NTRY(KW) = 1 
	    ELSE
C	reject: just increment the weight
	      NTRY(KW) = NTRY(KW)+1
	      WRITE (*,FMT='(F7.4,5F8.4,2F11.3,I3)') 
	1	P,CHSQ,CH2(KW),(KW-1)
	    ENDIF
	  ENDDO
	  TWALL = GMT2T68('NOW')-TSTART
	  WRITE (*,FMT='(I5,F10.0,2H P)') IMC,TWALL
	ENDDO
	STOP
	END

	SUBROUTINE NEOCAM_ORBIT(T68,E)
C	fake it up using 0.99 of Earth position
	REAL*8 T68
	REAL*4 E(3)
	CALL EPHEM(-1,T68,E)
	DO I=1,3
	  E(I) = -0.99*E(I)
	ENDDO
	RETURN
	END

	REAL*4 FUNCTION ERROR(P,QUIT)
C
	PARAMETER (NP=6)
	REAL*4 P(NP),QUIT
C       P(1..2) = RA,Dec in radians of pole
C	P(3) = ln(p_V)
C	P(4) = ln(Period) in hours
C	P(5) = ln(J) in MKS
C	P(6) = ln(dia) in km
C	if ERROR>QUIT then return without adding more terms
C
	PARAMETER (NDMX=9)
	PARAMETER (NBAND=2,NL=6)
	REAL*4 PHT(NDMX),DET(NDMX),RAT(NDMX),DECT(NDMX),DST(NDMX)
	REAL*4 SIGT(NBAND,NDMX),MAGT(NBAND,NDMX)
	REAL*8 T0(NDMX)
	COMMON /NRMLPTS/ T0,H,DH,PER,MAGT,SIGT,DST,DET,PHT,RAT,DECT,NDM
	COMMON /PREDICTED/ WMAGT(NBAND,NDMX)
	COMMON /BACKCHAN/ DIA,SDIA
	COMMON /WATERMARK/ BESTEVER,PB(NP)
C
	REAL*4 R(3),N(3),POLE(3)
	PARAMETER (DTR = 57.2957795, PI = 3.14159265)
	REAL*4 WMAG(NBAND)
	REAL*4 LAM(NL) /4.4130,5.0450,6.1758,7.2171,8.0235,9.6973/
	REAL*4 W(NL) /0.6778,0.3165,0.19955,0.18515,0.37066,0.20223/
        REAL*4 WVL0(NBAND) /4.6028,7.7542/
	PARAMETER (NLT=80)
	REAL*4 WVLT(NLT)
	PARAMETER (NALL=NBAND+NL+NLT)
	REAL*4 LN10OVER5
	PARAMETER (LN10OVER5 = 0.460517)
	REAL*4 WVALL(NALL)
	EQUIVALENCE (WVALL,LAM),(WVALL(NL+1),WVL0),
	1	(WVALL(NL+NBAND+1),WVLT)
	REAL*4 F(NALL)
	REAL*4 FALL(NALL,NDMX)
	REAL*4 F0(NBAND) /170.663,63.182/
	REAL*4 SIGMIN(NBAND) /NBAND*0.10/
	REAL*4 MICRONS,NUFNU
	DO I=1,NLT
	  WVLT(I) = 0.55*(28./0.55)**((I-1)/(NLT-1.))
	ENDDO
	CALL RADEC_TO_C(P(1),P(2),POLE)
C
	OMEGA = 2*PI/(3600*EXP(P(4)))
	INERTIA = EXP(P(5))
	EM = 0.95
C	T1 = 401 K
	T1 = (1367./(EM*5.67E-8))**0.25
	THETA1 = INERTIA*SQRT(OMEGA)/(EM*5.67E-8*T1**3)
C	range restrictions to protect the RC code
	TH1MAX = 10000
	THETA1 = TH1MAX*THETA1/(TH1MAX+THETA1)
	PVMAX = 1.5
	IF (P(3).GT.23.) THEN
	  PV = PVMAX
	ELSE
	  PV = EXP(P(3))
	  PV = PVMAX*PV/(PVMAX+PV)
	ENDIF
C
C	compute priors first
	ERROR = PJD_PV_PENALTY(P)
	IF (ERROR.GT.QUIT) RETURN
C
C
	ALB = 0.38399*PV
C
	ASYM = 1.
	ROT = 0.
C	write (*,fmt='(3f8.4,f10.4,1pe10.3)') pole,alb,theta1
	DO J=1,NDM
	  CALL EPHEM(-1,T0(J),R)
	  CALL RADEC_TO_C(RAT(J)/DTR,DECT(J)/DTR,N)
	  THETA = THETA1*DST(J)**1.5
	  DO I=1,3
	    R(I) = -R(I)
	  ENDDO
	  CALL PREDICT(R,N,DET(J),POLE,THETA,ALB,EM,ASYM,ROT,
	1	NALL,WVALL,F)
	  DO I=1,NALL
	    FALL(I,J) = F(I)
	  ENDDO
	  WMAGT(1,J)=-2.5*LOG10((W(1)*F(1)+W(2)*F(2))/F0(1))
	  WMAGT(2,J)=-2.5*LOG10((W(3)*F(3)+W(4)*F(4)+W(5)*F(5)+
	1	W(6)*F(6))/F0(2))
	ENDDO
C
	DIA = EXP(P(6))
	DM = -5*P(6)/ALOG(10.)
C
	DO J=1,NDM
	  DO I=1,NALL
	    FALL(I,J) = FALL(I,J)*DIA**2
	  ENDDO
	  DO I=1,NBAND
	    WMAGT(I,J) = WMAGT(I,J)+DM
	  ENDDO
	ENDDO
C
C	compute magnitude errors
	RCS = 0
	IF (DH.LT.2.) THEN
	  HCLC = 5*ALOG10(1329./(DIA*SQRT(PV)))
	  E = ABS(H-HCLC)/DH
	  EE = E*E
	  IF (E.GT.2.) EE=4+2*(E-2)
	  RCS = EE
	ENDIF
	DO J=1,NDM
	  DO I=1,NBAND
	    IF (SIGT(I,J).LT.8.99) THEN
C	input m = mzp-2.5*log_10(|F|) and 
C	sigma_m = sgn(F)*2.5*log_10(1+sigma(F)/|F|)
	      FOBS = 10.**(-0.4*MAGT(I,J))
	      IF (SIGT(I,J).LT.0.) FOBS = -FOBS
	      FCLC = 10.**(-0.4*WMAGT(I,J))
	      SIGF = ABS(FOBS)*(10.**(0.4*ABS(SIGT(I,J)))-1)
	      E = ABS(FOBS-FCLC)/SIGF
	      EE = E*E
	      IF (E.GT.2.) EE=4+2*(E-2)
	      RCS = RCS+EE
	    ENDIF
	  ENDDO
	ENDDO
C
C	add chi^2 from magnitude errors
	ERROR = ERROR+RCS
C
C	watermark code
C
	IF (ERROR.LT.BESTEVER) THEN
	  BESTEVER = ERROR
	  DO I=1,NP
	    PB(I) = P(I)
	  ENDDO
	  WRITE (22,FMT='(1H%,F11.6,5F12.6,2F12.4)') P,ERROR,RCS
	  WRITE (22,FMT='(1H%,F10.5,1P2E11.4)') PV,DIA,THETA1
	  WRITE (22,FMT='(A)') '/wvl ['
	  WRITE (22,FMT='(11F7.3)') (WVALL(I+NL+NBAND),I=1,NLT)
	  WRITE (22,FMT='(A)') '] def'
	  DO J=1,NDM
C	set up colors for different epochs
	    IF (MOD(J,4).EQ.1.AND.NDM.GT.1) 
	1	WRITE (22,FMT='(A)') '1 0 0 SRGB'
	    IF (MOD(J,4).EQ.2) WRITE (22,FMT='(A)') '0 0.7 0 SRGB'
	    IF (MOD(J,4).EQ.3) WRITE (22,FMT='(A)') '0 0 1 SRGB'
	    IF (MOD(J,4).EQ.0) WRITE (22,FMT='(A)') '0 setgray'
C
	    DO I=1,NBAND
	      DM = MAGT(I,J)-WMAGT(I,J)
	      FJY = FALL(I+NL,J)/10.**(0.4*DM)
	      MICRONS = WVALL(I+NL)
	      NUFNU = 1.E-23*FJY*(3.E14/MICRONS)
	      IF (SIGT(I,J).LT.9.00)
	1       WRITE (22,FMT='(2F7.3,1PE10.3,3H QQ)') 
	1		MICRONS,SIGT(I,J),NUFNU
	    ENDDO
	    WRITE (22,FMT='(A)') '/ft ['
	    DO I=1,NLT
	      MICRONS = WVALL(I+NL+NBAND)
	      NUFNU = 1.E-23*FALL(I+NL+NBAND,J)*(3.E14/MICRONS)
	      WRITE (22,FMT='(1PE10.3)') NUFNU
	    ENDDO
	    WRITE (22,FMT='(A)') '] def doit'
	  ENDDO
	  CLOSE (UNIT=22)
	ENDIF
C
	RETURN
	END

	FUNCTION PJD_PV_PENALTY(P)
	REAL*4 INERTIA,JLO,JHI
	REAL*4 ERROR
	PARAMETER (NP=6)
	REAL*4 P(NP),QUIT
C       P(1..2) = RA,Dec in radians of pole
C	P(3) = ln(p_V)
C	P(4) = ln(Period) in hours
C	P(5) = ln(J) in MKS
C	P(6) = ln(dia) in km
C	parameters of albedo prior
	PARAMETER (PK1=0.03,FD=0.251,PK2=0.167)
C	the PRE table shows the effect of the P, J, D prior
C	divide by PRE to get uniform in log(D) diameter prior
	INTEGER PRE(140)/6365,6371,6378,6384,6391,6399,6407,6415,6423,
	1  6432,6441,6451,6461,6472,6484,6496,6509,6522,6537,6552,6569,
	2  6586,6605,6625,6647,6670,6695,6722,6751,6783,6817,6855,6896,
	3  6941,6990,7044,7104,7170,7243,7324,7414,7514,7625,7747,7882,
	4  8030,8191,8365,8548,8740,8936,9133,9326,9512,7208,7267,7326,
	5  7383,7440,7497,7553,7608,7662,7715,7768,7820,7872,7922,7972,
	6  8021,8069,8117,8164,8210,8256,8300,8344,8387,8430,8472,8513,
	7  8553,8593,8632,8670,8708,8745,8782,8817,8853,8887,8921,8954,
	8  8987,9019,9051,9082,9112,9142,9171,9200,9229,9256,9284,9311,
	9  9337,9363,9388,9413,9438,9462,9485,9508,9531,9554,9575,9597,
	1  9618,9639,9659,9680,9699,9719,9738,9756,9775,9793,9810,9828,
	2  9845,9861,9878,9894,9910,9926,9941,9956,9971,9985,9999/
C
	PARAMETER (NDMX=9,NB=2)
	REAL*4 PHT(NDMX),DET(NDMX),RAT(NDMX),DECT(NDMX),DST(NDMX)
	REAL*4 SIGT(NB,NDMX),MAGT(NB,NDMX)
	REAL*8 T0(NDMX)
	COMMON /NRMLPTS/ T0,H,DH,PER,MAGT,SIGT,DST,DET,PHT,RAT,DECT,NDM
C
	OMEGA = 2*PI/(3600*EXP(P(4)))
	INERTIA = EXP(P(5))
C
C       albedo prior
C	d(p_V)/d(ln(P_v)) = p_V^2*[f*/pk1^2*exp(-0.5*(p_V/pk1)^2+(1-f)/pk2^2*exp(-0.5*(p_V/pk2)**2)]
C
	PVMAX = 1.5
	IF (P(3).GT.23.) THEN
	  PV = PVMAX
	ELSE
	  PV = EXP(P(3))
	  PV = PVMAX*PV/(PVMAX+PV)
	ENDIF
	ERROR = -2*ALOG(PV**2*(FD*EXP(-0.5*(PV/PK1)**2)/PK1**2
	1	+(1-FD)*EXP(-0.5*(PV/PK2)**2)/PK2**2))
C	diameter prior just to stay in bounds
	OMEGA = 2*PI/(3600*EXP(P(4)))
	DIA = EXP(P(6))
	IF (P(6).GT.ALOG(1000.)) ERROR=ERROR+10.*(P(6)-ALOG(1000.))**2
	IF (P(6).LT.ALOG(.001)) ERROR=ERROR+10.*(P(6)-ALOG(.001))**2
C	period prior
	IF (PER.LE.0.) THEN
	  TIPSPEED = OMEGA*DIA*500.
	  TPSP0 = 0.065
	  IF (DIA.GT.0.2) TPSP0 = TPSP0*(DIA/0.2)**(5./6.)
	  ERROR = ERROR+2*ALOG(1+ALOG(TIPSPEED/TPSP0)**2)
C	periods shorter than 2 hrs for D>200 m strongly discouraged
	  IF (DIA.GT.0.2 .AND. OMEGA.GT.0.00087266)
	1	ERROR=ERROR+100*ALOG(OMEGA/0.00087266)**2
	ELSE
	  ERROR = ERROR+1.E4*(P(4)-ALOG(PER))**2
	ENDIF
C	inertia prior
	JLO = 2.5
	JHI = 2500.
C	no regolith on fast spinning small bodies
	IF (OMEGA.GT.0.00087266) JLO=250.
	IF (INERTIA.GT.JHI) ERROR=ERROR+ALOG(INERTIA/JHI)**2
	IF (INERTIA.LT.JLO) ERROR=ERROR+ALOG(JLO/INERTIA)**2
C
	XD = 53*(ALOG(DIA/.001)/ALOG(199.99))+1
	IF (DIA.GT.0.2) XD = 55+85*ALOG(DIA/.20001)/ALOG(1000/.20001)
	ID = XD
	FDIA = XD-ID
	IF (ID.GE.140) ID=139
	IF (ID.LE.0) ID=1
	PRIOR = 0.0001*(FDIA*PRE(ID+1)+(1-FDIA)*PRE(ID))
	ERROR = ERROR+2*ALOG(PRIOR)
	PJD_PV_PENALTY = ERROR
	RETURN
	END

	SUBROUTINE PREDICT(R,N,DE,P,THRM,ALB,EM,ASYM,ROT,NL,LAM,FNU)
C
C	computes the fluxes at NL wavelengths in LAM(1..NL)
C	in units of microns.  Fluxes in FNU(1..NL) in Jy
C	R(1..3) is the position of the observer relative to the Sun in AU
C	N(1..3) is a unit vector toward the object from the Earth
C	DE is the distance of the object from the Earth in AU
C	P is a unit vector giving the rotation pole of the object
C	THRM is the dimensionless thermal inertia parameter
C	ALB is the bolometric Bond albedo
C	EM is the emissivity
C	ASYM is C/A for the assumed prolate ellipsoidal shape
C	ROT is the angle of rotation, ROT=0 for the long axis in the x-z plane
C	all vectors have to be in the same coordinate system
C
C	equivalent volume diameter (C*A*A)^(1/3) is fixed at 1 km
C
	REAL R(3),N(3),P(3),THRM,ASYM,ROT
	REAL DE		! distance to Earth in AU
	INTEGER NL
	REAL*4 LAM(NL),FNU(NL)
C
	REAL OBJ(3),Z(3),X(3),Y(3)
	INTEGER NFACE,NLAT,NTIME,NLIST,NPXL
	INTEGER I, K
	PARAMETER (NFACE=127,NLAT=16,NTIME=32)
	PARAMETER (NPXL=12*4**6)
	REAL AXIS(3,NPXL)
	REAL PFACE(3,NPXL)
	INTEGER*2 FCLST(NPXL)
	REAL INSOL(NTIME,NFACE,NLAT)		! solar input
	REAL T(NTIME,NFACE,NLAT)	! the temperatures
	REAL LATS(NLAT)	! the latitudes
	REAL CL,SL	! cos and sin of latitude
	REAL WLAT(NLAT)	! weight for latitude integration
	REAL SUNLAT,SSL,CSL	! latitude of Sun and its sine and cosine
	INTEGER ISL
	REAL ELONGATION,CELONG,SELONG	! elongation & cos & sin
	REAL THMAX	! max angle of crater walls
	REAL ALPHA	! phase angle for EVEC calculations
	REAL PHASE	! phase angle in degrees
	REAL ELHA	! local hour angle of EVEC
	REAL ELAT	! latitude of Earth
	REAL G(NTIME)	! normalized thermal response
	REAL TBRITE(12)
	REAL TMP	! temporary variable used to rotate vectors
	REAL THCC /1.19104272E-5/
	REAL OMEGA
	REAL HCK /1.4387752/	! h*c/k
	REAL XB		! h*nu/k*Tb
	REAL SOLCON /1367.E3/	! solar constant in erg/cm^2/sec
	PARAMETER (TSUN = 5600.)
	REAL SIGSB /5.6704E-5/	! Stefan-Boltzmann constant
	REAL TAT1AU /394.03854311248/
	REAL ALBEDO, EMISSIVITY
	REAL ETA,DIA
	REAL T0
	REAL DS		! distance to Sun in AU
	REAL AU /1.495979E13 /	! AU in cm
	REAL LXT(9) /-4.,-3.,-2.,-1.,0.,1.,2.,3.,4./
	INTEGER NLMX
	PARAMETER (NLMX=99)
	REAL TBT2(9),LX(NLMX)  ! variables for interpolating TB
	REAL FSUM(NLMX),F,FLAM(NLMX)
	REAL VSUM,FOPT
	INTEGER IDS,IA,IDE,ILAM
	REAL JY,MAB
	REAL SPLINE
	REAL PI,TWOPI,DTR
	PARAMETER (PI=3.1415927,TWOPI=6.2831853,DTR=57.29578)
	PARAMETER (Q=0.38399,GG=0.15,A1=3.33,B1=0.63,A2=1.87,B2=1.22)
	PARAMETER (PIR2PV = 1.6)
C
	INTEGER IMC,NEE
	REAL RMS
	REAL*8 SE,SSEE,SEE
	INTEGER IY,IZ,IR,JR
	REAL RAN
C
C
C	OBJ is the object wrt the Sun, R is the Earth wrt the Sun
C
	DO I=1,3
	  Z(I) = P(I)
	  OBJ(I) = DE*N(I)+R(I)
	ENDDO
	DS = SQRT(DOT(OBJ,OBJ))
	CALL NORM(OBJ)
	CALL NORM(Z)
	CSL = -DOT(Z,OBJ)
	SUNLAT = 90/DTR-ACOS(CSL)
C	set up object based unit vectors with z=P and Sun in x-z plane
C	OBJ is close to -X so do y = (-x) X (z)
	CALL CROSS(OBJ,P,Y)
	CALL NORM(Y)
C	Y is unit vector perp to POLE and OBJ
	CALL CROSS(Y,Z,X)
	XE = -DOT(X,N)
	YE = -DOT(Y,N)
	ZE = -DOT(Z,N)
	ELHA = ATAN2(YE,XE)
	ELAT = ATAN2(ZE,SQRT(XE**2+YE**2))
	ALPHA = ACOS(SIN(SUNLAT)*SIN(ELAT)+
	1	COS(SUNLAT)*COS(ELAT)*COS(ELHA))
	PHASE = DTR*ALPHA
C       compute the H-G formulation phase function
	TA2 = TAN(0.5*PH/DTR)
	PHI1 = EXP(-A1*TA2**B1)
	PHI2 = EXP(-A2*TA2**B2)
	PHI = (1-GG)*PHI1+GG*PHI2
C	WRITE (*,FMT='(3F9.4,3X,3F9.5)') OBJ,Z
C	WRITE (*,FMT='(F9.4,3F9.3,F9.1)') DS,SUNLAT,ELAT,ELHA,PHASE
C------------
	CALL MAKE_G(NTIME,G)
C
	THMAX = 45
	THMAX=THMAX/DTR
	ALBEDO = ALB
	EMISSIVITY = EM
	DO I=1,NLAT
	  SL = -1.+(2.*I-1.)/NLAT
	  CL = SQRT(1.-SL**2)
	  LATS(I) = ATAN(SL/CL)
	  WLAT(I) = 2./NLAT
	ENDDO
C	set up axes of the crater facets
	CALL SPHERE(NFACE,THMAX,AXIS,PFACE,FCLST)
	DO K=1,NLAT
C	compute the solar input
	  CALL SOLAR(SUNLAT,LATS(K),NFACE,NTIME,THMAX,PFACE,AXIS,
	1	FCLST,ALBEDO,INSOL(1,1,K))
	  CALL THERMAL(NFACE,NTIME,THRM,G,THMAX,INSOL(1,1,K),T(1,1,K))
	ENDDO
C
	CSL = COS(SUNLAT)
	SSL = SIN(SUNLAT)
C	solid angle of 1 km diameter at DE AU
	OMEGA = (PI/4)*(1.E5/(DE*AU))**2
	ALPHA = PHASE/DTR
	T0 = TAT1AU*((1-ALBEDO)/EMISSIVITY/DS**2)**0.25
	DO ILAM=1,NL
	  LX(ILAM) = ALOG(HCK*10000.0/(LAM(ILAM)*T0))/ALOG(2.)
	ENDDO
C	calculate the fluxes seen at Earth
	CALL FLUXES(ELAT,ELHA,NLAT,LATS,WLAT,NFACE,NTIME,AXIS,
	1	PFACE,FCLST,INSOL,T,TBRITE)
	DO I=1,12
	  IF (TBRITE(I).LE.0) THEN
	    WRITE (*,*) 'TBRITE < 0'
	    WRITE (*,*) ELAT,ELHA,SUNLAT,THRM
	    WRITE (*,*) TBRITE
	    IF(I.EQ.12) TBRITE(12)=0
	    IF (I.NE.12) STOP
	  ENDIF
	ENDDO
	VSUM = TBRITE(12)*(ALBEDO/(1-ALBEDO))
	CALL SPLSET(9,LXT,TBRITE,TBT2)
	DO ILAM = 1,NL
	  XB = HCK*10000./
	1	(LAM(ILAM)*T0*SPLINE(LX(ILAM),9,LXT,TBRITE,TBT2))
	  F = 0.
	  IF (XB.LT.70.) F = 1./(EXP(XB)-1.)
	  FSUM(ILAM) = THCC*OMEGA*(10000./LAM(ILAM))**4*EMISSIVITY*F
	  XS = HCK*10000./(TSUN*LAM(ILAM))
C       the solar flux, nu*F_\nu, in erg/cm^2/sec/oct
	  FSUN = (XS**4/(EXP(XS)-1))/((PI**4)/15)*SOLCON/DS**2
C	face-on white disk flux at 1 AU for 1 km is
C	(PI*R**2)*FSUN/(PI*AU**2) = FSUN*(R/AU)**2
	  PV = ALBEDO/Q
	  PIR=PV*EXP(MAX(0.,MIN(1.,LOG(LAM(ILAM)/0.55)/LOG(3.4/0.55)))*
	1	LOG(PIR2PV))
	  FHG = PHI*PIR*FSUN*(5.E4/(DE*AU))**2
	  FOPT = THCC*OMEGA*(10000./LAM(ILAM))**4*VSUM*
	1	2.5E-5/DS**2/(EXP(14387.752/LAM(ILAM)/5600.)-1)
	  FSUM(ILAM) = FSUM(ILAM)+FHG
	  JY = FSUM(ILAM)*LAM(ILAM)/3.E-9
	  FLAM(ILAM) = JY
	  FNU(ILAM) = JY
	ENDDO
C------------
C
	RETURN
	END

	SUBROUTINE MAKE_G(NTIME,G)
C	computes the flux into the surface for a unit triangle wave
C	temperature history
	IMPLICIT NONE
	INTEGER NTIME		! number of times day is divided into
	REAL G(NTIME)		! output array
	REAL PI,TWOPI,DTR
	PARAMETER (PI=3.1415927,TWOPI=6.2831853,DTR=57.29578)
	INTEGER N,I
	REAL TH,YN
	DO I=1,NTIME
	  G(I) = 0.
	ENDDO
	DO N=1,(20*NTIME)
	  YN =  NTIME*(1.0-COS(TWOPI*N/NTIME))*SQRT(N/2.)/(PI*N)**2
	  DO I=1,NTIME
	    TH = TWOPI*(I-1)/NTIME
	    G(I) = G(I) + YN*(SIN(N*TH)-COS(N*TH))
	  ENDDO
	ENDDO
C	WRITE (*,FMT='(8F9.4)') G
	RETURN
	END

	SUBROUTINE SPHERE(NFACE,THMAX,AXIS,PFACE,FCLST)
C
C	sets up AXIS and PFACE for a spherical cap crater
C
	IMPLICIT NONE
	INTEGER NFACE		! input number of facets
	REAL THMAX		! input max angle of wall, radians
	INTEGER NPXL
	PARAMETER (NPXL=12*4**6)
	REAL AXIS(3,NPXL)	! output: normal vectors for facets
	REAL PFACE(3,NPXL)	! output: centers of facets
	INTEGER*2 FCLST(NPXL)	! output, FACE for each pixel, 0 for none
	REAL CMAX,SMAX		! sine and cosine of THMAX
	INTEGER I,M,N		! counters
	INTEGER NRING
	INTEGER LONG
	REAL CLIM(6)		! cos for top of n'th ring
	INTEGER WF(127)
	REAL PI,TWOPI,DTR
	PARAMETER (PI=3.1415927,TWOPI=6.2831853,DTR=57.29578)
	CMAX = COS(THMAX)
	SMAX = SIN(THMAX)
	IF (NFACE.EQ.1) THEN
	  NRING = 0
	ELSE IF (NFACE.EQ.7) THEN
	  NRING = 1
	ELSE IF (NFACE.EQ.19) THEN
	  NRING = 2
	ELSE IF (NFACE.EQ.37) THEN
	  NRING = 3
	ELSE IF (NFACE.EQ.61) THEN
	  NRING = 4
	ELSE IF (NFACE.EQ.91) THEN
	  NRING = 5
	ELSE IF (NFACE.EQ.127) THEN
	  NRING = 6
	ELSE
	  STOP 'bad NFACE'
	ENDIF
C	WRITE (*,*) NRING,NFACE
	N = 1
	DO M=1,NRING
	  CLIM(M) = 1-(1-CMAX)*FLOAT(N)/FLOAT(NFACE)
	  N = N+6*M
	ENDDO
C	WRITE (*,FMT='(8F9.6)') (CLIM(M),M=1,NRING),CMAX
	DO I=1,NPXL
	  CALL NESTED_CENPIX(I-1,6,AXIS(1,I))
	  IF (AXIS(3,I).LT.CMAX) THEN
	    FCLST(I) = 0
	  ELSE
	    M = NRING
20	    IF (AXIS(3,I).LT.CLIM(M)) GO TO 30
	    M = M-1
	    IF (M.GE.1) GO TO 20
30	    IF (M.EQ.0) THEN
	      FCLST(I) = 1
	    ELSE
	      LONG = 6*M*(ATAN2(AXIS(2,I),AXIS(1,I))+PI)/TWOPI
	      IF (LONG.GE.6*M) LONG = 6*M-1
C	2+3*M*(M-1) = 2,8,20,38,62,92 ... for M=1,2,3,4,5,6
	      FCLST(I) = LONG + 2+3*M*(M-1)
	    ENDIF
	  ENDIF
	  CALL SET_VEC(-AXIS(1,I)/SMAX,-AXIS(2,I)/SMAX,
	1	(CMAX-AXIS(3,I))/SMAX,PFACE(1,I))
	ENDDO
	DO I=1,NPXL
	  IF (FCLST(I).GT.0) WF(FCLST(I)) = WF(FCLST(I))+1
	ENDDO
C	WRITE (*,FMT='(16I5)') (WF(I),I=1,NFACE)
C	WRITE (2) FCLST
	RETURN
	END

	SUBROUTINE SET_VEC(X,Y,Z,V)
	REAL X,Y,Z,V(3)
	V(1) = X
	V(2) = Y
	V(3) = Z
	RETURN
	END

	SUBROUTINE LWATE(ELAT,ELHA,LAT,NFACE,AXIS,PFACE,FCLST,WATE)
C	routine to compute Weights for the NFACE FACETS when the Earth is
C	at ELAT and ELHA.  Surface is at LAT
C	mouth of crater is circle of radius 1
	PARAMETER (NPXL=12*4**6)
	REAL AXIS(3,NPXL)
	REAL PFACE(3,NPXL)
	INTEGER*2 FCLST(NPXL)
	REAL WATE(NFACE)
	INTEGER WF(127)
	REAL E(3)
	REAL ELAT
	REAL ELHA
	REAL LAT
C
	CL = COS(ELAT)
	SL = SIN(ELAT)
	CALL SET_VEC(CL*COS(ELHA),CL*SIN(ELHA),SL,E)
C
	CL = COS(LAT)
	SL = SIN(LAT)
	TMP = CL*E(1)+SL*E(3)
	E(1) = CL*E(3)-SL*E(1)
	E(3) = TMP
	DO I=1,NFACE
	  WATE(I) = 0.
	  WF(I) = 0
	ENDDO
C	see if Earth is above horizon
	IF (E(3).LE.1.E-5) RETURN
C	above horizon so compute the weights
	SXZ = E(1)/E(3)
	SYZ = E(2)/E(3)
	DO I=1,NPXL
C	compute position where LOS crosses Z=0 surface
	  IF (FCLST(I).GT.0) THEN
	    WF(FCLST(I)) = WF(FCLST(I))+1
	    X = PFACE(1,I)-SXZ*PFACE(3,I)
	    Y = PFACE(2,I)-SYZ*PFACE(3,I)
	    RR = X*X+Y*Y
	    IF (RR.LE.1.) 
	1	WATE(FCLST(I)) = WATE(FCLST(I))+MAX(0.,DOT(E,AXIS(1,I)))
	  ENDIF
	ENDDO
C	normalize to cos(angle of incidence of surface) = E(3)
	DO I=1,NFACE
	  IF (WATE(I).LT.0.) GO TO 20
	  WATE(I) = WATE(I)/WF(I)
	ENDDO
	GO TO 40
20	WRITE (*,*) 'negative weights in LWATE'
	WRITE (*,FMT='(10F8.5)') WATE
40	SUM = 0
	DO I=1,NFACE
	  SUM = SUM+WATE(I)
	ENDDO
	IF (SUM.LE.0.) THEN
C	  WRITE (*,*) ' SUM LE 0 in LWATE'
	  RETURN
	ENDIF
	SUM = E(3)/SUM
	DO I=1,NFACE
	  WATE(I) = SUM*WATE(I)
	ENDDO
	RETURN
	END

	SUBROUTINE SOLAR(SUNLAT,LAT,NFACE,NTIME,THMAX,PFACE,AXIS,
	1	FCLST,ALBEDO,INSOL)
C	computes Solar input vs time and facet
	PARAMETER (NPXL=12*4**6)
	REAL AXIS(3,NPXL)
	REAL PFACE(3,NPXL)
	INTEGER*2 FCLST(NPXL)
	REAL ALBEDO
	REAL INSOL(NTIME,NFACE)		! solar input
	REAL PI,TWOPI,DTR
	PARAMETER (PI=3.1415927,TWOPI=6.2831853,DTR=57.29578)
	REAL LAT
	PARAMETER (NFMAX=127)
	REAL W(NFMAX)
C
	AFACE = NFACE*(1.+COS(THMAX))/2.
	D = (1.-COS(THMAX))/2./NFACE
C	weights from LWATE peak at 1/NFACE
C	(1+cos(thmax))/2 is (area of crater mouth)/(area of crater wall)
C
	DO I=1,NTIME
	  HA = ((I-0.5)/NTIME-0.5)*TWOPI
	  CALL LWATE(SUNLAT,HA,LAT,NFACE,AXIS,PFACE,FCLST,W)
	  TOTAL = 0
	  DO J=1,NFACE
	    SOL = AFACE*W(J)
	    TOTAL = TOTAL+SOL
	    INSOL(I,J) = SOL
	    IF (SOL.GT.1.25) THEN
	      WRITE (*,*) 'TOO BRIGHT IN SOLAR',I,J,SOL
	      INSOL(I,J) = 1.25
	    ENDIF
	  ENDDO
	  AVERAGE = TOTAL/NFACE
	  F = ALBEDO*(1-COS(THMAX))/2.
C	add in scattered sunlight
	  DO J=1,NFACE
	    INSOL(I,J) = (1-ALBEDO)*(INSOL(I,J) + F*AVERAGE/(1-F))
	  ENDDO
	ENDDO
	RETURN
	END

	SUBROUTINE THERMAL(NFACE,NTIME,X,G,THMAX,INSOL,T)
C	does the thermal inertia calculation
	REAL INSOL(NTIME,NFACE)		! solar input
	REAL T(NTIME,NFACE)	! the temperatures
	PARAMETER (NFMAX=127,NTMAX=32)
	PARAMETER (XMAX=10000.0)
	REAL A(NTMAX,NTMAX+1)
	REAL MUTUAL(NFMAX)
	REAL G(NTIME)
C
	D = (1.-COS(THMAX))/2./NFACE
	DO I=1,NFACE
	  MUTUAL(I) = D
	ENDDO
C	initial guess
	DO J=1,NFACE
	  S = 0.
	  DO I=1,NTIME
	    S = S+INSOL(I,J)
	  ENDDO
	  S = S/NTIME
	  S = S/16.
	  DO I=1,NTIME
	    T(I,J) = MAX(S,0.003,INSOL(I,J))**0.25
	  ENDDO
	ENDDO
C	iteration
	NIT = 1
100	EPSILON = 0
	DO JJ=1,NFACE
	  DO I=1,NTIME
C	infrared emission - solar input
	    S = T(I,JJ)**4-INSOL(I,JJ)
C
C	add in emission from other facets
C
	    DO K=1,NFACE
	      S = S-D*T(I,K)**4
	    ENDDO
C
C	set up partial derivatives
C
	    DO J=1,NTIME
	      K = I-J+1
	      IF (K.LE.0) K = K+NTIME
	      S = S-X*G(K)*T(J,JJ)
	      A(I,J) = -X*G(K)
	    ENDDO
	    A(I,NTIME+1) = S
	    A(I,I) = A(I,I)+4*T(I,JJ)**3
C	protect against large X
	    IF (X.GT.XMAX) A(I,I) = A(I,I)+10*X/XMAX
	  ENDDO
C	solve for corrections on this facet
	  CALL MATINV(A,NTMAX,NTIME,1,IFF)
	  IF (IFF.NE.0) STOP 'bad IFF'
	  ER = 0.
	  DO I=1,NTIME
	    ER = MAX(ER,ABS(A(I,NTIME+1))/T(I,JJ))
	  ENDDO
C	limit max change to 25%
	  F = 1.
	  IF (ER.GT.0.25) F = 0.25/ER
	  DO I=1,NTIME
	    DELTA = F*A(I,NTIME+1)
	    T(I,JJ) = T(I,JJ) - DELTA
	    EPSILON = MAX(EPSILON,ABS(DELTA))
	  ENDDO
	ENDDO	! JJ loop over facets
	NIT = NIT+1
C	put iteration limit at 500
	IF (EPSILON.GT.0.001 .AND. NIT.LT.500) GO TO 100
	IF (NIT.GE.500) WRITE (*,*) 'NIT,X=',NIT,X
	RETURN
	END

	SUBROUTINE FLUXES(ELAT,ELHA,NLAT,LATS,WLAT,NFACE,NTIME,AXIS,
	1	PFACE,FCLST,INSOL,T,TBRITE)
C	routine to compute the fluxes at Earth viewing point
C	ELAT and ELHA are the sub-Earth latitude and local hour angle
C	in radians.  +ELHA means afternoon, -ELHA is morning
C	NLAT is the number of latitudes, listed in LATS, with weights WLAT
C	weight should be COS(LAT)*DLAT
C	output is array TBRITE
C	1-9 is h*nu/k*To = 1/16,1/8,...8,16
C	10 is the bolometric or effective temperature
C	11 is the radio temperature for deep penetration and low dielectric
C	constant
C	12 is the result for optical light scaled so 
C	ALBEDO/(1-ALBEDO)*TBRITE(12) gives 1 for a white diffuse
C
	PARAMETER (NFMAX=127,NTMAX=32)
	REAL TBAR(NFMAX),WATE(NFMAX)
	PARAMETER (NPXL=12*4**6)
	REAL AXIS(3,NPXL)
	REAL PFACE(3,NPXL)
	INTEGER*2 FCLST(NPXL)
	REAL T(NTIME,NFACE,NLAT)	! the temperatures
	REAL TBRITE(12)
	REAL INSOL(NTIME,NFACE,NLAT)         ! solar input
	REAL LATS(NLAT)	! the latitudes
	REAL WLAT(NLAT)	! weight for latitude integration
	REAL*8 X
	REAL PI,TWOPI,DTR
	PARAMETER (PI=3.1415927,TWOPI=6.2831853,DTR=57.29578)
C
	WSUM = 0.	
	DO I=1,12
	  TBRITE(I) = 0.
	ENDDO
	DO K=1,NLAT
	  DO J=1,NFACE
	    S = 0.
	    DO I=1,NTIME
	      S = S+T(I,J,K)
	    ENDDO
	    TBAR(J) = S/NTIME
	  ENDDO
	  DO I=1,NTIME
	    HA = ((I-0.5)/NTIME - 0.5)*TWOPI-ELHA
	    CALL LWATE(ELAT,HA,LATS(K),NFACE,AXIS,PFACE,FCLST,WATE)
	    SW = 0.
	    DO J=1,NFACE
	      SW = SW+WATE(J)
	    ENDDO
	    IF (SW.GT.0.00001) THEN
	      DO J=1,NFACE
	        W = WLAT(K)*WATE(J)
	        IF (T(I,J,K).GT.0.) THEN
	          TBRITE(10) = TBRITE(10) + W*T(I,J,K)**4
	          TBRITE(11) = TBRITE(11) + W*TBAR(J)
	          TBRITE(12) = TBRITE(12) + W*INSOL(I,J,K)
	          X = 0.0625/T(I,J,K)
	          DO L=1,9
	            IF (X.LE.70.) TBRITE(L) = TBRITE(L) + W/(EXP(X)-1.)
	            X = X+X
	          ENDDO
	        ENDIF
	        WSUM = WSUM+W
	      ENDDO	! J=1,NFACE
	    ENDIF
	  ENDDO	! I=1,NTIME
	ENDDO	! K=1,NLAT
C	normalize
	DO I=1,12
	  TBRITE(I) = TBRITE(I)/WSUM
	ENDDO
C	convert from flux to brightness T
	TBRITE(10) = TBRITE(10)**0.25
	X = 0.0625
	DO I=1,9
	  IF (TBRITE(I).GT.0.) TBRITE(I) = X/ALOG(1.+1./TBRITE(I))
	  IF (I.GE.3 .AND. TBRITE(I).LE.0.) 
	1	TBRITE(I) = 2*TBRITE(I-1)-TBRITE(I-2)
	  X = X+X
	ENDDO
	RETURN
	END

	SUBROUTINE MATINV(A,NROW,N,NRHS,IFF)
C	single precision version of DMATINV
	DIMENSION A(NROW,*)
	M=N+NRHS
	NP=N+1
	NM=N-1
	IFF=0
	DO 40 I=1,NM
C	FIND PIVOT
	V=ABS(A(I,I))
	IP=I
	DO 10 J=I,NM
	IF(ABS(A(J+1,I)).LE.V) GO TO 10
	V=ABS(A(J+1,I))
	IP=J+1
10	CONTINUE
	IF(V.LE.0.) GO TO 666
	IF(IP.EQ.I) GO TO 25
C	INTERCHANGE ROWS TO PUT PIVOT ON DIAGONAL
	DO 20 J=I,M
	V=A(I,J)
	A(I,J)=A(IP,J)
20	A(IP,J)=V
25	IP=I+1
C	LOOP OVER ROWS BELOW THE DIAGONAL
	DO 40 IL=IP,N
	F=-A(IL,I)/A(I,I)
	A(IL,I)=0.
C	REDUCE THE ELEMENTS IN THE ROW
	DO 40 J=IP,M
40	A(IL,J)=A(IL,J)+A(I,J)*F
C	BACK SUBSTITUTION
	I=NP
	DO 60 II=1,N
	I=I-1
	IP=I+1
	DO 60 J=NP,M
	S=0.
	IF(I.EQ.N) GO TO 60
	DO 50 K=IP,N
50	S=S+A(K,J)*A(I,K)
60	A(I,J)=(A(I,J)-S)/A(I,I)
	RETURN
666	IFF=666
	RETURN
	END

      FUNCTION BRENT(AX,BX,CX,F,TOL,XMIN)
      PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0E-10)
      A=MIN(AX,CX)
      B=MAX(AX,CX)
      V=BX
      W=V
      X=V
      E=0.
      FX=F(X)
      FV=FX
      FW=FX
      DO 11 ITER=1,ITMAX
        XM=0.5*(A+B)
        TOL1=TOL*ABS(X)+ZEPS
        TOL2=2.*TOL1
        IF(ABS(X-XM).LE.(TOL2-.5*(B-A))) GOTO 3
        IF(ABS(E).GT.TOL1) THEN
          R=(X-W)*(FX-FV)
          Q=(X-V)*(FX-FW)
          P=(X-V)*Q-(X-W)*R
          Q=2.*(Q-R)
          IF(Q.GT.0.) P=-P
          Q=ABS(Q)
          ETEMP=E
          E=D
          IF(ABS(P).GE.ABS(.5*Q*ETEMP).OR.P.LE.Q*(A-X).OR. 
	1	P.GE.Q*(B-X)) GOTO 1
          D=P/Q
          U=X+D
          IF(U-A.LT.TOL2 .OR. B-U.LT.TOL2) D=SIGN(TOL1,XM-X)
          GOTO 2
        ENDIF
1       IF(X.GE.XM) THEN
          E=A-X
        ELSE
          E=B-X
        ENDIF
        D=CGOLD*E
2       IF(ABS(D).GE.TOL1) THEN
          U=X+D
        ELSE
          U=X+SIGN(TOL1,D)
        ENDIF
        FU=F(U)
        IF(FU.LE.FX) THEN
          IF(U.GE.X) THEN
            A=X
          ELSE
            B=X
          ENDIF
          V=W
          FV=FW
          W=X
          FW=FX
          X=U
          FX=FU
        ELSE
          IF(U.LT.X) THEN
            A=U
          ELSE
            B=U
          ENDIF
          IF(FU.LE.FW .OR. W.EQ.X) THEN
            V=W
            FV=FW
            W=U
            FW=FU
          ELSE IF(FU.LE.FV .OR. V.EQ.X .OR. V.EQ.W) THEN
            V=U
            FV=FU
          ENDIF
        ENDIF
11    CONTINUE
      WRITE (*,*) 'Brent exceed maximum iterations.'
3     XMIN=X
      BRENT=FX
      RETURN
      END
