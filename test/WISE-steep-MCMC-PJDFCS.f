C	WISE-steep-MCMC-PJDFCS.f
C	gfortran -o WISE-steep-MCMC-PJDFCS WISE-steep-MCMC-PJDFCS.f READ_TABLE-steep.f ~/Ned.a
C
C	modified 30-DEC-2017 to use 5933 for TSUN, the V-W1 color temperature
C	modified 21-DEC-2017 to properly account for NP in prob factor
C	modified 08-SEP-2017 to use 10 parameters, with P(9)=logit([b/a]^4)
C	and P(10) = logit([c/b]^4).  A uniform prior in [b/a]^4 and [c/b]^4
C	is used.
C
C	Lightcurve amplitudes are input if a period is specified
C	A second line with values for all 10 parameters is a new part of
C	the input.  Zero means not specified.  If all 10 parameters are
C	specified, a program just computes the model and stops - no MCMC. 
C
C	modified 11-MAY-2017 to use 8 parameters, with P(7)=logit(fc)
C	the crater fraction, and P(8)=ln(p_IR/p_V)
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
C	no writes to fort.23
C
C	READ_TABLE will read the tables made by make-tables-steep.f
C
	PARAMETER (NFACE=127)
	PARAMETER (NPXL=12*4**6)
	REAL AXIS(3,NPXL)
	REAL PFACE(3,NPXL)
	INTEGER*2 FCLST(NPXL)
	COMMON /SPHERE_TABLE/ AXIS,PFACE,FCLST
C
	INTEGER NP,NDMX,NB,NFP
	LOGICAL POLEKNOWN
	PARAMETER (NP=10,NDMX=99,NB=4)
	REAL*4 PHT(NDMX),DET(NDMX),RAT(NDMX),DECT(NDMX),DST(NDMX)
	REAL*4 SIGT(NB,NDMX),MAGT(NB,NDMX),AMPT(NB,NDMX)
	REAL*8 T0(NDMX)
	COMMON /NRMLPTS/ T0,H,DH,MAGT,SIGT,AMPT,DST,DET,PHT,RAT,DECT,NDM
	REAL PRIOR(NP)
	COMMON /PREDICTED/ WMAGT(NB,NDMX),WAMPT(NB,NDMX)
	REAL*8 MJD,JD(NDMX)
	REAL*8 T1JAN2010,T,JD1JAN2010
	REAL*8 TSTART,TWALL
	REAL*8 JULIAN_DAY,GMT2T68
	REAL*4 EA(3),A(3),E(3)
	REAL*4 WMAG(NB),SIG(NB),AMP(NB)
	REAL*4 P(NP)
	REAL*4 P1(3),P2(3),P3(3)
	CHARACTER*11 GMT
	CHARACTER*10 NAME
C
	COMMON /WATERMARK/ BESTEVER,PB(NP)
	PARAMETER (NW=200)
	REAL*4 PW(NP,NW),CH2(NW)
	INTEGER NTRY(NW)
	INTEGER PPMIL
	COMMON /WALKERS/ PW,CH2,NTRY
C
	REAL*4 SIGMIN(NB) /4*0.02/
	EXTERNAL RAN
C
	PARAMETER (MCTMX=72900)
	PARAMETER (DTR= 57.2957795)
	REAL MCT(NP+2,MCTMX)
C
	CALL RANDOMIZE
	CALL READ_TABLE
C	set up axes of the crater facets
	THMAX = 75/DTR
	CALL SPHERE(NFACE,THMAX,AXIS,PFACE,FCLST)
C
	BESTEVER = 1.E30
C
	NDM=0
	WRITE (*,FMT='(A,$)') 'Enter H,sig(H),Period,Name:'
	READ (*,FMT='(3F7.3,A)') H,DH,PER,NAME
	WRITE (*,FMT='(3F7.3,1X,A)') H,DH,PER,NAME
	WRITE  (*,FMT='(A,$)') 'Enter fixed:'
	READ (*,FMT='(10E8.3)') PRIOR
	IF (PER.NE.0.) PRIOR(4) = ALOG(PER)
	NPRIOR = 0
	DO I=1,NP
	  IF (PRIOR(I).NE.0.) NPRIOR = NPRIOR+1
	ENDDO
	POLEKNOWN = PRIOR(1).NE.0 .OR. PRIOR(2).NE.0
	NFP = NFP-NPRIOR
	WRITE (*,FMT='(I2,10F8.4)') NPRIOR,PRIOR
	WRITE (*,*) POLEKNOWN,NFP
C
	CALL TIME_TAG(2010,1,1,0,0,0.,T1JAN2010)
	JD1JAN2010 = JULIAN_DAY(T1JAN2010)
	WRITE (*,FMT='(A)') 
	1	'Enter MJD,RA,Dec,Delta,W1234 with sigmas & amps:'
	WRITE (*,FMT='(A)') '   MJD      RA      DEC     DE      DS'//
	1	'    PHASE W1234, sigmas & amps'
10	READ (*,900,END=100) MJD,RA,DEC,DELTA,
	1	(WMAG(I),SIG(I),AMP(I),I=1,4)
900	FORMAT(F11.5,2F11.7,F8.4,4(F7.3,F6.3,F5.2))
	IF (NDM.GE.NDMX) STOP '2 much data'
	NDM=NDM+1
	DO I=1,4
	  MAGT(I,NDM) = WMAG(I)
	  SIGT(I,NDM) = SQRT(SIG(I)**2+SIGMIN(I)**2)
	  AMPT(I,NDM) = AMP(I)
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
	1	PHT(NDM),(WMAG(I),SIG(I),AMP(I),I=1,4)
901	FORMAT(F9.3,2F8.3,2F8.5,F7.2,4(F7.3,F6.3,F5.2))
	GO TO 10
100	CLOSE (UNIT=1)
	WRITE (*,*) NDM
C	shortcut for just prediction
	IF (NPRIOR.EQ.NP) THEN
	  WRITE (*,*) 'prediction only short cut'
	  WRITE (*,FMT='(10F8.4)') PRIOR
	  DO I=1,NP
	    P(I) = PRIOR(I)
	  ENDDO
	  E = ERROR(P,1.E30)
	  STOP
	ENDIF
C
C	Monte Carlo Markov chain
C
	WRITE (21,*) ' Monte Carlo Markov Chain'
	CLOSE (UNIT=21)
	TSTART = GMT2T68('NOW')
	IMC = 0
	TWALL = GMT2T68('NOW')-TSTART
	WRITE (*,FMT='(I4,F10.0,2H P)') IMC,TWALL
	WRITE (*,*)
	DO I=1,NW
	  DO J=1,NP
	    PW(J,I) = PRIOR(J)
	  ENDDO
	  IF (PRIOR(1).EQ.0. .AND. PRIOR(2).EQ.0.) THEN
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
	  ENDIF
C	albedo uniform in log 0.01 to 1.01
	  IF (PRIOR(3).EQ.0.) PW(3,I) = ALOG(RAN(IR,JR)+0.01)
C	periods 0.0001 hours to 1000 hours
	  IF (PRIOR(4).EQ.0.) PW(4,I)=-4*ALOG(10.)+7*ALOG(10.)*RAN(IR,JR)
C	J = sqrt(C\rho\kappa) from 5 to 2500, uniform in log
	  IF (PRIOR(5).EQ.0.) PW(5,I) = ALOG(5.)+RAN(IR,JR)*ALOG(500.)
C	fc, color & shape priors
	  DO J=7,10
	    IF (PRIOR(J).EQ.0.) THEN
	      IF (J.NE.8) THEN
C	fc uniform in 0.0001 to 0.9999
	        FC = 0.0001+0.9998*RAN(IR,JR)
	        PW(J,I) = ALOG(FC/(1-FC))
	      ELSE
C	p_IR/p_V log normal
C	prior should be ln(PIR2PV)=0.288+/-0.340, see apj498_786_t1_mrt.f
C	also apj408731t1_mrt.txt Mainzer etal (2001) gives 0.563 for mean
	        PW(8,I) =   EXP(0.563+0.340*GAUSRN(IR,JR))
	      ENDIF
	    ENDIF
	  ENDDO
C	diameter from 1 meter to 1000 km
	  IF (PRIOR(6).EQ.0.) THEN
	    IF (I.EQ.1) THEN
	      PW(6,I) = ALOG(.001)+ALOG(1.E6)*RAN(IR,JR)
	    ELSE
	      PW(6,I) = PW(6,1)+0.5*GAUSRN(IR,JR) 
	    ENDIF
	    DO J=1,NP
	      P(J) = PW(J,I)
	    ENDDO
	    D1 = EXP(P(6))
	    CH2(I) = ERROR(P,1.E30)
	    IMC=IMC+1
	    IF (I.EQ.1) THEN
	      E1 = CH2(I)
C	reset diameter for better start
	      SWDD = 0
	      SW = 0
	      IF (DH.LT.8.) THEN
	        SW = SW+1
C	SDD is the ln[(D from H & PV)/D]
	        SWDD=SWDD+
	1	  ALOG(1329./10.**(H/5)/SQRT(EXP(PW(3,I))))-PW(6,I)
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
	      IF (E1.NE.E1) THEN
	        WRITE (*,*) 'NaN while setting D1'
	        P(6) = 0.
	      ELSE
	        P(6) = (SWDD/SW)+P(6)
	      ENDIF
	      PW(6,I) = P(6)
	      D2 = EXP(P(6))
	      CH2(I) = ERROR(P,1.E30)
	      IMC = IMC+1
	      E2 = CH2(I)
C	write (*,*) d1,e1,d2,e2
	    ENDIF
	  ENDIF
C
	  NTRY(I) = 0
	  WRITE (*,FMT='(I3,10F8.3,F15.2)') I,(PW(J,I),J=1,NP),CH2(I)
	ENDDO
	WRITE (*,*)
	NMCT = 0
	MCTARGET = 300
200	  DO KW=1,NW
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
	    DO I=3,NP
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
	    CHSQ = ERROR(P,CH2(KW)+100)
	    IMC = IMC+1
	    IF (MOD(IMC,1000).EQ.0) 
	1	WRITE (*,FMT='(I7,I6,F10.0)') IMC,NMCT,TWALL
	    Q = 0
	    XP = -0.5*(CHSQ-CH2(KW))
	    IF (XP.GT.50.) Q = 1
	    IF (POLEKNOWN) THEN
	      ZFACTOR = Z**(NFP-1)
	    ELSE
C	Z**(NP-2) was Z**4 until 21-DEC-2017
	      ZFACTOR = ABS(SZTH/SIN(TH))*Z**(NP-2)
	    ENDIF
	    IF (ABS(XP).LT.50.) 
	1	Q = MIN(1.0,ZFACTOR*EXP(XP))
	    U = RAN(IR,JR)
	    IF (U.LT.Q) THEN
C	accept
C
	      NMCT = NMCT+1
	      DO I=1,NP
	        MCT(I,NMCT) = PW(I,KW)
	      ENDDO
	      MCT(NP+1,NMCT) = CH2(KW)
	      MCT(NP+2,NMCT) = NTRY(KW)+0.5
	      IF (NMCT.GE.MCTARGET) THEN
	        NBURN = MCTARGET/3
	        WRITE (21,*) ' Monte Carlo Markov Chain -',NBURN,' burnt'
	        DO J=NBURN+1,MCTARGET
	          MTRY = MCT(NP+2,J)
	          WRITE (21,FMT='(F8.5,9F9.5,F17.5,I5)')
	1		(MCT(I,J),I=1,NP+1),MTRY
	        ENDDO
	        CLOSE (UNIT=21)
	        IF (MCTARGET.GE.(MCTMX/9)) THEN
	          MCTARGET = MCTARGET+MCTMX/9
	        ELSE
	          MCTARGET = 3*MCTARGET
	        ENDIF
	        IF (MCTARGET.GT.MCTMX) STOP
	      ENDIF
C
C	move to the new position
	      DO I=1,NP
	        PW(I,KW) = P(I)
	      ENDDO
	      CH2(KW) = CHSQ
	      NTRY(KW) = 1 
	    ELSE
C	reject: just increment the weight
	      NTRY(KW) = NTRY(KW)+1
	    ENDIF
	  ENDDO
	GO TO 200
	END

	REAL*4 FUNCTION ERROR(P,QUIT)
C
	PARAMETER (NP=10)
	REAL*4 P(NP),QUIT
C       P(1..2) = RA,Dec in radians of pole
C	P(3) = ln(p_V)
C	P(4) = ln(Period) in hours
C	P(5) = ln(J) in MKS
C	P(6) = ln(dia) in km
C	P(7) = ln(fc/(1-fc)), with fc = the crater coverage
C	P(8) = ln(p_IR/p_V), the IR to visible albedo ratio
C	P(9) = ln((B/A)**4/(1-(B/A)**4))
C	P(10) = ln((C/B)**4/(1-(C/B)**4))
C
C	if ERROR>QUIT then return without adding more terms
C
	PARAMETER (NDMX=99)
	PARAMETER (NBAND=4,NL=2*NBAND)
	REAL*4 PHT(NDMX),DET(NDMX),RAT(NDMX),DECT(NDMX),DST(NDMX)
	REAL*4 SIGT(NBAND,NDMX),MAGT(NBAND,NDMX),AMPT(NBAND,NDMX)
	REAL*8 T0(NDMX)
	COMMON /NRMLPTS/ T0,H,DH,MAGT,SIGT,AMPT,DST,DET,PHT,RAT,DECT,NDM
	COMMON /PREDICTED/ WMAGT(NBAND,NDMX),WAMPT(NBAND,NDMX)
	COMMON /BACKCHAN/ DIA,SDIA
	COMMON /ANGLES/ SUNLAT,ELAT,ELHA
	REAL SUNLATD(NDMX),ELATD(NDMX),ELHAD(NDMX)
	COMMON /WATERMARK/ BESTEVER,PB(NP)
C
	REAL*4 R(3),N(3),POLE(3)
	PARAMETER (DTR = 57.2957795, PI = 3.14159265)
	REAL*4 WMAG(NBAND)
	REAL*4 LAM(NL) /3.0974,3.6298,5.0450,4.4130,14.6388,10.0348,
	1	20.4156,23.5869/
	REAL*4 W(NL) /0.51167,0.47952,0.3165,0.6778,0.3775,0.5188,
	1	0.4554,0.5351/
	REAL*4 WVL0(NBAND) /3.3526,4.6028,11.5608,22.0883/
	PARAMETER (NLT=80)
	REAL*4 WVLT(NLT)
	PARAMETER (NALL=NBAND+NL+NLT)
	REAL*4 LN10OVER5
	PARAMETER (LN10OVER5 = 0.460517)
	REAL*4 WVALL(NALL)
	EQUIVALENCE (WVALL,LAM),(WVALL(NL+1),WVL0),
	1	(WVALL(NL+NBAND+1),WVLT)
	REAL*4 F(NALL),CNU(NALL),SNU(NALL)
	REAL*4 FALL(NALL,NDMX)
	REAL*4 F0(NBAND) /306.681,170.663,29.0448,8.2839/
	REAL*4 SIGMIN(NBAND) /4*0.10/
	REAL*4 MICRONS,NUFNU
	REAL*4 INERTIA,JLO,JHI
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
	DO I=1,NLT
	  WVLT(I) = 0.55*(28./0.55)**((I-1)/(NLT-1.))
	ENDDO
C
C	WRITE (*,FMT='(10F8.4)') P
	CALL RADEC_TO_C(P(1),P(2),POLE)
C
	OMEGA = 2*PI/(3600*EXP(P(4)))
	INERTIA = EXP(P(5))
	EM = 0.9
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
C
C       albedo prior
C	d(p_V)/d(ln(P_v)) = p_V^2*[f*/pk1^2*exp(-0.5*(p_V/pk1)^2+(1-f)/pk2^2*exp(-0.5*(p_V/pk2)**2)]
	ERROR = -2*ALOG(PV**2*(FD*EXP(-0.5*(PV/PK1)**2)/PK1**2
	1	+(1-FD)*EXP(-0.5*(PV/PK2)**2)/PK2**2))
C	diameter prior just to stay in bounds
	DIA = EXP(P(6))
	IF (P(6).GT.ALOG(1000.)) ERROR=ERROR+10.*(P(6)-ALOG(1000.))**2
	IF (P(6).LT.ALOG(.001)) ERROR=ERROR+10.*(P(6)-ALOG(.001))**2
C	period prior
C	boundary penalties
	IF (P(4).GT.10.) ERROR=ERROR+100*(P(4)-10)**2
	IF (P(4).LT.-10.) ERROR = ERROR+100*(P(4)+10)**2
C	tipspeed prior
	TIPSPEED = OMEGA*DIA*500.
	TPSP0 = 0.065
	IF (DIA.GT.0.2) TPSP0 = TPSP0*(DIA/0.2)**(5./6.)
	ERROR = ERROR+2*ALOG(1+ALOG(TIPSPEED/TPSP0)**2)
C	spin barrier: P < 2 hrs for D>200 m strongly discouraged
	IF (DIA.GT.0.2 .AND. OMEGA.GT.0.00087266)
	1	ERROR=ERROR+100*ALOG(OMEGA/0.00087266)**2
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
C	prior should be uniform in 0..1 for FC or p(x)dx=FC*(1-FC) with FC=e^x/(1+e^x)
C	X = ln(1/(1-fc)-1)
	DO J=7,10
	  IF (J.NE.8) THEN
	    IF (ABS(P(J)).GT.9.2) THEN
	      ERROR = ERROR+2*(ABS(P(J))+EXP(-ABS(P(J)))-ALOG(4.))
	    ELSE
	      FC = 1/(1+EXP(-P(J)))
C	1-1/(1+exp(-x)) = exp(-x)/(1+exp(-x)), so ln[]=-x-exp(-x)+ln4 for x<<0
C	ln[]=ln[4/(1+exp(-x))(1+exp(x))]
	      ERROR = ERROR-2*ALOG(4*FC*(1-FC))
	    ENDIF
	  ENDIF
	ENDDO
C	prior should be ln(PIR2PV)=0.563+/-0.340, see apj498_786_t1_mrt.f & apj408731t1_mrt.txt
	PIR2PV = EXP(P(8))
	ERROR = ERROR+((P(8)-0.563)/0.340)**2
C
	IF (ERROR.GT.QUIT) RETURN
C
C
	ALB = 0.38399*PV
C	shape
	B2A = (1/(1+EXP(-P(9))))**0.25
	C2B = (1/(1+EXP(-P(10))))**0.25
C	write (*,fmt='(3f8.4,f10.4,1pe10.3)') pole,alb,theta1
C	WRITE (*,FMT='(10F8.4)') P
C	WRITE (*,*) B2A,C2B
	DO J=1,NDM
	  CALL EPHEM(-1,T0(J),R)
	  CALL RADEC_TO_C(RAT(J)/DTR,DECT(J)/DTR,N)
	  THETA = THETA1*DST(J)**1.5
	  DO I=1,3
	    R(I) = -R(I)
	  ENDDO
C	  WRITE (*,*) 'starting predict'
	  CALL PREDICT(R,N,DET(J),POLE,THETA,ALB,EM,FC,PIR2PV,
	1	B2A,C2B,NALL,WVALL,F,CNU,SNU)
C	  WRITE (*,*) 'thru predict'
C	  WRITE (*,FMT='(12F6.2)') (WVALL(I),I=1,12)
C	  WRITE (*,FMT='(1P8E10.3)') (F(I),I=1,8)
C	save the angles
	  SUNLATD(J) = DTR*SUNLAT
	  ELATD(J) = DTR*ELAT
	  ELHAD(J) = DTR*ELHA
C
	  DO I=1,NALL
	    FALL(I,J) = F(I)
	  ENDDO
C	  WRITE (*,FMT='(1P7E11.3)') F
	  DO I=1,NBAND
	    BANDFLUX = W(2*I-1)*F(2*I-1)+W(2*I)*F(2*I)
	    BANDCOS = W(2*I-1)*CNU(2*I-1)+W(2*I)*CNU(2*I)
	    BANDSIN = W(2*I-1)*SNU(2*I-1)+W(2*I)*SNU(2*I)
	    WAMPT(I,J) = 2*1.086*SQRT(BANDCOS**2+BANDSIN**2)/BANDFLUX
	    WMAGT(I,J) = -2.5*LOG10(BANDFLUX/F0(I))
	  ENDDO
C	  WRITE (*,FMT='(4F6.2,2X,4F6.2)') (WMAGT(I,J),I=1,NBAND),
C	1	(WAMPT(I,J),I=1,NBAND)
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
C	amplitude errors
	      IF (AMPT(I,J).GT.0.) THEN
	        SIGAMP = 2.8*SIGT(I,J)
	        PRED = SQRT(SIGAMP**2+WAMPT(I,J)**2)
	        E = ABS(PRED-AMPT(I,J))/SIGAMP
	        EE = E*E
	        IF (E.GT.2.) EE = 4+2*(E-2)
	        RCS = RCS+EE
	      ENDIF
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
	  WRITE (22,FMT='(1H%,F11.6,9F12.6,2F12.4)') PB,ERROR,RCS
	  B2A = (1/(1+EXP(-PB(9))))**0.25
	  C2B = (1/(1+EXP(-PB(10))))**0.25
	  TWO_A = DIA/(C2B*B2A**2)**0.3333333
	  TWO_B = TWO_A*B2A
	  TWO_C = TWO_B*C2B
	  WRITE (22,FMT='(9H% 2(abc)=,F8.3,1Hx,F8.3,1Hx,F8.3,3H km)') 
	1	TWO_A,TWO_B,TWO_C
	  WRITE (22,FMT='(1H%,F10.5,1P2E11.4)') PV,DIA,THETA1
	  WRITE (22,FMT='(A)') '% DE,DS,SunLat,ELat,ELHA,Mags,Amps'
	  DO J=1,NDM
	    WRITE (22,FMT='(1H%,2F7.4,2F6.1,F7.1,4F6.2,1X,4F5.2)')
	1	DET(J),DST(J),SUNLATD(J),ELATD(J),ELHAD(J),
	2	(WMAGT(I,J),I=1,NBAND),(WAMPT(I,J),I=1,NBAND)
	  ENDDO
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

	SUBROUTINE PREDICT(R,N,DE,P,THRM,ALB,EM,
	1	CRATER,PIR2PV,B2A,C2B,NL,LAM,FNU,CNU,SNU)
C
C	computes the fluxes at NL wavelengths in LAM(1..NL)
C	in units of microns.  Fluxes in FNU(1..NL) in Jy
C	CNU & SNU are the cosine and sine amplitudes in 
C		F(t) = F+C*cos(4\pi*t/P)+S*sin(4\pi*t/P)
C	R(1..3) is the position of the observer relative to the Sun in AU
C	N(1..3) is a unit vector toward the object from the Earth
C	DE is the distance of the object from the Earth in AU
C	P is a unit vector giving the rotation pole of the object
C	THRM is the dimensionless thermal inertia parameter
C	ALB is the bolometric Bond albedo
C	EM is the emissivity
C	CRATER is coverage fraction of craters
C	PIR2PV is p_IR/p_V
C	B2A is B/A, C2B is C/B for assumed ellipsoidal shape
C	all vectors have to be in the same coordinate system
C
C	equivalent volume diameter (C*A*A)^(1/3) is fixed at 1 km
C
	REAL R(3),N(3),P(3),THRM,ASYM,ROT
	REAL DE		! distance to Earth in AU
	INTEGER NL
	REAL*4 LAM(NL),FNU(NL),CNU(NL),SNU(NL)
C
	REAL OBJ(3),Z(3),X(3),Y(3)
	INTEGER NFACE,NLAT,NTIME,NLIST,NPXL
	INTEGER I, K
	PARAMETER (NFACE=127,NLAT=16,NTIME=32)
C
	PARAMETER (NPXL=12*4**6)
	REAL AXIS(3,NPXL)
	REAL PFACE(3,NPXL)
	INTEGER*2 FCLST(NPXL)
	COMMON /SPHERE_TABLE/ AXIS,PFACE,FCLST
C
	REAL INSOL(NTIME,0:NFACE,NLAT)		! solar input
	REAL T(NTIME,0:NFACE,NLAT)	! the temperatures
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
	COMMON /ANGLES/ SUNLAT,ELAT,ELHA
	REAL G(NTIME)	! normalized thermal response
C	axes of ellipsoid A >= B >= C
	REAL A,B,C,WATES(NLAT,NTIME)
	INTEGER IWATES(NLAT,NTIME)
	REAL TBRITE(12,2,NTIME),AREA(NTIME)
	REAL TMP	! temporary variable used to rotate vectors
	REAL THCC /1.19104272E-5/
	REAL OMEGA
	REAL HCK /1.4387752/	! h*c/k
	REAL XB		! h*nu/k*Tb
	PARAMETER (TSUN=5933.)
	PARAMETER (SOLCON=1367.E3)	! solar constant in erg/cm^2/sec
	REAL SIGSB /5.6704E-5/	! Stefan-Boltzmann constant
	REAL TAT1AU /394.03854311248/
	REAL ALBEDO, EMISSIVITY
	REAL ETA,DIA
	REAL T0
	REAL DS		! distance to Sun in AU
	REAL AU /1.495979E13 /	! AU in cm
	REAL LXT(9) /-4.,-3.,-2.,-1.,0.,1.,2.,3.,4./
	INTEGER NLMX
	PARAMETER (NLMX=199)
	REAL FHG(NLMX)	! reflected light flux for 1 km sphere
	REAL TBT2(9,2,NTIME),LX(NLMX)  ! variables for interpolating TB
	REAL FSUM(NLMX,NTIME),F,FLAM(NLMX,NTIME)
	REAL VSUM,FOPT
	INTEGER IDS,IA,IDE,ILAM
	REAL JY,MAB
	REAL SPLINE
	REAL PI,TWOPI,DTR
	PARAMETER (PI=3.1415927,TWOPI=6.2831853,DTR=57.29578)
	PARAMETER (Q=0.38399,GG=0.25,A1=3.33,B1=0.63,A2=1.87,B2=1.22)
C
	INTEGER IMC,NEE
	REAL RMS
	REAL*8 SE,SSEE,SEE
	INTEGER IY,IZ,IR,JR
	REAL RAN
	LOGICAL TBNEG,TBNAN
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
	TA2 = TAN(0.5*PHASE/DTR)
	PHI1 = EXP(-A1*TA2**B1)
	PHI2 = EXP(-A2*TA2**B2)
	PHI = (1-GG)*PHI1+GG*PHI2
C	WRITE (*,FMT='(3F9.4,3X,3F9.5)') OBJ,Z
C	WRITE (*,FMT='(F9.4,3F9.3,F9.1)') DS,SUNLAT,ELAT,ELHA,PHASE
C------------
C
	THMAX = 75
	THMAX=THMAX/DTR
	ALBEDO = ALB
	EMISSIVITY = EM
	DO I=1,NLAT
	  SL = -1.+(2.*I-1.)/NLAT
	  CL = SQRT(1.-SL**2)
	  LATS(I) = ATAN(SL/CL)
	  WLAT(I) = 2./NLAT
	ENDDO
C
C	was:
C	DO K=1,NLAT
C	compute the solar input
C	  CALL SOLAR(SUNLAT,LATS(K),NFACE,NTIME,THMAX,PFACE,AXIS,
C	1	FCLST,ALBEDO,INSOL(1,1,K))
C	  CALL THERMAL(NFACE,NTIME,THRM,G,THMAX,INSOL(1,1,K),T(1,1,K))
C	ENDDO
C
C	now:
	CALL INTERP_TABLE(SUNLAT,THRM,INSOL,T)
	DO K=1,NLAT
	  DO J=0,NFACE
	    DO I=1,NTIME
	      IF (INSOL(I,J,K).NE.INSOL(I,J,K)
	1	.OR. T(I,J,K).NE.T(I,J,K)) THEN
	        WRITE (*,*) 'NaN in INSOL or T'
	        STOP
	      ENDIF
	    ENDDO
	  ENDDO
	ENDDO
C
	CSL = COS(SUNLAT)
	SSL = SIN(SUNLAT)
C	solid angle of 1 km diameter at DE AU
	OMEGA = (PI/4)*(1.E5/(DE*AU))**2
	ALPHA = PHASE/DTR
C	compute relected solar flux for 1 km object
	DO ILAM=1,NL
	  XS = HCK*10000./(TSUN*LAM(ILAM))
C       the solar flux, nu*F_\nu, in erg/cm^2/sec/oct.  1.056 is a fudge factor
C	to make the V and W1 fluxes correct
	  FSUN = 1.056*SOLCON*(XS**4/(EXP(XS)-1))/(DS**2*(PI**4)/15)
C	  IF (AL.LT.0.4) FSUN = 0.6*FSUN ! to make bolometric flux right
C	face-on white disk flux at 1 AU for 1 km is
C	(PI*R**2)*FSUN/(PI*AU**2) = FSUN*(R/AU)**2
	  PV = ALBEDO/Q
	  PIR=PV*EXP(MAX(0.,MIN(1.,LOG(LAM(ILAM)/0.55)/LOG(3.4/0.55)))*
	1	LOG(PIR2PV))
	  FHG(ILAM) = PHI*PIR*FSUN*(5.E4/(DE*AU))**2
	ENDDO
C	set diameters to give volume of 1 km diameter sphere
	A = 1/(C2B*B2A**2)**(1./3)
	B = A*B2A
	C = B*C2B
C	A*B*C = B2A*B2A*C2B/(C2B*B2A**2) = 1
	CALL ELLIPSOID_WEIGHTS(A,B,C,NLAT,NTIME,WATES)
C
	DO I=1,NLAT
	  DO J=1,NTIME
	    IWATES(I,J) = 100*WATES(I,J)
	  ENDDO
	ENDDO
C	WRITE (*,FMT='(16I5)') IWATES
C
	T0 = TAT1AU*((1-ALBEDO*(1+0.19*(PIR2PV-1)))/EMISSIVITY/DS**2)**0.25
	DO ILAM=1,NL
	  LX(ILAM) = ALOG(HCK*10000.0/(LAM(ILAM)*T0))/ALOG(2.)
	ENDDO
C	calculate the fluxes seen at Earth
	CALL FLUXES(ELAT,ELHA,NLAT,LATS,WLAT,NFACE,NTIME,AXIS,
	1	PFACE,FCLST,INSOL,T,WATES,TBRITE,AREA)
C	WRITE (*,FMT='(8H   Area:,8F9.4/(8X,8F9.4))') AREA
C	WRITE (*,FMT='(3HTB:,9F8.4)') (TBRITE(I,1,1),I=1,9)
C	WRITE (*,FMT='(3HTB:,9F8.4)') (TBRITE(I,1,5),I=1,9)
C	WRITE (*,FMT='(3HTB:,9F8.4)') (TBRITE(I,1,9),I=1,9)
C	WRITE (*,FMT='(3HTB:,9F8.4)') (TBRITE(I,1,13),I=1,9)
	TBNEG = .FALSE.
	TBNAN = .FALSE.
	DO J=1,NTIME
	 DO I=1,12
	  DO K=1,2
	    IF (TBRITE(I,K,J).NE.TBRITE(I,K,J)) TBNAN = .TRUE.
	    IF (TBRITE(I,K,J).LT.0.) TBNEG = .TRUE.
	  ENDDO
	 ENDDO
	ENDDO
	IF (TBNAN) THEN
	  WRITE (*,*) 'NaN in TB table'
	  STOP
	ENDIF
	IF (TBNEG) THEN
	  WRITE (*,*) 'TBRITE < 0'
	  WRITE (*,*) ELAT,ELHA,SUNLAT,THRM
	  WRITE (*,*) TBRITE
	  DO ILAM = 1,NL
	    FNU(ILAM) = 0
	  ENDDO
	  RETURN
	ENDIF
	DO K=1,NTIME
	  DO J=1,2
	    CALL SPLSET(9,LXT,TBRITE(1,J,K),TBT2(1,J,K))
	  ENDDO
	ENDDO
C	WRITE (*,*) 'thru SPLSET'
	TBMIN = 1.E4
	TBMAX = 0.
	DO K=1,NTIME
	 DO ILAM = 1,NL
	  TBC = T0*SPLINE(LX(ILAM),9,LXT,TBRITE(1,1,K),TBT2(1,1,K))
	  XB = HCK*10000./(LAM(ILAM)*TBC)
	  F1 = 0.
	  IF (XB.LT.70.) F1 = 1./(EXP(XB)-1.)
	  TBF = T0*SPLINE(LX(ILAM),9,LXT,TBRITE(1,2,K),TBT2(1,2,K))
	  TBMIN = MIN(TBMIN,TBC,TBF)
	  TBMAX = MAX(TBMAX,TBC,TBF)
	  XB = HCK*10000./(LAM(ILAM)*TBF)
	  F2 = 0.
	  IF (XB.LT.70.) F2 = 1./(EXP(XB)-1.)
	  F = CRATER*F1+(1-CRATER)*F2
	  FSUM(ILAM,K) = THCC*OMEGA*AREA(K)*
	1	(10000./LAM(ILAM))**4*EMISSIVITY*F
C	write (*,*) lam(ilam),fhg,phi,pir,fsun
C	stop
	  FSUM(ILAM,K) = FSUM(ILAM,K)+FHG(ILAM)*AREA(K)
C	convert to JY from nu*F_nu in erg/cm^2/sec
	  JY = FSUM(ILAM,K)*LAM(ILAM)/3.E-9
	  FLAM(ILAM,K) = JY
	 ENDDO
	ENDDO
C------------
C	summarize the lightcurves
C
C	WRITE (*,*) 'min,max TB:',TBMIN,TBMAX
C	WRITE (*,*) 'starting summary'
	DO I=1,NL
	  FNU(I) = 0
	  CNU(I) = 0
	  SNU(I) = 0
C	sum over lightcure
	  DO K=1,NTIME
	    FNU(I) = FNU(I)+FLAM(I,K)
	    CSN = COS((K-1)*PI*4./NTIME)
	    SSN = SIN((K-1)*PI*4./NTIME)
	    CNU(I) = CNU(I)+CSN*FLAM(I,K)
	    SNU(I) = SNU(I)+SSN*FLAM(I,K)
	  ENDDO
C	normalize
	  FNU(I) = FNU(I)/NTIME
	  CNU(I) = 2*CNU(I)/NTIME
	  SNU(I) = 2*SNU(I)/NTIME
	ENDDO
C
	RETURN
	END

	SUBROUTINE ELLIPSOID_WEIGHTS(A,B,C,NLAT,NTIME,WATES)
C
	REAL A,B,C
	INTEGER NLAT,NTIME
	REAL WATES(NLAT,NTIME)
	REAL MU,MUMIN,MUMAX
C
	REAL R(3)
	REAL N(3),E(3)
	REAL NORMAL(3)
	INTEGER I,J, LAT,LONG
	PARAMETER (IRES=7,NPIX=12*4**IRES)
	PARAMETER (PI = 3.1415926535898)
	PARAMETER (TWOPI = 2*PI)
	PARAMETER (DTR = 180/PI)
	DO I=1,NLAT
	  DO J=1,NTIME
	    WATES(I,J) = 0
	  ENDDO
	ENDDO
	W = FLOAT(NLAT)*FLOAT(NTIME)/NPIX
	MUMIN = 99.
	MUMAX = -99.
	DO I=0,NPIX-1
	  CALL NESTED_CENPIX(I,IRES,R)
	  RR = 1/((R(1)/A)**2+(R(2)/B)**2+(R(3)/C)**2)
	  E(1) = -R(2)*(A/B)**2
	  E(2) = R(1)
	  E(3) = 0
	  N(1) = -R(3)*R(1)
	  N(2) = -R(3)*R(2)
	  N(3) = C**2*((R(1)/A)**2+(R(2)/B)**2)
	  CALL CROSS(E,N,NORMAL)
	  CALL NORM(NORMAL)
	  MU = DOT(R,NORMAL)
	  MUMIN = MIN(MU,MUMIN)
	  MUMAX = MAX(MU,MUMAX)
	  DA = RR/DOT(R,NORMAL)
	  CALL C_TO_RADEC(NORMAL,RA,DEC)
	  LAT = NLAT*0.5*(NORMAL(3)+1.)+1.
	  LAT = MIN(NLAT,MAX(1,LAT))
	  LONG = 1+NTIME*RA/TWOPI
	  LONG = MIN(NTIME,MAX(1,LONG))
	  WATES(LAT,LONG) = WATES(LAT,LONG)+W*DA
	ENDDO
C	WRITE (*,*) MUMIN,MUMAX
	RETURN
	END

	SUBROUTINE SPHERE(NFACE,THMAX,AXIS,PFACE,FCLST)
C
C	sets up AXIS, PFACE & FCLST for a spherical cap crater
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
C
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
	PARAMETER (NPXL=12*4**6, NFMAX=127)
	REAL AXIS(3,NPXL)
	REAL PFACE(3,NPXL)
	INTEGER*2 FCLST(NPXL)
	REAL WATE(NFACE)
	INTEGER WF(NFMAX)
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
C	use average of cos factor over the healpix pixels in this facet
	  IF (WF(I).GE.2) WATE(I) = WATE(I)/WF(I)
	ENDDO
	GO TO 40
20	WRITE (*,*) 'negative weights in LWATE'
	WRITE (*,FMT='(10F8.5)') WATE
40	SUM = 0
	DO I=1,NFACE
	  SUM = SUM+WATE(I)
	ENDDO
	SUM1 = SUM
	IF (SUM.LE.0.) THEN
C	  WRITE (*,*) ' SUM LE 0 in LWATE'
	  RETURN
	ENDIF
	SUM = E(3)/SUM
	DO I=1,NFACE
	  WATE(I) = SUM*WATE(I)
	  IF (WATE(I).NE.WATE(I)) THEN
	    WRITE (*,*) ELAT,ELHA,LAT,NFACE
	    WRITE (*,*) SUM,SUM1,E
	    STOP 'NaN in LWATE'
	  ENDIF
	ENDDO
	RETURN
	END

	SUBROUTINE FLUXES(ELAT,ELHA,NLAT,LATS,WLAT,NFACE,NTIME,AXIS,
	1	PFACE,FCLST,INSOL,T,WATES,TBRITE,AREA)
C	routine to compute the fluxes at Earth viewing point
C	ELAT and ELHA are the sub-Earth latitude and local hour angle
C	in radians.  +ELHA means afternoon, -ELHA is morning
C	NLAT is the number of latitudes, listed in LATS, with weights WLAT
C	weight should be COS(LAT)*DLAT
C	output is array TBRITE(I,*,*)
C	I=1-9 is h*nu/k*To = 1/16,1/8,...8,16
C	10 is the bolometric or effective temperature
C	11 is the radio temperature for deep penetration and low dielectric
C	constant
C	12 is the result for optical light scaled so 
C	ALBEDO/(1-ALBEDO)*TBRITE(12) gives 1 for a white diffuse
C	TBRITE (*,1,*) is for the cratered model, (*,2,*) is for a flat model
C	TBRITE (*,*,K) is for the K'th time in the lightcurve
C
	PARAMETER (NFMAX=127,NTMAX=32)
	REAL TBAR(0:NFMAX),WATE(NFMAX)
	PARAMETER (NPXL=12*4**6)
	REAL AXIS(3,NPXL)
	REAL PFACE(3,NPXL)
	INTEGER*2 FCLST(NPXL)
	REAL T(NTIME,0:NFACE,NLAT)	! the temperatures
	REAL TBRITE(12,2,NTIME)
	REAL V(12)
	REAL WSUM(NTIME)
	REAL AREA(NTIME)
	REAL INSOL(NTIME,0:NFACE,NLAT)         ! solar input
	REAL LATS(NLAT)	! the latitudes
	REAL WLAT(NLAT)	! weight for latitude integration
	REAL WATES(NLAT,NTIME) ! area correction for aspherical shape
	REAL*8 X
	REAL PI,TWOPI,DTR
	PARAMETER (PI=3.1415927,TWOPI=6.2831853,DTR=57.29578)
C
	DO K=1,NTIME
	  WSUM(K) = 0.
	  DO J=1,2	
	    DO I=1,12
	      TBRITE(I,J,K) = 0.
	    ENDDO
	  ENDDO
	ENDDO
	saw = 0
	DO K=1,NLAT
	  DO J=0,NFACE
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
	      saw = saw+abs(wate(j))
	    ENDDO
	    IF (SW.GT.0.00001) THEN
	      DO J=1,NFACE
	        W = WLAT(K)*WATE(J)
	        IF (T(I,J,K).GT.0.) THEN
	          V(10) = W*T(I,J,K)**4
	          V(11) = W*TBAR(J)
	          V(12) = W*INSOL(I,J,K)
	          X = 0.0625/T(I,J,K)
	          DO L=1,9
	            IF (X.LE.70.) V(L) = W/(EXP(X)-1.)
	            X = X+X
	          ENDDO
	          DO IT = 1,NTIME
	            JT = I+IT-1
	            IF (JT.GT.NTIME) JT = JT-NTIME
	            DO L=1,12
	              TBRITE(L,1,IT) = TBRITE(L,1,IT)+V(L)*WATES(K,JT)
	            ENDDO
	          ENDDO
	        ENDIF
	        DO IT = 1,NTIME
	          JT = I+IT-1
	          IF (JT.GT.NTIME) JT = JT-NTIME
	          WSUM(IT) = WSUM(IT)+W*WATES(K,JT)
	        ENDDO
	      ENDDO	! J=1,NFACE
C	repeat for the flat model in FACE=0
	      W = WLAT(K)*SW
	      IF (T(I,0,K).GT.0.) THEN
	        V(10) = W*T(I,0,K)**4
	        V(11) = W*TBAR(0)
	        V(12) = W*INSOL(I,0,K)
	        X = 0.0625/T(I,0,K)
	        DO L=1,9
	          IF (X.LE.70.) V(L) = W/(EXP(X)-1.)
	          X = X+X
	        ENDDO
	        DO IT = 1,NTIME
	          JT = I+IT-1
	          IF (JT.GT.NTIME) JT = JT-NTIME
	          DO L=1,12
	            TBRITE(L,2,IT) = TBRITE(L,2,IT)+V(L)*WATES(K,JT)
	          ENDDO
	        ENDDO
	      ENDIF
	    ENDIF	! if SW>0
	  ENDDO	! I=1,NTIME
	ENDDO	! K=1,NLAT
C	normalize
	DO IT=1,NTIME
	  AREA(IT) = WSUM(IT)*(2./NTIME)
	  DO J=1,2
	    DO I=1,12
	      TBRITE(I,J,IT) = TBRITE(I,J,IT)/WSUM(IT)
	    ENDDO
C	convert from flux to brightness T
	    TBRITE(10,J,IT) = TBRITE(10,J,IT)**0.25
	    X = 0.0625
	    DO I=1,9
	      IF (TBRITE(I,J,IT).GT.0.) 
	1	TBRITE(I,J,IT) = X/ALOG(1.+1./TBRITE(I,J,IT))
	      IF (I.GE.3 .AND. TBRITE(I,J,IT).LE.0.) 
	1	TBRITE(I,J,IT) = 2*TBRITE(I-1,J,IT)-TBRITE(I-2,J,IT)
	      X = X+X
	    ENDDO
	  ENDDO
	ENDDO
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

	SUBROUTINE INTERP_TABLE(SUNLAT,THRM,INSOL,T)
C
C	interpolates in the tables read in by READ_TABLE
C
	REAL SUNLAT	! input subsolar latitude, in radians
	REAL THRM	! input dimensionless thermal inertia parameter Theta
	PARAMETER (NFACE=127,NLAT=16,NTIME=32,NSL=31,NTHRM=46)
	REAL INSOL(NTIME,0:NFACE,NLAT), T(NTIME,0:NFACE,NLAT)
C
	REAL INSOLT(NTIME,0:NFACE,NLAT,NSL)		! solar input
	REAL TT(NTIME,0:NFACE,NLAT,NTHRM,NSL)	! the temperatures
C
	COMMON /TABLES/ INSOLT,TT
C
	PARAMETER (DTR= 57.2957795)
	FX = 1+(NTHRM-1)*THRM/(THRM+1.6)
	IX = FX
	IF (IX.LE.0) IX=1
	IF (IX.GE.NTHRM) IX = NTHRM-1
	FX = FX-IX
	FSL = (SUNLAT*DTR+96)/6
	ISL = FSL
	IF (ISL.LE.0) ISL=1
	IF (ISL.GE.NSL) ISL=NSL-1
	FSL = FSL-ISL
C
	DO K=1,NLAT
	  DO J=0,NFACE
	    DO I=1,NTIME
	      INSOL(I,J,K) = (1-FSL)*INSOLT(I,J,K,ISL)+
	1	FSL*INSOLT(I,J,K,ISL+1)
	      T(I,J,K) = 
	1	(1-FSL)*((1-FX)*TT(I,J,K,IX,ISL)+FX*TT(I,J,K,IX+1,ISL))
	2	+FSL*((1-FX)*TT(I,J,K,IX,ISL+1)+FX*TT(I,J,K,IX+1,ISL+1))
	    ENDDO
	  ENDDO
	ENDDO
	RETURN
	END

