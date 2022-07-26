C	WISE-NEATM-MCMC.f
C	gfortran -o WISE-NEATM-MCMC WISE-NEATM-MCMC.f ~/Ned.a
C
C	created 21-DEC-2017
C
	PARAMETER (NP=4,NDMX=99)
	REAL*4 PHT(NDMX),DET(NDMX),RAT(NDMX),DECT(NDMX),DST(NDMX)
	REAL*4 SIGT(4,NDMX),MAGT(4,NDMX)
	REAL*8 T0(NDMX)
	COMMON /NRMLPTS/ T0,H,DH,PER,MAGT,SIGT,DST,DET,PHT,RAT,DECT,NDM
	COMMON /PREDICTED/ WMAGT(4,NDMX)
	REAL*8 MJD,JD(NDMX)
	REAL*8 T1JAN2010,T,JD1JAN2010
	REAL*8 TSTART,TWALL
	REAL*8 JULIAN_DAY,GMT2T68
	REAL*4 EA(3),A(3),E(3)
	REAL*4 WMAG(4),SIG(4)
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
	REAL*4 SIGMIN(4) /4*0.02/
	EXTERNAL RAN
C
	PARAMETER (MCTMX=72900)
	PARAMETER (DTR= 57.2957795)
	REAL MCT(NP+2,MCTMX)
C
	CALL RANDOMIZE
C
	BESTEVER = 1.E30
C
	NDM=0
C	period is just ignored.  Included for compatibility
	WRITE (*,FMT='(A,$)') 'Enter H,sig(H),Period,Name:'
	READ (*,FMT='(3F7.3,A)') H,DH,PER,NAME
	WRITE (*,FMT='(3F7.3,1X,A)') H,DH,PER,NAME
	CALL TIME_TAG(2010,1,1,0,0,0.,T1JAN2010)
	JD1JAN2010 = JULIAN_DAY(T1JAN2010)
	WRITE (*,FMT='(A)') 'Enter MJD,RA,Dec,Delta,W1234 with sigmas:'
	WRITE (*,FMT='(A)') '   MJD      RA      DEC     DE      DS'//
	1	'    PHASE W1234 & sigmas'
10	READ (*,900,END=100) MJD,RA,DEC,DELTA,(WMAG(I),SIG(I),I=1,4)
900	FORMAT(F11.5,2F11.7,F8.4,4(F7.3,F6.3))
	IF (NDM.GE.NDMX) STOP '2 much data'
	NDM=NDM+1
	DO I=1,4
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
	1	PHT(NDM),(WMAG(I),SIG(I),I=1,4)
901	FORMAT(F9.3,2F8.3,2F8.5,F7.2,4(F7.3,F6.3))
	GO TO 10
100	CLOSE (UNIT=1)
	WRITE (*,*) NDM
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
C	to use ln(p_V),ln(p_IR/p_V),ln(eta),ln(D) as parameters
C	albedo uniform in log 0.01 to 1.01
	  PW(1,I) = ALOG(RAN(IR,JR)+0.01)
C	periods 0.0001 hours to 1000 hours
C	p_IR/p_V log normal
C	prior should be ln(PIR2PV)=0.288+/-0.340, see apj498_786_t1_mrt.f
C	also apj408731t1_mrt.txt Mainzer etal (2011) gives 0.563 for mean
	  PW(2,I) =   EXP(0.563+0.340*GAUSRN(IR,JR))
	  PW(3,I) = 0.3365+0.4418*GAUSRN(IR,JR)
C	diameter from 1 meter to 1000 km
	  IF (I.EQ.1) THEN
	    PW(4,I) = ALOG(.001)+ALOG(1.E6)*RAN(IR,JR)
	  ELSE
	    PW(4,I) = PW(4,1)+0.5*GAUSRN(IR,JR) 
	  ENDIF
	  DO J=1,NP
	    P(J) = PW(J,I)
	  ENDDO
	  D1 = EXP(P(4))
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
	      SWDD=SWDD+ALOG(1329./10.**(H/5)/SQRT(EXP(PW(1,I))))-PW(4,I)
	    ENDIF
	    DO K=1,NDM
	      DO J=1,4
	        IF (SIGT(J,K).LT.8.) THEN
	          SW=SW+1
C	add 0.4605*(pred-obs) mag so if D is too big, get a negative
	          SWDD = SWDD+0.4605*(WMAGT(J,K)-MAGT(J,K))
	        ENDIF
	      ENDDO
	    ENDDO
	    IF (E1.NE.E1) THEN
	      P(4) = 0.
	    ELSE
	      P(4) = (SWDD/SW)+P(4)
	    ENDIF
	    PW(4,I) = P(4)
	    D2 = EXP(P(4))
	    CH2(I) = ERROR(P,1.E30)
	    IMC = IMC+1
	    E2 = CH2(I)
C	write (*,*) d1,e1,d2,e2
	  ENDIF
C
	  NTRY(I) = 0
	  WRITE (*,FMT='(I3,8F8.3,F15.2)') I,(PW(J,I),J=1,NP),CH2(I)
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
	    DO I=1,NP
	      P(I) = PW(I,JW)+Z*(PW(I,KW)-PW(I,JW))
	    ENDDO
	    TWALL = GMT2T68('NOW')-TSTART
	    CHSQ = ERROR(P,CH2(KW)+100)
	    IMC = IMC+1
	    IF (MOD(IMC,1000).EQ.0) 
	1	WRITE (*,FMT='(I7,I6,F10.0)') IMC,NMCT,TWALL
	    Q = 0
	    XP = -0.5*(CHSQ-CH2(KW))
	    IF (XP.GT.50.) Q = 1
	    IF (ABS(XP).LT.50.) 
	1	Q = MIN(1.0,ABS(Z**(NP-1)*EXP(XP)))
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
	          WRITE (21,FMT='(F8.5,3F9.5,F17.5,I5)')
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
	PARAMETER (NP=4)
	REAL*4 P(NP),QUIT
C	P(1) = ln(p_V)
C	P(2) = ln(p_IR/p_V), the IR to visible albedo ratio
C	P(3) = ln(eta)
C	P(4) = ln(dia) in km
C
	PARAMETER (NDMX=99)
	PARAMETER (NBAND=4,NL=2*NBAND)
	REAL*4 PHT(NDMX),DET(NDMX),RAT(NDMX),DECT(NDMX),DST(NDMX)
	REAL*4 SIGT(4,NDMX),MAGT(4,NDMX)
	REAL*8 T0(NDMX)
	COMMON /NRMLPTS/ T0,H,DH,PER,MAGT,SIGT,DST,DET,PHT,RAT,DECT,NDM
	COMMON /PREDICTED/ WMAGT(4,NDMX)
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
	REAL*4 F(NALL)
	REAL*4 FALL(NALL,NDMX)
	REAL*4 F0(NBAND) /306.681,170.663,29.0448,8.2839/
	REAL*4 SIGMIN(NBAND) /4*0.10/
	REAL*4 MICRONS,NUFNU
	REAL*4 INERTIA,JLO,JHI
C	parameters of albedo prior
	PARAMETER (PK1=0.03,FD=0.251,PK2=0.167)
C
	DO I=1,NLT
	  WVLT(I) = 0.55*(28./0.55)**((I-1)/(NLT-1.))
	ENDDO
C
C	T1 = 401 K
	EM = 0.90
	T1 = (1367./(EM*5.67E-8))**0.25
C	range restrictions to protect the RC code
	TH1MAX = 10000
	THETA1 = TH1MAX*THETA1/(TH1MAX+THETA1)
	PVMAX = 1.5
	IF (P(1).GT.23.) THEN
	  PV = PVMAX
	ELSE
	  PV = EXP(P(1))
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
	DIA = EXP(P(4))
	IF (P(4).GT.ALOG(1000.)) ERROR=ERROR+10.*(P(4)-ALOG(1000.))**2
	IF (P(4).LT.ALOG(.001)) ERROR=ERROR+10.*(P(4)-ALOG(.001))**2
C	prior should be ln(PIR2PV)=0.563+/-0.340, see apj498_786_t1_mrt.f & apj408731t1_mrt.txt
	PIR2PV = EXP(P(2))
	ERROR = ERROR+((P(2)-0.563)/0.340)**2
C	ETA prior, 1.4+/-0.5 for NEOs, but log normal here
	ETA = EXP(P(3))
	ERROR = ERROR+((P(3)-0.3365)/0.4418)**2
C
	IF (ERROR.GT.QUIT) RETURN
C
C
	ALB = 0.38399*PV
C	write (*,fmt='(3f8.4,f10.4,1pe10.3)') pole,alb,theta1
	DO J=1,NDM
	  CALL EPHEM(-1,T0(J),R)
	  CALL RADEC_TO_C(RAT(J)/DTR,DECT(J)/DTR,N)
	  DO I=1,3
	    R(I) = -R(I)
	  ENDDO
	  CALL PREDICT(R,N,DET(J),ALB,EM,PIR2PV,ETA,NALL,WVALL,F)
C
	  DO I=1,NALL
	    FALL(I,J) = F(I)
	  ENDDO
	  DO I=1,NBAND
	    WMAGT(I,J)=-2.5*LOG10((W(2*I-1)*F(2*I-1)+W(2*I)*F(2*I))
	1	/F0(I))
	  ENDDO
	ENDDO
C
	DIA = EXP(P(4))
	DM = -5*P(4)/ALOG(10.)
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
	  WRITE (22,FMT='(1H%,F11.6,7F12.6,2F12.4)') P,ERROR,RCS
	  WRITE (22,FMT='(1H%,F10.5,1P2E11.4)') PV,DIA,ETA
	  WRITE (22,FMT='(A)') '% DE,DS,Phase,Mags'
	  DO J=1,NDM
	    WRITE (22,FMT='(1H%,2F7.4,F7.1,4F6.2)') DET(J),
	1     DST(J),PHT(J),(WMAGT(I,J),I=1,NBAND)
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

	SUBROUTINE PREDICT(R,N,DE,ALB,EM,PIR2PV,ETA,NL,LAM,FNU)
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
C	FC is coverage fraction of crates
C	PIR2PV is p_IR/p_V
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
C
	INTEGER ISL
	REAL ELONGATION,CELONG,SELONG	! elongation & cos & sin
	REAL THMAX	! max angle of crater walls
	REAL ALPHA	! phase angle for EVEC calculations
	REAL PHASE	! phase angle in degrees
	REAL ELHA	! local hour angle of EVEC
	REAL ELAT	! latitude of Earth
	COMMON /ANGLES/ SUNLAT,ELAT,ELHA
	REAL TMP	! temporary variable used to rotate vectors
	REAL THCC /1.19104272E-5/
	REAL OMEGA
	REAL HCK /1.4387752/	! h*c/k
	REAL XB		! h*nu/k*Tb
	REAL LSUN
	PARAMETER (LSUN = 3.826E33)
	PARAMETER (TSUN=5933.)
	PARAMETER (SOLCON=1367.E3)	! solar constant in erg/cm^2/sec
	REAL SIGSB /5.6704E-5/	! Stefan-Boltzmann constant
	REAL TAT1AU /394.03854311248/
	REAL ALBEDO, EMISSIVITY
	REAL ETA,DIA
	REAL T0
	REAL DS		! distance to Sun in AU
	REAL AU /1.495979E13 /	! AU in cm
	INTEGER NLMX
	PARAMETER (NLMX=199)
	REAL FSUM(NLMX),F,FLAM(NLMX)
	REAL VSUM,FOPT
	INTEGER IDS,IA,IDE,ILAM
	REAL JY,MAB
	REAL SPLINE
	REAL PI,TWOPI,DTR
	PARAMETER (PI=3.1415927,TWOPI=6.2831853,DTR=57.29578)
	PARAMETER (Q=0.38399,GG=0.15,A1=3.33,B1=0.63,A2=1.87,B2=1.22)
C
	INTEGER IMC,NEE
	REAL RMS
	REAL*8 SE,SSEE,SEE
	INTEGER IY,IZ,IR,JR
	REAL RAN
	LOGICAL TBNEG
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
	ALBEDO = ALB
	PV = ALB/Q
	EMISSIVITY = EM
	OMEGA = (PI/4)*(1.E5/(DE*AU))**2
	ALPHA = PHASE/DTR
	T0=TAT1AU*((1-ALBEDO*(1+0.19*(PIR2PV-1)))/(EM*ETA*DS**2))**0.25
	DO ILAM = 1,NL
	  AL = LAM(ILAM)
	  AN = 10000./AL
C	x = h*nu/kT
	  XS = AN*HCK/TSUN
C	the solar flux, nu*F_\nu, in erg/cm^2/sec/oct, 1.056 is a fudge to make V and W1 right
	  FSUN = 1.056*(XS**4/(EXP(XS)-1))/
	1	((PI**4)/15)*LSUN/(DS*AU)**2/(4*PI)
C	  IF (AL.LT.0.4) FSUN = 0.6*FSUN ! to make bolometric flux right
	  XIR = AN*HCK/T0
C	compute average value of 1/(e^x-1) over disk for thermal emission
	  THERMAL = (4*PHASUM(XIR,PH,40)-PHASUM(XIR,PH,20))/3
C	flux is average surface brightness * solid angle, [erg/cm^2/sec/oct]
C	1.03 is a fudge factor
C	THCC has units erg-cm^2/sec, and times AN^4 gives erg/cm^2/sec
	  PIR=PV*EXP(MAX(0.,MIN(1.,LOG(LAM(ILAM)/0.55)/LOG(3.4/0.55)))*
	1	LOG(PIR2PV))
	  FLUX = (THERMAL*EM*THCC*AN**4 + 1.03*PIR*PHI*FSUN/PI)
	1	*OMEGA
	  FJY = 1.E23*FLUX/(3.E14/AL)
	  FNU(ILAM) = FJY
	ENDDO
C------------
C
	RETURN
	END

	REAL FUNCTION PHASUM(HNUKTSS,PHASE,NPOINT)
C
C	compute average of 1/(exp[h*nu/kT]-1) over the projected area of a
C	sphere with Sun at phase angle PHASE.
C	if HNUKTSS < 0 then do average of (T/TSS)^4 instead for bolometric
C
	IMPLICIT NONE
	REAL HNUKTSS	! h*nu/k*T at the sub solar point
	REAL PHASE	! phase angle in degrees
	INTEGER NPOINT	! # points to use in average in each of two directions
C
	REAL*8 SUM
	REAL XX
	INTEGER IMU
	REAL DMU
	REAL MU,ST,CPMN,PMX,PMN
	INTEGER IP
	REAL DP,P,CP,SP
	REAL Z
	REAL V
	REAL CA,SA
	REAL PI,DTR
	PARAMETER (PI = 3.1415926535898, DTR = 180/PI)
C
	CA = COS(PHASE/DTR)
	SA = SIN(PHASE/DTR)
	SUM = 0
	DMU = 1./NPOINT
	DO IMU = 1,NPOINT
	  MU = DMU*(IMU-0.5)
	  V = MU
	  IF (HNUKTSS.GT.0) THEN
	    XX = HNUKTSS/MU**0.25
	    V = 0.
	    IF (XX.LT.70.) V =1./(EXP(XX)-1.)
	  ENDIF
	  ST = SQRT(1.-MU**2)
	  IF (SA.NE.0.) THEN
	    CPMN = MU*CA/(ST*SA)
	  ELSE
	    CPMN = -1.
	    IF (MU.GT.0.) CPMN = 1.
	  ENDIF
	  PMX = PI
	  IF (CPMN.GT.-1.) THEN
	    IF (CPMN.LT.1.) THEN
	      PMN = ACOS(CPMN)
	    ELSE
	      PMN = 0.
	    ENDIF
	    DP = (PMX-PMN)/NPOINT
	    DO IP=1,NPOINT
	      P = PMN+DP*(IP-0.5)
	      SP = SIN(P)
	      CP = COS(P)
	      Z = MU*CA-ST*CP*SA
	      IF (Z.LT.0.) WRITE (*,*) Z,P,MU,PHASE
	      SUM = SUM+Z*V*DP*DMU*2
	    ENDDO
	  ENDIF
	ENDDO
	PHASUM = SUM/PI
	RETURN
	END
