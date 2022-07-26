C	NEOCAM-rc-predict.f
C	gfortran -o NEOCAM-rc-predict NEOCAM-rc-predict.f ~/Ned.a
C
C	23-APR-2016 modifed from NEOCAM-rc-MCMC.f
C	22-APR-2016 modified for NEOCAM bands and gfortran
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
	PARAMETER (NDMX=99)
	PARAMETER (NB=2)
	REAL*4 PHT(NDMX),DET(NDMX),RAT(NDMX),DECT(NDMX),DST(NDMX)
	REAL*4 SIGT(NB,NDMX),MAGT(NB,NDMX)
	REAL*8 MJDT(NDMX)
	COMMON /NRMLPTS/ MJDT,H,MAGT,DH,SIGT,DST,DET,PHT,RAT,DECT,NDM
	COMMON /PREDICTED/ WMAGT(NB,NDMX)
	COMMON /BACKCHAN/ DIA,SDIA
	REAL*8 MJD,JD(NDMX)
	REAL*8 T1JAN2010,T,JD1JAN2010
	REAL*8 TSTART,TWALL
	REAL*8 JULIAN_DAY,GMT2T68
	REAL*4 EA(3),A(3),E(3)
	REAL*4 TE(3)
	REAL*4 WMAG(4),SIG(4)
	REAL*4 P(4)
	REAL*4 P1(3),P2(3),P3(3)
	CHARACTER*11 GMT
	CHARACTER*7 NAME
C
	REAL*4 SIGMIN(NB) /NB*0.05/
	REAL*4 F0(NB) /170.663,63.182/

C	these are CBE 1 sigma per visit, so a quad could be 2x less
	REAL SNC1(9,5) 
	1	/9.41,6.75,5.25,4.36,3.80,3.43,3.17,2.99,2.87,
	2	7.62,5.88,4.79,4.09,3.63,3.32,3.09,2.93,2.82,
	3	5.69,4.79,4.16,3.71,3.39,3.16,2.99,2.86,2.76,
	4	4.57,4.08,3.70,3.41,3.19,3.02,2.89,2.79,2.72,
	5	3.87,3.61,3.38,3.19,3.04,2.92,2.82,2.74,2.68/

	REAL SNC2(9,5) 
	1	/31.52,25.02,20.77,17.92,15.90,14.43,13.33,12.51,11.89,
	2	25.57,21.58,18.65,16.49,14.89,13.67,12.75,12.04,11.50,
	3	19.38,17.36,15.74,14.46,13.45,12.63,11.96,11.43,11.00,
	4	15.85,14.76,13.80,12.99,12.31,11.75,11.28,10.88,10.56,
	5	13.65,13.05,12.48,11.96,11.50,11.11,10.77,10.48,10.24/

	REAL SIGF(NB,NDMX)
C
	EXTERNAL RAN
C
	CALL RANDOMIZE
	CALL INIT_ORBIT
	NDM=0
	CALL TIME_TAG(2010,1,1,0,0,0.,T1JAN2010)
	JD1JAN2010 = JULIAN_DAY(T1JAN2010)
	DTR = 45/ATAN(1.)
	READ (*,FMT='(3F7.3,1X,A,F8.3,F6.3,F10.5,2F8.5)')
	1	H,DH,PERIOD,NAME,DIA,ALBEDO,THETA1,RAP,DECP
	P(1) = RAP
	P(2) = DECP
	P(3) = LOG(THETA1)
	P(4) = LOG(ALBEDO)
	WRITE (*,FMT='(A)') 'Enter MJD,RA,Dec,Delta,Diameter,p_V,'//
	1	'Theta1,Pole RA,Dec:'
	WRITE (*,FMT='(A)') '   MJD      RA      DEC     DE      DS'//
	1	'    PHASE NC12 & sigmas'
10	READ (*,900,END=100) MJD,RA,DEC,DELTA
900	FORMAT(F11.5,2F11.7,F8.4,4(F7.3,F6.3))
	NDM=NDM+1
	RAT(NDM) = RA
	DECT(NDM) = DEC
	DET(NDM) = DELTA
	MJDT(NDM) = MJD
	CALL NEOCAM_ORBIT(MJDT(NDM),E)
	CALL CTOE(E,TE)
	CALL C_TO_RADEC(TE,SUNLONG,ELAT)
	CALL RADEC_TO_C(RAT(NDM)/DTR,DECT(NDM)/DTR,EA)
	DO J=1,3
	  EA(J) = EA(J)*DET(NDM)
	  A(J) = EA(J)+E(J)
	ENDDO
	CALL CTOE(EA,TE)
	CALL C_TO_RADEC(TE,ELONG,ELAT)
C	inteprolate to get sigmas
	DELONG = DTR*ABS(SUNLONG-ELONG)
	X = (DELONG-45)/10
	IX = X
	IX = MAX(0,MIN(IX,7))
	X = X-IX
	IX = IX+2
	Y = DTR*ABS(ELAT)/10
C	use quadratic near ecliptic
	IF (Y.LT.1) Y=Y*Y
	IY = Y
	IY = MIN(IY,3)
	Y = Y-IY
	IY = IY+2
	SIGF(1,NDM) = X*Y*SNC1(IX,IY)+X*(1-Y)*SNC1(IX,IY-1)+
	1	(1-X)*Y*SNC1(IX-1,IY)+(1-X)*(1-Y)*SNC1(IX-1,IY-1)
	SIGF(2,NDM) = X*Y*SNC2(IX,IY)+X*(1-Y)*SNC2(IX,IY-1)+
	1	(1-X)*Y*SNC2(IX-1,IY)+(1-X)*(1-Y)*SNC2(IX-1,IY-1)
C
	DST(NDM) = SQRT(DOT(A,A))
	PHT(NDM) = DTR*ACOS(DOT(EA,A)/(DST(NDM)*DET(NDM)))
	WRITE (*,901) MJD,RAT(NDM),DECT(NDM),DET(NDM),DST(NDM),
	1	PHT(NDM),SIGF(1,NDM),SIGF(2,NDM)
901	FORMAT(F9.3,2F8.3,2F8.5,F7.2,2F7.3,F10.2,F6.2)
	GO TO 10
100	E = ERROR(P)
	WRITE (*,FMT='(5F9.3)') DIA,P
	WRITE (*,FMT='(3F7.3,1X,A)') H,DH,PER,NAME
	DO I=1,NDM
	  MJD = MJDT(I)
	  DO J=1,NB
	    F = 1.E6*F0(J)/10.**(0.4*WMAGT(J,I))
	    SIGT(J,I) = SQRT(SIGMIN(J)**2+(2.5*LOG10(1+SIGF(J,I)/F))**2)
	  ENDDO
	  WRITE (*,902) MJD,RAT(I),DECT(I),DET(I),
	1	(WMAGT(J,I),SIGT(J,I),J=1,NB)
	ENDDO
902	FORMAT(F9.3,1H,,F8.4,1H,,F8.4,1H,,F7.5,2(1H,,F6.3,1H,,F4.2))
	STOP
	END

	SUBROUTINE INIT_ORBIT
C
	REAL*8 MJDT(256),OBST(3,256)
	COMMON /ORBITDATA/ MJDT,OBST
	OPEN (UNIT=1,FILE='/Users/wright/missions/NEOCAM/orbit.dat')
	DO I=1,256
	  READ (1,*) MJDT(I),(OBST(J,I),J=1,3)
	ENDDO
	CLOSE (UNIT=1)
	RETURN
	END

	SUBROUTINE NEOCAM_ORBIT(MJD,C)
	REAL*8 MJDT(256),OBST(3,256)
	COMMON /ORBITDATA/ MJDT,OBST
C	interpolate in Table
	REAL*4 E(3),C(3)
	REAL*8 MJD
	REAL*8 W(4)
C
	IF (MJD.LT.MJDT(1)-30 .OR. MJD.GT.MJDT(256)+30) STOP 'bad MJD'
	DO K=2,254
	  J=K
	  IF (MJD.LT.MJDT(J+1)) GO TO 10
	ENDDO
10	DO I=1,4
	  W(I)=1
	  DO K=1,4
	    IF (K.NE.I) 
	1	W(I)=W(I)*(MJD-MJDT(K+J-2))/(MJDT(I+J-2)-MJDT(K+J-2))
	  ENDDO
	ENDDO
	DO K=1,3
	  E(K)=0
	  DO I=1,4
	    E(K) = E(K)+OBST(K,I+J-2)*W(I)
	  ENDDO
	ENDDO
	CALL ETOC(E,C)
	RETURN
	END

	REAL*4 FUNCTION ERROR(P)
C
	REAL*4 P(4)
C       P(1..2) = RA,Dec in radians of pole
C	P(3) = ln(theta1)
C	P(4) = ln(p_V)
C
	PARAMETER (NDMX=99)
	PARAMETER (NBAND=2,NL=6)
	REAL*4 PHT(NDMX),DET(NDMX),RAT(NDMX),DECT(NDMX),DST(NDMX)
	REAL*4 SIGT(NBAND,NDMX),MAGT(NBAND,NDMX)
	REAL*8 MJDT(NDMX)
	REAL*8 DIAM(NBAND*NDMX+1),SDIAM(NBAND*NDMX+1),LND,SLND,RCS
	REAL*8 AREA(NBAND*NDMX+1),SAREA(NBAND*NDMX+1),DSQ,SDSQ
	COMMON /NRMLPTS/ MJDT,H,MAGT,DH,SIGT,DST,DET,PHT,RAT,DECT,NDM
	COMMON /PREDICTED/ WMAGT(NBAND,NDMX)
	COMMON /BACKCHAN/ DIA,SDIA
	COMMON /WATERMARK/ BESTEVER,PB(4)
C
	REAL*4 R(3),N(3),POLE(3)
	PARAMETER (DTR = 57.2957795)
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
	REAL*4 MICRONS,NUFNU
C	parameters of albedo prior
	PARAMETER (PK1=0.03,FD=0.251,PK2=0.167)
C
	DO I=1,NLT
	  WVLT(I) = 0.55*(28./0.55)**((I-1)/(NLT-1.))
	ENDDO
	CALL RADEC_TO_C(P(1),P(2),POLE)
C
C	range restrictions to protect the RC code
	TH1MAX = 10000
	IF (P(3).LT.23.) THEN
	  THETA1 = EXP(P(3))
	  THETA1 = TH1MAX*THETA1/(TH1MAX+THETA1)
	ELSE
	  THETA1 = TH1MAX
	ENDIF
	PVMAX = 1.5
	IF (P(4).GT.23.) THEN
	  PV = PVMAX
	ELSE
	  PV = EXP(P(4))
	  PV = PVMAX*PV/(PVMAX+PV)
	ENDIF
C
	ALB = 0.38399*PV
C
	EM = 0.95
	ASYM = 1.
	ROT = 0.
	DO J=1,NDM
	  CALL NEOCAM_ORBIT(MJDT(J),R)
	  CALL RADEC_TO_C(RAT(J)/DTR,DECT(J)/DTR,N)
	  THETA = THETA1*DST(J)**1.5
	  CALL PREDICT(R,N,DET(J),POLE,THETA,ALB,EM,ASYM,ROT,NALL,WVALL,F)
	  DO I=1,NALL
	    FALL(I,J) = F(I)
	  ENDDO
	  WMAGT(1,J)=-2.5*LOG10((W(1)*F(1)+W(2)*F(2))/F0(1))
	  WMAGT(2,J)=-2.5*LOG10((W(3)*F(3)+W(4)*F(4)+W(5)*F(5)+
	1	W(6)*F(6))/F0(2))
	ENDDO
C	scale to the diameter
	DM = -5.*ALOG10(DIA)
	DO J=1,NDM
	  DO I=1,NALL
	    FALL(I,J) = FALL(I,J)*DIA**2
	  ENDDO
	  DO I=1,NBAND
	    WMAGT(I,J) = WMAGT(I,J)+DM
	  ENDDO
	ENDDO
	ERROR=0
	RCS=0
C	always write fort.22
	WRITE (22,FMT='(1H%,F11.6,3F12.6,F12.4,F9.4)') P
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
	    DM = 0	! for PREDICT make it exact
	    FJY = FALL(I+NL,J)/10.**(0.4*DM)
	    MICRONS = WVALL(I+NL)
	    NUFNU = 1.E-23*FJY*(3.E14/MICRONS)
	    SGT = 0.1
	    WRITE (22,FMT='(2F7.3,1PE10.3,3H QQ)') MICRONS,SGT,NUFNU
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
C
	RETURN
	END

	SUBROUTINE ROBUST(NP,X,SG,MU,SIGMA,RCS)
C
C	NP input points in X with sigma's in SG
C	return estimated centroid MU and sigma SIGMA
C	with robust chi^2 in RCS
C
	REAL*8 X(NP),SG(NP),MU,SIGMA,RCS
	REAL*8 W,SW,SWX,CHI2
	PARAMETER (NPMX=500)
	REAL*4 DX(NPMX)
	REAL*4 Y(NPMX)
	REAL*8 XMIN,XMAX,XC
C
	IF (NP.LE.0) STOP 'No data'
	IF (NP.GT.NPMX) STOP 'Too much data'
	IF (NP.EQ.1) THEN
	  MU = X(1)
	  SIGMA = SG(1)
	  RCS = 0
	  RETURN
	ENDIF
	XMIN = X(1)
	XMAX = XMIN
	DO I=2,NP
	  XMIN = MIN(XMIN,X(I))
	  XMAX = MAX(XMAX,X(I))
	ENDDO
	XC = (XMIN+XMAX)/2
	DMAX = XMAX-XMIN
	IF (DMAX.LE.0) THEN
	  MU = XC
	  SIGMA = 0
	  DO I=1,NP
	    SIGMA = SIGMA+1/SG(I)**2
	  ENDDO
	  SIGMA = SQRT(1/SIGMA)
	  RCS = 0
	  RETURN
	ENDIF
	DO I=1,NP
	  DX(I) = X(I)-XC
	ENDDO
	CALL RSORT(NP,DX)
	NX = NPMX/NP
	NX = MAX(5,NX)
	N = NX*(NP-1)+1
C	WRITE(*,FMT='(10F8.4)') (DX(I),I=1,NP)
	DX(N) = DX(NP)
	DO II=1,(NP-1)
	  I = NP-II
	  DO JX=1,NX
	    F = (NX-JX+1.)/NX
	    J = NX*(I-1)+(NX+2-JX)
	    DX(J) = F*DX(I+1)+(1-F)*DX(I)
	  ENDDO
	ENDDO
C	WRITE(*,FMT='(10F8.4)') (DX(I),I=1,N)
	DO I=1,N
	  Y(I)=0
	  T = DX(I)
	  DO J=1,NP
	    E = (X(J)-XC-T)/SG(J)
	    IF (ABS(E).GT.2) THEN
	      CS = 4*ABS(E)-4
	    ELSE
	      CS = E*E
	    ENDIF
	    Y(I) = Y(I)+CS
	  ENDDO
	ENDDO
	BEST = Y(1)
	IB = 1
	DO I=2,N
	  IF (Y(I).LT.BEST) THEN
	    IB = I
	    BEST = Y(I)
	  ENDIF
	ENDDO
C
	IF (IB.EQ.1) THEN
	  D = (DX(3)-DX(1))/2
	ELSE IF (IB.EQ.N) THEN
	  D = (DX(N)-DX(N-2))/2
	ELSE
	  D = (DX(IB+1)-DX(IB-1))/2
	ENDIF
	D = MAX(D,DMAX/1024)
	DDY = 0
	DO WHILE (DDY.LE.0 .AND. D.LT.DMAX)
	  D = 2*D
	  DO I=1,3
	    Y(I) = 0
	    T = DX(IB)+(I-2)*D
	    DO J=1,NP
	      E = (X(J)-XC-T)/SG(J)
	      IF (ABS(E).GT.2) THEN
	        CS = 4*ABS(E)-4
	      ELSE
	        CS = E*E
	      ENDIF
	      Y(I) = Y(I)+CS
	    ENDDO
	  ENDDO
	  DDY = Y(3)-2*Y(2)+Y(1)
	  DY = (Y(3)-Y(1))/2
	ENDDO
	MU = XC+DX(IB)-D*DY/DDY
	SIGMA = D*SQRT(2/DDY)
	RCS = 0
	DO J=1,NP
	  E = (X(J)-MU)/SG(J)
	  IF (ABS(E).GT.2) THEN
	    CS = 4*ABS(E)-4
	  ELSE
	    CS = E*E
	  ENDIF
	  RCS = RCS+CS
	ENDDO
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
	PARAMETER (SOLCON = 1367.E3)	! solar constant in erg/cm^2/sec
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
	PARAMETER (PIR2PV = 1.3)
C
	INTEGER IMC,NEE
	REAL RMS
	REAL*8 SE,SSEE,SEE
	INTEGER IY,IZ,IR,JR
	REAL RAN
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
	T0 = TAT1AU*((1-ALBEDO*(1+0.19*(PIR2PV-1)))/EMISSIVITY/DS**2)**0.25
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
	1	2.5E-5/DS**2/(EXP(14387.752/LAM(ILAM)/TSUN)-1)
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
