C	rc-spherical-average.f -- checking RC code for energy conservation
C	gfortran rc-spherical-average.f ~/Ned.a
C
C	by Ned Wright
C	25-MAY-2017
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
	PARAMETER (NP3=12*4**3)
	REAL FBOLMAP(NP3)
	REAL R(3),N(3),P(3)
	PARAMETER (Q=0.38399)
	PARAMETER (NL=74,NC2=47)
	REAL FNU(NL),LAM(NL)
	REAL FBOLSUM(0:3)
	PARAMETER (AU=1.495979E13)	! AU in cm
	REAL PI,TWOPI,DTR
	PARAMETER (PI=3.1415927,TWOPI=6.2831853,DTR=57.29578)
C
	CALL READ_TABLE
	CALL SPHERE(NFACE,75./57.29577,AXIS,PFACE,FCLST)
	DO I=1,NL
	  LAM(I) = 0.1*1000.**((I-1)/73.)
	ENDDO
	DLNWV = ALOG(1000.)/73.
	WRITE (*,*) LAM(NC2)
C
	PV=0.2
	ALB=PV*Q
	THRM=0.2
	EM=0.9
	FC=0.5
	PIR2PV=2
	ASYM=1
	ROT=0
	DOMEGA = 4*3.14159/NP3
	DO ISL=1,5
	  SSL = 0.2*(ISL-0.5)
	  CSL = SQRT(1-SSL**2)
	  P(2)=0
	  P(3)=CSL
	  P(1)=SSL
C
	  DO J=0,3
	    FBOLSUM(J)=0
	  ENDDO
C
	  DO I=1,NP3
	    CALL NESTED_CENPIX(I-1,3,N)
	    DE=1
	    DO J=1,3
	      R(J)=-DE*N(J)
	    ENDDO
C	note: asteroid at (1,0,0), so for N=(1,0,0) & DE=1, need R=(0,0,0)
	    R(1)=R(1)+1
	    CALL PREDICT(R,N,DE,P,THRM,ALB,EM,FC,
	1	PIR2PV,ASYM,ROT,NL,LAM,FNU)
	    FBOL = 0
	    DO K=1,NL
C	Jy to nu*F_nu in erg/cm^2/sec: 3E14/lam(um)*F(Jy)*1E-23
	      FBOL = FBOL+FNU(K)*(3.E-9/LAM(K))
	    ENDDO
	    FBOL = FBOL*DLNWV
	    FNC2 = FNU(NC2)*3.E-9/LAM(NC2)
	    ALPHA = DTR*ACOS(N(1))
	    BC = FNC2/FBOL
	    PHI4Q = 4*PI*AU**2*FBOL/1.0736E16
	    WRITE (16,FMT='(F5.1,F7.4,F7.4,I2,2H P)') ALPHA,BC,PHI4Q,ISL
	    FBOLMAP(I) = FBOL
	    FBOLSUM(0) = FBOLSUM(0)+FBOL
	    DO J=1,3
	      FBOLSUM(J) = FBOLSUM(J)+FBOL*N(J)
	    ENDDO
	  ENDDO
	  DO I=1,NP3
	    FBOLMAP(I) = FBOLMAP(I)/(FBOLSUM(0)/NP3)
	  ENDDO
	  DO J=0,3
	    FBOLSUM(J) = FBOLSUM(J)*DOMEGA*AU**2
	  ENDDO
C	expected 0th moment: 1367*PI*500^2*1E7=1.0736E16
	  WRITE (*,FMT='(1P4E11.3)') FBOLSUM
	  WRITE (UNIT=10+ISL) FBOLMAP
	  CLOSE (UNIT=10+ISL)
	ENDDO
	STOP
	END

	SUBROUTINE PREDICT(R,N,DE,P,THRM,ALB,EM,
	1	FC,PIR2PV,ASYM,ROT,NL,LAM,FNU)
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
	REAL G(NTIME)	! normalized thermal response
	REAL TBRITE(12,2)
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
	REAL TBT2(9,2),LX(NLMX)  ! variables for interpolating TB
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
C	correction for IR albedo slope
	T0 = TAT1AU*((1-ALBEDO*(1+0.2*(PIR2PV-1)))/EMISSIVITY/DS**2)**0.25
	DO ILAM=1,NL
	  LX(ILAM) = ALOG(HCK*10000.0/(LAM(ILAM)*T0))/ALOG(2.)
	ENDDO
C	calculate the fluxes seen at Earth
	CALL FLUXES(ELAT,ELHA,NLAT,LATS,WLAT,NFACE,NTIME,AXIS,
	1	PFACE,FCLST,INSOL,T,TBRITE)
	TBNEG = .FALSE.
	DO I=1,12
	  IF (TBRITE(I,1).LT.0 .OR. TBRITE(I,2).LT.0) TBNEG = .TRUE.
	ENDDO
	IF (TBNEG) THEN
	  WRITE (*,*) 'TBRITE < 0'
	  WRITE (*,*) ELAT,ELHA,SUNLAT,THRM
	  WRITE (*,*) TBRITE
	  DO ILAM = 1,NL
	    FNU(ILAM) = 0
	  ENDDO
	  RETURN
	ENDIF
	VSUM = TBRITE(12,1)*(ALBEDO/(1-ALBEDO))
	DO J=1,2
	  CALL SPLSET(9,LXT,TBRITE(1,J),TBT2(1,J))
	ENDDO
	DO ILAM = 1,NL
	  XB = HCK*10000./
	1	(LAM(ILAM)*T0*SPLINE(LX(ILAM),9,LXT,TBRITE,TBT2))
	  F1 = 0.
	  IF (XB.LT.70.) F1 = 1./(EXP(XB)-1.)
	  XB = HCK*10000./(LAM(ILAM)*T0*
	1	SPLINE(LX(ILAM),9,LXT,TBRITE(1,2),TBT2(1,2)))
	  F2 = 0.
	  IF (XB.LT.70.) F2 = 1./(EXP(XB)-1.)
	  F = FC*F1+(1-FC)*F2
	  FSUM(ILAM) = THCC*OMEGA*(10000./LAM(ILAM))**4*EMISSIVITY*F
	  XS = HCK*10000./(TSUN*LAM(ILAM))
C       the solar flux, nu*F_\nu, in erg/cm^2/sec/oct
          FSUM(ILAM) = THCC*OMEGA*(10000./LAM(ILAM))**4*EMISSIVITY*F
          XS = HCK*10000./(TSUN*LAM(ILAM))
C       the solar flux, nu*F_\nu, in erg/cm^2/sec/oct. 1.056 is a fudge factor to make
C       the V and W1 fluxes correct
          FSUN = 1.056*SOLCON*(XS**4/(EXP(XS)-1))/(DS**2*(PI**4)/15)
	  IF (LAM(ILAM).LT.0.4) FSUN = 0.6*FSUN ! to make bolometric flux right
C	face-on white disk flux at 1 AU for 1 km is
C	(PI*R**2)*FSUN/(PI*AU**2) = FSUN*(R/AU)**2
	  PV = ALBEDO/Q
	  PIR=PV*EXP(MAX(0.,MIN(1.,LOG(LAM(ILAM)/0.55)/LOG(3.4/0.55)))*
	1	LOG(PIR2PV))
	  FHG = PHI*PIR*FSUN*(5.E4/(DE*AU))**2
C	write (*,*) lam(ilam),fhg,phi,pir,fsun
C	stop
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
C	TBRITE (*,1) is for the cratered model, (*,2) is for a flat model
C	ALBEDO/(1-ALBEDO)*TBRITE(12) gives 1 for a white diffuse
C
	PARAMETER (NFMAX=127,NTMAX=32)
	REAL TBAR(0:NFMAX),WATE(NFMAX)
	PARAMETER (NPXL=12*4**6)
	REAL AXIS(3,NPXL)
	REAL PFACE(3,NPXL)
	INTEGER*2 FCLST(NPXL)
	REAL T(NTIME,0:NFACE,NLAT)	! the temperatures
	REAL TBRITE(12,2)
	REAL INSOL(NTIME,0:NFACE,NLAT)         ! solar input
	REAL LATS(NLAT)	! the latitudes
	REAL WLAT(NLAT)	! weight for latitude integration
	REAL*8 X
	REAL PI,TWOPI,DTR
	PARAMETER (PI=3.1415927,TWOPI=6.2831853,DTR=57.29578)
C
	WSUM = 0.
	DO J=1,2	
	  DO I=1,12
	    TBRITE(I,J) = 0.
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
	          TBRITE(10,1) = TBRITE(10,1) + W*T(I,J,K)**4
	          TBRITE(11,1) = TBRITE(11,1) + W*TBAR(J)
	          TBRITE(12,1) = TBRITE(12,1) + W*INSOL(I,J,K)
	          X = 0.0625/T(I,J,K)
	          DO L=1,9
	            IF (X.LE.70.) 
	1		TBRITE(L,1) = TBRITE(L,1) + W/(EXP(X)-1.)
	            X = X+X
	          ENDDO
	        ENDIF
	        WSUM = WSUM+W
	      ENDDO	! J=1,NFACE
C	repeat for the flat model in FACE=0
	      W = WLAT(K)*SW
	      IF (T(I,0,K).GT.0.) THEN
	        TBRITE(10,2) = TBRITE(10,2) + W*T(I,0,K)**4
	        TBRITE(11,2) = TBRITE(11,2) + W*TBAR(0)
	        TBRITE(12,2) = TBRITE(12,2) + W*INSOL(I,0,K)
	        X = 0.0625/T(I,0,K)
	        DO L=1,9
	          IF (X.LE.70.) TBRITE(L,2) = TBRITE(L,2) + W/(EXP(X)-1.)
	          X = X+X
	        ENDDO
	      ENDIF
	    ENDIF
	  ENDDO	! I=1,NTIME
	ENDDO	! K=1,NLAT
C	normalize
	DO J=1,2
	  DO I=1,12
	    TBRITE(I,J) = TBRITE(I,J)/WSUM
	  ENDDO
C	convert from flux to brightness T
	  TBRITE(10,J) = TBRITE(10,J)**0.25
	  X = 0.0625
	  DO I=1,9
	    IF (TBRITE(I,J).GT.0.) TBRITE(I,J) = X/ALOG(1.+1./TBRITE(I,J))
	    IF (I.GE.3 .AND. TBRITE(I,J).LE.0.) 
	1	TBRITE(I,J) = 2*TBRITE(I-1,J)-TBRITE(I-2,J)
	    X = X+X
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

	SUBROUTINE READ_TABLE
	PARAMETER (NFACE=127,NLAT=16,NTIME=32,NSL=31,NTHRM=46)
	REAL INSOLT(NTIME,0:NFACE,NLAT,NSL)		! solar input
	REAL TT(NTIME,0:NFACE,NLAT,NTHRM,NSL)	! the temperatures
C
	COMMON /TABLES/ INSOLT,TT
C
	OPEN (UNIT=41,FILE='/Users/ned/asteroid/rc-tables-steep.dat',
	1	FORM='unformatted')
	DO ISL=1,NSL
	  READ (41) (((INSOLT(I,J,K,ISL),I=1,NTIME),J=0,NFACE),K=1,NLAT)
	  DO IT=1,NTHRM
	    READ (41) (((TT(I,J,K,IT,ISL),I=1,NTIME),J=0,NFACE),K=1,NLAT)
	  ENDDO
	ENDDO
	CLOSE (UNIT=41)
C
	SOLMIN = 1.E30
	SOLMAX = -SOLMIN
	DO ISL=1,NSL
	  DO K=1,NLAT
	    DO J=1,NFACE
	      DO I=1,NTIME
	        SOLMIN=MIN(SOLMIN,INSOLT(I,J,K,ISL))
	        SOLMAX=MAX(SOLMAX,INSOLT(I,J,K,ISL))
	      ENDDO
	    ENDDO
	  ENDDO
	ENDDO
	WRITE (*,*) SOLMIN,SOLMAX
	THMIN = 1.E30
	THMAX = -THMIN
	DO ISL=1,NSL
	  DO IT=1,NTHRM
	    DO K=1,NLAT
	      DO J=1,NFACE
	        DO I=1,NTIME
	          THMIN=MIN(THMIN,TT(I,J,K,IT,ISL))
	          THMAX=MAX(THMAX,TT(I,J,K,IT,ISL))
	        ENDDO
	      ENDDO
	    ENDDO
	  ENDDO
	ENDDO
	WRITE (*,*) THMIN,THMAX
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

