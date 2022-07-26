C	lightcurve.f -- modeling thermal emission of 
C	gfortran lightcurve.f ~/Ned.a
C	rotating cratered asteroids for NGSS viewing conditions
C	by Ned Wright
C
C       sub-solar temperature increased to 1.01620 1.03398
C       for 30 & 45 degrees THMAX
C
	IMPLICIT NONE
	INTEGER NFACE,NLAT,NTIME,NLIST,NPXL
	INTEGER I, J, K, IT
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
	INTEGER IRING	! counter for ALPHA loop
	REAL CA,SA	! cosine and sine of ALPHA
	REAL TH		! clock angle of EVEC around Sun vector
	INTEGER ITH,NTH	! counter and number of TH's
	REAL DTH	! delta TH
	REAL ELHA	! local hour angle of EVEC
	REAL ELAT	! latitude of Earth
	REAL EVEC(3)
	REAL G(NTIME)	! normalized thermal response
	REAL TBRITE(12,NTIME),AREA(NTIME)
	REAL XLIST(13),X	! thermal inertia parameter
	INTEGER IX
	INTEGER LIST	! counter for XLIST
	REAL TMP	! temporary variable used to rotate vectors
	REAL THCC /1.19104272E-5/
	REAL OMEGA
	REAL HCK /1.4387752/	! h*c/k
	REAL XB		! h*nu/k*Tb
	REAL SOLCON /1367.E3/	! solar constant in erg/cm^2/sec
	REAL SIGSB /5.6704E-5/	! Stefan-Boltzmann constant
	REAL TAT1AU /394.03854311248/
	REAL ALBEDO, EMISSIVITY
	REAL ETA,DIA
	REAL T0
	REAL DE		! distance to Earth in AU
	REAL DS		! distance to Sun in AU
	REAL AU /1.495979E13 /	! AU in cm
	REAL LXT(9) /-4.,-3.,-2.,-1.,0.,1.,2.,3.,4./
	INTEGER NLAM
	PARAMETER (NLAM=5)
	REAL LAMBDA(NLAM) /3.36,4.6,7.75,11.56,22.09/
	REAL NU
	REAL TBT2(9),LX(NLAM)  ! variables for interpolating TB
	REAL FSUM(NLAM),F,FLAM(NLAM)
	INTEGER IFNU(NLAM)
	REAL VSUM,FOPT
	INTEGER IDS,IA,IDE,ILAM
	REAL JY,MAB
	REAL SPLINE
	REAL PI,TWOPI,DTR
	PARAMETER (PI=3.1415927,TWOPI=6.2831853,DTR=57.29578)
C
	INTEGER IMC
	REAL RMS
	INTEGER IY,IZ,IR,JR
C	axes of ellipsoid A >= B >= C
	REAL A,B,C,WATES(NLAT,NTIME)
	INTEGER IWATES(NLAT)
	REAL TMIN,TMAX
C
	CALL MAKE_G(NTIME,G)
C
	THMAX = 45
	THMAX=THMAX/DTR
	ALBEDO = 0.06	! was 0.1 on 3/30/2007
	EMISSIVITY = 0.9
	DO I=1,NLAT
	  SL = -1.+(2.*I-1.)/NLAT
	  CL = SQRT(1.-SL**2)
	  LATS(I) = ATAN(SL/CL)
	  WLAT(I) = 2./NLAT
	ENDDO
C	set up axes of the crater facets
	CALL SPHERE(NFACE,THMAX,AXIS,PFACE,FCLST)
	WRITE (4,FMT='(A)') '% DS,ISL,PHASE,TH,X,ETA,DIA,FLAM'
	WRITE (3,*)
	CLOSE (UNIT=3)
C	WRITE (*,*) 'Enter Dsun,Sunlat,Phase,Clock,X,A,B,C'
	WRITE (*,*) 'Enter Dsun,Sunlat,Elat,ELHA,X,A,B,C'
C	READ (*,*) DS,SUNLAT,PHASE,TH,X,A,B,C
	READ (*,*) DS,SUNLAT,ELAT,ELHA,X,A,B,C
	WRITE (*,*) DS,SUNLAT,ELAT,ELHA,X,A,B,C
	ELAT = ELAT/DTR
	ELHA = ELHA/DTR
	SUNLAT = SUNLAT/DTR
	CSL = COS(SUNLAT)
	SSL = SIN(SUNLAT)
	PHASE = DTR*ACOS(SSL*SIN(ELAT)+CSL*COS(ELAT)*COS(ELHA))
	WRITE (*,*) 'Phase = ',PHASE
	CALL EGG_WEIGHTS(NLAT,NTIME,WATES)
	CALL ELLIPSOID_WEIGHTS(A,B,C,NLAT,NTIME,WATES)
	WRITE (*,FMT='(A)')
	DO I=1,NTIME
	  DO K=1,NLAT
	    IWATES(K) = 1000*WATES(K,I)+0.5
	  ENDDO
	  WRITE (*,FMT='(16I5)') IWATES
	ENDDO
	WRITE (*,FMT='(A)')
C
	  ISL = DTR*SUNLAT+0.5
	  X = 1.0
	  DO K=1,NLAT
C	compute the solar input
	    CALL SOLAR(SUNLAT,LATS(K),NFACE,NTIME,THMAX,PFACE,AXIS,
	1	FCLST,ALBEDO,INSOL(1,1,K))
C	    IF(K.EQ.NLAT/2) WRITE (*,FMT='(10F8.5)') (INSOL(I,1,K),I=1,NTIME)
C	compute the temperatures of the facets
	    CALL THERMAL(NFACE,NTIME,X,G,THMAX,INSOL(1,1,K),T(1,1,K))
C	    WRITE (*,FMT='(10F8.5)') (T(I,1,K),I=1,NTIME)
	  ENDDO
	  TMAX = T(1,1,1)
	  TMIN = TMAX
	  DO K=1,NLAT
	    DO J=1,NFACE
	      DO I=1,NTIME
	        TMIN = MIN(TMIN,T(I,J,K))
	        TMAX = MAX(TMAX,T(I,J,K))
	      ENDDO
	    ENDDO
	  ENDDO
	  WRITE (*,FMT='(1P2E14.6)') TMIN,TMAX
C	  WRITE (*,FMT='(A,3F10.4)') ' SUNLAT,THMAX,X:',SUNLAT,THMAX,X
C	generate view directions. 18 angles at correct phase angle for this DE
C	solid angle of 1 km diameter at 1 AU
	  OMEGA = (PI/4)*(1.E5/AU)**2
	  ALPHA = PHASE/DTR
C	fixed 15-MAY-2017:
C	was 30-MAR-2007:	  T0 = TAT1AU*(EMISSIVITY*DS**2)**0.25
	  T0 = TAT1AU/(EMISSIVITY*DS**2)**0.25
	  DO ILAM=1,NLAM
	    NU = 10000/LAMBDA(ILAM)
	    LX(ILAM) = ALOG(HCK*NU/T0)/ALOG(2.)
	  ENDDO
	GO TO 77
	  CA = COS(ALPHA)
	  SA = SIN(ALPHA)
C	set up vector ALPHA from Sun axis with clock angle TH
	  EVEC(1) = CA
	  EVEC(2) = COS(TH/DTR)*SA
	  EVEC(3) = SIN(TH/DTR)*SA
C	rotate by SUNLAT
	  TMP = CSL*EVEC(1)-SSL*EVEC(3)
	  EVEC(3) = CSL*EVEC(3)+SSL*EVEC(1)
	  EVEC(1) = TMP
C	for SUNLAT=90, SSL=1, CSL=0, ALPHA=90, EVEC=(-sin(TH),cos(TH),0)
	  ELHA = ATAN2(EVEC(2),EVEC(1))
	  ELAT = ATAN2(EVEC(3),SQRT(EVEC(1)**2+EVEC(2)**2))
C	calculate the fluxes seen at Earth
77	  CALL FLUXES(ELAT,ELHA,NLAT,LATS,WLAT,NFACE,NTIME,AXIS,
	1	PFACE,FCLST,INSOL,T,WATES,TBRITE,AREA)
C	  do it = 1,11,2
C	    write (*,fmt='(6f10.6)') (tbrite(i,it),i=1,11,2)
C	  enddo
C	  write (*,fmt='(16F5.2)') area
	  DO I=1,12
	    DO IT=1,NTIME
	      IF (TBRITE(I,IT).LE.0) THEN
	        WRITE (*,*) 'TBRITE < 0'
	        WRITE (*,*) ELAT,ELHA
	        WRITE (*,*) TBRITE
	        STOP
	      ENDIF
	    ENDDO
	  ENDDO
	  DO IT=1,NTIME
	    VSUM = TBRITE(12,IT)*(ALBEDO/(1-ALBEDO))
	    CALL SPLSET(9,LXT,TBRITE(1,IT),TBT2)
	    DO ILAM = 1,NLAM
	      NU = 10000/LAMBDA(ILAM)
	      XB = HCK*NU/(T0*SPLINE(LX(ILAM),9,LXT,TBRITE(1,IT),TBT2))
	      F = 0.
	      IF (XB.LT.70.) F = 1./(EXP(XB)-1.)
	      FSUM(ILAM) = THCC*OMEGA*NU**4*EMISSIVITY*F
	      FOPT = THCC*OMEGA*NU**4*VSUM*
	1	2.5E-5/DS**2/(EXP(HCK*NU/5600.)-1)
	      FSUM(ILAM) = AREA(IT)*(FSUM(ILAM)+FOPT)
	      JY = FSUM(ILAM)*LAMBDA(ILAM)/3.E-9
	      FLAM(ILAM) = JY
	      IFNU(ILAM) = 1.E6*JY+0.5
	    ENDDO
C	    CALL FIT_NEATM(ALBEDO,EMISSIVITY,DS,PHASE,2,LAMBDA(6),
C	1	FLAM(6),ETA,DIA)
	    WRITE (*,FMT='(I2,5I6,F7.3,2F7.4,2H P)') 
	1	IT,IFNU,AREA(IT)	!	,ETA,DIA
	  ENDDO
C
	STOP
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
C
C	routine to compute Weights for the NFACE FACETS when the Earth is
C	at ELAT and ELHA.  Surface is at LAT
C	mouth of crater is circle of radius 1
C
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
	IF (EPSILON.GT.0.001) GO TO 100
	RETURN
	END

	SUBROUTINE FLUXES(ELAT,ELHA,NLAT,LATS,WLAT,NFACE,NTIME,AXIS,
	1	PFACE,FCLST,INSOL,T,WATES,TBRITE,AREA)
C
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
	REAL TBRITE(12,NTIME)
	REAL V(12)
	REAL WATES(NLAT,NTIME)
	REAL AREA(NTIME)
	REAL WSUM(NTIME)
	REAL INSOL(NTIME,NFACE,NLAT)         ! solar input
	REAL LATS(NLAT)	! the latitudes
	REAL WLAT(NLAT)	! weight for latitude integration
	REAL*8 X
	REAL PI,TWOPI,DTR
	LOGICAL AOK
	PARAMETER (PI=3.1415927,TWOPI=6.2831853,DTR=57.29578)
C
	DO IT = 1,NTIME
	  WSUM(IT) = 0.	
	  DO I=1,12
	    TBRITE(I,IT) = 0.
	  ENDDO
	ENDDO
C
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
	      IF (WATE(J).LT.0.) WRITE (*,*) 'negative weight:',I,J,K
	      SW = SW+WATE(J)
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
		    V(L) = 0.
	            IF (X.LE.83.) V(L) = W/(EXP(X)-1.)
	            X = X+X
	          ENDDO
	          DO IT = 1,NTIME
	            JT = I+IT-1
	            IF (JT.GT.NTIME) JT = JT-NTIME
	            DO L=1,12
	              TBRITE(L,IT) = TBRITE(L,IT)+V(L)*WATES(K,JT)
	            ENDDO
	          ENDDO
	        ENDIF
	        DO IT = 1,NTIME
	          JT = I+IT-1
	          IF (JT.GT.NTIME) JT = JT-NTIME
	          WSUM(IT) = WSUM(IT)+W*WATES(K,JT)
	        ENDDO
	      ENDDO	! J=1,NFACE
	    ENDIF	! IF (SW>0)
	  ENDDO	! I=1,NTIME
	ENDDO	! K=1,NLAT
C	normalize
	WRITE (*,*)
C	WRITE (*,FMT='(11F7.3)') WSUM
	WRITE (*,*)
	DO IT=1,NTIME
	  DO I=1,12
	    TBRITE(I,IT) = TBRITE(I,IT)/WSUM(IT)
	  ENDDO
	  AREA(IT) = WSUM(IT)*(2./NTIME)
C	convert from flux to brightness T
	  TBRITE(10,IT) = TBRITE(10,IT)**0.25
	  AOK = TBRITE(1,IT).GT.0. 
	  DO I=2,9
	    AOK = AOK.AND.(TBRITE(I,IT).GT.0.)
	  ENDDO
	  IF (.NOT.AOK) WRITE (*,FMT='(1P9E8.1)') (TBRITE(I,IT),I=1,9)
	  X = 0.0625
	  DO I=1,9
	    IF (TBRITE(I,IT).GT.0.) 
	1	TBRITE(I,IT) = X/ALOG(1.+1./TBRITE(I,IT))
	    IF (I.LT.3 .AND. TBRITE(I,IT).LE.0.)  TBRITE (I,IT)=0.
	    IF (I.GE.3 .AND. TBRITE(I,IT).LE.0.) 
	1	TBRITE(I,IT) = 2*TBRITE(I-1,IT)-TBRITE(I-2,IT)
	    X = X+X
	  ENDDO
	  IF (.NOT.AOK) WRITE (*,FMT='(9F8.4)') (TBRITE(I,IT),I=1,9)
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

	SUBROUTINE NEATM(ALB,EPS,DS,PH,ETA,NL,LAM,FNU)
C
C	INPUT:	ALB, the albedo
C		EPS, the emissivity
C		DS, the Sun distance in AU
C		PH, the phase angle in degrees
C		ETA, the beaming parameter
C		NL, the number of wavelengths, must be >= 2
C		LAM(1..NL): the wavelengths in microns
C
C	OUTPUT:	FNU(1..NL): the fluxes in Jy for 1 km diameter at 1 AU
C
	IMPLICIT NONE
	REAL DS,PH,ALB,EPS,ETA,DIA
	INTEGER NL
	REAL LAM(NL),FNU(NL)
C
	REAL TSS,OPTICAL,THERMAL,XIR,XS,AL,AN,FSUN,FLUX,FJY
	INTEGER INU
	REAL*4 LSUN,AU,SB,R,TSUN,HCK,THCC,PI
C	LSUN = solar luminosity in erg/sec
C	AU = astronomical unit in cm
C	SB = Stefan-Boltzman constant in erg/cm^2/sec/K^4
C	R = radius of object in cm
C	TSUN = solar temperature for scattered light
C	HCK = h*c/k in cm*K
	PARAMETER (LSUN = 3.826E33, AU = 1.49E13, SB = 5.67E-5, R=5.E4)
	PARAMETER (TSUN = 5600.0, HCK = 1.4387752)
	PARAMETER (THCC=1.19104272E-5)
	PARAMETER (PI = 3.1415926535897)
C
	REAL PHASUM
C
	TSS = ((1-ALB)*LSUN/(4*PI*DS**2*AU**2)/EPS/SB/ETA)**0.25
C	compute average value of mu over disk for scattering
	OPTICAL = (4*PHASUM(-1.,PH,40)-PHASUM(-1.,PH,20))/3
	DO INU=1,NL
	  AL = LAM(INU)
	  AN = 10000./AL
C	x = h*nu/kT
	  XS = AN*HCK/TSUN
C	the solar flux, nu*F_\nu, in erg/cm^2/sec/oct
	  FSUN = (XS**4/(EXP(XS)-1))/(PI**4/15)*LSUN/(DS*AU)**2/(4*PI)
	  XIR = AN*HCK/TSS
C	compute average value of 1/(e^x-1) over disk for thermal emission
	  THERMAL = (4*PHASUM(XIR,PH,40)-PHASUM(XIR,PH,20))/3
C	flux is average surface brightness * solid angle, [erg/cm^2/sec/oct]
	  FLUX = (THERMAL*EPS*THCC*AN**4 + ALB*OPTICAL*FSUN/PI)
	1	*(PI*R**2/AU**2)
	  FJY = 1.E23*FLUX/(3.E14/AL)
	  FNU(INU) = FJY
	ENDDO
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
	  IF (HNUKTSS.GT.0) V =1./(EXP(HNUKTSS/MU**0.25)-1.)
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

	SUBROUTINE FIT_NEATM(ALB,EPS,DS,PH,NL,LAM,FNU,ETA,DIA)
C
C	INPUT:	ALB, the albedo
C		EPS, the emissivity
C		DS, the Sun distance in AU
C		PH, the phase angle in degrees
C		NL, the number of wavelengths, must be >= 2
C		LAM(1..NL): the wavelengths in microns
C		FNU(1..NL): the fluxes in Jy
C
C	OUTPUT:	ETA, the beaming parameter
C		DIA, the diameter in km if the observer distance is 1 AU
C
	REAL DS,PH,ALB,EPS,ETA,DIA
	INTEGER NL
	REAL LAM(NL),FNU(NL)
	INTEGER I,IT
	REAL EL,EH,X,L,H,SX,SXX,SXL,SL,SXH,SH
C
	REAL FL(20),FH(20)
	IF (NL.LT.2 .OR. NL.GT.20) STOP 'Bad NL'
	ALBEDO = ALB
10	EL = 0.5
	EH = 2.0
	NTRY = 0
	BL = DLNDIA_DLNWVL_NEATM(ALBEDO,EPS,DS,PH,NL,LAM,FNU,EL)
	NTRY = 1
	DO WHILE (BL.LT.0.)
	  EH = EL
	  EL = EL/1.25
	  BL = DLNDIA_DLNWVL_NEATM(ALBEDO,EPS,DS,PH,NL,LAM,FNU,EL)
	  NTRY = NTRY+1
	  IF (NTRY.GT.100) STOP 'in EL loop'
	ENDDO
	BH = DLNDIA_DLNWVL_NEATM(ALBEDO,EPS,DS,PH,NL,LAM,FNU,EH)
	DO WHILE (BH.GT.0.)
	  EH = 1.25*EH
	  BH = DLNDIA_DLNWVL_NEATM(ALBEDO,EPS,DS,PH,NL,LAM,FNU,EH)
	  NTRY = NTRY+1
	  IF (NTRY.GT.100) THEN
	    ALBEDO = 0.0
	    GO TO 10
	  ENDIF
	ENDDO
	DO IT = 1,9
	  ETA = EL*(EH/EL)**(BL/(BL-BH))
	  B = DLNDIA_DLNWVL_NEATM(ALBEDO,EPS,DS,PH,NL,LAM,FNU,ETA)
	  IF (B.GE.0.) THEN
	    EL = ETA
	    BL = B
	  ELSE
	    EH = ETA
	    BH = B
	  ENDIF
	ENDDO
	CALL NEATM(ALBEDO,EPS,DS,PH,ETA,NL,LAM,FL)
	SL = 0
	DO I=1,NL
	  L = ALOG(FNU(I)/FL(I))
	  SL = SL+L
	ENDDO
	DIA = EXP(0.5*SL/NL)
	RETURN
	END

	REAL FUNCTION DLNDIA_DLNWVL_NEATM(ALB,EPS,DS,PH,NL,LAM,FNU,ETA)
C
C	INPUT:	ALB, the albedo
C		EPS, the emissivity
C		DS, the Sun distance in AU
C		PH, the phase angle in degrees
C		NL, the number of wavelengths, must be >= 2
C		LAM(1..NL): the wavelengths in microns
C		FNU(1..NL): the fluxes in Jy
C		ETA, the beaming parameter
C
C	returned value: slope d\ln(Dia)/d\ln(\lambda) for these parameter
C
	REAL DS,PH,ALB,EPS,ETA,DIA
	INTEGER NL
	REAL LAM(NL),FNU(NL)
	INTEGER I,IT
	REAL X,L
	REAL*8 SX,SXX,SXL,SL
C
	REAL FL(20)
	IF (NL.LT.2 .OR. NL.GT.20) STOP 'Bad NL'
	CALL NEATM(ALB,EPS,DS,PH,ETA,NL,LAM,FL)
C	zero sums for fit
	SX = 0
	SXX = 0
	SL = 0
	SXL = 0
C	do sums for fit of ln flux ratio to ln lambda
	DO I=1,NL
	  X = ALOG(LAM(I))
	  L = 0.5*ALOG(FNU(I)/FL(I))
	  SX = SX+X
	  SXX = SXX+X*X
	  SL = SL+L
	  SXL = SXL+X*L
	ENDDO
C	compute slopes of the fits
	DLNDIA_DLNWVL_NEATM = (SXL-SX*SL/NL)/(SXX-SX*SX/NL)
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

	SUBROUTINE EGG_WEIGHTS(NLAT,NTIME,WATES)
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
	B = 1.0
	C = 1.0
	DO I=0,NPIX-1
	  CALL NESTED_CENPIX(I,IRES,R)
	  A = 2.5
	  IF (R(1).LT.0.) A = 1.25
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
