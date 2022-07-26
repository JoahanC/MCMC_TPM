C	make-table-steep.f
C	gfortran -o make-table-steep make-table-steep.f ~/Ned.a
C	computes tables of facet temperatures for 0 & 75^deg angle
C	writes to fort.11
	TSTART = SECNDS(0.)
	CALL MAKETABLESTEEP(0.06)
	TDONE = SECNDS(0.)
	DT = TDONE-TSTART
	WRITE (*,*) DT
	STOP
	END

	SUBROUTINE MAKETABLESTEEP(ALB)
C
C	ALB is the bolometric Bond albedo
C
	REAL THRM
C
	REAL OBJ(3),Z(3),X(3),Y(3)
	INTEGER NFACE,NLAT,NTIME,NLIST,NPXL
	INTEGER I, K
	PARAMETER (NFACE=127,NLAT=16,NTIME=32)
	PARAMETER (NPXL=12*4**6)
	REAL AXIS(3,NPXL)
	REAL PFACE(3,NPXL)
	INTEGER*2 FCLST(NPXL)
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
	REAL AU /1.495979E13 /	! AU in cm
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
	CALL MAKE_G(NTIME,G)
C
	THMAX = 75
	THMAX=THMAX/DTR
	ALBEDO = ALB
	EMISSIVITY = EM
C	setup the latitude table
	DO I=1,NLAT
	  SL = -1.+(2.*I-1.)/NLAT
	  CL = SQRT(1.-SL**2)
	  LATS(I) = ATAN(SL/CL)
	  WLAT(I) = 2./NLAT
	ENDDO
C	set up axes of the crater facets
	CALL SPHERE(NFACE,THMAX,AXIS,PFACE,FCLST)
	DO ISL=1,31
	  SUNLAT = PI*(-96+6*ISL)/180.
	  DO K=1,NLAT
C	compute the solar input
	    CALL SOLAR(SUNLAT,LATS(K),NFACE,NTIME,THMAX,PFACE,AXIS,
	1	FCLST,ALBEDO,INSOL(1,0,K))
	  ENDDO
	  WRITE (11) INSOL
	  DO ITHRM=1,46
C	x = THRM/(THRM+1.6) = (0..45)/45 so (x-1)THRM+1.6*x=0, THRM=-1.6*x/(x-1)
	    THRM = 3.E3
	    IF (ITHRM.NE.46) THEN
	      XXX = (ITHRM-1.)/45.
	      THRM = 1.6*XXX/(1-XXX)
	    ENDIF
	    DO K=1,NLAT
	      CALL THERMAL(NFACE,NTIME,THRM,G,THMAX,INSOL(1,0,K),
	1	T(1,0,K))
	    ENDDO
	    IF (ITHRM.EQ.46) THEN
	      DO K=1,NLAT
	        DO J=0,NFACE
	          T4SUM=0
	          DO I=1,NTIME
	            T4SUM=T4SUM+T(I,J,K)**4
	          ENDDO
	          TBAR = (T4SUM/NTIME)**0.25
	          DO I=1,NTIME
	            T(I,J,K)=TBAR
	          ENDDO
	        ENDDO
	      ENDDO
	    ENDIF
	    WRITE (11) T
	  ENDDO
	ENDDO
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
	REAL INSOL(NTIME,0:NFACE)		! solar input
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
	CL = COS(LAT)
	SL = SIN(LAT)
	CSL = COS(SUNLAT)
	SSL = SIN(SUNLAT)
	DO I=1,NTIME
	  HA = ((I-0.5)/NTIME-0.5)*TWOPI
	  CHA = COS(HA)
	  CALL LWATE(SUNLAT,HA,LAT,NFACE,AXIS,PFACE,FCLST,W)
	  TOTAL = 0
C	the flat surface case. NOTE - (1-ALBEDO) factor applied later
	  INSOL(I,0) = MAX(0.,SL*SSL+CL*CSL*CHA)
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
C	add in scattered sunlight. NOTE (1-ALBEDO) applied later
	  DO J=1,NFACE
	    INSOL(I,J) = (INSOL(I,J) + F*AVERAGE/(1-F))
	  ENDDO
	ENDDO
	RETURN
	END

	SUBROUTINE THERMAL(NFACE,NTIME,X,G,THMAX,INSOL,T)
C	does the thermal inertia calculation
	REAL INSOL(NTIME,0:NFACE)		! solar input
	REAL T(NTIME,0:NFACE)	! the temperatures
	PARAMETER (NFMAX=127,NTMAX=32)
	PARAMETER (XMAX=10000.0)
	REAL*8 A(NTMAX,NTMAX+1)
	REAL G(NTIME)
C
	D = (1.-COS(THMAX))/2./NFACE
C	initial guess
	DO J=0,NFACE
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
	DO JJ=0,NFACE
	  DO I=1,NTIME
C	infrared emission - solar input
	    S = T(I,JJ)**4-INSOL(I,JJ)
C
C	add in emission from other facets
C
	    IF (JJ.NE.0) THEN
	      DO K=1,NFACE
	        S = S-D*T(I,K)**4
	      ENDDO
	    ENDIF
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
	  CALL DMATINV(A,NTMAX,NTIME,1,IFF)
	  IF (IFF.NE.0) STOP 'bad IFF'
	  ER = 0.
	  DO I=1,NTIME
	    ERR = ABS(A(I,NTIME+1))/T(I,JJ)
	    ER = MAX(ER,ERR)
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
