C	read-WISE-rc-MCMC-PJDFC.f
C	gfortran -o read-WISE-rc-MCMC-PJDFC read-WISE-rc-MCMC-PJDFC.f ~/Ned.a
C
C	write D histogram on fort.2, "fine" on fort.32
C	writes D vs period scatter on fort.3
C	writes D vs albedo scatter on fort.4
C	writes hi-res pole position sky map on big-WISE-rc-MCMC.map
C
	PARAMETER (NMCX=1500000,NP=8)
	PARAMETER (NSKIP=1)
	PARAMETER (PVMAX = 1.5, TH1MAX = 1.0E4)
	REAL*4 P(NP,NMCX),D(NMCX),TH(NMCX),PV(NMCX)
	REAL*4 PIR2PV(NMCX),FC(NMCX)
	REAL*4 C(3),PEAK(3),CC(3),PPK(3),POLE(3)
	INTEGER WGT(NMCX)
	CHARACTER*80 FN
	PARAMETER (NDH=123,NDF=NDH*10)
	PARAMETER (DTR = 57.2957795, PI=3.14159265)
	INTEGER DHIST(NDH),DFINE(NDF)
	INTEGER NL,LIST(9999)
	REAL*4 WW(9999)
	PARAMETER (IRES=3,NPIX=12*4**IRES)
	PARAMETER (JRES=6,MPIX=12*4**JRES)
	REAL*4 BIGSKY(MPIX)
	REAL*4 SKY(NPIX)
	REAL*4 ABAR(NPIX)
	REAL*4 DBAR(NPIX)
	REAL*4 XBAR(NPIX)
	REAL*8 SW,SWD,SWDD,DIA,SDIA,SWA,SWAA,SWT,SWTT
	REAL*4 POLEBAR(3),MOIPOLE(3,3),EIGEN(3),VECTORS(3,3)
	REAL*4 INERTIA(NMCX),JMED,JM,JP
	REAL*4 PERIOD(NMCX),PMED,PM,PP
C
	WRITE (*,FMT='(A,$)') 'Input mcmc file:'
	READ (*,FMT='(A)') FN
	NC = LNBLNK(FN)
	OPEN (UNIT=21,FILE=FN(1:NC))
	DO I=1,NSKIP
	  READ (21,*)
	ENDDO
C	to isolate a region, put center RA,DEC below
	CALL RADEC_TO_C(160./DTR,-26./DTR,PPK)
C	and region radius below
	DVSMAX = (2*SIN(180./2./DTR))**2
	NMC = 1
	SWGTT = 0
	SWGT = 0
10	READ (21,FMT='(F8.5,7F9.5,F17.5,I5)',ERR=10,END=20) 
	1	(P(I,NMC),I=1,NP),CHSQ,WGT(NMC)
	SWGTT = SWGTT+WGT(NMC)
	CALL RADEC_TO_C(P(1,NMC),P(2,NMC),POLE)
	IF (DVS(PPK,POLE).GT.DVSMAX) GO TO 10
	SWGT = SWGT+WGT(NMC)
	DO I=3,NP
	  IF (ABS(P(I,NMC)).GT.23.) GO TO 10
	ENDDO
	NMC = NMC+1
	IF (NMC.LE.NMCX) GO TO 10
20	CLOSE (UNIT=21)
	NMC = NMC-1
C
	DO I=1,NDH
	  DHIST(I)=0
	ENDDO
	DO I=1,NDF
	  DFINE(I)=0
	ENDDO
	DO I=1,MPIX
	  BIGSKY(I) = 0
	ENDDO
	DO I=1,NPIX
	  SKY(I) = 0
	  ABAR(I) = 0
	  DBAR(I) = 0
	ENDDO
	SW = 0
	SWD = 0
	SWDD = 0
	SWA = 0
	SWAA = 0
	SWT = 0
	SWTT = 0
	DO J=1,3
	  DO I=1,3
	    MOIPOLE(I,J)=0
	  ENDDO
	  POLEBAR(J)=0
	ENDDO
C
	RADIUS = 7./DTR
	RR = (2*SIN(RADIUS/2))**2
C
	WRITE (*,*) 'entering I=1..NMC=',NMC,' loop'
	WRITE (*,*) 'patch & total weights:', SWGT,SWGTT
	DO J=1,NMC
	  D(J) = EXP(P(6,J))
	  SW = SW+WGT(J)
	  SWD = SWD+WGT(J)*D(J)
	  SWDD = SWDD+WGT(J)*D(J)**2
	  ID = 1.0+23.*ALOG10(D(J)*1000.)
	  IF (ID.GE.1 .AND. ID.LE.NDH) DHIST(ID)=DHIST(ID)+WGT(J)
	  ID = 1.0+230.*ALOG10(D(J)*1000.)
	  IF (ID.GE.1 .AND. ID.LE.NDF) DFINE(ID)=DFINE(ID)+WGT(J)
	  CALL RADEC_TO_C(P(1,J),P(2,J),C)
	  DO JJ=1,3
	    DO II=1,3
	      MOIPOLE(II,JJ) = MOIPOLE(II,JJ)+C(II)*C(JJ)*WGT(J)
	    ENDDO
	    POLEBAR(JJ) = POLEBAR(JJ) + C(JJ)*WGT(J)
	  ENDDO
	  CALL NESTED_PIXNO(C,IRES,IP)
	  SKY(IP+1) = SKY(IP+1)+WGT(J)
	  CALL NESTED_PIXLIST(C,RADIUS,JRES,NL,LIST)
	  SWW = 0
	  DO K=1,NL
	    CALL NESTED_CENPIX(LIST(K),JRES,CC)
	    X = DVS(CC,C)/RR
	    WW(K) = 0
	    IF (X.LT.0.98) WW(K) = EXP(-X/(1-X))
	    SWW = SWW+WW(K)
	  ENDDO
	  DO K=1,NL
	    BIGSKY(LIST(K)+1) = BIGSKY(LIST(K)+1)+64*WGT(J)*WW(K)/SWW
	  ENDDO
	  PIR2PV(J) = EXP(P(8,J))
	  FC(J) = 1/(1+EXP(-P(7,J)))
	  ALBEDO = EXP(P(3,J))
	  ALBEDO = ALBEDO*PVMAX/(ALBEDO+PVMAX)
	  PV(J) = ALBEDO
	  PERIOD(J) = EXP(MIN(11.5,P(4,J)))
	  INERTIA(J) = EXP(P(5,J))
	  OMEGA = 2*PI/(3600*EXP(P(4,J)))
	  T1 = 401 
	  EM = 0.95
          TH(J) = INERTIA(J)*SQRT(OMEGA)/(EM*5.67E-8*T1**3)
	  SWT = SWT+WGT(J)*TH(J)
	  SWTT = SWTT+WGT(J)*TH(J)**2
	  SWA = SWA+WGT(J)*ALBEDO
	  SWAA = SWAA+WGT(J)*ALBEDO**2
	  ABAR(IP+1) = ABAR(IP+1)+ALBEDO*WGT(J)
	  DBAR(IP+1) = DBAR(IP+1)+D(J)*WGT(J)
	  WRITE (3,FMT='(F10.4,1PE12.4,2H P)') D(J),PERIOD(J)
	  WRITE (4,FMT='(F10.4,F7.4,2H P)') D(J),ALBEDO
	ENDDO
C
	WRITE (*,*) 'out of 1..NMC loop'
	CLOSE (UNIT=3)
	CLOSE (UNIT=4)
C
	DO J=1,3
	  DO I=1,3
	    MOIPOLE(I,J) = MOIPOLE(I,J)/SW
	  ENDDO
	  POLEBAR(J) = POLEBAR(J)/SW
	ENDDO
	OPEN (UNIT=1,FILE='big-WISE-rc-MCMC.map',FORM='unformatted')
	WRITE (1) BIGSKY
	CLOSE (UNIT=1)
	OPEN (UNIT=1,FILE='pole-WISE-rc-MCMC.map',FORM='unformatted')
	WRITE (1) SKY
	CLOSE (UNIT=1)
	WRITE (2,FMT='(13I6)') DHIST
	WRITE (32,FMT='(13I6)') DFINE
	SKYMAX = 0
	DO IP=1,NPIX
	  IF (SKY(IP).GT.SKYMAX) THEN
	    CALL NESTED_CENPIX(IP-1,IRES,PEAK)
	    SKYMAX = SKY(IP)
	  ENDIF
	  IF (SKY(IP).GT.10.) THEN
	    ABAR(IP) = ABAR(IP)/SKY(IP)
	    DBAR(IP) = DBAR(IP)/SKY(IP)
	    XBAR(IP) = XBAR(IP)/SKY(IP)
	  ELSE
	    ABAR(IP) = -1.E30
	    DBAR(IP) = -1.E30
	    XBAR(IP) = -1.E30
	  ENDIF
	ENDDO
	OPEN (UNIT=1,FILE='albedo-WISE-rc-MCMC.map',FORM='unformatted')
	WRITE (1) ABAR
	CLOSE (UNIT=1)
	OPEN (UNIT=1,FILE='diameter-WISE-rc-MCMC.map',FORM='unformatted')
	WRITE (1) DBAR
	CLOSE (UNIT=1)
C
	DIA = SWD/SW
	SDIA = SQRT(SWDD/SW-DIA**2)
	CALL MED1SIGMA(NMC,D,WGT,DMED,DM,DP)
	SDP = 100*ALOG(DP/DMED)
	SDM = 100*ALOG(DMED/DM)
C
	ALB = SWA/SW
	SALB = SQRT(SWAA/SW-ALB**2)
	CALL MED1SIGMA(NMC,PV,WGT,PVMED,PVM,PVP)
	SPV = 50*ALOG(PVP/PVM)
C
	THETA1 = SWT/SW
	STH = SQRT(SWTT/SW-THETA1**2)
	CALL MED1SIGMA(NMC,TH,WGT,TH1MED,TH1M,TH1P)
	STH1 = 50*ALOG(TH1P/TH1M)
	STHP = 100*ALOG(TH1P/TH1MED)
	STHM = 100*ALOG(TH1MED/TH1M)
C
	WRITE (*,9001) 'dia=',DIA,SDIA,' median',DMED,SDP,SDM
900	FORMAT(A,F9.4,A,F9.4,A,F9.4,A,F5.1,1H%)
9001	FORMAT(A,F9.4,3H+/-,F9.4,A,F9.4,1H+,F6.1,1H-,F6.1,1H%)
	WRITE (*,900) 
	1	'p_V = ',ALB,' +/-',SALB,' median',PVMED,' +/- ',SPV
	WRITE (*,9001) 'theta1=',THETA1,STH,' median',TH1MED,STHP,STHM
901	FORMAT(A,1PE11.3,A,1PE11.3,A,0PF7.3,A,F7.3,A)
	CALL MED1SIGMA(NMC,PERIOD,WGT,PMED,PM,PP)
	SPP = 100*ALOG(PP/PMED)
	SPM = 100*ALOG(PM/PMED)
	WRITE (*,FMT='(A,F8.2,1H+,F5.1,F6.1,1H%)') 'Period [h] =',
	1	PMED,SPP,SPM
	CALL MED1SIGMA(NMC,INERTIA,WGT,JMED,JM,JP)
	SJP=100*ALOG(JP/JMED)
	SJM=100*ALOG(JM/JMED)
	WRITE (*,FMT='(A,F8.1,1H+,F5.1,F6.1,1H%)') 'sqrt(kappa*rho*C)=',
	1	JMED,SJP,SJM
	CALL MED1SIGMA(NMC,FC,WGT,FCMED,FCM,FCP)
	SFCP = FCP-FCMED
	SFCM = FCM-FCMED
	WRITE (*,FMT='(A,F6.3,1H+,F5.3,F6.3)') 
	1	'crater fraction=',FCMED,SFCP,SFCM
	CALL MED1SIGMA(NMC,PIR2PV,WGT,CMED,CM,CP)
	SCP = 100*ALOG(CP/CMED)
	SCM = 100*ALOG(CM/CMED)
	WRITE (*,FMT='(A,F6.3,1H+,F4.1,F5.1,1H%)') 
	1	'p_IR/p_V=',CMED,SCP,SCM
	CALL C_TO_RADEC(PEAK,RA,DEC)
	RA = DTR*RA
	DEC = DTR*DEC
	WRITE (*,*) 'pole peak at =',RA,DEC
	CALL C_TO_RADEC(POLEBAR,RA,DEC)
	RA = DTR*RA
	DEC = DTR*DEC
	PB = SQRT(DOT(POLEBAR,POLEBAR))
	WRITE (*,*) 'mean pole at =',RA,DEC,' |<p>|=',PB
	CALL JACOBI(MOIPOLE,3,3,EIGEN,VECTORS,NROT)
	DO J=1,3
	  CALL C_TO_RADEC(VECTORS(J,1),RA,DEC)
	  RA = DTR*RA
	  DEC = DTR*DEC
	  WRITE (*,FMT='(A,F9.6,A,2F7.1)') 'Moment eigenvalue=',EIGEN(J)
	1	,' at RA,Dec=',RA,DEC
	ENDDO
	STOP
	END

	SUBROUTINE MED1SIGMA(NMC,X,W,XMED,XM,XP)
	REAL X(*)
	INTEGER W(*)
	PARAMETER (NMCX=1500000)
	COMPLEX*8 C(NMCX)
	REAL V(3)
	SW = 0
	DO J=1,NMC
	  SW = SW+W(J)
	  C(J) = CMPLX(X(J),FLOAT(W(J)))
	ENDDO
	CALL CSORT(NMC,C)
	DO I=1,3
	  T = SW*(0.50+0.34*(I-2))
	  S = 0
	  DO J=1,NMC
	    S = S+AIMAG(C(J))
	    IF (S.GT.T) THEN
	      F = (T-(S-AIMAG(C(J))))/AIMAG(C(J))
	      V(I) = F*C(J)+(1-F)*C(J-1)
	      GO TO 100
	    ENDIF
	  ENDDO
100	ENDDO
	XMED = V(2)
	XM = V(1)
	XP = V(3)
	RETURN
	END
