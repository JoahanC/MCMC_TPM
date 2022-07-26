	SUBROUTINE ROBUST(NP,X,SG,MU,SIGMA,RCS)
C
C	NP input points in X with sigma's in SG
C	return estimated centroid in MU and sigma of MU in SIGMA
C	with robust chi^2 in RCS
C	deviations larger than EMX sigma are Winsorized
C
	PARAMETER (EMX = 2.0)
	REAL*8 X(NP),SG(NP),MU,SIGMA,RCS
	REAL*8 W,SW,SWX,CHI2
	PARAMETER (NPMX=500)
	REAL*4 DX(2*NPMX)
	REAL*4 Y(2*NPMX)
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
C	make sorted list of deviations
	DO I=1,NP
	  DX(I) = X(I)-XC
	ENDDO
	CALL RSORT(NP,DX)
C	This has a forced assumption that there
C	are fewer than 100 data points
	NX = NPMX/NP
	NX = MAX(5,NX)
	N = NX*(NP-1)+1
C	WRITE(*,FMT='(10F8.4)') (DX(I),I=1,NP)
C	make table of values to try looking for minimum

C	This is now using the initialized array and building
C	sub points from back to front, so by the time you
C	overwrite the original data, it doesn't matter
	
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
C	calculate robust chi^2 for each value in DX
	DO I=1,N
	  Y(I)=0
	  T = DX(I)
	  DO J=1,NP
	    E = ABS((X(J)-XC-T)/SG(J))
	    IF (E.GT.EMX) THEN
	      CS = 2*EMX*E-EMX**2
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
	      E = ABS((X(J)-XC-T)/SG(J))
	      IF (E.GT.EMX) THEN
	        CS = 2*EMX*E-EMX**2
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
	  E = ABS((X(J)-MU)/SG(J))
	  IF (E.GT.EMX) THEN
	    CS = 2*EMX*E-EMX**2
	  ELSE
	    CS = E*E
	  ENDIF
	  RCS = RCS+CS
	ENDDO
	RETURN
	END
