C	SHOWPXL.FOR

	SUBROUTINE WORLD(G,X,Y)
	REAL*4 G(3),X,Y
	COMMON /CNTRL/ ISW,ICD
	COMMON /GNOMONIC/ FOCAL
	PARAMETER (PIOVER2=1.5707963)
	IF (ISW.LT.1 .OR. ISW.GT.7) STOP 'BAD ISW'
	GO TO (10,20,30,40,50,60,70), ISW
C	AITOFF
10	CALL AITOFF(G,X,Y)
	GO TO 100
C	MOLLWEIDE
20	CALL MOLLWEIDE(G,X,Y)
	GO TO 100
C	SINUSOIDAL
30	B=ATAN2(G(3),SQRT(G(2)**2+G(1)**2))
	Y=B/PIOVER2
	X=-ATAN2(G(2),G(1))*COS(B)/PIOVER2
	GO TO 100
C	POLAR
40	DEN=SQRT(1.+ABS(G(3)))
	Y=-G(2)/DEN
	X=-1.+G(1)/DEN
	IF(G(3).LT.0.) X=-X
	GO TO 100
C	LAT/LONG
50	X=-ATAN2(G(2),G(1))/PIOVER2
	Y=ATAN2(G(3),SQRT(G(2)**2+G(1)**2))/PIOVER2
	GO TO 100
C	gnomonic
60	X = -FOCAL*G(2)/G(1)
	Y = FOCAL*G(3)/G(1)
	GO TO 100
C	stereographic
70	X = -2*FOCAL*G(2)/(1+G(1))
	Y = 2*FOCAL*G(3)/(1+G(1))
	GO TO 100
100	RETURN
	END

	SUBROUTINE AWORLD(X,Y,G)
	REAL*4 X,Y,G(3)
	COMMON /CNTRL/ ISW,ICD
	COMMON /GNOMONIC/ FOCAL
	PARAMETER (PIOVER2=1.5707963)
	IF (ISW.LT.1 .OR. ISW.GT.7) STOP 'BAD ISW'
	GO TO (10,20,30,40,50,60,70), ISW
C	AITOFF
10	CALL AAITOFF(X,Y,G)
	GO TO 100
C	MOLLWEIDE
20	CALL AMOLLWEIDE(X,Y,G)
	GO TO 100
C	SINUSOIDAL
30	G(3)=SIN(Y*PIOVER2)
	CB=COS(Y*PIOVER2)
	IF (CB.LE.0.0001) THEN
	  G(3) = 1.
	  IF (Y.LT.0.) G(3) = -1.
	  G(1) = 0.
	  G(2) = 0.
	ELSE
	  G(1)=COS(X*PIOVER2/COS(Y*PIOVER2))*CB
	  G(2)=-SIN(X*PIOVER2/COS(Y*PIOVER2))*CB
	ENDIF
	GO TO 100
C	POLAR
40	G(3)=1.-(Y*Y+(1.-ABS(X))**2)
	DEN = 1.+G(3)
	IF (DEN.GT.0.) THEN
	  DEN=SQRT(1.+G(3))
	  G(2)=-DEN*Y
	  G(1)=DEN*(1.-ABS(X))
	ELSE
	  G(1) = 0.
	  G(2) = 0.
	ENDIF
	IF (X.GT.0.) G(3)=-G(3)
	GO TO 100
C	LAT/LONG
50	G(3)=SIN(Y*PIOVER2)
	CB=COS(Y*PIOVER2)
	G(1)=COS(X*PIOVER2)*CB
	G(2)=-SIN(X*PIOVER2)*CB
	GO TO 100
C	gnomonic
60	G(1) = FOCAL
	G(2) = -X
	G(3) = Y
	CALL NORM(G)
	GO TO 100
C	stereographic
C	r = s/(1+c) = s(1-c)/(1-c^2) = (1-c)/s
C	r^2 = (1-c)/(1+c)
C	c = (1-r^2)/(1+r^2)
C	for FOCAL=1
C	Y=2, RR=1, G(1)=0, G(3)=1
C	Y=1, X=0 => RR=0.25, G(1) = 0.6, G(3) = 0.8
C	Y=0.1, X=0 => RR=0.0025, G(1) = 0.995 , G(3) = 0.1
70	G(2) = -X/(2*FOCAL)
	G(3) = Y/(2*FOCAL)
	RR = G(2)**2+G(3)**2
	G(1) = (1-RR)/(1+RR)
	G(2) = (1+G(1))*G(2)
	G(3) = (1+G(1))*G(3)
	CALL NORM(G)
100	RETURN
	END

	SUBROUTINE AITOFF(GAL,X,Y)
C	CONVERTS DIRECTION COSINES IN GAL INTO X,Y POSITION
C	ON AITOFF EQUAL AREA MAP.  GAL(3)=SIN(B)
C	GAL(1)=COS(L)*COS(B), GAL(2)=SIN(L)*COS(B)
C	-1<Y<1, -2<X<2
	DIMENSION GAL(3)
	  R=SQRT(GAL(1)**2+GAL(2)**2)
C	  CHECK FOR POLES AND EXIT IF AT POLE
	  IF(R.GT.1.E-5) GO TO 10
	    X=0.
	    Y=GAL(3)
	    RETURN
C	CUT LONGITUDE IN HALF
10	  C=GAL(1)/R
	  S=GAL(2)/R
	  IF(C.GT.0.0) GO TO 20
	    S2=SQRT((1.-C)/2.)
	    IF(S.LT.0.) S2=-S2
	    C2=S/2./S2
	    GO TO 30
20	    C2=SQRT((1.+C)/2.)
	    S2=S/2./C2
30	C2=R*C2
	S2=R*S2
C
C	THE UNIT VECTOR WITH HALF LONGITUDE IS NOW (C2,S2,GAL(3))
C
	DEN=SQRT(1.+C2)
	Y=GAL(3)/DEN
	X=-2.*S2/DEN
	RETURN
	END

	SUBROUTINE AAITOFF(X,Y,GAL)
C	INVERTS AITOFF EQUAL AREA MAP
C	-2<X<2, -1<Y<1 GIVES DIRECTION COSINES IN GAL
	DIMENSION GAL(3)
	DATA R2/1.414214/
	YY=Y
	XX=X/2.
	R=SQRT(XX**2+YY**2)
C	R=SIN(THETA/2)
	IF (R.LT.0.01) THEN
C	  EXPANSION FOR REGION NEAR THE CENTER
	  GAL(3)=R2*YY
	  GAL(1)=1.0
	  GAL(2)=-R2*X
	  CALL NORM(GAL)
	  RETURN
	ENDIF
	ST2=R/R2
	ONEMRR=1.-ST2**2
	CT2=SQRT(AMAX1(0.,ONEMRR))
	ST=2.*ST2*CT2
	CT=CT2**2-ST2**2
	GAL(3)=ST*YY/R
	S=-ST*XX/R
	R=SQRT(CT**2+S**2)
	C=CT/R
	S=S/R
	GAL(1)=(C*C-S*S)*R
	GAL(2)=2.*C*S*R
	RETURN
	END

	SUBROUTINE MOLLWEIDE(GAL,X,Y)
C
C	MOLLWEIDE PROJECTION - EQUAL AREA INTO ELLIPSE
C	-2 LE X LE 2, -1 LE Y LE 1
C
C	output has EAST to the left (Negative X) for SKYMAPS
C
C	UNLIKE AITOFF, HAS HORIZONTAL LATITUDE LINES
C
C	THE RULE IS:
C	Given t such that: 2t + sin(2t) = pi*sin(lat)
C		then	   Y = sin(t)
C			   X = -2*long*cos(t)/pi
C
	REAL*4 GAL(3),LONG,X,Y
	PARAMETER (PI=3.1415927,TWOPI=6.2831853,DTR=57.29578)
	ASL=ABS(GAL(3))
	IF(ASL.LT..707) THEN
	  CTHETA=ASL
	  F=PI*CTHETA/2.0
	  FF=F*F
	  T=F*(0.5+FF*(0.04166666667+FF*(.008333333333+FF*.0021329365)))
	  E=T+0.5*SIN(2.0*T)-F
	  DEDT=1.0+COS(2.0*T)
	  T=T-E/DEDT
	  ST=SIN(T)
	  CT=COS(T)
	ELSE
C	  SOLVE X-SIN(X)=2*PI*SIN(THETA/2)**2
C	  THETA = PI/2 - ALAT
	  RHS=0.5*PI*((1.-ASL)**2+GAL(1)**2+GAL(2)**2)
	  X=(6.*RHS)**0.33333333
	  DX=1.
	  DO WHILE (DX.GT.1.E-5)
	    XX=X*X
	    IF (X.LT.0.8) THEN
	      F=X**5*(0.0083333333-XX*0.0001984127)
	    ELSE
	      F=SIN(X)-X*(1.-.166666667*XX)
	    ENDIF
	    XX=(6.*(RHS+F))**0.33333333
	    DX=ABS(XX-X)
	    X=XX
	  ENDDO
C	  T=PI/2-X/2
	  ST=COS(0.5*X)
	  CT=SIN(0.5*X)
	ENDIF
	IF(GAL(3).LT.0.) ST=-ST
	Y=ST
	LONG=0.
	IF(GAL(1).NE.0. .OR. GAL(2).NE.0.) LONG=ATAN2(GAL(2),GAL(1))
	X=-2.*LONG*CT/PI
	RETURN
	END

	SUBROUTINE AMOLLWEIDE(X,Y,GAL)
	REAL*4 X,Y,GAL(3)
C
C	INVERTS THE MOLLWEIDE PROJECTION WITH -2 LE X LE 2
C		-1 LE Y LE 1 -> UNIT VECTOR IN GAL
C
C	ASSUMES EAST IS TO THE LEFT FOR SKYMAPS
C
C	X,Y = (+/-2,0) -> (-1,0,0)
C	    = (0,0)    -> (1,0,0)
C	    = (0,+/-1) -> (0,0,+/-1)
C	    = (+/-1,0) -> (0,-/+1,0)
C
	REAL*4 LONG
	PARAMETER (PI=3.1415927,TWOPI=6.2831853,DTR=57.29578)
	ST=Y
	CT=SQRT(MAX(0.,1.-ST**2))
	LONG=0.
	IF (CT.GT.0.) LONG=-0.5*PI*X/CT
	T=ATAN2(ST,CT)
	GAL(3)=(2.*T+SIN(2.*T))/PI
	CB=SQRT(MAX(0.,1.-GAL(3)**2))
	GAL(1)=COS(LONG)*CB
	GAL(2)=SIN(LONG)*CB
	RETURN
	END
              
