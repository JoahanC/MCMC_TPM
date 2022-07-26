	REAL*8 FUNCTION TDN(T)
C
C	routine to find the next descending node time given T
C	T and TDN are seconds since 5/24/68 [JD=2440000.5]
C
C	ORBIT (T,XYZ) returns vector XYZ in km relative to Earth's center
C	in current of date celestial corrdinates
C
C	but this routine converts into ecliptic coordinates to give the
C	time the satellite descends through the ecliptic
C
	REAL*8 T,T1(2),TB,TT
	REAL*4 XYZ(3,3),C(3)
	T1(1) = T
	T1(2) = T+400
	DO I=1,25
	  DO J=1,2
	    CALL ORBIT(T1(J),C)
	    CALL CTOE(C,XYZ(1,J))
	  ENDDO
	  IF (XYZ(3,1).GE.0. .AND. XYZ(3,2).LE.0) GO TO 100
	  T1(1) = T1(2)
	  T1(2) = T1(1)+400
	ENDDO
	STOP 'DN not found'
100	E = 100.
	DO WHILE (E.GT.0.01)
C	if XYZ(3,1)=0 return T1.  If XYZ(3,2)=0 return T2
	  TT = T1(1) - (T1(2)-T1(1))*XYZ(3,1)/(XYZ(3,2)-XYZ(3,1))
	  CALL ORBIT(TT,C)
	  CALL CTOE(C,XYZ(1,3))
C	if South of Equator replace T2 and XYZ(,2), otherwise replace T1...
	  J = 2
	  IF (XYZ(3,1).GE.0.) J=1
	  T1(J) = TT
	  DO I=1,3
	    XYZ(I,J) = XYZ(I,3)
	  ENDDO
	  E = ABS(XYZ(3,1))
	  TB = T1(1)
	  IF (ABS(XYZ(3,2)).LT.E) THEN
	    E = ABS(XYZ(3,2))
	    TB = T1(2)
	  ENDIF
	ENDDO
	TDN = TB
	RETURN
	END
