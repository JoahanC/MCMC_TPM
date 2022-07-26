	LOGICAL FUNCTION SAA(T68)
C
C	returns true if position is within 1E-3 pixels hit contour
C	for WISE
C
	IMPLICIT NONE
	REAL*8 T68	! input, time in seconds since 5/24/68
C	table read from map
	REAL*4 MAXLONG,MINLONG
	PARAMETER (MINLONG= -90.,MAXLONG=10.)
	INTEGER NPT
	PARAMETER (NPT= 11)
	REAL*4 SAA_NL(NPT) /-20.,-10.,-7.,-4.,-1.,0.,-2.,-10.,-18.,
	1  -22.,-26./
	REAL*4 SAA_SL(NPT) /-10.,-30.,-42.,2*-46.,-45.,-42.,-38.,-36.,
	1  -32.,-30./
	REAL*4 R(3)
	REAL*4 LAT,LONG
	INTEGER J,N
	REAL*4 X,Y
	REAL*8 DTR
	PARAMETER (DTR=57.2957795130824D0)
C	get satellite position in geocentric km
	CALL TERRESTIAL(T68,R)
	LONG = DTR*ATAN2(R(2),R(1))
	SAA = .FALSE.
	IF (LONG.LT.MINLONG .OR. LONG.GT.MAXLONG) RETURN
	LAT = DTR*ATAN2(R(3),SQRT(R(1)**2+R(2)**2))
C	if latitude not in range return false
C	else - in the SAA, return true
	X = (NPT-1)*(LONG-MINLONG)/(MAXLONG-MINLONG)
	J = X
	X = X-J
	J = J+1
	IF (J.LE.0 .OR. J.GE.NPT) RETURN
	SAA = (LAT.LT.(1-X)*SAA_NL(J)+X*SAA_NL(J+1)) .AND.
	1	(LAT.GT.(1-X)*SAA_SL(J)+X*SAA_SL(J+1))
	RETURN
	END
