	LOGICAL FUNCTION SAA(T68)
C
C	returns true if position is within the 100 /cm^2/sec contour
C	for 600 km
C
	IMPLICIT NONE
	REAL*8 T68	! input, time in seconds since 5/24/68
C	table read from map
	REAL*4 MINLONG
	PARAMETER (MINLONG= -90.)
	INTEGER NLONG
	PARAMETER (NLONG= 8)
	REAL*4 NLIM(NLONG) /-21.,-20.,-6.,-4.,-4.,-6.,-9.,-23./
	REAL*4 SLIM(NLONG) /-40.,-39.,-48.,-50.,-49.,-42.,-37.,-30./
	REAL*4 R(3)
	REAL*4 LAT,LONG
	INTEGER J
	REAL*4 X
	REAL*8 DTR
	PARAMETER (DTR=57.2957795130824)
C	get satellite position in geocentric km
	CALL TERRESTIAL(T68,R)
	LONG = DTR*ATAN2(R(2),R(1))
	SAA = .FALSE.
	IF (LONG.LT.MINLONG) RETURN
	  X = 1+(LONG-MINLONG)/15
	  J = X
C	  if longitude not in range return false
	  IF (J.LT.1 .OR. J.GE.NLONG) RETURN
C	    interpolate in table
	    X=X-J
	    LAT = DTR*ATAN2(R(3),SQRT(R(1)**2+R(2)**2))
C	    if latitude not in range return false
	    IF (LAT.GT.NLIM(J)+X*(NLIM(J+1)-NLIM(J))) RETURN
	    IF (LAT.LT.SLIM(J)+X*(SLIM(J+1)-SLIM(J))) RETURN
C	else - in the SAA, return true
	SAA = .TRUE.
	RETURN
	END
