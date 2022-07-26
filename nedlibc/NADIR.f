	SUBROUTINE NADIR(T,E50)
C
C	GIVES UNIT VECTOR IN DIRECTION OF APPARENT NADIR = EARTH
C	ALLOWING FOR OBLATENESS
C
C	BY NED WRIGHT 7/28/81
C
C	PRECESSION PUT IN 11/27/84
C
C	INPUT: T IS REAL*8 SECONDS SINCE 5/24/68
C	OUTPUT: E50 IS UNIT VECTOR TOWARDS CENTER OF HORIZON
C	IN J2000 CELESTIAL COORDINATES
C
	REAL*8 T,TOLDEPOCH
	REAL*4 EARTH(3),E50(3),R(3,3)
	DATA TOLDEPOCH/-9999.D10/
	SAVE
	IF (ABS(T-TOLDEPOCH).GT.3.D6) THEN
	  TOLDEPOCH=T
	  CALL GET_EPOCH(T,EPOCH)
	  CALL PREMAT(EPOCH,R)
	  DO I=1,3
	    R(I,I)=R(I,I)+1.0
	  ENDDO
	ENDIF
	CALL ORBIT(T,EARTH)
C	OBLATENESS FACTOR
	HH=EARTH(1)**2+EARTH(2)**2+EARTH(3)**2
	F=1.+0.0066*(6371.**2/HH)
	EARTH(1)=-EARTH(1)
	EARTH(2)=-EARTH(2)
	EARTH(3)=-F*EARTH(3)
	CALL NORM(EARTH)
	CALL VMAT(R,EARTH,E50)
	RETURN
	END
