	SUBROUTINE SPITZER(T68,XYZ)
C	Spitzer-ephemeris.f
C	from T68 time in seconds since 24 May 1968 returns
C	Spitzer position in J2000 celestial XYZ
C	for earlier than 1-OCT-2003 returns the position of the Earth
	REAL*8 T68
	REAL*8 T01OCT2003 /1115683225.52D0/
	REAL*8 JULIAN_DAY,JDIN,JD1,JD2
	REAL*4 XYZ(3),XYZ1(3),XYZ2(3),V1(3),V2(3)
	REAL*4 E(3)
C
	IF (T68.LT.T01OCT2003) THEN
	  CALL EPHEM(-1,T68,E)
	  DO I=1,3
	    XYZ(I) = -E(I)
	  ENDDO
	  RETURN
	ENDIF
C
	OPEN (UNIT=22,
	1  FILE='/Users/wright/missions/WISE/Spitzer-ephemeris.txt')
	READ (22,*)
	JDIN = JULIAN_DAY(T68)
	READ (22,*,END=666) JD1
	READ (22,*) XYZ1
	READ (22,*) V1
	READ (22,*)
10	READ (22,*,END=666) JD2
	READ (22,*) XYZ2
	READ (22,*) V2
	READ (22,*)
C	WRITE (*,FMT='(F11.0,6F11.5)') JD2,XYZ2,V2
	IF (JDIN.LT.JD2) THEN
	  X = (JDIN-JD1)/(JD2-JD1)
	  W = (1-X*X*(3-2*X))
	  DO I=1,3
	    E(I) = (XYZ1(I)+(JDIN-JD1)*V1(I))*W
	1	+    (XYZ2(I)+(JDIN-JD2)*V2(I))*(1-W)
	  ENDDO
	  CLOSE (UNIT=22)
	  CALL ETOC(E,XYZ)
	  RETURN 
	ENDIF
	JD1 = JD2
	DO I=1,3
	  XYZ1(I) = XYZ2(I)
	  V1(I) = V2(I)
	ENDDO
	GO TO 10
666	STOP 'T68 2 big'
	END
