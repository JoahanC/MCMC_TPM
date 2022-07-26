	REAL*4 FUNCTION DELTA_T(T68)
C
C	T68 is REAL*8 seconds since 24 May 1968, UT 00:00:00
C	which was JD 2440000.5
C
C	DELTA_T is the difference between IAT (or ET or TDT) and UTC
C	since 5/24/68
C
C	the value returned includes the 32.18 seconds between terrestrial
C	time and TAI
C
C	The seconds to be used are ephemeris seconds, or IAT seconds
C	Therefore the conversion from UT date and time must allow for the
C	rates adjustments used before 1972 and the leap seconds used later
C	This routine first calculates UTC ignoring these corrections, then
C	looks in a table to interpolate for DT = ET - UT
C	Linear interpolation is used, so two entries are made for each
C	leap second.  Entry for JDATE should be the JD at noon the day 
C	after the leap second.
C
	IMPLICIT NONE
	INTEGER*4 NPT
	PARAMETER (NPT=34)
	INTEGER*4 JDATE(NPT),I,NL,NH,NTRY
	REAL*8 DT(NPT),DTNOW,SLOPE
	LOGICAL SAVED/.FALSE./
	INTEGER*4 JD68/2440001/
	REAL*8 DT68/38.66/
	REAL*8 T68,TAB68(NPT)
C
	SAVE SAVED,TAB68
C
	DATA (JDATE(I),DT(I),I=1,NPT)/
	1		2440001,	38.66,	      	!	5/24/1968
	1		2445517,	53.18,		!	7/1/1983
	1		2445517,	54.18,		!	7/1/1983
	1		2446248,	54.18,		!	7/1/1985
	1		2446248,	55.18,		!	7/1/1985
	1		2447162,	55.18,		!       1/1/1988
	1		2447162,	56.18,		!	1/1/1988
	1		2447893,	56.18,		!	1/1/1990
	1		2447893,	57.18,		!	1/1/1990
	1		2448258,	57.18,		!	1/1/1991
	1		2448258,	58.18,		!	1/1/1991
	1		2448805,	58.18,		!	7/1/1992
	1		2448805,	59.18,		!	7/1/1992
	1		2449170,	59.18,		!     	7/1/1993
	1		2449170,	60.18,		!     	7/1/1993
	1		2449535,	60.18,		!	7/1/1994
	1		2449535,	61.18,		!	7/1/1994
	1		2450084,	61.18,		!	1/1/1996
	1		2450084,	62.18,		!	1/1/1996
	1		2450631,	62.18,		!	7/1/1997
	1		2450631,	63.18,		!	7/1/1997
	1		2451180,	63.18,		!	1/1/1999
	1		2451180,	64.18,		!	1/1/1999
	1		2453737,	64.18,		!	1/1/2006
	1		2453737,	65.18,		!	1/1/2006
	1		2454833,	65.18,		!	1/1/2009
	1		2454833,	66.18,		!	1/1/2009
	1		2456110,	66.18,		!	7/1/2012
	1		2456110,	67.18,		!	7/1/2012
	1		2457205,	67.18,		!	7/1/2015
	1		2457205,	68.18,		!	7/1/2015
	1		2457755,	68.18,		!	1/1/2017
	1		2457755,	69.18,		!	1/1/2017
	1		2488070,	69.18/		!	1/1/2100
C
C---------------------------------------------------
	IF (.NOT. SAVED) THEN
C				 YYDDDHHMMSSFFF
	  SAVED = .TRUE.
	  DO I=1,NPT
	    TAB68(I)=(JDATE(I)-JD68)*864.D2 + (DT(I)-DT68)
	  ENDDO
	ENDIF
C
C	Find position in Table
C
	NL=1
	NH=NPT-1
	DO WHILE (NL.LT.NH)
	  NTRY=(NL+NH+1)/2
	  IF(TAB68(NTRY).GT.T68) THEN
	    NH=NTRY-1
	  ELSE
	    NL=NTRY
	  ENDIF
	ENDDO
	SLOPE = (DT(NL+1)-DT(NL))/(TAB68(NL+1)-TAB68(NL))
	DTNOW=DT(NL)+ (T68-TAB68(NL)) * SLOPE
	DELTA_T = (DTNOW-DT68)	! time correction ET-UTC in seconds
	RETURN
	END
