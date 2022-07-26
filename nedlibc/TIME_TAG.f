	SUBROUTINE TIME_TAG(IY,IM,ID,IH,MINUTE,SECOND,T)
C
C	CONVERTS YEAR, MONTH, DAY, HOUR, MINUTE, SECOND INTO
C	REAL*8 TIMETAG T, SECONDS SINCE 5/24/68 00:00:00 UT
C	The time tag T is REAL*8 seconds since 24 May 1968, UT 00:00:00
C	which was JD 2440000.5
C
C	The seconds to be used are ephemris seconds, or IAT seconds
C	Therefore the conversion from UT date and time must allow for the
C	rates adjustments used before 19?? and the leap seconds used later
C	This routine first calculates UTC ignoring these corrections, then
C	looks in a table to interpolate for DT = ET - UT
C	Linear interpolation is used, so two entries are made for each
C	leap second.  Entry for JDATE should be the JD at noon the day 
C	after the leap second.
C
	REAL*8 T,DT
	REAL*4 SECOND,DELTA_T
	INTEGER*4 JD,I367,JULIAN,JDZERO
	DATA JDZERO/2440001/			!JD at 5/24/68 UT noon
	DATA I367/367/
	JD=-7*(IY+(IM+9)/12)/4-3*((IY+(IM-9)/7)/100+1)/4+275*IM/9+ID
	JULIAN=JD+I367*IY+1721029
C	lookup DT at UT noon to avoid leap seconds
	T=(JULIAN-JDZERO)*864.D2+432.D2
	DT = DELTA_T(T)
	T=(JULIAN-JDZERO)*864.D2+IH*36.D2+MINUTE*6.D1+SECOND+DT
	RETURN
	END
