	REAL*4 FUNCTION DMS(R)
C
C	CONVERTS INPUT REAL*4 R FROM "DD.MMSS" INTO RADIANS
C
	REAL*4 R
	PARAMETER (PI=3.1415927,TWOPI=6.2831853,DTR=57.29577951)
	V=ABS(R)
	ID=V
	V=100.*(V-FLOAT(ID))
	IM=V
	V=100.*(V-FLOAT(IM))
	V=FLOAT(ID)+FLOAT(IM)/60.+V/3600.
	IF(R.LT.0.) V=-V
	DMS=V/DTR
	RETURN
	END
