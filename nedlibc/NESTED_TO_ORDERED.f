	SUBROUTINE NESTED_TO_ORDERED(IPIX,IRES,ILAT,ILONG,IHALF,NLONG,
	1	JPIX)
	INTEGER IPIX,IRES	! input
	INTEGER ILAT,ILONG,IHALF,JPIX	! output
C
C	convert nested pixel number into lat, long ordered pixel number
C	order of JPIX: 0-3 for northernmost row, 4-11 for next row...
C	N = 2**IRES is the number of pixels for side of each rhombus
C	NPIX = 12*4**IRES = 12*N**2	range of IPIX and JPIX 0..NPIX-1
C	range of IRES 0..13
C	ILAT = 0 for northernmost row to 4*N-2 for southernmost row
C	NLONG is number of longitudes in this row: range 4 to 4*N
C	ILONG = 0..NLONG-1
C	IHALF = 0 or 1
C	longitude = 360*(ILONG+0.5*IHALF)/NLONG degrees
C
	INTEGER*2 KX,KY
	COMMON /SPRCNPX/ KX(0:1023),KY(0:1023)
	REAL*8 TWOPI,PI,PIOVER2,PIOVER4
	PARAMETER (TWOPI = 6.28318530717958D0)
	PARAMETER (PI = 3.14159265358979D0)
	PARAMETER (PIOVER2 = 1.5707963267949D0)
	PARAMETER (PIOVER4 = 0.785398163397448D0)
C	lower left corner X coordinates in units of pi/4
	INTEGER IXLLC(12)/1,3,5,7,0,2,4,6,1,3,5,7/
C	lower left corner Y coordinates in units of 2/3
	INTEGER IYLLC(12)/0,0,0,0,-1,-1,-1,-1,-2,-2,-2,-2/
C	table of powers of 2
	INTEGER*4 IT14(14)/1,2,4,8,16,32,64,128,256,512,1024,2048,4096,
	1	8192/
C	table of powers of 4
	INTEGER*4 IT28(14)/1,4,16,64,256,1024,4096,16384,65536,262144,
	1	1048576,4194304,16777216,67108864/
	IF (KY(1023).NE.31) CALL SPRCNPX_SET
	IF (IRES.LT.0 .OR. IRES.GT.13) 
	1	STOP 'IRES out of range in NESTED_TO_ORDERED'
	NPIX = 12*IT28(IRES+1)
	IF (IPIX.LT.0 .OR. IPIX.GE.NPIX) STOP
	1	'IPIX out of range in IPIX_2_JPIX'
	NFACE = IPIX/IT28(IRES+1)
	N=MOD(IPIX,IT28(IRES+1))
	I=MOD(N,1024)
	N=N/1024
	J=MOD(N,1024)
	K=N/1024
	JX=1024*KX(K)+32*KX(J)+KX(I)
	JY=1024*KY(K)+32*KY(J)+KY(I)
C
	N = IT14(IRES+1)
C	range of JX and JY is 0..N-1
	IY = JX+JY
	IX = JX-JY
	IF (NFACE.LE.3) THEN
	  ILAT = 2*N-2-IY
	  IF (IY.GT.N-1) THEN
	    ILONG = (ILAT+1)*NFACE + (IX+ILAT)/2
	    IHALF = 1
	    NLONG = 4*(ILAT+1)
	  ELSE
	    ILONG = N*NFACE + (N+IX)/2
	    IHALF = MOD(ILAT,2)
	    NLONG = 4*N
	  ENDIF
	ELSE IF (NFACE.LE.7) THEN
	  NLONG = 4*N
	  ILAT = 3*N-2-IY
	  ILONG = N*(NFACE-4) + (2*N+IX)/2 - N
	  IF (ILONG.LE.0) ILONG = ILONG+4*N
	  IF (ILONG.GE.4*N) ILONG = ILONG-4*N
	  IHALF = MOD(ILAT,2)
	ELSE
	  ILAT = 4*N-2-IY
	  IF (IY.LE.N-1) THEN
	    ILONG = (IY+1)*(NFACE-8) + (IX+IY)/2
	    IHALF = 1
	    NLONG = 4*(IY+1)
	  ELSE
	    ILONG = N*(NFACE-8) + (IX+N)/2
	    IHALF = MOD(ILAT,2)
	    NLONG = 4*N
	  ENDIF
	ENDIF
C	compute JPIX
	IF (ILAT.LE.N) THEN
	  JPIX = 2*ILAT*(ILAT+1)+ILONG
	ELSE IF (ILAT.LE.3*N) THEN
	  JPIX = 4*N*(ILAT-N) + 2*N*(N+1) + ILONG
	ELSE
	  JPIX = N*(2+10*N) + 2*(ILAT-3*N)*(N-1 + (N-(ILAT-3*N))) + ILONG
	ENDIF
C
	RETURN
	END
