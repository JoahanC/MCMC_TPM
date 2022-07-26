	SUBROUTINE NESTED_EDGES(IPIX,IRES,NP,E)
	INTEGER IPIX	! nested scheme pixel number, input, 0..(12*4**IRES -1)
	INTEGER IRES	! resolution, range 0..13
C
C	convert "nested" scheme pixel number into unit vectors
C	in array E which trace the edges of the pixel. NP points
C	per edge. E must have second dimension GE 4*NP
C
	REAL*4 E(3,*)
	INTEGER*2 IX,IY
	COMMON /SPRCNPX/ IX(0:1023),IY(0:1023)
	REAL*8 TWOPI,PI,PIOVER2,PIOVER4
	PARAMETER (TWOPI = 6.28318530717958D0)
	PARAMETER (PI = 3.14159265358979D0)
	PARAMETER (PIOVER2 = 1.5707963267949D0)
	PARAMETER (PIOVER4 = 0.785398163397448D0)
C	lower left corner X coordinates in units of pi/4
	INTEGER IXLLC(12)/1,3,5,7,0,2,4,6,1,3,5,7/
C	lower left corner Y coordinates in units of 2/3
	INTEGER IYLLC(12)/0,0,0,0,-1,-1,-1,-1,-2,-2,-2,-2/
C	central longitude index for the triangular pieces
	INTEGER NOFFACE(12)/0,1,2,3,-1,-1,-1,-1,0,1,2,3/
	INTEGER IDX(4) /0,1,1,0/
	INTEGER IDY(4) /0,0,1,1/
C	table of powers of 2
	INTEGER*4 IT14(14)/1,2,4,8,16,32,64,128,256,512,1024,2048,4096,
	1	8192/
C	table of powers of 4
	INTEGER*4 IT28(14)/1,4,16,64,256,1024,4096,16384,65536,262144,
	1	1048576,4194304,16777216,67108864/
	IF (IY(1023).NE.31) CALL SPRCNPX_SET
	IF (IRES.LT.0 .OR. IRES.GT.13) STOP 'IRES out of range in EDGES'
	NPIX = 12*IT28(IRES+1)
	IF (IPIX.LT.0 .OR. IPIX.GE.NPIX) STOP 'IPIX out of range in EDGES'
	NFACE = IPIX/IT28(IRES+1)
	N=MOD(IPIX,IT28(IRES+1))
	I=MOD(N,1024)
	N=N/1024
	J=MOD(N,1024)
	K=N/1024
	JX=1024*IX(K)+32*IX(J)+IX(I)
	JY=1024*IY(K)+32*IY(J)+IY(I)
	XO = FLOAT(JX+IDX(1))/IT14(IRES+1)
	YO = FLOAT(JY+IDY(1))/IT14(IRES+1)
	JP = 0
	DO JC=1,4
	  IC = MOD(JC,4)+1
	  XN = FLOAT(JX+IDX(IC))/IT14(IRES+1)
	  YN = FLOAT(JY+IDY(IC))/IT14(IRES+1)
	  DO IP=1,NP
	    F = (IP-1.)/NP
	    XP = (1.-F)*XO + F*XN
	    YP = (1.-F)*YO + F*YN
	    X = PIOVER4*(XP-YP)
	    Y = (XP+YP)/1.5
	    YY = Y + IYLLC(NFACE+1)/1.5
	    XX = X + IXLLC(NFACE+1)*PIOVER4
	    IF (XX.LT.0.) XX = XX+TWOPI
	    N = NOFFACE(NFACE+1)
	    YA = ABS(YY)
	    IF (YA.GT.(2./3.)) THEN
	      SINLAT = 2*YA-0.75*YA*YA-1./3.
	      COSLAT = SQRT((1.+SINLAT)/12.)*(4-3.*YA)
	      SINLAT = SIGN(SINLAT,YY)
	      IF (N.GE.0) THEN
	        DIV = 2.-1.5*YA
	        IF (DIV.GT.0) THEN
	          ALONG = (XX-(N+0.5)*PIOVER2*(1.5*YA-1.))/(2.-1.5*YA)
	        ELSE
	          ALONG = (N+0.5)*PIOVER2
	        ENDIF
	      ELSE
	        ALONG = XX
	      ENDIF
	    ELSE
	      SINLAT = YY
	      COSLAT = SQRT(1.-SINLAT**2)
	      ALONG = XX
	    ENDIF
	    JP = JP+1
	    E(3,JP) = SINLAT
	    E(1,JP) = COSLAT*COS(ALONG)
	    E(2,JP) = COSLAT*SIN(ALONG)
	  ENDDO
	  XO = XN
	  YO = YN
	ENDDO
	RETURN
	END
