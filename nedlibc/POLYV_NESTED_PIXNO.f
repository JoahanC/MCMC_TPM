	SUBROUTINE POLYV_NESTED_PIXNO(EE,NE,IRES,IPIX)
C
	REAL*4 EE(4,NE),E(3)	! input unit vector
	INTEGER NE		! number of input vectors
	INTEGER IRES		! resolution
	INTEGER IPIX(NE)	! pixel number in `nested' scheme
C
C	written by Ned Wright 27-Feb-01 based on NESTED_PIXNO
C
C	written by Ned Wright 22-Jul-97 from Hinshaw's verbal description
C	of Gorski's scheme
C
C	NPIX = number of pixels = 12*4**IRES
C	range of IRES 0..13
C	range of IPIX 0..(NPIX-1)
C
C	sphere is cut into an equatorial band and four triangles on the top
C	and 4 on the bottom.  Band is cut into 4 squares and each is cut into
C	four triangles.  total 24 triangles now reassembled into 12 diamonds
C	for equal area need edge of equatorial band at sin(lat) = 2/3 or
C	lat = 41.81^o
C	y = sin(lat) up to this point
C	y_max = 4/3
C	in polar regions, y > 2/3
C	and	d(area) = dy*2*pi*(4/3 -y)*(3/2) = 2*pi*cos(lat)*d(lat)
C	thus 2*(y-2/3)-0.75*(y^2-4/9) = [sin(lat)-2/3]
C	and y = (4/3) - sqrt((4/3)*(1-sin(lat)))
C
	REAL*8 TWOPI,PI,PIOVER2,PIOVER4
	PARAMETER (TWOPI = 6.28318530717958D0)
	PARAMETER (PI = 3.14159265358979D0)
	PARAMETER (PIOVER2 = 1.5707963267949D0)
	PARAMETER (PIOVER4 = 0.785398163397448D0)
	PARAMETER (R2V3 = 0.66666666666667)
C	table of powers of 2
	INTEGER*4 IT14(14)/1,2,4,8,16,32,64,128,256,512,1024,2048,4096,
	1	8192/
C	table of powers of 4
	INTEGER*4 IT28(14)/1,4,16,64,256,1024,4096,16384,65536,262144,
	1	1048576,4194304,16777216,67108864/
C	lower left corner X coordinates in units of pi/4
	INTEGER IXLLC(12)/1,3,5,7,0,2,4,6,1,3,5,7/
C	lower left corner Y coordinates in units of 2/3
	INTEGER IYLLC(12)/0,0,0,0,-1,-1,-1,-1,-2,-2,-2,-2/
C	table of integers with bits separated made by SPRCNPX_SET
	INTEGER*4 IX,IY
	COMMON /SPRBIT/ IX(128),IY(128)
	IF (IX(128).EQ.0) CALL SPRBIT_SET
	IF (IRES.LT.0 .OR. IRES.GT.13) STOP 'IRES out of range in PIXEL'
	DO IE=1,NE
C	need E to be normalized
C	EM = EE(1)**2 + EE(2)**2 + EE(3)**2
C	IF(ABS(EM-1).GT.5.E-3) THEN
C	  F = 1./SQRT(EM)
C	ELSE
C	  F = 8./(3.+EM*(6.-EM))
C	ENDIF
C	DO I=1,3
C	  E(I) = F*EE(I)
C	ENDDO
	DO I=1,3
	  E(I) = EE(I,IE)
	ENDDO
C
	Z = ABS(E(3))
	X = 0.
	IF (E(2).NE.0. .OR. E(1).NE.0.) X = ATAN2(E(2),E(1))
	IF (X.LT.0.) X=X+TWOPI
	N = X/PIOVER2
	IF (N.GT.3) N=3
	IF (Z.LE.R2V3) THEN
	  YY = E(3)
	  XX = X
	ELSE
	  Y = (4.-SQRT(12.*(E(1)**2+E(2)**2)/(1+Z)))/3.
	  YY = SIGN(Y,E(3))
	  XX = X*(2.-1.5*Y)+(N+0.5)*PIOVER2*(1.5*Y-1.)
	ENDIF
C	face finding logic
	IF (E(3).GT.R2V3) THEN
	  NFACE = N
	ELSE IF (E(3).LT.-R2V3) THEN
	  NFACE = N+8
	ELSE
	  M = X/PIOVER2 + 0.5
	  DX = X/PIOVER2 - M
C	  IF (DX.GT.2.) DX = DX-4.
	  IF (Z.LT.(R2V3-1.33333333*ABS(DX))) THEN
	    NFACE = X/PIOVER2+4.5
	    IF (NFACE.GE.8) NFACE = NFACE-4
	  ELSE
	    NFACE = N
	    IF (E(3).LT.0.) NFACE = NFACE+8
	  ENDIF
	ENDIF
	X = XX - IXLLC(NFACE+1)*PIOVER4
	IF (X.GT.PI) X = X-TWOPI
	Y = YY - IYLLC(NFACE+1)/1.5
C	resolve the diamond into a square
	XP = (X/PIOVER4+1.5*Y)/2.0
	YP = (1.5*Y-X/PIOVER4)/2.0
C	write (*,*) xx,yy,nface
C	write (*,*) x,y
C	write (*,*) xp,yp
	I=IT14(IRES+1)*XP
	J=IT14(IRES+1)*YP
	I=MIN(IT14(IRES+1)-1,I)			!Roundoff protection for cases very
	J=MIN(IT14(IRES+1)-1,J)			!near an edge of the cube
	IL=MOD(I,128)
	IH=I/128
	JL=MOD(J,128)
	JH=J/128
	NINFACE = IX(IL+1)+IY(JL+1)+16384*(IX(IH+1)+IY(JH+1))
	IPIX(IE) = NFACE*IT28(IRES+1) + NINFACE
	ENDDO	! IE=1,NE
	RETURN
	END
