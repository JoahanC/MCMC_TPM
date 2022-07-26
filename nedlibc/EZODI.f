	SUBROUTINE EZODI(E,R,NNU,NU,INU)
	INTEGER NNU
	REAL*4 E(3),R(3),NU(NNU),INU(NNU)
C
C	zodiacal light model for arbitrary direction and point of view
C	by Ned Wright, UCLA
C
C	E is a unit vector direction of view
C	R is the viewing location, in AU.  Both E and R are in the 
C	ecliptic J2000 frame.  Using 1950 ecliptic is close enough...
C
C	NNU is the number of wavenumbers
C	NU is the list of wavenumbers for which the ZODI light should be found
C	INU is the OUTPUT list of intensites in erg/cm**2/sec/sr/cm**-1
C
C	EXTERNALS: 	DOT	dot product of two 3 vectors
C			NORM	normalized a 3 vector
C			MAV	3x3 matrix times 3 vector -> 3 vector
C			PLANCK	blackbody of NU (/cm) and T (K) in
C				erg/cm^2/sec/sr/(/cm)	(Eplee's)
C
	REAL*4 TSUN,FSCAT,FDUST(2),TDUST(2)
C
C	ZODI LIGHT PARAMETERS: DENSITY LAW
	DATA ALPHA,BETA/1.3,2.1/
C	TEMPERATURES AND EMISSIVITIES
	DATA TSUN,FSCAT/5500.,12.E-14/
C	DUST
	DATA FDUST,TDUST/1.5E-10,9.E-10,265.2,148.1/
C
	REAL*4 A(3,3),X(3),Y(3),B(3)
C
C	MATRIX TO TRANSFORM TO ZODIACAL CLOUD FRAME
	DATA A/.9975,0.,+.05,0.,1.,0.,-.05,0.,.9975/
	CALL MAV(A,E,Y)
	CALL MAV(A,R,X)
	CALL NORM(Y)
C	GET IMPACT PARAMETER
	XY=DOT(X,Y)
	DO 100 I=1,3
100	B(I)=X(I)-XY*Y(I)
	BB=SQRT(DOT(B,B))
	T1=XY/BB
	S=SQRT(1.+T1**2)
	W1=1./(S*(T1+S))
C
C	INTEGRATION VARIABLE IS w = 1 - sin(theta)
C		where X = B + b tan(theta) Y
C		Range of w is from W1 to 0.
C
	NW=4
	DW=W1/NW
	W=DW/2.
C	ZERO THE SUMS
	DO 110 I=1,NNU
110	  INU(I)=0.
C	DO THE INTEGRATION OVER LINE OF SIGHT
	DO 150 IW=1,NW
	  C=SQRT(W*(2.-W))		!cos(theta)
	  S=1.-W
	  T=BB*S/C
C	FIND THE VECTOR POSITION
	  DO 120 I=1,3
120	    X(I)=B(I)+T*Y(I)
	  RR=SQRT(DOT(X,X))
	  Z=X(3)
C	SOLAR FLUX VARIES LIKE R**-2
	  FSUN=1./RR**2
C	OVER 1/LAMBDA DUST T GOES LIKE FSUN**1/5
	  TF=FSUN**0.2
C	DENSITY LAW FROM HELIOS MODEL
	  DENSE=EXP(-BETA*ABS(Z)/RR)/RR**ALPHA
	  DS=BB*DW/C**3
	  PHASE=DOT(Y,X)/RR
C	TRY RAYLEIGH PHASE LAW ON SCATTERED LIGHT, BUT ALBEDO IS GRAY
	  SCAT=FSCAT*DENSE*(1.+PHASE**2)*FSUN*DS
C	ADD SCATTERED LIGHT TO INU
	  DO 130 I=1,NNU
130	    INU(I)=INU(I)+SCAT*PLANCK(NU(I),TSUN)
C	TWO TYPES OF DUST
	  DO 140 ID=1,2
	  THERM=FDUST(ID)*DENSE*DS
	  TD=TDUST(ID)*TF
C	ADD DUST TO INU
	    DO 140 I=1,NNU
140	      INU(I)=INU(I)+THERM*PLANCK(NU(I),TD)*NU(I)
150	  W=W+DW
	RETURN
	END
