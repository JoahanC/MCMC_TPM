	SUBROUTINE PROTONS(T,NE,E,DNDE)
C
C	GIVES PROTON FLUX AT COBE
C	BY NED WRIGHT 30-MAR-83
C
C	INPUT:	T IS REAL*8 SECONDS SINCE 5/24/68
C		NE IS NUMBER OF ENERGIES
C		E IS A REAL*4 ARRAY OF ENERGIES IN MEV
C	OUTPUT:	DNDE IS A REAL*4 ARRAY OF DIFFERENTIAL SPECTRUM
C
	REAL*8 T
	DIMENSION E(NE),DNDE(NE),R(3)
	DIMENSION PEAKS(3,3)
C
C	PEAKS CONTAINS POSITIONS OF THREE GAUSSIANS THAT DEFINE THE
C		SOUTH ATLANTIC ANOMALY
	DATA PEAKS/ .3830,-.8214,-.4226,
	1	    .9830,-.1294,-.1305,
	2	    .7789,-.5421,-.3154/
C
	PARAMETER (PI=3.1415927,TWOPI=6.2831853,DTR=57.29578)
C
	CALL TERRESTIAL(T,R)
	CALL NORM(R)
C
C	COMPUTE RATE WITH E GT 50 MEV
	E50=0.
	DO 20 I=1,3
20	E50=1000.*EXP(-8.5934*DVS(R,PEAKS(1,I)))+E50
C
C	ASSUME 1/E**2 SPECTRUM SO E50*50./E**2
	DO 40 I=1,NE
40	DNDE(I)=50.*E50/E(I)**2
	RETURN
	END
