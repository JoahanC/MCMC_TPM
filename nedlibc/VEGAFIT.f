	FUNCTION VEGAFIT(WVL)
C       better fit: 1.0158E-16*(1-0.0083*LOG(W/LC)**2)*B_nu(14454)
C	returns Vega flux in Jy for WVL in microns
C
C	table of correction factors for 437 nm to 2015 nm based on
C	Vega.dat file from 
C	http://www.astro.yale.edu/dokkum/newfirm/instrument/filters.html
	INTEGER IFC(32) /2345,1658,1871,2830,2210,2353,2435,
	1       2515,2629,3116,2608,2604,2595,2540,2414,1453,
	2       0733,0708,0433,0743,0464,0459,0695,0463,0466,
	3       0447,0483,0429,0152,0175,0296,0023/
C
	REAL*4 LC /8.891/
	DATA OMEGA,CBEST,TVEGA/1.0158E-16,-0.0083,14454./
	VEGAFIT = OMEGA*(1+CBEST*LOG(WVL/LC)**2)*
	1	PLANCK(10000/WVL,TVEGA)
	2	/2.99792458E-13
	IF (WVL.LT.2.0153) THEN
	  X = MAX(0.,46*LOG10(WVL/0.4270))
	  I = X
	  X = X-I
	  RATIO = 1+0.0001*((1-X)*IFC(I+1)+X*IFC(I+2))
	  VEGAFIT = VEGAFIT/RATIO
	ENDIF
	RETURN
	END
