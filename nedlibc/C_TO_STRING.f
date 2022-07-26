	SUBROUTINE C_TO_STRING(DCEL,NDIG,STRING)
C
C	CONVERTS RA+DEC UNIT VECTOR DCEL INTO A FORMATTED STRING FOR OUTPUT
C
C	DCEL IS REAL*4 (3) INPUT ONLY
C	NDIG IS AN INTEGER INPUT THAT CONTROLS PRECISION
C	FOR NDIG = 0, THE OUTPUT IS HHhMMmSS -DDdMM.M = 17 characters
C	           1,               HHhMMmSS.S -DDdMM'SS = 20 chars
C		   2,		    HHhMMmSS.SS -DDdMM'SS.S = 23 chars
C
C	ROUNDING CAN GIVE 60.0 for the seconds PART
C
	CHARACTER*(*) STRING
	CHARACTER*1 IDS
	PARAMETER (PI=3.1415927,TWOPI=6.2831853,DTR=57.295778)
	REAL*4 DCEL(3)
	RA=ATAN2(DCEL(2),DCEL(1))
	IF(RA.LT.0.) RA=RA+TWOPI
	DEC=ATAN2(DCEL(3),SQRT(DCEL(1)**2+DCEL(2)**2))
	IF(DEC.LT.0.) THEN
	  IDS='-'
	  DEC=-DEC
	ELSE
	  IDS=' '
	ENDIF
	CALL R TO DMS(RA/15.,IH,IRM,SR)
	CALL R TO DMS(DEC,ID,IDM,SD)
	IF (NDIG.EQ.2) THEN
	  WRITE (STRING,900) IH,IRM,SR,IDS,ID,IDM,SD
	  IS=13
900	FORMAT(I2,1Hh,I2,1Hm,F5.2,1X,A,I2,1Hd,I2,1H',F4.1)
	ELSE IF(NDIG.EQ.1) THEN
	  ISD=SD+.5
	  IF(ISD.EQ.60) THEN
	    ISD=0
	    IDM=IDM+1
	    IF(IDM.EQ.60) THEN
	      IDM=0
	      ID=ID+1
	    ENDIF
	  ENDIF
	  WRITE (STRING,901) IH,IRM,SR,IDS,ID,IDM,ISD
	  IS=12
901	FORMAT(I2,1Hh,I2,1Hm,F4.1,1X,A,I2,1Hd,I2,1H',I2)
	ELSE
	  RM=IDM+SD/60.
	  ISR=SR+.5
	  IF(ISR.EQ.60) THEN
	    ISR=0
	    IRM=IRM+1
	    IF(IRM.EQ.60) THEN
	      IRM=0
	      IH=IH+1
	      IF(IH.EQ.24) IH=0
	    ENDIF
	  ENDIF
	  WRITE (STRING,902) IH,IRM,ISR,IDS,ID,RM
	  IS=10
902	FORMAT(I2,1Hh,I2,1Hm,I2,1X,A,I2,1Hd,F4.1)
	ENDIF
	IF(STRING(IS:IS+1).EQ.'- ') STRING(IS:IS+1) = ' -'
	RETURN
	END
