C	robust-average-singles.f
C	gfortran robust-average-singles.f ~/missions/WISE/TLEs-orbit.f ~/Ned.a
C	to process the IRSA table on a text view of a NEOWISER .tbl
C	from a moving object search
C	with RA,DEC,W1,W2,ccflags,MJD,w1flux,sw1,w2flux,sw2 selected
C	first 5 columns from IRSA are cntr_u,dist_x,pang_x,ra_u,dec_u  
C	file on standard input
C
	CHARACTER*400 INL
	CHARACTER*80 FN
	CHARACTER*7 PACKED
	CHARACTER*1 CODE
	INTEGER L(100,2)
	PARAMETER (NLMX=199)
	INTEGER NC(NLMX),NT(NLMX)
	REAL*8 VALUE
	REAL*8 MJD(NLMX),RA(NLMX),DEC(NLMX),FLUX(4,NLMX),
	1	FSIG(4,NLMX)
	REAL*8 MAG(4),SMAG(4)
	REAL*8 SZP(4),NZP(4),SZZ(4),ZP(4)
	REAL*8 SRA,SDEC,SMJD,RAC,DECC,MJDC,SF,SFF
	REAL*8 FLX(NLMX),SIG(NLMX),F,SIGF,RCHI2,MBAR,SIGM
	LOGICAL IFNUM
C
C	see http://www.minorplanetcenter.net/iau/info/ObsNote.html for codes
C	'c' = crowded field; 
	CODE = ' '
C	packed designation
	PACKED = 'N00bzsy'
	CALL INIT_ORBIT
	WRITE (*,FMT='(A,$)') 'input file:'
	READ (*,FMT='(A)') FN
	NCF = LNBLNK(FN)
	OPEN (UNIT=1,FILE=FN(1:NCF))
	WRITE (*,FMT='(A,$)') 'enter NSKIP,DMAX,PACKED:'
	READ (*,*) NSKIP,DMAX,PACKED
	DO ISKIP=1,NSKIP
	  READ (1,FMT='(A)') INL
	ENDDO
C	no nulls in the table
	DO IB = 1,2
	  SZP(IB) = 0
	  NZP(IB) = 0
	  SZZ(IB) = 0
	ENDDO
	IR = 0
10	  READ (1,FMT='(A)',END=20) INL
	  NC(IR+1) = LNBLNK(INL)
	  CALL TOKENIZE(INL,NTOKEN,L)
	  IF (IFNUM(INL(L(2,1):L(2,2)),VALUE)) DIST = VALUE
C	  if the distance from the predicted position is larger than DMAX, skip the line
	  IF (DIST.GT.DMAX) GO TO 10
	  IF (IFNUM(INL(L(11,1):L(11,2)),VALUE)) MJD(IR+1) = VALUE
	  IR = IR+1
	  IF (IFNUM(INL(L(6,1):L(6,2)),VALUE)) RA(IR) = VALUE
	  IF (IFNUM(INL(L(7,1):L(7,2)),VALUE)) DEC(IR) = VALUE
	  DO IB=1,2
	    IF (IFNUM(INL(L(IB+7,1):L(IB+7,2)),VALUE)) MAG(IB) = VALUE
	    IF (IFNUM(INL(L(2*IB+10,1):L(2*IB+10,2)),VALUE)) 
	1	FLUX(IB,IR) = VALUE
	    IF (IFNUM(INL(L(2*IB+11,1):L(2*IB+11,2)),VALUE)) 
	1	FSIG(IB,IR) = VALUE
	    IF (FLUX(IB,IR).GT.2.0D0*FSIG(IB,IR)) THEN
	      ZP(IB) = MAG(IB)+2.5D0*LOG10(FLUX(IB,IR))
	    ELSE
	      ZP(IB) = MAG(IB)+2.5D0*LOG10(MAX(0.D0,FLUX(IB,IR))+
	1	2*FSIG(IB,IR))
	    ENDIF
	    NZP(IB) = NZP(IB)+1
	    SZP(IB) = SZP(IB)+ZP(IB)
	    SZZ(IB) = SZZ(IB)+ZP(IB)**2
	  ENDDO
	  NT(IR) = NTOKEN
	  CALL MPC_FORMAT(PACKED,CODE,MJD(IR),RA(IR),DEC(IR))
	  GO TO 10
20	NLINE = IR
	WRITE (*,FMT='(18I4)') (NC(I),I=1,NLINE)
	WRITE (*,FMT='(18I4)') (NT(I),I=1,NLINE)
	DO IB = 1,2
	  WRITE (*,*) NZP(IB),SZP(IB),SZZ(IB)
	  ZP(IB) = SZP(IB)/NZP(IB)
	  SZP(IB) = SQRT((SZZ(IB)/NZP(IB)-ZP(IB)**2))
	  WRITE (*,*) ZP(IB),SZP(IB)
	ENDDO
	SRA = 0
	SDEC = 0
	SMJD = 0
	DO I=1,NLINE
	  SRA = SRA+RA(I)
	  SDEC = SDEC+DEC(I)
	  SMJD = SMJD+MJD(I)
	ENDDO
	DO IB=1,2
	  SF = 0
	  SW = 0
	  DO I=1,NLINE
	    W = 1/FSIG(IB,I)**2
	    SW = SW+W
	    SF = SF + W*FLUX(IB,I)
	  ENDDO
	  SF = SF/SW
	  SFF = 0
	  DO I=1,NLINE
	    SFF = SFF+((FLUX(IB,I)-SF)/FSIG(IB,I))**2
	  ENDDO
	  WRITE (*,*) IB,' chi^2= ',SFF
	  SFF = MAX(NLINE-1.D0,SFF)
	  MAG(IB) = ZP(IB)-2.5D0*LOG10(SF)
	  SMAG(IB) = 1.0857362*SQRT(SFF/(NLINE-1))/(SQRT(SW)*SF)
C	  write (*,*) IB,MAG(IB),SFF
	ENDDO
	RAC = SRA/NLINE
	DECC = SDEC/NLINE
	MJDC = SMJD/NLINE
	DO IB = 1,2
C	  get weighted mean of fluxes and sigmas
	  SW=0
	  SWF=0
	  DO IR=1,NLINE
	    FLX(IR) = FLUX(IB,IR)
	    SIG(IR) = FSIG(IB,IR)
	    SW=SW+1/FSIG(IB,IR)**2
	    SWF=SWF+FLUX(IB,IR)/FSIG(IB,IR)**2
	  ENDDO
	  F = SWF/SW
	  SIGF = 1/SQRT(SW)
	  CHI2 = 0
	  DO IR=1,NLINE
	    CHI2 = CHI2+((FLUX(IB,IR)-F)/FSIG(IB,IR))**2
	  ENDDO
	  MBAR = ZP(IB)-2.5*LOG10(ABS(F))
	  SIGM = 2.5*LOG10(1+SIGF/ABS(F))
	  IF (F.LT.0.) SIGM = -SIGM
C	  get robust mean of fluxes and sigmas
	  CALL ROBUST(NLINE,FLX,SIG,F,SIGF,RCHI2)
	  RMAG = ZP(IB)-2.5*LOG10(ABS(F))
	  MAG(IB) = RMAG
	  RSIG = 2.5*LOG10(1+SIGF/ABS(F))
	  IF (F.LT.0.) RSIG = -RSIG
	  SMAG(IB) = RSIG*SQRT(MAX(1.,RCHI2/(NLINE-1)))
	  WRITE (*,FMT='(I1,2(F8.3,F6.3,F8.2),1P2E11.3)') 
	1	IB,MBAR,SIGM,CHI2,RMAG,RSIG,RCHI2,F,SIGF
	ENDDO
	WRITE (*,FMT='(F9.3,2F11.6,6F7.3)') 
	1	MJDC,RAC,DECC,(MAG(IB),SMAG(IB),IB=1,2)
	STOP
	END

	SUBROUTINE MPC_FORMAT(PACKED,CODE,MJD,TRA,TDEC)
	REAL*8 XYZ(3),AXYZ(3)
	REAL*4 RXYZ(3)
	REAL*8 MJD,RA,DEC,TRA,TDEC
	REAL*8 JD,JD1JAN2010,T68,T1JAN2010,TUT
	REAL*8 JULIAN_DAY
	REAL*8 DAY,FDAY
	INTEGER IMJD,IH,MON,MINUTE,IDAY
	CHARACTER*36 CXYZ
	CHARACTER*3 CMON
	CHARACTER*9 CDAY
	CHARACTER*12 CRA,CDEC
	CHARACTER*1 CDS,CODE
	CHARACTER*7 PACKED
C
	CALL MJD_2_T68(MJD,T68)
C
	CALL T68_TO_UTC(T68,IYEAR,MON,IDAY,IH,MINUTE,SECOND)
C	FDAY = (SECOND+60*(MINUTE+60*IH))/86400
C	write (*,FMT='(F11.5,F8.5,F15.3,I5,4I3,F8.4)') 
C	1	MJD,FDAY,T68,IYEAR,MON,IDAY,IH,MINUTE,SECOND
	IMJD = MJD
	FDAY = MJD-IMJD
	IH = 24*FDAY
	MINUTE = 1440*(FDAY-IH/24.D0)
	SECOND = 86400*(FDAY-IH/24.D0-MINUTE/1440.D0)
	DAY = FDAY+IDAY
	CALL ORBIT(T68,RXYZ)
	CALL PRECESS2J2000(T68,RXYZ)
	DO I=1,3
	  XYZ(I) = RXYZ(I)
	ENDDO
	RA = TRA
	DEC = TDEC
	IY=IYEAR
	IMON=MON
	IH = RA/15
	RA = 60*(RA/15-IH)
	IRM = RA
	RA = 60*(RA-IRM)
	CDS = '+'
	IF (DEC.LT.0.) CDS = '-'
	DEC = ABS(DEC)
	ID = DEC
	DEC = 60*(DEC-ID)
	IDM = DEC
	DEC = 60*(DEC-IDM)
	DO J=1,3
	  AXYZ(J) = ABS(XYZ(J))
	ENDDO
	WRITE (CXYZ,FMT='(3F12.4)') AXYZ
	WRITE (CMON,FMT='(I3)') IMON
	IF (CMON(2:2).EQ.' ') CMON(2:2) = '0'
	WRITE (CDAY,FMT='(F9.5)') DAY
	IF (CDAY(2:2).EQ.' ') CDAY(2:2) = '0'
	WRITE (CRA,FMT='(2I3,F6.2)') IH,IRM,RA
	IF (CRA(2:2).EQ.' ') CRA(2:2) = '0'
	IF (CRA(5:5).EQ.' ') CRA(5:5) = '0'
	IF (CRA(8:8).EQ.' ') CRA(8:8) = '0'
	WRITE (CDEC,FMT='(I4,I3,F5.1)') ID,IDM,DEC
	IF (CDEC(3:3).EQ.' ') CDEC(3:3) = '0'
	IF (CDEC(6:6).EQ.' ') CDEC(6:6) = '0'
	CDEC(2:2) = CDS
	DO J=1,3
	  IC = 3+12*(J-1)
	  IF(XYZ(J).LT.0) CXYZ(IC:IC) = '-'
	  IF(XYZ(J).GE.0) CXYZ(IC:IC) = '+'
	ENDDO
	WRITE (*,900) PACKED,CODE,IY,CMON,CDAY,CRA,CDEC
	WRITE (*,901) PACKED,CODE,IY,CMON,CDAY,CXYZ
C
	RETURN
900	FORMAT(5X,A,1X,A,1HS,I4,A,A,A,A,22X,3HC51)
901	FORMAT(5X,A,1X,A,1Hs,I4,A,A,2H 1,A,8X,3HC51)
	END
