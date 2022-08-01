C	Gorski_map2ps.f - to make a Postscript file from a HEALPIX map, any
C	resolution up to 12582912 pixels (IRES=10)
C	uses Ned Wright's scheme for indexed color to reduce file size
C	link with SHOWPXL and CSORT
C	to make:
C gfortran Gorski_map2ps.f showpxl.f ~/Ned.a
C
C	the input file is unformatted binary with a single record
C	containing an array of size NPIX=12*4^IRES entries.
C	either REAL*4 or INTEGER*2 are allowed.  The size of the file
C	determines IRES and whether it is I*2 or R*4
C
C	output file is a Postscript file with the map
C	NX is the number of pixels in the Postscript file across the width
C	of the map.  NX=800 is usually good enough
C
C	externals in showpxl.f:
C	CALL AWORLD(X,Y,G) returns unit vector on map for position X,Y
C	CALL WORLD(G,X,Y) returns X and Y point on map corresponding to
C		unit vector G
C	SET_SCALE converts map into 0..255
C
	PARAMETER (IRESMAX=10,NPMAX=12*4**IRESMAX,NPM1 = NPMAX-1)
	REAL*4 BIG(0:NPM1)
	INTEGER*2 IBIG(0:NPM1)
	REAL*8 SSCL(33:118)
	INTEGER*4 NSCL(33:118)
	REAL*4 G(3),C(3),E(3),AMAT(3,3),TMAT(3,3),VMIN,VMAX,DTR,AL,B
	PARAMETER (DTR = 57.2957795130824)
	REAL BLINE(4)/-45.,0.001,45.,-0.001/
	CHARACTER*132 OUT,TITLE
	CHARACTER*4 CMD
	CHARACTER*14 FMT
	CHARACTER*64 FN
	INTEGER I,NC,NT,J,IX,IY,IV
C
	INTEGER ISW,ICD
	COMMON /CNTRL/ ISW,ICD
	COMMON /SCALEBAR/ SCALE(6)
	INTEGER IR(0:255),IG(0:255),IB(0:255)
	COMMON /COLORS/ IR,IG,IB
	COMMON /GNOMONIC/ FOCAL
C
C ISW = type of MAP: 1..5, Aitoff, Molleweide,sinusoidal,polar,lat/long
C ICD = display coords: 1..3, galactic, celestial, ecliptic
C	but note that input is supposed to be ecliptic
	ISW = 2		! Molleweide
	ICD = 1		! same as input, means galactic for DMR
C
C	try to read a BIG map
C
	NPX = NPMAX
1	WRITE (*,900) 'ENTER FILE NAME:'
	READ (*,898) FN
	NC = LNBLNK(FN)
	OPEN (UNIT=1,FILE=FN(1:NC),STATUS='OLD',
	1	FORM='UNFORMATTED',ERR=3)
	GO TO 5
3	WRITE (*,FMT='(A)') 'File not found, try again.'
	GO TO 1
5	READ (1,ERR=10,END=10) (BIG(I),I=0,NPX-1)
	CLOSE (UNIT=1)
	GO TO 20
C
C	BIG read failed so try I*2
C
10	REWIND 1
	READ (1,ERR=15,END=15) (IBIG(I),I=0,NPX-1)
	CLOSE (UNIT=1)
	DO I=0,NPX-1
	  BIG(I) = IBIG(I)
	ENDDO
	GO TO 20
C
C	IBIG read failed so try again for lower resolution
C
15	REWIND 1
	NPX = NPX/4
	GO TO 5
20	IRES = ALOG(NPX/12.)/ALOG(4.) + 0.5
	WRITE (*,*) 'NPX, IRES = ',NPX,IRES
	CALL SET_SCALE(NPX,BIG,IBIG)
	WRITE (*,FMT='(1P6E13.5)') SCALE
        DO IV=33,118
          SSCL(IV) = 0
          NSCL(IV) = 0
        ENDDO
        DO NP=0,(NPX-1)
          IV = IBIG(NP)/3.0 + 33.5
          IV = MAX(33,MIN(IV,118))
          NSCL(IV) = NSCL(IV)+1
          SSCL(IV) = SSCL(IV)+BIG(NP)
        ENDDO
        DO IV=33,118
          IF (NSCL(IV).GT.0) SSCL(IV) = SSCL(IV)/NSCL(IV)
        ENDDO
30	WRITE (*,900) 'DISPLAY COORDINATES: 1=GALACTIC, 2=CELESTIAL, '//
	1	'3=ECLIPTIC:'
	READ (*,*) IC
	IF (IC.LT.1 .OR.IC.GT.3) GO TO 30
40	WRITE (*,900) 'INPUT COORDINATES: 1=GALACTIC, 2=CELESTIAL, '//
	1	'3=ECLIPTIC:'
	READ (*,*) JC
	IF (JC.LT.1 .OR.JC.GT.3) GO TO 40
50	WRITE (*,*) 
	1	'Enter 1 for Aitoff, 2 for Mollweide, 3 for Sinusoidal,'
	WRITE (*,900) '      4 for Polar, 5 for lat/long, 6 for gnomonic:'
	READ (*,*) ISW
	IF (ISW.LT.1 .OR. ISW.GT.7) GO TO 50
C
C       MAKE THE TRANSFORMATION MATRIX
C
        DO I=1,3
C       make I'th unit vector -> G
          DO J=1,3
            G(J)=0.
          ENDDO
          G(I) = 1.0
          write (*,*) 'display ',g
          IF (IC.EQ.2) THEN
            DO J=1,3
              C(J)=G(J)
            ENDDO
            CALL CTOG(C,G)
          ELSE IF (IC.EQ.3) THEN
            CALL ETOC(G,C)
            CALL CTOG(C,G)
          ENDIF
          write (*,*) 'galactic',g
          IF (JC.EQ.1) THEN
            DO J=1,3
              E(J) = G(J)
            ENDDO
          ELSE IF (JC.EQ.2) THEN
            CALL GTOC(G,E)
          ELSE IF(JC.EQ.3) THEN
            CALL GTOC(G,C)
            CALL CTOE(C,E)
          ENDIF
          write (*,*) 'ecliptic',e
          DO J=1,3
            AMAT(J,I) = E(J)
          ENDDO
        ENDDO
	IF (ISW.GE.6) THEN
	  WRITE (*,900) 'Enter center long,lat in degrees:'
	  READ (*,*) ALONG,ALAT
C	multiply A*(long mat)*(lat mat)
	  SL = SIN(ALONG/DTR)
	  CL = COS(ALONG/DTR)
C	|  c -s  0 |
C	|  s  c  0 |
C	|  0  0  1 |
	  DO I=1,3
	    TMAT(I,1) = AMAT(I,1)*CL + AMAT(I,2)*SL
	    TMAT(I,2) = AMAT(I,2)*CL - AMAT(I,1)*SL
	    TMAT(I,3) = AMAT(I,3)
	  ENDDO
	  SL = SIN(ALAT/DTR)
	  CL = COS(ALAT/DTR)
C	map (1,0,0) to (cos(lat),0,sin(lat)) and (.995,0,.1) to 
C	(cos(lat+0.1),0,sin(lat+0.1))
C	| c  0 -s |
C	| 0  1  0 |
C	| s  0  c |
	  DO I=1,3
	    AMAT(I,1) = CL*TMAT(I,1) + SL*TMAT(I,3)
	    AMAT(I,2) = TMAT(I,2)
	    AMAT(I,3) = CL*TMAT(I,3) - SL*TMAT(I,1)
	  ENDDO
	  WRITE (*,900) 'Enter expansion factor:'
	  READ (*,*) FOCAL
	ENDIF
        write (*,fmt='(3f12.5)') ((amat(i,j),j=1,3),i=1,3)
	WRITE (*,900) 'ENTER NX:'
	READ (*,*) NX
	NX = 16*((NX+8)/16)
	WRITE (*,*) 'NX = ',NX
	NY = NX/2
	WRITE (*,900) 'ENTER TITLE:'
	READ (*,898) TITLE
	NT = LNBLNK(TITLE)
898	FORMAT(A)
900	FORMAT(A,$)
	CALL GET_COLOR_BAR
        DO I=0,85
          J = 3*I + 0.5
          IR(I) = IR(J)
          IG(I) = IG(J)
          IB(I) = IB(J)
        ENDDO
C
C	write Postscript header
C
	WRITE (*,900) 'ENTER OUTPUT FILE NAME:'
	READ (*,898) FN
	NC = LNBLNK(FN)
	OPEN (UNIT=2,FILE=FN(1:NC),STATUS='UNKNOWN')
	WRITE (2,902) '%!PS-Adobe-1.0'
	WRITE (2,902) '%%Title: '//FN(1:NC)//', '//TITLE(1:NT)
	WRITE (2,902) '%%Creator: gorski_map2ps by Ned Wright'
	WRITE (2,902) '%%DocumentFonts: Times-Roman'
	WRITE (2,902) '%%BoundingBox: 17 57 594 418'
	WRITE (2,902) '%%Pages: 1'
	WRITE (2,902) '%%EndComments'
	WRITE (2,902) '/M {moveto} def /mt {2 div M} def'
	WRITE (2,902) '/L {lineto} def /lnt {2 div L} def'
	WRITE (2,902) '/S {stroke} def /NP {newpath} def'
	WRITE (2,902) '/SRGB {setrgbcolor} def'
	WRITE (2,902) '/labelat {dup stringwidth pop -2 div 0 rmoveto'
	1	//' show} def'
C       output the color table
        WRITE (2,902) 'currentfile 258 string readhexstring'
        WRITE (2,FMT='(24Z2.2)') (IR(I),IG(I),IB(I),I=0,85)
	WRITE (2,902) 'pop /ctbl exch def'
C	define strings for storage.  picstr is one line of file
	WRITE (2,902) '/picstr 128 string def /cpicstr 384 string def'
C	routine to convert string with values 0..255 into a 3x longer
C	string with R,G,B values
	WRITE (2,902) '/ctcvrt {/outpos 0 def {33 sub 3 mul ctbl exch 3'
	WRITE (2,902) 'getinterval cpicstr outpos 3 -1 roll putinterval'
	WRITE (2,902) 'outpos 3 add /outpos exch def} forall cpicstr} def'
	WRITE (2,902) '% detailed scalebar:'
	WRITE (2,FMT='(2H% ,1P7E11.3)') SSCL
	WRITE (2,902) '% number of pixels at each level:'
	WRITE (2,FMT='(2H% ,7I11)') NSCL
C
	WRITE (2,902) 
	WRITE (2,902) '%%Page: 1 1'
	WRITE (2,902) 'gsave'
	WRITE (2,902) '18 100 translate 0.8 dup scale'
	WRITE (2,902) '10 72 mul 10000 div 10 72 mul 10000 div scale'
	WRITE (2,902) 'gsave 10000 5000 scale'
	IF (ISW.LE.2) THEN
	  WRITE (2,902) 'NP 0.5 0.5 0.50 0 360 arc clip'
	ELSE IF (ISW.EQ.3) THEN
C	clipping path for sinusoidal map
	  WRITE (2,902) 'NP 0.5 1 M'
	  DO IY = 1,30
	    X = 0.5-0.5*SIN(6.*IY/DTR)
	    Y = 1-IY/30.
	    WRITE (2,FMT='(2F7.4,2H L)') X, Y
	  ENDDO
	  DO IY = 1,29
	    X = 0.5+0.5*SIN(6.*IY/DTR)
	    Y = 0.00+1.00*IY/30
	    WRITE (2,FMT='(2F7.4,2H L)') X, Y
	  ENDDO
	  WRITE (2,902) 'closepath clip'
	ELSE IF (ISW.EQ.4) THEN
C	clipping path for polar eyeball map
	  WRITE (2,902) 'NP 1 0.5 M'
	  DO IY = 1,30
	    X = 0.75+0.25*COS(6.*IY/DTR)
	    Y = 0.5+0.5*SIN(5.999*IY/DTR)
	    WRITE (2,FMT='(2F7.4,2H L)') X, Y
	  ENDDO
	  DO IY = 1,30
	    X = 0.25+0.25*COS(6.*IY/DTR)
	    Y = 0.5+0.5*SIN(5.999*IY/DTR)
	    WRITE (2,FMT='(2F7.4,2H L)') X, Y
	  ENDDO
	  DO IY = 1,30
	    X = 0.25-0.25*COS(6.*IY/DTR)
	    Y = 0.5-0.5*SIN(5.999*IY/DTR)
	    WRITE (2,FMT='(2F7.4,2H L)') X, Y
	  ENDDO
	  DO IY = 1,29
	    X = 0.75-0.25*COS(6.*IY/DTR)
	    Y = 0.5-0.5*SIN(5.999*IY/DTR)
	    WRITE (2,FMT='(2F7.4,2H L)') X, Y
	  ENDDO
	  WRITE (2,902) 'closepath clip'
	ELSE IF (ISW.EQ.5) THEN
	  WRITE (2,902) 'NP 0 0 M 1 0 L 1 1 L'
	  WRITE (2,902) '0 1 L closepath clip'
	ELSE
	  WRITE (2,902) 'NP 0.125 0 M 0.875 0 L 0.875 1 L'
	  WRITE (2,902) '0.125 1 L closepath clip'
	ENDIF
C	WRITE (2,902) '512 256 8 [512 0 0 256 0 0]'
	WRITE (2,FMT='(I4,I5,4H 8 [,I4,4H 0 0,I5,5H 0 0])') NX,NY,NX,NY
	WRITE (2,902) '{currentfile picstr readline pop ctcvrt}'
	WRITE (2,902) 'false 3 colorimage'
902	FORMAT(A)
C
C	do the ellipse
C
	NC = 0
	DO IY=1,NY
	  Y = (IY-0.5*(1+NY))*(2.0/NY)
	  DO IX = 1,NX
	    X = (IX-0.5*(1+NX))*(4.0/NX)
	    YY = Y
	    CALL AWORLD(X,YY,G)
	    CALL MAV(AMAT,G,E)
	    CALL NESTED_PIXNO(E,IRES,NP)
	    IV = IBIG(NP)/3.0 + 33.5
	    IV = MAX(33,MIN(IV,118))
	    OUT(NC+1:NC+1) = CHAR(IV)
	    NC = NC+1
	    IF (NC.GE.128) THEN
C       fudge the lines starting with %% to avoid dvips bug
	      IF (OUT(1:2).EQ.'%%') OUT(1:1) = '$'
	      IF (OUT(1:2).EQ.'%!') OUT(1:1) = '&'
	      WRITE (2,902) OUT(1:NC)
	      NC = 0
	    ENDIF
	  ENDDO
	ENDDO
C
C	draw in the grid
C
	WRITE (2,902) 'grestore'	! stop clipping, back to 10000 scale
	WRITE (2,902) '10 setlinewidth 0 setgray'
	WRITE (*,900) 'Grid?:'
	READ (*,898) FN
	NC = LNBLNK(FN)
	IF (NC.LE.0 .OR. (FN(1:1).NE.'Y' .AND. FN(1:1).NE.'y')) GO TO 160
	DO I=-4,4
	  AL = 44.999*I/DTR
	  WRITE (2,902) 'NP'
	  CMD = '  mt'
	  DO J=-90,0,3
	    B = (0.99999*J-0.0002)/DTR
	    G(3) = SIN(B)
	    G(2) = COS(B)*SIN(AL)
	    G(1) = COS(B)*COS(AL)
	    CALL WORLD(G,X,Y)
	    IX = 5000.5+2500.0*X
	    IF (IX.GT.9999) IX=9999
	    IY = 5000.5+5000.0*Y
	    IF (IY.GT.9999) IY=9999
	    WRITE (2,904) IX,IY,CMD
904	FORMAT(2I5,A)
	    CMD = ' lnt'
	  ENDDO
	  WRITE (2,902) 'S NP'
	  CMD = '  mt'
	  DO J=0,90,3
	    B = (0.99999*J+0.0002)/DTR
	    G(3) = SIN(B)
	    G(2) = COS(B)*SIN(AL)
	    G(1) = COS(B)*COS(AL)
	    CALL WORLD(G,X,Y)
	    IX = 5000.5+2500.0*X
	    IF (IX.GT.9999) IX=9999
	    IY = 5000.5+5000.0*Y
	    IF (IY.GT.9999) IY=9999
	    WRITE (2,904) IX,IY,CMD
	    CMD = ' lnt'
	  ENDDO
	  WRITE (2,902) 'S'
	ENDDO
	NB = 3
	IF (ISW.EQ.4) NB = 4
	DO J=1,NB
	  B = BLINE(J)/DTR
	  WRITE (2,902) 'NP'
	  CMD = '  mt'
	  DO I=-90,90,2
	    AL = 1.9999*I/DTR
	    G(3) =SIN(B)
	    G(2) = COS(B)*SIN(AL)
	    G(1) = COS(B)*COS(AL)
	    CALL WORLD(G,X,Y)
	    IX = 5000.5+2500.0*X
	    IF (IX.GT.9999) IX=9999
	    IY = 5000.5+5000.0*Y
	    IF (IY.GT.9999) IY=9999
	    WRITE (2,904) IX,IY,CMD
	    CMD = ' lnt'
	  ENDDO
	  WRITE (2,902) 'S'
	ENDDO
C
160	WRITE (*,900) 'Scale Bar?:'
	READ (*,898) FN
	NC = LNBLNK(FN)
	IF (NC.LE.0 .OR. (FN(1:1).NE.'Y'.AND. FN(1:1).NE.'y')) GO TO 180
	WRITE (2,902) '/Times-Roman findfont 250 scalefont setfont'
C
C	label scale bar
C
        VMAX = MAX(SCALE(1),SCALE(6))
        VMIN = MIN(SCALE(1),SCALE(6))
        VMAX = AMAX1(VMAX,-10.*VMIN)
        FMT = 'F6.0'
        IF (VMAX.LT.9999.) FMT = 'F6.1'
        IF (VMAX.LT.999.9) FMT = 'F6.2'
        IF (VMAX.LT.99.99) FMT = 'F6.3'
        IF (VMAX.LT.9.999) FMT = 'F6.4'
	FMT = '(1H(,'//FMT(1:4)//',1H))'
	DO I=1,6
	  IX = 1000+1600*(I-1)
	  WRITE (2,905) IX
905	FORMAT(I4,' -700 moveto')
	  WRITE (2,FMT) SCALE(I)
	  WRITE (2,902) 'labelat'
	ENDDO
C
C	draw scale bar and print
C
	WRITE (2,902) 'gsave 1000 -400 translate 8000 200 scale'
	WRITE (2,902) '86 1 8 [86 0 0 1 0 0]'
	WRITE (2,902) '{ ctbl }'
	WRITE (2,902) 'false 3 colorimage'
	WRITE (2,902) 'grestore'
C
C	do TITLE
C
180	WRITE (2,902) '/Times-Roman findfont 500 scalefont setfont'
	WRITE (2,902) '5000 5200 moveto'
	WRITE (2,902) '('//TITLE(1:NT)//')'
	WRITE (2,902) 'labelat'
	WRITE (2,902) 'grestore'
	WRITE (2,902) 'showpage'
	WRITE (2,902) '%%Trailer'
	STOP
	END

	SUBROUTINE SET_SCALE(NPX,R,IMAP)
C
C	CONVERTS REAL MAP IN ARRAY R INTO DISPLAY RANGE IN
C	IMAP.  VALUES IN IMAP RANGE FROM 1 TO IMAX-1
C	Values less than THRESHOLD are assumed to be missing data
C
	PARAMETER (THRESHOLD=-1.0E10)
	PARAMETER (IRESMAX=10,NPIX=12*4**IRESMAX)
	PARAMETER (IMAX=255)
	INTEGER*4 NPX,I,J
	REAL*4 R(NPIX)
	INTEGER*2  IMAP(NPIX)
	COMPLEX CM(NPIX)
	REAL*4 RM(2,NPIX)
	INTEGER*4 IM(2,NPIX)
	EQUIVALENCE (CM,RM,IM)
	LOGICAL WRAP
	COMMON /SCALEBAR/ SCALE(6)
	DO I=1,NPX
	  RM(1,I) = R(I)
	  IM(2,I) = I
	ENDDO
	CALL CSORT(NPX,CM)
	VMAX = REAL(CM(NPX))
	I = 1
1	VMIN = REAL(CM(I))
	IF (VMIN.LE.THRESHOLD) THEN
	  I=I+1
	  IF (I.LE.NPX) GO TO 1
	ENDIF
	ILOW = I
	PENULTIMATE = REAL(CM(NPX-1))
	WRITE (*,*) 'SECOND HIGHEST VALUE = ',PENULTIMATE
	WRITE (*,*) 'RANGE OF VALUES: ',VMIN,' TO ',VMAX
	WRITE (*,*) 'NEGATIVE SUPPRESSES WRAPPING'
	WRITE (*,901) 'ENTER 1 FOR LIN, 2 FOR LOG, 3 FOR HIST'//
	1	', 4 FOR quasiLOG, 5 for pt src:'
901	FORMAT(1X,A,$)
	READ (*,902) IANS
902	FORMAT(8I10)
	WRAP = (IANS.GT.0)
	IANS = ABS(IANS)
	VMN = 0.
	VMX = 0.
	IF (IANS.EQ.1) THEN
C		LINEAR
	  WRITE (*,901) 'ENTER NEW VALUES FOR VMIN AND VMAX:'
	  READ (*,903) VMN,VMX
903	FORMAT(8E10.3)
	  IF (VMN.NE.0. .OR. VMX.NE.0.) THEN
	    VMIN = VMN
	    VMAX = VMX
	  ENDIF
	  DO I=1,NPX
	    IF (R(I).LE.THRESHOLD) THEN
	      IMAP(I) = 0
	    ELSE
	      IMAP(I) = 1.0 + (IMAX-0.001)*(R(I)-VMIN)/(VMAX-VMIN)
	      IF (.NOT.WRAP) IMAP(I) = MIN(IMAX,MAX(0,IMAP(I)))
	      IF (WRAP) IMAP(I) = MOD(IMAP(I),IMAX+1)
	    ENDIF
	  ENDDO
	  DO I=1,6
	    SCALE(I) = VMIN + (I-1)*(VMAX-VMIN)/5.0
	  ENDDO
	ELSE IF (IANS.EQ.2) THEN
C		LOGARITHMIC
10	  WRITE (*,901) 'ENTER NEW VALUES FOR VMIN AND VMAX:'
	  READ (*,903) VMN,VMX
	  IF (VMN.NE.0. .OR. VMX.NE.0.) THEN
	    VMIN = VMN
	    VMAX = VMX
	  ENDIF
	  IF (VMIN.LE.0. .OR. VMAX.LE.0.) GO TO 10
	  DO I=1,NPX
	    IMAP(I) = 0
	    IF (R(I).GT.THRESHOLD) IMAP(I) = 1
	    IF (R(I).GT.VMIN)
	1     IMAP(I) = 1.+(IMAX-.001)*ALOG(R(I)/VMIN)/ALOG(VMAX/VMIN)
	    IF (.NOT.WRAP) IMAP(I) = MIN(IMAX,MAX(0,IMAP(I)))
	    IF (WRAP) IMAP(I) = MOD(IMAP(I),IMAX+1)
	  ENDDO
	  DO I=1,6
	    SCALE(I) = VMIN*EXP((I-1)*ALOG(VMAX/VMIN)/5.0)
	  ENDDO
        ELSE IF (IANS.EQ.4) THEN
C		QUASI-LOGARITHMIC
	  WRITE (*,901) 'ENTER NEW VALUES FOR VMIN AND VMAX:'
	  READ (*,903) VMN,VMX
	  IF (VMN.NE.0. .OR. VMX.NE.0.) THEN
	    VMIN = VMN
	    VMAX = VMX
	  ENDIF
	  DV = MAX(VMIN,-2.0*VMIN)
	  DO I=1,NPX
	    IMAP(I) = 0
	    IF (R(I).GT.THRESHOLD) IMAP(I) = 1
	    IF (R(I).GT.VMIN)
	1     IMAP(I) = 1.5 + (IMAX-1)*
	2     SQRT(ALOG((R(I)+DV)/(VMIN+DV))/ALOG((VMAX+DV)/(VMIN+DV)))
	    IF (.NOT.WRAP) IMAP(I) = MIN(IMAX,MAX(0,IMAP(I)))
	    IF (WRAP) IMAP(I) = MOD(IMAP(I),IMAX+1)
	  ENDDO
	  DO I=1,6
	    SCALE(I) = (VMIN+DV)*EXP(((I-1.)/5.)**2*
	1	ALOG((VMAX+DV)/(VMIN+DV)))-DV
          ENDDO
        ELSE IF (IANS.EQ.4) THEN
C	MODIFIED HISTOGRAM EQUALIZED TO SHOW STARS WELL
	  DO I=1,NPX
	    J = IM(2,I)
	    IMAP(J) = 0
	    IF (I.GE.ILOW) THEN
	      X = FLOAT(I-ILOW)/FLOAT(NPX-ILOW)
C	x=0 -> 0; 0.5 -> 0.167, 0.9 -> 0.45
	      XX = 0.1*X/(1.1-X)
	      IF (.NOT.WRAP) XX = 0.05*X/(1.05-X)
	      IMAP(J) = 1.0+(IMAX-.001)*XX
	    ENDIF
	  ENDDO
	  DO I=1,6
	    XX = (I-1.0)/5.0
C	1.1*XX - X*XX = 0.1*X, X=1.1*XX/(XX+0.1)
	    X = 1.1*XX/(XX+0.1)
	    IF (.NOT.WRAP) X = 1.05*XX/(XX+0.05)
	    J = ILOW+0.5+X*(NPX-ILOW)
	    SCALE(I) = REAL(CM(J))
	  ENDDO
	ELSE
C		HISTOGRAM EQUALIZED
	  DO I=1,NPX
	    J = IM(2,I)
	    IMAP(J) = 0
	    IF (I.GE.ILOW) THEN
	      X = FLOAT(I-ILOW)/FLOAT(NPX-ILOW)
	      IF (.NOT.WRAP) THEN
	        X = 2*(X-0.5)
	        X = X*(1+X*X)/2.
	        X = (1+X)/2.
	      ENDIF
	      IMAP(J) = 1.0+(IMAX-.001)*X
	    ENDIF
	  ENDDO
	  DO I=1,6
	    X = (I-1.0)/5.0
	    IF (.NOT.WRAP) THEN
	      Y = 2*(X-0.5)
	      X = Y
	      DO J=1,10
	        E = X*(1+X*X)/2. - Y
	        DEDX = 0.5*(1+3*X*X)
	        X = X-E/DEDX
	      ENDDO
	      X = (1+X)/2
	    ENDIF
	    J = ILOW+0.5+X*(NPX-ILOW)
	    SCALE(I) = REAL(CM(J))
	  ENDDO
	ENDIF
	RETURN
	END

	SUBROUTINE GET_COLOR_BAR
	CHARACTER*80 FN
	INTEGER IFAKE,IR(256),IG(256),IB(256),JR(256),JG(256),
	1	JB(256)
	COMMON /COLORS/ IR,IG,IB
	BYTE BFAKE
	EQUIVALENCE (BFAKE,IFAKE)
	PARAMETER (NR = 4, NG = 5, NB = 5)
	PARAMETER (WHITEN = 0.25)
	REAL*4 XR(NR),YR(NR),XG(NG),YG(NG),XB(NB),YB(NB)
	LOGICAL USEFILE /.FALSE./
	DATA XR/0., 0.20, 0.68, 1.0/
	DATA YR/0., 0.0, 1.0, 1.0/
	DATA XG/0., 0.005, 0.47, 0.995, 1.0/
	DATA YG/0.7, 0.7,   0.0,  0.0,   0.0/
	DATA XB/0., 0.005, 0.68, 0.995, 1.0/
	DATA YB/1.0, 1.0,   1.0,  0.0,   0.0/
10	WRITE (*,FMT='(1X,A,$)') 'Enter Color File name:'
	READ (*,FMT='(A)') FN
	NC = LNBLNK(FN)
	USEFILE = (NC.GT.0)
	IF (USEFILE) THEN
	  OPEN (UNIT=2,FILE=FN(1:NC),STATUS='OLD',ERR=20)
	  GO TO 30
20	  WRITE (*,FMT='(A)') 'File not found, try again.'
	  GO TO 10
30	  SR = 0.
	  SG = 0.
	  SB = 0.
	  DO I=1,256
	    READ (2,FMT='(3Z3)') JR(I),JG(I),JB(I)
	    SR = SR+JR(I)
	    SG = SG+JG(I)
	    SB = SB+JB(I)
	  ENDDO
	  CLOSE (UNIT=2)
	  WRITE (*,*) 'Sum of RGB = ',SR,SG,SB
	  DO I=1,256
	    J = I
	    IR(I) = JR(J)
	    IG(I) = JG(J)
	    IB(I) = JB(J)
	  ENDDO
	ELSE
	  DO I=1,256
	    X = (I-1.)/255.
	    IR(I) = 255.0*((1.-WHITEN)*TERP(X,NR,XR,YR)+WHITEN)
	    IG(I) = 255.0*((1.-WHITEN)*TERP(X,NG,XG,YG)+WHITEN)
	    IB(I) = 255.0*((1.-WHITEN)*TERP(X,NB,XB,YB)+WHITEN)
	  ENDDO
	ENDIF
	RETURN
	END

	FUNCTION TERP(X,N,XT,YT)
	REAL*4 X, XT(N), YT(N)
	I = 1
	DO WHILE (I.LT.N-1 .AND. X.GT.XT(I+1))
	  I = I+1
	ENDDO
	F = (X-XT(I))/(XT(I+1)-XT(I))
	TERP = YT(I) + (YT(I+1)-YT(I))*F
	RETURN
	END
