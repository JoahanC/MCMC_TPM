	SUBROUTINE RPOLYLINE(NIN,X,Y,XRANGE,YRANGE,XE,NOUT,XOUT,YOUT)
C
C	Scans array of (X,Y) 1..NIN to find set of points
C	(XOUT,YOUT) 1..NOUT that such that straight lines thru
C	the output array adequately represents the input, such that
C	SQRT((dX/XRANGE)**2+(dY/YRANGE)**2) < XE
C

	REAL*4 X(*),Y(*),XOUT(*),YOUT(*)
	REAL*4 KX,KY,JX,JY,LX,LY
	INTEGER NOW,DONE,IPV(3,3)
	REAL*4 XE

	IF(NIN.EQ.0) RETURN
	N=NIN
C
	XOUT(1)=X(1)
	YOUT(1)=Y(1)
	KX=X(1)/XRANGE
	KY=Y(1)/YRANGE
	NOUT=1
	DONE=1
	IWORST=0
	DO 200 NOW=2,N
	    NN=NOW
	    JX=X(NOW)/XRANGE
	    JY=Y(NOW)/YRANGE
C
C	IF CURRENT POINT EQUALS LAST PLOTTED POINT, CHOOSE ARBITRARY DIRECTION
C
	    IF(JX.EQ.KX .AND. JY.EQ.KY) THEN
		DX=1.
		DY=0.
		D=0.
	    ELSE
		DX=JX-KX
		DY=JY-KY
		DD=DX*DX+DY*DY
		D=SQRT(DD)
		DX=DX/D
		DY=DY/D
	    END IF
	    NTEST=NOW-DONE-1
	    NSKIP=1
	    NHIGH=NOW-1
	    NLOW=DONE+1
	    IF(NTEST.GT.5) THEN
		NSKIP=SQRT(FLOAT(NTEST))
		NSKIP=MAX(NSKIP,NTEST/4+1)
		NLOW=NHIGH-NSKIP*(NTEST/NSKIP)
		IF(IWORST.NE.0) THEN
		    LX=X(IWORST)/XRANGE
		    LY=Y(IWORST)/YRANGE
		    EX=LX-KX
		    EY=LY-KY
		    EPS=ABS(EX*DY-EY*DX)
		    IF(EPS .GE. XE) GO TO 70
		    PROGRESS=EX*DX+EY*DY
		    IF(PROGRESS.LT.-XE .OR. PROGRESS.GT.D+XE) GO TO 70
		ENDIF
	    ENDIF
	    DO I=NLOW,NHIGH,NSKIP
		LX=X(I)/XRANGE
		LY=Y(I)/YRANGE
		EX=LX-KX
		EY=LY-KY
		EPS=ABS(EX*DY-EY*DX)
		IF(EPS.GT.WORST) THEN
		    WORST=EPS
		    IWORST=I
		ENDIF
		IF(EPS .GE. XE) GO TO 70
		PROGRESS=EX*DX+EY*DY
		IF(PROGRESS.LT.-XE .OR. PROGRESS.GT.D+XE) GO TO 70
	    END DO
	    GO TO 100
C
C	LINEARITY TEST FAILED, SO PLOT UP TO THE PREVIOUS POINT (NOW-1)
C
70	    LX=X(NN-1)/XRANGE
	    LY=Y(NN-1)/YRANGE
	    IWORST=0
	    WORST=0.
	    NOUT=NOUT+1
	    XOUT(NOUT)=X(NN-1)
	    YOUT(NOUT)=Y(NN-1)
	    KX=LX
	    KY=LY
	    DONE=NN-1
C
C	CHECK TO FORCE PLOTTING OF LAST POINT
C
100	    IF(NOW.EQ.N .AND. DONE.LT.N) THEN
		NOUT=NOUT+1
		XOUT(NOUT)=X(N)
		YOUT(NOUT)=Y(N)
		DONE=N
	    END IF
200	CONTINUE
	RETURN
	END
