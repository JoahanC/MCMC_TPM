      SUBROUTINE SPLINSET(X,Y,N,YP1,YP2,Y2)
      PARAMETER (NMAX=2000)
      DIMENSION Y2(N),X(N),Y(N),U(NMAX)
      IF (N.GT.NMAX) WRITE (*,*) 'N too big in SPLINSET'
      DO 5 I=2,N
        IF(X(I).LE.X(I-1)) WRITE (*,666)
666   FORMAT(' X ARRAY NOT IN ORDER FOR SPLINSET')
5     CONTINUE
      IF(YP1.GT.0.99E30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))/
     1     (X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF(YP2.GT.0.99E30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YP2-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END
