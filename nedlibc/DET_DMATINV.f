	REAL*8 FUNCTION DET_DMATINV(A,NROW,N,NRHS,IFF)
C
C	solve the linear equations A*X=B for the
C	N x N matrix A dimensioned (NROW,N+NRHS)
C	the NRHS right hand sides B are stored in A(*,N+1..N+NRHS)
C	on entry.  On exit they are replaced by the solutions X
C
	REAL*8 A(NROW,*),V,F,S
	REAL*8 DET
	INTEGER NROW,N,NRHS,IFF,M,NP,NM,I,IP,J,II,K
	M=N+NRHS
	NP=N+1
	NM=N-1
	IFF=0
	DET = 1.D0
	DO 40 I=1,NM
C	FIND PIVOT
	  V = ABS(A(I,I))
	  IP=I
	  DO 10 J=I,NM
	    IF(DABS(A(J+1,I)).LE.V) GO TO 10
	      V=DABS(A(J+1,I))
	      IP=J+1
10	  CONTINUE
	  IF(V.LE.0.D0) GO TO 666
	  IF(IP.EQ.I) GO TO 25
C	INTERCHANGE ROWS TO PUT PIVOT ON DIAGONAL
	  DO 20 J=I,M
	    V=A(I,J)
	    A(I,J)=A(IP,J)
20	    A(IP,J)=V
	  DET = -DET
25	  IP=I+1
C	UPDATE DETERMINANT
	  DET = DET*A(I,I)
C	LOOP OVER ROWS BELOW THE DIAGONAL
	  DO 40 IL=IP,N
	    F=-A(IL,I)/A(I,I)
	    A(IL,I)=0.D0
C	REDUCE THE ELEMENTS IN THE ROW
	    DO 40 J=IP,M
40	      A(IL,J)=A(IL,J)+A(I,J)*F
	DET = DET*A(N,N)
C	BACK SUBSTITUTION
	I=NP
	DO 60 II=1,N
	  I=I-1
	  IP=I+1
	  DO 60 J=NP,M
	    S=0.D0
	    IF(I.EQ.N) GO TO 60
	    DO 50 K=IP,N
50	      S=S+A(K,J)*A(I,K)
60	    A(I,J)=(A(I,J)-S)/A(I,I)
	DET_DMATINV = DET
	RETURN
666	IFF=666
	DMATINV = 0
	RETURN
	END
