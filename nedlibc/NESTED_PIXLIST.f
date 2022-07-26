	SUBROUTINE NESTED_PIXLIST(G,RADIUS,IRES,NL,LIST)
C
C	G = unit vector
C	RADIUS in radians
C	NL is the number of IRES pixels within RADIUS of G
C	and LIST is the list of pixel numbers
C
	REAL*4 G(3),RADIUS
	INTEGER IRES,NL,LIST(*)
C
	REAL*4 E(3),N(3),C(3)
	INTEGER JRES,NJ,LJ(9999)
	JRES = ALOG(2./RADIUS)/ALOG(2.)
	JRES = MIN(JRES,IRES)
	SPACING = 0.35/2**JRES
	NDO = RADIUS/SPACING + 2.0
	N(1) = 0
	N(2) = 0
	N(3) = 1
	NJ = 0
	CALL ORTH2(G,N,E)
	DO IX = -NDO,NDO
	  DO IY = -NDO,NDO
	    DO K = 1,3
	      C(K) = G(K)+SPACING*(IX*E(K)+IY*N(K))
	    ENDDO
	    CALL NORM(C)
	    CALL NESTED_PIXNO(C,JRES,JP)
	    DO J=1,NJ
	      IF (JP.EQ.LJ(J)) GO TO 50
	    ENDDO
	    NJ = NJ+1
	    LJ(NJ) = JP
50	  ENDDO
	ENDDO
C
	RR = (2.*SIN(RADIUS/2.))**2
	IMUL = 4**(IRES-JRES)
	NL = 0
	DO J=1,NJ
	  KL = IMUL*LJ(J)
	  KH = KL+IMUL-1
	  DO K = KL,KH
	    CALL NESTED_CENPIX(K,IRES,C)
	    IF (DVS(G,C).LT.RR) THEN
	      NL = NL+1
	      LIST(NL) = K
	    ENDIF
	  ENDDO
	ENDDO
	RETURN
	END

