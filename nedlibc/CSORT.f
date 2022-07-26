      SUBROUTINE CSORT(N,RA)
C	CSORT.FOR
C	heap sort of N complex numbers in RA, with real part as the key
      COMPLEX RA(N),RRA
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(REAL(RA(J)).LT.REAL(RA(J+1)))J=J+1
          ENDIF
          IF(REAL(RRA).LT.REAL(RA(J)))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END

