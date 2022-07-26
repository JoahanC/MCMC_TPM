	INTEGER*4 FUNCTION I4FITS(BBBB)
C
C	converts 4 blocks from fits file into a real*4 value
C	on Power PC Macs this does not need  a byte order reversal
C
	BYTE BBBB(4),LOCAL(4)
	INTEGER*4 VALUE
	EQUIVALENCE (LOCAL,VALUE)
C       byte order testing setup
	BYTE BTEST(4) /1,2,3,4/
	INTEGER*4 JTEST
	EQUIVALENCE (JTEST,BTEST)
C
        IF (JTEST.EQ.16909060) THEN
C       Power PC CPU
C
	  DO I=1,4
	    LOCAL(I) = BBBB(I)
	  ENDDO
	ELSE IF (JTEST.EQ.67305985) THEN
C       Intel CPU
	  DO I=1,4
	    LOCAL(I) = BBBB(5-I)
	  ENDDO
	ELSE
	  STOP 'Bad JTEST'
	ENDIF
	I4FITS = VALUE
	RETURN
	END
