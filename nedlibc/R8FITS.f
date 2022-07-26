	REAL*8 FUNCTION R8FITS(BBBB)
C
C	converts 8 bytes from fits file into a real*8 value
C	on Power PC Macs this does not need  a byte order reversal
C
	BYTE BBBB(8),LOCAL(8)
	REAL*8 VALUE
	EQUIVALENCE (LOCAL,VALUE)
C       byte order testing setup
	BYTE BTEST(4) /1,2,3,4/
	INTEGER*4 JTEST
	EQUIVALENCE (JTEST,BTEST)
C
        IF (JTEST.EQ.16909060) THEN
C       Power PC CPU
C
	  DO I=1,8
	    LOCAL(I) = BBBB(I)
	  ENDDO
	ELSE IF (JTEST.EQ.67305985) THEN
C       Intel CPU
	  DO I=1,8
	    LOCAL(I) = BBBB(9-I)
	  ENDDO
	ELSE
	  STOP 'Bad JTEST'
	ENDIF
	R8FITS = VALUE
	RETURN
	END
