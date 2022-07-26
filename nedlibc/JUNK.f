	SUBROUTINE JUNK(X)
C	force the &*%$# compiler to flush the registers!
	REAL*8 X,Y
	COMMON /JNKCMN/ Y
	Y = X+1
	RETURN
	END
