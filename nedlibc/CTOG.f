	SUBROUTINE CTOG(C2000,G)
C
C	CONVERTS VECTOR C FROM CELESTIAL COORDINATES (2000.0)
C	TO GALACTIC COORDINATE VECTOR G
C
	DIMENSION G(3),C2000(3),R(3,3)
	DATA R/ -0.0548809,  0.4941080, -0.8676667,
	1       -0.8734369, -0.4448323, -0.1980717,
 	2       -0.4838349,  0.7469817,  0.4559848/
	CALL MAV(R,C2000,G)
	RETURN
	ENTRY GTOC(G,C2000)
C
C	CONVERTS GALACTIC IN G TO CELESTIAL IN C
C
	CALL VMAT(R,G,C2000)
	RETURN
	END
