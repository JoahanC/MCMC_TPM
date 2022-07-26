      FUNCTION SPLINE(Z,NP,X,Y,TAB)
C
C	COMPUTES SPLINE INTERPOLATION AFTER SPLSET HAS SETUP TAB
C	OR USE THE Y2 ARRAY RETURNED BY SPLINSET
C
      DIMENSION X(NP),Y(NP),TAB(NP)
      DATA NT/0/
      IF(Z.LE.X(1)) GO TO 40
      IF(Z.GE.X(NP)) GO TO 50
      NL=1
      NH=NP-1
      IF(NT.LT.NL .OR. NT.GT.NH) GO TO 17
C
C	BINARY SEARCH OF X FOR X(NT)<Z<X(NT+1) BUT START AT OLD NT
C
14    IF(Z.GT.X(NT+1)) GO TO 15
      IF(Z.LT.X(NT)) GO TO 16
      GO TO 30
15    NL=NT+1
      GO TO 17
16    NH=NT-1
17    NT=(NL+NH)/2
      IF(NL.LT.NH) GO TO 14
30    H=X(NT+1)-X(NT)
      A=(X(NT+1)-Z)/H
      B=(Z-X(NT))/H
      SPLINE=A*Y(NT)+B*Y(NT+1)+
     1   ((A**3-A)*TAB(NT)+(B**3-B)*TAB(NT+1))*(H**2)/6.0
      RETURN
40    NT=1
      D=(Y(2)-Y(1))/(X(2)-X(1))-(X(2)-X(1))*(2.0*TAB(1)+TAB(2))/6.0
      SPLINE=Y(1)+(Z-X(1))*D
      RETURN
50    NT=NP-1
      D=(Y(NP)-Y(NT))/(X(NP)-X(NT))+(X(NP)-X(NT))*
     1   (2.0*TAB(NP)+TAB(NT))/6.0
      SPLINE=Y(NP)+(Z-X(NP))*D
      RETURN
      END
