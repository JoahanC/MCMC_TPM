      SUBROUTINE SPLOPE(NP,X,Y,TAB)
C
C	SETS UP FOR SPLINE EVALUATION USING OLD FORM OF DATA FORMAT
C	TAB IS DIM. (3,NP) AND TAB(1,1) HAS Y'(1)
C	TAB(1,NP) HAS Y'(NP)
C
      DIMENSION TAB(3,NP),X(NP),Y(NP)
      YP1=TAB(1,1)
      YPN=TAB(1,NP)
      CALL SPLINSET(X,Y,NP,YP1,YPN,TAB)
      RETURN
      END
