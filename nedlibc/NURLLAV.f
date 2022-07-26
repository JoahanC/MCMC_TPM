C   IMSL ROUTINE NAME   - RLLAV
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - VAX/SINGLE
C
C   LATEST REVISION     - JUNE 1, 1980
C
C   PURPOSE             - PERFORM LINEAR REGRESSION USING
C                           THE LEAST ABSOLUTE VALUES CRITERION
C
C   USAGE               - CALL RLLAV (XY,IXY,N,M,IOPT,BETA,
C                           SUMRE,ITER,IRANK,IWK,WK,IER)
C
C   ARGUMENTS    XY     - INPUT/OUTPUT.  ON INPUT XY IS AN N BY
C                           (M+5) MATRIX CONTAINING IN THE FIRST
C                           M COLUMNS THE VALUES OF THE INDEPENDENT
C                           VARIABLES AND IN THE (M+1)TH COLUMN
C                           THE VALUES OF THE DEPENDENT VARIABLE.
C                           THE LAST FOUR COLUMNS ARE WORKSPACE.
C                           ON OUTPUT, THE RESIDUALS ARE STORED IN
C                           THE (M+5)TH COLUMN OF XY.  THE REST OF
C                           XY IS DESTROYED.
C                IXY    - INPUT ROW DIMENSION OF XY EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           IN THE CALLING PROGRAM.
C                N      - INPUT NUMBER OF OBSERVATIONS.
C                M      - INPUT NUMBER OF INDEPENDENT VARIABLES.
C                IOPT   - INPUT OPTION PARAMETER SPECIFYING, IF
C                           IOPT = 0, THAT AN INTERCEPT TERM IS
C                           TO BE INCLUDED IN THE MODEL, OTHERWISE
C                           NO INTERCEPT TERM IS INCLUDED.
C                BETA   - OUTPUT VECTOR CONTAINING THE ESTIMATED
C                           COEFFICIENTS.  IF IOPT = 0, BETA IS
C                           OF LENGTH M+1 WITH THE INTERCEPT
C                           ESTIMATE IN THE (M+1)TH POSITION,
C                           OTHERWISE BETA IS OF LENGTH M.
C                SUMRE  - OUTPUT SUM OF THE ABSOLUTE VALUES OF
C                           THE RESIDUALS.
C                ITER   - OUTPUT NUMBER OF ITERATIONS REQUIRED
C                           FOR SOLUTION.
C                IRANK  - OUTPUT RANK OF THE MATRIX OF INDEPENDENT
C                           VARIABLES (INCLUDING THE CONSTANT COLUMN
C                           FOR THE INTERCEPT IF IOPT=0).
C                IWK    - WORK VECTOR OF LENGTH N.
C                WK     - WORK VECTOR OF LENGTH 2*M+4.
C                IER    - ERROR PARAMETER.  (OUTPUT)
C                           WARNING ERROR
C                             IER = 33 INDICATES THE SOLUTION
C                             OBTAINED IS LIKELY TO BE NONUNIQUE.
C                           TERMINAL ERROR
C                             IER = 129 INDICATES THAT CALCULATION
C                             TERMINATED PREMATURELY DUE TO
C                             ROUNDING ERRORS.
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   REQD. IMSL ROUTINES - SINGLE/UERTST,UGETIO
C                       - DOUBLE/UERTST,UGETIO,VXADD,VXSTO
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   COPYRIGHT           - 1980 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C        
C-----------------------------------------------------------------------
C
      SUBROUTINE NURLLAV  (XY,IXY,N,M,IOPT,BETA,SUMRE,ITER,IRANK,IWK,
     1                   WK,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            IXY,N,M,IOPT,ITER,IRANK,IWK(1),IER
      REAL*4             XY(IXY,1),BETA(1),SUMRE,WK(1)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            IN,IOUT,I,J,KL,KOUNT,KR,K,L,M1,M2,M3,M5,MM,N1,
     1                   N2
      REAL*4             AMAX,AMIN,BIG,D,FMM,ONE,PIVOT,TOLER,TWO,ZERO
      LOGICAL            STAGE,TEST
      REAL*8             SUM
      DATA               TOLER/1.1921E-07/
      DATA               BIG/1.7E+38/
      DATA               ONE/1.0/,ZERO/0.0/,TWO/2.0/
C                                  FIRST EXECUTABLE STATEMENT
      MM = M
      N1 = N+1
      N2 = N+2
      M1 = M+1
      M2 = M+2
      M5 = M+5
      IF (IOPT.NE.0) GO TO 10
      MM = M+1
      DO 5 I=1,N            
         XY(I,M2) = XY(I,M1)
         XY(I,M1) = ONE
    5 CONTINUE
      M1 = MM+1
      M2 = MM+2
   10 M3 = MM+3
      FMM = MM
      DO 15 J=1,MM
         WK(M1+J) = J
         BETA(J) = ZERO
   15 CONTINUE
      DO 30 I=1,N
         XY(I,M2) = MM+I
         XY(I,M3) = XY(I,M1)
         IF (XY(I,M3).GE.ZERO) GO TO 25
         DO 20 J=1,M2
            XY(I,J) = -XY(I,J)
   20    CONTINUE
   25    XY(I,M5) = ZERO
   30 CONTINUE
C                                  COMPUTE THE MARGINAL COSTS
      DO 40 J=1,M1
         SUM = 0.D0
         DO 35 I=1,N
            SUM = SUM+XY(I,J)
   35    CONTINUE
         WK(J) = SUM
   40 CONTINUE
C                                  STAGE I
C                                    DETERMINE THE VECTOR TO
C                                    ENTER THE BASIS
      STAGE = .TRUE.
      KOUNT = 0
      KR = 1
      KL = 1
   45 AMAX = -1
      DO 50 J=KR,MM
         IF (ABS(WK(M1+J)).GT.FMM) GO TO 50
         D = ABS(WK(J))
         IF (D.LE.AMAX) GO TO 50
         AMAX = D
         IN = J
   50 CONTINUE
      IF (WK(IN).GE.ZERO) GO TO 60
      DO 55 I=1,N
         XY(I,IN) = -XY(I,IN)
   55 CONTINUE
      WK(IN) = -WK(IN)
      WK(M1+IN) = -WK(M1+IN)
C                                  DETERMINE THE VECTOR TO LEAVE THE
C                                    BASIS
   60 K = 0
      DO 65 I=KL,N
         D = XY(I,IN)
         IF (D.LE.TOLER) GO TO 65
         K = K+1
         XY(K,M3) = XY(I,M1)/D
         IWK(K) = I
         TEST = .TRUE.
   65 CONTINUE
   70 IF (K.GT.0) GO TO 75
      TEST = .FALSE.
      GO TO 85
   75 AMIN = BIG
      DO 80 I=1,K
         IF (XY(I,M3).GE.AMIN) GO TO 80
         J = I
         AMIN = XY(I,M3)
         IOUT = IWK(I)
   80 CONTINUE
      XY(J,M3) = XY(K,M3)
      IWK(J) = IWK(K)
      K = K-1
C                                  CHECK FOR LINEAR DEPENDENCE IN
C                                    STAGE I
   85 IF (TEST.OR..NOT.STAGE) GO TO 95
C
C	SWAP KR AND IN, BOTH LE MM
C
      DO 90 I=1,N
         D = XY(I,KR)
         XY(I,KR) = XY(I,IN)
         XY(I,IN) = D
   90 CONTINUE
      D = WK(KR)
      WK(KR) = WK(IN)
      WK(IN) = D
      D = WK(M1+KR)
      WK(M1+KR) = WK(M1+IN)
      WK(M1+IN) = D
      KR = KR+1
      GO TO 145
   95 IF (TEST) GO TO 100
      IER = 129
      GO TO 9000
  100 PIVOT = XY(IOUT,IN)
      IF (WK(IN)-PIVOT-PIVOT.LE.TOLER) GO TO 110
C
C	MODIFY ELEMENTS IN THE IOUT'TH ROW
C
      DO 105 J=KR,M1
         D = XY(IOUT,J)
         WK(J) = WK(J)-D-D
         XY(IOUT,J) = -D
  105 CONTINUE
      XY(IOUT,M2) = -XY(IOUT,M2)
      GO TO 70
C                                  PIVOT ON XY(IOUT,IN)
  110 DO 115 J=KR,M1
         IF (J.EQ.IN) GO TO 115
         XY(IOUT,J) = XY(IOUT,J)/PIVOT
  115 CONTINUE
      DO 125 I=1,N
         IF (I.EQ.IOUT) GO TO 125
         D = XY(I,IN)
         DO 120 J=KR,M1
            IF (J.EQ.IN) GO TO 120
C
C	MIXES I'TH AND IOUT'TH ROWS IN COLUMNS KR - M1
C
            XY(I,J) = XY(I,J)-D*XY(IOUT,J)
  120    CONTINUE
  125 CONTINUE
      DO 130 J=KR,M1
         IF (J.EQ.IN) GO TO 130
         WK(J) = WK(J)-WK(IN)*XY(IOUT,J)
  130 CONTINUE
      DO 135 I=1,N
         IF (I.EQ.IOUT) GO TO 135
         XY(I,IN) = -XY(I,IN)/PIVOT
  135 CONTINUE
      WK(IN) = -WK(IN)/PIVOT
      XY(IOUT,IN) = ONE/PIVOT
      D = XY(IOUT,M2)
      XY(IOUT,M2) = WK(M1+IN)
      WK(M1+IN) = D
      KOUNT = KOUNT+1
      IF (.NOT.STAGE) GO TO 150
C                                  INTERCHANGE ROWS IN STAGE I
      KL = KL+1
C
C	INTERCHANGE ROWS IOUT AND KOUNT
C
      DO 140 J=KR,M2
         D = XY(IOUT,J)
         XY(IOUT,J) = XY(KOUNT,J)
         XY(KOUNT,J) = D
  140 CONTINUE
  145 IF (KOUNT+KR.NE.M1) GO TO 45
C                                  STAGE II
      STAGE = .FALSE.
C                                  DETERMINE THE VECTOR TO ENTER THE
C                                    BASIS
  150 AMAX = -BIG
      DO 160 J=KR,MM
         D = WK(J)
         IF (D.GE.ZERO) GO TO 155
         IF (D.GT.(-TWO)) GO TO 160
         D = -D-TWO
  155    IF (D.LE.AMAX) GO TO 160
         AMAX = D
         IN = J
  160 CONTINUE
      IF (AMAX.LE.TOLER) GO TO 170
      IF (WK(IN).GT.ZERO) GO TO 60
C
C	NEGATE IN'TH COLUMN
C
      DO 165 I=1,N
         XY(I,IN) = -XY(I,IN)
  165 CONTINUE
      WK(IN) = -WK(IN)-TWO
      WK(M1+IN) = -WK(M1+IN)
      GO TO 60
C                                  PREPARE OUTPUT
  170 L = KL-1
      DO 180 I=1,L
         IF (XY(I,M1).GE.ZERO) GO TO 180
C
C	NEGATE I'TH ROW
C
         DO 175 J=KR,M2
            XY(I,J) = -XY(I,J)
  175    CONTINUE
  180 CONTINUE
      IER = 33
      IF (KR.NE.1) GO TO 190
      DO 185 J=1,MM
         D = ABS(WK(J))
         IF (D.LE.TOLER.OR.TWO-D.LE.TOLER) GO TO 190
  185 CONTINUE
      IER = 0
  190 DO 205 I=1,N
C
C	RIGHT HAND SIDE MACHINATIONS
C
         K = XY(I,M2)
         D = XY(I,M1)
         IF (K.GT.0) GO TO 195
         K = -K
         D = -D
  195    IF (I.GE.KL) GO TO 200
         BETA(K) = D
         GO TO 205
  200    K = K-MM
         XY(K,M5) = D
  205 CONTINUE
      ITER = KOUNT
      IRANK = M1-KR
      SUM = 0.D0
      DO 210 I=KL,N
         SUM = SUM+XY(I,M1)
  210 CONTINUE
      SUMRE = SUM
      IF (IER.EQ.0) GO TO 9005
 9000 CONTINUE

C     CALL UERTST (IER,6HRLLAV )	taken out by MBK. Don't want
C  errors logged into a stupid file.
C

 9005 RETURN
      END
