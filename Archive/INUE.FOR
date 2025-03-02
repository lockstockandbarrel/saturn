C                                                                       INUE  10
C     ..................................................................INUE  20
C                                                                       INUE  30
C        SUBROUTINE INUE                                                INUE  40
C                                                                       INUE  50
C        PURPOSE                                                        INUE  60
C           COMPUTE THE MODIFIED BESSEL FUNCTIONS I FOR ORDERS 1 TO N   INUE  70
C                                                                       INUE  80
C        USAGE                                                          INUE  90
C           CALL INUE(X,N,ZI,RI)                                        INUE 100
C                                                                       INUE 110
C        DESCRIPTION OF PARAMETERS                                      INUE 120
C           X     -GIVEN ARGUMENT OF THE BESSEL FUNCTIONS I             INUE 130
C           N     -GIVEN MAXIMUM ORDER OF BESSEL FUNCTIONS I            INUE 140
C           ZI    -GIVEN VALUE OF BESSEL FUNCTION I OF ORDER ZERO       INUE 150
C                  FOR ARGUMENT X                                       INUE 160
C           RI    -RESULTANT VECTOR OF DIMENSION N, CONTAINING THE      INUE 170
C                  VALUES OF THE FUNCTIONS I FOR ORDERS 1 TO N          INUE 180
C                                                                       INUE 190
C        REMARKS                                                        INUE 200
C           THE VALUE OF ZI MAY BE CALCULATED USING SUBROUTINE I0.      INUE 210
C           USING A DIFFERENT VALUE HAS THE EFFECT THAT ALL VALUES OF   INUE 220
C           BESSEL FUNCTIONS I ARE MULTIPLIED BY THE  FACTOR ZI/I(0,X)  INUE 230
C           WHERE I(0,X) IS THE VALUE OF I FOR ORDER 0 AND ARGUMENT X.  INUE 240
C           THIS MAY BE USED DISADVANTAGEOUSLY IF ONLY THE RATIOS OF I  INUE 250
C           FOR DIFFERENT ORDERS ARE REQUIRED.                          INUE 260
C                                                                       INUE 270
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  INUE 280
C           NONE                                                        INUE 290
C                                                                       INUE 300
C        METHOD                                                         INUE 310
C           THE VALUES ARE OBTAINED USING BACKWARD RECURRENCE RELATION  INUE 320
C           TECHNIQUE. THE RATIO I(N+1,X)/I(N,X) IS OBTAINED FROM A     INUE 330
C           CONTINUED FRACTION.                                         INUE 340
C           FOR REFERENCE SEE                                           INUE 350
C           G. BLANCH,'NUMERICAL EVALUATION OF CONTINUED FRACTIONS',    INUE 360
C           SIAM REVIEW, VOL.6,NO.4,1964,PP.383-421.                    INUE 370
C                                                                       INUE 380
C     ..................................................................INUE 390
C                                                                       INUE 400
      SUBROUTINE INUE(X,N,ZI,RI)                                        INUE 410
      DIMENSION RI(1)                                                   INUE 420
      IF(N)10,10,1                                                      INUE 430
    1 FN=N+N                                                            INUE 440
      Q1=X/FN                                                           INUE 450
      IF(ABS(X)-5.E-4)6,6,2                                             INUE 460
    2 A0=1.                                                             INUE 470
      A1=0.                                                             INUE 480
      B0=0.                                                             INUE 490
      B1=1.                                                             INUE 500
      FI=FN                                                             INUE 510
    3 FI=FI+2.                                                          INUE 520
      AN=FI/ABS(X)                                                      INUE 530
      A=AN*A1+A0                                                        INUE 540
      B=AN*B1+B0                                                        INUE 550
      A0=A1                                                             INUE 560
      B0=B1                                                             INUE 570
      A1=A                                                              INUE 580
      B1=B                                                              INUE 590
      Q0=Q1                                                             INUE 600
      Q1=A/B                                                            INUE 610
      IF(ABS((Q1-Q0)/Q1)-1.E-6)4,4,3                                    INUE 620
    4 IF(X)5,6,6                                                        INUE 630
    5 Q1=-Q1                                                            INUE 640
    6 K=N                                                               INUE 650
    7 Q1=X/(FN+X*Q1)                                                    INUE 660
      RI(K)=Q1                                                          INUE 670
      FN=FN-2.                                                          INUE 680
      K=K-1                                                             INUE 690
      IF(K)8,8,7                                                        INUE 700
    8 FI=ZI                                                             INUE 710
      DO 9 I=1,N                                                        INUE 720
      FI=FI*RI(I)                                                       INUE 730
    9 RI(I)=FI                                                          INUE 740
   10 RETURN                                                            INUE 750
      END                                                               INUE 760
