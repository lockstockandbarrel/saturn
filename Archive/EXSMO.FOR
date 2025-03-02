C                                                                       EXSM  10
C     ..................................................................EXSM  20
C                                                                       EXSM  30
C        SUBROUTINE EXSMO                                               EXSM  40
C                                                                       EXSM  50
C        PURPOSE                                                        EXSM  60
C           TO FIND THE TRIPLE EXPONENTIAL SMOOTHED SERIES S OF THE     EXSM  70
C           GIVEN SERIES X.                                             EXSM  80
C                                                                       EXSM  90
C        USAGE                                                          EXSM 100
C           CALL EXSMO (X,NX,AL,A,B,C,S)                                EXSM 110
C                                                                       EXSM 120
C        DESCRIPTION OF PARAMETERS                                      EXSM 130
C           X     - INPUT VECTOR OF LENGTH NX CONTAINING TIME SERIES    EXSM 140
C                   DATA WHICH IS TO BE EXPONENTIALLY SMOOTHED.         EXSM 150
C           NX    - THE NUMBER OF ELEMENTS IN X.                        EXSM 160
C           AL    - SMOOTHING CONSTANT, ALPHA.  AL MUST BE GREATER THAN EXSM 170
C                   ZERO AND LESS THAN ONE.                             EXSM 180
C           A,B,C - COEFFICIENTS OF THE PREDICTION EQUATION WHERE S IS  EXSM 190
C                   PREDICTED T PERIODS HENCE BY                        EXSM 200
C                                 A + B*T + C*T*T/2.                    EXSM 210
C                   AS INPUT-- IF A=B=C=0, PROGRAM WILL PROVIDE INITIAL EXSM 220
C                   VALUES.  IF AT LEAST ONE OF A,B,C IS NOT ZERO,      EXSM 230
C                   PROGRAM WILL TAKE GIVEN VALUES AS INITIAL VALUES.   EXSM 240
C                   AS OUTPUT-- A,B,C CONTAIN LATEST, UPDATED COEFFI-   EXSM 250
C                   CIENTS OF PREDICTION.                               EXSM 260
C           S     - OUTPUT VECTOR OF LENGTH NX CONTAINING TRIPLE        EXSM 270
C                   EXPONENTIALLY SMOOTHED TIME SERIES.                 EXSM 280
C                                                                       EXSM 290
C        REMARKS                                                        EXSM 300
C           NONE                                                        EXSM 310
C                                                                       EXSM 320
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  EXSM 330
C           NONE                                                        EXSM 340
C                                                                       EXSM 350
C        METHOD                                                         EXSM 360
C           REFER TO R. G. BROWN, 'SMOOTHING, FORECASTING AND PREDICTIONEXSM 370
C           OF DISCRETE TIME SERIES', PRENTICE-HALL, N.J., 1963,        EXSM 380
C           PP. 140 TO 144.                                             EXSM 390
C                                                                       EXSM 400
C     ..................................................................EXSM 410
C                                                                       EXSM 420
      SUBROUTINE EXSMO (X,NX,AL,A,B,C,S)                                EXSM 430
      DIMENSION X(1),S(1)                                               EXSM 440
C                                                                       EXSM 450
C     IF A=B=C=0.0, GENERATE INITIAL VALUES OF A, B, AND C              EXSM 460
C                                                                       EXSM 470
      IF(A) 140, 110, 140                                               EXSM 480
  110 IF(B) 140, 120, 140                                               EXSM 490
  120 IF(C) 140, 130, 140                                               EXSM 500
  130 C=X(1)-2.0*X(2)+X(3)                                              EXSM 510
      B=X(2)-X(1)-1.5*C                                                 EXSM 520
      A=X(1)-B-0.5*C                                                    EXSM 530
C                                                                       EXSM 540
  140 BE=1.0-AL                                                         EXSM 550
      BECUB=BE*BE*BE                                                    EXSM 560
      ALCUB=AL*AL*AL                                                    EXSM 570
C                                                                       EXSM 580
C     DO THE FOLLOWING FOR I=1 TO NX                                    EXSM 590
C                                                                       EXSM 600
      DO 150 I=1,NX                                                     EXSM 610
C                                                                       EXSM 620
C        FIND S(I) FOR ONE PERIOD AHEAD                                 EXSM 630
C                                                                       EXSM 640
      S(I)=A+B+0.5*C                                                    EXSM 650
C                                                                       EXSM 660
C        UPDATE COEFFICIENTS A, B, AND C                                EXSM 670
C                                                                       EXSM 680
      DIF=S(I)-X(I)                                                     EXSM 690
      A=X(I)+BECUB*DIF                                                  EXSM 700
      B=B+C-1.5*AL*AL*(2.0-AL)*DIF                                      EXSM 710
  150 C=C-ALCUB*DIF                                                     EXSM 720
      RETURN                                                            EXSM 730
      END                                                               EXSM 740
