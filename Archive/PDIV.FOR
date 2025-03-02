C                                                                       PDIV  10
C     ..................................................................PDIV  20
C                                                                       PDIV  30
C        SUBROUTINE PDIV                                                PDIV  40
C                                                                       PDIV  50
C        PURPOSE                                                        PDIV  60
C           DIVIDE ONE POLYNOMIAL BY ANOTHER                            PDIV  70
C                                                                       PDIV  80
C        USAGE                                                          PDIV  90
C           CALL PDIV(P,IDIMP,X,IDIMX,Y,IDIMY,TOL,IER)                  PDIV 100
C                                                                       PDIV 110
C        DESCRIPTION OF PARAMETERS                                      PDIV 120
C           P     - RESULTANT VECTOR OF INTEGRAL PART                   PDIV 130
C           IDIMP - DIMENSION OF P                                      PDIV 140
C           X     - VECTOR OF COEFFICIENTS FOR DIVIDEND POLYNOMIAL,     PDIV 150
C                   ORDERED FROM SMALLEST TO LARGEST POWER. IT IS       PDIV 160
C                   REPLACED BY REMAINDER AFTER DIVISION.               PDIV 170
C           IDIMX - DIMENSION OF X                                      PDIV 180
C           Y     - VECTOR OF COEFFICIENTS FOR DIVISOR POLYNOMIAL,      PDIV 190
C                   ORDERED FROM SMALLEST TO LARGEST POWER              PDIV 200
C           IDIMY - DIMENSION OF Y                                      PDIV 210
C           TOL   - TOLERANCE VALUE BELOW WHICH COEFFICIENTS ARE        PDIV 220
C                   ELIMINATED DURING NORMALIZATION                     PDIV 230
C           IER   - ERROR CODE. 0 IS NORMAL, 1 IS FOR ZERO DIVISOR      PDIV 240
C                                                                       PDIV 250
C        REMARKS                                                        PDIV 260
C           THE REMAINDER R REPLACES X.                                 PDIV 270
C           THE DIVISOR Y REMAINS UNCHANGED.                            PDIV 280
C           IF DIMENSION OF Y EXCEEDS DIMENSION OF X, IDIMP IS SET TO   PDIV 290
C           ZERO AND CALCULATION IS BYPASSED                            PDIV 300
C                                                                       PDIV 310
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PDIV 320
C           PNORM                                                       PDIV 330
C                                                                       PDIV 340
C        METHOD                                                         PDIV 350
C           POLYNOMIAL X IS DIVIDED BY POLYNOMIAL Y GIVING INTEGER PART PDIV 360
C           P AND REMAINDER R SUCH THAT X = P*Y + R.                    PDIV 370
C           DIVISOR Y AND REMAINDER  VECTOR GET NORMALIZED.             PDIV 380
C                                                                       PDIV 390
C     ..................................................................PDIV 400
C                                                                       PDIV 410
      SUBROUTINE PDIV(P,IDIMP,X,IDIMX,Y,IDIMY,TOL,IER)                  PDIV 420
      DIMENSION P(1),X(1),Y(1)                                          PDIV 430
C                                                                       PDIV 440
      CALL PNORM (Y,IDIMY,TOL)                                          PDIV 450
      IF(IDIMY) 50,50,10                                                PDIV 460
   10 IDIMP=IDIMX-IDIMY+1                                               PDIV 470
      IF(IDIMP) 20,30,60                                                PDIV 480
C                                                                       PDIV 490
C     DEGREE OF DIVISOR WAS GREATER THAN DEGREE OF DIVIDEND             PDIV 500
C                                                                       PDIV 510
   20 IDIMP=0                                                           PDIV 520
   30 IER=0                                                             PDIV 530
   40 RETURN                                                            PDIV 540
C                                                                       PDIV 550
C     Y IS ZERO POLYNOMIAL                                              PDIV 560
C                                                                       PDIV 570
   50 IER=1                                                             PDIV 580
      GO TO 40                                                          PDIV 590
C                                                                       PDIV 600
C     START REDUCTION                                                   PDIV 610
C                                                                       PDIV 620
   60 IDIMX=IDIMY-1                                                     PDIV 630
      I=IDIMP                                                           PDIV 640
   70 II=I+IDIMX                                                        PDIV 650
      P(I)=X(II)/Y(IDIMY)                                               PDIV 660
C                                                                       PDIV 670
C     SUBTRACT MULTIPLE OF DIVISOR                                      PDIV 680
C                                                                       PDIV 690
      DO 80 K=1,IDIMX                                                   PDIV 700
      J=K-1+I                                                           PDIV 710
      X(J)=X(J)-P(I)*Y(K)                                               PDIV 720
   80 CONTINUE                                                          PDIV 730
      I=I-1                                                             PDIV 740
      IF(I) 90,90,70                                                    PDIV 750
C                                                                       PDIV 760
C     NORMALIZE REMAINDER POLYNOMIAL                                    PDIV 770
C                                                                       PDIV 780
   90 CALL PNORM(X,IDIMX,TOL)                                           PDIV 790
      GO TO 30                                                          PDIV 800
      END                                                               PDIV 810
