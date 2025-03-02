C                                                                       SMIR  10
C     ..................................................................SMIR  20
C                                                                       SMIR  30
C        SUBROUTINE SMIRN                                               SMIR  40
C                                                                       SMIR  50
C        PURPOSE                                                        SMIR  60
C           COMPUTES VALUES OF THE LIMITING DISTRIBUTION FUNCTION FOR   SMIR  70
C           THE KOLMOGOROV-SMIRNOV STATISTIC.                           SMIR  80
C                                                                       SMIR  90
C        USAGE                                                          SMIR 100
C           CALL SMIRN(X,Y)                                             SMIR 110
C                                                                       SMIR 120
C        DESCRIPTION OF PARAMETERS                                      SMIR 130
C           X    - THE ARGUMENT OF THE SMIRN FUNCTION                   SMIR 140
C           Y    - THE RESULTANT SMIRN FUNCTION VALUE                   SMIR 150
C                                                                       SMIR 160
C        REMARKS                                                        SMIR 170
C           Y IS SET TO ZERO IF X IS NOT GREATER THAN 0.27, AND IS SET  SMIR 180
C           TO ONE IF X IS NOT LESS THAN 3.1.  ACCURACY TESTS WERE MADE SMIR 190
C           REFERRING TO THE TABLE GIVEN IN THE REFERENCE BELOW.        SMIR 200
C           TWO ARGUMENTS, X= 0.62, AND X = 1.87 GAVE RESULTS WHICH     SMIR 210
C           DIFFER FROM THE SMIRNOV TABLES BY 2.9 AND 1.9 IN THE 5TH    SMIR 220
C           DECIMAL PLACE.  ALL OTHER RESULTS SHOWED SMALLER ERRORS,    SMIR 230
C           AND ERROR SPECIFICATIONS ARE GIVEN IN THE ACCURACY TABLES   SMIR 240
C           IN THIS MANUAL.  IN DOUBLE PRECISION MODE, THESE SAME       SMIR 250
C           ARGUMENTS RESULTED IN DIFFERENCES FROM TABLED VALUES BY 3   SMIR 260
C           AND 2 IN THE 5TH  DECIMAL PLACE.  IT IS NOTED IN            SMIR 270
C           LINDGREN (REFERENCE BELOW) THAT FOR HIGH SIGNIFICANCE LEVELSSMIR 280
C           (SAY, .01 AND .05) ASYMPTOTIC FORMULAS GIVE VALUES WHICH ARESMIR 290
C           TOO HIGH ( BY 1.5 PER CENT WHEN N = 80).  THAT IS, AT HIGH  SMIR 300
C           SIGNIFICANCE LEVELS, THE HYPOTHESIS OF NO DIFFERENCE WILL BESMIR 310
C           REJECTED TOO SELDOM USING ASYMPTOTIC FORMULAS.              SMIR 320
C                                                                       SMIR 330
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  SMIR 340
C           NONE                                                        SMIR 350
C                                                                       SMIR 360
C        METHOD                                                         SMIR 370
C           THE METHOD IS DESCRIBED BY W. FELLER-ON THE KOLMOGOROV-     SMIR 380
C           SMIRNOV LIMIT THEOREMS FOR EMPIRICAL DISTRIBUTIONS- ANNALS  SMIR 390
C           OF MATH. STAT., 19, 1948, 177-189, BY N. SMIRNOV--TABLE     SMIR 400
C           FOR ESTIMATING THE GOODNESS OF FIT OF EMPIRICAL             SMIR 410
C           DISTRIBUTIONS- ANNALS OF MATH. STAT., 19, 1948, 279-281,    SMIR 420
C           AND GIVEN IN LINDGREN, STATISTICAL THEORY, THE MACMILLAN    SMIR 430
C           COMPANY, N. Y., 1962.                                       SMIR 440
C                                                                       SMIR 450
C     ..................................................................SMIR 460
C                                                                       SMIR 470
      SUBROUTINE SMIRN(X,Y)                                             SMIR 480
C     DOUBLE PRECISION X,Q1,Q2,Q4,Q8,Y                                  SMIR 490
C                                                                       SMIR 500
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE CSMIR 510
C        IN COLUMN ONE OF THE DOUBLE PRECISION CARD ABOVE SHOULD BE     SMIR 520
C        REMOVED, AND THE C IN COLUMN ONE OF THE STATEMENTS NUMBERED    SMIR 530
C        C   3, C   5, AND C   8 SHOULD BE REMOVED, AND THESE CARDS     SMIR 540
C        SHOULD REPLACE THE STATEMENTS NUMBERED 3, 5, AND 8,            SMIR 550
C        RESPECTIVELY.  ALL ROUTINES CALLING THIS ROUTINE MUST ALSO     SMIR 560
C        PROVIDE DOUBLE PRECISION ARGUMENTS TO THIS ROUTINE.            SMIR 570
C                                                                       SMIR 580
C     ..................................................................SMIR 590
C                                                                       SMIR 600
      IF(X-.27)1,1,2                                                    SMIR 610
    1 Y=0.0                                                             SMIR 620
      GO TO 9                                                           SMIR 630
    2 IF(X-1.0)3,6,6                                                    SMIR 640
    3 Q1=EXP(-1.233701/X**2)                                            SMIR 650
C   3 Q1=DEXP(-1.233700550136170/X**2)                                  SMIR 660
      Q2=Q1*Q1                                                          SMIR 670
      Q4=Q2*Q2                                                          SMIR 680
      Q8=Q4*Q4                                                          SMIR 690
      IF(Q8-1.0E-25)4,5,5                                               SMIR 700
    4 Q8=0.0                                                            SMIR 710
    5 Y=(2.506628/X)*Q1*(1.0+Q8*(1.0+Q8*Q8))                            SMIR 720
C   5 Y=(2.506628274631001/X)*Q1*(1.0D0+Q8*(1.0D0+Q8*Q8))               SMIR 730
      GO TO 9                                                           SMIR 740
    6 IF(X-3.1)8,7,7                                                    SMIR 750
    7 Y=1.0                                                             SMIR 760
      GO TO 9                                                           SMIR 770
    8 Q1=EXP(-2.0*X*X)                                                  SMIR 780
C   8 Q1=DEXP(-2.0D0*X*X)                                               SMIR 790
      Q2=Q1*Q1                                                          SMIR 800
      Q4=Q2*Q2                                                          SMIR 810
      Q8=Q4*Q4                                                          SMIR 820
      Y=1.0-2.0*(Q1-Q4+Q8*(Q1-Q8))                                      SMIR 830
    9 RETURN                                                            SMIR 840
      END                                                               SMIR 850
