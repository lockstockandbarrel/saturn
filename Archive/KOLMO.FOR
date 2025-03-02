C                                                                       KLMO  10
C     ..................................................................KLMO  20
C                                                                       KLMO  30
C        SUBROUTINE KOLMO                                               KLMO  40
C                                                                       KLMO  50
C        PURPOSE                                                        KLMO  60
C           TESTS THE DIFFERENCE BETWEEN EMPIRICAL AND THEORETICAL      KLMO  70
C           DISTRIBUTIONS  USING THE KOLMOGOROV-SMIRNOV TEST            KLMO  80
C                                                                       KLMO  90
C        USAGE                                                          KLMO 100
C           CALL KOLMO(X,N,Z,PROB,IFCOD,U,S,IER)                        KLMO 110
C                                                                       KLMO 120
C        DESCRIPTION OF PARAMETERS                                      KLMO 130
C           X    - INPUT VECTOR OF N INDEPENDENT OBSERVATIONS.  ON      KLMO 140
C                  RETURN FROM KOLMO, X HAS BEEN SORTED INTO A          KLMO 150
C                  MONOTONIC NON-DECREASING SEQUENCE.                   KLMO 160
C           N    - NUMBER OF OBSERVATIONS IN X                          KLMO 170
C           Z    - OUTPUT VARIABLE CONTAINING THE GREATEST VALUE WITH   KLMO 180
C                  RESPECT TO X OF  SQRT(N)*ABS(FN(X)-F(X)) WHERE       KLMO 190
C                  F(X) IS A  THEORETICAL DISTRIBUTION FUNCTION AND     KLMO 200
C                  FN(X) AN EMPIRICAL DISTRIBUTION FUNCTION.            KLMO 210
C           PROB - OUTPUT VARIABLE CONTAINING THE PROBABILITY OF        KLMO 220
C                  THE STATISTIC BEING GREATER THAN OR EQUAL TO Z IF    KLMO 230
C                  THE HYPOTHESIS THAT X IS FROM THE DENSITY UNDER      KLMO 240
C                  CONSIDERATION IS TRUE.  E.G., PROB = 0.05 IMPLIES    KLMO 250
C                  THAT ONE CAN REJECT THE NULL HYPOTHESIS THAT THE SET KLMO 260
C                  X IS FROM THE DENSITY UNDER CONSIDERATION WITH 5 PER KLMO 270
C                  CENT PROBABILITY OF BEING INCORRECT.  PROB = 1. -    KLMO 280
C                  SMIRN(Z).                                            KLMO 290
C           IFCOD- A CODE DENOTING THE PARTICULAR THEORETICAL           KLMO 300
C                  PROBABILITY DISTRIBUTION FUNCTION BEING CONSIDERED.  KLMO 310
C                  = 1---F(X) IS THE NORMAL PDF.                        KLMO 320
C                  = 2---F(X) IS THE EXPONENTIAL PDF.                   KLMO 330
C                  = 3---F(X) IS THE CAUCHY PDF.                        KLMO 340
C                  = 4---F(X) IS THE UNIFORM PDF.                       KLMO 350
C                  = 5---F(X) IS USER SUPPLIED.                         KLMO 360
C           U    - WHEN IFCOD IS 1 OR 2, U IS THE MEAN OF THE DENSITY   KLMO 370
C                  GIVEN ABOVE.                                         KLMO 380
C                  WHEN IFCOD IS 3, U IS THE MEDIAN OF THE CAUCHY       KLMO 390
C                  DENSITY.                                             KLMO 400
C                  WHEN IFCOD IS 4, U IS THE LEFT ENDPOINT OF THE       KLMO 410
C                  UNIFORM DENSITY.                                     KLMO 420
C                  WHEN IFCOD IS 5, U IS USER SPECIFIED.                KLMO 430
C           S    - WHEN IFCOD IS 1 OR 2, S IS THE STANDARD DEVIATION OF KLMO 440
C                  DENSITY GIVEN ABOVE, AND SHOULD BE POSITIVE.         KLMO 450
C                  WHEN IFCOD IS 3, U - S SPECIFIES THE FIRST QUARTILE  KLMO 460
C                  OF THE CAUCHY DENSITY.  S SHOULD BE NON-ZERO.        KLMO 470
C                  IF IFCOD IS 4, S IS THE RIGHT ENDPOINT OF THE UNIFORMKLMO 480
C                  DENSITY.  S SHOULD BE GREATER THAN U.                KLMO 490
C                  IF IFCOD IS 5, S IS USER SPECIFIED.                  KLMO 500
C           IER  - ERROR INDICATOR WHICH IS NON-ZERO IF S VIOLATES ABOVEKLMO 510
C                  CONVENTIONS.  ON RETURN NO TEST HAS BEEN MADE, AND X KLMO 520
C                  AND Y HAVE BEEN SORTED INTO MONOTONIC NON-DECREASING KLMO 530
C                  SEQUENCES.  IER IS SET TO ZERO ON ENTRY TO KOLMO.    KLMO 540
C                  IER IS CURRENTLY SET TO ONE IF THE USER-SUPPLIED PDF KLMO 550
C                  IS REQUESTED FOR TESTING.  THIS SHOULD BE CHANGED    KLMO 560
C                  (SEE REMARKS) WHEN SOME PDF IS SUPPLIED BY THE USER. KLMO 570
C                                                                       KLMO 580
C        REMARKS                                                        KLMO 590
C           N SHOULD BE GREATER THAN OR EQUAL TO 100.  (SEE THE         KLMO 600
C           MATHEMATICAL DESCRIPTION GIVEN FOR THE PROGRAM SMIRN,       KLMO 610
C           CONCERNING ASYMPTOTIC FORMULAE)  ALSO, PROBABILITY LEVELS   KLMO 620
C           DETERMINED BY THIS PROGRAM WILL NOT BE CORRECT IF THE       KLMO 630
C           SAME SAMPLES ARE USED TO ESTIMATE PARAMETERS FOR THE        KLMO 640
C           CONTINUOUS DISTRIBUTIONS WHICH ARE USED IN THIS TEST.       KLMO 650
C           (SEE THE MATHEMATICAL DESCRIPTION FOR THIS PROGRAM)         KLMO 660
C           F(X) SHOULD BE A CONTINUOUS FUNCTION.                       KLMO 670
C           ANY USER SUPPLIED CUMULATIVE PROBABILITY DISTRIBUTION       KLMO 680
C           FUNCTION SHOULD BE CODED BEGINNING WITH STATEMENT 26 BELOW, KLMO 690
C           AND SHOULD RETURN TO STATEMENT 27.                          KLMO 700
C                                                                       KLMO 710
C           DOUBLE PRECISION USAGE---IT IS DOUBTFUL THAT THE USER WILL  KLMO 720
C           WISH TO PERFORM THIS TEST USING DOUBLE PRECISION ACCURACY.  KLMO 730
C           IF ONE WISHES TO COMMUNICATE WITH KOLMO IN A DOUBLE         KLMO 740
C           PRECISION PROGRAM, HE SHOULD CALL THE FORTRAN SUPPLIED      KLMO 750
C           PROGRAM SNGL(X) PRIOR TO CALLING KOLMO, AND CALL THE        KLMO 760
C           FORTRAN SUPPLIED PROGRAM DBLE(X) AFTER EXITING FROM KOLMO.  KLMO 770
C           (NOTE THAT SUBROUTINE SMIRN DOES HAVE DOUBLE PRECISION      KLMO 780
C           CAPABILITY AS SUPPLIED BY THIS PACKAGE.)                    KLMO 790
C                                                                       KLMO 800
C                                                                       KLMO 810
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  KLMO 820
C           SMIRN, NDTR, AND ANY USER SUPPLIED SUBROUTINES REQUIRED.    KLMO 830
C                                                                       KLMO 840
C        METHOD                                                         KLMO 850
C           FOR REFERENCE, SEE (1) W. FELLER--ON THE KOLMOGOROV-SMIRNOV KLMO 860
C           LIMIT THEOREMS FOR EMPIRICAL DISTRIBUTIONS--                KLMO 870
C           ANNALS OF MATH. STAT., 19, 1948.  177-189,                  KLMO 880
C           (2) N. SMIRNOV--TABLE FOR ESTIMATING THE GOODNESS OF FIT    KLMO 890
C           OF EMPIRICAL DISTRIBUTIONS--ANNALS OF MATH. STAT., 19,      KLMO 900
C           1948.  279-281.                                             KLMO 910
C           (3) R. VON MISES--MATHEMATICAL THEORY OF PROBABILITY AND    KLMO 920
C           STATISTICS--ACADEMIC PRESS, NEW YORK, 1964.  490-493,       KLMO 930
C           (4) B.V. GNEDENKO--THE THEORY OF PROBABILITY--CHELSEA       KLMO 940
C           PUBLISHING COMPANY, NEW YORK, 1962.  384-401.               KLMO 950
C                                                                       KLMO 960
C     ..................................................................KLMO 970
C                                                                       KLMO 980
      SUBROUTINE KOLMO(X,N,Z,PROB,IFCOD,U,S,IER)                        KLMO 990
      DIMENSION X(1)                                                    KLMO1000
C                                                                       KLMO1010
C          NON DECREASING ORDERING OF X(I)'S  (DUBY METHOD)             KLMO1020
C                                                                       KLMO1030
      IER=0                                                             KLMO1040
      DO 5 I=2,N                                                        KLMO1050
      IF(X(I)-X(I-1))1,5,5                                              KLMO1060
    1 TEMP=X(I)                                                         KLMO1070
      IM=I-1                                                            KLMO1080
      DO 3 J=1,IM                                                       KLMO1090
      L=I-J                                                             KLMO1100
      IF(TEMP-X(L))2,4,4                                                KLMO1110
    2 X(L+1)=X(L)                                                       KLMO1120
    3 CONTINUE                                                          KLMO1130
      X(1)=TEMP                                                         KLMO1140
      GO TO 5                                                           KLMO1150
    4 X(L+1)=TEMP                                                       KLMO1160
    5 CONTINUE                                                          KLMO1170
C                                                                       KLMO1180
C           COMPUTES MAXIMUM DEVIATION DN IN ABSOLUTE VALUE BETWEEN     KLMO1190
C           EMPIRICAL AND THEORETICAL DISTRIBUTIONS                     KLMO1200
C                                                                       KLMO1210
      NM1=N-1                                                           KLMO1220
      XN=N                                                              KLMO1230
      DN=0.0                                                            KLMO1240
      FS=0.0                                                            KLMO1250
      IL=1                                                              KLMO1260
    6 DO 7  I=IL,NM1                                                    KLMO1270
      J=I                                                               KLMO1280
      IF(X(J)-X(J+1))9,7,9                                              KLMO1290
    7 CONTINUE                                                          KLMO1300
    8 J=N                                                               KLMO1310
    9 IL=J+1                                                            KLMO1320
      FI=FS                                                             KLMO1330
      FS=FLOAT(J)/XN                                                    KLMO1340
      IF(IFCOD-2)10,13,17                                               KLMO1350
   10 IF(S)11,11,12                                                     KLMO1360
   11 IER=1                                                             KLMO1370
      GO TO 29                                                          KLMO1380
   12 Z =(X(J)-U)/S                                                     KLMO1390
      CALL NDTR(Z,Y,D)                                                  KLMO1400
      GO TO 27                                                          KLMO1410
   13 IF(S)11,11,14                                                     KLMO1420
   14 Z=(X(J)-U)/S+1.0                                                  KLMO1430
      IF(Z)15,15,16                                                     KLMO1440
   15 Y=0.0                                                             KLMO1450
      GO TO 27                                                          KLMO1460
   16 Y=1.-EXP(-Z)                                                      KLMO1470
      GO TO 27                                                          KLMO1480
   17 IF(IFCOD-4)18,20,26                                               KLMO1490
   18 IF(S)19,11,19                                                     KLMO1500
   19 Y=ATAN((X(J)-U)/S)*0.3183099+0.5                                  KLMO1510
      GO TO 27                                                          KLMO1520
   20 IF(S-U)11,11,21                                                   KLMO1530
   21 IF(X(J)-U)22,22,23                                                KLMO1540
   22 Y=0.0                                                             KLMO1550
      GO TO 27                                                          KLMO1560
   23 IF(X(J)-S)25,25,24                                                KLMO1570
   24 Y=1.0                                                             KLMO1580
      GO TO 27                                                          KLMO1590
   25 Y=(X(J)-U)/(S-U)                                                  KLMO1600
      GO TO 27                                                          KLMO1610
   26 IER=1                                                             KLMO1620
      GO TO 29                                                          KLMO1630
   27 EI=ABS(Y-FI)                                                      KLMO1640
      ES=ABS(Y-FS)                                                      KLMO1650
      DN=AMAX1(DN,EI,ES)                                                KLMO1660
      IF(IL-N)6,8,28                                                    KLMO1670
C                                                                       KLMO1680
C           COMPUTES Z=DN*SQRT(N)  AND  PROBABILITY                     KLMO1690
C                                                                       KLMO1700
   28 Z=DN*SQRT(XN)                                                     KLMO1710
      CALL SMIRN(Z,PROB)                                                KLMO1720
      PROB=1.0-PROB                                                     KLMO1730
   29 RETURN                                                            KLMO1740
      END                                                               KLMO1750
