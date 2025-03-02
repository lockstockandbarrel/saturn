C                                                                       TETR  10
C     ..................................................................TETR  20
C                                                                       TETR  30
C        SUBROUTINE TETRA                                               TETR  40
C                                                                       TETR  50
C        PURPOSE                                                        TETR  60
C           COMPUTE A TETRACHORIC CORRELATION COEFFICIENT BETWEEN TWO   TETR  70
C           VARIABLES WHERE DATA IN BOTH VARIABLES HAVE BEEN REDUCED    TETR  80
C           ARTIFICIALLY TO TWO CATEGORIES.                             TETR  90
C                                                                       TETR 100
C        USAGE                                                          TETR 110
C           CALL TETRA (N,U,V,HU,HV,R,RS,IE)                            TETR 120
C                                                                       TETR 130
C        DESCRIPTION OF PARAMETERS                                      TETR 140
C           N  - NUMBER OF OBSERVATIONS                                 TETR 150
C           U  - INPUT VECTOR OF LENGTH N CONTAINING THE FIRST VARIABLE TETR 160
C                REDUCED TO TWO CATEGORIES                              TETR 170
C           V  - INPUT VECTOR OF LENGTH N CONTAINING THE SECOND VARIABLETETR 180
C                REDUCED TO TWO CATEGORIES                              TETR 190
C           HU - INPUT NUMERICAL CODE INDICATING THE HIGHER CATEGORY OF TETR 200
C                THE FIRST VARIABLE.  IF ANY VALUE OF VARIABLE U IS     TETR 210
C                EQUAL TO OR GREATER THAN HU, IT WILL BE CLASSIFIED AS  TETR 220
C                THE HIGHER CATEGORY, OTHERWISE AS THE LOWER CATEGORY.  TETR 230
C           HV - SAME AS HU EXCEPT THAT HV IS FOR THE SECOND VARIABLE.  TETR 240
C           R  - TETRACHORIC CORRELATION COMPUTED                       TETR 250
C           RS - STANDARD ERROR OF TETRACHORIC CORRELATION COMPUTED     TETR 260
C           IE - ERROR CODE                                             TETR 270
C                0 - NO ERROR                                           TETR 280
C                1 - UNABLE TO COMPUTE A TETRACHORIC CORRELATION DUE TO TETR 290
C                    THE FACT THAT AT LEAST ONE CELL SHOWS ZERO FRE-    TETR 300
C                    QUENCY IN THE 2X2 CONTINGENCY TABLE CONSTRUCTED    TETR 310
C                    FROM INPUT DATA.  IN THIS CASE, R AND RS ARE SET   TETR 320
C                    TO 10**75.  (SEE GUILFORD, 1956)                   TETR 330
C                2 - THE ROOT SOLVER GIVES MULTIPLE ROOTS, OR NO ROOTS, TETR 340
C                    R, IN THE INTERVAL (-1,1) INCLUSIVE. R AND RS ARE  TETR 350
C                    SET TO 10**75.                                     TETR 360
C                3 - UNABLE TO COMPUTE A SATISFACTORY VALUE OF TETRA-   TETR 370
C                    CHORIC CORRELATION USING NEWTON-RAPHSON METHOD OF  TETR 380
C                    APPROXIMATION TO THE ROOT OF THE EQUATION.  R AND  TETR 390
C                    RS ARE SET TO 10**75.  SEE SUBROUTINE POLRT ERROR  TETR 400
C                    INDICATORS.                                        TETR 410
C                4 - HIGH ORDER COEFFICIENT OF THE POLYNOMIAL IS ZERO.  TETR 420
C                    SEE SUBROUTINE POLRT ERROR INDICATORS.             TETR 430
C                                                                       TETR 440
C        REMARKS                                                        TETR 450
C           VALUES OF VARIABLES U AND V MUST BE NUMERICAL, AND          TETR 460
C           ALPHABETIC AND SPECIAL CHARACTERS MUST NOT BE USED.         TETR 470
C           FOR A DEPENDABLE RESULT FOR TETRACHORIC CORRELATION,        TETR 480
C           IT IS RECOMMENDED THAT N BE AT LEAST 200 OR GREATER.        TETR 490
C                                                                       TETR 500
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  TETR 510
C           NDTRI                                                       TETR 520
C           POLRT--THIS POLYNOMIAL ROOT ROUTINE WAS SELECTED BECAUSE OF TETR 530
C                  ITS SMALL STORAGE REQUIREMENT.  OTHER SSP ROUTINES   TETR 540
C                  WHICH COULD REPLACE POLRT ARE PRQD AND PRBM.  THEIR  TETR 550
C                  USE WOULD REQUIRE MODIFICATION OF TETRA.             TETR 560
C                                                                       TETR 570
C        METHOD                                                         TETR 580
C           REFER TO J. P. GUILFORD, 'FUNDAMENTAL STATISTICS IN PSYCHO- TETR 590
C           LOGY AND EDUCATION', MCGRAW-HILL, NEW YORK, 1956, CHAPTER 13TETR 600
C           AND W. P. ELDERTON, 'FREQUENCY CURVES AND CORRELATION' 4-TH TETR 610
C           ED., CAMBRIDGE UNIVERSITY PRESS, 1953, CHAPTER 9.           TETR 620
C                                                                       TETR 630
C     ..................................................................TETR 640
C                                                                       TETR 650
      SUBROUTINE TETRA (N,U,V,HU,HV,R,RS,IE)                            TETR 660
C                                                                       TETR 670
      DIMENSION XCOF(8),COF(8),ROOTR(7),ROOTI(7)                        TETR 680
      DIMENSION U(1),V(1)                                               TETR 690
      DOUBLE PRECISION X31,X32,X312,X322                                TETR 700
C                                                                       TETR 710
C        CONSTRUCT A 2X2 CONTINGENCY TABLE                              TETR 720
C                                                                       TETR 730
      A=0.0                                                             TETR 740
      B=0.0                                                             TETR 750
      C=0.0                                                             TETR 760
      D=0.0                                                             TETR 770
      DO 40 I=1,N                                                       TETR 780
      IF(U(I)-HU) 10, 25, 25                                            TETR 790
   10 IF(V(I)-HV) 15, 20, 20                                            TETR 800
   15 D=D+1.0                                                           TETR 810
      GO TO 40                                                          TETR 820
   20 B=B+1.0                                                           TETR 830
      GO TO 40                                                          TETR 840
   25 IF(V(I)-HV) 30, 35, 35                                            TETR 850
   30 C=C+1.0                                                           TETR 860
      GO TO 40                                                          TETR 870
   35 A=A+1.0                                                           TETR 880
   40 CONTINUE                                                          TETR 890
C                                                                       TETR 900
C        TEST WHETHER ANY CELL IN THE CONTINGENCY TABLE IS ZERO.        TETR 910
C        IF SO, RETURN TO THE CALLING ROUTINE WITH R=0.0 AND IE=1.      TETR 920
C                                                                       TETR 930
      IE=0                                                              TETR 940
      IF(A) 60, 60, 45                                                  TETR 950
   45 IF(B) 60, 60, 50                                                  TETR 960
   50 IF(C) 60, 60, 55                                                  TETR 970
   55 IF(D) 60, 60, 70                                                  TETR 980
   60 IE=1                                                              TETR 990
      GO TO 86                                                          TETR1000
C                                                                       TETR1010
C        COMPUTE P1, Q1, P2, AND Q2                                     TETR1020
C                                                                       TETR1030
   70 FN=N                                                              TETR1040
      P1=(A+C)/FN                                                       TETR1050
      Q1=(B+D)/FN                                                       TETR1060
      P2=(A+B)/FN                                                       TETR1070
      Q2=(C+D)/FN                                                       TETR1080
C                                                                       TETR1090
C        FIND THE STANDARD NORMAL DEVIATES AT Q1 AND Q2, AND THE        TETR1100
C        ORDINATES AT THOSE POINTS                                      TETR1110
C                                                                       TETR1120
      CALL NDTRI (Q1,X1,Y1,ER)                                          TETR1130
      CALL NDTRI (Q2,X2,Y2,ER)                                          TETR1140
C                                                                       TETR1150
C        COMPUTE THE TETRACHORIC CORRELATION COEFFICIENT                TETR1160
C                                                                       TETR1170
      IF(X1) 76, 72, 76                                                 TETR1180
   72 IF(X2) 76, 74, 76                                                 TETR1190
   74 R=0.0                                                             TETR1200
      GO TO 90                                                          TETR1210
   76 XCOF(1)=-((A*D-B*C)/(Y1*Y2*FN*FN))                                TETR1220
      XCOF(2)=1.0                                                       TETR1230
      XCOF(3)=X1*X2/2.0                                                 TETR1240
      XCOF(4)=(X1*X1-1.0)*(X2*X2-1.0)/6.0                               TETR1250
      X31=DBLE(X1)                                                      TETR1260
      X32=DBLE(X2)                                                      TETR1270
      X312=X31**2                                                       TETR1280
      X322=X32**2                                                       TETR1290
      XCOF(5)=SNGL(X31*(X312-3.0D0)*X32*(X322-3.0D0)/24.0D0)            TETR1300
      XCOF(6)=SNGL((X312*(X312-6.0D0)+3.0D0)*(X322*(X322-6.0D0)+3.0D0)  TETR1310
     1        /120.0D0)                                                 TETR1320
      XCOF(7)=SNGL(X31*(X312*(X312-10.0D0)+15.0D0)*X32*(X322*(X322-10.0 TETR1330
     1        D0)+15.0D0)/720.0D0)                                      TETR1340
      XCOF(8)=SNGL((((X312-15.0D0)*X312+45.0D0)*X312-15.0D0)*(((X322-   TETR1350
     1        15.0D0)*X322+45.0D0)*X322-15.0D0)/5040.0D0)               TETR1360
C                                                                       TETR1370
      CALL POLRT (XCOF,COF,7,ROOTR,ROOTI,IER)                           TETR1380
C                                                                       TETR1390
      J=0                                                               TETR1400
      IF(IER) 78, 78, 84                                                TETR1410
   78 DO 82 I=1,7                                                       TETR1420
      IF(ABS(ROOTI(I))-.5*ABS(ROOTR(I))*1.0E-6)79,79,82                 TETR1430
   79 R=ROOTR(I)                                                        TETR1440
      IF(ABS(R)-1.0)81,81,80                                            TETR1450
   80 R=1.E75                                                           TETR1460
      GO TO 82                                                          TETR1470
   81 J=J+1                                                             TETR1480
   82 CONTINUE                                                          TETR1490
      IF(J-1)83,88,83                                                   TETR1500
   83 IE=2                                                              TETR1510
      GO TO 86                                                          TETR1520
C                                                                       TETR1530
C        UNABLE TO COMPUTE R                                            TETR1540
C                                                                       TETR1550
   84 IE=IER                                                            TETR1560
   86 R=1.0E75                                                          TETR1570
      RS=R                                                              TETR1580
      GO TO 100                                                         TETR1590
   88 IF(R-1.0E75)90,83,83                                              TETR1600
C                                                                       TETR1610
C        STANDARD ERROR OF R=0.0                                        TETR1620
C                                                                       TETR1630
   90 RS= SQRT(P1*P2*Q1*Q2)/(Y1*Y2* SQRT(FN))                           TETR1640
C                                                                       TETR1650
  100 RETURN                                                            TETR1660
      END                                                               TETR1670
