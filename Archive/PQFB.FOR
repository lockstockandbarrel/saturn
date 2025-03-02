C                                                                       PQFB  10
C     ..................................................................PQFB  20
C                                                                       PQFB  30
C        SUBROUTINE PQFB                                                PQFB  40
C                                                                       PQFB  50
C        PURPOSE                                                        PQFB  60
C           TO FIND AN APPROXIMATION Q(X)=Q1+Q2*X+X*X TO A QUADRATIC    PQFB  70
C           FACTOR OF A GIVEN POLYNOMIAL P(X) WITH REAL COEFFICIENTS.   PQFB  80
C                                                                       PQFB  90
C        USAGE                                                          PQFB 100
C           CALL PQFB(C,IC,Q,LIM,IER)                                   PQFB 110
C                                                                       PQFB 120
C        DESCRIPTION OF PARAMETERS                                      PQFB 130
C           C   - INPUT VECTOR CONTAINING THE COEFFICIENTS OF P(X) -    PQFB 140
C                 C(1) IS THE CONSTANT TERM (DIMENSION IC)              PQFB 150
C           IC  - DIMENSION OF C                                        PQFB 160
C           Q   - VECTOR OF DIMENSION 4 - ON INPUT Q(1) AND Q(2) MUST   PQFB 170
C                 CONTAIN INITIAL GUESSES FOR Q1 AND Q2 - ON RETURN Q(1)PQFB 180
C                 AND Q(2) CONTAIN THE REFINED COEFFICIENTS Q1 AND Q2 OFPQFB 190
C                 Q(X), WHILE Q(3) AND Q(4) CONTAIN THE COEFFICIENTS A  PQFB 200
C                 AND B OF A+B*X, WHICH IS THE REMAINDER OF THE QUOTIENTPQFB 210
C                 OF P(X) BY Q(X)                                       PQFB 220
C           LIM - INPUT VALUE SPECIFYING THE MAXIMUM NUMBER OF          PQFB 230
C                 ITERATIONS TO BE PERFORMED                            PQFB 240
C           IER - RESULTING ERROR PARAMETER (SEE REMARKS)               PQFB 250
C                 IER= 0 - NO ERROR                                     PQFB 260
C                 IER= 1 - NO CONVERGENCE WITHIN LIM ITERATIONS         PQFB 270
C                 IER=-1 - THE POLYNOMIAL P(X) IS CONSTANT OR UNDEFINED PQFB 280
C                          - OR OVERFLOW OCCURRED IN NORMALIZING P(X)   PQFB 290
C                 IER=-2 - THE POLYNOMIAL P(X) IS OF DEGREE 1           PQFB 300
C                 IER=-3 - NO FURTHER REFINEMENT OF THE APPROXIMATION TOPQFB 310
C                          A QUADRATIC FACTOR IS FEASIBLE, DUE TO EITHERPQFB 320
C                          DIVISION BY 0, OVERFLOW OR AN INITIAL GUESS  PQFB 330
C                          THAT IS NOT SUFFICIENTLY CLOSE TO A FACTOR OFPQFB 340
C                          P(X)                                         PQFB 350
C                                                                       PQFB 360
C        REMARKS                                                        PQFB 370
C           (1)  IF IER=-1 THERE IS NO COMPUTATION OTHER THAN THE       PQFB 380
C                POSSIBLE NORMALIZATION OF C.                           PQFB 390
C           (2)  IF IER=-2 THERE IS NO COMPUTATION OTHER THAN THE       PQFB 400
C                NORMALIZATION OF C.                                    PQFB 410
C           (3)  IF IER =-3  IT IS SUGGESTED THAT A NEW INITIAL GUESS BEPQFB 420
C                MADE FOR A QUADRATIC FACTOR.  Q, HOWEVER, WILL CONTAIN PQFB 430
C                THE VALUES ASSOCIATED WITH THE ITERATION THAT YIELDED  PQFB 440
C                THE SMALLEST NORM OF THE MODIFIED LINEAR REMAINDER.    PQFB 450
C           (4)  IF IER=1, THEN, ALTHOUGH THE NUMBER OF ITERATIONS LIM  PQFB 460
C                WAS TOO SMALL TO INDICATE CONVERGENCE, NO OTHER PROB-  PQFB 470
C                LEMS HAVE BEEN DETECTED, AND Q WILL CONTAIN THE VALUES PQFB 480
C                ASSOCIATED WITH THE ITERATION THAT YIELDED THE SMALLESTPQFB 490
C                NORM OF THE MODIFIED LINEAR REMAINDER.                 PQFB 500
C           (5)  FOR COMPLETE DETAIL SEE THE DOCUMENTATION FOR          PQFB 510
C                SUBROUTINES PQFB AND DPQFB.                            PQFB 520
C                                                                       PQFB 530
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PQFB 540
C           NONE                                                        PQFB 550
C                                                                       PQFB 560
C        METHOD                                                         PQFB 570
C           COMPUTATION IS BASED ON BAIRSTOW'S ITERATIVE METHOD.  (SEE  PQFB 580
C           WILKINSON, J.H., THE EVALUATION OF THE ZEROS OF ILL-CON-    PQFB 590
C           DITIONED POLYNOMIALS (PART ONE AND TWO), NUMERISCHE MATHE-  PQFB 600
C           MATIK, VOL.1 (1959), PP. 150-180, OR HILDEBRAND, F.B.,      PQFB 610
C           INTRODUCTION TO NUMERICAL ANALYSIS, MC GRAW-HILL, NEW YORK/ PQFB 620
C           TORONTO/LONDON, 1956, PP. 472-476.)                         PQFB 630
C                                                                       PQFB 640
C     ..................................................................PQFB 650
C                                                                       PQFB 660
      SUBROUTINE PQFB(C,IC,Q,LIM,IER)                                   PQFB 670
C                                                                       PQFB 680
C                                                                       PQFB 690
      DIMENSION C(1),Q(1)                                               PQFB 700
C                                                                       PQFB 710
C        TEST ON LEADING ZERO COEFFICIENTS                              PQFB 720
      IER=0                                                             PQFB 730
      J=IC+1                                                            PQFB 740
    1 J=J-1                                                             PQFB 750
      IF(J-1)40,40,2                                                    PQFB 760
    2 IF(C(J))3,1,3                                                     PQFB 770
C                                                                       PQFB 780
C        NORMALIZATION OF REMAINING COEFFICIENTS                        PQFB 790
    3 A=C(J)                                                            PQFB 800
      IF(A-1.)4,6,4                                                     PQFB 810
    4 DO 5 I=1,J                                                        PQFB 820
      C(I)=C(I)/A                                                       PQFB 830
      CALL OVERFL(N)                                                    PQFB 840
      IF(N-2)40,5,5                                                     PQFB 850
    5 CONTINUE                                                          PQFB 860
C                                                                       PQFB 870
C        TEST ON NECESSITY OF BAIRSTOW ITERATION                        PQFB 880
    6 IF(J-3)41,38,7                                                    PQFB 890
C                                                                       PQFB 900
C        PREPARE BAIRSTOW ITERATION                                     PQFB 910
    7 EPS=1.E-6                                                         PQFB 920
      EPS1=1.E-3                                                        PQFB 930
      L=0                                                               PQFB 940
      LL=0                                                              PQFB 950
      Q1=Q(1)                                                           PQFB 960
      Q2=Q(2)                                                           PQFB 970
      QQ1=0.                                                            PQFB 980
      QQ2=0.                                                            PQFB 990
      AA=C(1)                                                           PQFB1000
      BB=C(2)                                                           PQFB1010
      CB=ABS(AA)                                                        PQFB1020
      CA=ABS(BB)                                                        PQFB1030
      IF(CB-CA)8,9,10                                                   PQFB1040
    8 CC=CB+CB                                                          PQFB1050
      CB=CB/CA                                                          PQFB1060
      CA=1.                                                             PQFB1070
      GO TO 11                                                          PQFB1080
    9 CC=CA+CA                                                          PQFB1090
      CA=1.                                                             PQFB1100
      CB=1.                                                             PQFB1110
      GO TO 11                                                          PQFB1120
   10 CC=CA+CA                                                          PQFB1130
      CA=CA/CB                                                          PQFB1140
      CB=1.                                                             PQFB1150
   11 CD=CC*.1                                                          PQFB1160
C                                                                       PQFB1170
C        START BAIRSTOW ITERATION                                       PQFB1180
C        PREPARE NESTED MULTIPLICATION                                  PQFB1190
   12 A=0.                                                              PQFB1200
      B=A                                                               PQFB1210
      A1=A                                                              PQFB1220
      B1=A                                                              PQFB1230
      I=J                                                               PQFB1240
      QQQ1=Q1                                                           PQFB1250
      QQQ2=Q2                                                           PQFB1260
      DQ1=HH                                                            PQFB1270
      DQ2=H                                                             PQFB1280
C                                                                       PQFB1290
C        START NESTED MULTIPLICATION                                    PQFB1300
   13 H=-Q1*B-Q2*A+C(I)                                                 PQFB1310
      CALL OVERFL(N)                                                    PQFB1320
      IF(N-2)42,14,14                                                   PQFB1330
   14 B=A                                                               PQFB1340
      A=H                                                               PQFB1350
      I=I-1                                                             PQFB1360
      IF(I-1)18,15,16                                                   PQFB1370
   15 H=0.                                                              PQFB1380
   16 H=-Q1*B1-Q2*A1+H                                                  PQFB1390
      CALL OVERFL(N)                                                    PQFB1400
      IF(N-2)42,17,17                                                   PQFB1410
   17 C1=B1                                                             PQFB1420
      B1=A1                                                             PQFB1430
      A1=H                                                              PQFB1440
      GO TO 13                                                          PQFB1450
C        END OF NESTED MULTIPLICATION                                   PQFB1460
C                                                                       PQFB1470
C        TEST ON SATISFACTORY ACCURACY                                  PQFB1480
   18 H=CA*ABS(A)+CB*ABS(B)                                             PQFB1490
      IF(LL)19,19,39                                                    PQFB1500
   19 L=L+1                                                             PQFB1510
      IF(ABS(A)-EPS*ABS(C(1)))20,20,21                                  PQFB1520
   20 IF(ABS(B)-EPS*ABS(C(2)))39,39,21                                  PQFB1530
C                                                                       PQFB1540
C        TEST ON LINEAR REMAINDER OF MINIMUM NORM                       PQFB1550
   21 IF(H-CC)22,22,23                                                  PQFB1560
   22 AA=A                                                              PQFB1570
      BB=B                                                              PQFB1580
      CC=H                                                              PQFB1590
      QQ1=Q1                                                            PQFB1600
      QQ2=Q2                                                            PQFB1610
C                                                                       PQFB1620
C        TEST ON LAST ITERATION STEP                                    PQFB1630
   23 IF(L-LIM)28,28,24                                                 PQFB1640
C                                                                       PQFB1650
C        TEST ON RESTART OF BAIRSTOW ITERATION WITH ZERO INITIAL GUESS  PQFB1660
   24 IF(H-CD)43,43,25                                                  PQFB1670
   25 IF(Q(1))27,26,27                                                  PQFB1680
   26 IF(Q(2))27,42,27                                                  PQFB1690
   27 Q(1)=0.                                                           PQFB1700
      Q(2)=0.                                                           PQFB1710
      GO TO 7                                                           PQFB1720
C                                                                       PQFB1730
C        PERFORM ITERATION STEP                                         PQFB1740
   28 HH=AMAX1(ABS(A1),ABS(B1),ABS(C1))                                 PQFB1750
      IF(HH)42,42,29                                                    PQFB1760
   29 A1=A1/HH                                                          PQFB1770
      B1=B1/HH                                                          PQFB1780
      C1=C1/HH                                                          PQFB1790
      H=A1*C1-B1*B1                                                     PQFB1800
      IF(H)30,42,30                                                     PQFB1810
   30 A=A/HH                                                            PQFB1820
      B=B/HH                                                            PQFB1830
      HH=(B*A1-A*B1)/H                                                  PQFB1840
      H=(A*C1-B*B1)/H                                                   PQFB1850
      Q1=Q1+HH                                                          PQFB1860
      Q2=Q2+H                                                           PQFB1870
C        END OF ITERATION STEP                                          PQFB1880
C                                                                       PQFB1890
C        TEST ON SATISFACTORY RELATIVE ERROR OF ITERATED VALUES         PQFB1900
      IF(ABS(HH)-EPS*ABS(Q1))31,31,33                                   PQFB1910
   31 IF(ABS(H)-EPS*ABS(Q2))32,32,33                                    PQFB1920
   32 LL=1                                                              PQFB1930
      GO TO 12                                                          PQFB1940
C                                                                       PQFB1950
C        TEST ON DECREASING RELATIVE ERRORS                             PQFB1960
   33 IF(L-1)12,12,34                                                   PQFB1970
   34 IF(ABS(HH)-EPS1*ABS(Q1))35,35,12                                  PQFB1980
   35 IF(ABS(H)-EPS1*ABS(Q2))36,36,12                                   PQFB1990
   36 IF(ABS(QQQ1*HH)-ABS(Q1*DQ1))37,44,44                              PQFB2000
   37 IF(ABS(QQQ2*H)-ABS(Q2*DQ2))12,44,44                               PQFB2010
C        END OF BAIRSTOW ITERATION                                      PQFB2020
C                                                                       PQFB2030
C        EXIT IN CASE OF QUADRATIC POLYNOMIAL                           PQFB2040
   38 Q(1)=C(1)                                                         PQFB2050
      Q(2)=C(2)                                                         PQFB2060
      Q(3)=0.                                                           PQFB2070
      Q(4)=0.                                                           PQFB2080
      RETURN                                                            PQFB2090
C                                                                       PQFB2100
C        EXIT IN CASE OF SUFFICIENT ACCURACY                            PQFB2110
   39 Q(1)=Q1                                                           PQFB2120
      Q(2)=Q2                                                           PQFB2130
      Q(3)=A                                                            PQFB2140
      Q(4)=B                                                            PQFB2150
      RETURN                                                            PQFB2160
C                                                                       PQFB2170
C        ERROR EXIT IN CASE OF ZERO OR CONSTANT POLYNOMIAL              PQFB2180
   40 IER=-1                                                            PQFB2190
      RETURN                                                            PQFB2200
C                                                                       PQFB2210
C        ERROR EXIT IN CASE OF LINEAR POLYNOMIAL                        PQFB2220
   41 IER=-2                                                            PQFB2230
      RETURN                                                            PQFB2240
C                                                                       PQFB2250
C        ERROR EXIT IN CASE OF NONREFINED QUADRATIC FACTOR              PQFB2260
   42 IER=-3                                                            PQFB2270
      GO TO 44                                                          PQFB2280
C                                                                       PQFB2290
C        ERROR EXIT IN CASE OF UNSATISFACTORY ACCURACY                  PQFB2300
   43 IER=1                                                             PQFB2310
   44 Q(1)=QQ1                                                          PQFB2320
      Q(2)=QQ2                                                          PQFB2330
      Q(3)=AA                                                           PQFB2340
      Q(4)=BB                                                           PQFB2350
      RETURN                                                            PQFB2360
      END                                                               PQFB2370
