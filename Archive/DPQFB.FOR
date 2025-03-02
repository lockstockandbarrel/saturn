C                                                                       DPQF  10
C     ..................................................................DPQF  20
C                                                                       DPQF  30
C        SUBROUTINE DPQFB                                               DPQF  40
C                                                                       DPQF  50
C        PURPOSE                                                        DPQF  60
C           TO FIND AN APPROXIMATION Q(X)=Q1+Q2*X+X*X TO A QUADRATIC    DPQF  70
C           FACTOR OF A GIVEN POLYNOMIAL P(X) WITH REAL COEFFICIENTS.   DPQF  80
C                                                                       DPQF  90
C        USAGE                                                          DPQF 100
C           CALL DPQFB(C,IC,Q,LIM,IER)                                  DPQF 110
C                                                                       DPQF 120
C        DESCRIPTION OF PARAMETERS                                      DPQF 130
C           C   - DOUBLE PRECISION INPUT VECTOR CONTAINING THE          DPQF 140
C                 COEFFICIENTS OF P(X) - C(1) IS THE CONSTANT TERM      DPQF 150
C                 (DIMENSION IC)                                        DPQF 160
C           IC  - DIMENSION OF C                                        DPQF 170
C           Q   - DOUBLE PRECISION VECTOR OF DIMENSION 4 - ON INPUT Q(1)DPQF 180
C                 AND Q(2) CONTAIN INITIAL GUESSES FOR Q1 AND Q2 - ON   DPQF 190
C                 RETURN Q(1) AND Q(2) CONTAIN THE REFINED COEFFICIENTS DPQF 200
C                 Q1 AND Q2 OF Q(X), WHILE Q(3) AND Q(4) CONTAIN THE    DPQF 210
C                 COEFFICIENTS A AND B OF A+B*X, WHICH IS THE REMAINDER DPQF 220
C                 OF THE QUOTIENT OF P(X) BY Q(X)                       DPQF 230
C           LIM - INPUT VALUE SPECIFYING THE MAXIMUM NUMBER OF          DPQF 240
C                 ITERATIONS TO BE PERFORMED                            DPQF 250
C           IER - RESULTING ERROR PARAMETER (SEE REMARKS)               DPQF 260
C                 IER= 0 - NO ERROR                                     DPQF 270
C                 IER= 1 - NO CONVERGENCE WITHIN LIM ITERATIONS         DPQF 280
C                 IER=-1 - THE POLYNOMIAL P(X) IS CONSTANT OR UNDEFINED DPQF 290
C                          - OR OVERFLOW OCCURRED IN NORMALIZING P(X)   DPQF 300
C                 IER=-2 - THE POLYNOMIAL P(X) IS OF DEGREE 1           DPQF 310
C                 IER=-3 - NO FURTHER REFINEMENT OF THE APPROXIMATION TODPQF 320
C                          A QUADRATIC FACTOR IS FEASIBLE, DUE TO EITHERDPQF 330
C                          DIVISION BY 0, OVERFLOW OR AN INITIAL GUESS  DPQF 340
C                          THAT IS NOT SUFFICIENTLY CLOSE TO A FACTOR OFDPQF 350
C                          P(X)                                         DPQF 360
C                                                                       DPQF 370
C        REMARKS                                                        DPQF 380
C           (1)  IF IER=-1 THERE IS NO COMPUTATION OTHER THAN THE       DPQF 390
C                POSSIBLE NORMALIZATION OF C.                           DPQF 400
C           (2)  IF IER=-2 THERE IS NO COMPUTATION OTHER THAN THE       DPQF 410
C                NORMALIZATION OF C.                                    DPQF 420
C           (3)  IF IER =-3  IT IS SUGGESTED THAT A NEW INITIAL GUESS BEDPQF 430
C                MADE FOR A QUADRATIC FACTOR.  Q, HOWEVER, WILL CONTAIN DPQF 440
C                THE VALUES ASSOCIATED WITH THE ITERATION THAT YIELDED  DPQF 450
C                THE SMALLEST NORM OF THE MODIFIED LINEAR REMAINDER.    DPQF 460
C           (4)  IF IER=1, THEN, ALTHOUGH THE NUMBER OF ITERATIONS LIM  DPQF 470
C                WAS TOO SMALL TO INDICATE CONVERGENCE, NO OTHER PROB-  DPQF 480
C                LEMS HAVE BEEN DETECTED, AND Q WILL CONTAIN THE VALUES DPQF 490
C                ASSOCIATED WITH THE ITERATION THAT YIELDED THE SMALLESTDPQF 500
C                NORM OF THE MODIFIED LINEAR REMAINDER.                 DPQF 510
C           (5)  FOR COMPLETE DETAIL SEE THE DOCUMENTATION FOR          DPQF 520
C                SUBROUTINES PQFB AND DPQFB.                            DPQF 530
C                                                                       DPQF 540
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DPQF 550
C           NONE                                                        DPQF 560
C                                                                       DPQF 570
C        METHOD                                                         DPQF 580
C           COMPUTATION IS BASED ON BAIRSTOW'S ITERATIVE METHOD.  (SEE  DPQF 590
C           WILKINSON, J.H., THE EVALUATION OF THE ZEROS OF ILL-CON-    DPQF 600
C           DITIONED POLYNOMIALS (PART ONE AND TWO), NUMERISCHE MATHE-  DPQF 610
C           MATIK, VOL.1 (1959), PP. 150-180, OR HILDEBRAND, F.B.,      DPQF 620
C           INTRODUCTION TO NUMERICAL ANALYSIS, MC GRAW-HILL, NEW YORK/ DPQF 630
C           TORONTO/LONDON, 1956, PP. 472-476.)                         DPQF 640
C                                                                       DPQF 650
C     ..................................................................DPQF 660
C                                                                       DPQF 670
      SUBROUTINE DPQFB(C,IC,Q,LIM,IER)                                  DPQF 680
C                                                                       DPQF 690
C                                                                       DPQF 700
      DIMENSION C(1),Q(1)                                               DPQF 710
      DOUBLE PRECISION A,B,AA,BB,CA,CB,CC,CD,A1,B1,C1,H,HH,Q1,Q2,QQ1,   DPQF 720
     1                 QQ2,QQQ1,QQQ2,DQ1,DQ2,EPS,EPS1,C,Q               DPQF 730
C                                                                       DPQF 740
C        TEST ON LEADING ZERO COEFFICIENTS                              DPQF 750
      IER=0                                                             DPQF 760
      J=IC+1                                                            DPQF 770
    1 J=J-1                                                             DPQF 780
      IF(J-1)40,40,2                                                    DPQF 790
    2 IF(C(J))3,1,3                                                     DPQF 800
C                                                                       DPQF 810
C        NORMALIZATION OF REMAINING COEFFICIENTS                        DPQF 820
    3 A=C(J)                                                            DPQF 830
      IF(A-1.D0)4,6,4                                                   DPQF 840
    4 DO 5 I=1,J                                                        DPQF 850
      C(I)=C(I)/A                                                       DPQF 860
      CALL OVERFL(N)                                                    DPQF 870
      IF(N-2)40,5,5                                                     DPQF 880
    5 CONTINUE                                                          DPQF 890
C                                                                       DPQF 900
C        TEST ON NECESSITY OF BAIRSTOW ITERATION                        DPQF 910
    6 IF(J-3)41,38,7                                                    DPQF 920
C                                                                       DPQF 930
C        PREPARE BAIRSTOW ITERATION                                     DPQF 940
    7 EPS=1.D-14                                                        DPQF 950
      EPS1=1.D-6                                                        DPQF 960
      L=0                                                               DPQF 970
      LL=0                                                              DPQF 980
      Q1=Q(1)                                                           DPQF 990
      Q2=Q(2)                                                           DPQF1000
      QQ1=0.D0                                                          DPQF1010
      QQ2=0.D0                                                          DPQF1020
      AA=C(1)                                                           DPQF1030
      BB=C(2)                                                           DPQF1040
      CB=DABS(AA)                                                       DPQF1050
      CA=DABS(BB)                                                       DPQF1060
      IF(CB-CA)8,9,10                                                   DPQF1070
    8 CC=CB+CB                                                          DPQF1080
      CB=CB/CA                                                          DPQF1090
      CA=1.D0                                                           DPQF1100
      GO TO 11                                                          DPQF1110
    9 CC=CA+CA                                                          DPQF1120
      CA=1.D0                                                           DPQF1130
      CB=1.D0                                                           DPQF1140
      GO TO 11                                                          DPQF1150
   10 CC=CA+CA                                                          DPQF1160
      CA=CA/CB                                                          DPQF1170
      CB=1.D0                                                           DPQF1180
   11 CD=CC*.1D0                                                        DPQF1190
C                                                                       DPQF1200
C        START BAIRSTOW ITERATION                                       DPQF1210
C        PREPARE NESTED MULTIPLICATION                                  DPQF1220
   12 A=0.D0                                                            DPQF1230
      B=A                                                               DPQF1240
      A1=A                                                              DPQF1250
      B1=A                                                              DPQF1260
      I=J                                                               DPQF1270
      QQQ1=Q1                                                           DPQF1280
      QQQ2=Q2                                                           DPQF1290
      DQ1=HH                                                            DPQF1300
      DQ2=H                                                             DPQF1310
C                                                                       DPQF1320
C        START NESTED MULTIPLICATION                                    DPQF1330
   13 H=-Q1*B-Q2*A+C(I)                                                 DPQF1340
      CALL OVERFL(N)                                                    DPQF1350
      IF(N-2)42,14,14                                                   DPQF1360
   14 B=A                                                               DPQF1370
      A=H                                                               DPQF1380
      I=I-1                                                             DPQF1390
      IF(I-1)18,15,16                                                   DPQF1400
   15 H=0.D0                                                            DPQF1410
   16 H=-Q1*B1-Q2*A1+H                                                  DPQF1420
      CALL OVERFL(N)                                                    DPQF1430
      IF(N-2)42,17,17                                                   DPQF1440
   17 C1=B1                                                             DPQF1450
      B1=A1                                                             DPQF1460
      A1=H                                                              DPQF1470
      GO TO 13                                                          DPQF1480
C        END OF NESTED MULTIPLICATION                                   DPQF1490
C                                                                       DPQF1500
C        TEST ON SATISFACTORY ACCURACY                                  DPQF1510
   18 H=CA*DABS(A)+CB*DABS(B)                                           DPQF1520
      IF(LL)19,19,39                                                    DPQF1530
   19 L=L+1                                                             DPQF1540
      IF(DABS(A)-EPS*DABS(C(1)))20,20,21                                DPQF1550
   20 IF(DABS(B)-EPS*DABS(C(2)))39,39,21                                DPQF1560
C                                                                       DPQF1570
C        TEST ON LINEAR REMAINDER OF MINIMUM NORM                       DPQF1580
   21 IF(H-CC)22,22,23                                                  DPQF1590
   22 AA=A                                                              DPQF1600
      BB=B                                                              DPQF1610
      CC=H                                                              DPQF1620
      QQ1=Q1                                                            DPQF1630
      QQ2=Q2                                                            DPQF1640
C                                                                       DPQF1650
C        TEST ON LAST ITERATION STEP                                    DPQF1660
   23 IF(L-LIM)28,28,24                                                 DPQF1670
C                                                                       DPQF1680
C        TEST ON RESTART OF BAIRSTOW ITERATION WITH ZERO INITIAL GUESS  DPQF1690
   24 IF(H-CD)43,43,25                                                  DPQF1700
   25 IF(Q(1))27,26,27                                                  DPQF1710
   26 IF(Q(2))27,42,27                                                  DPQF1720
   27 Q(1)=0.D0                                                         DPQF1730
      Q(2)=0.D0                                                         DPQF1740
      GO TO 7                                                           DPQF1750
C                                                                       DPQF1760
C        PERFORM ITERATION STEP                                         DPQF1770
   28 HH=DMAX1(DABS(A1),DABS(B1),DABS(C1))                              DPQF1780
      IF(HH)42,42,29                                                    DPQF1790
   29 A1=A1/HH                                                          DPQF1800
      B1=B1/HH                                                          DPQF1810
      C1=C1/HH                                                          DPQF1820
      H=A1*C1-B1*B1                                                     DPQF1830
      IF(H)30,42,30                                                     DPQF1840
   30 A=A/HH                                                            DPQF1850
      B=B/HH                                                            DPQF1860
      HH=(B*A1-A*B1)/H                                                  DPQF1870
      H=(A*C1-B*B1)/H                                                   DPQF1880
      Q1=Q1+HH                                                          DPQF1890
      Q2=Q2+H                                                           DPQF1900
C        END OF ITERATION STEP                                          DPQF1910
C                                                                       DPQF1920
C        TEST ON SATISFACTORY RELATIVE ERROR OF ITERATED VALUES         DPQF1930
      IF(DABS(HH)-EPS*DABS(Q1))31,31,33                                 DPQF1940
   31 IF(DABS(H)-EPS*DABS(Q2))32,32,33                                  DPQF1950
   32 LL=1                                                              DPQF1960
      GO TO 12                                                          DPQF1970
C                                                                       DPQF1980
C        TEST ON DECREASING RELATIVE ERRORS                             DPQF1990
   33 IF(L-1)12,12,34                                                   DPQF2000
   34 IF(DABS(HH)-EPS1*DABS(Q1))35,35,12                                DPQF2010
   35 IF(DABS(H)-EPS1*DABS(Q2))36,36,12                                 DPQF2020
   36 IF(DABS(QQQ1*HH)-DABS(Q1*DQ1))37,44,44                            DPQF2030
   37 IF(DABS(QQQ2*H)-DABS(Q2*DQ2))12,44,44                             DPQF2040
C        END OF BAIRSTOW ITERATION                                      DPQF2050
C                                                                       DPQF2060
C        EXIT IN CASE OF QUADRATIC POLYNOMIAL                           DPQF2070
   38 Q(1)=C(1)                                                         DPQF2080
      Q(2)=C(2)                                                         DPQF2090
      Q(3)=0.D0                                                         DPQF2100
      Q(4)=0.D0                                                         DPQF2110
      RETURN                                                            DPQF2120
C                                                                       DPQF2130
C        EXIT IN CASE OF SUFFICIENT ACCURACY                            DPQF2140
   39 Q(1)=Q1                                                           DPQF2150
      Q(2)=Q2                                                           DPQF2160
      Q(3)=A                                                            DPQF2170
      Q(4)=B                                                            DPQF2180
      RETURN                                                            DPQF2190
C                                                                       DPQF2200
C        ERROR EXIT IN CASE OF ZERO OR CONSTANT POLYNOMIAL              DPQF2210
   40 IER=-1                                                            DPQF2220
      RETURN                                                            DPQF2230
C                                                                       DPQF2240
C        ERROR EXIT IN CASE OF LINEAR POLYNOMIAL                        DPQF2250
   41 IER=-2                                                            DPQF2260
      RETURN                                                            DPQF2270
C                                                                       DPQF2280
C        ERROR EXIT IN CASE OF NONREFINED QUADRATIC FACTOR              DPQF2290
   42 IER=-3                                                            DPQF2300
      GO TO 44                                                          DPQF2310
C                                                                       DPQF2320
C        ERROR EXIT IN CASE OF UNSATISFACTORY ACCURACY                  DPQF2330
   43 IER=1                                                             DPQF2340
   44 Q(1)=QQ1                                                          DPQF2350
      Q(2)=QQ2                                                          DPQF2360
      Q(3)=AA                                                           DPQF2370
      Q(4)=BB                                                           DPQF2380
      RETURN                                                            DPQF2390
      END                                                               DPQF2400
