C                                                                       DPRB  10
C     ..................................................................DPRB  20
C                                                                       DPRB  30
C        SUBROUTINE DPRBM                                               DPRB  40
C                                                                       DPRB  50
C        PURPOSE                                                        DPRB  60
C           TO CALCULATE ALL REAL AND COMPLEX ROOTS OF A GIVEN          DPRB  70
C           POLYNOMIAL WITH REAL COEFFICIENTS.                          DPRB  80
C                                                                       DPRB  90
C        USAGE                                                          DPRB 100
C           CALL DPRBM (C,IC,RR,RC,POL,IR,IER)                          DPRB 110
C                                                                       DPRB 120
C        DESCRIPTION OF PARAMETERS                                      DPRB 130
C           C      - DOUBLE PRECISION INPUT VECTOR CONTAINING THE       DPRB 140
C                    COEFFICIENTS OF THE GIVEN POLYNOMIAL. COEFFICIENTS DPRB 150
C                    ARE ORDERED FROM LOW TO HIGH. ON RETURN COEFFI-    DPRB 160
C                    CIENTS ARE DIVIDED BY THE LAST NONZERO TERM.       DPRB 170
C           IC     - DIMENSION OF VECTORS C, RR, RC, AND POL.           DPRB 180
C           RR     - RESULTANT DOUBLE PRECISION VECTOR OF REAL PARTS    DPRB 190
C                    OF THE ROOTS.                                      DPRB 200
C           RC     - RESULTANT DOUBLE PRECISION VECTOR OF COMPLEX PARTS DPRB 210
C                    OF THE ROOTS.                                      DPRB 220
C           POL    - RESULTANT DOUBLE PRECISION VECTOR OF COEFFICIENTS  DPRB 230
C                    OF THE POLYNOMIAL WITH CALCULATED ROOTS.           DPRB 240
C                    COEFFICIENTS ARE ORDERED FROM LOW TO HIGH (SEE     DPRB 250
C                    REMARK 4).                                         DPRB 260
C           IR     - OUTPUT VALUE SPECIFYING THE NUMBER OF CALCULATED   DPRB 270
C                    ROOTS. NORMALLY IR IS EQUAL TO IC-1.               DPRB 280
C           IER    - RESULTANT ERROR PARAMETER CODED AS FOLLOWS         DPRB 290
C                     IER=0  - NO ERROR,                                DPRB 300
C                     IER=1  - SUBROUTINE DPQFB RECORDS POOR CONVERGENCEDPRB 310
C                              AT SOME QUADRATIC FACTORIZATION WITHIN   DPRB 320
C                              100 ITERATION STEPS,                     DPRB 330
C                     IER=2  - POLYNOMIAL IS DEGENERATE, I.E. ZERO OR   DPRB 340
C                              CONSTANT,                                DPRB 350
C                              OR OVERFLOW IN NORMALIZATION OF GIVEN    DPRB 360
C                              POLYNOMIAL,                              DPRB 370
C                     IER=3  - THE SUBROUTINE IS BYPASSED DUE TO        DPRB 380
C                              SUCCESSIVE ZERO DIVISORS OR OVERFLOWS    DPRB 390
C                              IN QUADRATIC FACTORIZATION OR DUE TO     DPRB 400
C                              COMPLETELY UNSATISFACTORY ACCURACY,      DPRB 410
C                     IER=-1 - CALCULATED COEFFICIENT VECTOR HAS LESS   DPRB 420
C                              THAN SIX CORRECT SIGNIFICANT DIGITS.     DPRB 430
C                              THIS REVEALS POOR ACCURACY OF CALCULATED DPRB 440
C                              ROOTS.                                   DPRB 450
C                                                                       DPRB 460
C        REMARKS                                                        DPRB 470
C           (1) REAL PARTS OF THE ROOTS ARE STORED IN RR(1) UP TO RR(IR)DPRB 480
C               AND CORRESPONDING COMPLEX PARTS IN RC(1) UP TO RC(IR).  DPRB 490
C           (2) ERROR MESSAGE IER=1 INDICATES POOR CONVERGENCE WITHIN   DPRB 500
C               100 ITERATION STEPS AT SOME QUADRATIC FACTORIZATION     DPRB 510
C               PERFORMED BY SUBROUTINE DPQFB.                          DPRB 520
C           (3) NO ACTION BESIDES ERROR MESSAGE IER=2 IN CASE OF A ZERO DPRB 530
C               OR CONSTANT POLYNOMIAL. THE SAME ERROR MESSAGE IS GIVEN DPRB 540
C               IN CASE OF AN OVERFLOW IN NORMALIZATION OF GIVEN        DPRB 550
C               POLYNOMIAL.                                             DPRB 560
C           (4) ERROR MESSAGE IER=3 INDICATES SUCCESSIVE ZERO DIVISORS  DPRB 570
C               OR OVERFLOWS OR COMPLETELY UNSATISFACTORY ACCURACY AT   DPRB 580
C               ANY QUADRATIC FACTORIZATION PERFORMED BY                DPRB 590
C               SUBROUTINE DPQFB. IN THIS CASE CALCULATION IS BYPASSED. DPRB 600
C               IR RECORDS THE NUMBER OF CALCULATED ROOTS.              DPRB 610
C               POL(1),...,POL(J-IR) ARE THE COEFFICIENTS OF THE        DPRB 620
C               REMAINING POLYNOMIAL, WHERE J IS THE ACTUAL NUMBER OF   DPRB 630
C               COEFFICIENTS IN VECTOR C (NORMALLY J=IC).               DPRB 640
C           (5) IF CALCULATED COEFFICIENT VECTOR HAS LESS THAN SIX      DPRB 650
C               CORRECT SIGNIFICANT DIGITS THOUGH ALL QUADRATIC         DPRB 660
C               FACTORIZATIONS SHOWED SATISFACTORY ACCURACY, THE ERROR  DPRB 670
C               MESSAGE IER=-1 IS GIVEN.                                DPRB 680
C           (6) THE FINAL COMPARISON BETWEEN GIVEN AND CALCULATED       DPRB 690
C               COEFFICIENT VECTOR IS PERFORMED ONLY IF ALL ROOTS HAVE  DPRB 700
C               BEEN CALCULATED. IN THIS CASE THE NUMBER OF ROOTS IR IS DPRB 710
C               EQUAL TO THE ACTUAL DEGREE OF THE POLYNOMIAL (NORMALLY  DPRB 720
C               IR=IC-1). THE MAXIMAL RELATIVE ERROR OF THE COEFFICIENT DPRB 730
C               VECTOR IS RECORDED IN RR(IR+1).                         DPRB 740
C                                                                       DPRB 750
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DPRB 760
C           SUBROUTINE DPQFB    QUADRATIC FACTORIZATION OF A POLYNOMIAL DPRB 770
C                               BY BAIRSTOW ITERATION.                  DPRB 780
C                                                                       DPRB 790
C        METHOD                                                         DPRB 800
C           THE ROOTS OF THE POLYNOMIAL ARE CALCULATED BY MEANS OF      DPRB 810
C           SUCCESSIVE QUADRATIC FACTORIZATION PERFORMED BY BAIRSTOW    DPRB 820
C           ITERATION. X**2 IS USED AS INITIAL GUESS FOR THE FIRST      DPRB 830
C           QUADRATIC FACTOR, AND FURTHER EACH CALCULATED QUADRATIC     DPRB 840
C           FACTOR IS USED AS INITIAL GUESS FOR THE NEXT ONE. AFTER     DPRB 850
C           COMPUTATION OF ALL ROOTS THE COEFFICIENT VECTOR IS          DPRB 860
C           CALCULATED AND COMPARED WITH THE GIVEN ONE.                 DPRB 870
C           FOR REFERENCE, SEE J. H. WILKINSON, THE EVALUATION OF THE   DPRB 880
C           ZEROS OF ILL-CONDITIONED POLYNOMIALS (PART ONE AND TWO),    DPRB 890
C           NUMERISCHE MATHEMATIK, VOL.1 (1959), PP.150-180.            DPRB 900
C                                                                       DPRB 910
C     ..................................................................DPRB 920
C                                                                       DPRB 930
      SUBROUTINE DPRBM(C,IC,RR,RC,POL,IR,IER)                           DPRB 940
C                                                                       DPRB 950
C                                                                       DPRB 960
      DIMENSION C(1),RR(1),RC(1),POL(1),Q(4)                            DPRB 970
      DOUBLE PRECISION C,RR,RC,POL,Q,EPS,A,B,H,Q1,Q2                    DPRB 980
C                                                                       DPRB 990
C        TEST ON LEADING ZERO COEFFICIENTS                              DPRB1000
      EPS=1.D-6                                                         DPRB1010
      LIM=100                                                           DPRB1020
      IR=IC+1                                                           DPRB1030
    1 IR=IR-1                                                           DPRB1040
      IF(IR-1)42,42,2                                                   DPRB1050
    2 IF(C(IR))3,1,3                                                    DPRB1060
C                                                                       DPRB1070
C        WORK UP ZERO ROOTS AND NORMALIZE REMAINING POLYNOMIAL          DPRB1080
    3 IER=0                                                             DPRB1090
      J=IR                                                              DPRB1100
      L=0                                                               DPRB1110
      A=C(IR)                                                           DPRB1120
      DO 8 I=1,IR                                                       DPRB1130
      IF(L)4,4,7                                                        DPRB1140
    4 IF(C(I))6,5,6                                                     DPRB1150
    5 RR(I)=0.D0                                                        DPRB1160
      RC(I)=0.D0                                                        DPRB1170
      POL(J)=0.D0                                                       DPRB1180
      J=J-1                                                             DPRB1190
      GO TO 8                                                           DPRB1200
    6 L=1                                                               DPRB1210
      IST=I                                                             DPRB1220
      J=0                                                               DPRB1230
    7 J=J+1                                                             DPRB1240
      C(I)=C(I)/A                                                       DPRB1250
      POL(J)=C(I)                                                       DPRB1260
      CALL OVERFL(N)                                                    DPRB1270
      IF(N-2)42,8,8                                                     DPRB1280
    8 CONTINUE                                                          DPRB1290
C                                                                       DPRB1300
C        START BAIRSTOW ITERATION                                       DPRB1310
      Q1=0.D0                                                           DPRB1320
      Q2=0.D0                                                           DPRB1330
    9 IF(J-2)33,10,14                                                   DPRB1340
C                                                                       DPRB1350
C        DEGREE OF RESTPOLYNOMIAL IS EQUAL TO ONE                       DPRB1360
   10 A=POL(1)                                                          DPRB1370
      RR(IST)=-A                                                        DPRB1380
      RC(IST)=0.D0                                                      DPRB1390
      IR=IR-1                                                           DPRB1400
      Q2=0.D0                                                           DPRB1410
      IF(IR-1)13,13,11                                                  DPRB1420
   11 DO 12 I=2,IR                                                      DPRB1430
      Q1=Q2                                                             DPRB1440
      Q2=POL(I+1)                                                       DPRB1450
   12 POL(I)=A*Q2+Q1                                                    DPRB1460
   13 POL(IR+1)=A+Q2                                                    DPRB1470
      GO TO 34                                                          DPRB1480
C        THIS IS BRANCH TO COMPARISON OF COEFFICIENT VECTORS C AND POL  DPRB1490
C                                                                       DPRB1500
C        DEGREE OF RESTPOLYNOMIAL IS GREATER THAN ONE                   DPRB1510
   14 DO 22 L=1,10                                                      DPRB1520
      N=1                                                               DPRB1530
   15 Q(1)=Q1                                                           DPRB1540
      Q(2)=Q2                                                           DPRB1550
      CALL DPQFB(POL,J,Q,LIM,I)                                         DPRB1560
      IF(I)16,24,23                                                     DPRB1570
   16 IF(Q1)18,17,18                                                    DPRB1580
   17 IF(Q2)18,21,18                                                    DPRB1590
   18 GO TO (19,20,19,21),N                                             DPRB1600
   19 Q1=-Q1                                                            DPRB1610
      N=N+1                                                             DPRB1620
      GO TO 15                                                          DPRB1630
   20 Q2=-Q2                                                            DPRB1640
      N=N+1                                                             DPRB1650
      GO TO 15                                                          DPRB1660
   21 Q1=1.D0+Q1                                                        DPRB1670
   22 Q2=1.D0-Q2                                                        DPRB1680
C                                                                       DPRB1690
C        ERROR EXIT DUE TO UNSATISFACTORY RESULTS OF FACTORIZATION      DPRB1700
      IER=3                                                             DPRB1710
      IR=IR-J                                                           DPRB1720
      RETURN                                                            DPRB1730
C                                                                       DPRB1740
C        WORK UP RESULTS OF QUADRATIC FACTORIZATION                     DPRB1750
   23 IER=1                                                             DPRB1760
   24 Q1=Q(1)                                                           DPRB1770
      Q2=Q(2)                                                           DPRB1780
C                                                                       DPRB1790
C        PERFORM DIVISION OF FACTORIZED POLYNOMIAL BY QUADRATIC FACTOR  DPRB1800
      B=0.D0                                                            DPRB1810
      A=0.D0                                                            DPRB1820
      I=J                                                               DPRB1830
   25 H=-Q1*B-Q2*A+POL(I)                                               DPRB1840
      POL(I)=B                                                          DPRB1850
      B=A                                                               DPRB1860
      A=H                                                               DPRB1870
      I=I-1                                                             DPRB1880
      IF(I-2)26,26,25                                                   DPRB1890
   26 POL(2)=B                                                          DPRB1900
      POL(1)=A                                                          DPRB1910
C                                                                       DPRB1920
C        MULTIPLY POLYNOMIAL WITH CALCULATED ROOTS BY QUADRATIC FACTOR  DPRB1930
      L=IR-1                                                            DPRB1940
      IF(J-L)27,27,29                                                   DPRB1950
   27 DO 28 I=J,L                                                       DPRB1960
   28 POL(I-1)=POL(I-1)+POL(I)*Q2+POL(I+1)*Q1                           DPRB1970
   29 POL(L)=POL(L)+POL(L+1)*Q2+Q1                                      DPRB1980
      POL(IR)=POL(IR)+Q2                                                DPRB1990
C                                                                       DPRB2000
C        CALCULATE ROOT-PAIR FROM QUADRATIC FACTOR X*X+Q2*X+Q1          DPRB2010
      H=-.5D0*Q2                                                        DPRB2020
      A=H*H-Q1                                                          DPRB2030
      B=DSQRT(DABS(A))                                                  DPRB2040
      IF(A)30,30,31                                                     DPRB2050
   30 RR(IST)=H                                                         DPRB2060
      RC(IST)=B                                                         DPRB2070
      IST=IST+1                                                         DPRB2080
      RR(IST)=H                                                         DPRB2090
      RC(IST)=-B                                                        DPRB2100
      GO TO 32                                                          DPRB2110
   31 B=H+DSIGN(B,H)                                                    DPRB2120
      RR(IST)=Q1/B                                                      DPRB2130
      RC(IST)=0.D0                                                      DPRB2140
      IST=IST+1                                                         DPRB2150
      RR(IST)=B                                                         DPRB2160
      RC(IST)=0.D0                                                      DPRB2170
   32 IST=IST+1                                                         DPRB2180
      J=J-2                                                             DPRB2190
      GO TO 9                                                           DPRB2200
C                                                                       DPRB2210
C        SHIFT BACK ELEMENTS OF POL BY 1 AND COMPARE VECTORS POL AND C  DPRB2220
   33 IR=IR-1                                                           DPRB2230
   34 A=0.D0                                                            DPRB2240
      DO 38 I=1,IR                                                      DPRB2250
      Q1=C(I)                                                           DPRB2260
      Q2=POL(I+1)                                                       DPRB2270
      POL(I)=Q2                                                         DPRB2280
      IF(Q1)35,36,35                                                    DPRB2290
   35 Q2=(Q1-Q2)/Q1                                                     DPRB2300
   36 Q2=DABS(Q2)                                                       DPRB2310
      IF(Q2-A)38,38,37                                                  DPRB2320
   37 A=Q2                                                              DPRB2330
   38 CONTINUE                                                          DPRB2340
      I=IR+1                                                            DPRB2350
      POL(I)=1.D0                                                       DPRB2360
      RR(I)=A                                                           DPRB2370
      RC(I)=0.D0                                                        DPRB2380
      IF(IER)39,39,41                                                   DPRB2390
   39 IF(A-EPS)41,41,40                                                 DPRB2400
C                                                                       DPRB2410
C        WARNING DUE TO POOR ACCURACY OF CALCULATED COEFFICIENT VECTOR  DPRB2420
   40 IER=-1                                                            DPRB2430
   41 RETURN                                                            DPRB2440
C                                                                       DPRB2450
C        ERROR EXIT DUE TO DEGENERATE POLYNOMIAL OR OVERFLOW IN         DPRB2460
C        NORMALIZATION                                                  DPRB2470
   42 IER=2                                                             DPRB2480
      IR=0                                                              DPRB2490
      RETURN                                                            DPRB2500
      END                                                               DPRB2510
