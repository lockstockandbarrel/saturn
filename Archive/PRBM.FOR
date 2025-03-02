C                                                                       PRBM  10
C     ..................................................................PRBM  20
C                                                                       PRBM  30
C        SUBROUTINE PRBM                                                PRBM  40
C                                                                       PRBM  50
C        PURPOSE                                                        PRBM  60
C           TO CALCULATE ALL REAL AND COMPLEX ROOTS OF A GIVEN          PRBM  70
C           POLYNOMIAL WITH REAL COEFFICIENTS.                          PRBM  80
C                                                                       PRBM  90
C        USAGE                                                          PRBM 100
C           CALL PRBM (C,IC,RR,RC,POL,IR,IER)                           PRBM 110
C                                                                       PRBM 120
C        DESCRIPTION OF PARAMETERS                                      PRBM 130
C           C      - INPUT VECTOR CONTAINING THE COEFFICIENTS OF THE    PRBM 140
C                    GIVEN POLYNOMIAL. COEFFICIENTS ARE ORDERED FROM    PRBM 150
C                    LOW TO HIGH. ON RETURN COEFFICIENTS ARE DIVIDED    PRBM 160
C                    BY THE LAST NONZERO TERM.                          PRBM 170
C           IC     - DIMENSION OF VECTORS C, RR, RC, AND POL.           PRBM 180
C           RR     - RESULTANT VECTOR OF REAL PARTS OF THE ROOTS.       PRBM 190
C           RC     - RESULTANT VECTOR OF COMPLEX PARTS OF THE ROOTS.    PRBM 200
C           POL    - RESULTANT VECTOR OF COEFFICIENTS OF THE POLYNOMIAL PRBM 210
C                    WITH CALCULATED ROOTS. COEFFICIENTS ARE ORDERED    PRBM 220
C                    FROM LOW TO HIGH (SEE REMARK 4).                   PRBM 230
C           IR     - OUTPUT VALUE SPECIFYING THE NUMBER OF CALCULATED   PRBM 240
C                    ROOTS. NORMALLY IR IS EQUAL TO IC-1.               PRBM 250
C           IER    - RESULTANT ERROR PARAMETER CODED AS FOLLOWS         PRBM 260
C                     IER=0  - NO ERROR,                                PRBM 270
C                     IER=1  - SUBROUTINE PQFB RECORDS POOR CONVERGENCE PRBM 280
C                              AT SOME QUADRATIC FACTORIZATION WITHIN   PRBM 290
C                              50 ITERATION STEPS,                      PRBM 300
C                     IER=2  - POLYNOMIAL IS DEGENERATE, I.E. ZERO OR   PRBM 310
C                              CONSTANT,                                PRBM 320
C                              OR OVERFLOW IN NORMALIZATION OF GIVEN    PRBM 330
C                              POLYNOMIAL,                              PRBM 340
C                     IER=3  - THE SUBROUTINE IS BYPASSED DUE TO        PRBM 350
C                              SUCCESSIVE ZERO DIVISORS OR OVERFLOWS    PRBM 360
C                              IN QUADRATIC FACTORIZATION OR DUE TO     PRBM 370
C                              COMPLETELY UNSATISFACTORY ACCURACY,      PRBM 380
C                     IER=-1 - CALCULATED COEFFICIENT VECTOR HAS LESS   PRBM 390
C                              THAN THREE CORRECT SIGNIFICANT DIGITS.   PRBM 400
C                              THIS REVEALS POOR ACCURACY OF CALCULATED PRBM 410
C                              ROOTS.                                   PRBM 420
C                                                                       PRBM 430
C        REMARKS                                                        PRBM 440
C           (1) REAL PARTS OF THE ROOTS ARE STORED IN RR(1) UP TO RR(IR)PRBM 450
C               AND CORRESPONDING COMPLEX PARTS IN RC(1) UP TO RC(IR).  PRBM 460
C           (2) ERROR MESSAGE IER=1 INDICATES POOR CONVERGENCE WITHIN   PRBM 470
C               50 ITERATION STEPS AT SOME QUADRQTIC FACTORIZATION      PRBM 480
C               PERFORMED BY SUBROUTINE PQFB.                           PRBM 490
C           (3) NO ACTION BESIDES ERROR MESSAGE IER=2 IN CASE OF A ZERO PRBM 500
C               OR CONSTANT POLYNOMIAL. THE SAME ERROR MESSAGE IS GIVEN PRBM 510
C               IN CASE OF AN OVERFLOW IN NORMALIZATION OF GIVEN        PRBM 520
C               POLYNOMIAL.                                             PRBM 530
C           (4) ERROR MESSAGE IER=3 INDICATES SUCCESSIVE ZERO DIVISORS  PRBM 540
C               OR OVERFLOWS OR COMPLETELY UNSATISFACTORY ACCURACY AT   PRBM 550
C               ANY QUADRATIC FACTORIZATION PERFORMED BY                PRBM 560
C               SUBROUTINE PQFB. IN THIS CASE CALCULATION IS BYPASSED.  PRBM 570
C               IR RECORDS THE NUMBER OF CALCULATED ROOTS.              PRBM 580
C               POL(1),...,POL(J-IR) ARE THE COEFFICIENTS OF THE        PRBM 590
C               REMAINING POLYNOMIAL, WHERE J IS THE ACTUAL NUMBER OF   PRBM 600
C               COEFFICIENTS IN VECTOR C (NORMALLY J=IC).               PRBM 610
C           (5) IF CALCULATED COEFFICIENT VECTOR HAS LESS THAN THREE    PRBM 620
C               CORRECT SIGNIFICANT DIGITS THOUGH ALL QUADRATIC         PRBM 630
C               FACTORIZATIONS SHOWED SATISFACTORY ACCURACY, THE ERROR  PRBM 640
C               MESSAGE IER=-1 IS GIVEN.                                PRBM 650
C           (6) THE FINAL COMPARISON BETWEEN GIVEN AND CALCULATED       PRBM 660
C               COEFFICIENT VECTOR IS PERFORMED ONLY IF ALL ROOTS HAVE  PRBM 670
C               BEEN CALCULATED. IN THIS CASE THE NUMBER OF ROOTS IR IS PRBM 680
C               EQUAL TO THE ACTUAL DEGREE OF THE POLYNOMIAL (NORMALLY  PRBM 690
C               IR=IC-1). THE MAXIMAL RELATIVE ERROR OF THE COEFFICIENT PRBM 700
C               VECTOR IS RECORDED IN RR(IR+1).                         PRBM 710
C                                                                       PRBM 720
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PRBM 730
C           SUBROUTINE PQFB     QUADRATIC FACTORIZATION OF A POLYNOMIAL PRBM 740
C                               BY BAIRSTOW ITERATION.                  PRBM 750
C                                                                       PRBM 760
C        METHOD                                                         PRBM 770
C           THE ROOTS OF THE POLYNOMIAL ARE CALCULATED BY MEANS OF      PRBM 780
C           SUCCESSIVE QUADRATIC FACTORIZATION PERFORMED BY BAIRSTOW    PRBM 790
C           ITERATION. X**2 IS USED AS INITIAL GUESS FOR THE FIRST      PRBM 800
C           QUADRATIC FACTOR, AND FURTHER EACH CALCULATED QUADRATIC     PRBM 810
C           FACTOR IS USED AS INITIAL GUESS FOR THE NEXT ONE. AFTER     PRBM 820
C           COMPUTATION OF ALL ROOTS THE COEFFICIENT VECTOR IS          PRBM 830
C           CALCULATED AND COMPARED WITH THE GIVEN ONE.                 PRBM 840
C           FOR REFERENCE, SEE J. H. WILKINSON, THE EVALUATION OF THE   PRBM 850
C           ZEROS OF ILL-CONDITIONED POLYNOMIALS (PART ONE AND TWO),    PRBM 860
C           NUMERISCHE MATHEMATIK, VOL.1 (1959), PP.150-180.            PRBM 870
C                                                                       PRBM 880
C     ..................................................................PRBM 890
C                                                                       PRBM 900
      SUBROUTINE PRBM(C,IC,RR,RC,POL,IR,IER)                            PRBM 910
C                                                                       PRBM 920
C                                                                       PRBM 930
      DIMENSION C(1),RR(1),RC(1),POL(1),Q(4)                            PRBM 940
C                                                                       PRBM 950
C        TEST ON LEADING ZERO COEFFICIENTS                              PRBM 960
      EPS=1.E-3                                                         PRBM 970
      LIM=50                                                            PRBM 980
      IR=IC+1                                                           PRBM 990
    1 IR=IR-1                                                           PRBM1000
      IF(IR-1)42,42,2                                                   PRBM1010
    2 IF(C(IR))3,1,3                                                    PRBM1020
C                                                                       PRBM1030
C        WORK UP ZERO ROOTS AND NORMALIZE REMAINING POLYNOMIAL          PRBM1040
    3 IER=0                                                             PRBM1050
      J=IR                                                              PRBM1060
      L=0                                                               PRBM1070
      A=C(IR)                                                           PRBM1080
      DO 8 I=1,IR                                                       PRBM1090
      IF(L)4,4,7                                                        PRBM1100
    4 IF(C(I))6,5,6                                                     PRBM1110
    5 RR(I)=0.                                                          PRBM1120
      RC(I)=0.                                                          PRBM1130
      POL(J)=0.                                                         PRBM1140
      J=J-1                                                             PRBM1150
      GO TO 8                                                           PRBM1160
    6 L=1                                                               PRBM1170
      IST=I                                                             PRBM1180
      J=0                                                               PRBM1190
    7 J=J+1                                                             PRBM1200
      C(I)=C(I)/A                                                       PRBM1210
      POL(J)=C(I)                                                       PRBM1220
      CALL OVERFL(N)                                                    PRBM1230
      IF(N-2)42,8,8                                                     PRBM1240
    8 CONTINUE                                                          PRBM1250
C                                                                       PRBM1260
C        START BAIRSTOW ITERATION                                       PRBM1270
      Q1=0.                                                             PRBM1280
      Q2=0.                                                             PRBM1290
    9 IF(J-2)33,10,14                                                   PRBM1300
C                                                                       PRBM1310
C        DEGREE OF RESTPOLYNOMIAL IS EQUAL TO ONE                       PRBM1320
   10 A=POL(1)                                                          PRBM1330
      RR(IST)=-A                                                        PRBM1340
      RC(IST)=0.                                                        PRBM1350
      IR=IR-1                                                           PRBM1360
      Q2=0.                                                             PRBM1370
      IF(IR-1)13,13,11                                                  PRBM1380
   11 DO 12 I=2,IR                                                      PRBM1390
      Q1=Q2                                                             PRBM1400
      Q2=POL(I+1)                                                       PRBM1410
   12 POL(I)=A*Q2+Q1                                                    PRBM1420
   13 POL(IR+1)=A+Q2                                                    PRBM1430
      GO TO 34                                                          PRBM1440
C        THIS IS BRANCH TO COMPARISON OF COEFFICIENT VECTORS C AND POL  PRBM1450
C                                                                       PRBM1460
C        DEGREE OF RESTPOLYNOMIAL IS GREATER THAN ONE                   PRBM1470
   14 DO 22 L=1,10                                                      PRBM1480
      N=1                                                               PRBM1490
   15 Q(1)=Q1                                                           PRBM1500
      Q(2)=Q2                                                           PRBM1510
      CALL PQFB(POL,J,Q,LIM,I)                                          PRBM1520
      IF(I)16,24,23                                                     PRBM1530
   16 IF(Q1)18,17,18                                                    PRBM1540
   17 IF(Q2)18,21,18                                                    PRBM1550
   18 GO TO (19,20,19,21),N                                             PRBM1560
   19 Q1=-Q1                                                            PRBM1570
      N=N+1                                                             PRBM1580
      GO TO 15                                                          PRBM1590
   20 Q2=-Q2                                                            PRBM1600
      N=N+1                                                             PRBM1610
      GO TO 15                                                          PRBM1620
   21 Q1=1.+Q1                                                          PRBM1630
   22 Q2=1.-Q2                                                          PRBM1640
C                                                                       PRBM1650
C        ERROR EXIT DUE TO UNSATISFACTORY RESULTS OF FACTORIZATION      PRBM1660
      IER=3                                                             PRBM1670
      IR=IR-J                                                           PRBM1680
      RETURN                                                            PRBM1690
C                                                                       PRBM1700
C        WORK UP RESULTS OF QUADRATIC FACTORIZATION                     PRBM1710
   23 IER=1                                                             PRBM1720
   24 Q1=Q(1)                                                           PRBM1730
      Q2=Q(2)                                                           PRBM1740
C                                                                       PRBM1750
C        PERFORM DIVISION OF FACTORIZED POLYNOMIAL BY QUADRATIC FACTOR  PRBM1760
      B=0.                                                              PRBM1770
      A=0.                                                              PRBM1780
      I=J                                                               PRBM1790
   25 H=-Q1*B-Q2*A+POL(I)                                               PRBM1800
      POL(I)=B                                                          PRBM1810
      B=A                                                               PRBM1820
      A=H                                                               PRBM1830
      I=I-1                                                             PRBM1840
      IF(I-2)26,26,25                                                   PRBM1850
   26 POL(2)=B                                                          PRBM1860
      POL(1)=A                                                          PRBM1870
C                                                                       PRBM1880
C        MULTIPLY POLYNOMIAL WITH CALCULATED ROOTS BY QUADRATIC FACTOR  PRBM1890
      L=IR-1                                                            PRBM1900
      IF(J-L)27,27,29                                                   PRBM1910
   27 DO 28 I=J,L                                                       PRBM1920
   28 POL(I-1)=POL(I-1)+POL(I)*Q2+POL(I+1)*Q1                           PRBM1930
   29 POL(L)=POL(L)+POL(L+1)*Q2+Q1                                      PRBM1940
      POL(IR)=POL(IR)+Q2                                                PRBM1950
C                                                                       PRBM1960
C        CALCULATE ROOT-PAIR FROM QUADRATIC FACTOR X*X+Q2*X+Q1          PRBM1970
      H=-.5*Q2                                                          PRBM1980
      A=H*H-Q1                                                          PRBM1990
      B=SQRT(ABS(A))                                                    PRBM2000
      IF(A)30,30,31                                                     PRBM2010
   30 RR(IST)=H                                                         PRBM2020
      RC(IST)=B                                                         PRBM2030
      IST=IST+1                                                         PRBM2040
      RR(IST)=H                                                         PRBM2050
      RC(IST)=-B                                                        PRBM2060
      GO TO 32                                                          PRBM2070
   31 B=H+SIGN(B,H)                                                     PRBM2080
      RR(IST)=Q1/B                                                      PRBM2090
      RC(IST)=0.                                                        PRBM2100
      IST=IST+1                                                         PRBM2110
      RR(IST)=B                                                         PRBM2120
      RC(IST)=0.                                                        PRBM2130
   32 IST=IST+1                                                         PRBM2140
      J=J-2                                                             PRBM2150
      GO TO 9                                                           PRBM2160
C                                                                       PRBM2170
C        SHIFT BACK ELEMENTS OF POL BY 1 AND COMPARE VECTORS POL AND C  PRBM2180
   33 IR=IR-1                                                           PRBM2190
   34 A=0.                                                              PRBM2200
      DO 38 I=1,IR                                                      PRBM2210
      Q1=C(I)                                                           PRBM2220
      Q2=POL(I+1)                                                       PRBM2230
      POL(I)=Q2                                                         PRBM2240
      IF(Q1)35,36,35                                                    PRBM2250
   35 Q2=(Q1-Q2)/Q1                                                     PRBM2260
   36 Q2=ABS(Q2)                                                        PRBM2270
      IF(Q2-A)38,38,37                                                  PRBM2280
   37 A=Q2                                                              PRBM2290
   38 CONTINUE                                                          PRBM2300
      I=IR+1                                                            PRBM2310
      POL(I)=1.                                                         PRBM2320
      RR(I)=A                                                           PRBM2330
      RC(I)=0.                                                          PRBM2340
      IF(IER)39,39,41                                                   PRBM2350
   39 IF(A-EPS)41,41,40                                                 PRBM2360
C                                                                       PRBM2370
C        WARNING DUE TO POOR ACCURACY OF CALCULATED COEFFICIENT VECTOR  PRBM2380
   40 IER=-1                                                            PRBM2390
   41 RETURN                                                            PRBM2400
C                                                                       PRBM2410
C        ERROR EXIT DUE TO DEGENERATE POLYNOMIAL OR OVERFLOW IN         PRBM2420
C        NORMALIZATION                                                  PRBM2430
   42 IER=2                                                             PRBM2440
      IR=0                                                              PRBM2450
      RETURN                                                            PRBM2460
      END                                                               PRBM2470
