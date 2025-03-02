C                                                                       PRQD  10
C     ..................................................................PRQD  20
C                                                                       PRQD  30
C        SUBROUTINE PRQD                                                PRQD  40
C                                                                       PRQD  50
C        PURPOSE                                                        PRQD  60
C           CALCULATE ALL REAL AND COMPLEX ROOTS OF A GIVEN POLYNOMIAL  PRQD  70
C           WITH REAL COEFFICIENTS.                                     PRQD  80
C                                                                       PRQD  90
C        USAGE                                                          PRQD 100
C           CALL PRQD(C,IC,Q,E,POL,IR,IER)                              PRQD 110
C                                                                       PRQD 120
C        DESCRIPTION OF PARAMETERS                                      PRQD 130
C           C     - COEFFICIENT VECTOR OF GIVEN POLYNOMIAL              PRQD 140
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           PRQD 150
C                   THE GIVEN COEFFICIENT VECTOR GETS DIVIDED BY THE    PRQD 160
C                   LAST NONZERO TERM                                   PRQD 170
C           IC    - DIMENSION OF VECTOR C                               PRQD 180
C           Q     - WORKING STORAGE OF DIMENSION IC                     PRQD 190
C                   ON RETURN Q CONTAINS REAL PARTS OF ROOTS            PRQD 200
C           E     - WORKING STORAGE OF DIMENSION IC                     PRQD 210
C                   ON RETURN E CONTAINS COMPLEX PARTS OF ROOTS         PRQD 220
C           POL   - WORKING STORAGE OF DIMENSION IC                     PRQD 230
C                   ON RETURN POL CONTAINS THE COEFFICIENTS OF THE      PRQD 240
C                   POLYNOMIAL WITH CALCULATED ROOTS                    PRQD 250
C                   THIS RESULTING COEFFICIENT VECTOR HAS DIMENSION IR+1PRQD 260
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           PRQD 270
C           IR    - NUMBER OF CALCULATED ROOTS                          PRQD 280
C                   NORMALLY IR IS EQUAL TO DIMENSION IC MINUS ONE      PRQD 290
C           IER   - RESULTING ERROR PARAMETER. SEE REMARKS              PRQD 300
C                                                                       PRQD 310
C        REMARKS                                                        PRQD 320
C           THE REAL PART OF THE ROOTS IS STORED IN Q(1) UP TO Q(IR)    PRQD 330
C           CORRESPONDING COMPLEX PARTS ARE STORED IN E(1) UP TO E(IR). PRQD 340
C           IER = 0 MEANS NO ERRORS                                     PRQD 350
C           IER = 1 MEANS NO CONVERGENCE WITH FEASIBLE TOLERANCE        PRQD 360
C           IER = 2 MEANS POLYNOMIAL IS DEGENERATE (CONSTANT OR ZERO)   PRQD 370
C           IER = 3 MEANS SUBROUTINE WAS ABANDONED DUE TO ZERO DIVISOR  PRQD 380
C           IER = 4 MEANS THERE EXISTS NO S-FRACTION                    PRQD 390
C           IER =-1 MEANS CALCULATED COEFFICIENT VECTOR REVEALS POOR    PRQD 400
C                   ACCURACY OF THE CALCULATED ROOTS.                   PRQD 410
C                   THE CALCULATED COEFFICIENT VECTOR HAS LESS THAN     PRQD 420
C                   3 CORRECT DIGITS.                                   PRQD 430
C           THE FINAL COMPARISON BETWEEN GIVEN AND CALCULATED           PRQD 440
C           COEFFICIENT VECTOR IS PERFORMED ONLY IF ALL ROOTS HAVE BEEN PRQD 450
C           CALCULATED.                                                 PRQD 460
C           THE MAXIMAL RELATIVE ERROR OF THE COEFFICIENT VECTOR IS     PRQD 470
C           RECORDED IN Q(IR+1).                                        PRQD 480
C                                                                       PRQD 490
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PRQD 500
C           NONE                                                        PRQD 510
C                                                                       PRQD 520
C        METHOD                                                         PRQD 530
C           THE ROOTS OF THE POLYNOMIAL ARE CALCULATED BY MEANS OF      PRQD 540
C           THE QUOTIENT-DIFFERENCE ALGORITHM WITH DISPLACEMENT.        PRQD 550
C           REFERENCE                                                   PRQD 560
C           H.RUTISHAUSER, DER QUOTIENTEN-DIFFERENZEN-ALGORITHMUS,      PRQD 570
C           BIRKHAEUSER, BASEL/STUTTGART, 1957.                         PRQD 580
C                                                                       PRQD 590
C     ..................................................................PRQD 600
C                                                                       PRQD 610
      SUBROUTINE PRQD(C,IC,Q,E,POL,IR,IER)                              PRQD 620
C                                                                       PRQD 630
C      DIMENSIONED DUMMY VARIABLES                                      PRQD 640
      DIMENSION E(1),Q(1),C(1),POL(1)                                   PRQD 650
C                                                                       PRQD 660
C        NORMALIZATION OF GIVEN POLYNOMIAL                              PRQD 670
C           TEST OF DIMENSION                                           PRQD 680
C        IR CONTAINS INDEX OF HIGHEST COEFFICIENT                       PRQD 690
      IER=0                                                             PRQD 700
      IR=IC                                                             PRQD 710
      EPS=1.E-6                                                         PRQD 720
      TOL=1.E-3                                                         PRQD 730
      LIMIT=10*IC                                                       PRQD 740
      KOUNT=0                                                           PRQD 750
    1 IF(IR-1)79,79,2                                                   PRQD 760
C                                                                       PRQD 770
C        DROP TRAILING ZERO COEFFICIENTS                                PRQD 780
    2 IF(C(IR))4,3,4                                                    PRQD 790
    3 IR=IR-1                                                           PRQD 800
      GOTO 1                                                            PRQD 810
C                                                                       PRQD 820
C           REARRANGEMENT OF GIVEN POLYNOMIAL                           PRQD 830
C        EXTRACTION OF ZERO ROOTS                                       PRQD 840
    4 O=1./C(IR)                                                        PRQD 850
      IEND=IR-1                                                         PRQD 860
      ISTA=1                                                            PRQD 870
      NSAV=IR+1                                                         PRQD 880
      JBEG=1                                                            PRQD 890
C                                                                       PRQD 900
C        Q(J)=1.                                                        PRQD 910
C        Q(J+I)=C(IR-I)/C(IR)                                           PRQD 920
C        Q(IR)=C(J)/C(IR)                                               PRQD 930
C        WHERE J IS THE INDEX OF THE LOWEST NONZERO COEFFICIENT         PRQD 940
      DO 9 I=1,IR                                                       PRQD 950
      J=NSAV-I                                                          PRQD 960
      IF(C(I))7,5,7                                                     PRQD 970
    5 GOTO(6,8),JBEG                                                    PRQD 980
    6 NSAV=NSAV+1                                                       PRQD 990
      Q(ISTA)=0.                                                        PRQD1000
      E(ISTA)=0.                                                        PRQD1010
      ISTA=ISTA+1                                                       PRQD1020
      GOTO 9                                                            PRQD1030
    7 JBEG=2                                                            PRQD1040
    8 Q(J)=C(I)*O                                                       PRQD1050
      C(I)=Q(J)                                                         PRQD1060
    9 CONTINUE                                                          PRQD1070
C                                                                       PRQD1080
C           INITIALIZATION                                              PRQD1090
      ESAV=0.                                                           PRQD1100
      Q(ISTA)=0.                                                        PRQD1110
   10 NSAV=IR                                                           PRQD1120
C                                                                       PRQD1130
C        COMPUTATION OF DERIVATIVE                                      PRQD1140
      EXPT=IR-ISTA                                                      PRQD1150
      E(ISTA)=EXPT                                                      PRQD1160
      DO 11 I=ISTA,IEND                                                 PRQD1170
      EXPT=EXPT-1.0                                                     PRQD1180
      POL(I+1)=EPS*ABS(Q(I+1))+EPS                                      PRQD1190
   11 E(I+1)=Q(I+1)*EXPT                                                PRQD1200
C                                                                       PRQD1210
C        TEST OF REMAINING DIMENSION                                    PRQD1220
      IF(ISTA-IEND)12,20,60                                             PRQD1230
   12 JEND=IEND-1                                                       PRQD1240
C                                                                       PRQD1250
C        COMPUTATION OF S-FRACTION                                      PRQD1260
      DO 19 I=ISTA,JEND                                                 PRQD1270
      IF(I-ISTA)13,16,13                                                PRQD1280
   13 IF(ABS(E(I))-POL(I+1))14,14,16                                    PRQD1290
C                                                                       PRQD1300
C        THE GIVEN POLYNOMIAL HAS MULTIPLE ROOTS, THE COEFFICIENTS OF   PRQD1310
C        THE COMMON FACTOR ARE STORED FROM Q(NSAV) UP TO Q(IR)          PRQD1320
   14 NSAV=I                                                            PRQD1330
      DO 15 K=I,JEND                                                    PRQD1340
      IF(ABS(E(K))-POL(K+1))15,15,80                                    PRQD1350
   15 CONTINUE                                                          PRQD1360
      GOTO 21                                                           PRQD1370
C                                                                       PRQD1380
C           EUCLIDEAN ALGORITHM                                         PRQD1390
   16 DO 19 K=I,IEND                                                    PRQD1400
      E(K+1)=E(K+1)/E(I)                                                PRQD1410
      Q(K+1)=E(K+1)-Q(K+1)                                              PRQD1420
      IF(K-I)18,17,18                                                   PRQD1430
C                                                                       PRQD1440
C        TEST FOR SMALL DIVISOR                                         PRQD1450
   17 IF(ABS(Q(I+1))-POL(I+1))80,80,19                                  PRQD1460
   18 Q(K+1)=Q(K+1)/Q(I+1)                                              PRQD1470
      POL(K+1)=POL(K+1)/ABS(Q(I+1))                                     PRQD1480
      E(K)=Q(K+1)-E(K)                                                  PRQD1490
   19 CONTINUE                                                          PRQD1500
   20 Q(IR)=-Q(IR)                                                      PRQD1510
C                                                                       PRQD1520
C           THE DISPLACEMENT EXPT IS SET TO 0 AUTOMATICALLY.            PRQD1530
C           E(ISTA)=0.,Q(ISTA+1),...,E(NSAV-1),Q(NSAV),E(NSAV)=0.,      PRQD1540
C           FORM A DIAGONAL OF THE QD-ARRAY.                            PRQD1550
C        INITIALIZATION OF BOUNDARY VALUES                              PRQD1560
   21 E(ISTA)=0.                                                        PRQD1570
      NRAN=NSAV-1                                                       PRQD1580
   22 E(NRAN+1)=0.                                                      PRQD1590
C                                                                       PRQD1600
C           TEST FOR LINEAR OR CONSTANT FACTOR                          PRQD1610
C        NRAN-ISTA IS DEGREE-1                                          PRQD1620
      IF(NRAN-ISTA)24,23,31                                             PRQD1630
C                                                                       PRQD1640
C        LINEAR FACTOR                                                  PRQD1650
   23 Q(ISTA+1)=Q(ISTA+1)+EXPT                                          PRQD1660
      E(ISTA+1)=0.                                                      PRQD1670
C                                                                       PRQD1680
C        TEST FOR UNFACTORED COMMON DIVISOR                             PRQD1690
   24 E(ISTA)=ESAV                                                      PRQD1700
      IF(IR-NSAV)60,60,25                                               PRQD1710
C                                                                       PRQD1720
C        INITIALIZE QD-ALGORITHM FOR COMMON DIVISOR                     PRQD1730
   25 ISTA=NSAV                                                         PRQD1740
      ESAV=E(ISTA)                                                      PRQD1750
      GOTO 10                                                           PRQD1760
C                                                                       PRQD1770
C        COMPUTATION OF ROOT PAIR                                       PRQD1780
   26 P=P+EXPT                                                          PRQD1790
C                                                                       PRQD1800
C           TEST FOR REALITY                                            PRQD1810
      IF(O)27,28,28                                                     PRQD1820
C                                                                       PRQD1830
C           COMPLEX ROOT PAIR                                           PRQD1840
   27 Q(NRAN)=P                                                         PRQD1850
      Q(NRAN+1)=P                                                       PRQD1860
      E(NRAN)=T                                                         PRQD1870
      E(NRAN+1)=-T                                                      PRQD1880
      GOTO 29                                                           PRQD1890
C                                                                       PRQD1900
C           REAL ROOT PAIR                                              PRQD1910
   28 Q(NRAN)=P-T                                                       PRQD1920
      Q(NRAN+1)=P+T                                                     PRQD1930
      E(NRAN)=0.                                                        PRQD1940
C                                                                       PRQD1950
C           REDUCTION OF DEGREE BY 2 (DEFLATION)                        PRQD1960
   29 NRAN=NRAN-2                                                       PRQD1970
      GOTO 22                                                           PRQD1980
C                                                                       PRQD1990
C        COMPUTATION OF REAL ROOT                                       PRQD2000
   30 Q(NRAN+1)=EXPT+P                                                  PRQD2010
C                                                                       PRQD2020
C           REDUCTION OF DEGREE BY 1 (DEFLATION)                        PRQD2030
      NRAN=NRAN-1                                                       PRQD2040
      GOTO 22                                                           PRQD2050
C                                                                       PRQD2060
C        START QD-ITERATION                                             PRQD2070
   31 JBEG=ISTA+1                                                       PRQD2080
      JEND=NRAN-1                                                       PRQD2090
      TEPS=EPS                                                          PRQD2100
      TDELT=1.E-2                                                       PRQD2110
   32 KOUNT=KOUNT+1                                                     PRQD2120
      P=Q(NRAN+1)                                                       PRQD2130
      R=ABS(E(NRAN))                                                    PRQD2140
C                                                                       PRQD2150
C           TEST FOR CONVERGENCE                                        PRQD2160
      IF(R-TEPS)30,30,33                                                PRQD2170
   33 S=ABS(E(JEND))                                                    PRQD2180
C                                                                       PRQD2190
C        IS THERE A REAL ROOT NEXT                                      PRQD2200
      IF(S-R)38,38,34                                                   PRQD2210
C                                                                       PRQD2220
C        IS DISPLACEMENT SMALL ENOUGH                                   PRQD2230
   34 IF(R-TDELT)36,35,35                                               PRQD2240
   35 P=0.                                                              PRQD2250
   36 O=P                                                               PRQD2260
      DO 37 J=JBEG,NRAN                                                 PRQD2270
      Q(J)=Q(J)+E(J)-E(J-1)-O                                           PRQD2280
C                                                                       PRQD2290
C           TEST FOR SMALL DIVISOR                                      PRQD2300
      IF(ABS(Q(J))-POL(J))81,81,37                                      PRQD2310
   37 E(J)=Q(J+1)*E(J)/Q(J)                                             PRQD2320
      Q(NRAN+1)=-E(NRAN)+Q(NRAN+1)-O                                    PRQD2330
      GOTO 54                                                           PRQD2340
C                                                                       PRQD2350
C        CALCULATE DISPLACEMENT FOR DOUBLE ROOTS                        PRQD2360
C           QUADRATIC EQUATION FOR DOUBLE ROOTS                         PRQD2370
C           X**2-(Q(NRAN)+Q(NRAN+1)+E(NRAN))*X+Q(NRAN)*Q(NRAN+1)=0      PRQD2380
   38 P=0.5*(Q(NRAN)+E(NRAN)+Q(NRAN+1))                                 PRQD2390
      O=P*P-Q(NRAN)*Q(NRAN+1)                                           PRQD2400
      T=SQRT(ABS(O))                                                    PRQD2410
C                                                                       PRQD2420
C        TEST FOR CONVERGENCE                                           PRQD2430
      IF(S-TEPS)26,26,39                                                PRQD2440
C                                                                       PRQD2450
C        ARE THERE COMPLEX ROOTS                                        PRQD2460
   39 IF(O)43,40,40                                                     PRQD2470
   40 IF(P)42,41,41                                                     PRQD2480
   41 T=-T                                                              PRQD2490
   42 P=P+T                                                             PRQD2500
      R=S                                                               PRQD2510
      GOTO 34                                                           PRQD2520
C                                                                       PRQD2530
C        MODIFICATION FOR COMPLEX ROOTS                                 PRQD2540
C        IS DISPLACEMENT SMALL ENOUGH                                   PRQD2550
   43 IF(S-TDELT)44,35,35                                               PRQD2560
C                                                                       PRQD2570
C        INITIALIZATION                                                 PRQD2580
   44 O=Q(JBEG)+E(JBEG)-P                                               PRQD2590
C                                                                       PRQD2600
C        TEST FOR SMALL DIVISOR                                         PRQD2610
      IF(ABS(O)-POL(JBEG))81,81,45                                      PRQD2620
   45 T=(T/O)**2                                                        PRQD2630
      U=E(JBEG)*Q(JBEG+1)/(O*(1.+T))                                    PRQD2640
      V=O+U                                                             PRQD2650
      KOUNT=KOUNT+2                                                     PRQD2660
C                                                                       PRQD2670
C        THREEFOLD LOOP FOR COMPLEX DISPLACEMENT                        PRQD2680
      DO 53 J=JBEG,NRAN                                                 PRQD2690
      O=Q(J+1)+E(J+1)-U-P                                               PRQD2700
C                                                                       PRQD2710
C        TEST FOR SMALL DIVISOR                                         PRQD2720
      IF(ABS(V)-POL(J))46,46,49                                         PRQD2730
   46 IF(J-NRAN)81,47,81                                                PRQD2740
   47 EXPT=EXPT+P                                                       PRQD2750
      IF(ABS(E(JEND))-TOL)48,48,81                                      PRQD2760
   48 P=0.5*(V+O-E(JEND))                                               PRQD2770
      O=P*P-(V-U)*(O-U*T-O*W*(1.+T)/Q(JEND))                            PRQD2780
      T=SQRT(ABS(O))                                                    PRQD2790
      GOTO 26                                                           PRQD2800
C                                                                       PRQD2810
C           TEST FOR SMALL DIVISOR                                      PRQD2820
   49 IF(ABS(O)-POL(J+1))46,46,50                                       PRQD2830
   50 W=U*O/V                                                           PRQD2840
      T=T*(V/O)**2                                                      PRQD2850
      Q(J)=V+W-E(J-1)                                                   PRQD2860
      U=0.                                                              PRQD2870
      IF(J-NRAN)51,52,52                                                PRQD2880
   51 U=Q(J+2)*E(J+1)/(O*(1.+T))                                        PRQD2890
   52 V=O+U-W                                                           PRQD2900
C                                                                       PRQD2910
C        TEST FOR SMALL DIVISOR                                         PRQD2920
      IF(ABS(Q(J))-POL(J))81,81,53                                      PRQD2930
   53 E(J)=W*V*(1.+T)/Q(J)                                              PRQD2940
      Q(NRAN+1)=V-E(NRAN)                                               PRQD2950
   54 EXPT=EXPT+P                                                       PRQD2960
      TEPS=TEPS*1.1                                                     PRQD2970
      TDELT=TDELT*1.1                                                   PRQD2980
      IF(KOUNT-LIMIT)32,55,55                                           PRQD2990
C                                                                       PRQD3000
C        NO CONVERGENCE WITH FEASIBLE TOLERANCE                         PRQD3010
C           ERROR RETURN IN CASE OF UNSATISFACTORY CONVERGENCE          PRQD3020
   55 IER=1                                                             PRQD3030
C                                                                       PRQD3040
C        REARRANGE CALCULATED ROOTS                                     PRQD3050
   56 IEND=NSAV-NRAN-1                                                  PRQD3060
      E(ISTA)=ESAV                                                      PRQD3070
      IF(IEND)59,59,57                                                  PRQD3080
   57 DO 58 I=1,IEND                                                    PRQD3090
      J=ISTA+I                                                          PRQD3100
      K=NRAN+1+I                                                        PRQD3110
      E(J)=E(K)                                                         PRQD3120
   58 Q(J)=Q(K)                                                         PRQD3130
   59 IR=ISTA+IEND                                                      PRQD3140
C                                                                       PRQD3150
C        NORMAL RETURN                                                  PRQD3160
   60 IR=IR-1                                                           PRQD3170
      IF(IR)78,78,61                                                    PRQD3180
C                                                                       PRQD3190
C        REARRANGE CALCULATED ROOTS                                     PRQD3200
   61 DO 62 I=1,IR                                                      PRQD3210
      Q(I)=Q(I+1)                                                       PRQD3220
   62 E(I)=E(I+1)                                                       PRQD3230
C                                                                       PRQD3240
C        CALCULATE COEFFICIENT VECTOR FROM ROOTS                        PRQD3250
      POL(IR+1)=1.                                                      PRQD3260
      IEND=IR-1                                                         PRQD3270
      JBEG=1                                                            PRQD3280
      DO 69 J=1,IR                                                      PRQD3290
      ISTA=IR+1-J                                                       PRQD3300
      O=0.                                                              PRQD3310
      P=Q(ISTA)                                                         PRQD3320
      T=E(ISTA)                                                         PRQD3330
      IF(T)65,63,65                                                     PRQD3340
C                                                                       PRQD3350
C        MULTIPLY WITH LINEAR FACTOR                                    PRQD3360
   63 DO 64 I=ISTA,IR                                                   PRQD3370
      POL(I)=O-P*POL(I+1)                                               PRQD3380
   64 O=POL(I+1)                                                        PRQD3390
      GOTO 69                                                           PRQD3400
   65 GOTO(66,67),JBEG                                                  PRQD3410
   66 JBEG=2                                                            PRQD3420
      POL(ISTA)=0.                                                      PRQD3430
      GOTO 69                                                           PRQD3440
C                                                                       PRQD3450
C        MULTIPLY WITH QUADRATIC FACTOR                                 PRQD3460
   67 JBEG=1                                                            PRQD3470
      U=P*P+T*T                                                         PRQD3480
      P=P+P                                                             PRQD3490
      DO 68 I=ISTA,IEND                                                 PRQD3500
      POL(I)=O-P*POL(I+1)+U*POL(I+2)                                    PRQD3510
   68 O=POL(I+1)                                                        PRQD3520
      POL(IR)=O-P                                                       PRQD3530
   69 CONTINUE                                                          PRQD3540
      IF(IER)78,70,78                                                   PRQD3550
C                                                                       PRQD3560
C        COMPARISON OF COEFFICIENT VECTORS, IE. TEST OF ACCURACY        PRQD3570
   70 P=0.                                                              PRQD3580
      DO 75 I=1,IR                                                      PRQD3590
      IF(C(I))72,71,72                                                  PRQD3600
   71 O=ABS(POL(I))                                                     PRQD3610
      GOTO 73                                                           PRQD3620
   72 O=ABS((POL(I)-C(I))/C(I))                                         PRQD3630
   73 IF(P-O)74,75,75                                                   PRQD3640
   74 P=O                                                               PRQD3650
   75 CONTINUE                                                          PRQD3660
      IF(P-TOL)77,76,76                                                 PRQD3670
   76 IER=-1                                                            PRQD3680
   77 Q(IR+1)=P                                                         PRQD3690
      E(IR+1)=0.                                                        PRQD3700
   78 RETURN                                                            PRQD3710
C                                                                       PRQD3720
C        ERROR RETURNS                                                  PRQD3730
C           ERROR RETURN FOR POLYNOMIALS OF DEGREE LESS THAN 1          PRQD3740
   79 IER=2                                                             PRQD3750
      IR=0                                                              PRQD3760
      RETURN                                                            PRQD3770
C                                                                       PRQD3780
C           ERROR RETURN IF THERE EXISTS NO S-FRACTION                  PRQD3790
   80 IER=4                                                             PRQD3800
      IR=ISTA                                                           PRQD3810
      GOTO 60                                                           PRQD3820
C                                                                       PRQD3830
C           ERROR RETURN IN CASE OF INSTABLE QD-ALGORITHM               PRQD3840
   81 IER=3                                                             PRQD3850
      GOTO 56                                                           PRQD3860
      END                                                               PRQD3870
