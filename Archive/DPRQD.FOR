C                                                                       DPRQ  10
C     ..................................................................DPRQ  20
C                                                                       DPRQ  30
C        SUBROUTINE DPRQD                                               DPRQ  40
C                                                                       DPRQ  50
C        PURPOSE                                                        DPRQ  60
C           CALCULATE ALL REAL AND COMPLEX ROOTS OF A GIVEN POLYNOMIAL  DPRQ  70
C           WITH REAL COEFFICIENTS.                                     DPRQ  80
C                                                                       DPRQ  90
C        USAGE                                                          DPRQ 100
C           CALL DPRQD(C,IC,Q,E,POL,IR,IER)                             DPRQ 110
C                                                                       DPRQ 120
C        DESCRIPTION OF PARAMETERS                                      DPRQ 130
C           C     - COEFFICIENT VECTOR OF GIVEN POLYNOMIAL              DPRQ 140
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           DPRQ 150
C                   THE GIVEN COEFFICIENT VECTOR GETS DIVIDED BY THE    DPRQ 160
C                   LAST NONZERO TERM                                   DPRQ 170
C                   DOUBLE PRECISION ARRAY                              DPRQ 180
C           IC    - DIMENSION OF VECTOR C                               DPRQ 190
C           Q     - WORKING STORAGE OF DIMENSION IC                     DPRQ 200
C                   ON RETURN Q CONTAINS REAL PARTS OF ROOTS            DPRQ 210
C                   DOUBLE PRECISION ARRAY                              DPRQ 220
C           E     - WORKING STORAGE OF DIMENSION IC                     DPRQ 230
C                   ON RETURN E CONTAINS COMPLEX PARTS OF ROOTS         DPRQ 240
C                   DOUBLE PRECISION ARRAY                              DPRQ 250
C           POL   - WORKING STORAGE OF DIMENSION IC                     DPRQ 260
C                   ON RETURN POL CONTAINS THE COEFFICIENTS OF THE      DPRQ 270
C                   POLYNOMIAL WITH CALCULATED ROOTS                    DPRQ 280
C                   THIS RESULTING COEFFICIENT VECTOR HAS DIMENSION IR+1DPRQ 290
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           DPRQ 300
C                   DOUBLE PRECISION ARRAY                              DPRQ 310
C           IR    - NUMBER OF CALCULATED ROOTS                          DPRQ 320
C                   NORMALLY IR IS EQUAL TO DIMENSION IC MINUS ONE      DPRQ 330
C           IER   - RESULTING ERROR PARAMETER. SEE REMARKS              DPRQ 340
C                                                                       DPRQ 350
C        REMARKS                                                        DPRQ 360
C           THE REAL PART OF THE ROOTS IS STORED IN Q(1) UP TO Q(IR)    DPRQ 370
C           CORRESPONDING COMPLEX PARTS ARE STORED IN E(1) UP TO E(IR). DPRQ 380
C           IER = 0 MEANS NO ERRORS                                     DPRQ 390
C           IER = 1 MEANS NO CONVERGENCE WITH FEASIBLE TOLERANCE        DPRQ 400
C           IER = 2 MEANS POLYNOMIAL IS DEGENERATE (CONSTANT OR ZERO)   DPRQ 410
C           IER = 3 MEANS SUBROUTINE WAS ABANDONED DUE TO ZERO DIVISOR  DPRQ 420
C           IER = 4 MEANS THERE EXISTS NO S-FRACTION                    DPRQ 430
C           IER =-1 MEANS CALCULATED COEFFICIENT VECTOR REVEALS POOR    DPRQ 440
C                   ACCURACY OF THE CALCULATED ROOTS.                   DPRQ 450
C                   THE CALCULATED COEFFICIENT VECTOR HAS LESS THAN     DPRQ 460
C                   6 CORRECT DIGITS.                                   DPRQ 470
C           THE FINAL COMPARISON BETWEEN GIVEN AND CALCULATED           DPRQ 480
C           COEFFICIENT VECTOR IS PERFORMED ONLY IF ALL ROOTS HAVE BEEN DPRQ 490
C           CALCULATED.                                                 DPRQ 500
C           THE MAXIMAL RELATIVE ERROR OF THE COEFFICIENT VECTOR IS     DPRQ 510
C           RECORDED IN Q(IR+1).                                        DPRQ 520
C                                                                       DPRQ 530
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DPRQ 540
C           NONE                                                        DPRQ 550
C                                                                       DPRQ 560
C        METHOD                                                         DPRQ 570
C           THE ROOTS OF THE POLYNOMIAL ARE CALCULATED BY MEANS OF      DPRQ 580
C           THE QUOTIENT-DIFFERENCE ALGORITHM WITH DISPLACEMENT.        DPRQ 590
C           REFERENCE                                                   DPRQ 600
C           H.RUTISHAUSER, DER QUOTIENTEN-DIFFERENZEN-ALGORITHMUS,      DPRQ 610
C           BIRKHAEUSER, BASEL/STUTTGART, 1957.                         DPRQ 620
C                                                                       DPRQ 630
C     ..................................................................DPRQ 640
C                                                                       DPRQ 650
      SUBROUTINE DPRQD(C,IC,Q,E,POL,IR,IER)                             DPRQ 660
C                                                                       DPRQ 670
C      DIMENSIONED DUMMY VARIABLES                                      DPRQ 680
      DIMENSION E(1),Q(1),C(1),POL(1)                                   DPRQ 690
      DOUBLE PRECISION Q,E,O,P,T,EXPT,ESAV,U,V,W,C,POL,EPS              DPRQ 700
C                                                                       DPRQ 710
C        NORMALIZATION OF GIVEN POLYNOMIAL                              DPRQ 720
C           TEST OF DIMENSION                                           DPRQ 730
C        IR CONTAINS INDEX OF HIGHEST COEFFICIENT                       DPRQ 740
      IR=IC                                                             DPRQ 750
      IER=0                                                             DPRQ 760
      EPS=1.D-16                                                        DPRQ 770
      TOL=1.E-6                                                         DPRQ 780
      LIMIT=10*IC                                                       DPRQ 790
      KOUNT=0                                                           DPRQ 800
    1 IF(IR-1)79,79,2                                                   DPRQ 810
C                                                                       DPRQ 820
C        DROP TRAILING ZERO COEFFICIENTS                                DPRQ 830
    2 IF(C(IR))4,3,4                                                    DPRQ 840
    3 IR=IR-1                                                           DPRQ 850
      GOTO 1                                                            DPRQ 860
C                                                                       DPRQ 870
C           REARRANGEMENT OF GIVEN POLYNOMIAL                           DPRQ 880
C        EXTRACTION OF ZERO ROOTS                                       DPRQ 890
    4 O=1.0D0/C(IR)                                                     DPRQ 900
      IEND=IR-1                                                         DPRQ 910
      ISTA=1                                                            DPRQ 920
      NSAV=IR+1                                                         DPRQ 930
      JBEG=1                                                            DPRQ 940
C                                                                       DPRQ 950
C        Q(J)=1.                                                        DPRQ 960
C        Q(J+I)=C(IR-I)/C(IR)                                           DPRQ 970
C        Q(IR)=C(J)/C(IR)                                               DPRQ 980
C        WHERE J IS THE INDEX OF THE LOWEST NONZERO COEFFICIENT         DPRQ 990
      DO 9 I=1,IR                                                       DPRQ1000
      J=NSAV-I                                                          DPRQ1010
      IF(C(I))7,5,7                                                     DPRQ1020
    5 GOTO(6,8),JBEG                                                    DPRQ1030
    6 NSAV=NSAV+1                                                       DPRQ1040
      Q(ISTA)=0.D0                                                      DPRQ1050
      E(ISTA)=0.D0                                                      DPRQ1060
      ISTA=ISTA+1                                                       DPRQ1070
      GOTO 9                                                            DPRQ1080
    7 JBEG=2                                                            DPRQ1090
    8 Q(J)=C(I)*O                                                       DPRQ1100
      C(I)=Q(J)                                                         DPRQ1110
    9 CONTINUE                                                          DPRQ1120
C                                                                       DPRQ1130
C           INITIALIZATION                                              DPRQ1140
      ESAV=0.D0                                                         DPRQ1150
      Q(ISTA)=0.D0                                                      DPRQ1160
   10 NSAV=IR                                                           DPRQ1170
C                                                                       DPRQ1180
C        COMPUTATION OF DERIVATIVE                                      DPRQ1190
      EXPT=IR-ISTA                                                      DPRQ1200
      E(ISTA)=EXPT                                                      DPRQ1210
      DO 11 I=ISTA,IEND                                                 DPRQ1220
      EXPT=EXPT-1.0D0                                                   DPRQ1230
      POL(I+1)=EPS*DABS(Q(I+1))+EPS                                     DPRQ1240
   11 E(I+1)=Q(I+1)*EXPT                                                DPRQ1250
C                                                                       DPRQ1260
C        TEST OF REMAINING DIMENSION                                    DPRQ1270
      IF(ISTA-IEND)12,20,60                                             DPRQ1280
   12 JEND=IEND-1                                                       DPRQ1290
C                                                                       DPRQ1300
C        COMPUTATION OF S-FRACTION                                      DPRQ1310
      DO 19 I=ISTA,JEND                                                 DPRQ1320
      IF(I-ISTA)13,16,13                                                DPRQ1330
   13 IF(DABS(E(I))-POL(I+1))14,14,16                                   DPRQ1340
C                                                                       DPRQ1350
C        THE GIVEN POLYNOMIAL HAS MULTIPLE ROOTS, THE COEFFICIENTS OF   DPRQ1360
C        THE COMMON FACTOR ARE STORED FROM Q(NSAV) UP TO Q(IR)          DPRQ1370
   14 NSAV=I                                                            DPRQ1380
      DO 15 K=I,JEND                                                    DPRQ1390
      IF(DABS(E(K))-POL(K+1))15,15,80                                   DPRQ1400
   15 CONTINUE                                                          DPRQ1410
      GOTO 21                                                           DPRQ1420
C                                                                       DPRQ1430
C           EUCLIDEAN ALGORITHM                                         DPRQ1440
   16 DO 19 K=I,IEND                                                    DPRQ1450
      E(K+1)=E(K+1)/E(I)                                                DPRQ1460
      Q(K+1)=E(K+1)-Q(K+1)                                              DPRQ1470
      IF(K-I)18,17,18                                                   DPRQ1480
C                                                                       DPRQ1490
C        TEST FOR SMALL DIVISOR                                         DPRQ1500
   17 IF(DABS(Q(I+1))-POL(I+1))80,80,19                                 DPRQ1510
   18 Q(K+1)=Q(K+1)/Q(I+1)                                              DPRQ1520
      POL(K+1)=POL(K+1)/DABS(Q(I+1))                                    DPRQ1530
      E(K)=Q(K+1)-E(K)                                                  DPRQ1540
   19 CONTINUE                                                          DPRQ1550
   20 Q(IR)=-Q(IR)                                                      DPRQ1560
C                                                                       DPRQ1570
C           THE DISPLACEMENT EXPT IS SET TO 0 AUTOMATICALLY.            DPRQ1580
C           E(ISTA)=0.,Q(ISTA+1),...,E(NSAV-1),Q(NSAV),E(NSAV)=0.,      DPRQ1590
C           FORM A DIAGONAL OF THE QD-ARRAY.                            DPRQ1600
C        INITIALIZATION OF BOUNDARY VALUES                              DPRQ1610
   21 E(ISTA)=0.D0                                                      DPRQ1620
      NRAN=NSAV-1                                                       DPRQ1630
   22 E(NRAN+1)=0.D0                                                    DPRQ1640
C                                                                       DPRQ1650
C           TEST FOR LINEAR OR CONSTANT FACTOR                          DPRQ1660
C        NRAN-ISTA IS DEGREE-1                                          DPRQ1670
      IF(NRAN-ISTA)24,23,31                                             DPRQ1680
C                                                                       DPRQ1690
C        LINEAR FACTOR                                                  DPRQ1700
   23 Q(ISTA+1)=Q(ISTA+1)+EXPT                                          DPRQ1710
      E(ISTA+1)=0.D0                                                    DPRQ1720
C                                                                       DPRQ1730
C        TEST FOR UNFACTORED COMMON DIVISOR                             DPRQ1740
   24 E(ISTA)=ESAV                                                      DPRQ1750
      IF(IR-NSAV)60,60,25                                               DPRQ1760
C                                                                       DPRQ1770
C        INITIALIZE QD-ALGORITHM FOR COMMON DIVISOR                     DPRQ1780
   25 ISTA=NSAV                                                         DPRQ1790
      ESAV=E(ISTA)                                                      DPRQ1800
      GOTO 10                                                           DPRQ1810
C                                                                       DPRQ1820
C        COMPUTATION OF ROOT PAIR                                       DPRQ1830
   26 P=P+EXPT                                                          DPRQ1840
C                                                                       DPRQ1850
C        TEST FOR REALITY                                               DPRQ1860
      IF(O)27,28,28                                                     DPRQ1870
C                                                                       DPRQ1880
C        COMPLEX ROOT PAIR                                              DPRQ1890
   27 Q(NRAN)=P                                                         DPRQ1900
      Q(NRAN+1)=P                                                       DPRQ1910
      E(NRAN)=T                                                         DPRQ1920
      E(NRAN+1)=-T                                                      DPRQ1930
      GOTO 29                                                           DPRQ1940
C                                                                       DPRQ1950
C        REAL ROOT PAIR                                                 DPRQ1960
   28 Q(NRAN)=P-T                                                       DPRQ1970
      Q(NRAN+1)=P+T                                                     DPRQ1980
      E(NRAN)=0.D0                                                      DPRQ1990
C                                                                       DPRQ2000
C           REDUCTION OF DEGREE BY 2 (DEFLATION)                        DPRQ2010
   29 NRAN=NRAN-2                                                       DPRQ2020
      GOTO 22                                                           DPRQ2030
C                                                                       DPRQ2040
C        COMPUTATION OF REAL ROOT                                       DPRQ2050
   30 Q(NRAN+1)=EXPT+P                                                  DPRQ2060
C                                                                       DPRQ2070
C           REDUCTION OF DEGREE BY 1 (DEFLATION)                        DPRQ2080
      NRAN=NRAN-1                                                       DPRQ2090
      GOTO 22                                                           DPRQ2100
C                                                                       DPRQ2110
C        START QD-ITERATION                                             DPRQ2120
   31 JBEG=ISTA+1                                                       DPRQ2130
      JEND=NRAN-1                                                       DPRQ2140
      TEPS=EPS                                                          DPRQ2150
      TDELT=1.E-2                                                       DPRQ2160
   32 KOUNT=KOUNT+1                                                     DPRQ2170
      P=Q(NRAN+1)                                                       DPRQ2180
      R=ABS(SNGL(E(NRAN)))                                              DPRQ2190
C                                                                       DPRQ2200
C           TEST FOR CONVERGENCE                                        DPRQ2210
      IF(R-TEPS)30,30,33                                                DPRQ2220
   33 S=ABS(SNGL(E(JEND)))                                              DPRQ2230
C                                                                       DPRQ2240
C        IS THERE A REAL ROOT NEXT                                      DPRQ2250
      IF(S-R)38,38,34                                                   DPRQ2260
C                                                                       DPRQ2270
C        IS DISPLACEMENT SMALL ENOUGH                                   DPRQ2280
   34 IF(R-TDELT)36,35,35                                               DPRQ2290
   35 P=0.D0                                                            DPRQ2300
   36 O=P                                                               DPRQ2310
      DO 37 J=JBEG,NRAN                                                 DPRQ2320
      Q(J)=Q(J)+E(J)-E(J-1)-O                                           DPRQ2330
C                                                                       DPRQ2340
C           TEST FOR SMALL DIVISOR                                      DPRQ2350
      IF(DABS(Q(J))-POL(J))81,81,37                                     DPRQ2360
   37 E(J)=Q(J+1)*E(J)/Q(J)                                             DPRQ2370
      Q(NRAN+1)=-E(NRAN)+Q(NRAN+1)-O                                    DPRQ2380
      GOTO 54                                                           DPRQ2390
C                                                                       DPRQ2400
C        CALCULATE DISPLACEMENT FOR DOUBLE ROOTS                        DPRQ2410
C           QUADRATIC EQUATION FOR DOUBLE ROOTS                         DPRQ2420
C           X**2-(Q(NRAN)+Q(NRAN+1)+E(NRAN))*X+Q(NRAN)*Q(NRAN+1)=0      DPRQ2430
   38 P=0.5D0*(Q(NRAN)+E(NRAN)+Q(NRAN+1))                               DPRQ2440
      O=P*P-Q(NRAN)*Q(NRAN+1)                                           DPRQ2450
      T=DSQRT(DABS(O))                                                  DPRQ2460
C                                                                       DPRQ2470
C        TEST FOR CONVERGENCE                                           DPRQ2480
      IF(S-TEPS)26,26,39                                                DPRQ2490
C                                                                       DPRQ2500
C        ARE THERE COMPLEX ROOTS                                        DPRQ2510
   39 IF(O)43,40,40                                                     DPRQ2520
   40 IF(P)42,41,41                                                     DPRQ2530
   41 T=-T                                                              DPRQ2540
   42 P=P+T                                                             DPRQ2550
      R=S                                                               DPRQ2560
      GOTO 34                                                           DPRQ2570
C                                                                       DPRQ2580
C        MODIFICATION FOR COMPLEX ROOTS                                 DPRQ2590
C        IS DISPLACEMENT SMALL ENOUGH                                   DPRQ2600
   43 IF(S-TDELT)44,35,35                                               DPRQ2610
C                                                                       DPRQ2620
C           INITIALIZATION                                              DPRQ2630
   44 O=Q(JBEG)+E(JBEG)-P                                               DPRQ2640
C                                                                       DPRQ2650
C           TEST FOR SMALL DIVISOR                                      DPRQ2660
      IF(DABS(O)-POL(JBEG))81,81,45                                     DPRQ2670
   45 T=(T/O)**2                                                        DPRQ2680
      U=E(JBEG)*Q(JBEG+1)/(O*(1.0D0+T))                                 DPRQ2690
      V=O+U                                                             DPRQ2700
C                                                                       DPRQ2710
C        THREEFOLD LOOP FOR COMPLEX DISPLACEMENT                        DPRQ2720
      KOUNT=KOUNT+2                                                     DPRQ2730
      DO 53 J=JBEG,NRAN                                                 DPRQ2740
      O=Q(J+1)+E(J+1)-U-P                                               DPRQ2750
C                                                                       DPRQ2760
C           TEST FOR SMALL DIVISOR                                      DPRQ2770
      IF(DABS(V)-POL(J))46,46,49                                        DPRQ2780
   46 IF(J-NRAN)81,47,81                                                DPRQ2790
   47 EXPT=EXPT+P                                                       DPRQ2800
      IF(ABS(SNGL(E(JEND)))-TOL)48,48,81                                DPRQ2810
   48 P=0.5D0*(V+O-E(JEND))                                             DPRQ2820
      O=P*P-(V-U)*(O-U*T-O*W*(1.D0+T)/Q(JEND))                          DPRQ2830
      T=DSQRT(DABS(O))                                                  DPRQ2840
      GOTO 26                                                           DPRQ2850
C                                                                       DPRQ2860
C           TEST FOR SMALL DIVISOR                                      DPRQ2870
   49 IF(DABS(O)-POL(J+1))46,46,50                                      DPRQ2880
   50 W=U*O/V                                                           DPRQ2890
      T=T*(V/O)**2                                                      DPRQ2900
      Q(J)=V+W-E(J-1)                                                   DPRQ2910
      U=0.D0                                                            DPRQ2920
      IF(J-NRAN)51,52,52                                                DPRQ2930
   51 U=Q(J+2)*E(J+1)/(O*(1.D0+T))                                      DPRQ2940
   52 V=O+U-W                                                           DPRQ2950
C                                                                       DPRQ2960
C           TEST FOR SMALL DIVISOR                                      DPRQ2970
      IF(DABS(Q(J))-POL(J))81,81,53                                     DPRQ2980
   53 E(J)=W*V*(1.0D0+T)/Q(J)                                           DPRQ2990
      Q(NRAN+1)=V-E(NRAN)                                               DPRQ3000
   54 EXPT=EXPT+P                                                       DPRQ3010
      TEPS=TEPS*1.1                                                     DPRQ3020
      TDELT=TDELT*1.1                                                   DPRQ3030
      IF(KOUNT-LIMIT)32,55,55                                           DPRQ3040
C                                                                       DPRQ3050
C        NO CONVERGENCE WITH FEASIBLE TOLERANCE                         DPRQ3060
C           ERROR RETURN IN CASE OF UNSATISFACTORY CONVERGENCE          DPRQ3070
   55 IER=1                                                             DPRQ3080
C                                                                       DPRQ3090
C        REARRANGE CALCULATED ROOTS                                     DPRQ3100
   56 IEND=NSAV-NRAN-1                                                  DPRQ3110
      E(ISTA)=ESAV                                                      DPRQ3120
      IF(IEND)59,59,57                                                  DPRQ3130
   57 DO 58 I=1,IEND                                                    DPRQ3140
      J=ISTA+I                                                          DPRQ3150
      K=NRAN+1+I                                                        DPRQ3160
      E(J)=E(K)                                                         DPRQ3170
   58 Q(J)=Q(K)                                                         DPRQ3180
   59 IR=ISTA+IEND                                                      DPRQ3190
C                                                                       DPRQ3200
C        NORMAL RETURN                                                  DPRQ3210
   60 IR=IR-1                                                           DPRQ3220
      IF(IR)78,78,61                                                    DPRQ3230
C                                                                       DPRQ3240
C        REARRANGE CALCULATED ROOTS                                     DPRQ3250
   61 DO 62 I=1,IR                                                      DPRQ3260
      Q(I)=Q(I+1)                                                       DPRQ3270
   62 E(I)=E(I+1)                                                       DPRQ3280
C                                                                       DPRQ3290
C        CALCULATE COEFFICIENT VECTOR FROM ROOTS                        DPRQ3300
      POL(IR+1)=1.D0                                                    DPRQ3310
      IEND=IR-1                                                         DPRQ3320
      JBEG=1                                                            DPRQ3330
      DO 69 J=1,IR                                                      DPRQ3340
      ISTA=IR+1-J                                                       DPRQ3350
      O=0.D0                                                            DPRQ3360
      P=Q(ISTA)                                                         DPRQ3370
      T=E(ISTA)                                                         DPRQ3380
      IF(T)65,63,65                                                     DPRQ3390
C                                                                       DPRQ3400
C        MULTIPLY WITH LINEAR FACTOR                                    DPRQ3410
   63 DO 64 I=ISTA,IR                                                   DPRQ3420
      POL(I)=O-P*POL(I+1)                                               DPRQ3430
   64 O=POL(I+1)                                                        DPRQ3440
      GOTO 69                                                           DPRQ3450
   65 GOTO(66,67),JBEG                                                  DPRQ3460
   66 JBEG=2                                                            DPRQ3470
      POL(ISTA)=0.D0                                                    DPRQ3480
      GOTO 69                                                           DPRQ3490
C                                                                       DPRQ3500
C        MULTIPLY WITH QUADRATIC FACTOR                                 DPRQ3510
   67 JBEG=1                                                            DPRQ3520
      U=P*P+T*T                                                         DPRQ3530
      P=P+P                                                             DPRQ3540
      DO 68 I=ISTA,IEND                                                 DPRQ3550
      POL(I)=O-P*POL(I+1)+U*POL(I+2)                                    DPRQ3560
   68 O=POL(I+1)                                                        DPRQ3570
      POL(IR)=O-P                                                       DPRQ3580
   69 CONTINUE                                                          DPRQ3590
      IF(IER)78,70,78                                                   DPRQ3600
C                                                                       DPRQ3610
C        COMPARISON OF COEFFICIENT VECTORS, IE. TEST OF ACCURACY        DPRQ3620
   70 P=0.D0                                                            DPRQ3630
      DO 75 I=1,IR                                                      DPRQ3640
      IF(C(I))72,71,72                                                  DPRQ3650
   71 O=DABS(POL(I))                                                    DPRQ3660
      GOTO 73                                                           DPRQ3670
   72 O=DABS((POL(I)-C(I))/C(I))                                        DPRQ3680
   73 IF(P-O)74,75,75                                                   DPRQ3690
   74 P=O                                                               DPRQ3700
   75 CONTINUE                                                          DPRQ3710
      IF(SNGL(P)-TOL)77,76,76                                           DPRQ3720
   76 IER=-1                                                            DPRQ3730
   77 Q(IR+1)=P                                                         DPRQ3740
      E(IR+1)=0.D0                                                      DPRQ3750
   78 RETURN                                                            DPRQ3760
C                                                                       DPRQ3770
C        ERROR RETURNS                                                  DPRQ3780
C           ERROR RETURN FOR POLYNOMIALS OF DEGREE LESS THAN 1          DPRQ3790
   79 IER=2                                                             DPRQ3800
      IR=0                                                              DPRQ3810
      RETURN                                                            DPRQ3820
C                                                                       DPRQ3830
C           ERROR RETURN IF THERE EXISTS NO S-FRACTION                  DPRQ3840
   80 IER=4                                                             DPRQ3850
      IR=ISTA                                                           DPRQ3860
      GOTO 60                                                           DPRQ3870
C                                                                       DPRQ3880
C           ERROR RETURN IN CASE OF INSTABLE QD-ALGORITHM               DPRQ3890
   81 IER=3                                                             DPRQ3900
      GOTO 56                                                           DPRQ3910
      END                                                               DPRQ3920
