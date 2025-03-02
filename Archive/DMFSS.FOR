C                                                                       DMSS  10
C     ..................................................................DMSS  20
C                                                                       DMSS  30
C        SUBROUTINE DMFSS                                               DMSS  40
C                                                                       DMSS  50
C        PURPOSE                                                        DMSS  60
C           GIVEN A SYMMETRIC POSITIVE SEMI DEFINITE MATRIX ,DMFSS WILL DMSS  70
C           (1) DETERMINE THE RANK AND LINEARLY INDEPENDENT ROWS AND    DMSS  80
C               COLUMNS                                                 DMSS  90
C           (2) FACTOR A SYMMETRIC SUBMATRIX OF MAXIMAL RANK            DMSS 100
C           (3) EXPRESS NONBASIC ROWS IN TERMS OF BASIC ONES,           DMSS 110
C               EXPRESS NONBASIC COLUMNS IN TERMS OF BASIC ONES         DMSS 120
C               EXPRESS BASIC VARIABLES IN TERMS OF FREE ONES           DMSS 130
C           SUBROUTINE DMFSS MAY BE USED AS A PREPARATORY STEP FOR THE  DMSS 140
C           CALCULATION OF THE LEAST SQUARES SOLUTION OF MINIMAL        DMSS 150
C           LENGTH OF A SYSTEM OF LINEAR EQUATIONS WITH SYMMETRIC       DMSS 160
C           POSITIVE SEMI-DEFINITE COEFFICIENT MATRIX                   DMSS 170
C                                                                       DMSS 180
C        USAGE                                                          DMSS 190
C           CALL DMFSS(A,N,EPS,IRANK,TRAC)                              DMSS 200
C                                                                       DMSS 210
C        DESCRIPTION OF PARAMETERS                                      DMSS 220
C           A     - UPPER TRIANGULAR PART OF GIVEN SYMMETRIC SEMI-      DMSS 230
C                   DEFINITE MATRIX STORED COLUMNWISE IN COMPRESSED FORMDMSS 240
C                   ON RETURN A CONTAINS THE MATRIX T AND, IF IRANK IS  DMSS 250
C                   LESS THAN N, THE MATRICES U AND TU                  DMSS 260
C                   A MUST BE OF DOUBLE PRECISION                       DMSS 270
C           N     - DIMENSION OF GIVEN MATRIX A                         DMSS 280
C           EPS   - TESTVALUE FOR ZERO AFFECTED BY ROUND-OFF NOISE      DMSS 290
C           IRANK - RESULTANT VARIABLE, CONTAINING THE RANK OF GIVEN    DMSS 300
C                   MATRIX A IF A IS SEMI-DEFINITE                      DMSS 310
C                   IRANK = 0 MEANS A HAS NO POSITIVE DIAGONAL ELEMENT  DMSS 320
C                             AND/OR EPS IS NOT ABSOLUTELY LESS THAN ONEDMSS 330
C                   IRANK =-1 MEANS DIMENSION N IS NOT POSITIVE         DMSS 340
C                   IRANK =-2 MEANS COMPLETE FAILURE, POSSIBLY DUE TO   DMSS 350
C                             INADEQUATE RELATIVE TOLERANCE EPS         DMSS 360
C           TRAC  - VECTOR OF DIMENSION N CONTAINING THE                DMSS 370
C                   SOURCE INDEX OF THE I-TH PIVOT ROW IN ITS I-TH      DMSS 380
C                   LOCATION, THIS MEANS THAT TRAC CONTAINS THE         DMSS 390
C                   PRODUCT REPRESENTATION OF THE PERMUTATION WHICH     DMSS 400
C                   IS APPLIED TO ROWS AND COLUMNS OF A IN TERMS OF     DMSS 410
C                   TRANSPOSITIONS                                      DMSS 420
C                   TRAC MUST BE OF DOUBLE PRECISION                    DMSS 430
C                                                                       DMSS 440
C        REMARKS                                                        DMSS 450
C           EPS MUST BE ABSOLUTELY LESS THAN ONE. A SENSIBLE VALUE IS   DMSS 460
C           SOMEWHERE IN BETWEEN 10**(-4) AND 10**(-6)                  DMSS 470
C           THE ABSOLUTE VALUE OF INPUT PARAMETER EPS IS USED AS        DMSS 480
C           RELATIVE TOLERANCE.                                         DMSS 490
C           IN ORDER TO PRESERVE SYMMETRY ONLY PIVOTING ALONG THE       DMSS 500
C           DIAGONAL IS BUILT IN.                                       DMSS 510
C           ALL PIVOTELEMENTS MUST BE GREATER THAN THE ABSOLUTE VALUE   DMSS 520
C           OF EPS TIMES ORIGINAL DIAGONAL ELEMENT                      DMSS 530
C           OTHERWISE THEY ARE TREATED AS IF THEY WERE ZERO             DMSS 540
C           MATRIX A REMAINS UNCHANGED IF THE RESULTANT VALUE IRANK     DMSS 550
C           EQUALS ZERO                                                 DMSS 560
C                                                                       DMSS 570
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DMSS 580
C           NONE                                                        DMSS 590
C                                                                       DMSS 600
C        METHOD                                                         DMSS 610
C           THE SQUARE ROOT METHOD WITH DIAGONAL PIVOTING IS USED FOR   DMSS 620
C           CALCULATION OF THE RIGHT HAND TRIANGULAR FACTOR.            DMSS 630
C           IN CASE OF AN ONLY SEMI-DEFINITE MATRIX THE SUBROUTINE      DMSS 640
C           RETURNS THE IRANK X IRANK UPPER TRIANGULAR FACTOR T OF A    DMSS 650
C           SUBMATRIX OF MAXIMAL RANK, THE IRANK X (N-IRANK) MATRIX U   DMSS 660
C           AND THE (N-IRANK) X (N-IRANK) UPPER TRIANGULAR TU SUCH      DMSS 670
C           THAT TRANSPOSE(TU)*TU=I+TRANSPOSE(U)*U                      DMSS 680
C                                                                       DMSS 690
C     ..................................................................DMSS 700
C                                                                       DMSS 710
      SUBROUTINE DMFSS(A,N,EPS,IRANK,TRAC)                              DMSS 720
C                                                                       DMSS 730
C                                                                       DMSS 740
C        DIMENSIONED DUMMY VARIABLES                                    DMSS 750
      DIMENSION A(1),TRAC(1)                                            DMSS 760
      DOUBLE PRECISION SUM,A,TRAC,PIV,HOLD                              DMSS 770
C                                                                       DMSS 780
C        TEST OF SPECIFIED DIMENSION                                    DMSS 790
      IF(N)36,36,1                                                      DMSS 800
C                                                                       DMSS 810
C        INITIALIZE TRIANGULAR FACTORIZATION                            DMSS 820
    1 IRANK=0                                                           DMSS 830
      ISUB=0                                                            DMSS 840
      KPIV=0                                                            DMSS 850
      J=0                                                               DMSS 860
      PIV=0.D0                                                          DMSS 870
C                                                                       DMSS 880
C        SEARCH FIRST PIVOT ELEMENT                                     DMSS 890
      DO 3 K=1,N                                                        DMSS 900
      J=J+K                                                             DMSS 910
      TRAC(K)=A(J)                                                      DMSS 920
      IF(A(J)-PIV)3,3,2                                                 DMSS 930
    2 PIV=A(J)                                                          DMSS 940
      KSUB=J                                                            DMSS 950
      KPIV=K                                                            DMSS 960
    3 CONTINUE                                                          DMSS 970
C                                                                       DMSS 980
C        START LOOP OVER ALL ROWS OF A                                  DMSS 990
      DO 32 I=1,N                                                       DMSS1000
      ISUB=ISUB+I                                                       DMSS1010
      IM1=I-1                                                           DMSS1020
    4 KMI=KPIV-I                                                        DMSS1030
      IF(KMI)35,9,5                                                     DMSS1040
C                                                                       DMSS1050
C        PERFORM PARTIAL COLUMN INTERCHANGE                             DMSS1060
    5 JI=KSUB-KMI                                                       DMSS1070
      IDC=JI-ISUB                                                       DMSS1080
      JJ=ISUB-IM1                                                       DMSS1090
      DO 6 K=JJ,ISUB                                                    DMSS1100
      KK=K+IDC                                                          DMSS1110
      HOLD=A(K)                                                         DMSS1120
      A(K)=A(KK)                                                        DMSS1130
    6 A(KK)=HOLD                                                        DMSS1140
C                                                                       DMSS1150
C        PERFORM PARTIAL ROW INTERCHANGE                                DMSS1160
      KK=KSUB                                                           DMSS1170
      DO 7 K=KPIV,N                                                     DMSS1180
      II=KK-KMI                                                         DMSS1190
      HOLD=A(KK)                                                        DMSS1200
      A(KK)=A(II)                                                       DMSS1210
      A(II)=HOLD                                                        DMSS1220
    7 KK=KK+K                                                           DMSS1230
C                                                                       DMSS1240
C        PERFORM REMAINING INTERCHANGE                                  DMSS1250
      JJ=KPIV-1                                                         DMSS1260
      II=ISUB                                                           DMSS1270
      DO 8 K=I,JJ                                                       DMSS1280
      HOLD=A(II)                                                        DMSS1290
      A(II)=A(JI)                                                       DMSS1300
      A(JI)=HOLD                                                        DMSS1310
      II=II+K                                                           DMSS1320
    8 JI=JI+1                                                           DMSS1330
    9 IF(IRANK)22,10,10                                                 DMSS1340
C                                                                       DMSS1350
C        RECORD INTERCHANGE IN TRANSPOSITION VECTOR                     DMSS1360
   10 TRAC(KPIV)=TRAC(I)                                                DMSS1370
      TRAC(I)=KPIV                                                      DMSS1380
C                                                                       DMSS1390
C        MODIFY CURRENT PIVOT ROW                                       DMSS1400
      KK=IM1-IRANK                                                      DMSS1410
      KMI=ISUB-KK                                                       DMSS1420
      PIV=0.D0                                                          DMSS1430
      IDC=IRANK+1                                                       DMSS1440
      JI=ISUB-1                                                         DMSS1450
      JK=KMI                                                            DMSS1460
      JJ=ISUB-I                                                         DMSS1470
      DO 19 K=I,N                                                       DMSS1480
      SUM=0.D0                                                          DMSS1490
C                                                                       DMSS1500
C        BUILD UP SCALAR PRODUCT IF NECESSARY                           DMSS1510
      IF(KK)13,13,11                                                    DMSS1520
   11 DO 12 J=KMI,JI                                                    DMSS1530
      SUM=SUM-A(J)*A(JK)                                                DMSS1540
   12 JK=JK+1                                                           DMSS1550
   13 JJ=JJ+K                                                           DMSS1560
      IF(K-I)14,14,16                                                   DMSS1570
   14 SUM=A(ISUB)+SUM                                                   DMSS1580
C                                                                       DMSS1590
C        TEST RADICAND FOR LOSS OF SIGNIFICANCE                         DMSS1600
      IF(SUM-DABS(A(ISUB)*DBLE(EPS)))20,20,15                           DMSS1610
   15 A(ISUB)=DSQRT(SUM)                                                DMSS1620
      KPIV=I+1                                                          DMSS1630
      GOTO 19                                                           DMSS1640
   16 SUM=(A(JK)+SUM)/A(ISUB)                                           DMSS1650
      A(JK)=SUM                                                         DMSS1660
C                                                                       DMSS1670
C        SEARCH FOR NEXT PIVOT ROW                                      DMSS1680
      IF(A(JJ))19,19,17                                                 DMSS1690
   17 TRAC(K)=TRAC(K)-SUM*SUM                                           DMSS1700
      HOLD=TRAC(K)/A(JJ)                                                DMSS1710
      IF(PIV-HOLD)18,19,19                                              DMSS1720
   18 PIV=HOLD                                                          DMSS1730
      KPIV=K                                                            DMSS1740
      KSUB=JJ                                                           DMSS1750
   19 JK=JJ+IDC                                                         DMSS1760
      GOTO 32                                                           DMSS1770
C                                                                       DMSS1780
C        CALCULATE MATRIX OF DEPENDENCIES U                             DMSS1790
   20 IF(IRANK)21,21,37                                                 DMSS1800
   21 IRANK=-1                                                          DMSS1810
      GOTO 4                                                            DMSS1820
   22 IRANK=IM1                                                         DMSS1830
      II=ISUB-IRANK                                                     DMSS1840
      JI=II                                                             DMSS1850
      DO 26 K=1,IRANK                                                   DMSS1860
      JI=JI-1                                                           DMSS1870
      JK=ISUB-1                                                         DMSS1880
      JJ=K-1                                                            DMSS1890
      DO 26 J=I,N                                                       DMSS1900
      IDC=IRANK                                                         DMSS1910
      SUM=0.D0                                                          DMSS1920
      KMI=JI                                                            DMSS1930
      KK=JK                                                             DMSS1940
      IF(JJ)25,25,23                                                    DMSS1950
   23 DO 24 L=1,JJ                                                      DMSS1960
      IDC=IDC-1                                                         DMSS1970
      SUM=SUM-A(KMI)*A(KK)                                              DMSS1980
      KMI=KMI-IDC                                                       DMSS1990
   24 KK=KK-1                                                           DMSS2000
   25 A(KK)=(SUM+A(KK))/A(KMI)                                          DMSS2010
   26 JK=JK+J                                                           DMSS2020
C                                                                       DMSS2030
C        CALCULATE I+TRANSPOSE(U)*U                                     DMSS2040
      JJ=ISUB-I                                                         DMSS2050
      PIV=0.D0                                                          DMSS2060
      KK=ISUB-1                                                         DMSS2070
      DO 31 K=I,N                                                       DMSS2080
      JJ=JJ+K                                                           DMSS2090
      IDC=0                                                             DMSS2100
      DO 28 J=K,N                                                       DMSS2110
      SUM=0.D0                                                          DMSS2120
      KMI=JJ+IDC                                                        DMSS2130
      DO 27 L=II,KK                                                     DMSS2140
      JK=L+IDC                                                          DMSS2150
   27 SUM=SUM+A(L)*A(JK)                                                DMSS2160
      A(KMI)=SUM                                                        DMSS2170
   28 IDC=IDC+J                                                         DMSS2180
      A(JJ)=A(JJ)+1.D0                                                  DMSS2190
      TRAC(K)=A(JJ)                                                     DMSS2200
C                                                                       DMSS2210
C        SEARCH NEXT DIAGONAL ELEMENT                                   DMSS2220
      IF(PIV-A(JJ))29,30,30                                             DMSS2230
   29 KPIV=K                                                            DMSS2240
      KSUB=JJ                                                           DMSS2250
      PIV=A(JJ)                                                         DMSS2260
   30 II=II+K                                                           DMSS2270
      KK=KK+K                                                           DMSS2280
   31 CONTINUE                                                          DMSS2290
      GOTO 4                                                            DMSS2300
   32 CONTINUE                                                          DMSS2310
   33 IF(IRANK)35,34,35                                                 DMSS2320
   34 IRANK=N                                                           DMSS2330
   35 RETURN                                                            DMSS2340
C                                                                       DMSS2350
C        ERROR RETURNS                                                  DMSS2360
C                                                                       DMSS2370
C        RETURN IN CASE OF ILLEGAL DIMENSION                            DMSS2380
   36 IRANK=-1                                                          DMSS2390
      RETURN                                                            DMSS2400
C                                                                       DMSS2410
C        INSTABLE FACTORIZATION OF I+TRANSPOSE(U)*U                     DMSS2420
   37 IRANK=-2                                                          DMSS2430
      RETURN                                                            DMSS2440
      END                                                               DMSS2450
