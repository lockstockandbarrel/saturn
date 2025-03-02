C                                                                       MFSS  10
C     ..................................................................MFSS  20
C                                                                       MFSS  30
C        SUBROUTINE MFSS                                                MFSS  40
C                                                                       MFSS  50
C        PURPOSE                                                        MFSS  60
C           GIVEN A SYMMETRIC POSITIVE SEMI DEFINITE MATRIX , MFSS WILL MFSS  70
C           (1) DETERMINE THE RANK AND LINEARLY INDEPENDENT ROWS AND    MFSS  80
C               COLUMNS                                                 MFSS  90
C           (2) FACTOR A SYMMETRIC SUBMATRIX OF MAXIMAL RANK            MFSS 100
C           (3) EXPRESS NONBASIC ROWS IN TERMS OF BASIC ONES,           MFSS 110
C               EXPRESS NONBASIC COLUMNS IN TERMS OF BASIC ONES         MFSS 120
C               EXPRESS BASIC VARIABLES IN TERMS OF FREE ONES           MFSS 130
C           SUBROUTINE MFSS MAY BE USED AS A PREPARATORY STEP FOR THE   MFSS 140
C           CALCULATION OF THE LEAST SQUARES SOLUTION OF MINIMAL        MFSS 150
C           LENGTH OF A SYSTEM OF LINEAR EQUATIONS WITH SYMMETRIC       MFSS 160
C           POSITIVE SEMI-DEFINITE COEFFICIENT MATRIX                   MFSS 170
C                                                                       MFSS 180
C        USAGE                                                          MFSS 190
C           CALL MFSS(A,N,EPS,IRANK,TRAC)                               MFSS 200
C                                                                       MFSS 210
C        DESCRIPTION OF PARAMETERS                                      MFSS 220
C           A     - UPPER TRIANGULAR PART OF GIVEN SYMMETRIC SEMI-      MFSS 230
C                   DEFINITE MATRIX STORED COLUMNWISE IN COMPRESSED FORMMFSS 240
C                   ON RETURN A CONTAINS THE MATRIX T AND, IF IRANK IS  MFSS 250
C                   LESS THAN N, THE MATRICES U AND TU                  MFSS 260
C           N     - DIMENSION OF GIVEN MATRIX A                         MFSS 270
C           EPS   - TESTVALUE FOR ZERO AFFECTED BY ROUND-OFF NOISE      MFSS 280
C           IRANK - RESULTANT VARIABLE, CONTAINING THE RANK OF GIVEN    MFSS 290
C                   MATRIX A IF A IS SEMI-DEFINITE                      MFSS 300
C                   IRANK = 0 MEANS A HAS NO POSITIVE DIAGONAL ELEMENT  MFSS 310
C                             AND/OR EPS IS NOT ABSOLUTELY LESS THAN ONEMFSS 320
C                   IRANK =-1 MEANS DIMENSION N IS NOT POSITIVE         MFSS 330
C                   IRANK =-2 MEANS COMPLETE FAILURE, POSSIBLY DUE TO   MFSS 340
C                             INADEQUATE RELATIVE TOLERANCE EPS         MFSS 350
C           TRAC  - VECTOR OF DIMENSION N CONTAINING THE                MFSS 360
C                   SOURCE INDEX OF THE I-TH PIVOT ROW IN ITS I-TH      MFSS 370
C                   LOCATION, THIS MEANS THAT TRAC CONTAINS THE         MFSS 380
C                   PRODUCT REPRESENTATION OF THE PERMUTATION WHICH     MFSS 390
C                   IS APPLIED TO ROWS AND COLUMNS OF A IN TERMS OF     MFSS 400
C                   TRANSPOSITIONS                                      MFSS 410
C                                                                       MFSS 420
C        REMARKS                                                        MFSS 430
C           EPS MUST BE ABSOLUTELY LESS THAN ONE. A SENSIBLE VALUE IS   MFSS 440
C           SOMEWHERE IN BETWEEN 10**(-4) AND 10**(-6)                  MFSS 450
C           THE ABSOLUTE VALUE OF INPUT PARAMETER EPS IS USED AS        MFSS 460
C           RELATIVE TOLERANCE.                                         MFSS 470
C           IN ORDER TO PRESERVE SYMMETRY ONLY PIVOTING ALONG THE       MFSS 480
C           DIAGONAL IS BUILT IN.                                       MFSS 490
C           ALL PIVOTELEMENTS MUST BE GREATER THAN THE ABSOLUTE VALUE   MFSS 500
C           OF EPS TIMES ORIGINAL DIAGONAL ELEMENT                      MFSS 510
C           OTHERWISE THEY ARE TREATED AS IF THEY WERE ZERO             MFSS 520
C           MATRIX A REMAINS UNCHANGED IF THE RESULTANT VALUE IRANK     MFSS 530
C           EQUALS ZERO                                                 MFSS 540
C                                                                       MFSS 550
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  MFSS 560
C           NONE                                                        MFSS 570
C                                                                       MFSS 580
C        METHOD                                                         MFSS 590
C           THE SQUARE ROOT METHOD WITH DIAGONAL PIVOTING IS USED FOR   MFSS 600
C           CALCULATION OF THE RIGHT HAND TRIANGULAR FACTOR.            MFSS 610
C           IN CASE OF AN ONLY SEMI-DEFINITE MATRIX THE SUBROUTINE      MFSS 620
C           RETURNS THE IRANK X IRANK UPPER TRIANGULAR FACTOR T OF A    MFSS 630
C           SUBMATRIX OF MAXIMAL RANK, THE IRANK X (N-IRANK) MATRIX U   MFSS 640
C           AND THE (N-IRANK) X (N-IRANK) UPPER TRIANGULAR TU SUCH      MFSS 650
C           THAT TRANSPOSE(TU)*TU=I+TRANSPOSE(U)*U                      MFSS 660
C                                                                       MFSS 670
C     ..................................................................MFSS 680
C                                                                       MFSS 690
      SUBROUTINE MFSS(A,N,EPS,IRANK,TRAC)                               MFSS 700
C                                                                       MFSS 710
C                                                                       MFSS 720
C        DIMENSIONED DUMMY VARIABLES                                    MFSS 730
      DIMENSION A(1),TRAC(1)                                            MFSS 740
      DOUBLE PRECISION SUM                                              MFSS 750
C                                                                       MFSS 760
C        TEST OF SPECIFIED DIMENSION                                    MFSS 770
      IF(N)36,36,1                                                      MFSS 780
C                                                                       MFSS 790
C        INITIALIZE TRIANGULAR FACTORIZATION                            MFSS 800
    1 IRANK=0                                                           MFSS 810
      ISUB=0                                                            MFSS 820
      KPIV=0                                                            MFSS 830
      J=0                                                               MFSS 840
      PIV=0.                                                            MFSS 850
C                                                                       MFSS 860
C        SEARCH FIRST PIVOT ELEMENT                                     MFSS 870
      DO 3 K=1,N                                                        MFSS 880
      J=J+K                                                             MFSS 890
      TRAC(K)=A(J)                                                      MFSS 900
      IF(A(J)-PIV)3,3,2                                                 MFSS 910
    2 PIV=A(J)                                                          MFSS 920
      KSUB=J                                                            MFSS 930
      KPIV=K                                                            MFSS 940
    3 CONTINUE                                                          MFSS 950
C                                                                       MFSS 960
C        START LOOP OVER ALL ROWS OF A                                  MFSS 970
      DO 32 I=1,N                                                       MFSS 980
      ISUB=ISUB+I                                                       MFSS 990
      IM1=I-1                                                           MFSS1000
    4 KMI=KPIV-I                                                        MFSS1010
      IF(KMI)35,9,5                                                     MFSS1020
C                                                                       MFSS1030
C        PERFORM PARTIAL COLUMN INTERCHANGE                             MFSS1040
    5 JI=KSUB-KMI                                                       MFSS1050
      IDC=JI-ISUB                                                       MFSS1060
      JJ=ISUB-IM1                                                       MFSS1070
      DO 6 K=JJ,ISUB                                                    MFSS1080
      KK=K+IDC                                                          MFSS1090
      HOLD=A(K)                                                         MFSS1100
      A(K)=A(KK)                                                        MFSS1110
    6 A(KK)=HOLD                                                        MFSS1120
C                                                                       MFSS1130
C        PERFORM PARTIAL ROW INTERCHANGE                                MFSS1140
      KK=KSUB                                                           MFSS1150
      DO 7 K=KPIV,N                                                     MFSS1160
      II=KK-KMI                                                         MFSS1170
      HOLD=A(KK)                                                        MFSS1180
      A(KK)=A(II)                                                       MFSS1190
      A(II)=HOLD                                                        MFSS1200
    7 KK=KK+K                                                           MFSS1210
C                                                                       MFSS1220
C        PERFORM REMAINING INTERCHANGE                                  MFSS1230
      JJ=KPIV-1                                                         MFSS1240
      II=ISUB                                                           MFSS1250
      DO 8 K=I,JJ                                                       MFSS1260
      HOLD=A(II)                                                        MFSS1270
      A(II)=A(JI)                                                       MFSS1280
      A(JI)=HOLD                                                        MFSS1290
      II=II+K                                                           MFSS1300
    8 JI=JI+1                                                           MFSS1310
    9 IF(IRANK)22,10,10                                                 MFSS1320
C                                                                       MFSS1330
C        RECORD INTERCHANGE IN TRANSPOSITION VECTOR                     MFSS1340
   10 TRAC(KPIV)=TRAC(I)                                                MFSS1350
      TRAC(I)=KPIV                                                      MFSS1360
C                                                                       MFSS1370
C        MODIFY CURRENT PIVOT ROW                                       MFSS1380
      KK=IM1-IRANK                                                      MFSS1390
      KMI=ISUB-KK                                                       MFSS1400
      PIV=0.                                                            MFSS1410
      IDC=IRANK+1                                                       MFSS1420
      JI=ISUB-1                                                         MFSS1430
      JK=KMI                                                            MFSS1440
      JJ=ISUB-I                                                         MFSS1450
      DO 19 K=I,N                                                       MFSS1460
      SUM=0.D0                                                          MFSS1470
C                                                                       MFSS1480
C        BUILD UP SCALAR PRODUCT IF NECESSARY                           MFSS1490
      IF(KK)13,13,11                                                    MFSS1500
   11 DO 12 J=KMI,JI                                                    MFSS1510
      SUM=SUM-A(J)*A(JK)                                                MFSS1520
   12 JK=JK+1                                                           MFSS1530
   13 JJ=JJ+K                                                           MFSS1540
      IF(K-I)14,14,16                                                   MFSS1550
   14 SUM=A(ISUB)+SUM                                                   MFSS1560
C                                                                       MFSS1570
C        TEST RADICAND FOR LOSS OF SIGNIFICANCE                         MFSS1580
      IF(SUM-ABS(A(ISUB)*EPS))20,20,15                                  MFSS1590
   15 A(ISUB)=DSQRT(SUM)                                                MFSS1600
      KPIV=I+1                                                          MFSS1610
      GOTO 19                                                           MFSS1620
   16 SUM=(A(JK)+SUM)/A(ISUB)                                           MFSS1630
      A(JK)=SUM                                                         MFSS1640
C                                                                       MFSS1650
C        SEARCH FOR NEXT PIVOT ROW                                      MFSS1660
      IF(A(JJ))19,19,17                                                 MFSS1670
   17 TRAC(K)=TRAC(K)-SUM*SUM                                           MFSS1680
      HOLD=TRAC(K)/A(JJ)                                                MFSS1690
      IF(PIV-HOLD)18,19,19                                              MFSS1700
   18 PIV=HOLD                                                          MFSS1710
      KPIV=K                                                            MFSS1720
      KSUB=JJ                                                           MFSS1730
   19 JK=JJ+IDC                                                         MFSS1740
      GOTO 32                                                           MFSS1750
C                                                                       MFSS1760
C        CALCULATE MATRIX OF DEPENDENCIES U                             MFSS1770
   20 IF(IRANK)21,21,37                                                 MFSS1780
   21 IRANK=-1                                                          MFSS1790
      GOTO 4                                                            MFSS1800
   22 IRANK=IM1                                                         MFSS1810
      II=ISUB-IRANK                                                     MFSS1820
      JI=II                                                             MFSS1830
      DO 26 K=1,IRANK                                                   MFSS1840
      JI=JI-1                                                           MFSS1850
      JK=ISUB-1                                                         MFSS1860
      JJ=K-1                                                            MFSS1870
      DO 26 J=I,N                                                       MFSS1880
      IDC=IRANK                                                         MFSS1890
      SUM=0.D0                                                          MFSS1900
      KMI=JI                                                            MFSS1910
      KK=JK                                                             MFSS1920
      IF(JJ)25,25,23                                                    MFSS1930
   23 DO 24 L=1,JJ                                                      MFSS1940
      IDC=IDC-1                                                         MFSS1950
      SUM=SUM-A(KMI)*A(KK)                                              MFSS1960
      KMI=KMI-IDC                                                       MFSS1970
   24 KK=KK-1                                                           MFSS1980
   25 A(KK)=(SUM+A(KK))/A(KMI)                                          MFSS1990
   26 JK=JK+J                                                           MFSS2000
C                                                                       MFSS2010
C        CALCULATE I+TRANSPOSE(U)*U                                     MFSS2020
      JJ=ISUB-I                                                         MFSS2030
      PIV=0.                                                            MFSS2040
      KK=ISUB-1                                                         MFSS2050
      DO 31 K=I,N                                                       MFSS2060
      JJ=JJ+K                                                           MFSS2070
      IDC=0                                                             MFSS2080
      DO 28 J=K,N                                                       MFSS2090
      SUM=0.D0                                                          MFSS2100
      KMI=JJ+IDC                                                        MFSS2110
      DO 27 L=II,KK                                                     MFSS2120
      JK=L+IDC                                                          MFSS2130
   27 SUM=SUM+A(L)*A(JK)                                                MFSS2140
      A(KMI)=SUM                                                        MFSS2150
   28 IDC=IDC+J                                                         MFSS2160
      A(JJ)=A(JJ)+1.D0                                                  MFSS2170
      TRAC(K)=A(JJ)                                                     MFSS2180
C                                                                       MFSS2190
C        SEARCH NEXT DIAGONAL ELEMENT                                   MFSS2200
      IF(PIV-A(JJ))29,30,30                                             MFSS2210
   29 KPIV=K                                                            MFSS2220
      KSUB=JJ                                                           MFSS2230
      PIV=A(JJ)                                                         MFSS2240
   30 II=II+K                                                           MFSS2250
      KK=KK+K                                                           MFSS2260
   31 CONTINUE                                                          MFSS2270
      GOTO 4                                                            MFSS2280
   32 CONTINUE                                                          MFSS2290
   33 IF(IRANK)35,34,35                                                 MFSS2300
   34 IRANK=N                                                           MFSS2310
   35 RETURN                                                            MFSS2320
C                                                                       MFSS2330
C        ERROR RETURNS                                                  MFSS2340
C                                                                       MFSS2350
C        RETURN IN CASE OF ILLEGAL DIMENSION                            MFSS2360
   36 IRANK=-1                                                          MFSS2370
      RETURN                                                            MFSS2380
C                                                                       MFSS2390
C        INSTABLE FACTORIZATION OF I+TRANSPOSE(U)*U                     MFSS2400
   37 IRANK=-2                                                          MFSS2410
      RETURN                                                            MFSS2420
      END                                                               MFSS2430
