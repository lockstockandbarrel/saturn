C                                                                       MCHB  10
C     ..................................................................MCHB  20
C                                                                       MCHB  30
C        SUBROUTINE MCHB                                                MCHB  40
C                                                                       MCHB  50
C        PURPOSE                                                        MCHB  60
C           FOR A GIVEN POSITIVE-DEFINITE M BY M MATRIX A WITH SYMMETRICMCHB  70
C           BAND STRUCTURE AND - IF NECESSARY - A GIVEN GENERAL M BY N  MCHB  80
C           MATRIX R, THE FOLLOWING CALCULATIONS (DEPENDENT ON THE      MCHB  90
C           VALUE OF THE DECISION PARAMETER IOP) ARE PERFORMED          MCHB 100
C           (1) MATRIX A IS FACTORIZED (IF IOP IS NOT NEGATIVE), THAT   MCHB 110
C               MEANS BAND MATRIX TU WITH UPPER CODIAGONALS ONLY IS     MCHB 120
C               GENERATED ON THE LOCATIONS OF A SUCH THAT               MCHB 130
C               TRANSPOSE(TU)*TU=A.                                     MCHB 140
C           (2) MATRIX R IS MULTIPLIED ON THE LEFT BY INVERSE(TU)       MCHB 150
C               AND/OR INVERSE(TRANSPOSE(TU)) AND THE RESULT IS STORED  MCHB 160
C               IN THE LOCATIONS OF R.                                  MCHB 170
C           THIS SUBROUTINE ESPECIALLY CAN BE USED TO SOLVE THE SYSTEM  MCHB 180
C           OF SIMULTANEOUS LINEAR EQUATIONS A*X=R WITH POSITIVE-       MCHB 190
C           DEFINITE COEFFICIENT MATRIX A OF SYMMETRIC BAND STRUCTURE.  MCHB 200
C                                                                       MCHB 210
C        USAGE                                                          MCHB 220
C           CALL MCHB (R,A,M,N,MUD,IOP,EPS,IER)                         MCHB 230
C                                                                       MCHB 240
C        DESCRIPTION OF PARAMETERS                                      MCHB 250
C           R      - INPUT IN CASES IOP=-3,-2,-1,1,2,3  M BY N RIGHT    MCHB 260
C                          HAND SIDE MATRIX,                            MCHB 270
C                          IN CASE IOP=0  IRRELEVANT.                   MCHB 280
C                    OUTPUT IN CASES IOP=1,-1  INVERSE(A)*R,            MCHB 290
C                           IN CASES IOP=2,-2  INVERSE(TU)*R,           MCHB 300
C                           IN CASES IOP=3,-3  INVERSE(TRANSPOSE(TU))*R,MCHB 310
C                           IN CASE  IOP=0     UNCHANGED.               MCHB 320
C           A      - INPUT IN CASES IOP=0,1,2,3 M BY M POSITIVE-DEFINITEMCHB 330
C                          COEFFICIENT MATRIX OF SYMMETRIC BAND STRUC-  MCHB 340
C                          TURE STORED IN COMPRESSED FORM (SEE REMARKS),MCHB 350
C                          IN CASES IOP=-1,-2,-3  M BY M BAND MATRIX TU MCHB 360
C                          WITH UPPER CODIAGONALS ONLY, STORED IN       MCHB 370
C                          COMPRESSED FORM (SEE REMARKS).               MCHB 380
C                    OUTPUT IN ALL CASES  BAND MATRIX TU WITH UPPER     MCHB 390
C                           CODIAGONALS ONLY, STORED IN COMPRESSED FORM MCHB 400
C                           (THAT MEANS UNCHANGED IF IOP=-1,-2,-3).     MCHB 410
C           M      - INPUT VALUE SPECIFYING THE NUMBER OF ROWS AND      MCHB 420
C                    COLUMNS OF A AND THE NUMBER OF ROWS OF R.          MCHB 430
C           N      - INPUT VALUE SPECIFYING THE NUMBER OF COLUMNS OF R  MCHB 440
C                    (IRRELEVANT IN CASE IOP=0).                        MCHB 450
C           MUD    - INPUT VALUE SPECIFYING THE NUMBER OF UPPER         MCHB 460
C                    CODIAGONALS OF A.                                  MCHB 470
C           IOP    - ONE OF THE VALUES -3,-2,-1,0,1,2,3 GIVEN AS INPUT  MCHB 480
C                    AND USED AS DECISION PARAMETER.                    MCHB 490
C           EPS    - INPUT VALUE USED AS RELATIVE TOLERANCE FOR TEST ON MCHB 500
C                    LOSS OF SIGNIFICANT DIGITS.                        MCHB 510
C           IER    - RESULTING ERROR PARAMETER CODED AS FOLLOWS         MCHB 520
C                     IER=0  - NO ERROR,                                MCHB 530
C                     IER=-1 - NO RESULT BECAUSE OF WRONG INPUT         MCHB 540
C                              PARAMETERS M,MUD,IOP (SEE REMARKS),      MCHB 550
C                              OR BECAUSE OF A NONPOSITIVE RADICAND AT  MCHB 560
C                              SOME FACTORIZATION STEP,                 MCHB 570
C                              OR BECAUSE OF A ZERO DIAGONAL ELEMENT    MCHB 580
C                              AT SOME DIVISION STEP.                   MCHB 590
C                     IER=K  - WARNING DUE TO POSSIBLE LOSS OF SIGNIFI- MCHB 600
C                              CANCE INDICATED AT FACTORIZATION STEP K+1MCHB 610
C                              WHERE RADICAND WAS NO LONGER GREATER     MCHB 620
C                              THAN EPS*A(K+1,K+1).                     MCHB 630
C                                                                       MCHB 640
C        REMARKS                                                        MCHB 650
C           UPPER PART OF SYMMETRIC BAND MATRIX A CONSISTING OF MAIN    MCHB 660
C           DIAGONAL AND MUD UPPER CODIAGONALS (RESP. BAND MATRIX TU    MCHB 670
C           CONSISTING OF MAIN DIAGONAL AND MUD UPPER CODIAGONALS)      MCHB 680
C           IS ASSUMED TO BE STORED IN COMPRESSED FORM, I.E. ROWWISE    MCHB 690
C           IN TOTALLY NEEDED M+MUD*(2M-MUD-1)/2 SUCCESSIVE STORAGE     MCHB 700
C           LOCATIONS. ON RETURN UPPER BAND FACTOR TU (ON THE LOCATIONS MCHB 710
C           OF A) IS STORED IN THE SAME WAY.                            MCHB 720
C           RIGHT HAND SIDE MATRIX R IS ASSUMED TO BE STORED COLUMNWISE MCHB 730
C           IN N*M SUCCESSIVE STORAGE LOCATIONS. ON RETURN RESULT MATRIXMCHB 740
C           INVERSE(A)*R OR INVERSE(TU)*R OR INVERSE(TRANSPOSE(TU))*R   MCHB 750
C           IS STORED COLUMNWISE TOO ON THE LOCATIONS OF R.             MCHB 760
C           INPUT PARAMETERS M, MUD, IOP SHOULD SATISFY THE FOLLOWING   MCHB 770
C           RESTRICTIONS     MUD NOT LESS THAN ZERO,                    MCHB 780
C                            1+MUD NOT GREATER THAN M,                  MCHB 790
C                            ABS(IOP) NOT GREATER THAN 3.               MCHB 800
C           NO ACTION BESIDES ERROR MESSAGE IER=-1 TAKES PLACE IF THESE MCHB 810
C           RESTRICTIONS ARE NOT SATISFIED.                             MCHB 820
C           THE PROCEDURE GIVES RESULTS IF THE RESTRICTIONS ON INPUT    MCHB 830
C           PARAMETERS ARE SATISFIED, IF RADICANDS AT ALL FACTORIZATION MCHB 840
C           STEPS ARE POSITIVE AND/OR IF ALL DIAGONAL ELEMENTS OF       MCHB 850
C           UPPER BAND FACTOR TU ARE NONZERO.                           MCHB 860
C                                                                       MCHB 870
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  MCHB 880
C           NONE                                                        MCHB 890
C                                                                       MCHB 900
C        METHOD                                                         MCHB 910
C           FACTORIZATION IS DONE USING CHOLESKY-S SQUARE-ROOT METHOD,  MCHB 920
C           WHICH GENERATES THE UPPER BAND MATRIX TU SUCH THAT          MCHB 930
C           TRANSPOSE(TU)*TU=A. TU IS RETURNED AS RESULT ON THE         MCHB 940
C           LOCATIONS OF A. FURTHER, DEPENDENT ON THE ACTUAL VALUE OF   MCHB 950
C           IOP, DIVISION OF R BY TRANSPOSE(TU) AND/OR TU IS PERFORMED  MCHB 960
C           AND THE RESULT IS RETURNED ON THE LOCATIONS OF R.           MCHB 970
C           FOR REFERENCE, SEE H. RUTISHAUSER, ALGORITHMUS 1 - LINEARES MCHB 980
C           GLEICHUNGSSYSTEM MIT SYMMETRISCHER POSITIV-DEFINITER        MCHB 990
C           BANDMATRIX NACH CHOLESKY - , COMPUTING (ARCHIVES FOR        MCHB1000
C           ELECTRONIC COMPUTING), VOL.1, ISS.1 (1966), PP.77-78.       MCHB1010
C                                                                       MCHB1020
C     ..................................................................MCHB1030
C                                                                       MCHB1040
      SUBROUTINE MCHB(R,A,M,N,MUD,IOP,EPS,IER)                          MCHB1050
C                                                                       MCHB1060
C                                                                       MCHB1070
      DIMENSION R(1),A(1)                                               MCHB1080
      DOUBLE PRECISION TOL,SUM,PIV                                      MCHB1090
C                                                                       MCHB1100
C        TEST ON WRONG INPUT PARAMETERS                                 MCHB1110
      IF(IABS(IOP)-3)1,1,43                                             MCHB1120
    1 IF(MUD)43,2,2                                                     MCHB1130
    2 MC=MUD+1                                                          MCHB1140
      IF(M-MC)43,3,3                                                    MCHB1150
    3 MR=M-MUD                                                          MCHB1160
      IER=0                                                             MCHB1170
C                                                                       MCHB1180
C        MC IS THE MAXIMUM NUMBER OF ELEMENTS IN THE ROWS OF ARRAY A    MCHB1190
C        MR IS THE INDEX OF THE LAST ROW IN ARRAY A WITH MC ELEMENTS    MCHB1200
C                                                                       MCHB1210
C     ******************************************************************MCHB1220
C                                                                       MCHB1230
C        START FACTORIZATION OF MATRIX A                                MCHB1240
      IF(IOP)24,4,4                                                     MCHB1250
    4 IEND=0                                                            MCHB1260
      LLDST=MUD                                                         MCHB1270
      DO 23 K=1,M                                                       MCHB1280
      IST=IEND+1                                                        MCHB1290
      IEND=IST+MUD                                                      MCHB1300
      J=K-MR                                                            MCHB1310
      IF(J)6,6,5                                                        MCHB1320
    5 IEND=IEND-J                                                       MCHB1330
    6 IF(J-1)8,8,7                                                      MCHB1340
    7 LLDST=LLDST-1                                                     MCHB1350
    8 LMAX=MUD                                                          MCHB1360
      J=MC-K                                                            MCHB1370
      IF(J)10,10,9                                                      MCHB1380
    9 LMAX=LMAX-J                                                       MCHB1390
   10 ID=0                                                              MCHB1400
      TOL=A(IST)*EPS                                                    MCHB1410
C                                                                       MCHB1420
C        START FACTORIZATION-LOOP OVER K-TH ROW                         MCHB1430
      DO 23 I=IST,IEND                                                  MCHB1440
      SUM=0.D0                                                          MCHB1450
      IF(LMAX)14,14,11                                                  MCHB1460
C                                                                       MCHB1470
C        PREPARE INNER LOOP                                             MCHB1480
   11 LL=IST                                                            MCHB1490
      LLD=LLDST                                                         MCHB1500
C                                                                       MCHB1510
C        START INNER LOOP                                               MCHB1520
      DO 13 L=1,LMAX                                                    MCHB1530
      LL=LL-LLD                                                         MCHB1540
      LLL=LL+ID                                                         MCHB1550
      SUM=SUM+A(LL)*A(LLL)                                              MCHB1560
      IF(LLD-MUD)12,13,13                                               MCHB1570
   12 LLD=LLD+1                                                         MCHB1580
   13 CONTINUE                                                          MCHB1590
C        END OF INNER LOOP                                              MCHB1600
C                                                                       MCHB1610
C        TRANSFORM ELEMENT A(I)                                         MCHB1620
   14 SUM=DBLE(A(I))-SUM                                                MCHB1630
      IF(I-IST)15,15,20                                                 MCHB1640
C                                                                       MCHB1650
C        A(I) IS DIAGONAL ELEMENT. ERROR TEST.                          MCHB1660
   15 IF(SUM)43,43,16                                                   MCHB1670
C                                                                       MCHB1680
C        TEST ON LOSS OF SIGNIFICANT DIGITS AND WARNING                 MCHB1690
   16 IF(SUM-TOL)17,17,19                                               MCHB1700
   17 IF(IER)18,18,19                                                   MCHB1710
   18 IER=K-1                                                           MCHB1720
C                                                                       MCHB1730
C        COMPUTATION OF PIVOT ELEMENT                                   MCHB1740
   19 PIV=DSQRT(SUM)                                                    MCHB1750
      A(I)=PIV                                                          MCHB1760
      PIV=1.D0/PIV                                                      MCHB1770
      GO TO 21                                                          MCHB1780
C                                                                       MCHB1790
C        A(I) IS NOT DIAGONAL ELEMENT                                   MCHB1800
   20 A(I)=SUM*PIV                                                      MCHB1810
C                                                                       MCHB1820
C        UPDATE ID AND LMAX                                             MCHB1830
   21 ID=ID+1                                                           MCHB1840
      IF(ID-J)23,23,22                                                  MCHB1850
   22 LMAX=LMAX-1                                                       MCHB1860
   23 CONTINUE                                                          MCHB1870
C                                                                       MCHB1880
C        END OF FACTORIZATION-LOOP OVER K-TH ROW                        MCHB1890
C        END OF FACTORIZATION OF MATRIX A                               MCHB1900
C                                                                       MCHB1910
C     ******************************************************************MCHB1920
C                                                                       MCHB1930
C        PREPARE MATRIX DIVISIONS                                       MCHB1940
      IF(IOP)24,44,24                                                   MCHB1950
   24 ID=N*M                                                            MCHB1960
      IEND=IABS(IOP)-2                                                  MCHB1970
      IF(IEND)25,35,25                                                  MCHB1980
C                                                                       MCHB1990
C     ******************************************************************MCHB2000
C                                                                       MCHB2010
C        START DIVISION BY TRANSPOSE OF MATRIX TU (TU IS STORED IN      MCHB2020
C        LOCATIONS OF A)                                                MCHB2030
   25 IST=1                                                             MCHB2040
      LMAX=0                                                            MCHB2050
      J=-MR                                                             MCHB2060
      LLDST=MUD                                                         MCHB2070
      DO 34 K=1,M                                                       MCHB2080
      PIV=A(IST)                                                        MCHB2090
      IF(PIV)26,43,26                                                   MCHB2100
   26 PIV=1.D0/PIV                                                      MCHB2110
C                                                                       MCHB2120
C        START BACKSUBSTITUTION-LOOP FOR K-TH ROW OF MATRIX R           MCHB2130
      DO 30 I=K,ID,M                                                    MCHB2140
      SUM=0.D0                                                          MCHB2150
      IF(LMAX)30,30,27                                                  MCHB2160
C                                                                       MCHB2170
C        PREPARE INNER LOOP                                             MCHB2180
   27 LL=IST                                                            MCHB2190
      LLL=I                                                             MCHB2200
      LLD=LLDST                                                         MCHB2210
C                                                                       MCHB2220
C        START INNER LOOP                                               MCHB2230
      DO 29 L=1,LMAX                                                    MCHB2240
      LL=LL-LLD                                                         MCHB2250
      LLL=LLL-1                                                         MCHB2260
      SUM=SUM+A(LL)*R(LLL)                                              MCHB2270
      IF(LLD-MUD)28,29,29                                               MCHB2280
   28 LLD=LLD+1                                                         MCHB2290
   29 CONTINUE                                                          MCHB2300
C        END OF INNER LOOP                                              MCHB2310
C                                                                       MCHB2320
C        TRANSFORM ELEMENT R(I)                                         MCHB2330
   30 R(I)=PIV*(DBLE(R(I))-SUM)                                         MCHB2340
C        END OF BACKSUBSTITUTION-LOOP FOR K-TH ROW OF MATRIX R          MCHB2350
C                                                                       MCHB2360
C        UPDATE PARAMETERS LMAX, IST AND LLDST                          MCHB2370
      IF(MC-K)32,32,31                                                  MCHB2380
   31 LMAX=K                                                            MCHB2390
   32 IST=IST+MC                                                        MCHB2400
      J=J+1                                                             MCHB2410
      IF(J)34,34,33                                                     MCHB2420
   33 IST=IST-J                                                         MCHB2430
      LLDST=LLDST-1                                                     MCHB2440
   34 CONTINUE                                                          MCHB2450
C                                                                       MCHB2460
C        END OF DIVISION BY TRANSPOSE OF MATRIX TU                      MCHB2470
C                                                                       MCHB2480
C     ******************************************************************MCHB2490
C                                                                       MCHB2500
C        START DIVISION BY MATRIX TU (TU IS STORED ON LOCATIONS OF A)   MCHB2510
      IF(IEND)35,35,44                                                  MCHB2520
   35 IST=M+(MUD*(M+M-MC))/2+1                                          MCHB2530
      LMAX=0                                                            MCHB2540
      K=M                                                               MCHB2550
   36 IEND=IST-1                                                        MCHB2560
      IST=IEND-LMAX                                                     MCHB2570
      PIV=A(IST)                                                        MCHB2580
      IF(PIV)37,43,37                                                   MCHB2590
   37 PIV=1.D0/PIV                                                      MCHB2600
      L=IST+1                                                           MCHB2610
C                                                                       MCHB2620
C        START BACKSUBSTITUTION-LOOP FOR K-TH ROW OF MATRIX R           MCHB2630
      DO 40 I=K,ID,M                                                    MCHB2640
      SUM=0.D0                                                          MCHB2650
      IF(LMAX)40,40,38                                                  MCHB2660
   38 LLL=I                                                             MCHB2670
C                                                                       MCHB2680
C        START INNER LOOP                                               MCHB2690
      DO 39 LL=L,IEND                                                   MCHB2700
      LLL=LLL+1                                                         MCHB2710
   39 SUM=SUM+A(LL)*R(LLL)                                              MCHB2720
C        END OF INNER LOOP                                              MCHB2730
C                                                                       MCHB2740
C        TRANSFORM ELEMENT R(I)                                         MCHB2750
   40 R(I)=PIV*(DBLE(R(I))-SUM)                                         MCHB2760
C        END OF BACKSUBSTITUTION-LOOP FOR K-TH ROW OF MATRIX R          MCHB2770
C                                                                       MCHB2780
C                                                                       MCHB2790
C        UPDATE PARAMETERS LMAX AND K                                   MCHB2800
      IF(K-MR)42,42,41                                                  MCHB2810
   41 LMAX=LMAX+1                                                       MCHB2820
   42 K=K-1                                                             MCHB2830
      IF(K)44,44,36                                                     MCHB2840
C                                                                       MCHB2850
C        END OF DIVISION BY MATRIX TU                                   MCHB2860
C                                                                       MCHB2870
C     ******************************************************************MCHB2880
C                                                                       MCHB2890
C        ERROR EXIT IN CASE OF WRONG INPUT PARAMETERS OR PIVOT ELEMENT  MCHB2900
C        LESS THAN OR EQUAL TO ZERO                                     MCHB2910
   43 IER=-1                                                            MCHB2920
   44 RETURN                                                            MCHB2930
      END                                                               MCHB2940
