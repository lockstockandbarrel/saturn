C                                                                       DMCH  10
C     ..................................................................DMCH  20
C                                                                       DMCH  30
C        SUBROUTINE DMCHB                                               DMCH  40
C                                                                       DMCH  50
C        PURPOSE                                                        DMCH  60
C           FOR A GIVEN POSITIVE-DEFINITE M BY M MATRIX A WITH SYMMETRICDMCH  70
C           BAND STRUCTURE AND - IF NECESSARY - A GIVEN GENERAL M BY N  DMCH  80
C           MATRIX R, THE FOLLOWING CALCULATIONS (DEPENDENT ON THE      DMCH  90
C           VALUE OF THE DECISION PARAMETER IOP) ARE PERFORMED          DMCH 100
C           (1) MATRIX A IS FACTORIZED (IF IOP IS NOT NEGATIVE), THAT   DMCH 110
C               MEANS BAND MATRIX TU WITH UPPER CODIAGONALS ONLY IS     DMCH 120
C               GENERATED ON THE LOCATIONS OF A SUCH THAT               DMCH 130
C               TRANSPOSE(TU)*TU=A.                                     DMCH 140
C           (2) MATRIX R IS MULTIPLIED ON THE LEFT BY INVERSE(TU)       DMCH 150
C               AND/OR INVERSE(TRANSPOSE(TU)) AND THE RESULT IS STORED  DMCH 160
C               IN THE LOCATIONS OF R.                                  DMCH 170
C           THIS SUBROUTINE ESPECIALLY CAN BE USED TO SOLVE THE SYSTEM  DMCH 180
C           OF SIMULTANEOUS LINEAR EQUATIONS A*X=R WITH POSITIVE-       DMCH 190
C           DEFINITE COEFFICIENT MATRIX A OF SYMMETRIC BAND STRUCTURE.  DMCH 200
C                                                                       DMCH 210
C        USAGE                                                          DMCH 220
C           CALL DMCHB (R,A,M,N,MUD,IOP,EPS,IER)                        DMCH 230
C                                                                       DMCH 240
C        DESCRIPTION OF PARAMETERS                                      DMCH 250
C           R      - INPUT IN CASES IOP=-3,-2,-1,1,2,3  DOUBLE PRECISIONDMCH 260
C                          M BY N RIGHT HAND SIDE MATRIX,               DMCH 270
C                          IN CASE IOP=0  IRRELEVANT.                   DMCH 280
C                    OUTPUT IN CASES IOP=1,-1  INVERSE(A)*R,            DMCH 290
C                           IN CASES IOP=2,-2  INVERSE(TU)*R,           DMCH 300
C                           IN CASES IOP=3,-3  INVERSE(TRANSPOSE(TU))*R,DMCH 310
C                           IN CASE  IOP=0     UNCHANGED.               DMCH 320
C           A      - INPUT IN CASES IOP=0,1,2,3  DOUBLE PRECISION M BY MDMCH 330
C                          POSITIVE-DEFINITE COEFFICIENT MATRIX OF      DMCH 340
C                          SYMMETRIC BAND STRUCTURE STORED IN           DMCH 350
C                          COMPRESSED FORM (SEE REMARKS),               DMCH 360
C                          IN CASES IOP=-1,-2,-3 DOUBLE PRECISION M BY MDMCH 370
C                          BAND MATRIX TU WITH UPPER CODIAGONALS ONLY,  DMCH 380
C                          STORED IN COMPRESSED FORM (SEE REMARKS).     DMCH 390
C                    OUTPUT IN ALL CASES  BAND MATRIX TU WITH UPPER     DMCH 400
C                           CODIAGONALS ONLY, STORED IN COMPRESSED FORM DMCH 410
C                           (THAT MEANS UNCHANGED IF IOP=-1,-2,-3).     DMCH 420
C           M      - INPUT VALUE SPECIFYING THE NUMBER OF ROWS AND      DMCH 430
C                    COLUMNS OF A AND THE NUMBER OF ROWS OF R.          DMCH 440
C           N      - INPUT VALUE SPECIFYING THE NUMBER OF COLUMNS OF R  DMCH 450
C                    (IRRELEVANT IN CASE IOP=0).                        DMCH 460
C           MUD    - INPUT VALUE SPECIFYING THE NUMBER OF UPPER         DMCH 470
C                    CODIAGONALS OF A.                                  DMCH 480
C           IOP    - ONE OF THE VALUES -3,-2,-1,0,1,2,3 GIVEN AS INPUT  DMCH 490
C                    AND USED AS DECISION PARAMETER.                    DMCH 500
C           EPS    - SINGLE PRECISION INPUT VALUE USED AS RELATIVE      DMCH 510
C                    TOLERANCE FOR TEST ON LOSS OF SIGNIFICANT DIGITS.  DMCH 520
C           IER    - RESULTING ERROR PARAMETER CODED AS FOLLOWS         DMCH 530
C                     IER=0  - NO ERROR,                                DMCH 540
C                     IER=-1 - NO RESULT BECAUSE OF WRONG INPUT         DMCH 550
C                              PARAMETERS M,MUD,IOP (SEE REMARKS),      DMCH 560
C                              OR BECAUSE OF A NONPOSITIVE RADICAND AT  DMCH 570
C                              SOME FACTORIZATION STEP,                 DMCH 580
C                              OR BECAUSE OF A ZERO DIAGONAL ELEMENT    DMCH 590
C                              AT SOME DIVISION STEP.                   DMCH 600
C                     IER=K  - WARNING DUE TO POSSIBLE LOSS OF SIGNIFI- DMCH 610
C                              CANCE INDICATED AT FACTORIZATION STEP K+1DMCH 620
C                              WHERE RADICAND WAS NO LONGER GREATER     DMCH 630
C                              THAN EPS*A(K+1,K+1).                     DMCH 640
C                                                                       DMCH 650
C        REMARKS                                                        DMCH 660
C           UPPER PART OF SYMMETRIC BAND MATRIX A CONSISTING OF MAIN    DMCH 670
C           DIAGONAL AND MUD UPPER CODIAGONALS (RESP. BAND MATRIX TU    DMCH 680
C           CONSISTING OF MAIN DIAGONAL AND MUD UPPER CODIAGONALS)      DMCH 690
C           IS ASSUMED TO BE STORED IN COMPRESSED FORM, I.E. ROWWISE    DMCH 700
C           IN TOTALLY NEEDED M+MUD*(2M-MUD-1)/2 SUCCESSIVE STORAGE     DMCH 710
C           LOCATIONS. ON RETURN UPPER BAND FACTOR TU (ON THE LOCATIONS DMCH 720
C           OF A) IS STORED IN THE SAME WAY.                            DMCH 730
C           RIGHT HAND SIDE MATRIX R IS ASSUMED TO BE STORED COLUMNWISE DMCH 740
C           IN N*M SUCCESSIVE STORAGE LOCATIONS. ON RETURN RESULT MATRIXDMCH 750
C           INVERSE(A)*R OR INVERSE(TU)*R OR INVERSE(TRANSPOSE(TU))*R   DMCH 760
C           IS STORED COLUMNWISE TOO ON THE LOCATIONS OF R.             DMCH 770
C           INPUT PARAMETERS M, MUD, IOP SHOULD SATISFY THE FOLLOWING   DMCH 780
C           RESTRICTIONS     MUD NOT LESS THAN ZERO,                    DMCH 790
C                            1+MUD NOT GREATER THAN M,                  DMCH 800
C                            ABS(IOP) NOT GREATER THAN 3.               DMCH 810
C           NO ACTION BESIDES ERROR MESSAGE IER=-1 TAKES PLACE IF THESE DMCH 820
C           RESTRICTIONS ARE NOT SATISFIED.                             DMCH 830
C           THE PROCEDURE GIVES RESULTS IF THE RESTRICTIONS ON INPUT    DMCH 840
C           PARAMETERS ARE SATISFIED, IF RADICANDS AT ALL FACTORIZATION DMCH 850
C           STEPS ARE POSITIVE AND/OR IF ALL DIAGONAL ELEMENTS OF       DMCH 860
C           UPPER BAND FACTOR TU ARE NONZERO.                           DMCH 870
C                                                                       DMCH 880
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DMCH 890
C           NONE                                                        DMCH 900
C                                                                       DMCH 910
C        METHOD                                                         DMCH 920
C           FACTORIZATION IS DONE USING CHOLESKY-S SQUARE-ROOT METHOD,  DMCH 930
C           WHICH GENERATES THE UPPER BAND MATRIX TU SUCH THAT          DMCH 940
C           TRANSPOSE(TU)*TU=A. TU IS RETURNED AS RESULT ON THE         DMCH 950
C           LOCATIONS OF A. FURTHER, DEPENDENT ON THE ACTUAL VALUE OF   DMCH 960
C           IOP, DIVISION OF R BY TRANSPOSE(TU) AND/OR TU IS PERFORMED  DMCH 970
C           AND THE RESULT IS RETURNED ON THE LOCATIONS OF R.           DMCH 980
C           FOR REFERENCE, SEE H. RUTISHAUSER, ALGORITHMUS 1 - LINEARES DMCH 990
C           GLEICHUNGSSYSTEM MIT SYMMETRISCHER POSITIV-DEFINITER        DMCH1000
C           BANDMATRIX NACH CHOLESKY - , COMPUTING (ARCHIVES FOR        DMCH1010
C           ELECTRONIC COMPUTING), VOL.1, ISS.1 (1966), PP.77-78.       DMCH1020
C                                                                       DMCH1030
C     ..................................................................DMCH1040
C                                                                       DMCH1050
      SUBROUTINE DMCHB(R,A,M,N,MUD,IOP,EPS,IER)                         DMCH1060
C                                                                       DMCH1070
C                                                                       DMCH1080
      DIMENSION R(1),A(1)                                               DMCH1090
      DOUBLE PRECISION TOL,SUM,PIV,R,A                                  DMCH1100
C                                                                       DMCH1110
C        TEST ON WRONG INPUT PARAMETERS                                 DMCH1120
      IF(IABS(IOP)-3)1,1,43                                             DMCH1130
    1 IF(MUD)43,2,2                                                     DMCH1140
    2 MC=MUD+1                                                          DMCH1150
      IF(M-MC)43,3,3                                                    DMCH1160
    3 MR=M-MUD                                                          DMCH1170
      IER=0                                                             DMCH1180
C                                                                       DMCH1190
C        MC IS THE MAXIMUM NUMBER OF ELEMENTS IN THE ROWS OF ARRAY A    DMCH1200
C        MR IS THE INDEX OF THE LAST ROW IN ARRAY A WITH MC ELEMENTS    DMCH1210
C                                                                       DMCH1220
C     ******************************************************************DMCH1230
C                                                                       DMCH1240
C        START FACTORIZATION OF MATRIX A                                DMCH1250
      IF(IOP)24,4,4                                                     DMCH1260
    4 IEND=0                                                            DMCH1270
      LLDST=MUD                                                         DMCH1280
      DO 23 K=1,M                                                       DMCH1290
      IST=IEND+1                                                        DMCH1300
      IEND=IST+MUD                                                      DMCH1310
      J=K-MR                                                            DMCH1320
      IF(J)6,6,5                                                        DMCH1330
    5 IEND=IEND-J                                                       DMCH1340
    6 IF(J-1)8,8,7                                                      DMCH1350
    7 LLDST=LLDST-1                                                     DMCH1360
    8 LMAX=MUD                                                          DMCH1370
      J=MC-K                                                            DMCH1380
      IF(J)10,10,9                                                      DMCH1390
    9 LMAX=LMAX-J                                                       DMCH1400
   10 ID=0                                                              DMCH1410
      TOL=A(IST)*EPS                                                    DMCH1420
C                                                                       DMCH1430
C        START FACTORIZATION-LOOP OVER K-TH ROW                         DMCH1440
      DO 23 I=IST,IEND                                                  DMCH1450
      SUM=0.D0                                                          DMCH1460
      IF(LMAX)14,14,11                                                  DMCH1470
C                                                                       DMCH1480
C        PREPARE INNER LOOP                                             DMCH1490
   11 LL=IST                                                            DMCH1500
      LLD=LLDST                                                         DMCH1510
C                                                                       DMCH1520
C        START INNER LOOP                                               DMCH1530
      DO 13 L=1,LMAX                                                    DMCH1540
      LL=LL-LLD                                                         DMCH1550
      LLL=LL+ID                                                         DMCH1560
      SUM=SUM+A(LL)*A(LLL)                                              DMCH1570
      IF(LLD-MUD)12,13,13                                               DMCH1580
   12 LLD=LLD+1                                                         DMCH1590
   13 CONTINUE                                                          DMCH1600
C        END OF INNER LOOP                                              DMCH1610
C                                                                       DMCH1620
C        TRANSFORM ELEMENT A(I)                                         DMCH1630
   14 SUM=A(I)-SUM                                                      DMCH1640
      IF(I-IST)15,15,20                                                 DMCH1650
C                                                                       DMCH1660
C        A(I) IS DIAGONAL ELEMENT. ERROR TEST.                          DMCH1670
   15 IF(SUM)43,43,16                                                   DMCH1680
C                                                                       DMCH1690
C        TEST ON LOSS OF SIGNIFICANT DIGITS AND WARNING                 DMCH1700
   16 IF(SUM-TOL)17,17,19                                               DMCH1710
   17 IF(IER)18,18,19                                                   DMCH1720
   18 IER=K-1                                                           DMCH1730
C                                                                       DMCH1740
C        COMPUTATION OF PIVOT ELEMENT                                   DMCH1750
   19 PIV=DSQRT(SUM)                                                    DMCH1760
      A(I)=PIV                                                          DMCH1770
      PIV=1.D0/PIV                                                      DMCH1780
      GO TO 21                                                          DMCH1790
C                                                                       DMCH1800
C        A(I) IS NOT DIAGONAL ELEMENT                                   DMCH1810
   20 A(I)=SUM*PIV                                                      DMCH1820
C                                                                       DMCH1830
C        UPDATE ID AND LMAX                                             DMCH1840
   21 ID=ID+1                                                           DMCH1850
      IF(ID-J)23,23,22                                                  DMCH1860
   22 LMAX=LMAX-1                                                       DMCH1870
   23 CONTINUE                                                          DMCH1880
C                                                                       DMCH1890
C        END OF FACTORIZATION-LOOP OVER K-TH ROW                        DMCH1900
C        END OF FACTORIZATION OF MATRIX A                               DMCH1910
C                                                                       DMCH1920
C     ******************************************************************DMCH1930
C                                                                       DMCH1940
C        PREPARE MATRIX DIVISIONS                                       DMCH1950
      IF(IOP)24,44,24                                                   DMCH1960
   24 ID=N*M                                                            DMCH1970
      IEND=IABS(IOP)-2                                                  DMCH1980
      IF(IEND)25,35,25                                                  DMCH1990
C                                                                       DMCH2000
C     ******************************************************************DMCH2010
C                                                                       DMCH2020
C        START DIVISION BY TRANSPOSE OF MATRIX TU (TU IS STORED IN      DMCH2030
C        LOCATIONS OF A)                                                DMCH2040
   25 IST=1                                                             DMCH2050
      LMAX=0                                                            DMCH2060
      J=-MR                                                             DMCH2070
      LLDST=MUD                                                         DMCH2080
      DO 34 K=1,M                                                       DMCH2090
      PIV=A(IST)                                                        DMCH2100
      IF(PIV)26,43,26                                                   DMCH2110
   26 PIV=1.D0/PIV                                                      DMCH2120
C                                                                       DMCH2130
C        STA-T BACKSUBSTITUTION-LOOP FOR K-TH ROW OF MATRIX R           DMCH2140
      DO 30 I=K,ID,M                                                    DMCH2150
      SUM=0.D0                                                          DMCH2160
      IF(LMAX)30,30,27                                                  DMCH2170
C                                                                       DMCH2180
C        PREPARE INNER LOOP                                             DMCH2190
   27 LL=IST                                                            DMCH2200
      LLL=I                                                             DMCH2210
      LLD=LLDST                                                         DMCH2220
C                                                                       DMCH2230
C        START INNER LOOP                                               DMCH2240
      DO 29 L=1,LMAX                                                    DMCH2250
      LL=LL-LLD                                                         DMCH2260
      LLL=LLL-1                                                         DMCH2270
      SUM=SUM+A(LL)*R(LLL)                                              DMCH2280
      IF(LLD-MUD)28,29,29                                               DMCH2290
   28 LLD=LLD+1                                                         DMCH2300
   29 CONTINUE                                                          DMCH2310
C        END OF INNER LOOP                                              DMCH2320
C                                                                       DMCH2330
C        TRANSFORM ELEMENT R(I)                                         DMCH2340
   30 R(I)=PIV*(R(I)-SUM)                                               DMCH2350
C        END OF BACKSUBSTITUTION-LOOP FOR K-TH ROW OF MATRIX R          DMCH2360
C                                                                       DMCH2370
C        UPDATE PARAMETERS LMAX, IST AND LLDST                          DMCH2380
      IF(MC-K)32,32,31                                                  DMCH2390
   31 LMAX=K                                                            DMCH2400
   32 IST=IST+MC                                                        DMCH2410
      J=J+1                                                             DMCH2420
      IF(J)34,34,33                                                     DMCH2430
   33 IST=IST-J                                                         DMCH2440
      LLDST=LLDST-1                                                     DMCH2450
   34 CONTINUE                                                          DMCH2460
C                                                                       DMCH2470
C        END OF DIVISION BY TRANSPOSE OF MATRIX TU                      DMCH2480
C                                                                       DMCH2490
C     ******************************************************************DMCH2500
C                                                                       DMCH2510
C        START DIVISION BY MATRIX TU (TU IS STORED ON LOCATIONS OF A)   DMCH2520
      IF(IEND)35,35,44                                                  DMCH2530
   35 IST=M+(MUD*(M+M-MC))/2+1                                          DMCH2540
      LMAX=0                                                            DMCH2550
      K=M                                                               DMCH2560
   36 IEND=IST-1                                                        DMCH2570
      IST=IEND-LMAX                                                     DMCH2580
      PIV=A(IST)                                                        DMCH2590
      IF(PIV)37,43,37                                                   DMCH2600
   37 PIV=1.D0/PIV                                                      DMCH2610
      L=IST+1                                                           DMCH2620
C                                                                       DMCH2630
C        START BACKSUBSTITUTION-LOOP FOR K-TH ROW OF MATRIX R           DMCH2640
      DO 40 I=K,ID,M                                                    DMCH2650
      SUM=0.D0                                                          DMCH2660
      IF(LMAX)40,40,38                                                  DMCH2670
   38 LLL=I                                                             DMCH2680
C                                                                       DMCH2690
C        START INNER LOOP                                               DMCH2700
      DO 39 LL=L,IEND                                                   DMCH2710
      LLL=LLL+1                                                         DMCH2720
   39 SUM=SUM+A(LL)*R(LLL)                                              DMCH2730
C        END OF INNER LOOP                                              DMCH2740
C                                                                       DMCH2750
C        TRANSFORM ELEMENT R(I)                                         DMCH2760
   40 R(I)=PIV*(R(I)-SUM)                                               DMCH2770
C        END OF BACKSUBSTITUTION-LOOP FOR K-TH ROW OF MATRIX R          DMCH2780
C                                                                       DMCH2790
C        UPDATE PARAMETERS LMAX AND K                                   DMCH2800
      IF(K-MR)42,42,41                                                  DMCH2810
   41 LMAX=LMAX+1                                                       DMCH2820
   42 K=K-1                                                             DMCH2830
      IF(K)44,44,36                                                     DMCH2840
C                                                                       DMCH2850
C        END OF DIVISION BY MATRIX TU                                   DMCH2860
C                                                                       DMCH2870
C     ******************************************************************DMCH2880
C                                                                       DMCH2890
C        ERROR EXIT IN CASE OF WRONG INPUT PARAMETERS OR PIVOT ELEMENT  DMCH2900
C        LESS THAN OR EQUAL TO ZERO                                     DMCH2910
   43 IER=-1                                                            DMCH2920
   44 RETURN                                                            DMCH2930
      END                                                               DMCH2940
