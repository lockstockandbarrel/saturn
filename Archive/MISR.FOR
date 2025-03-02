C                                                                       MISR  10
C     ..................................................................MISR  20
C                                                                       MISR  30
C        SUBROUTINE MISR                                                MISR  40
C                                                                       MISR  50
C        PURPOSE                                                        MISR  60
C           COMPUTE MEANS, STANDARD DEVIATIONS, SKEWNESS AND KURTOSIS,  MISR  70
C           CORRELATION COEFFICIENTS, REGRESSION COEFFICIENTS, AND      MISR  80
C           STANDARD ERRORS OF REGRESSION COEFFICIENTS WHEN THERE ARE   MISR  90
C           MISSING DATA POINTS.  THE USER IDENTIFIES THE MISSING DATA  MISR 100
C           BY MEANS OF A NUMERIC CODE.  THOSE VALUES HAVING THIS CODE  MISR 110
C           ARE SKIPPED IN COMPUTING THE STATISTICS.  IN THE CASE OF THEMISR 120
C           CORRELATION COEFFICIENTS, ANY PAIR OF VALUES ARE SKIPPED IF MISR 130
C           EITHER ONE OF THEM ARE MISSING.                             MISR 140
C                                                                       MISR 150
C        USAGE                                                          MISR 160
C           CALL MISR (NO,M,X,CODE,XBAR,STD,SKEW,CURT,R,N,A,B,S,IER)    MISR 170
C                                                                       MISR 180
C        DESCRIPTION OF PARAMETERS                                      MISR 190
C           NO   - NUMBER OF OBSERVATIONS                               MISR 200
C           M    - NUMBER OF VARIABLES                                  MISR 210
C           X    - INPUT DATA MATRIX OF SIZE NO X M.                    MISR 220
C           CODE - INPUT VECTOR OF LENGTH M, WHICH CONTAINS A NUMERIC   MISR 230
C                  MISSING DATA CODE FOR EACH VARIABLE. ANY OBSERVATION MISR 240
C                  FOR A GIVEN VARIABLE HAVING A VALUE EQUAL TO THE CODEMISR 250
C                  WILL BE DROPPED FOR THE COMPUTATIONS.                MISR 260
C           XBAR - OUTPUT VECTOR OF LENGTH M CONTAINING MEANS           MISR 270
C           STD  - OUTPUT VECTOR OF LENGTH M CONTAINING STANDARD DEVI-  MISR 280
C                  ATIONS                                               MISR 290
C           SKEW - OUTPUT VECTOR OF LENGTH M CONTAINING SKEWNESS        MISR 300
C           CURT - OUTPUT VECTOR OF LENGTH M CONTAINING KURTOSIS        MISR 310
C           R    - OUTPUT MATRIX OF PRODUCT-MOMENT CORRELATION          MISR 320
C                  COEFFICIENTS.  THIS WILL BE THE UPPER TRIANGULAR     MISR 330
C                  MATRIX ONLY, SINCE THE M X M MATRIX OF COEFFICIENTS  MISR 340
C                  IS SYMMETRIC. (STORAGE MODE 1)                       MISR 350
C           N    - OUTPUT MATRIX OF NUMBER OF PAIRS OF OBSERVATIONS USEDMISR 360
C                  IN COMPUTING THE CORRELATION COEFFICIENTS.  ONLY THE MISR 370
C                  UPPER TRIANGULAR PORTION OF THE MATRIX IS GIVEN.     MISR 380
C                  (STORAGE MODE 1)                                     MISR 390
C           A    - OUTPUT MATRIX (M BY M)  CONTAINING INTERCEPTS OF     MISR 400
C                  REGRESSION LINES (A) OF THE FORM Y=A+BX.  THE FIRST  MISR 410
C                  SUBSCRIPT OF THIS MATRIX REFERS TO THE INDEPENDENT   MISR 420
C                  VARIABLE AND THE SECOND TO THE DEPENDENT VARIABLE.   MISR 430
C                  FOR EXAMPLE, A(1,3) CONTAINS THE INTERCEPT OF THE    MISR 440
C                  REGRESSION LINE FOR TWO VARIABLES WHERE VARIABLE 1   MISR 450
C                  IS INDEPENDENT AND VARIABLE 3 IS DEPENDENT.  NOTE    MISR 460
C                  THAT MATRIX A IS STORED IN A VECTOR FORM.            MISR 470
C           B    - OUTPUT MATRIX (M BY M)  CONTAINING REGRESSION        MISR 480
C                  COEFFICIENTS (B) CORRESPONDING TO THE VALUES OF      MISR 490
C                  INTERCEPTS CONTAINED IN THE OUTPUT MATRIX A.         MISR 500
C           S    - OUTPUT MATRIX (M BY M)  CONTAINING STANDARD ERRORS   MISR 510
C                  OF REGRESSION COEFFICIENTS CORRESPONDING TO THE      MISR 520
C                  COEFFICIENTS CONTAINED IN THE OUTPUT MATRIX B.       MISR 530
C           IER  - 0, NO ERROR.                                         MISR 540
C                  1, IF NUMBER OF NON-MISSING DATA ELEMENTS FOR J-TH   MISR 550
C                     VARIABLE IS TWO OR LESS.  IN THIS CASE, STD(J),   MISR 560
C                     SKEW(J), AND CURT(J) ARE SET TO 10**75.  ALL      MISR 570
C                     VALUES OF R, A, B, AND S RELATED TO THIS VARIABLE MISR 580
C                     ARE ALSO SET TO 10**75.                           MISR 590
C                  2, IF VARIANCE OF J-TH VARIABLE IS LESS THAN         MISR 600
C                     10**(-20).  IN THIS CASE, STD(J), SKEW(J), AND    MISR 610
C                     CURT(J) ARE SET TO 10**75.  ALL VALUES OF R, A,   MISR 620
C                     B, AND S RELATED TO THIS VARIABLE ARE ALSO SET TO MISR 630
C                     10**75.                                           MISR 640
C                                                                       MISR 650
C        REMARKS                                                        MISR 660
C           THIS SUBROUTINE CANNOT DISTINGUISH A BLANK AND A ZERO.      MISR 670
C           THEREFORE, IF A BLANK IS SPECIFIED AS A MISSING DATA CODE INMISR 680
C           INPUT CARDS, IT WILL BE TREATED AS 0 (ZERO).                MISR 690
C                                                                       MISR 700
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  MISR 710
C           NONE                                                        MISR 720
C                                                                       MISR 730
C        METHOD                                                         MISR 740
C           LEAST SQUARES REGRESSION LINES AND PRODUCT-MOMENT CORRE-    MISR 750
C           LATION COEFFICIENTS ARE COMPUTED.                           MISR 760
C                                                                       MISR 770
C     ..................................................................MISR 780
C                                                                       MISR 790
      SUBROUTINE MISR (NO,M,X,CODE,XBAR,STD,SKEW,CURT,R,N,A,B,S,IER)    MISR 800
C                                                                       MISR 810
      DIMENSION X(1),CODE(1),XBAR(1),STD(1),SKEW(1),CURT(1),R(1),N(1)   MISR 820
      DIMENSION A(1),B(1),S(1)                                          MISR 830
C                                                                       MISR 840
C        COMPUTE MEANS                                                  MISR 850
C                                                                       MISR 860
      IER=0                                                             MISR 870
      L=0                                                               MISR 880
      DO 20 J=1,M                                                       MISR 890
      FN=0.0                                                            MISR 900
      XBAR(J)=0.0                                                       MISR 910
      DO 15 I=1,NO                                                      MISR 920
      L=L+1                                                             MISR 930
      IF(X(L)-CODE(J)) 12, 15, 12                                       MISR 940
   12 FN=FN+1.0                                                         MISR 950
      XBAR(J)=XBAR(J)+X(L)                                              MISR 960
   15 CONTINUE                                                          MISR 970
      IF(FN) 16, 16, 17                                                 MISR 980
   16 XBAR(J)=0.0                                                       MISR 990
      GO TO 20                                                          MISR1000
   17 XBAR(J)=XBAR(J)/FN                                                MISR1010
   20 CONTINUE                                                          MISR1020
C                                                                       MISR1030
C        SET-UP WORK AREAS AND TEST WHETHER DATA IS MISSING             MISR1040
C                                                                       MISR1050
      L=0                                                               MISR1060
      DO 55 J=1,M                                                       MISR1070
      LJJ=NO*(J-1)                                                      MISR1080
      SKEW(J)=0.0                                                       MISR1090
      CURT(J)=0.0                                                       MISR1100
      KI=M*(J-1)                                                        MISR1110
      KJ=J-M                                                            MISR1120
      DO 54 I=1,J                                                       MISR1130
      KI=KI+1                                                           MISR1140
      KJ=KJ+M                                                           MISR1150
      SUMX=0.0                                                          MISR1160
      SUMY=0.0                                                          MISR1170
      TI=0.0                                                            MISR1180
      TJ=0.0                                                            MISR1190
      TII=0.0                                                           MISR1200
      TJJ=0.0                                                           MISR1210
      TIJ=0.0                                                           MISR1220
      NIJ=0                                                             MISR1230
      LI=NO*(I-1)                                                       MISR1240
      LJ=LJJ                                                            MISR1250
      L=L+1                                                             MISR1260
      DO 38 K=1,NO                                                      MISR1270
      LI=LI+1                                                           MISR1280
      LJ=LJ+1                                                           MISR1290
      IF(X(LI)-CODE(I)) 30, 38, 30                                      MISR1300
   30 IF(X(LJ)-CODE(J)) 35, 38, 35                                      MISR1310
C                                                                       MISR1320
C        BOTH DATA ARE PRESENT                                          MISR1330
C                                                                       MISR1340
   35 XX=X(LI)-XBAR(I)                                                  MISR1350
      YY=X(LJ)-XBAR(J)                                                  MISR1360
      TI=TI+XX                                                          MISR1370
      TII=TII+XX**2                                                     MISR1380
      TJ=TJ+YY                                                          MISR1390
      TJJ=TJJ+YY**2                                                     MISR1400
      TIJ=TIJ+XX*YY                                                     MISR1410
      NIJ=NIJ+1                                                         MISR1420
      SUMX=SUMX+X(LI)                                                   MISR1430
      SUMY=SUMY+X(LJ)                                                   MISR1440
      IF(I-J) 38, 37, 37                                                MISR1450
   37 SKEW(J)=SKEW(J)+YY**3                                             MISR1460
      CURT(J)=CURT(J)+YY**4                                             MISR1470
   38 CONTINUE                                                          MISR1480
C                                                                       MISR1490
C        COMPUTE SUM OF CROSS-PRODUCTS OF DEVIATIONS                    MISR1500
C                                                                       MISR1510
      IF(NIJ) 40, 40, 39                                                MISR1520
   39 FN=NIJ                                                            MISR1530
      R(L)=TIJ-TI*TJ/FN                                                 MISR1540
      N(L)=NIJ                                                          MISR1550
      TII=TII-TI*TI/FN                                                  MISR1560
      TJJ=TJJ-TJ*TJ/FN                                                  MISR1570
C                                                                       MISR1580
C        COMPUTE STANDARD DEVIATION, SKEWNESS, AND KURTOSIS             MISR1590
C                                                                       MISR1600
   40 IF(I-J) 47, 41, 47                                                MISR1610
   41 IF(NIJ-2) 42,42,43                                                MISR1620
   42 IER=1                                                             MISR1630
      R(L)=1.0E75                                                       MISR1640
      A(KI)=1.0E75                                                      MISR1650
      B(KI)=1.0E75                                                      MISR1660
      S(KI)=1.0E75                                                      MISR1670
      GO TO 45                                                          MISR1680
C                                                                       MISR1690
   43 STD(J)=R(L)                                                       MISR1700
      R(L)=1.0                                                          MISR1710
      A(KI)=0.0                                                         MISR1720
      B(KI)=1.0                                                         MISR1730
      S(KI)=0.0                                                         MISR1740
C                                                                       MISR1750
      IF(STD(J)-(1.0E-20)) 44,44,46                                     MISR1760
   44 IER=2                                                             MISR1770
   45 STD(J)=1.0E75                                                     MISR1780
      SKEW(J)=1.0E75                                                    MISR1790
      CURT(J)=1.0E75                                                    MISR1800
      GO TO 55                                                          MISR1810
C                                                                       MISR1820
   46 WORK=STD(J)/FN                                                    MISR1830
      SKEW(J)=(SKEW(J)/FN)/(WORK*SQRT(WORK))                            MISR1840
      CURT(J)=((CURT(J)/FN)/WORK**2)-3.0                                MISR1850
      STD(J)=SQRT(STD(J)/(FN-1.0))                                      MISR1860
      GO TO 55                                                          MISR1870
C                                                                       MISR1880
C        COMPUTE REGRESSION COEFFICIENTS                                MISR1890
C                                                                       MISR1900
   47 IF(NIJ-2) 48,48,50                                                MISR1910
   48 IER=1                                                             MISR1920
   49 R(L)=1.0E75                                                       MISR1930
      A(KI)=1.0E75                                                      MISR1940
      B(KI)=1.0E75                                                      MISR1950
      S(KI)=1.0E75                                                      MISR1960
      A(KJ)=1.0E75                                                      MISR1970
      B(KJ)=1.0E75                                                      MISR1980
      S(KJ)=1.0E75                                                      MISR1990
      GO TO 54                                                          MISR2000
C                                                                       MISR2010
   50 IF(TII-(1.0E-20)) 52,52,51                                        MISR2020
   51 IF(TJJ-(1.0E-20)) 52,52,53                                        MISR2030
   52 IER=2                                                             MISR2040
      GO TO 49                                                          MISR2050
C                                                                       MISR2060
   53 SUMX=SUMX/FN                                                      MISR2070
      SUMY=SUMY/FN                                                      MISR2080
      B(KI)=R(L)/TII                                                    MISR2090
      A(KI)=SUMY-B(KI)*SUMX                                             MISR2100
      B(KJ)=R(L)/TJJ                                                    MISR2110
      A(KJ)=SUMX-B(KJ)*SUMY                                             MISR2120
C                                                                       MISR2130
C        COMPUTE CORRELATION COEFFICIENTS                               MISR2140
C                                                                       MISR2150
      R(L)=R(L)/(SQRT(TII)*SQRT(TJJ))                                   MISR2160
C                                                                       MISR2170
C        COMPUTE STANDARD ERRORS OF REGRESSION COEFFICIENTS             MISR2180
C                                                                       MISR2190
      RR=R(L)**2                                                        MISR2200
      SUMX=(TJJ-TJJ*RR)/(FN-2)                                          MISR2210
      S(KI)=SQRT(SUMX/TII)                                              MISR2220
      SUMY=(TII-TII*RR)/(FN-2)                                          MISR2230
      S(KJ)=SQRT(SUMY/TJJ)                                              MISR2240
C                                                                       MISR2250
   54 CONTINUE                                                          MISR2260
   55 CONTINUE                                                          MISR2270
C                                                                       MISR2280
      RETURN                                                            MISR2290
      END                                                               MISR2300
