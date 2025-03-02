C                                                                       CORR  10
C     ..................................................................CORR  20
C                                                                       CORR  30
C        SUBROUTINE CORRE                                               CORR  40
C                                                                       CORR  50
C        PURPOSE                                                        CORR  60
C           COMPUTE MEANS, STANDARD DEVIATIONS, SUMS OF CROSS-PRODUCTS  CORR  70
C           OF DEVIATIONS, AND CORRELATION COEFFICIENTS.                CORR  80
C                                                                       CORR  90
C        USAGE                                                          CORR 100
C           CALL CORRE (N,M,IO,X,XBAR,STD,RX,R,B,D,T)                   CORR 110
C                                                                       CORR 120
C        DESCRIPTION OF PARAMETERS                                      CORR 130
C           N     - NUMBER OF OBSERVATIONS. N MUST BE > OR = TO 2.      CORR 140
C           M     - NUMBER OF VARIABLES. M MUST BE > OR = TO 1.         CORR 150
C           IO    - OPTION CODE FOR INPUT DATA                          CORR 160
C                   0 IF DATA ARE TO BE READ IN FROM INPUT DEVICE IN THECORR 170
C                     SPECIAL SUBROUTINE NAMED DATA.  (SEE SUBROUTINES  CORR 180
C                     USED BY THIS SUBROUTINE BELOW.)                   CORR 190
C                   1 IF ALL DATA ARE ALREADY IN CORE.                  CORR 200
C           X     - IF IO=0, THE VALUE OF X IS 0.0.                     CORR 210
C                   IF IO=1, X IS THE INPUT MATRIX (N BY M) CONTAINING  CORR 220
C                            DATA.                                      CORR 230
C           XBAR  - OUTPUT VECTOR OF LENGTH M CONTAINING MEANS.         CORR 240
C           STD   - OUTPUT VECTOR OF LENGTH M CONTAINING STANDARD       CORR 250
C                   DEVIATIONS.                                         CORR 260
C           RX    - OUTPUT MATRIX (M X M) CONTAINING SUMS OF CROSS-     CORR 270
C                   PRODUCTS OF DEVIATIONS FROM MEANS.                  CORR 280
C           R     - OUTPUT MATRIX (ONLY UPPER TRIANGULAR PORTION OF THE CORR 290
C                   SYMMETRIC MATRIX OF M BY M) CONTAINING CORRELATION  CORR 300
C                   COEFFICIENTS.  (STORAGE MODE OF 1)                  CORR 310
C           B     - OUTPUT VECTOR OF LENGTH M CONTAINING THE DIAGONAL   CORR 320
C                   OF THE MATRIX OF SUMS OF CROSS-PRODUCTS OF          CORR 330
C                   DEVIATIONS FROM MEANS.                              CORR 340
C           D     - WORKING VECTOR OF LENGTH M.                         CORR 350
C           T     - WORKING VECTOR OF LENGTH M.                         CORR 360
C                                                                       CORR 370
C        REMARKS                                                        CORR 380
C           CORRE WILL NOT ACCEPT A CONSTANT VECTOR.                    CORR 390
C                                                                       CORR 400
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  CORR 410
C           DATA(M,D) - THIS SUBROUTINE MUST BE PROVIDED BY THE USER.   CORR 420
C                       (1) IF IO=0, THIS SUBROUTINE IS EXPECTED TO     CORR 430
C                           FURNISH AN OBSERVATION IN VECTOR D FROM AN  CORR 440
C                           EXTERNAL INPUT DEVICE.                      CORR 450
C                       (2) IF IO=1, THIS SUBROUTINE IS NOT USED BY     CORR 460
C                           CORRE BUT MUST EXIST IN JOB DECK. IF USER   CORR 470
C                           HAS NOT SUPPLIED A SUBROUTINE NAMED DATA,   CORR 480
C                           THE FOLLOWING IS SUGGESTED.                 CORR 490
C                                SUBROUTINE DATA                        CORR 500
C                                RETURN                                 CORR 510
C                                END                                    CORR 520
C                                                                       CORR 530
C        METHOD                                                         CORR 540
C           PRODUCT-MOMENT CORRELATION COEFFICIENTS ARE COMPUTED.       CORR 550
C                                                                       CORR 560
C     ..................................................................CORR 570
C                                                                       CORR 580
      SUBROUTINE CORRE (N,M,IO,X,XBAR,STD,RX,R,B,D,T)                   CORR 590
      DIMENSION X(1),XBAR(1),STD(1),RX(1),R(1),B(1),D(1),T(1)           CORR 600
C                                                                       CORR 610
C        ...............................................................CORR 620
C                                                                       CORR 630
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE  CORR 640
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION      CORR 650
C        STATEMENT WHICH FOLLOWS.                                       CORR 660
C                                                                       CORR 670
C     DOUBLE PRECISION XBAR,STD,RX,R,B,T                                CORR 680
C                                                                       CORR 690
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS    CORR 700
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS      CORR 710
C        ROUTINE.                                                       CORR 720
C                                                                       CORR 730
C        THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO      CORR 740
C        CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  SQRT AND ABS IN   CORR 750
C        STATEMENT 220 MUST BE CHANGED TO DSQRT AND DABS.               CORR 760
C                                                                       CORR 770
C        ...............................................................CORR 780
C                                                                       CORR 790
C     INITIALIZATION                                                    CORR 800
C                                                                       CORR 810
      DO 100 J=1,M                                                      CORR 820
      B(J)=0.0                                                          CORR 830
  100 T(J)=0.0                                                          CORR 840
      K=(M*M+M)/2                                                       CORR 850
      DO 102 I=1,K                                                      CORR 860
  102 R(I)=0.0                                                          CORR 870
      FN=N                                                              CORR 880
      L=0                                                               CORR 890
C                                                                       CORR 900
      IF(IO) 105, 127, 105                                              CORR 910
C                                                                       CORR 920
C     DATA ARE ALREADY IN CORE                                          CORR 930
C                                                                       CORR 940
  105 DO 108 J=1,M                                                      CORR 950
      DO 107 I=1,N                                                      CORR 960
      L=L+1                                                             CORR 970
  107 T(J)=T(J)+X(L)                                                    CORR 980
      XBAR(J)=T(J)                                                      CORR 990
  108 T(J)=T(J)/FN                                                      CORR1000
C                                                                       CORR1010
      DO 115 I=1,N                                                      CORR1020
      JK=0                                                              CORR1030
      L=I-N                                                             CORR1040
      DO 110 J=1,M                                                      CORR1050
      L=L+N                                                             CORR1060
      D(J)=X(L)-T(J)                                                    CORR1070
  110 B(J)=B(J)+D(J)                                                    CORR1080
      DO 115 J=1,M                                                      CORR1090
      DO 115 K=1,J                                                      CORR1100
      JK=JK+1                                                           CORR1110
  115 R(JK)=R(JK)+D(J)*D(K)                                             CORR1120
      GO TO 205                                                         CORR1130
C                                                                       CORR1140
C     READ OBSERVATIONS AND CALCULATE TEMPORARY                         CORR1150
C     MEANS FROM THESE DATA IN T(J)                                     CORR1160
C                                                                       CORR1170
  127 IF(N-M) 130, 130, 135                                             CORR1180
  130 KK=N                                                              CORR1190
      GO TO 137                                                         CORR1200
  135 KK=M                                                              CORR1210
  137 DO 140 I=1,KK                                                     CORR1220
      CALL DATA (M,D)                                                   CORR1230
      DO 140 J=1,M                                                      CORR1240
      T(J)=T(J)+D(J)                                                    CORR1250
      L=L+1                                                             CORR1260
  140 RX(L)=D(J)                                                        CORR1270
      FKK=KK                                                            CORR1280
      DO 150 J=1,M                                                      CORR1290
      XBAR(J)=T(J)                                                      CORR1300
  150 T(J)=T(J)/FKK                                                     CORR1310
C                                                                       CORR1320
C     CALCULATE SUMS OF CROSS-PRODUCTS OF DEVIATIONS                    CORR1330
C     FROM TEMPORARY MEANS FOR M OBSERVATIONS                           CORR1340
C                                                                       CORR1350
      L=0                                                               CORR1360
      DO 180 I=1,KK                                                     CORR1370
      JK=0                                                              CORR1380
      DO 170 J=1,M                                                      CORR1390
      L=L+1                                                             CORR1400
  170 D(J)=RX(L)-T(J)                                                   CORR1410
      DO 180 J=1,M                                                      CORR1420
      B(J)=B(J)+D(J)                                                    CORR1430
      DO 180 K=1,J                                                      CORR1440
      JK=JK+1                                                           CORR1450
  180 R(JK)=R(JK)+D(J)*D(K)                                             CORR1460
C                                                                       CORR1470
      IF(N-KK) 205, 205, 185                                            CORR1480
C                                                                       CORR1490
C     READ THE REST OF OBSERVATIONS ONE AT A TIME, SUM                  CORR1500
C     THE OBSERVATION, AND CALCULATE SUMS OF CROSS-                     CORR1510
C     PRODUCTS OF DEVIATIONS FROM TEMPORARY MEANS                       CORR1520
C                                                                       CORR1530
  185 KK=N-KK                                                           CORR1540
      DO 200 I=1,KK                                                     CORR1550
      JK=0                                                              CORR1560
      CALL DATA (M,D)                                                   CORR1570
      DO 190 J=1,M                                                      CORR1580
      XBAR(J)=XBAR(J)+D(J)                                              CORR1590
      D(J)=D(J)-T(J)                                                    CORR1600
  190 B(J)=B(J)+D(J)                                                    CORR1610
      DO 200 J=1,M                                                      CORR1620
      DO 200 K=1,J                                                      CORR1630
      JK=JK+1                                                           CORR1640
  200 R(JK)=R(JK)+D(J)*D(K)                                             CORR1650
C                                                                       CORR1660
C     CALCULATE MEANS                                                   CORR1670
C                                                                       CORR1680
  205 JK=0                                                              CORR1690
      DO 210 J=1,M                                                      CORR1700
      XBAR(J)=XBAR(J)/FN                                                CORR1710
C                                                                       CORR1720
C     ADJUST SUMS OF CROSS-PRODUCTS OF DEVIATIONS                       CORR1730
C     FROM TEMPORARY MEANS                                              CORR1740
C                                                                       CORR1750
      DO 210 K=1,J                                                      CORR1760
      JK=JK+1                                                           CORR1770
  210 R(JK)=R(JK)-B(J)*B(K)/FN                                          CORR1780
C                                                                       CORR1790
C     CALCULATE CORRELATION COEFFICIENTS                                CORR1800
C                                                                       CORR1810
      JK=0                                                              CORR1820
      DO 220 J=1,M                                                      CORR1830
      JK=JK+J                                                           CORR1840
  220 STD(J)= SQRT( ABS(R(JK)))                                         CORR1850
      DO 230 J=1,M                                                      CORR1860
      DO 230 K=J,M                                                      CORR1870
      JK=J+(K*K-K)/2                                                    CORR1880
      L=M*(J-1)+K                                                       CORR1890
      RX(L)=R(JK)                                                       CORR1900
      L=M*(K-1)+J                                                       CORR1910
      RX(L)=R(JK)                                                       CORR1920
      IF(STD(J)*STD(K)) 225, 222, 225                                   CORR1930
  222 R(JK)=0.0                                                         CORR1940
      GO TO 230                                                         CORR1950
  225 R(JK)=R(JK)/(STD(J)*STD(K))                                       CORR1960
  230 CONTINUE                                                          CORR1970
C                                                                       CORR1980
C     CALCULATE STANDARD DEVIATIONS                                     CORR1990
C                                                                       CORR2000
      FN=SQRT(FN-1.0)                                                   CORR2010
      DO 240 J=1,M                                                      CORR2020
  240 STD(J)=STD(J)/FN                                                  CORR2030
C                                                                       CORR2040
C     COPY THE DIAGONAL OF THE MATRIX OF SUMS OF CROSS-PRODUCTS OF      CORR2050
C     DEVIATIONS FROM MEANS.                                            CORR2060
C                                                                       CORR2070
      L=-M                                                              CORR2080
      DO 250 I=1,M                                                      CORR2090
      L=L+M+1                                                           CORR2100
  250 B(I)=RX(L)                                                        CORR2110
      RETURN                                                            CORR2120
      END                                                               CORR2130
