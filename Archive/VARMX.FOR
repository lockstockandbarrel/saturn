C                                                                       VARM  10
C     ..................................................................VARM  20
C                                                                       VARM  30
C        SUBROUTINE VARMX                                               VARM  40
C                                                                       VARM  50
C        PURPOSE                                                        VARM  60
C           PERFORM ORTHOGONAL ROTATIONS OF A FACTOR MATRIX.  THIS      VARM  70
C           SUBROUTINE NORMALLY OCCURS IN A SEQUENCE OF CALLS TO SUB-   VARM  80
C           ROUTINES CORRE, EIGEN, TRACE, LOAD, VARMX IN THE PERFORMANCEVARM  90
C           OF A FACTOR ANALYSIS.                                       VARM 100
C                                                                       VARM 110
C        USAGE                                                          VARM 120
C           CALL VARMX (M,K,A,NC,TV,H,F,D,IER)                          VARM 130
C                                                                       VARM 140
C        DESCRIPTION OF PARAMETERS                                      VARM 150
C           M     - NUMBER OF VARIABLES AND NUMBER OF ROWS OF MATRIX A. VARM 160
C           K     - NUMBER OF FACTORS.                                  VARM 170
C           A     - INPUT IS THE ORIGINAL FACTOR MATRIX, AND OUTPUT IS  VARM 180
C                   THE ROTATED FACTOR MATRIX.  THE ORDER OF MATRIX A   VARM 190
C                   IS M X K.                                           VARM 200
C           NC    - OUTPUT VARIABLE CONTAINING THE NUMBER OF ITERATION  VARM 210
C                   CYCLES PERFORMED.                                   VARM 220
C           TV    - OUTPUT VECTOR CONTAINING THE VARIANCE OF THE FACTOR VARM 230
C                   MATRIX FOR EACH ITERATION CYCLE.  THE VARIANCE PRIORVARM 240
C                   TO THE FIRST ITERATION CYCLE IS ALSO CALCULATED.    VARM 250
C                   THIS MEANS THAT NC+1 VARIANCES ARE STORED IN VECTOR VARM 260
C                   TV.  MAXIMUM NUMBER OF ITERATION CYCLES ALLOWED IN  VARM 270
C                   THIS SUBROUTINE IS 50.  THEREFORE, THE LENGTH OF    VARM 280
C                   VECTOR TV IS 51.                                    VARM 290
C           H     - OUTPUT VECTOR OF LENGTH M CONTAINING THE ORIGINAL   VARM 300
C                   COMMUNALITIES.                                      VARM 310
C           F     - OUTPUT VECTOR OF LENGTH M CONTAINING THE FINAL      VARM 320
C                   COMMUNALITIES.                                      VARM 330
C           D     - OUTPUT VECTOR OF LENGTH M CONTAINING THE DIFFERENCESVARM 340
C                   BETWEEN THE ORIGINAL AND FINAL COMMUNALITIES.       VARM 350
C           IER   - ERROR INDICATOR                                     VARM 351
C                   IER=0 - NO ERROR                                    VARM 352
C                   IER=1 - CONVERGENCE WAS NOT ACHIEVED IN 50 CYCLES   VARM 353
C                           OF ROTATION                                 VARM 354
C                                                                       VARM 360
C        REMARKS                                                        VARM 370
C           IF VARIANCE COMPUTED AFTER EACH ITERATION CYCLE DOES NOT    VARM 380
C           INCREASE FOR FOUR SUCCESSIVE TIMES, THE SUBROUTINE STOPS    VARM 390
C           ROTATION.                                                   VARM 400
C                                                                       VARM 410
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  VARM 420
C           NONE                                                        VARM 430
C                                                                       VARM 440
C        METHOD                                                         VARM 450
C           KAISER'S VARIMAX ROTATION AS DESCRIBED IN 'COMPUTER PROGRAM VARM 460
C           FOR VARIMAX ROTATION IN FACTOR ANALYSIS' BY THE SAME AUTHOR,VARM 470
C           EDUCATIONAL AND PSYCHOLOGICAL MEASUREMENT, VOL XIX, NO. 3,  VARM 480
C           1959.                                                       VARM 490
C                                                                       VARM 500
C     ..................................................................VARM 510
C                                                                       VARM 520
      SUBROUTINE VARMX (M,K,A,NC,TV,H,F,D,IER)                          VARM 530
      DIMENSION A(1),TV(1),H(1),F(1),D(1)                               VARM 540
C                                                                       VARM 550
C        ...............................................................VARM 560
C                                                                       VARM 570
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE  VARM 580
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION      VARM 590
C                                                                       VARM 600
C     DOUBLE PRECISION A,TV,H,F,D,TVLT,CONS,AA,BB,CC,DD,U,T,B,COS4T,    VARM 610
C    1                 SIN4T,TAN4T,SINP,COSP,CTN4T,COS2T,SIN2T,COST,SINTVARM 620
C                                                                       VARM 630
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS    VARM 640
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS      VARM 650
C        ROUTINE.                                                       VARM 660
C                                                                       VARM 670
C        THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO      VARM 680
C        CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  SQRT IN STATEMENTSVARM 690
C        115, 290, 330, 350, AND 355 MUST BE CHANGED TO DSQRT.  ABS IN  VARM 700
C        STATEMENTS 280, 320, AND 375 MUST BE CHANGED TO DABS.          VARM 710
C                                                                       VARM 720
C        ...............................................................VARM 730
C                                                                       VARM 740
C     INITIALIZATION                                                    VARM 750
C                                                                       VARM 760
      IER=0                                                             VARM 761
      EPS=0.00116                                                       VARM 770
      TVLT=0.0                                                          VARM 780
      LL=K-1                                                            VARM 790
      NV=1                                                              VARM 800
      NC=0                                                              VARM 810
      FN=M                                                              VARM 820
      FFN=FN*FN                                                         VARM 830
      CONS=0.7071066                                                    VARM 840
C                                                                       VARM 850
C     CALCULATE ORIGINAL COMMUNALITIES                                  VARM 860
C                                                                       VARM 870
      DO 110 I=1,M                                                      VARM 880
      H(I)=0.0                                                          VARM 890
      DO 110 J=1,K                                                      VARM 900
      L=M*(J-1)+I                                                       VARM 910
  110 H(I)=H(I)+A(L)*A(L)                                               VARM 920
C                                                                       VARM 930
C     CALCULATE NORMALIZED FACTOR MATRIX                                VARM 940
C                                                                       VARM 950
      DO 120 I=1,M                                                      VARM 960
  115 H(I)= SQRT(H(I))                                                  VARM 970
      DO 120 J=1,K                                                      VARM 980
      L=M*(J-1)+I                                                       VARM 990
  120 A(L)=A(L)/H(I)                                                    VARM1000
      GO TO 132                                                         VARM1010
C                                                                       VARM1020
C     CALCULATE VARIANCE FOR FACTOR MATRIX                              VARM1030
C                                                                       VARM1040
  130 NV=NV+1                                                           VARM1050
      TVLT=TV(NV-1)                                                     VARM1060
  132 TV(NV)=0.0                                                        VARM1070
      DO 150 J=1,K                                                      VARM1080
      AA=0.0                                                            VARM1090
      BB=0.0                                                            VARM1100
      LB=M*(J-1)                                                        VARM1110
      DO 140 I=1,M                                                      VARM1120
      L=LB+I                                                            VARM1130
      CC=A(L)*A(L)                                                      VARM1140
      AA=AA+CC                                                          VARM1150
  140 BB=BB+CC*CC                                                       VARM1160
  150 TV(NV)=TV(NV)+(FN*BB-AA*AA)/FFN                                   VARM1170
      IF(NV-51)160,155,155                                              VARM1180
  155 IER=1                                                             VARM1181
      GO TO 430                                                         VARM1182
C                                                                       VARM1190
C     PERFORM CONVERGENCE TEST                                          VARM1200
C                                                                       VARM1210
  160 IF((TV(NV)-TVLT)-(1.E-7)) 170, 170, 190                           VARM1220
  170 NC=NC+1                                                           VARM1230
      IF(NC-3) 190, 190, 430                                            VARM1240
C                                                                       VARM1250
C     ROTATION OF TWO FACTORS CONTINUES UP TO                           VARM1260
C     THE STATEMENT 120.                                                VARM1270
C                                                                       VARM1280
  190 DO 420 J=1,LL                                                     VARM1290
      L1=M*(J-1)                                                        VARM1300
      II=J+1                                                            VARM1310
C                                                                       VARM1320
C        CALCULATE NUM AND DEN                                          VARM1330
C                                                                       VARM1340
      DO 420 K1=II,K                                                    VARM1350
      L2=M*(K1-1)                                                       VARM1360
      AA=0.0                                                            VARM1370
      BB=0.0                                                            VARM1380
      CC=0.0                                                            VARM1390
      DD=0.0                                                            VARM1400
      DO 230 I=1,M                                                      VARM1410
      L3=L1+I                                                           VARM1420
      L4=L2+I                                                           VARM1430
      U=(A(L3)+A(L4))*(A(L3)-A(L4))                                     VARM1440
      T=A(L3)*A(L4)                                                     VARM1450
      T=T+T                                                             VARM1460
      CC=CC+(U+T)*(U-T)                                                 VARM1470
      DD=DD+2.0*U*T                                                     VARM1480
      AA=AA+U                                                           VARM1490
  230 BB=BB+T                                                           VARM1500
      T=DD-2.0*AA*BB/FN                                                 VARM1510
      B=CC-(AA*AA-BB*BB)/FN                                             VARM1520
C                                                                       VARM1530
C        COMPARISON OF NUM AND DEN                                      VARM1540
C                                                                       VARM1550
      IF(T-B) 280, 240, 320                                             VARM1560
  240 IF((T+B)-EPS) 420, 250, 250                                       VARM1570
C                                                                       VARM1580
C        NUM + DEN IS GREATER THAN OR EQUAL TO THE                      VARM1590
C        TOLERANCE FACTOR                                               VARM1600
C                                                                       VARM1610
  250 COS4T=CONS                                                        VARM1620
      SIN4T=CONS                                                        VARM1630
      GO TO 350                                                         VARM1640
C                                                                       VARM1650
C        NUM IS LESS THAN DEN                                           VARM1660
C                                                                       VARM1670
  280 TAN4T= ABS(T)/ ABS(B)                                             VARM1680
      IF(TAN4T-EPS) 300, 290, 290                                       VARM1690
  290 COS4T=1.0/ SQRT(1.0+TAN4T*TAN4T)                                  VARM1700
      SIN4T=TAN4T*COS4T                                                 VARM1710
      GO TO 350                                                         VARM1720
  300 IF(B) 310, 420, 420                                               VARM1730
  310 SINP=CONS                                                         VARM1740
      COSP=CONS                                                         VARM1750
      GO TO 400                                                         VARM1760
C                                                                       VARM1770
C        NUM IS GREATER THAN DEN                                        VARM1780
C                                                                       VARM1790
  320 CTN4T= ABS(T/B)                                                   VARM1800
      IF(CTN4T-EPS) 340, 330, 330                                       VARM1810
  330 SIN4T=1.0/ SQRT(1.0+CTN4T*CTN4T)                                  VARM1820
      COS4T=CTN4T*SIN4T                                                 VARM1830
      GO TO 350                                                         VARM1840
  340 COS4T=0.0                                                         VARM1850
      SIN4T=1.0                                                         VARM1860
C                                                                       VARM1870
C        DETERMINE COS THETA AND SIN THETA                              VARM1880
C                                                                       VARM1890
  350 COS2T= SQRT((1.0+COS4T)/2.0)                                      VARM1900
      SIN2T=SIN4T/(2.0*COS2T)                                           VARM1910
  355 COST= SQRT((1.0+COS2T)/2.0)                                       VARM1920
      SINT=SIN2T/(2.0*COST)                                             VARM1930
C                                                                       VARM1940
C        DETERMINE COS PHI AND SIN PHI                                  VARM1950
C                                                                       VARM1960
      IF(B) 370, 370, 360                                               VARM1970
  360 COSP=COST                                                         VARM1980
      SINP=SINT                                                         VARM1990
      GO TO 380                                                         VARM2000
  370 COSP=CONS*COST+CONS*SINT                                          VARM2010
  375 SINP= ABS(CONS*COST-CONS*SINT)                                    VARM2020
  380 IF(T) 390, 390, 400                                               VARM2030
  390 SINP=-SINP                                                        VARM2040
C                                                                       VARM2050
C        PERFORM ROTATION                                               VARM2060
C                                                                       VARM2070
  400 DO 410 I=1,M                                                      VARM2080
      L3=L1+I                                                           VARM2090
      L4=L2+I                                                           VARM2100
      AA=A(L3)*COSP+A(L4)*SINP                                          VARM2110
      A(L4)=-A(L3)*SINP+A(L4)*COSP                                      VARM2120
  410 A(L3)=AA                                                          VARM2130
  420 CONTINUE                                                          VARM2140
      GO TO 130                                                         VARM2150
C                                                                       VARM2160
C     DENORMALIZE VARIMAX LOADINGS                                      VARM2170
C                                                                       VARM2180
  430 DO 440 I=1,M                                                      VARM2190
      DO 440 J=1,K                                                      VARM2200
      L=M*(J-1)+I                                                       VARM2210
  440 A(L)=A(L)*H(I)                                                    VARM2220
C                                                                       VARM2230
C     CHECK ON COMMUNALITIES                                            VARM2240
C                                                                       VARM2250
      NC=NV-1                                                           VARM2260
      DO 450 I=1,M                                                      VARM2270
  450 H(I)=H(I)*H(I)                                                    VARM2280
      DO 470 I=1,M                                                      VARM2290
      F(I)=0.0                                                          VARM2300
      DO 460 J=1,K                                                      VARM2310
      L=M*(J-1)+I                                                       VARM2320
  460 F(I)=F(I)+A(L)*A(L)                                               VARM2330
  470 D(I)=H(I)-F(I)                                                    VARM2340
      RETURN                                                            VARM2350
      END                                                               VARM2360
