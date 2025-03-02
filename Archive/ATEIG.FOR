C                                                                       ATEI  10
C     ..................................................................ATEI  20
C                                                                       ATEI  30
C        SUBROUTINE ATEIG                                               ATEI  40
C                                                                       ATEI  50
C        PURPOSE                                                        ATEI  60
C           COMPUTE THE EIGENVALUES OF A REAL ALMOST TRIANGULAR MATRIX  ATEI  70
C                                                                       ATEI  80
C        USAGE                                                          ATEI  90
C           CALL ATEIG(M,A,RR,RI,IANA,IA)                               ATEI 100
C                                                                       ATEI 110
C        DESCRIPTION OF THE PARAMETERS                                  ATEI 120
C           M      ORDER OF THE MATRIX                                  ATEI 130
C           A      THE INPUT MATRIX, M BY M                             ATEI 140
C           RR     VECTOR CONTAINING THE REAL PARTS OF THE EIGENVALUES  ATEI 150
C                  ON RETURN                                            ATEI 160
C           RI     VECTOR CONTAINING THE IMAGINARY PARTS OF THE EIGEN-  ATEI 170
C                  VALUES ON RETURN                                     ATEI 180
C           IANA   VECTOR WHOSE DIMENSION MUST BE GREATER THAN OR EQUAL ATEI 190
C                  TO M, CONTAINING ON RETURN INDICATIONS ABOUT THE WAY ATEI 200
C                  THE EIGENVALUES APPEARED (SEE MATH. DESCRIPTION)     ATEI 210
C           IA     SIZE OF THE FIRST DIMENSION ASSIGNED TO THE ARRAY A  ATEI 220
C                  IN THE CALLING PROGRAM WHEN THE MATRIX IS IN DOUBLE  ATEI 230
C                  SUBSCRIPTED DATA STORAGE MODE.                       ATEI 240
C                  IA=M WHEN THE MATRIX IS IN SSP VECTOR STORAGE MODE.  ATEI 250
C                                                                       ATEI 260
C        REMARKS                                                        ATEI 270
C           THE ORIGINAL MATRIX IS DESTROYED                            ATEI 280
C           THE DIMENSION OF RR AND RI MUST BE GREATER OR EQUAL TO M    ATEI 290
C                                                                       ATEI 300
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  ATEI 310
C           NONE                                                        ATEI 320
C                                                                       ATEI 330
C        METHOD                                                         ATEI 340
C           QR DOUBLE ITERATION                                         ATEI 350
C                                                                       ATEI 360
C        REFERENCES                                                     ATEI 370
C           J.G.F. FRANCIS - THE QR TRANSFORMATION---THE COMPUTER       ATEI 380
C           JOURNAL, VOL. 4, NO. 3, OCTOBER 1961, VOL. 4, NO. 4, JANUARYATEI 390
C           1962.  J. H. WILKINSON - THE ALGEBRAIC EIGENVALUE PROBLEM - ATEI 400
C           CLARENDON PRESS, OXFORD, 1965.                              ATEI 410
C                                                                       ATEI 420
C     ..................................................................ATEI 430
C                                                                       ATEI 440
      SUBROUTINE ATEIG(M,A,RR,RI,IANA,IA)                               ATEI 450
      DIMENSION A(1),RR(1),RI(1),PRR(2),PRI(2),IANA(1)                  ATEI 460
      INTEGER P,P1,Q                                                    ATEI 470
C                                                                       ATEI 480
      E7=1.0E-8                                                         ATEI 490
      E6=1.0E-6                                                         ATEI 500
      E10=1.0E-10                                                       ATEI 510
      DELTA=0.5                                                         ATEI 520
      MAXIT=30                                                          ATEI 530
C                                                                       ATEI 540
C        INITIALIZATION                                                 ATEI 550
C                                                                       ATEI 560
      N=M                                                               ATEI 570
   20 N1=N-1                                                            ATEI 580
      IN=N1*IA                                                          ATEI 590
      NN=IN+N                                                           ATEI 600
      IF(N1) 30,1300,30                                                 ATEI 610
   30 NP=N+1                                                            ATEI 620
C                                                                       ATEI 630
C        ITERATION COUNTER                                              ATEI 640
C                                                                       ATEI 650
      IT=0                                                              ATEI 660
C                                                                       ATEI 670
C        ROOTS OF THE 2ND ORDER MAIN SUBMATRIX AT THE PREVIOUS          ATEI 680
C        ITERATION                                                      ATEI 690
C                                                                       ATEI 700
      DO 40 I=1,2                                                       ATEI 710
      PRR(I)=0.0                                                        ATEI 720
   40 PRI(I)=0.0                                                        ATEI 730
C                                                                       ATEI 740
C        LAST TWO SUBDIAGONAL ELEMENTS AT THE PREVIOUS ITERATION        ATEI 750
C                                                                       ATEI 760
      PAN=0.0                                                           ATEI 770
      PAN1=0.0                                                          ATEI 780
C                                                                       ATEI 790
C        ORIGIN SHIFT                                                   ATEI 800
C                                                                       ATEI 810
      R=0.0                                                             ATEI 820
      S=0.0                                                             ATEI 830
C                                                                       ATEI 840
C        ROOTS OF THE LOWER MAIN 2 BY 2 SUBMATRIX                       ATEI 850
C                                                                       ATEI 860
      N2=N1-1                                                           ATEI 870
      IN1=IN-IA                                                         ATEI 880
      NN1=IN1+N                                                         ATEI 890
      N1N=IN+N1                                                         ATEI 900
      N1N1=IN1+N1                                                       ATEI 910
   60 T=A(N1N1)-A(NN)                                                   ATEI 920
      U=T*T                                                             ATEI 930
      V=4.0*A(N1N)*A(NN1)                                               ATEI 940
      IF(ABS(V)-U*E7) 100,100,65                                        ATEI 950
   65 T=U+V                                                             ATEI 960
      IF(ABS(T)-AMAX1(U,ABS(V))*E6) 67,67,68                            ATEI 970
   67 T=0.0                                                             ATEI 980
   68 U=(A(N1N1)+A(NN))/2.0                                             ATEI 990
      V=SQRT(ABS(T))/2.0                                                ATEI1000
      IF(T)140,70,70                                                    ATEI1010
   70 IF(U) 80,75,75                                                    ATEI1020
   75 RR(N1)=U+V                                                        ATEI1030
      RR(N)=U-V                                                         ATEI1040
      GO TO 130                                                         ATEI1050
   80 RR(N1)=U-V                                                        ATEI1060
      RR(N)=U+V                                                         ATEI1070
      GO TO 130                                                         ATEI1080
  100 IF(T)120,110,110                                                  ATEI1090
  110 RR(N1)=A(N1N1)                                                    ATEI1100
      RR(N)=A(NN)                                                       ATEI1110
      GO TO 130                                                         ATEI1120
  120 RR(N1)=A(NN)                                                      ATEI1130
      RR(N)=A(N1N1)                                                     ATEI1140
  130 RI(N)=0.0                                                         ATEI1150
      RI(N1)=0.0                                                        ATEI1160
      GO TO 160                                                         ATEI1170
  140 RR(N1)=U                                                          ATEI1180
      RR(N)=U                                                           ATEI1190
      RI(N1)=V                                                          ATEI1200
      RI(N)=-V                                                          ATEI1210
  160 IF(N2)1280,1280,180                                               ATEI1220
C                                                                       ATEI1230
C        TESTS OF CONVERGENCE                                           ATEI1240
C                                                                       ATEI1250
  180 N1N2=N1N1-IA                                                      ATEI1260
      RMOD=RR(N1)*RR(N1)+RI(N1)*RI(N1)                                  ATEI1270
      EPS=E10*SQRT(RMOD)                                                ATEI1280
      IF(ABS(A(N1N2))-EPS)1280,1280,240                                 ATEI1290
  240 IF(ABS(A(NN1))-E10*ABS(A(NN))) 1300,1300,250                      ATEI1300
  250 IF(ABS(PAN1-A(N1N2))-ABS(A(N1N2))*E6) 1240,1240,260               ATEI1310
  260 IF(ABS(PAN-A(NN1))-ABS(A(NN1))*E6)1240,1240,300                   ATEI1320
  300 IF(IT-MAXIT) 320,1240,1240                                        ATEI1330
C                                                                       ATEI1340
C        COMPUTE THE SHIFT                                              ATEI1350
C                                                                       ATEI1360
  320 J=1                                                               ATEI1370
      DO 360 I=1,2                                                      ATEI1380
      K=NP-I                                                            ATEI1390
      IF(ABS(RR(K)-PRR(I))+ABS(RI(K)-PRI(I))-DELTA*(ABS(RR(K))          ATEI1400
     1    +ABS(RI(K)))) 340,360,360                                     ATEI1410
  340 J=J+I                                                             ATEI1420
  360 CONTINUE                                                          ATEI1430
      GO TO (440,460,460,480),J                                         ATEI1440
  440 R=0.0                                                             ATEI1450
      S=0.0                                                             ATEI1460
      GO TO 500                                                         ATEI1470
  460 J=N+2-J                                                           ATEI1480
      R=RR(J)*RR(J)                                                     ATEI1490
      S=RR(J)+RR(J)                                                     ATEI1500
      GO TO 500                                                         ATEI1510
  480 R=RR(N)*RR(N1)-RI(N)*RI(N1)                                       ATEI1520
      S=RR(N)+RR(N1)                                                    ATEI1530
C                                                                       ATEI1540
C        SAVE THE LAST TWO SUBDIAGONAL TERMS AND THE ROOTS OF THE       ATEI1550
C        SUBMATRIX BEFORE ITERATION                                     ATEI1560
C                                                                       ATEI1570
  500 PAN=A(NN1)                                                        ATEI1580
      PAN1=A(N1N2)                                                      ATEI1590
      DO 520 I=1,2                                                      ATEI1600
      K=NP-I                                                            ATEI1610
      PRR(I)=RR(K)                                                      ATEI1620
  520 PRI(I)=RI(K)                                                      ATEI1630
C                                                                       ATEI1640
C        SEARCH FOR A PARTITION OF THE MATRIX, DEFINED BY P AND Q       ATEI1650
C                                                                       ATEI1660
      P=N2                                                              ATEI1670
      IF (N-3)600,600,525                                               ATEI1675
  525 IPI=N1N2                                                          ATEI1680
      DO 580 J=2,N2                                                     ATEI1690
      IPI=IPI-IA-1                                                      ATEI1700
      IF(ABS(A(IPI))-EPS) 600,600,530                                   ATEI1710
  530 IPIP=IPI+IA                                                       ATEI1720
      IPIP2=IPIP+IA                                                     ATEI1730
      D=A(IPIP)*(A(IPIP)-S)+A(IPIP2)*A(IPIP+1)+R                        ATEI1740
      IF(D)540,560,540                                                  ATEI1750
  540 IF(ABS(A(IPI)*A(IPIP+1))*(ABS(A(IPIP)+A(IPIP2+1)-S)+ABS(A(IPIP2+2)ATEI1760
     1 )) -ABS(D)*EPS) 620,620,560                                      ATEI1770
  560 P=N1-J                                                            ATEI1780
  580 CONTINUE                                                          ATEI1790
  600 Q=P                                                               ATEI1800
      GO TO 680                                                         ATEI1810
  620 P1=P-1                                                            ATEI1820
      Q=P1                                                              ATEI1830
      IF (P1-1) 680,680,650                                             ATEI1835
  650 DO 660 I=2, P1                                                    ATEI1840
      IPI=IPI-IA-1                                                      ATEI1850
      IF(ABS(A(IPI))-EPS)680,680,660                                    ATEI1860
  660 Q=Q-1                                                             ATEI1870
C                                                                       ATEI1880
C        QR DOUBLE ITERATION                                            ATEI1890
C                                                                       ATEI1900
  680 II=(P-1)*IA+P                                                     ATEI1910
      DO 1220 I=P,N1                                                    ATEI1920
      II1=II-IA                                                         ATEI1930
      IIP=II+IA                                                         ATEI1940
      IF(I-P)720,700,720                                                ATEI1950
  700 IPI=II+1                                                          ATEI1960
      IPIP=IIP+1                                                        ATEI1970
C                                                                       ATEI1980
C        INITIALIZATION OF THE TRANSFORMATION                           ATEI1990
C                                                                       ATEI2000
      G1=A(II)*(A(II)-S)+A(IIP)*A(IPI)+R                                ATEI2010
      G2=A(IPI)*(A(IPIP)+A(II)-S)                                       ATEI2020
      G3=A(IPI)*A(IPIP+1)                                               ATEI2030
      A(IPI+1)=0.0                                                      ATEI2040
      GO TO 780                                                         ATEI2050
  720 G1=A(II1)                                                         ATEI2060
      G2=A(II1+1)                                                       ATEI2070
      IF(I-N2)740,740,760                                               ATEI2080
  740 G3=A(II1+2)                                                       ATEI2090
      GO TO 780                                                         ATEI2100
  760 G3=0.0                                                            ATEI2110
  780 CAP=SQRT(G1*G1+G2*G2+G3*G3)                                       ATEI2120
      IF(CAP)800,860,800                                                ATEI2130
  800 IF(G1)820,840,840                                                 ATEI2140
  820 CAP=-CAP                                                          ATEI2150
  840 T=G1+CAP                                                          ATEI2160
      PSI1=G2/T                                                         ATEI2170
      PSI2=G3/T                                                         ATEI2180
      ALPHA=2.0/(1.0+PSI1*PSI1+PSI2*PSI2)                               ATEI2190
      GO TO 880                                                         ATEI2200
  860 ALPHA=2.0                                                         ATEI2210
      PSI1=0.0                                                          ATEI2220
      PSI2=0.0                                                          ATEI2230
  880 IF(I-Q)900,960,900                                                ATEI2240
  900 IF(I-P)920,940,920                                                ATEI2250
  920 A(II1)=-CAP                                                       ATEI2260
      GO TO 960                                                         ATEI2270
  940 A(II1)=-A(II1)                                                    ATEI2280
C                                                                       ATEI2290
C        ROW OPERATION                                                  ATEI2300
C                                                                       ATEI2310
  960 IJ=II                                                             ATEI2320
      DO 1040 J=I,N                                                     ATEI2330
      T=PSI1*A(IJ+1)                                                    ATEI2340
      IF(I-N1)980,1000,1000                                             ATEI2350
  980 IP2J=IJ+2                                                         ATEI2360
      T=T+PSI2*A(IP2J)                                                  ATEI2370
 1000 ETA=ALPHA*(T+A(IJ))                                               ATEI2380
      A(IJ)=A(IJ)-ETA                                                   ATEI2390
      A(IJ+1)=A(IJ+1)-PSI1*ETA                                          ATEI2400
      IF(I-N1)1020,1040,1040                                            ATEI2410
 1020 A(IP2J)=A(IP2J)-PSI2*ETA                                          ATEI2420
 1040 IJ=IJ+IA                                                          ATEI2430
C                                                                       ATEI2440
C        COLUMN OPERATION                                               ATEI2450
C                                                                       ATEI2460
      IF(I-N1)1080,1060,1060                                            ATEI2470
 1060 K=N                                                               ATEI2480
      GO TO 1100                                                        ATEI2490
 1080 K=I+2                                                             ATEI2500
 1100 IP=IIP-I                                                          ATEI2510
      DO 1180 J=Q,K                                                     ATEI2520
      JIP=IP+J                                                          ATEI2530
      JI=JIP-IA                                                         ATEI2540
      T=PSI1*A(JIP)                                                     ATEI2550
      IF(I-N1)1120,1140,1140                                            ATEI2560
 1120 JIP2=JIP+IA                                                       ATEI2570
      T=T+PSI2*A(JIP2)                                                  ATEI2580
 1140 ETA=ALPHA*(T+A(JI))                                               ATEI2590
      A(JI)=A(JI)-ETA                                                   ATEI2600
      A(JIP)=A(JIP)-ETA*PSI1                                            ATEI2610
      IF(I-N1)1160,1180,1180                                            ATEI2620
 1160 A(JIP2)=A(JIP2)-ETA*PSI2                                          ATEI2630
 1180 CONTINUE                                                          ATEI2640
      IF(I-N2)1200,1220,1220                                            ATEI2650
 1200 JI=II+3                                                           ATEI2660
      JIP=JI+IA                                                         ATEI2670
      JIP2=JIP+IA                                                       ATEI2680
      ETA=ALPHA*PSI2*A(JIP2)                                            ATEI2690
      A(JI)=-ETA                                                        ATEI2700
      A(JIP)=-ETA*PSI1                                                  ATEI2710
      A(JIP2)=A(JIP2)-ETA*PSI2                                          ATEI2720
 1220 II=IIP+1                                                          ATEI2730
      IT=IT+1                                                           ATEI2740
      GO TO 60                                                          ATEI2750
C                                                                       ATEI2760
C        END OF ITERATION                                               ATEI2770
C                                                                       ATEI2780
 1240 IF(ABS(A(NN1))-ABS(A(N1N2))) 1300,1280,1280                       ATEI2790
C                                                                       ATEI2800
C        TWO EIGENVALUES HAVE BEEN FOUND                                ATEI2810
C                                                                       ATEI2820
 1280 IANA(N)=0                                                         ATEI2830
      IANA(N1)=2                                                        ATEI2840
      N=N2                                                              ATEI2850
      IF(N2)1400,1400,20                                                ATEI2860
C                                                                       ATEI2870
C        ONE EIGENVALUE HAS BEEN FOUND                                  ATEI2880
C                                                                       ATEI2890
 1300 RR(N)=A(NN)                                                       ATEI2900
      RI(N)=0.0                                                         ATEI2910
      IANA(N)=1                                                         ATEI2920
      IF(N1)1400,1400,1320                                              ATEI2930
 1320 N=N1                                                              ATEI2940
      GO TO 20                                                          ATEI2950
 1400 RETURN                                                            ATEI2960
      END                                                               ATEI2970
