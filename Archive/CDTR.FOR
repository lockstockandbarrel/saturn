C                                                                       CDTR  10
C     ..................................................................CDTR  20
C                                                                       CDTR  30
C        SUBROUTINE CDTR                                                CDTR  40
C                                                                       CDTR  50
C        PURPOSE                                                        CDTR  60
C           COMPUTES P(X) = PROBABILITY THAT THE RANDOM VARIABLE U,     CDTR  70
C           DISTRIBUTED ACCORDING TO THE CHI-SQUARE DISTRIBUTION WITH G CDTR  80
C           DEGREES OF FREEDOM, IS LESS THAN OR EQUAL TO X.  F(G,X), THECDTR  90
C           ORDINATE OF THE CHI-SQUARE DENSITY AT X, IS ALSO COMPUTED.  CDTR 100
C                                                                       CDTR 110
C        USAGE                                                          CDTR 120
C           CALL CDTR(X,G,P,D,IER)                                      CDTR 130
C                                                                       CDTR 140
C        DESCRIPTION OF PARAMETERS                                      CDTR 150
C           X   - INPUT SCALAR FOR WHICH P(X) IS COMPUTED.              CDTR 160
C           G   - NUMBER OF DEGREES OF FREEDOM OF THE CHI-SQUARE        CDTR 170
C                 DISTRIBUTION.  G IS A CONTINUOUS PARAMETER.           CDTR 180
C           P   - OUTPUT PROBABILITY.                                   CDTR 190
C           D   - OUTPUT DENSITY.                                       CDTR 200
C           IER - RESULTANT ERROR CODE WHERE                            CDTR 210
C               IER= 0 --- NO ERROR                                     CDTR 220
C               IER=-1 --- AN INPUT PARAMETER IS INVALID.  X IS LESS    CDTR 230
C                          THAN 0.0, OR G IS LESS THAN 0.5 OR GREATER   CDTR 240
C                          THAN 2*10**(+5).  P AND D ARE SET TO -1.E75. CDTR 250
C               IER=+1 --- INVALID OUTPUT.  P IS LESS THAN ZERO OR      CDTR 260
C                          GREATER THAN ONE, OR SERIES FOR T1 (SEE      CDTR 270
C                          MATHEMATICAL DESCRIPTION) HAS FAILED TO      CDTR 280
C                          CONVERGE.  P IS SET TO 1.E75.                CDTR 290
C                                                                       CDTR 300
C        REMARKS                                                        CDTR 310
C           SEE MATHEMATICAL DESCRIPTION.                               CDTR 320
C                                                                       CDTR 330
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  CDTR 340
C           DLGAM                                                       CDTR 350
C           NDTR                                                        CDTR 360
C                                                                       CDTR 370
C        METHOD                                                         CDTR 380
C           REFER TO R.E. BARGMANN AND S.P. GHOSH, STATISTICAL          CDTR 390
C           DISTRIBUTION PROGRAMS FOR A COMPUTER LANGUAGE,              CDTR 400
C           IBM RESEARCH REPORT RC-1094, 1963.                          CDTR 410
C                                                                       CDTR 420
C     ..................................................................CDTR 430
C                                                                       CDTR 440
      SUBROUTINE CDTR(X,G,P,D,IER)                                      CDTR 450
      DOUBLE PRECISION XX,DLXX,X2,DLX2,GG,G2,DLT3,THETA,THP1,           CDTR 460
     1GLG2,DD,T11,SER,CC,XI,FAC,TLOG,TERM,GTH,A2,A,B,C,DT2,DT3,THPI     CDTR 470
C                                                                       CDTR 480
C        TEST FOR VALID INPUT DATA                                      CDTR 490
C                                                                       CDTR 500
      IF(G-(.5-1.E-5)) 590,10,10                                        CDTR 510
   10 IF(G-2.E+5) 20,20,590                                             CDTR 520
   20 IF(X) 590,30,30                                                   CDTR 530
C                                                                       CDTR 540
C        TEST FOR X NEAR 0.0                                            CDTR 550
C                                                                       CDTR 560
   30 IF(X-1.E-8) 40,40,80                                              CDTR 570
   40 P=0.0                                                             CDTR 580
      IF(G-2.) 50,60,70                                                 CDTR 590
   50 D=1.E75                                                           CDTR 600
      GO TO 610                                                         CDTR 610
   60 D=0.5                                                             CDTR 620
      GO TO 610                                                         CDTR 630
   70 D=0.0                                                             CDTR 640
      GO TO 610                                                         CDTR 650
C                                                                       CDTR 660
C        TEST FOR X GREATER THAN 1.E+6                                  CDTR 670
C                                                                       CDTR 680
   80 IF(X-1.E+6) 100,100,90                                            CDTR 690
   90 D=0.0                                                             CDTR 700
      P=1.0                                                             CDTR 710
      GO TO 610                                                         CDTR 720
C                                                                       CDTR 730
C        SET PROGRAM PARAMETERS                                         CDTR 740
C                                                                       CDTR 750
  100 XX=DBLE(X)                                                        CDTR 760
      DLXX=DLOG(XX)                                                     CDTR 770
      X2=XX/2.D0                                                        CDTR 780
      DLX2=DLOG(X2)                                                     CDTR 790
      GG=DBLE(G)                                                        CDTR 800
      G2=GG/2.D0                                                        CDTR 810
C                                                                       CDTR 820
C        COMPUTE ORDINATE                                               CDTR 830
C                                                                       CDTR 840
      CALL DLGAM(G2,GLG2,IOK)                                           CDTR 850
      DD=(G2-1.D0)*DLXX-X2-G2*.6931471805599453 -GLG2                   CDTR 860
      IF(DD-1.68D02) 110,110,120                                        CDTR 870
  110 IF(DD+1.68D02) 130,130,140                                        CDTR 880
  120 D=1.E75                                                           CDTR 890
      GO TO 150                                                         CDTR 900
  130 D=0.0                                                             CDTR 910
      GO TO 150                                                         CDTR 920
  140 DD=DEXP(DD)                                                       CDTR 930
      D=SNGL(DD)                                                        CDTR 940
C                                                                       CDTR 950
C        TEST FOR G GREATER THAN 1000.0                                 CDTR 960
C        TEST FOR X GREATER THAN 2000.0                                 CDTR 970
C                                                                       CDTR 980
  150 IF(G-1000.) 160,160,180                                           CDTR 990
  160 IF(X-2000.) 190,190,170                                           CDTR1000
  170 P=1.0                                                             CDTR1010
      GO TO 610                                                         CDTR1020
  180 A=DLOG(XX/GG)/3.D0                                                CDTR1030
      A=DEXP(A)                                                         CDTR1040
      B=2.D0/(9.D0*GG)                                                  CDTR1050
      C=(A-1.D0+B)/DSQRT(B)                                             CDTR1060
      SC=SNGL(C)                                                        CDTR1070
      CALL NDTR(SC,P,DUMMY)                                             CDTR1080
      GO TO 490                                                         CDTR1090
C                                                                       CDTR1100
C        COMPUTE THETA                                                  CDTR1110
C                                                                       CDTR1120
  190 K= IDINT(G2)                                                      CDTR1130
      THETA=G2-DFLOAT(K)                                                CDTR1140
      IF(THETA-1.D-8) 200,200,210                                       CDTR1150
  200 THETA=0.D0                                                        CDTR1160
  210 THP1=THETA+1.D0                                                   CDTR1170
C                                                                       CDTR1180
C        SELECT METHOD OF COMPUTING T1                                  CDTR1190
C                                                                       CDTR1200
      IF(THETA) 230,230,220                                             CDTR1210
  220 IF(XX-10.D0) 260,260,320                                          CDTR1220
C                                                                       CDTR1230
C        COMPUTE T1 FOR THETA EQUALS 0.0                                CDTR1240
C                                                                       CDTR1250
  230 IF(X2-1.68D02) 250,240,240                                        CDTR1260
  240 T1=1.0                                                            CDTR1270
      GO TO 400                                                         CDTR1280
  250 T11=1.D0-DEXP(-X2)                                                CDTR1290
      T1=SNGL(T11)                                                      CDTR1300
      GO TO 400                                                         CDTR1310
C                                                                       CDTR1320
C        COMPUTE T1 FOR THETA GREATER THAN 0.0 AND                      CDTR1330
C        X LESS THAN OR EQUAL TO 10.0                                   CDTR1340
C                                                                       CDTR1350
  260 SER=X2*(1.D0/THP1 -X2/(THP1+1.D0))                                CDTR1360
      J=+1                                                              CDTR1370
      CC=DFLOAT(J)                                                      CDTR1380
      DO 270 IT1=3,30                                                   CDTR1390
      XI=DFLOAT(IT1)                                                    CDTR1400
      CALL DLGAM(XI,FAC,IOK)                                            CDTR1410
      TLOG= XI*DLX2-FAC-DLOG(XI+THETA)                                  CDTR1420
      TERM=DEXP(TLOG)                                                   CDTR1430
      TERM=DSIGN(TERM,CC)                                               CDTR1440
      SER=SER+TERM                                                      CDTR1450
      CC=-CC                                                            CDTR1460
      IF(DABS(TERM)-1.D-9) 280,270,270                                  CDTR1470
  270 CONTINUE                                                          CDTR1480
      GO TO 600                                                         CDTR1490
  280 IF(SER) 600,600,290                                               CDTR1500
  290 CALL DLGAM(THP1,GTH,IOK)                                          CDTR1510
      TLOG=THETA*DLX2+DLOG(SER)-GTH                                     CDTR1520
      IF(TLOG+1.68D02) 300,300,310                                      CDTR1530
  300 T1=0.0                                                            CDTR1540
      GO TO 400                                                         CDTR1550
  310 T11=DEXP(TLOG)                                                    CDTR1560
      T1=SNGL(T11)                                                      CDTR1570
      GO TO 400                                                         CDTR1580
C                                                                       CDTR1590
C        COMPUTE T1 FOR THETA GREATER THAN 0.0 AND                      CDTR1600
C        X GREATER THAN 10.0 AND LESS THAN 2000.0                       CDTR1610
C                                                                       CDTR1620
  320 A2=0.D0                                                           CDTR1630
      DO 340 I=1,25                                                     CDTR1640
      XI=DFLOAT(I)                                                      CDTR1650
      CALL DLGAM(THP1,GTH,IOK)                                          CDTR1660
      T11=-(13.D0*XX)/XI +THP1*DLOG(13.D0*XX/XI) -GTH-DLOG(XI)          CDTR1670
      IF(T11+1.68D02) 340,340,330                                       CDTR1680
  330 T11=DEXP(T11)                                                     CDTR1690
      A2=A2+T11                                                         CDTR1700
  340 CONTINUE                                                          CDTR1710
      A=1.01282051+THETA/156.D0-XX/312.D0                               CDTR1720
      B=DABS(A)                                                         CDTR1730
      C= -X2+THP1*DLX2+DLOG(B)-GTH-3.951243718581427                    CDTR1740
      IF(C+1.68D02) 370,370,350                                         CDTR1750
  350 IF (A) 360,370,380                                                CDTR1760
  360 C=-DEXP(C)                                                        CDTR1770
      GO TO 390                                                         CDTR1780
  370 C=0.D0                                                            CDTR1790
      GO TO 390                                                         CDTR1800
  380 C=DEXP(C)                                                         CDTR1810
  390 C=A2+C                                                            CDTR1820
      T11=1.D0-C                                                        CDTR1830
      T1=SNGL(T11)                                                      CDTR1840
C                                                                       CDTR1850
C        SELECT PROPER EXPRESSION FOR P                                 CDTR1860
C                                                                       CDTR1870
  400 IF(G-2.) 420,410,410                                              CDTR1880
  410 IF(G-4.) 450,460,460                                              CDTR1890
C                                                                       CDTR1900
C        COMPUTE P FOR G GREATER THAN ZERO AND LESS THAN 2.0            CDTR1910
C                                                                       CDTR1920
  420 CALL DLGAM(THP1,GTH,IOK)                                          CDTR1930
      DT2=THETA*DLXX-X2-THP1*.6931471805599453 -GTH                     CDTR1940
      IF(DT2+1.68D02) 430,430,440                                       CDTR1950
  430 P=T1                                                              CDTR1960
      GO TO 490                                                         CDTR1970
  440 DT2=DEXP(DT2)                                                     CDTR1980
      T2=SNGL(DT2)                                                      CDTR1990
      P=T1+T2+T2                                                        CDTR2000
      GO TO 490                                                         CDTR2010
C                                                                       CDTR2020
C        COMPUTE P FOR G GREATER THAN OR EQUAL TO 2.0                   CDTR2030
C        AND LESS THAN 4.0                                              CDTR2040
C                                                                       CDTR2050
  450 P=T1                                                              CDTR2060
      GO TO 490                                                         CDTR2070
C                                                                       CDTR2080
C        COMPUTE P FOR G GREATER THAN OR EQUAL TO 4.0                   CDTR2090
C        AND LESS THAN OR EQUAL TO 1000.0                               CDTR2100
C                                                                       CDTR2110
  460 DT3=0.D0                                                          CDTR2120
      DO 480 I3=2,K                                                     CDTR2130
      THPI=DFLOAT(I3)+THETA                                             CDTR2140
      CALL DLGAM(THPI,GTH,IOK)                                          CDTR2150
      DLT3=THPI*DLX2-DLXX-X2-GTH                                        CDTR2160
      IF(DLT3+1.68D02) 480,480,470                                      CDTR2170
  470 DT3=DT3+DEXP(DLT3)                                                CDTR2180
  480 CONTINUE                                                          CDTR2190
      T3=SNGL(DT3)                                                      CDTR2200
      P=T1-T3-T3                                                        CDTR2210
C                                                                       CDTR2220
C        SET ERROR INDICATOR                                            CDTR2230
C                                                                       CDTR2240
  490 IF(P) 500,520,520                                                 CDTR2250
  500 IF(ABS(P)-1.E-7) 510,510,600                                      CDTR2260
  510 P=0.0                                                             CDTR2270
      GO TO 610                                                         CDTR2280
  520 IF(1.-P) 530,550,550                                              CDTR2290
  530 IF(ABS(1.-P)-1.E-7) 540,540,600                                   CDTR2300
  540 P=1.0                                                             CDTR2310
      GO TO 610                                                         CDTR2320
  550 IF(P-1.E-8) 560,560,570                                           CDTR2330
  560 P=0.0                                                             CDTR2340
      GO TO 610                                                         CDTR2350
  570 IF((1.0-P)-1.E-8) 580,580,610                                     CDTR2360
  580 P=1.0                                                             CDTR2370
      GO TO 610                                                         CDTR2380
  590 IER=-1                                                            CDTR2390
      D=-1.E75                                                          CDTR2400
      P=-1.E75                                                          CDTR2410
      GO TO 620                                                         CDTR2420
  600 IER=+1                                                            CDTR2430
      P= 1.E75                                                          CDTR2440
      GO TO 620                                                         CDTR2450
  610 IER=0                                                             CDTR2460
  620 RETURN                                                            CDTR2470
      END                                                               CDTR2480
