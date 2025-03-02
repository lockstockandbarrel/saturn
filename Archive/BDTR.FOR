C                                                                       BDTR  10
C     ..................................................................BDTR  20
C                                                                       BDTR  30
C        SUBROUTINE BDTR                                                BDTR  40
C                                                                       BDTR  50
C        PURPOSE                                                        BDTR  60
C           COMPUTES P(X) = PROBABILITY THAT THE RANDOM VARIABLE U,     BDTR  70
C           DISTRIBUTED ACCORDING TO THE BETA DISTRIBUTION WITH         BDTR  80
C           PARAMETERS A AND B, IS LESS THAN OR EQUAL TO X.  F(A,B,X),  BDTR  90
C           THE ORDINATE OF THE BETA DENSITY AT X, IS ALSO COMPUTED.    BDTR 100
C                                                                       BDTR 110
C        USAGE                                                          BDTR 120
C           CALL BDTR(X,A,B,P,D,IER)                                    BDTR 130
C                                                                       BDTR 140
C        DESCRIPTION OF PARAMETERS                                      BDTR 150
C           X   - INPUT SCALAR FOR WHICH P(X) IS COMPUTED.              BDTR 160
C           A   - BETA DISTRIBUTION PARAMETER (CONTINUOUS).             BDTR 170
C           B   - BETA DISTRIBUTION PARAMETER (CONTINUOUS).             BDTR 180
C           P   - OUTPUT PROBABILITY.                                   BDTR 190
C           D   - OUTPUT DENSITY.                                       BDTR 200
C           IER - RESULTANT ERROR CODE WHERE                            BDTR 210
C               IER= 0 --- NO ERROR                                     BDTR 220
C               IER=-1,+1  CDTR HAS BEEN CALLED AND AN ERROR HAS        BDTR 230
C                          OCCURRED.  SEE CDTR.                         BDTR 240
C               IER=-2 --- AN INPUT PARAMETER IS INVALID.  X IS LESS    BDTR 250
C                          THAN 0.0 OR GREATER THAN 1.0, OR EITHER A OR BDTR 260
C                          B IS LESS THAN 0.5 OR GREATER THAN 10**(+5). BDTR 270
C                          P AND D ARE SET TO -1.E75.                   BDTR 280
C               IER=+2 --- INVALID OUTPUT.  P IS LESS THAN ZERO OR      BDTR 290
C                          GREATER THAN ONE.  P IS SET TO 1.E75.        BDTR 300
C                                                                       BDTR 310
C        REMARKS                                                        BDTR 320
C           SEE MATHEMATICAL DESCRIPTION.                               BDTR 330
C                                                                       BDTR 340
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  BDTR 350
C           DLGAM                                                       BDTR 360
C           NDTR                                                        BDTR 370
C           CDTR                                                        BDTR 380
C                                                                       BDTR 390
C        METHOD                                                         BDTR 400
C           REFER TO R.E. BARGMANN AND S.P. GHOSH, STATISTICAL          BDTR 410
C           DISTRIBUTION PROGRAMS FOR A COMPUTER LANGUAGE,              BDTR 420
C           IBM RESEARCH REPORT RC-1094, 1963.                          BDTR 430
C                                                                       BDTR 440
C     ..................................................................BDTR 450
C                                                                       BDTR 460
      SUBROUTINE BDTR(X,A,B,P,D,IER)                                    BDTR 470
      DOUBLE PRECISION XX,DLXX,DL1X,AA,BB,G1,G2,G3,G4,DD,PP,XO,FF,FN,   BDTR 480
     1XI,SS,CC,RR,DLBETA                                                BDTR 490
C                                                                       BDTR 500
C        TEST FOR VALID INPUT DATA                                      BDTR 510
C                                                                       BDTR 520
      IF(A-(.5-1.E-5)) 640,10,10                                        BDTR 530
   10 IF(B-(.5-1.E-5)) 640,20,20                                        BDTR 540
   20 IF(A-1.E+5) 30,30,640                                             BDTR 550
   30 IF(B-1.E+5) 40,40,640                                             BDTR 560
   40 IF(X) 640,50,50                                                   BDTR 570
   50 IF(1.-X) 640,60,60                                                BDTR 580
C                                                                       BDTR 590
C        COMPUTE LOG(BETA(A,B))                                         BDTR 600
C                                                                       BDTR 610
   60 AA=DBLE(A)                                                        BDTR 620
      BB=DBLE(B)                                                        BDTR 630
      CALL DLGAM(AA,G1,IOK)                                             BDTR 640
      CALL DLGAM(BB,G2,IOK)                                             BDTR 650
      CALL DLGAM(AA+BB,G3,IOK)                                          BDTR 660
      DLBETA=G1+G2-G3                                                   BDTR 670
C                                                                       BDTR 680
C        TEST FOR X NEAR 0.0 OR 1.0                                     BDTR 690
C                                                                       BDTR 700
      IF(X-1.E-8) 80,80,70                                              BDTR 710
   70 IF((1.-X)-1.E-8) 130,130,140                                      BDTR 720
   80 P=0.0                                                             BDTR 730
      IF(A-1.) 90,100,120                                               BDTR 740
   90 D=1.E+75                                                          BDTR 750
      GO TO 660                                                         BDTR 760
  100 DD=-DLBETA                                                        BDTR 770
      IF(DD+1.68D02)  120,120,110                                       BDTR 780
  110 DD=DEXP(DD)                                                       BDTR 790
      D=SNGL(DD)                                                        BDTR 800
      GO TO 660                                                         BDTR 810
  120 D=0.0                                                             BDTR 820
      GO TO 660                                                         BDTR 830
  130 P=1.0                                                             BDTR 840
      IF(B-1.) 90,100,120                                               BDTR 850
C                                                                       BDTR 860
C        SET PROGRAM PARAMETERS                                         BDTR 870
C                                                                       BDTR 880
  140 XX=DBLE(X)                                                        BDTR 890
      DLXX=DLOG(XX)                                                     BDTR 900
      DL1X=DLOG(1.D0-XX)                                                BDTR 910
      XO=XX/(1.D0-XX)                                                   BDTR 920
      ID=0                                                              BDTR 930
C                                                                       BDTR 940
C        COMPUTE ORDINATE                                               BDTR 950
C                                                                       BDTR 960
      DD=(AA-1.D0)*DLXX+(BB-1.D0)*DL1X-DLBETA                           BDTR 970
      IF(DD-1.68D02) 150,150,160                                        BDTR 980
  150 IF(DD+1.68D02) 170,170,180                                        BDTR 990
  160 D=1.E75                                                           BDTR1000
      GO TO 190                                                         BDTR1010
  170 D=0.0                                                             BDTR1020
      GO TO 190                                                         BDTR1030
  180 DD=DEXP(DD)                                                       BDTR1040
      D=SNGL(DD)                                                        BDTR1050
C                                                                       BDTR1060
C        A OR B OR BOTH WITHIN 1.E-8 OF 1.0                             BDTR1070
C                                                                       BDTR1080
  190 IF(ABS(A-1.)-1.E-8)  200,200,210                                  BDTR1090
  200 IF(ABS(B-1.)-1.E-8)  220,220,230                                  BDTR1100
  210 IF(ABS(B-1.)-1.E-8)  260,260,290                                  BDTR1110
  220 P=X                                                               BDTR1120
      GO TO 660                                                         BDTR1130
  230 PP=BB*DL1X                                                        BDTR1140
      IF(PP+1.68D02) 240,240,250                                        BDTR1150
  240 P=1.0                                                             BDTR1160
      GO TO 660                                                         BDTR1170
  250 PP=DEXP(PP)                                                       BDTR1180
      PP=1.D0-PP                                                        BDTR1190
      P=SNGL(PP)                                                        BDTR1200
      GO TO 600                                                         BDTR1210
  260 PP=AA*DLXX                                                        BDTR1220
      IF(PP+1.68D02) 270,270,280                                        BDTR1230
  270 P=0.0                                                             BDTR1240
      GO TO 660                                                         BDTR1250
  280 PP=DEXP(PP)                                                       BDTR1260
      P=SNGL(PP)                                                        BDTR1270
      GO TO 600                                                         BDTR1280
C                                                                       BDTR1290
C        TEST FOR A OR B GREATER THAN 1000.0                            BDTR1300
C                                                                       BDTR1310
  290 IF(A-1000.) 300,300,310                                           BDTR1320
  300 IF(B-1000.) 330,330,320                                           BDTR1330
  310 XX=2.D0*AA/XO                                                     BDTR1340
      XS=SNGL(XX)                                                       BDTR1350
      AA=2.D0*BB                                                        BDTR1360
      DF=SNGL(AA)                                                       BDTR1370
      CALL CDTR(XS,DF,P,DUMMY,IER)                                      BDTR1380
      P=1.0-P                                                           BDTR1390
      GO TO 670                                                         BDTR1400
  320 XX=2.D0*BB*XO                                                     BDTR1410
      XS=SNGL(XX)                                                       BDTR1420
      AA=2.D0*AA                                                        BDTR1430
      DF=SNGL(AA)                                                       BDTR1440
      CALL CDTR(XS,DF,P,DUMMY,IER)                                      BDTR1450
      GO TO 670                                                         BDTR1460
C                                                                       BDTR1470
C        SELECT PARAMETERS FOR CONTINUED FRACTION COMPUTATION           BDTR1480
C                                                                       BDTR1490
  330 IF(X-.5) 340,340,380                                              BDTR1500
  340 IF(AA-1.D0) 350,350,360                                           BDTR1510
  350 RR=AA+1.D0                                                        BDTR1520
      GO TO 370                                                         BDTR1530
  360 RR=AA                                                             BDTR1540
  370 DD=DLXX/5.D0                                                      BDTR1550
      DD=DEXP(DD)                                                       BDTR1560
      DD=(RR-1.D0)-(RR+BB-1.D0)*XX*DD +2.D0                             BDTR1570
      IF(DD) 420,420,430                                                BDTR1580
  380 IF(BB-1.D0) 390,390,400                                           BDTR1590
  390 RR=BB+1.D0                                                        BDTR1600
      GO TO 410                                                         BDTR1610
  400 RR=BB                                                             BDTR1620
  410 DD=DL1X/5.D0                                                      BDTR1630
      DD=DEXP(DD)                                                       BDTR1640
      DD=(RR-1.D0)-(AA+RR-1.D0)*(1.D0-XX)*DD +2.D0                      BDTR1650
      IF(DD) 430,430,420                                                BDTR1660
  420 ID=1                                                              BDTR1670
      FF=DL1X                                                           BDTR1680
      DL1X=DLXX                                                         BDTR1690
      DLXX=FF                                                           BDTR1700
      XO=1.D0/XO                                                        BDTR1710
      FF=AA                                                             BDTR1720
      AA=BB                                                             BDTR1730
      BB=FF                                                             BDTR1740
      G2=G1                                                             BDTR1750
C                                                                       BDTR1760
C        TEST FOR A LESS THAN 1.0                                       BDTR1770
C                                                                       BDTR1780
  430 FF=0.D0                                                           BDTR1790
      IF(AA-1.D0) 440,440,470                                           BDTR1800
  440 CALL DLGAM(AA+1.D0,G4,IOK)                                        BDTR1810
      DD=AA*DLXX+BB*DL1X+G3-G2-G4                                       BDTR1820
      IF(DD+1.68D02) 460,460,450                                        BDTR1830
  450 FF=FF+DEXP(DD)                                                    BDTR1840
  460 AA=AA+1.D0                                                        BDTR1850
C                                                                       BDTR1860
C        COMPUTE P USING CONTINUED FRACTION EXPANSION                   BDTR1870
C                                                                       BDTR1880
  470 FN=AA+BB-1.D0                                                     BDTR1890
      RR=AA-1.D0                                                        BDTR1900
      II=80                                                             BDTR1910
      XI=DFLOAT(II)                                                     BDTR1920
      SS=((BB-XI)*(RR+XI))/((RR+2.D0*XI-1.D0)*(RR+2.D0*XI))             BDTR1930
      SS=SS*XO                                                          BDTR1940
      DO 480 I=1,79                                                     BDTR1950
      II=80-I                                                           BDTR1960
      XI=DFLOAT(II)                                                     BDTR1970
      DD=(XI*(FN+XI))/((RR+2.D0*XI+1.D0)*(RR+2.D0*XI))                  BDTR1980
      DD=DD*XO                                                          BDTR1990
      CC=((BB-XI)*(RR+XI))/((RR+2.D0*XI-1.D0)*(RR+2.D0*XI))             BDTR2000
      CC=CC*XO                                                          BDTR2010
      SS=CC/(1.D0+DD/(1.D0-SS))                                         BDTR2020
  480 CONTINUE                                                          BDTR2030
      SS=1.D0/(1.D0-SS)                                                 BDTR2040
      IF(SS) 650,650,490                                                BDTR2050
  490 CALL DLGAM(AA+BB,G1,IOK)                                          BDTR2060
      CALL DLGAM(AA+1.D0,G4,IOK)                                        BDTR2070
      CC=G1-G2-G4+AA*DLXX+(BB-1.D0)*DL1X                                BDTR2080
      PP=CC+DLOG(SS)                                                    BDTR2090
      IF(PP+1.68D02) 500,500,510                                        BDTR2100
  500 PP=FF                                                             BDTR2110
      GO TO 520                                                         BDTR2120
  510 PP=DEXP(PP)+FF                                                    BDTR2130
  520 IF(ID) 540,540,530                                                BDTR2140
  530 PP=1.D0-PP                                                        BDTR2150
  540 P=SNGL(PP)                                                        BDTR2160
C                                                                       BDTR2170
C        SET ERROR INDICATOR                                            BDTR2180
C                                                                       BDTR2190
      IF(P) 550,570,570                                                 BDTR2200
  550 IF(ABS(P)-1.E-7) 560,560,650                                      BDTR2210
  560 P=0.0                                                             BDTR2220
      GO TO 660                                                         BDTR2230
  570 IF(1.-P) 580,600,600                                              BDTR2240
  580 IF(ABS(1.-P)-1.E-7) 590,590,650                                   BDTR2250
  590 P=1.0                                                             BDTR2260
      GO TO 660                                                         BDTR2270
  600 IF(P-1.E-8) 610,610,620                                           BDTR2280
  610 P=0.0                                                             BDTR2290
      GO TO 660                                                         BDTR2300
  620 IF((1.0-P)-1.E-8) 630,630,660                                     BDTR2310
  630 P=1.0                                                             BDTR2320
      GO TO 660                                                         BDTR2330
  640 IER=-2                                                            BDTR2340
      D=-1.E75                                                          BDTR2350
      P=-1.E75                                                          BDTR2360
      GO TO 670                                                         BDTR2370
  650 IER=+2                                                            BDTR2380
      P= 1.E75                                                          BDTR2390
      GO TO 670                                                         BDTR2400
  660 IER=0                                                             BDTR2410
  670 RETURN                                                            BDTR2420
      END                                                               BDTR2430
