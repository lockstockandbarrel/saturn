C                                                                       PROB  10
C     ..................................................................PROB  20
C                                                                       PROB  30
C        SUBROUTINE PROBT                                               PROB  40
C                                                                       PROB  50
C        PURPOSE                                                        PROB  60
C           TO OBTAIN MAXIMUM LIKELIHOOD ESTIMATES FOR THE PARAMETERS A PROB  70
C           AND B IN THE PROBIT EQUATION  Y = A + BX.  AN ITERATIVE     PROB  80
C           SCHEME IS USED.  THE INPUT TO THE SUBROUTINE CONSISTS OF K  PROB  90
C           DIFFERENT DOSAGE LEVELS APPLIED TO K GROUPS OF SUBJECTS, ANDPROB 100
C           THE NUMBER OF SUBJECTS IN EACH GROUP RESPONDING TO THE      PROB 110
C           RESPECTIVE DOSAGE OF THE DRUG.                              PROB 120
C                                                                       PROB 130
C        USAGE                                                          PROB 140
C           CALL PROBT (K,X,S,R,LOG,ANS,W1,W2,IER)                      PROB 150
C                                                                       PROB 160
C        DESCRIPTION OF PARAMETERS                                      PROB 170
C           K   - NUMBER OF DIFFERENT DOSE LEVELS OF THE DRUG.  K SHOULDPROB 180
C                 BE GREATER THAN 2.                                    PROB 190
C           X   - INPUT VECTOR OF LENGTH K CONTAINING THE DOSE LEVEL OF PROB 200
C                 THE DRUG TESTED.  X MUST BE NON-NEGATIVE.             PROB 210
C           S   - INPUT VECTOR OF LENGTH K CONTAINING THE NUMBER OF     PROB 220
C                 SUBJECTS TESTED AT EACH DOSE LEVEL                    PROB 230
C           R   - INPUT VECTOR OF LENGTH K CONTAINING THE NUMBER OF     PROB 240
C                 SUBJECTS AT EACH LEVEL RESPONDING TO THE DRUG         PROB 250
C           LOG - INPUT OPTION CODE                                     PROB 260
C                 1- IF IT IS DESIRED TO CONVERT THE DOSE LEVELS TO     PROB 270
C                    COMMON LOGARITHMS.  THE DOSAGE LEVELS SHOULD BE    PROB 280
C                    NON-NULL IN THIS CASE.                             PROB 290
C                 0- IF NO CONVERSION IS DESIRED                        PROB 300
C           ANS - OUTPUT VECTOR OF LENGTH 4 CONTAINING THE FOLLOWING    PROB 310
C                 RESULTS                                               PROB 320
C                 ANS(1)- ESTIMATE OF THE INTERCEPT CONSTANT A          PROB 330
C                 ANS(2)- ESTIMATE OF THE PROBIT REGRESSION COEFFICIENT PROB 340
C                         B                                             PROB 350
C                 ANS(3)- CHI-SQUARED VALUE FOR A TEST OF SIGNIFICANCE  PROB 360
C                         OF THE FINAL PROBIT EQUATION                  PROB 370
C                 ANS(4)- DEGREES OF FREEDOM FOR THE CHI-SQUARE         PROB 380
C                         STATISTIC                                     PROB 390
C           W1  - OUTPUT VECTOR OF LENGTH K CONTAINING THE PROPORTIONS  PROB 400
C                 OF SUBJECTS RESPONDING TO THE VARIOUS DOSE LEVELS OF  PROB 410
C                 THE DRUG                                              PROB 420
C           W2  - OUTPUT VECTOR OF LENGTH K CONTAINING THE VALUES OF THEPROB 430
C                 EXPECTED PROBIT FOR THE VARIOUS LEVELS OF A DRUG      PROB 440
C           IER - 1 IF K IS NOT GREATER THAN 2.                         PROB 450
C                 2 IF SOME DOSAGE LEVEL IS NEGATIVE, OR IF THE INPUT   PROB 460
C                   OPTION CODE LOG IS 1 AND SOME DOSAGE LEVEL IS ZERO. PROB 470
C                 3 IF SOME ELEMENT OF S IS NOT POSITIVE.               PROB 480
C                 4 IF NUMBER OF SUBJECTS RESPONDING IS GREATER THAN    PROB 490
C                 NUMBER OF SUBJECTS TESTED.                            PROB 500
C                 ONLY IF IER IS ZERO IS A PROBIT ANALYSIS PERFORMED.   PROB 510
C                 OTHERWISE, ANS, W1, AND W2 ARE SET TO ZERO.           PROB 520
C                                                                       PROB 530
C        REMARKS                                                        PROB 540
C           THE PROGRAM WILL ITERATE ON THE PROBIT EQUATION UNTIL TWO   PROB 550
C           SUCCESSIVE SOLUTIONS PRODUCE CHANGES OF LESS THAN 10**(-7). PROB 560
C                                                                       PROB 570
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PROB 580
C           NDTR                                                        PROB 590
C           NDTRI                                                       PROB 600
C                                                                       PROB 610
C        METHOD                                                         PROB 620
C           REFER TO D. J. FINNEY, PROBIT ANALYSIS. 2ND ED. (CAMBRIDGE, PROB 630
C           1952)                                                       PROB 640
C                                                                       PROB 650
C     ..................................................................PROB 660
C                                                                       PROB 670
      SUBROUTINE PROBT (K,X,S,R,LOG,ANS,W1,W2,IER)                      PROB 680
C                                                                       PROB 690
      DIMENSION X(1),S(1),R(1),ANS(1),W1(1),W2(1)                       PROB 700
C                                                                       PROB 710
C        TEST WHETHER LOG CONVERSION IS NEEDED                          PROB 720
C                                                                       PROB 730
      IER=0                                                             PROB 740
      IF(K-2)5,5,7                                                      PROB 750
    5 IER = 1                                                           PROB 760
      GO TO 90                                                          PROB 770
    7 DO 8 I=1,K                                                        PROB 780
      IF(X(I))12,8,8                                                    PROB 790
    8 CONTINUE                                                          PROB 800
      IF(LOG-1) 16,10,16                                                PROB 810
   10 DO 15 I=1,K                                                       PROB 820
      IF(X(I))12,12,14                                                  PROB 830
   12 IER=2                                                             PROB 840
      GO TO 90                                                          PROB 850
   14 X(I)= ALOG10(X(I))                                                PROB 860
   15 CONTINUE                                                          PROB 870
C                                                                       PROB 880
C        COMPUTE PROPORTIONS OF OBJECTS RESPONDING                      PROB 890
C                                                                       PROB 900
   16 DO 18 I=1,K                                                       PROB 910
      IF(S(I)-R(I)) 17,18,18                                            PROB 920
   17 IER=4                                                             PROB 930
      GO TO 90                                                          PROB 940
   18 CONTINUE                                                          PROB 950
   20 DO 23 I=1,K                                                       PROB 960
      IF(S(I))21,21,22                                                  PROB 970
   21 IER=3                                                             PROB 980
      GO TO 90                                                          PROB 990
   22 W1(I)=R(I)/S(I)                                                   PROB1000
   23 CONTINUE                                                          PROB1010
C                                                                       PROB1020
C        COMPUTE INITIAL ESTIMATES OF INTERCEPT AND PROBIT REGRESSION   PROB1030
C        COEFFICIENT                                                    PROB1040
C                                                                       PROB1050
      WN=0.0                                                            PROB1060
      XBAR=0.0                                                          PROB1070
      SNWY=0.0                                                          PROB1080
      SXX=0.0                                                           PROB1090
      SXY=0.0                                                           PROB1100
C                                                                       PROB1110
      DO 30 I=1,K                                                       PROB1120
      P=W1(I)                                                           PROB1130
      IF(P) 30, 30, 24                                                  PROB1140
   24 IF(P-1.0) 25, 30, 30                                              PROB1150
   25 WN=WN+1.0                                                         PROB1160
C                                                                       PROB1170
      CALL NDTRI (P,Z,D,IER)                                            PROB1180
C                                                                       PROB1190
      Z=Z+5.0                                                           PROB1200
      XBAR=XBAR+X(I)                                                    PROB1210
      SNWY=SNWY+Z                                                       PROB1220
      SXX=SXX+X(I)**2                                                   PROB1230
      SXY=SXY+X(I)*Z                                                    PROB1240
   30 CONTINUE                                                          PROB1250
C                                                                       PROB1260
      B=(SXY-(XBAR*SNWY)/WN)/(SXX-(XBAR*XBAR)/WN)                       PROB1270
      XBAR=XBAR/WN                                                      PROB1280
      SNWY=SNWY/WN                                                      PROB1290
      A=SNWY-B*XBAR                                                     PROB1300
      DD=0.0                                                            PROB1310
C                                                                       PROB1320
C        COMPUTE EXPECTED PROBIT                                        PROB1330
C                                                                       PROB1340
      DO 31 I=1,K                                                       PROB1350
   31 W2(I)=A+B*X(I)                                                    PROB1360
C                                                                       PROB1370
   33 SNW=0.0                                                           PROB1380
      SNWX=0.0                                                          PROB1390
      SNWY=0.0                                                          PROB1400
      SNWXX=0.0                                                         PROB1410
      SNWXY=0.0                                                         PROB1420
      DO 50 I=1,K                                                       PROB1430
      Y=W2(I)                                                           PROB1440
C                                                                       PROB1450
C        FIND A WEIGHTING COEFFICIENT FOR PROBIT ANALYSIS               PROB1460
C                                                                       PROB1470
      D=Y-5.0                                                           PROB1480
C                                                                       PROB1490
      CALL NDTR (D,P,Z)                                                 PROB1500
C                                                                       PROB1510
      Q=1.0-P                                                           PROB1520
      W=(Z*Z)/(P*Q)                                                     PROB1530
C                                                                       PROB1540
C        COMPUTE WORKING PROBIT                                         PROB1550
C                                                                       PROB1560
      IF(Y-5.0) 35, 35, 40                                              PROB1570
   35 WP=(Y-P/Z)+W1(I)/Z                                                PROB1580
      GO TO 45                                                          PROB1590
   40 WP=(Y+Q/Z)-(1.0-W1(I))/Z                                          PROB1600
C                                                                       PROB1610
C        SUM INTERMEDIATE RESULTS                                       PROB1620
C                                                                       PROB1630
   45 WN=W*S(I)                                                         PROB1640
      SNW=SNW+WN                                                        PROB1650
      SNWX=SNWX+WN*X(I)                                                 PROB1660
      SNWY=SNWY+WN*WP                                                   PROB1670
      SNWXX=SNWXX+WN*X(I)**2                                            PROB1680
   50 SNWXY=SNWXY+WN*X(I)*WP                                            PROB1690
C                                                                       PROB1700
C        COMPUTE NEW ESTIMATES OF INTERCEPT AND COEFFICIENT             PROB1710
C                                                                       PROB1720
      XBAR=SNWX/SNW                                                     PROB1730
C                                                                       PROB1740
      SXX=SNWXX-(SNWX)*(SNWX)/SNW                                       PROB1750
      SXY=SNWXY-(SNWX)*(SNWY)/SNW                                       PROB1760
      B=SXY/SXX                                                         PROB1770
C                                                                       PROB1780
      A=SNWY/SNW-B*XBAR                                                 PROB1790
C                                                                       PROB1800
C        EXAMINE THE CHANGES IN Y                                       PROB1810
C                                                                       PROB1820
      SXX=0.0                                                           PROB1830
      DO 60 I=1,K                                                       PROB1840
      Y=A+B*X(I)                                                        PROB1850
      D=W2(I)-Y                                                         PROB1860
      SXX=SXX+D*D                                                       PROB1870
   60 W2(I)=Y                                                           PROB1880
      IF(( ABS(DD-SXX))-(1.0E-7)) 65, 65, 63                            PROB1890
   63 DD=SXX                                                            PROB1900
      GO TO 33                                                          PROB1910
C                                                                       PROB1920
C        STORE INTERCEPT AND COEFFICIENT                                PROB1930
C                                                                       PROB1940
   65 ANS(1)=A                                                          PROB1950
      ANS(2)=B                                                          PROB1960
C                                                                       PROB1970
C        COMPUTE CHI-SQUARE                                             PROB1980
C                                                                       PROB1990
      ANS(3)=0.0                                                        PROB2000
      DO 70 I=1,K                                                       PROB2010
      Y=W2(I)-5.0                                                       PROB2020
C                                                                       PROB2030
      CALL NDTR (Y,P,D)                                                 PROB2040
C                                                                       PROB2050
      AA=R(I)-S(I)*P                                                    PROB2060
      DD=S(I)*P*(1.0-P)                                                 PROB2070
   70 ANS(3)=ANS(3)+AA*AA/DD                                            PROB2080
C                                                                       PROB2090
C        DEGREES OF FREEDOM FOR CHI-SQUARE                              PROB2100
C                                                                       PROB2110
      ANS(4)=K-2                                                        PROB2120
C                                                                       PROB2130
   80 RETURN                                                            PROB2140
   90 DO 100 I=1,K                                                      PROB2150
      W1(I)=0.0                                                         PROB2160
  100 W2(I)=0.0                                                         PROB2170
      DO 110 I=1,4                                                      PROB2180
  110 ANS(I)=0.0                                                        PROB2190
      GO TO 80                                                          PROB2200
      END                                                               PROB2210
