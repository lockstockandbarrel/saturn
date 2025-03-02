C                                                                       DISC  10
C     ..................................................................DISC  20
C                                                                       DISC  30
C        SUBROUTINE DISCR                                               DISC  40
C                                                                       DISC  50
C        PURPOSE                                                        DISC  60
C           COMPUTE A SET OF LINEAR FUNCTIONS WHICH SERVE AS INDICES    DISC  70
C           FOR CLASSIFYING AN INDIVIDUAL INTO ONE OF SEVERAL GROUPS.   DISC  80
C           NORMALLY THIS SUBROUTINE IS USED IN THE PERFORMANCE OF      DISC  90
C           DISCRIMINANT ANALYSIS.                                      DISC 100
C                                                                       DISC 110
C        USAGE                                                          DISC 120
C           CALL DISCR (K,M,N,X,XBAR,D,CMEAN,V,C,P,LG)                  DISC 130
C                                                                       DISC 140
C        DESCRIPTION OF PARAMETERS                                      DISC 150
C           K     - NUMBER OF GROUPS. K MUST BE GREATER THAN ONE.       DISC 160
C           M     - NUMBER OF VARIABLES                                 DISC 170
C           N     - INPUT VECTOR OF LENGTH K CONTAINING SAMPLE SIZES OF DISC 180
C                   GROUPS.                                             DISC 190
C           X     - INPUT VECTOR CONTAINING DATA IN THE MANNER EQUIVA-  DISC 200
C                   LENT TO A 3-DIMENSIONAL FORTRAN ARRAY, X(1,1,1),    DISC 210
C                   X(2,1,1), X(3,1,1), ETC.  THE FIRST SUBSCRIPT IS    DISC 220
C                   CASE NUMBER, THE SECOND SUBSCRIPT IS VARIABLE NUMBERDISC 230
C                   AND THE THIRD SUBSCRIPT IS GROUP NUMBER.  THE       DISC 240
C                   LENGTH OF VECTOR X IS EQUAL TO THE TOTAL NUMBER OF  DISC 250
C                   DATA POINTS, T*M, WHERE T = N(1)+N(2)+...+N(K).     DISC 260
C           XBAR  - INPUT MATRIX (M X K) CONTAINING MEANS OF M VARIABLESDISC 270
C                   IN K GROUPS                                         DISC 280
C           D     - INPUT MATRIX (M X M) CONTAINING THE INVERSE OF      DISC 290
C                   POOLED DISPERSION MATRIX.                           DISC 300
C           CMEAN - OUTPUT VECTOR OF LENGTH M CONTAINING COMMON MEANS.  DISC 310
C           V     - OUTPUT VARIABLE CONTAINING GENERALIZED MAHALANOBIS  DISC 320
C                   D-SQUARE.                                           DISC 330
C           C     - OUTPUT MATRIX (M+1 X K) CONTAINING THE COEFFICIENTS DISC 340
C                   OF DISCRIMINANT FUNCTIONS.  THE FIRST POSITION OF   DISC 350
C                   EACH COLUMN (FUNCTION) CONTAINS THE VALUE OF THE    DISC 360
C                   CONSTANT FOR THAT FUNCTION.                         DISC 370
C           P     - OUTPUT VECTOR CONTAINING THE PROBABILITY ASSOCIATED DISC 380
C                   WITH THE LARGEST DISCRIMINANT FUNCTIONS OF ALL CASESDISC 390
C                   IN ALL GROUPS.  CALCULATED RESULTS ARE STORED IN THEDISC 400
C                   MANNER EQUIVALENT TO A 2-DIMENSIONAL AREA (THE      DISC 410
C                   FIRST SUBSCRIPT IS CASE NUMBER, AND THE SECOND      DISC 420
C                   SUBSCRIPT IS GROUP NUMBER).  VECTOR P HAS LENGTH    DISC 430
C                   EQUAL TO THE TOTAL NUMBER OF CASES, T (T = N(1)+N(2)DISC 440
C                   +...+N(K)).                                         DISC 450
C           LG    - OUTPUT VECTOR CONTAINING THE SUBSCRIPTS OF THE      DISC 460
C                   LARGEST DISCRIMINANT FUNCTIONS STORED IN VECTOR P.  DISC 470
C                   THE LENGTH OF VECTOR LG IS THE SAME AS THE LENGTH   DISC 480
C                   OF VECTOR P.                                        DISC 490
C                                                                       DISC 500
C        REMARKS                                                        DISC 510
C           THE NUMBER OF VARIABLES MUST BE GREATER THAN OR EQUAL TO    DISC 520
C           THE NUMBER OF GROUPS.                                       DISC 530
C                                                                       DISC 540
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DISC 550
C           NONE                                                        DISC 560
C                                                                       DISC 570
C        METHOD                                                         DISC 580
C           REFER TO 'BMD COMPUTER PROGRAMS MANUAL', EDITED BY W. J.    DISC 590
C           DIXON, UCLA, 1964, AND T. W. ANDERSON, 'INTRODUCTION TO     DISC 600
C           MULTIVARIATE STATISTICAL ANALYSIS', JOHN WILEY AND SONS,    DISC 610
C           1958, SECTION 6.6-6.8.                                      DISC 620
C                                                                       DISC 630
C     ..................................................................DISC 640
C                                                                       DISC 650
      SUBROUTINE DISCR (K,M,N,X,XBAR,D,CMEAN,V,C,P,LG)                  DISC 660
      DIMENSION N(1),X(1),XBAR(1),D(1),CMEAN(1),C(1),P(1),LG(1)         DISC 670
C                                                                       DISC 680
C        ...............................................................DISC 690
C                                                                       DISC 700
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE  DISC 710
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION      DISC 720
C        STATEMENT WHICH FOLLOWS.                                       DISC 730
C                                                                       DISC 740
C     DOUBLE PRECISION XBAR,D,CMEAN,V,C,SUM,P,PL                        DISC 750
C                                                                       DISC 760
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS    DISC 770
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS      DISC 780
C        ROUTINE.                                                       DISC 790
C                                                                       DISC 800
C        THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO      DISC 810
C        CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  EXP IN STATEMENT  DISC 820
C        250 MUST BE CHANGED TO DEXP.                                   DISC 830
C                                                                       DISC 840
C        ...............................................................DISC 850
C                                                                       DISC 860
C     CALCULATE COMMON MEANS                                            DISC 870
C                                                                       DISC 880
      N1=N(1)                                                           DISC 890
      DO 100 I=2,K                                                      DISC 900
  100 N1=N1+N(I)                                                        DISC 910
      FNT=N1                                                            DISC 920
      DO 110 I=1,K                                                      DISC 930
  110 P(I)=N(I)                                                         DISC 940
      DO 130 I=1,M                                                      DISC 950
      CMEAN(I)=0                                                        DISC 960
      N1=I-M                                                            DISC 970
      DO 120 J=1,K                                                      DISC 980
      N1=N1+M                                                           DISC 990
  120 CMEAN(I)=CMEAN(I)+P(J)*XBAR(N1)                                   DISC1000
  130 CMEAN(I)=CMEAN(I)/FNT                                             DISC1010
C                                                                       DISC1020
C     CALCULATE GENERALIZED MAHALANOBIS D SQUARE                        DISC1030
C                                                                       DISC1040
      L=0                                                               DISC1050
      DO 140 I=1,K                                                      DISC1060
      DO 140 J=1,M                                                      DISC1070
      L=L+1                                                             DISC1080
  140 C(L)=XBAR(L)-CMEAN(J)                                             DISC1090
      V=0.0                                                             DISC1100
      L=0                                                               DISC1110
      DO 160 J=1,M                                                      DISC1120
      DO 160 I=1,M                                                      DISC1130
      N1=I-M                                                            DISC1140
      N2=J-M                                                            DISC1150
      SUM=0.0                                                           DISC1160
      DO 150 IJ=1,K                                                     DISC1170
      N1=N1+M                                                           DISC1180
      N2=N2+M                                                           DISC1190
  150 SUM=SUM+P(IJ)*C(N1)*C(N2)                                         DISC1200
      L=L+1                                                             DISC1210
  160 V=V+D(L)*SUM                                                      DISC1220
C                                                                       DISC1230
C     CALCULATE THE COEFFICIENTS OF DISCRIMINANT FUNCTIONS              DISC1240
C                                                                       DISC1250
      N2=0                                                              DISC1260
      DO 190 KA=1,K                                                     DISC1270
      DO 170 I=1,M                                                      DISC1280
      N2=N2+1                                                           DISC1290
  170 P(I)=XBAR(N2)                                                     DISC1300
      IQ=(M+1)*(KA-1)+1                                                 DISC1310
      SUM=0.0                                                           DISC1320
      DO 180 J=1,M                                                      DISC1330
      N1=J-M                                                            DISC1340
      DO 180 L=1,M                                                      DISC1350
      N1=N1+M                                                           DISC1360
  180 SUM=SUM+D(N1)*P(J)*P(L)                                           DISC1370
      C(IQ)=-(SUM/2.0)                                                  DISC1380
      DO 190 I=1,M                                                      DISC1390
      N1=I-M                                                            DISC1400
      IQ=IQ+1                                                           DISC1410
      C(IQ)=0.0                                                         DISC1420
      DO 190 J=1,M                                                      DISC1430
      N1=N1+M                                                           DISC1440
  190 C(IQ)=C(IQ)+D(N1)*P(J)                                            DISC1450
C                                                                       DISC1460
C     FOR EACH CASE IN EACH GROUP, CALCULATE..                          DISC1470
C                                                                       DISC1480
C        DISCRIMINANT FUNCTIONS                                         DISC1490
C                                                                       DISC1500
      LBASE=0                                                           DISC1510
      N1=0                                                              DISC1520
      DO 270 KG=1,K                                                     DISC1530
      NN=N(KG)                                                          DISC1540
      DO 260 I=1,NN                                                     DISC1550
      L=I-NN+LBASE                                                      DISC1560
      DO 200 J=1,M                                                      DISC1570
      L=L+NN                                                            DISC1580
  200 D(J)=X(L)                                                         DISC1590
      N2=0                                                              DISC1600
      DO 220 KA=1,K                                                     DISC1610
      N2=N2+1                                                           DISC1620
      SUM=C(N2)                                                         DISC1630
      DO 210 J=1,M                                                      DISC1640
      N2=N2+1                                                           DISC1650
  210 SUM=SUM+C(N2)*D(J)                                                DISC1660
  220 XBAR(KA)=SUM                                                      DISC1670
C                                                                       DISC1680
C        THE LARGEST DISCRIMINANT FUNCTION                              DISC1690
C                                                                       DISC1700
      L=1                                                               DISC1710
      SUM=XBAR(1)                                                       DISC1720
      DO 240 J=2,K                                                      DISC1730
      IF(SUM-XBAR(J)) 230, 240, 240                                     DISC1740
  230 L=J                                                               DISC1750
      SUM=XBAR(J)                                                       DISC1760
  240 CONTINUE                                                          DISC1770
C                                                                       DISC1780
C        PROBABILITY ASSOCIATED WITH THE LARGEST DISCRIMINANT FUNCTION  DISC1790
C                                                                       DISC1800
      PL=0.0                                                            DISC1810
      DO 250 J=1,K                                                      DISC1820
  250 PL=PL+ EXP(XBAR(J)-SUM)                                           DISC1830
      N1=N1+1                                                           DISC1840
      LG(N1)=L                                                          DISC1850
  260 P(N1)=1.0/PL                                                      DISC1860
  270 LBASE=LBASE+NN*M                                                  DISC1870
C                                                                       DISC1880
      RETURN                                                            DISC1890
      END                                                               DISC1900
