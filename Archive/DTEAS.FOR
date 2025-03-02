C                                                                       DTEA  10
C     ..................................................................DTEA  20
C                                                                       DTEA  30
C        SUBROUTINE DTEAS                                               DTEA  40
C                                                                       DTEA  50
C        PURPOSE                                                        DTEA  60
C           CALCULATE THE LIMIT OF A GIVEN SEQUENCE BY MEANS OF THE     DTEA  70
C           EPSILON-ALGORITHM.                                          DTEA  80
C                                                                       DTEA  90
C        USAGE                                                          DTEA 100
C           CALL DTEAS(X,N,FIN,EPS,IER)                                 DTEA 110
C                                                                       DTEA 120
C        DESCRIPTION OF PARAMETERS                                      DTEA 130
C           X      - DOUBLE PRECISION VECTOR WHOSE COMPONENTS ARE TERMS DTEA 140
C                    OF THE GIVEN SEQUENCE. ON RETURN THE COMPONENTS OF DTEA 150
C                    VECTOR X ARE DESTROYED.                            DTEA 160
C           N      - DIMENSION OF INPUT VECTOR X.                       DTEA 170
C           FIN    - RESULTANT SCALAR IN DOUBLE PRECISION CONTAINING ON DTEA 180
C                    RETURN THE LIMIT OF THE GIVEN SEQUENCE.            DTEA 190
C           EPS    - SINGLE PRECISION INPUT VALUE, WHICH SPECIFIES THE  DTEA 200
C                    UPPER BOUND OF THE RELATIVE (ABSOLUTE) ERROR IF THEDTEA 210
C                    COMPONENTS OF X ARE ABSOLUTELY GREATER (LESS) THAN DTEA 220
C                    ONE.                                               DTEA 230
C                    CALCULATION IS TERMINATED AS SOON AS THREE TIMES INDTEA 240
C                    SUCCESSION THE RELATIVE (ABSOLUTE) DIFFERENCE      DTEA 250
C                    BETWEEN NEIGHBOURING TERMS IS NOT GREATER THAN EPS.DTEA 260
C           IER    - RESULTANT ERROR PARAMETER CODED IN THE FOLLOWING   DTEA 270
C                    FORM                                               DTEA 280
C                     IER=0  - NO ERROR                                 DTEA 290
C                     IER=1  - REQUIRED ACCURACY NOT REACHED WITH       DTEA 300
C                              MAXIMAL NUMBER OF ITERATIONS             DTEA 310
C                     IER=-1 - INTEGER N IS LESS THAN TEN.              DTEA 320
C                                                                       DTEA 330
C        REMARKS                                                        DTEA 340
C           NO ACTION BESIDES ERROR MESSAGE IN CASE N LESS THAN TEN.    DTEA 350
C           THE CHARACTER OF THE GIVEN INFINITE SEQUENCE MUST BE        DTEA 360
C           RECOGNIZABLE BY THOSE N COMPONENTS OF THE INPUT VECTOR X.   DTEA 370
C                                                                       DTEA 380
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DTEA 390
C           NONE                                                        DTEA 400
C                                                                       DTEA 410
C        METHOD                                                         DTEA 420
C           THE CONVERGENCE OF THE GIVEN SEQUENCE IS ACCELERATED BY     DTEA 430
C           MEANS OF THE E(2)-TRANSFORMATION, USED IN AN ITERATIVE WAY. DTEA 440
C           FOR REFERENCE, SEE                                          DTEA 450
C           ALGORITHM 215,SHANKS, CACM 1963, NO. 11, PP. 662. AND       DTEA 460
C           P. WYNN, SINGULAR RULES FOR CERTAIN NON-LINEAR ALGORITHMS   DTEA 470
C           BIT VOL. 3, 1963, PP. 175-195.                              DTEA 480
C                                                                       DTEA 490
C     ..................................................................DTEA 500
C                                                                       DTEA 510
      SUBROUTINE DTEAS(X,N,FIN,EPS,IER)                                 DTEA 520
C                                                                       DTEA 530
      DIMENSION X(1)                                                    DTEA 540
      DOUBLE PRECISION X,FIN,W1,W2,W3,W4,W5,W6,W7,T                     DTEA 550
C                                                                       DTEA 560
C        TEST ON WRONG INPUT PARAMETER N                                DTEA 570
C                                                                       DTEA 580
      NEW=N                                                             DTEA 590
      IF(NEW-10)1,2,2                                                   DTEA 600
    1 IER=-1                                                            DTEA 610
      RETURN                                                            DTEA 620
C                                                                       DTEA 630
C        CALCULATE INITIAL VALUES FOR THE EPSILON ARRAY                 DTEA 640
C                                                                       DTEA 650
    2 ISW1=0                                                            DTEA 660
      ISW2=0                                                            DTEA 670
      W1=1.D38                                                          DTEA 680
      W7=X(4)-X(3)                                                      DTEA 690
      IF(W7)3,4,3                                                       DTEA 700
    3 W1=1.D0/W7                                                        DTEA 710
C                                                                       DTEA 720
    4 W5=1.D38                                                          DTEA 730
      W7=X(2)-X(1)                                                      DTEA 740
      IF(W7)5,6,5                                                       DTEA 750
    5 W5=1.D0/W7                                                        DTEA 760
C                                                                       DTEA 770
    6 W4=X(3)-X(2)                                                      DTEA 780
      IF(W4)9,7,9                                                       DTEA 790
    7 W4=1.D38                                                          DTEA 800
      T=X(2)                                                            DTEA 810
      W2=X(3)                                                           DTEA 820
    8 W3=1.D38                                                          DTEA 830
      GO TO 17                                                          DTEA 840
C                                                                       DTEA 850
    9 W4=1.D0/W4                                                        DTEA 860
C                                                                       DTEA 870
      T=1.D38                                                           DTEA 880
      W7=W4-W5                                                          DTEA 890
      IF(W7)10,11,10                                                    DTEA 900
   10 T=X(2)+1.D0/W7                                                    DTEA 910
C                                                                       DTEA 920
   11 W2=W1-W4                                                          DTEA 930
      IF(W2)15,12,15                                                    DTEA 940
   12 W2=1.D38                                                          DTEA 950
      IF(T-1.D38)13,14,14                                               DTEA 960
   13 ISW2=1                                                            DTEA 970
   14 W3=W4                                                             DTEA 980
      GO TO 17                                                          DTEA 990
C                                                                       DTEA1000
   15 W2=X(3)+1.D0/W2                                                   DTEA1010
      W7=W2-T                                                           DTEA1020
      IF(W7)16,8,16                                                     DTEA1030
   16 W3=W4+1.D0/W7                                                     DTEA1040
C                                                                       DTEA1050
   17 ISW1=ISW2                                                         DTEA1060
      ISW2=0                                                            DTEA1070
      IMIN=4                                                            DTEA1080
C                                                                       DTEA1090
C        CALCULATE DIAGONALS OF THE EPSILON ARRAY IN A DO-LOOP          DTEA1100
C                                                                       DTEA1110
      DO 40 I=5,NEW                                                     DTEA1120
      IAUS=I-IMIN                                                       DTEA1130
      W4=1.D38                                                          DTEA1140
      W5=X(I-1)                                                         DTEA1150
      W7=X(I)-X(I-1)                                                    DTEA1160
      IF(W7)18,24,18                                                    DTEA1170
   18 W4=1.D0/W7                                                        DTEA1180
C                                                                       DTEA1190
      IF(W1-1.D38)19,25,25                                              DTEA1200
   19 W6=W4-W1                                                          DTEA1210
C                                                                       DTEA1220
C        TEST FOR NECESSITY OF A SINGULAR RULE                          DTEA1230
C                                                                       DTEA1240
      IF(DABS(W6)-DABS(W4)*1.D-12)20,20,22                              DTEA1250
   20 ISW2=1                                                            DTEA1260
      IF(W6)22,21,22                                                    DTEA1270
   21 W5=1.D38                                                          DTEA1280
      W6=W1                                                             DTEA1290
      IF(W2-1.D38)28,26,26                                              DTEA1300
   22 W5=X(I-1)+1.D0/W6                                                 DTEA1310
C                                                                       DTEA1320
C        FIRST TEST FOR LOSS OF SIGNIFICANCE                            DTEA1330
C                                                                       DTEA1340
      IF(DABS(W5)-DABS(X(I-1))*1.D-10)23,24,24                          DTEA1350
   23 IF(W5)36,24,36                                                    DTEA1360
C                                                                       DTEA1370
   24 W7=W5-W2                                                          DTEA1380
      IF(W7)27,25,27                                                    DTEA1390
   25 W6=1.D38                                                          DTEA1400
   26 ISW2=0                                                            DTEA1410
      X(IAUS)=W2                                                        DTEA1420
      GO TO 37                                                          DTEA1430
   27 W6=W1+1.D0/W7                                                     DTEA1440
   28 IF(ISW1-1)33,29,29                                                DTEA1450
C                                                                       DTEA1460
C        CALCULATE X(IAUS) WITH HELP OF SINGULAR RULE                   DTEA1470
C                                                                       DTEA1480
   29 IF(W2-1.D38)30,32,32                                              DTEA1490
   30 W7=W5/(W2-W5)+T/(W2-T)+X(I-2)/(X(I-2)-W2)                         DTEA1500
      IF(1.D0+W7)31,38,31                                               DTEA1510
   31 X(IAUS)=W7*W2/(1.D0+W7)                                           DTEA1520
      GO TO 39                                                          DTEA1530
C                                                                       DTEA1540
   32 X(IAUS)=W5+T-X(I-2)                                               DTEA1550
      GO TO 39                                                          DTEA1560
C                                                                       DTEA1570
   33 W7=W6-W3                                                          DTEA1580
      IF(W7)34,38,34                                                    DTEA1590
   34 X(IAUS)=W2+1.D0/W7                                                DTEA1600
C                                                                       DTEA1610
C        SECOND TEST FOR LOSS OF SIGNIFICANCE                           DTEA1620
C                                                                       DTEA1630
      IF(DABS(X(IAUS))-DABS(W2)*1.D-10)35,37,37                         DTEA1640
   35 IF(X(IAUS))36,37,36                                               DTEA1650
C                                                                       DTEA1660
   36 NEW=IAUS-1                                                        DTEA1670
      ISW2=0                                                            DTEA1680
      GO TO 41                                                          DTEA1690
C                                                                       DTEA1700
   37 IF(W2-1.D38)39,38,38                                              DTEA1710
   38 X(IAUS)=1.D38                                                     DTEA1720
      IMIN=I                                                            DTEA1730
C                                                                       DTEA1740
   39 W1=W4                                                             DTEA1750
      T=W2                                                              DTEA1760
      W2=W5                                                             DTEA1770
      W3=W6                                                             DTEA1780
      ISW1=ISW2                                                         DTEA1790
   40 ISW2=0                                                            DTEA1800
C                                                                       DTEA1810
      NEW=NEW-IMIN                                                      DTEA1820
C                                                                       DTEA1830
C        TEST FOR ACCURACY                                              DTEA1840
C                                                                       DTEA1850
   41 IEND=NEW-1                                                        DTEA1860
      DO 47 I=1,IEND                                                    DTEA1870
      HE1=DABS(X(I)-X(I+1))                                             DTEA1880
      HE2=DABS(X(I+1))                                                  DTEA1890
      IF(HE1-EPS)44,44,42                                               DTEA1900
   42 IF(HE2-1.)46,46,43                                                DTEA1910
   43 IF(HE1-EPS*HE2)44,44,46                                           DTEA1920
   44 ISW2=ISW2+1                                                       DTEA1930
      IF(3-ISW2)45,45,47                                                DTEA1940
   45 FIN=X(I)                                                          DTEA1950
      IER=0                                                             DTEA1960
      RETURN                                                            DTEA1970
C                                                                       DTEA1980
   46 ISW2=0                                                            DTEA1990
   47 CONTINUE                                                          DTEA2000
C                                                                       DTEA2010
      IF(NEW-6)48,2,2                                                   DTEA2020
   48 FIN=X(NEW)                                                        DTEA2030
      IER=1                                                             DTEA2040
      RETURN                                                            DTEA2050
      END                                                               DTEA2060
