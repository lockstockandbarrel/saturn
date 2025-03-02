C                                                                       TEAS  10
C     ..................................................................TEAS  20
C                                                                       TEAS  30
C        SUBROUTINE TEAS                                                TEAS  40
C                                                                       TEAS  50
C        PURPOSE                                                        TEAS  60
C           CALCULATE THE LIMIT OF A GIVEN SEQUENCE BY MEANS OF THE     TEAS  70
C           EPSILON-ALGORITHM.                                          TEAS  80
C                                                                       TEAS  90
C        USAGE                                                          TEAS 100
C           CALL TEAS(X,N,FIN,EPS,IER)                                  TEAS 110
C                                                                       TEAS 120
C        DESCRIPTION OF PARAMETERS                                      TEAS 130
C           X      - VECTOR WHOSE COMPONENTS ARE TERMS OF THE GIVEN     TEAS 140
C                    SEQUENCE. ON RETURN THE COMPONENTS OF VECTOR X     TEAS 150
C                    ARE DESTROYED.                                     TEAS 160
C           N      - DIMENSION OF INPUT VECTOR X.                       TEAS 170
C           FIN    - RESULTANT SCALAR CONTAINING ON RETURN THE LIMIT    TEAS 180
C                    OF THE GIVEN SEQUENCE.                             TEAS 190
C           EPS    - AN INPUT VALUE, WHICH SPECIFIES THE UPPER BOUND    TEAS 200
C                    OF THE RELATIVE (ABSOLUTE) ERROR IF THE COMPONENTS TEAS 210
C                    OF X ARE ABSOLUTELY GREATER (LESS) THAN ONE.       TEAS 220
C                    CALCULATION IS TERMINATED AS SOON AS THREE TIMES INTEAS 230
C                    SUCCESSION THE RELATIVE (ABSOLUTE) DIFFERENCE      TEAS 240
C                    BETWEEN NEIGHBOURING TERMS IS NOT GREATER THAN EPS.TEAS 250
C           IER    - RESULTANT ERROR PARAMETER CODED IN THE FOLLOWING   TEAS 260
C                    FORM                                               TEAS 270
C                     IER=0  - NO ERROR                                 TEAS 280
C                     IER=1  - REQUIRED ACCURACY NOT REACHED WITH       TEAS 290
C                              MAXIMAL NUMBER OF ITERATIONS             TEAS 300
C                     IER=-1 - INTEGER N IS LESS THAN TEN.              TEAS 310
C                                                                       TEAS 320
C        REMARKS                                                        TEAS 330
C           NO ACTION BESIDES ERROR MESSAGE IN CASE N LESS THAN TEN.    TEAS 340
C           THE CHARACTER OF THE GIVEN INFINITE SEQUENCE MUST BE        TEAS 350
C           RECOGNIZABLE BY THOSE N COMPONENTS OF THE INPUT VECTOR X.   TEAS 360
C                                                                       TEAS 370
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  TEAS 380
C           NONE                                                        TEAS 390
C                                                                       TEAS 400
C        METHOD                                                         TEAS 410
C           THE CONVERGENCE OF THE GIVEN SEQUENCE IS ACCELERATED BY     TEAS 420
C           MEANS OF THE E(2)-TRANSFORMATION, USED IN AN ITERATIVE WAY. TEAS 430
C           FOR REFERENCE, SEE                                          TEAS 440
C           ALGORITHM 215,SHANKS, CACM 1963, NO. 11, PP. 662. AND       TEAS 450
C           P. WYNN, SINGULAR RULES FOR CERTAIN NON-LINEAR ALGORITHMS   TEAS 460
C           BIT VOL. 3, 1963, PP. 175-195.                              TEAS 470
C                                                                       TEAS 480
C     ..................................................................TEAS 490
C                                                                       TEAS 500
      SUBROUTINE TEAS(X,N,FIN,EPS,IER)                                  TEAS 510
C                                                                       TEAS 520
      DIMENSION X(1)                                                    TEAS 530
C                                                                       TEAS 540
C        TEST ON WRONG INPUT PARAMETER N                                TEAS 550
C                                                                       TEAS 560
      NEW=N                                                             TEAS 570
      IF(NEW-10)1,2,2                                                   TEAS 580
    1 IER=-1                                                            TEAS 590
      RETURN                                                            TEAS 600
C                                                                       TEAS 610
C        CALCULATE INITIAL VALUES FOR THE EPSILON ARRAY                 TEAS 620
C                                                                       TEAS 630
    2 ISW1=0                                                            TEAS 640
      ISW2=0                                                            TEAS 650
      W1=1.E38                                                          TEAS 660
      W7=X(4)-X(3)                                                      TEAS 670
      IF(W7)3,4,3                                                       TEAS 680
    3 W1=1./W7                                                          TEAS 690
C                                                                       TEAS 700
    4 W5=1.E38                                                          TEAS 710
      W7=X(2)-X(1)                                                      TEAS 720
      IF(W7)5,6,5                                                       TEAS 730
    5 W5=1./W7                                                          TEAS 740
C                                                                       TEAS 750
    6 W4=X(3)-X(2)                                                      TEAS 760
      IF(W4)9,7,9                                                       TEAS 770
    7 W4=1.E38                                                          TEAS 780
      T=X(2)                                                            TEAS 790
      W2=X(3)                                                           TEAS 800
    8 W3=1.E38                                                          TEAS 810
      GO TO 17                                                          TEAS 820
C                                                                       TEAS 830
    9 W4=1./W4                                                          TEAS 840
C                                                                       TEAS 850
      T=1.E38                                                           TEAS 860
      W7=W4-W5                                                          TEAS 870
      IF(W7)10,11,10                                                    TEAS 880
   10 T=X(2)+1./W7                                                      TEAS 890
C                                                                       TEAS 900
   11 W2=W1-W4                                                          TEAS 910
      IF(W2)15,12,15                                                    TEAS 920
   12 W2=1.E38                                                          TEAS 930
      IF(T-1.E38)13,14,14                                               TEAS 940
   13 ISW2=1                                                            TEAS 950
   14 W3=W4                                                             TEAS 960
      GO TO 17                                                          TEAS 970
C                                                                       TEAS 980
   15 W2=X(3)+1./W2                                                     TEAS 990
      W7=W2-T                                                           TEAS1000
      IF(W7)16,8,16                                                     TEAS1010
   16 W3=W4+1./W7                                                       TEAS1020
C                                                                       TEAS1030
   17 ISW1=ISW2                                                         TEAS1040
      ISW2=0                                                            TEAS1050
      IMIN=4                                                            TEAS1060
C                                                                       TEAS1070
C        CALCULATE DIAGONALS OF THE EPSILON ARRAY IN A DO-LOOP          TEAS1080
C                                                                       TEAS1090
      DO 40 I=5,NEW                                                     TEAS1100
      IAUS=I-IMIN                                                       TEAS1110
      W4=1.E38                                                          TEAS1120
      W5=X(I-1)                                                         TEAS1130
      W7=X(I)-X(I-1)                                                    TEAS1140
      IF(W7)18,24,18                                                    TEAS1150
   18 W4=1./W7                                                          TEAS1160
C                                                                       TEAS1170
      IF(W1-1.E38)19,25,25                                              TEAS1180
   19 W6=W4-W1                                                          TEAS1190
C                                                                       TEAS1200
C        TEST FOR NECESSITY OF A SINGULAR RULE                          TEAS1210
C                                                                       TEAS1220
      IF(ABS(W6)-ABS(W4)*1.E-4)20,20,22                                 TEAS1230
   20 ISW2=1                                                            TEAS1240
      IF(W6)22,21,22                                                    TEAS1250
   21 W5=1.E38                                                          TEAS1260
      W6=W1                                                             TEAS1270
      IF(W2-1.E38)28,26,26                                              TEAS1280
   22 W5=X(I-1)+1./W6                                                   TEAS1290
C                                                                       TEAS1300
C        FIRST TEST FOR LOSS OF SIGNIFICANCE                            TEAS1310
C                                                                       TEAS1320
      IF(ABS(W5)-ABS(X(I-1))*1.E-5)23,24,24                             TEAS1330
   23 IF(W5)36,24,36                                                    TEAS1340
C                                                                       TEAS1350
   24 W7=W5-W2                                                          TEAS1360
      IF(W7)27,25,27                                                    TEAS1370
   25 W6=1.E38                                                          TEAS1380
   26 ISW2=0                                                            TEAS1390
      X(IAUS)=W2                                                        TEAS1400
      GO TO 37                                                          TEAS1410
   27 W6=W1+1./W7                                                       TEAS1420
   28 IF(ISW1-1)33,29,29                                                TEAS1430
C                                                                       TEAS1440
C        CALCULATE X(IAUS) WITH HELP OF SINGULAR RULE                   TEAS1450
C                                                                       TEAS1460
   29 IF(W2-1.E38)30,32,32                                              TEAS1470
   30 W7=W5/(W2-W5)+T/(W2-T)+X(I-2)/(X(I-2)-W2)                         TEAS1480
      IF(1.+W7)31,38,31                                                 TEAS1490
   31 X(IAUS)=W7*W2/(1.+W7)                                             TEAS1500
      GO TO 39                                                          TEAS1510
C                                                                       TEAS1520
   32 X(IAUS)=W5+T-X(I-2)                                               TEAS1530
      GO TO 39                                                          TEAS1540
C                                                                       TEAS1550
   33 W7=W6-W3                                                          TEAS1560
      IF(W7)34,38,34                                                    TEAS1570
   34 X(IAUS)=W2+1./W7                                                  TEAS1580
C                                                                       TEAS1590
C        SECOND TEST FOR LOSS OF SIGNIFICANCE                           TEAS1600
C                                                                       TEAS1610
      IF(ABS(X(IAUS))-ABS(W2)*1.E-5)35,37,37                            TEAS1620
   35 IF(X(IAUS))36,37,36                                               TEAS1630
C                                                                       TEAS1640
   36 NEW=IAUS-1                                                        TEAS1650
      ISW2=0                                                            TEAS1660
      GO TO 41                                                          TEAS1670
C                                                                       TEAS1680
   37 IF(W2-1.E38)39,38,38                                              TEAS1690
   38 X(IAUS)=1.E38                                                     TEAS1700
      IMIN=I                                                            TEAS1710
C                                                                       TEAS1720
   39 W1=W4                                                             TEAS1730
      T=W2                                                              TEAS1740
      W2=W5                                                             TEAS1750
      W3=W6                                                             TEAS1760
      ISW1=ISW2                                                         TEAS1770
   40 ISW2=0                                                            TEAS1780
C                                                                       TEAS1790
      NEW=NEW-IMIN                                                      TEAS1800
C                                                                       TEAS1810
C        TEST FOR ACCURACY                                              TEAS1820
C                                                                       TEAS1830
   41 IEND=NEW-1                                                        TEAS1840
      DO 47 I=1,IEND                                                    TEAS1850
      W1=ABS(X(I)-X(I+1))                                               TEAS1860
      W2=ABS(X(I+1))                                                    TEAS1870
      IF(W1-EPS)44,44,42                                                TEAS1880
   42 IF(W2-1.)46,46,43                                                 TEAS1890
   43 IF(W1-EPS*W2)44,44,46                                             TEAS1900
   44 ISW2=ISW2+1                                                       TEAS1910
      IF(3-ISW2)45,45,47                                                TEAS1920
   45 FIN=X(I)                                                          TEAS1930
      IER=0                                                             TEAS1940
      RETURN                                                            TEAS1950
C                                                                       TEAS1960
   46 ISW2=0                                                            TEAS1970
   47 CONTINUE                                                          TEAS1980
C                                                                       TEAS1990
      IF(NEW-6)48,2,2                                                   TEAS2000
   48 FIN=X(NEW)                                                        TEAS2010
      IER=1                                                             TEAS2020
      RETURN                                                            TEAS2030
      END                                                               TEAS2040
