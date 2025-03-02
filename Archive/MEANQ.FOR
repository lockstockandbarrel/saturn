C                                                                       MEAN  10
C     ..................................................................MEAN  20
C                                                                       MEAN  30
C        SUBROUTINE MEANQ                                               MEAN  40
C                                                                       MEAN  50
C        PURPOSE                                                        MEAN  60
C           COMPUTE SUM OF SQUARES, DEGREES OF FREEDOM, AND MEAN SQUARE MEAN  70
C           USING THE MEAN SQUARE OPERATOR.  THIS SUBROUTINE NORMALLY   MEAN  80
C           FOLLOWS CALLS TO AVDAT AND AVCAL SUBROUTINES IN THE PER-    MEAN  90
C           FORMANCE OF ANALYSIS OF VARIANCE FOR A COMPLETE FACTORIAL   MEAN 100
C           DESIGN.                                                     MEAN 110
C                                                                       MEAN 120
C        USAGE                                                          MEAN 130
C           CALL MEANQ (K,LEVEL,X,GMEAN,SUMSQ,NDF,SMEAN,MSTEP,KOUNT,    MEAN 140
C                        LASTS)                                         MEAN 150
C                                                                       MEAN 160
C        DESCRIPTION OF PARAMETERS                                      MEAN 170
C           K     - NUMBER OF VARIABLES (FACTORS). K MUST BE .GT. ONE.  MEAN 180
C           LEVEL - INPUT VECTOR OF LENGTH K CONTAINING LEVELS (CATE-   MEAN 190
C                   GORIES) WITHIN EACH VARIABLE.                       MEAN 200
C           X     - INPUT VECTOR CONTAINING THE RESULT OF THE SIGMA AND MEAN 210
C                   DELTA OPERATORS. THE LENGTH OF X IS                 MEAN 220
C                   (LEVEL(1)+1)*(LEVEL(2)+1)*...*(LEVEL(K)+1).         MEAN 230
C           GMEAN - OUTPUT VARIABLE CONTAINING GRAND MEAN.              MEAN 240
C           SUMSQ - OUTPUT VECTOR CONTAINING SUMS OF SQUARES.  THE      MEAN 250
C                   LENGTH OF SUMSQ IS 2 TO THE K-TH POWER MINUS ONE,   MEAN 260
C                   (2**K)-1.                                           MEAN 270
C           NDF   - OUTPUT VECTOR CONTAINING DEGREES OF FREEDOM.  THE   MEAN 280
C                   LENGTH OF NDF IS 2 TO THE K-TH POWER MINUS ONE,     MEAN 290
C                   (2**K)-1.                                           MEAN 300
C           SMEAN - OUTPUT VECTOR CONTAINING MEAN SQUARES.  THE         MEAN 310
C                   LENGTH OF SMEAN IS 2 TO THE K-TH POWER MINUS ONE,   MEAN 320
C                   (2**K)-1.                                           MEAN 330
C           MSTEP - WORKING VECTOR OF LENGTH K.                         MEAN 340
C           KOUNT - WORKING VECTOR OF LENGTH K.                         MEAN 350
C           LASTS - WORKING VECTOR OF LENGTH K.                         MEAN 360
C                                                                       MEAN 370
C        REMARKS                                                        MEAN 380
C           THIS SUBROUTINE MUST FOLLOW SUBROUTINE AVCAL                MEAN 390
C                                                                       MEAN 400
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  MEAN 410
C           NONE                                                        MEAN 420
C                                                                       MEAN 430
C        METHOD                                                         MEAN 440
C           THE METHOD IS BASED ON THE TECHNIQUE DISCUSSED BY H. O.     MEAN 450
C           HARTLEY IN 'MATHEMATICAL METHODS FOR DIGITAL COMPUTERS',    MEAN 460
C           EDITED BY A. RALSTON AND H. WILF, JOHN WILEY AND SONS,      MEAN 470
C           1962, CHAPTER 20.                                           MEAN 480
C                                                                       MEAN 490
C     ..................................................................MEAN 500
C                                                                       MEAN 510
      SUBROUTINE MEANQ (K,LEVEL,X,GMEAN,SUMSQ,NDF,SMEAN,MSTEP,KOUNT,    MEAN 520
     1                  LASTS)                                          MEAN 530
      DIMENSION LEVEL(1),X(1),SUMSQ(1),NDF(1),SMEAN(1),MSTEP(1),        MEAN 540
     1          KOUNT(1),LASTS(1)                                       MEAN 550
C                                                                       MEAN 560
C        ...............................................................MEAN 570
C                                                                       MEAN 580
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE  MEAN 590
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION      MEAN 600
C        STATEMENT WHICH FOLLOWS.                                       MEAN 610
C                                                                       MEAN 620
C     DOUBLE PRECISION X,GMEAN,SUMSQ,SMEAN,FN1                          MEAN 630
C                                                                       MEAN 640
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS    MEAN 650
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS      MEAN 660
C        ROUTINE.                                                       MEAN 670
C                                                                       MEAN 680
C        ...............................................................MEAN 690
C                                                                       MEAN 700
C     CALCULATE TOTAL NUMBER OF DATA                                    MEAN 710
C                                                                       MEAN 720
      N=LEVEL(1)                                                        MEAN 730
      DO 150 I=2,K                                                      MEAN 740
  150 N=N*LEVEL(I)                                                      MEAN 750
C                                                                       MEAN 760
C     SET UP CONTROL FOR MEAN SQUARE OPERATOR                           MEAN 770
C                                                                       MEAN 780
      LASTS(1)=LEVEL(1)                                                 MEAN 790
      DO 178 I=2,K                                                      MEAN 800
  178 LASTS(I)=LEVEL(I)+1                                               MEAN 810
      NN=1                                                              MEAN 820
C                                                                       MEAN 830
C     CLEAR THE AREA TO STORE SUMS OF SQUARES                           MEAN 840
C                                                                       MEAN 850
      LL=(2**K)-1                                                       MEAN 860
      MSTEP(1)=1                                                        MEAN 870
      DO 180 I=2,K                                                      MEAN 880
  180 MSTEP(I)=MSTEP(I-1)*2                                             MEAN 890
      DO 185 I=1,LL                                                     MEAN 900
  185 SUMSQ(I)=0.0                                                      MEAN 910
C                                                                       MEAN 920
C     PERFORM MEAN SQUARE OPERATOR                                      MEAN 930
C                                                                       MEAN 940
      DO 190 I=1,K                                                      MEAN 950
  190 KOUNT(I)=0                                                        MEAN 960
  200 L=0                                                               MEAN 970
      DO 260 I=1,K                                                      MEAN 980
      IF(KOUNT(I)-LASTS(I)) 210, 250, 210                               MEAN 990
  210 IF(L) 220, 220, 240                                               MEAN1000
  220 KOUNT(I)=KOUNT(I)+1                                               MEAN1010
      IF(KOUNT(I)-LEVEL(I)) 230, 230, 250                               MEAN1020
  230 L=L+MSTEP(I)                                                      MEAN1030
      GO TO 260                                                         MEAN1040
  240 IF(KOUNT(I)-LEVEL(I)) 230, 260, 230                               MEAN1050
  250 KOUNT(I)=0                                                        MEAN1060
  260 CONTINUE                                                          MEAN1070
      IF(L) 285, 285, 270                                               MEAN1080
  270 SUMSQ(L)=SUMSQ(L)+X(NN)*X(NN)                                     MEAN1090
      NN=NN+1                                                           MEAN1100
      GO TO 200                                                         MEAN1110
C                                                                       MEAN1120
C     CALCULATE THE GRAND MEAN                                          MEAN1130
C                                                                       MEAN1140
  285 FN=N                                                              MEAN1150
      GMEAN=X(NN)/FN                                                    MEAN1160
C                                                                       MEAN1170
C     CALCULATE FIRST DIVISOR REQUIRED TO FORM SUM OF SQUARES AND SECONDMEAN1180
C     DIVISOR, WHICH IS EQUAL TO DEGREES OF FREEDOM, REQUIRED TO FORM   MEAN1190
C     MEAN SQUARES                                                      MEAN1200
C                                                                       MEAN1210
      DO 310 I=2,K                                                      MEAN1220
  310 MSTEP(I)=0                                                        MEAN1230
      NN=0                                                              MEAN1240
      MSTEP(1)=1                                                        MEAN1250
  320 ND1=1                                                             MEAN1260
      ND2=1                                                             MEAN1270
      DO 340 I=1,K                                                      MEAN1280
      IF(MSTEP(I)) 330, 340, 330                                        MEAN1290
  330 ND1=ND1*LEVEL(I)                                                  MEAN1300
      ND2=ND2*(LEVEL(I)-1)                                              MEAN1310
  340 CONTINUE                                                          MEAN1320
      FN1=N*ND1                                                         MEAN1330
      FN2=ND2                                                           MEAN1340
      NN=NN+1                                                           MEAN1350
      SUMSQ(NN)=SUMSQ(NN)/FN1                                           MEAN1360
      NDF(NN)=ND2                                                       MEAN1370
      SMEAN(NN)=SUMSQ(NN)/FN2                                           MEAN1380
      IF(NN-LL) 345, 370, 370                                           MEAN1390
  345 DO 360 I=1,K                                                      MEAN1400
      IF(MSTEP(I)) 347, 350, 347                                        MEAN1410
  347 MSTEP(I)=0                                                        MEAN1420
      GO TO 360                                                         MEAN1430
  350 MSTEP(I)=1                                                        MEAN1440
      GO TO 320                                                         MEAN1450
  360 CONTINUE                                                          MEAN1460
  370 RETURN                                                            MEAN1470
      END                                                               MEAN1480
