C                                                                       DELS  10
C     ..................................................................DELS  20
C                                                                       DELS  30
C        SUBROUTINE DGELS                                               DELS  40
C                                                                       DELS  50
C        PURPOSE                                                        DELS  60
C           TO SOLVE A SYSTEM OF SIMULTANEOUS LINEAR EQUATIONS WITH     DELS  70
C           SYMMETRIC COEFFICIENT MATRIX UPPER TRIANGULAR PART OF WHICH DELS  80
C           IS ASSUMED TO BE STORED COLUMNWISE.                         DELS  90
C                                                                       DELS 100
C        USAGE                                                          DELS 110
C           CALL DGELS(R,A,M,N,EPS,IER,AUX)                             DELS 120
C                                                                       DELS 130
C        DESCRIPTION OF PARAMETERS                                      DELS 140
C           R      - DOUBLE PRECISION M BY N RIGHT HAND SIDE MATRIX     DELS 150
C                    (DESTROYED). ON RETURN R CONTAINS THE SOLUTION OF  DELS 160
C                    THE EQUATIONS.                                     DELS 170
C           A      - UPPER TRIANGULAR PART OF THE SYMMETRIC DOUBLE      DELS 180
C                    PRECISION M BY M COEFFICIENT MATRIX.  (DESTROYED)  DELS 190
C           M      - THE NUMBER OF EQUATIONS IN THE SYSTEM.             DELS 200
C           N      - THE NUMBER OF RIGHT HAND SIDE VECTORS.             DELS 210
C           EPS    - SINGLE PRECISION INPUT CONSTANT WHICH IS USED AS   DELS 220
C                    RELATIVE TOLERANCE FOR TEST ON LOSS OF             DELS 230
C                    SIGNIFICANCE.                                      DELS 240
C           IER    - RESULTING ERROR PARAMETER CODED AS FOLLOWS         DELS 250
C                    IER=0  - NO ERROR,                                 DELS 260
C                    IER=-1 - NO RESULT BECAUSE OF M LESS THAN 1 OR     DELS 270
C                             PIVOT ELEMENT AT ANY ELIMINATION STEP     DELS 280
C                             EQUAL TO 0,                               DELS 290
C                    IER=K  - WARNING DUE TO POSSIBLE LOSS OF SIGNIFI-  DELS 300
C                             CANCE INDICATED AT ELIMINATION STEP K+1,  DELS 310
C                             WHERE PIVOT ELEMENT WAS LESS THAN OR      DELS 320
C                             EQUAL TO THE INTERNAL TOLERANCE EPS TIMES DELS 330
C                             ABSOLUTELY GREATEST MAIN DIAGONAL         DELS 340
C                             ELEMENT OF MATRIX A.                      DELS 350
C           AUX    - DOUBLE PRECISION AUXILIARY STORAGE ARRAY           DELS 360
C                    WITH DIMENSION M-1.                                DELS 370
C                                                                       DELS 380
C        REMARKS                                                        DELS 390
C           UPPER TRIANGULAR PART OF MATRIX A IS ASSUMED TO BE STORED   DELS 400
C           COLUMNWISE IN M*(M+1)/2 SUCCESSIVE STORAGE LOCATIONS, RIGHT DELS 410
C           HAND SIDE MATRIX R COLUMNWISE IN N*M SUCCESSIVE STORAGE     DELS 420
C           LOCATIONS. ON RETURN SOLUTION MATRIX R IS STORED COLUMNWISE DELS 430
C           TOO.                                                        DELS 440
C           THE PROCEDURE GIVES RESULTS IF THE NUMBER OF EQUATIONS M IS DELS 450
C           GREATER THAN 0 AND PIVOT ELEMENTS AT ALL ELIMINATION STEPS  DELS 460
C           ARE DIFFERENT FROM 0. HOWEVER WARNING IER=K - IF GIVEN -    DELS 470
C           INDICATES POSSIBLE LOSS OF SIGNIFICANCE. IN CASE OF A WELL  DELS 480
C           SCALED MATRIX A AND APPROPRIATE TOLERANCE EPS, IER=K MAY BE DELS 490
C           INTERPRETED THAT MATRIX A HAS THE RANK K. NO WARNING IS     DELS 500
C           GIVEN IN CASE M=1.                                          DELS 510
C           ERROR PARAMETER IER=-1 DOES NOT NECESSARILY MEAN THAT       DELS 520
C           MATRIX A IS SINGULAR, AS ONLY MAIN DIAGONAL ELEMENTS        DELS 530
C           ARE USED AS PIVOT ELEMENTS. POSSIBLY SUBROUTINE DGELG (WHICHDELS 540
C           WORKS WITH TOTAL PIVOTING) WOULD BE ABLE TO FIND A SOLUTION.DELS 550
C                                                                       DELS 560
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DELS 570
C           NONE                                                        DELS 580
C                                                                       DELS 590
C        METHOD                                                         DELS 600
C           SOLUTION IS DONE BY MEANS OF GAUSS-ELIMINATION WITH         DELS 610
C           PIVOTING IN MAIN DIAGONAL, IN ORDER TO PRESERVE             DELS 620
C           SYMMETRY IN REMAINING COEFFICIENT MATRICES.                 DELS 630
C                                                                       DELS 640
C     ..................................................................DELS 650
C                                                                       DELS 660
      SUBROUTINE DGELS(R,A,M,N,EPS,IER,AUX)                             DELS 670
C                                                                       DELS 680
C                                                                       DELS 690
      DIMENSION A(1),R(1),AUX(1)                                        DELS 700
      DOUBLE PRECISION R,A,AUX,PIV,TB,TOL,PIVI                          DELS 710
      IF(M)24,24,1                                                      DELS 720
C                                                                       DELS 730
C     SEARCH FOR GREATEST MAIN DIAGONAL ELEMENT                         DELS 740
    1 IER=0                                                             DELS 750
      PIV=0.D0                                                          DELS 760
      L=0                                                               DELS 770
      DO 3 K=1,M                                                        DELS 780
      L=L+K                                                             DELS 790
      TB=DABS(A(L))                                                     DELS 800
      IF(TB-PIV)3,3,2                                                   DELS 810
    2 PIV=TB                                                            DELS 820
      I=L                                                               DELS 830
      J=K                                                               DELS 840
    3 CONTINUE                                                          DELS 850
      TOL=EPS*PIV                                                       DELS 860
C     MAIN DIAGONAL ELEMENT A(I)=A(J,J) IS FIRST PIVOT ELEMENT.         DELS 870
C     PIV CONTAINS THE ABSOLUTE VALUE OF A(I).                          DELS 880
C                                                                       DELS 890
C                                                                       DELS 900
C     START ELIMINATION LOOP                                            DELS 910
      LST=0                                                             DELS 920
      NM=N*M                                                            DELS 930
      LEND=M-1                                                          DELS 940
      DO 18 K=1,M                                                       DELS 950
C                                                                       DELS 960
C     TEST ON USEFULNESS OF SYMMETRIC ALGORITHM                         DELS 970
      IF(PIV)24,24,4                                                    DELS 980
    4 IF(IER)7,5,7                                                      DELS 990
    5 IF(PIV-TOL)6,6,7                                                  DELS1000
    6 IER=K-1                                                           DELS1010
    7 LT=J-K                                                            DELS1020
      LST=LST+K                                                         DELS1030
C                                                                       DELS1040
C     PIVOT ROW REDUCTION AND ROW INTERCHANGE IN RIGHT HAND SIDE R      DELS1050
      PIVI=1.D0/A(I)                                                    DELS1060
      DO 8 L=K,NM,M                                                     DELS1070
      LL=L+LT                                                           DELS1080
      TB=PIVI*R(LL)                                                     DELS1090
      R(LL)=R(L)                                                        DELS1100
    8 R(L)=TB                                                           DELS1110
C                                                                       DELS1120
C     IS ELIMINATION TERMINATED                                         DELS1130
      IF(K-M)9,19,19                                                    DELS1140
C                                                                       DELS1150
C     ROW AND COLUMN INTERCHANGE AND PIVOT ROW REDUCTION IN MATRIX A.   DELS1160
C     ELEMENTS OF PIVOT COLUMN ARE SAVED IN AUXILIARY VECTOR AUX.       DELS1170
    9 LR=LST+(LT*(K+J-1))/2                                             DELS1180
      LL=LR                                                             DELS1190
      L=LST                                                             DELS1200
      DO 14 II=K,LEND                                                   DELS1210
      L=L+II                                                            DELS1220
      LL=LL+1                                                           DELS1230
      IF(L-LR)12,10,11                                                  DELS1240
   10 A(LL)=A(LST)                                                      DELS1250
      TB=A(L)                                                           DELS1260
      GO TO 13                                                          DELS1270
   11 LL=L+LT                                                           DELS1280
   12 TB=A(LL)                                                          DELS1290
      A(LL)=A(L)                                                        DELS1300
   13 AUX(II)=TB                                                        DELS1310
   14 A(L)=PIVI*TB                                                      DELS1320
C                                                                       DELS1330
C     SAVE COLUMN INTERCHANGE INFORMATION                               DELS1340
      A(LST)=LT                                                         DELS1350
C                                                                       DELS1360
C     ELEMENT REDUCTION AND SEARCH FOR NEXT PIVOT                       DELS1370
      PIV=0.D0                                                          DELS1380
      LLST=LST                                                          DELS1390
      LT=0                                                              DELS1400
      DO 18 II=K,LEND                                                   DELS1410
      PIVI=-AUX(II)                                                     DELS1420
      LL=LLST                                                           DELS1430
      LT=LT+1                                                           DELS1440
      DO 15 LLD=II,LEND                                                 DELS1450
      LL=LL+LLD                                                         DELS1460
      L=LL+LT                                                           DELS1470
   15 A(L)=A(L)+PIVI*A(LL)                                              DELS1480
      LLST=LLST+II                                                      DELS1490
      LR=LLST+LT                                                        DELS1500
      TB=DABS(A(LR))                                                    DELS1510
      IF(TB-PIV)17,17,16                                                DELS1520
   16 PIV=TB                                                            DELS1530
      I=LR                                                              DELS1540
      J=II+1                                                            DELS1550
   17 DO 18 LR=K,NM,M                                                   DELS1560
      LL=LR+LT                                                          DELS1570
   18 R(LL)=R(LL)+PIVI*R(LR)                                            DELS1580
C     END OF ELIMINATION LOOP                                           DELS1590
C                                                                       DELS1600
C                                                                       DELS1610
C     BACK SUBSTITUTION AND BACK INTERCHANGE                            DELS1620
   19 IF(LEND)24,23,20                                                  DELS1630
   20 II=M                                                              DELS1640
      DO 22 I=2,M                                                       DELS1650
      LST=LST-II                                                        DELS1660
      II=II-1                                                           DELS1670
      L=A(LST)+.5D0                                                     DELS1680
      DO 22 J=II,NM,M                                                   DELS1690
      TB=R(J)                                                           DELS1700
      LL=J                                                              DELS1710
      K=LST                                                             DELS1720
      DO 21 LT=II,LEND                                                  DELS1730
      LL=LL+1                                                           DELS1740
      K=K+LT                                                            DELS1750
   21 TB=TB-A(K)*R(LL)                                                  DELS1760
      K=J+L                                                             DELS1770
      R(J)=R(K)                                                         DELS1780
   22 R(K)=TB                                                           DELS1790
   23 RETURN                                                            DELS1800
C                                                                       DELS1810
C                                                                       DELS1820
C     ERROR RETURN                                                      DELS1830
   24 IER=-1                                                            DELS1840
      RETURN                                                            DELS1850
      END                                                               DELS1860
