C                                                                       GELS  10
C     ..................................................................GELS  20
C                                                                       GELS  30
C        SUBROUTINE GELS                                                GELS  40
C                                                                       GELS  50
C        PURPOSE                                                        GELS  60
C           TO SOLVE A SYSTEM OF SIMULTANEOUS LINEAR EQUATIONS WITH     GELS  70
C           SYMMETRIC COEFFICIENT MATRIX UPPER TRIANGULAR PART OF WHICH GELS  80
C           IS ASSUMED TO BE STORED COLUMNWISE.                         GELS  90
C                                                                       GELS 100
C        USAGE                                                          GELS 110
C           CALL GELS(R,A,M,N,EPS,IER,AUX)                              GELS 120
C                                                                       GELS 130
C        DESCRIPTION OF PARAMETERS                                      GELS 140
C           R      - M BY N RIGHT HAND SIDE MATRIX.  (DESTROYED)        GELS 150
C                    ON RETURN R CONTAINS THE SOLUTION OF THE EQUATIONS.GELS 160
C           A      - UPPER TRIANGULAR PART OF THE SYMMETRIC             GELS 170
C                    M BY M COEFFICIENT MATRIX.  (DESTROYED)            GELS 180
C           M      - THE NUMBER OF EQUATIONS IN THE SYSTEM.             GELS 190
C           N      - THE NUMBER OF RIGHT HAND SIDE VECTORS.             GELS 200
C           EPS    - AN INPUT CONSTANT WHICH IS USED AS RELATIVE        GELS 210
C                    TOLERANCE FOR TEST ON LOSS OF SIGNIFICANCE.        GELS 220
C           IER    - RESULTING ERROR PARAMETER CODED AS FOLLOWS         GELS 230
C                    IER=0  - NO ERROR,                                 GELS 240
C                    IER=-1 - NO RESULT BECAUSE OF M LESS THAN 1 OR     GELS 250
C                             PIVOT ELEMENT AT ANY ELIMINATION STEP     GELS 260
C                             EQUAL TO 0,                               GELS 270
C                    IER=K  - WARNING DUE TO POSSIBLE LOSS OF SIGNIFI-  GELS 280
C                             CANCE INDICATED AT ELIMINATION STEP K+1,  GELS 290
C                             WHERE PIVOT ELEMENT WAS LESS THAN OR      GELS 300
C                             EQUAL TO THE INTERNAL TOLERANCE EPS TIMES GELS 310
C                             ABSOLUTELY GREATEST MAIN DIAGONAL         GELS 320
C                             ELEMENT OF MATRIX A.                      GELS 330
C           AUX    - AN AUXILIARY STORAGE ARRAY WITH DIMENSION M-1.     GELS 340
C                                                                       GELS 350
C        REMARKS                                                        GELS 360
C           UPPER TRIANGULAR PART OF MATRIX A IS ASSUMED TO BE STORED   GELS 370
C           COLUMNWISE IN M*(M+1)/2 SUCCESSIVE STORAGE LOCATIONS, RIGHT GELS 380
C           HAND SIDE MATRIX R COLUMNWISE IN N*M SUCCESSIVE STORAGE     GELS 390
C           LOCATIONS. ON RETURN SOLUTION MATRIX R IS STORED COLUMNWISE GELS 400
C           TOO.                                                        GELS 410
C           THE PROCEDURE GIVES RESULTS IF THE NUMBER OF EQUATIONS M IS GELS 420
C           GREATER THAN 0 AND PIVOT ELEMENTS AT ALL ELIMINATION STEPS  GELS 430
C           ARE DIFFERENT FROM 0. HOWEVER WARNING IER=K - IF GIVEN -    GELS 440
C           INDICATES POSSIBLE LOSS OF SIGNIFICANCE. IN CASE OF A WELL  GELS 450
C           SCALED MATRIX A AND APPROPRIATE TOLERANCE EPS, IER=K MAY BE GELS 460
C           INTERPRETED THAT MATRIX A HAS THE RANK K. NO WARNING IS     GELS 470
C           GIVEN IN CASE M=1.                                          GELS 480
C           ERROR PARAMETER IER=-1 DOES NOT NECESSARILY MEAN THAT       GELS 490
C           MATRIX A IS SINGULAR, AS ONLY MAIN DIAGONAL ELEMENTS        GELS 500
C           ARE USED AS PIVOT ELEMENTS. POSSIBLY SUBROUTINE GELG (WHICH GELS 510
C           WORKS WITH TOTAL PIVOTING) WOULD BE ABLE TO FIND A SOLUTION.GELS 520
C                                                                       GELS 530
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  GELS 540
C           NONE                                                        GELS 550
C                                                                       GELS 560
C        METHOD                                                         GELS 570
C           SOLUTION IS DONE BY MEANS OF GAUSS-ELIMINATION WITH         GELS 580
C           PIVOTING IN MAIN DIAGONAL, IN ORDER TO PRESERVE             GELS 590
C           SYMMETRY IN REMAINING COEFFICIENT MATRICES.                 GELS 600
C                                                                       GELS 610
C     ..................................................................GELS 620
C                                                                       GELS 630
      SUBROUTINE GELS(R,A,M,N,EPS,IER,AUX)                              GELS 640
C                                                                       GELS 650
C                                                                       GELS 660
      DIMENSION A(1),R(1),AUX(1)                                        GELS 670
      IF(M)24,24,1                                                      GELS 680
C                                                                       GELS 690
C     SEARCH FOR GREATEST MAIN DIAGONAL ELEMENT                         GELS 700
    1 IER=0                                                             GELS 710
      PIV=0.                                                            GELS 720
      L=0                                                               GELS 730
      DO 3 K=1,M                                                        GELS 740
      L=L+K                                                             GELS 750
      TB=ABS(A(L))                                                      GELS 760
      IF(TB-PIV)3,3,2                                                   GELS 770
    2 PIV=TB                                                            GELS 780
      I=L                                                               GELS 790
      J=K                                                               GELS 800
    3 CONTINUE                                                          GELS 810
      TOL=EPS*PIV                                                       GELS 820
C     MAIN DIAGONAL ELEMENT A(I)=A(J,J) IS FIRST PIVOT ELEMENT.         GELS 830
C     PIV CONTAINS THE ABSOLUTE VALUE OF A(I).                          GELS 840
C                                                                       GELS 850
C                                                                       GELS 860
C     START ELIMINATION LOOP                                            GELS 870
      LST=0                                                             GELS 880
      NM=N*M                                                            GELS 890
      LEND=M-1                                                          GELS 900
      DO 18 K=1,M                                                       GELS 910
C                                                                       GELS 920
C     TEST ON USEFULNESS OF SYMMETRIC ALGORITHM                         GELS 930
      IF(PIV)24,24,4                                                    GELS 940
    4 IF(IER)7,5,7                                                      GELS 950
    5 IF(PIV-TOL)6,6,7                                                  GELS 960
    6 IER=K-1                                                           GELS 970
    7 LT=J-K                                                            GELS 980
      LST=LST+K                                                         GELS 990
C                                                                       GELS1000
C     PIVOT ROW REDUCTION AND ROW INTERCHANGE IN RIGHT HAND SIDE R      GELS1010
      PIVI=1./A(I)                                                      GELS1020
      DO 8 L=K,NM,M                                                     GELS1030
      LL=L+LT                                                           GELS1040
      TB=PIVI*R(LL)                                                     GELS1050
      R(LL)=R(L)                                                        GELS1060
    8 R(L)=TB                                                           GELS1070
C                                                                       GELS1080
C     IS ELIMINATION TERMINATED                                         GELS1090
      IF(K-M)9,19,19                                                    GELS1100
C                                                                       GELS1110
C     ROW AND COLUMN INTERCHANGE AND PIVOT ROW REDUCTION IN MATRIX A.   GELS1120
C     ELEMENTS OF PIVOT COLUMN ARE SAVED IN AUXILIARY VECTOR AUX.       GELS1130
    9 LR=LST+(LT*(K+J-1))/2                                             GELS1140
      LL=LR                                                             GELS1150
      L=LST                                                             GELS1160
      DO 14 II=K,LEND                                                   GELS1170
      L=L+II                                                            GELS1180
      LL=LL+1                                                           GELS1190
      IF(L-LR)12,10,11                                                  GELS1200
   10 A(LL)=A(LST)                                                      GELS1210
      TB=A(L)                                                           GELS1220
      GO TO 13                                                          GELS1230
   11 LL=L+LT                                                           GELS1240
   12 TB=A(LL)                                                          GELS1250
      A(LL)=A(L)                                                        GELS1260
   13 AUX(II)=TB                                                        GELS1270
   14 A(L)=PIVI*TB                                                      GELS1280
C                                                                       GELS1290
C     SAVE COLUMN INTERCHANGE INFORMATION                               GELS1300
      A(LST)=LT                                                         GELS1310
C                                                                       GELS1320
C     ELEMENT REDUCTION AND SEARCH FOR NEXT PIVOT                       GELS1330
      PIV=0.                                                            GELS1340
      LLST=LST                                                          GELS1350
      LT=0                                                              GELS1360
      DO 18 II=K,LEND                                                   GELS1370
      PIVI=-AUX(II)                                                     GELS1380
      LL=LLST                                                           GELS1390
      LT=LT+1                                                           GELS1400
      DO 15 LLD=II,LEND                                                 GELS1410
      LL=LL+LLD                                                         GELS1420
      L=LL+LT                                                           GELS1430
   15 A(L)=A(L)+PIVI*A(LL)                                              GELS1440
      LLST=LLST+II                                                      GELS1450
      LR=LLST+LT                                                        GELS1460
      TB=ABS(A(LR))                                                     GELS1470
      IF(TB-PIV)17,17,16                                                GELS1480
   16 PIV=TB                                                            GELS1490
      I=LR                                                              GELS1500
      J=II+1                                                            GELS1510
   17 DO 18 LR=K,NM,M                                                   GELS1520
      LL=LR+LT                                                          GELS1530
   18 R(LL)=R(LL)+PIVI*R(LR)                                            GELS1540
C     END OF ELIMINATION LOOP                                           GELS1550
C                                                                       GELS1560
C                                                                       GELS1570
C     BACK SUBSTITUTION AND BACK INTERCHANGE                            GELS1580
   19 IF(LEND)24,23,20                                                  GELS1590
   20 II=M                                                              GELS1600
      DO 22 I=2,M                                                       GELS1610
      LST=LST-II                                                        GELS1620
      II=II-1                                                           GELS1630
      L=A(LST)+.5                                                       GELS1640
      DO 22 J=II,NM,M                                                   GELS1650
      TB=R(J)                                                           GELS1660
      LL=J                                                              GELS1670
      K=LST                                                             GELS1680
      DO 21 LT=II,LEND                                                  GELS1690
      LL=LL+1                                                           GELS1700
      K=K+LT                                                            GELS1710
   21 TB=TB-A(K)*R(LL)                                                  GELS1720
      K=J+L                                                             GELS1730
      R(J)=R(K)                                                         GELS1740
   22 R(K)=TB                                                           GELS1750
   23 RETURN                                                            GELS1760
C                                                                       GELS1770
C                                                                       GELS1780
C     ERROR RETURN                                                      GELS1790
   24 IER=-1                                                            GELS1800
      RETURN                                                            GELS1810
      END                                                               GELS1820
