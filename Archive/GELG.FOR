C                                                                       GELG  10
C     ..................................................................GELG  20
C                                                                       GELG  30
C        SUBROUTINE GELG                                                GELG  40
C                                                                       GELG  50
C        PURPOSE                                                        GELG  60
C           TO SOLVE A GENERAL SYSTEM OF SIMULTANEOUS LINEAR EQUATIONS. GELG  70
C                                                                       GELG  80
C        USAGE                                                          GELG  90
C           CALL GELG(R,A,M,N,EPS,IER)                                  GELG 100
C                                                                       GELG 110
C        DESCRIPTION OF PARAMETERS                                      GELG 120
C           R      - THE M BY N MATRIX OF RIGHT HAND SIDES.  (DESTROYED)GELG 130
C                    ON RETURN R CONTAINS THE SOLUTION OF THE EQUATIONS.GELG 140
C           A      - THE M BY M COEFFICIENT MATRIX.  (DESTROYED)        GELG 150
C           M      - THE NUMBER OF EQUATIONS IN THE SYSTEM.             GELG 160
C           N      - THE NUMBER OF RIGHT HAND SIDE VECTORS.             GELG 170
C           EPS    - AN INPUT CONSTANT WHICH IS USED AS RELATIVE        GELG 180
C                    TOLERANCE FOR TEST ON LOSS OF SIGNIFICANCE.        GELG 190
C           IER    - RESULTING ERROR PARAMETER CODED AS FOLLOWS         GELG 200
C                    IER=0  - NO ERROR,                                 GELG 210
C                    IER=-1 - NO RESULT BECAUSE OF M LESS THAN 1 OR     GELG 220
C                             PIVOT ELEMENT AT ANY ELIMINATION STEP     GELG 230
C                             EQUAL TO 0,                               GELG 240
C                    IER=K  - WARNING DUE TO POSSIBLE LOSS OF SIGNIFI-  GELG 250
C                             CANCE INDICATED AT ELIMINATION STEP K+1,  GELG 260
C                             WHERE PIVOT ELEMENT WAS LESS THAN OR      GELG 270
C                             EQUAL TO THE INTERNAL TOLERANCE EPS TIMES GELG 280
C                             ABSOLUTELY GREATEST ELEMENT OF MATRIX A.  GELG 290
C                                                                       GELG 300
C        REMARKS                                                        GELG 310
C           INPUT MATRICES R AND A ARE ASSUMED TO BE STORED COLUMNWISE  GELG 320
C           IN M*N RESP. M*M SUCCESSIVE STORAGE LOCATIONS. ON RETURN    GELG 330
C           SOLUTION MATRIX R IS STORED COLUMNWISE TOO.                 GELG 340
C           THE PROCEDURE GIVES RESULTS IF THE NUMBER OF EQUATIONS M IS GELG 350
C           GREATER THAN 0 AND PIVOT ELEMENTS AT ALL ELIMINATION STEPS  GELG 360
C           ARE DIFFERENT FROM 0. HOWEVER WARNING IER=K - IF GIVEN -    GELG 370
C           INDICATES POSSIBLE LOSS OF SIGNIFICANCE. IN CASE OF A WELL  GELG 380
C           SCALED MATRIX A AND APPROPRIATE TOLERANCE EPS, IER=K MAY BE GELG 390
C           INTERPRETED THAT MATRIX A HAS THE RANK K. NO WARNING IS     GELG 400
C           GIVEN IN CASE M=1.                                          GELG 410
C                                                                       GELG 420
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  GELG 430
C           NONE                                                        GELG 440
C                                                                       GELG 450
C        METHOD                                                         GELG 460
C           SOLUTION IS DONE BY MEANS OF GAUSS-ELIMINATION WITH         GELG 470
C           COMPLETE PIVOTING.                                          GELG 480
C                                                                       GELG 490
C     ..................................................................GELG 500
C                                                                       GELG 510
      SUBROUTINE GELG(R,A,M,N,EPS,IER)                                  GELG 520
C                                                                       GELG 530
C                                                                       GELG 540
      DIMENSION A(1),R(1)                                               GELG 550
      IF(M)23,23,1                                                      GELG 560
C                                                                       GELG 570
C     SEARCH FOR GREATEST ELEMENT IN MATRIX A                           GELG 580
    1 IER=0                                                             GELG 590
      PIV=0.                                                            GELG 600
      MM=M*M                                                            GELG 610
      NM=N*M                                                            GELG 620
      DO 3 L=1,MM                                                       GELG 630
      TB=ABS(A(L))                                                      GELG 640
      IF(TB-PIV)3,3,2                                                   GELG 650
    2 PIV=TB                                                            GELG 660
      I=L                                                               GELG 670
    3 CONTINUE                                                          GELG 680
      TOL=EPS*PIV                                                       GELG 690
C     A(I) IS PIVOT ELEMENT. PIV CONTAINS THE ABSOLUTE VALUE OF A(I).   GELG 700
C                                                                       GELG 710
C                                                                       GELG 720
C     START ELIMINATION LOOP                                            GELG 730
      LST=1                                                             GELG 740
      DO 17 K=1,M                                                       GELG 750
C                                                                       GELG 760
C     TEST ON SINGULARITY                                               GELG 770
      IF(PIV)23,23,4                                                    GELG 780
    4 IF(IER)7,5,7                                                      GELG 790
    5 IF(PIV-TOL)6,6,7                                                  GELG 800
    6 IER=K-1                                                           GELG 810
    7 PIVI=1./A(I)                                                      GELG 820
      J=(I-1)/M                                                         GELG 830
      I=I-J*M-K                                                         GELG 840
      J=J+1-K                                                           GELG 850
C     I+K IS ROW-INDEX, J+K COLUMN-INDEX OF PIVOT ELEMENT               GELG 860
C                                                                       GELG 870
C     PIVOT ROW REDUCTION AND ROW INTERCHANGE IN RIGHT HAND SIDE R      GELG 880
      DO 8 L=K,NM,M                                                     GELG 890
      LL=L+I                                                            GELG 900
      TB=PIVI*R(LL)                                                     GELG 910
      R(LL)=R(L)                                                        GELG 920
    8 R(L)=TB                                                           GELG 930
C                                                                       GELG 940
C     IS ELIMINATION TERMINATED                                         GELG 950
      IF(K-M)9,18,18                                                    GELG 960
C                                                                       GELG 970
C     COLUMN INTERCHANGE IN MATRIX A                                    GELG 980
    9 LEND=LST+M-K                                                      GELG 990
      IF(J)12,12,10                                                     GELG1000
   10 II=J*M                                                            GELG1010
      DO 11 L=LST,LEND                                                  GELG1020
      TB=A(L)                                                           GELG1030
      LL=L+II                                                           GELG1040
      A(L)=A(LL)                                                        GELG1050
   11 A(LL)=TB                                                          GELG1060
C                                                                       GELG1070
C     ROW INTERCHANGE AND PIVOT ROW REDUCTION IN MATRIX A               GELG1080
   12 DO 13 L=LST,MM,M                                                  GELG1090
      LL=L+I                                                            GELG1100
      TB=PIVI*A(LL)                                                     GELG1110
      A(LL)=A(L)                                                        GELG1120
   13 A(L)=TB                                                           GELG1130
C                                                                       GELG1140
C     SAVE COLUMN INTERCHANGE INFORMATION                               GELG1150
      A(LST)=J                                                          GELG1160
C                                                                       GELG1170
C     ELEMENT REDUCTION AND NEXT PIVOT SEARCH                           GELG1180
      PIV=0.                                                            GELG1190
      LST=LST+1                                                         GELG1200
      J=0                                                               GELG1210
      DO 16 II=LST,LEND                                                 GELG1220
      PIVI=-A(II)                                                       GELG1230
      IST=II+M                                                          GELG1240
      J=J+1                                                             GELG1250
      DO 15 L=IST,MM,M                                                  GELG1260
      LL=L-J                                                            GELG1270
      A(L)=A(L)+PIVI*A(LL)                                              GELG1280
      TB=ABS(A(L))                                                      GELG1290
      IF(TB-PIV)15,15,14                                                GELG1300
   14 PIV=TB                                                            GELG1310
      I=L                                                               GELG1320
   15 CONTINUE                                                          GELG1330
      DO 16 L=K,NM,M                                                    GELG1340
      LL=L+J                                                            GELG1350
   16 R(LL)=R(LL)+PIVI*R(L)                                             GELG1360
   17 LST=LST+M                                                         GELG1370
C     END OF ELIMINATION LOOP                                           GELG1380
C                                                                       GELG1390
C                                                                       GELG1400
C     BACK SUBSTITUTION AND BACK INTERCHANGE                            GELG1410
   18 IF(M-1)23,22,19                                                   GELG1420
   19 IST=MM+M                                                          GELG1430
      LST=M+1                                                           GELG1440
      DO 21 I=2,M                                                       GELG1450
      II=LST-I                                                          GELG1460
      IST=IST-LST                                                       GELG1470
      L=IST-M                                                           GELG1480
      L=A(L)+.5                                                         GELG1490
      DO 21 J=II,NM,M                                                   GELG1500
      TB=R(J)                                                           GELG1510
      LL=J                                                              GELG1520
      DO 20 K=IST,MM,M                                                  GELG1530
      LL=LL+1                                                           GELG1540
   20 TB=TB-A(K)*R(LL)                                                  GELG1550
      K=J+L                                                             GELG1560
      R(J)=R(K)                                                         GELG1570
   21 R(K)=TB                                                           GELG1580
   22 RETURN                                                            GELG1590
C                                                                       GELG1600
C                                                                       GELG1610
C     ERROR RETURN                                                      GELG1620
   23 IER=-1                                                            GELG1630
      RETURN                                                            GELG1640
      END                                                               GELG1650
