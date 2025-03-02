C                                                                       DELG  10
C     ..................................................................DELG  20
C                                                                       DELG  30
C        SUBROUTINE DGELG                                               DELG  40
C                                                                       DELG  50
C        PURPOSE                                                        DELG  60
C           TO SOLVE A GENERAL SYSTEM OF SIMULTANEOUS LINEAR EQUATIONS. DELG  70
C                                                                       DELG  80
C        USAGE                                                          DELG  90
C           CALL DGELG(R,A,M,N,EPS,IER)                                 DELG 100
C                                                                       DELG 110
C        DESCRIPTION OF PARAMETERS                                      DELG 120
C           R      - DOUBLE PRECISION M BY N RIGHT HAND SIDE MATRIX     DELG 130
C                    (DESTROYED). ON RETURN R CONTAINS THE SOLUTIONS    DELG 140
C                    OF THE EQUATIONS.                                  DELG 150
C           A      - DOUBLE PRECISION M BY M COEFFICIENT MATRIX         DELG 160
C                    (DESTROYED).                                       DELG 170
C           M      - THE NUMBER OF EQUATIONS IN THE SYSTEM.             DELG 180
C           N      - THE NUMBER OF RIGHT HAND SIDE VECTORS.             DELG 190
C           EPS    - SINGLE PRECISION INPUT CONSTANT WHICH IS USED AS   DELG 200
C                    RELATIVE TOLERANCE FOR TEST ON LOSS OF             DELG 210
C                    SIGNIFICANCE.                                      DELG 220
C           IER    - RESULTING ERROR PARAMETER CODED AS FOLLOWS         DELG 230
C                    IER=0  - NO ERROR,                                 DELG 240
C                    IER=-1 - NO RESULT BECAUSE OF M LESS THAN 1 OR     DELG 250
C                             PIVOT ELEMENT AT ANY ELIMINATION STEP     DELG 260
C                             EQUAL TO 0,                               DELG 270
C                    IER=K  - WARNING DUE TO POSSIBLE LOSS OF SIGNIFI-  DELG 280
C                             CANCE INDICATED AT ELIMINATION STEP K+1,  DELG 290
C                             WHERE PIVOT ELEMENT WAS LESS THAN OR      DELG 300
C                             EQUAL TO THE INTERNAL TOLERANCE EPS TIMES DELG 310
C                             ABSOLUTELY GREATEST ELEMENT OF MATRIX A.  DELG 320
C                                                                       DELG 330
C        REMARKS                                                        DELG 340
C           INPUT MATRICES R AND A ARE ASSUMED TO BE STORED COLUMNWISE  DELG 350
C           IN M*N RESP. M*M SUCCESSIVE STORAGE LOCATIONS. ON RETURN    DELG 360
C           SOLUTION MATRIX R IS STORED COLUMNWISE TOO.                 DELG 370
C           THE PROCEDURE GIVES RESULTS IF THE NUMBER OF EQUATIONS M IS DELG 380
C           GREATER THAN 0 AND PIVOT ELEMENTS AT ALL ELIMINATION STEPS  DELG 390
C           ARE DIFFERENT FROM 0. HOWEVER WARNING IER=K - IF GIVEN -    DELG 400
C           INDICATES POSSIBLE LOSS OF SIGNIFICANCE. IN CASE OF A WELL  DELG 410
C           SCALED MATRIX A AND APPROPRIATE TOLERANCE EPS, IER=K MAY BE DELG 420
C           INTERPRETED THAT MATRIX A HAS THE RANK K. NO WARNING IS     DELG 430
C           GIVEN IN CASE M=1.                                          DELG 440
C                                                                       DELG 450
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DELG 460
C           NONE                                                        DELG 470
C                                                                       DELG 480
C        METHOD                                                         DELG 490
C           SOLUTION IS DONE BY MEANS OF GAUSS-ELIMINATION WITH         DELG 500
C           COMPLETE PIVOTING.                                          DELG 510
C                                                                       DELG 520
C     ..................................................................DELG 530
C                                                                       DELG 540
      SUBROUTINE DGELG(R,A,M,N,EPS,IER)                                 DELG 550
C                                                                       DELG 560
C                                                                       DELG 570
      DIMENSION A(1),R(1)                                               DELG 580
      DOUBLE PRECISION R,A,PIV,TB,TOL,PIVI                              DELG 590
      IF(M)23,23,1                                                      DELG 600
C                                                                       DELG 610
C     SEARCH FOR GREATEST ELEMENT IN MATRIX A                           DELG 620
    1 IER=0                                                             DELG 630
      PIV=0.D0                                                          DELG 640
      MM=M*M                                                            DELG 650
      NM=N*M                                                            DELG 660
      DO 3 L=1,MM                                                       DELG 670
      TB=DABS(A(L))                                                     DELG 680
      IF(TB-PIV)3,3,2                                                   DELG 690
    2 PIV=TB                                                            DELG 700
      I=L                                                               DELG 710
    3 CONTINUE                                                          DELG 720
      TOL=EPS*PIV                                                       DELG 730
C     A(I) IS PIVOT ELEMENT. PIV CONTAINS THE ABSOLUTE VALUE OF A(I).   DELG 740
C                                                                       DELG 750
C                                                                       DELG 760
C     START ELIMINATION LOOP                                            DELG 770
      LST=1                                                             DELG 780
      DO 17 K=1,M                                                       DELG 790
C                                                                       DELG 800
C     TEST ON SINGULARITY                                               DELG 810
      IF(PIV)23,23,4                                                    DELG 820
    4 IF(IER)7,5,7                                                      DELG 830
    5 IF(PIV-TOL)6,6,7                                                  DELG 840
    6 IER=K-1                                                           DELG 850
    7 PIVI=1.D0/A(I)                                                    DELG 860
      J=(I-1)/M                                                         DELG 870
      I=I-J*M-K                                                         DELG 880
      J=J+1-K                                                           DELG 890
C     I+K IS ROW-INDEX, J+K COLUMN-INDEX OF PIVOT ELEMENT               DELG 900
C                                                                       DELG 910
C     PIVOT ROW REDUCTION AND ROW INTERCHANGE IN RIGHT HAND SIDE R      DELG 920
      DO 8 L=K,NM,M                                                     DELG 930
      LL=L+I                                                            DELG 940
      TB=PIVI*R(LL)                                                     DELG 950
      R(LL)=R(L)                                                        DELG 960
    8 R(L)=TB                                                           DELG 970
C                                                                       DELG 980
C     IS ELIMINATION TERMINATED                                         DELG 990
      IF(K-M)9,18,18                                                    DELG1000
C                                                                       DELG1010
C     COLUMN INTERCHANGE IN MATRIX A                                    DELG1020
    9 LEND=LST+M-K                                                      DELG1030
      IF(J)12,12,10                                                     DELG1040
   10 II=J*M                                                            DELG1050
      DO 11 L=LST,LEND                                                  DELG1060
      TB=A(L)                                                           DELG1070
      LL=L+II                                                           DELG1080
      A(L)=A(LL)                                                        DELG1090
   11 A(LL)=TB                                                          DELG1100
C                                                                       DELG1110
C     ROW INTERCHANGE AND PIVOT ROW REDUCTION IN MATRIX A               DELG1120
   12 DO 13 L=LST,MM,M                                                  DELG1130
      LL=L+I                                                            DELG1140
      TB=PIVI*A(LL)                                                     DELG1150
      A(LL)=A(L)                                                        DELG1160
   13 A(L)=TB                                                           DELG1170
C                                                                       DELG1180
C     SAVE COLUMN INTERCHANGE INFORMATION                               DELG1190
      A(LST)=J                                                          DELG1200
C                                                                       DELG1210
C     ELEMENT REDUCTION AND NEXT PIVOT SEARCH                           DELG1220
      PIV=0.D0                                                          DELG1230
      LST=LST+1                                                         DELG1240
      J=0                                                               DELG1250
      DO 16 II=LST,LEND                                                 DELG1260
      PIVI=-A(II)                                                       DELG1270
      IST=II+M                                                          DELG1280
      J=J+1                                                             DELG1290
      DO 15 L=IST,MM,M                                                  DELG1300
      LL=L-J                                                            DELG1310
      A(L)=A(L)+PIVI*A(LL)                                              DELG1320
      TB=DABS(A(L))                                                     DELG1330
      IF(TB-PIV)15,15,14                                                DELG1340
   14 PIV=TB                                                            DELG1350
      I=L                                                               DELG1360
   15 CONTINUE                                                          DELG1370
      DO 16 L=K,NM,M                                                    DELG1380
      LL=L+J                                                            DELG1390
   16 R(LL)=R(LL)+PIVI*R(L)                                             DELG1400
   17 LST=LST+M                                                         DELG1410
C     END OF ELIMINATION LOOP                                           DELG1420
C                                                                       DELG1430
C                                                                       DELG1440
C     BACK SUBSTITUTION AND BACK INTERCHANGE                            DELG1450
   18 IF(M-1)23,22,19                                                   DELG1460
   19 IST=MM+M                                                          DELG1470
      LST=M+1                                                           DELG1480
      DO 21 I=2,M                                                       DELG1490
      II=LST-I                                                          DELG1500
      IST=IST-LST                                                       DELG1510
      L=IST-M                                                           DELG1520
      L=A(L)+.5D0                                                       DELG1530
      DO 21 J=II,NM,M                                                   DELG1540
      TB=R(J)                                                           DELG1550
      LL=J                                                              DELG1560
      DO 20 K=IST,MM,M                                                  DELG1570
      LL=LL+1                                                           DELG1580
   20 TB=TB-A(K)*R(LL)                                                  DELG1590
      K=J+L                                                             DELG1600
      R(J)=R(K)                                                         DELG1610
   21 R(K)=TB                                                           DELG1620
   22 RETURN                                                            DELG1630
C                                                                       DELG1640
C                                                                       DELG1650
C     ERROR RETURN                                                      DELG1660
   23 IER=-1                                                            DELG1670
      RETURN                                                            DELG1680
      END                                                               DELG1690
