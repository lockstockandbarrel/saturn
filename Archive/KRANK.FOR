C                                                                       KRAN  10
C     ..................................................................KRAN  20
C                                                                       KRAN  30
C        SUBROUTINE KRANK                                               KRAN  40
C                                                                       KRAN  50
C        PURPOSE                                                        KRAN  60
C           TEST CORRELATION BETWEEN TWO VARIABLES BY MEANS OF KENDALL  KRAN  70
C           RANK CORRELATION COEFFICIENT                                KRAN  80
C                                                                       KRAN  90
C        USAGE                                                          KRAN 100
C           CALL KRANK(A,B,R,N,TAU,SD,Z,NR)                             KRAN 110
C                                                                       KRAN 120
C        DESCRIPTION OF PARAMETERS                                      KRAN 130
C           A   - INPUT VECTOR OF N OBSERVATIONS FOR FIRST VARIABLE     KRAN 140
C           B   - INPUT VECTOR OF N OBSERVATIONS FOR SECOND VARIABLE    KRAN 150
C           R   - OUTPUT VECTOR OF RANKED DATA OF LENGTH 2*N. SMALLEST  KRAN 160
C                 OBSERVATION IS RANKED 1, LARGEST IS RANKED N. TIES    KRAN 170
C                 ARE ASSIGNED AVERAGE OF TIED RANKS.                   KRAN 180
C           N   - NUMBER OF OBSERVATIONS                                KRAN 190
C           TAU - KENDALL RANK CORRELATION COEFFICIENT (OUTPUT)         KRAN 200
C           SD  - STANDARD DEVIATION (OUTPUT)                           KRAN 210
C           Z   - TEST OF SIGNIFICANCE OF TAU IN TERMS OF NORMAL        KRAN 220
C                 DISTRIBUTION (OUTPUT)                                 KRAN 230
C           NR  - CODE, 0 FOR UNRANKED DATA IN A AND B, 1 FOR RANKED    KRAN 240
C                 DATA IN A AND B (INPUT)                               KRAN 250
C                                                                       KRAN 260
C        REMARKS                                                        KRAN 270
C           SD AND Z ARE SET TO ZERO IF N IS LESS THAN TEN              KRAN 280
C                                                                       KRAN 290
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  KRAN 300
C           RANK                                                        KRAN 310
C           TIE                                                         KRAN 320
C                                                                       KRAN 330
C        METHOD                                                         KRAN 340
C           DESCRIBED IN S. SIEGEL, 'NONPARAMETRIC STATISTICS FOR THE   KRAN 350
C           BEHAVIORAL SCIENCES', MCGRAW-HILL, NEW YORK, 1956,          KRAN 360
C           CHAPTER 9                                                   KRAN 370
C                                                                       KRAN 380
C     ..................................................................KRAN 390
C                                                                       KRAN 400
      SUBROUTINE KRANK(A,B,R,N,TAU,SD,Z,NR)                             KRAN 410
      DIMENSION A(1),B(1),R(1)                                          KRAN 420
C                                                                       KRAN 430
      SD=0.0                                                            KRAN 440
      Z=0.0                                                             KRAN 450
      FN=N                                                              KRAN 460
      FN1=N*(N-1)                                                       KRAN 470
C                                                                       KRAN 480
C        DETERMINE WHETHER DATA IS RANKED                               KRAN 490
C                                                                       KRAN 500
      IF(NR-1) 5, 10, 5                                                 KRAN 510
C                                                                       KRAN 520
C        RANK DATA IN A AND B VECTORS AND ASSIGN TIED OBSERVATIONS      KRAN 530
C        AVERAGE OF TIED RANKS                                          KRAN 540
C                                                                       KRAN 550
    5 CALL RANK (A,R,N)                                                 KRAN 560
      CALL RANK (B,R(N+1),N)                                            KRAN 570
      GO TO 40                                                          KRAN 580
C                                                                       KRAN 590
C        MOVE RANKED DATA TO R VECTOR                                   KRAN 600
C                                                                       KRAN 610
   10 DO 20 I=1,N                                                       KRAN 620
   20 R(I)=A(I)                                                         KRAN 630
      DO 30 I=1,N                                                       KRAN 640
      J=I+N                                                             KRAN 650
   30 R(J)=B(I)                                                         KRAN 660
C                                                                       KRAN 670
C        SORT RANK VECTOR R IN SEQUENCE OF VARIABLE A                   KRAN 680
C                                                                       KRAN 690
   40 ISORT=0                                                           KRAN 700
      DO 50 I=2,N                                                       KRAN 710
      IF(R(I)-R(I-1)) 45,50,50                                          KRAN 720
   45 ISORT=ISORT+1                                                     KRAN 730
      RSAVE=R(I)                                                        KRAN 740
      R(I)=R(I-1)                                                       KRAN 750
      R(I-1)=RSAVE                                                      KRAN 760
      I2=I+N                                                            KRAN 770
      SAVER=R(I2)                                                       KRAN 780
      R(I2)=R(I2-1)                                                     KRAN 790
      R(I2-1)=SAVER                                                     KRAN 800
   50 CONTINUE                                                          KRAN 810
      IF(ISORT) 40,55,40                                                KRAN 820
C                                                                       KRAN 830
C        COMPUTE S ON VARIABLE B. STARTING WITH THE FIRST RANK, ADD 1   KRAN 840
C        TO S FOR EACH LARGER RANK TO ITS RIGHT AND SUBTRACT 1 FOR EACH KRAN 850
C        SMALLER RANK.  REPEAT FOR ALL RANKS.                           KRAN 860
C                                                                       KRAN 870
   55 S=0.0                                                             KRAN 880
      NM=N-1                                                            KRAN 890
      DO 60 I=1,NM                                                      KRAN 900
      J=N+I                                                             KRAN 910
      DO 60 L=I,N                                                       KRAN 920
      K=N+L                                                             KRAN 930
      IF(R(K)-R(J)) 56,60,57                                            KRAN 940
   56 S=S-1.0                                                           KRAN 950
      GO TO 60                                                          KRAN 960
   57 S=S+1.0                                                           KRAN 970
   60 CONTINUE                                                          KRAN 980
C                                                                       KRAN 990
C        COMPUTE TIED SCORE INDEX FOR BOTH VARIABLES                    KRAN1000
C                                                                       KRAN1010
      KT=2                                                              KRAN1020
      CALL TIE(R,N,KT,TA)                                               KRAN1030
      CALL TIE(R(N+1),N,KT,TB)                                          KRAN1040
C                                                                       KRAN1050
C        COMPUTE TAU                                                    KRAN1060
C                                                                       KRAN1070
      IF(TA) 70,65,70                                                   KRAN1080
   65 IF(TB) 70,67,70                                                   KRAN1090
   67 TAU=S/(0.5*FN1)                                                   KRAN1100
      GO TO 80                                                          KRAN1110
   70 TAU=S/((SQRT(0.5*FN1-TA))*(SQRT(0.5*FN1-TB)))                     KRAN1120
C                                                                       KRAN1130
C     COMPUTE STANDARD DEVIATION AND Z IF N IS 10 OR LARGER             KRAN1140
C                                                                       KRAN1150
   80 IF(N-10) 90,85,85                                                 KRAN1160
   85 SD=(SQRT((2.0*(FN+FN+5.0))/(9.0*FN1)))                            KRAN1170
      Z=TAU/SD                                                          KRAN1180
   90 RETURN                                                            KRAN1190
      END                                                               KRAN1200
