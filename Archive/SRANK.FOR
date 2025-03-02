C                                                                       SRAN  10
C     ..................................................................SRAN  20
C                                                                       SRAN  30
C        SUBROUTINE SRANK                                               SRAN  40
C                                                                       SRAN  50
C        PURPOSE                                                        SRAN  60
C           TEST CORRELATION BETWEEN TWO VARIABLES BY MEANS OF SPEARMAN SRAN  70
C           RANK CORRELATION COEFFICIENT                                SRAN  80
C                                                                       SRAN  90
C        USAGE                                                          SRAN 100
C           CALL SRANK(A,B,R,N,RS,T,NDF,NR)                             SRAN 110
C                                                                       SRAN 120
C        DESCRIPTION OF PARAMETERS                                      SRAN 130
C           A   - INPUT VECTOR OF N OBSERVATIONS FOR FIRST VARIABLE     SRAN 140
C           B   - INPUT VECTOR OF N OBSERVATIONS FOR SECOND VARIABLE    SRAN 150
C           R   - OUTPUT VECTOR FOR RANKED DATA, LENGTH IS 2*N. SMALLESTSRAN 160
C                 OBSERVATION IS RANKED 1, LARGEST IS RANKED N. TIES    SRAN 170
C                 ARE ASSIGNED AVERAGE OF TIED RANKS.                   SRAN 180
C           N   - NUMBER OF OBSERVATIONS                                SRAN 190
C           RS  - SPEARMAN RANK CORRELATION COEFFICIENT (OUTPUT)        SRAN 200
C           T   - TEST OF SIGNIFICANCE OF RS (OUTPUT)                   SRAN 210
C           NDF - NUMBER OF DEGREES OF FREEDOM (OUTPUT)                 SRAN 220
C           NR  - CODE, 0 FOR UNRANKED DATA IN A AND B, 1 FOR RANKED    SRAN 230
C                 DATA IN A AND B (INPUT)                               SRAN 240
C                                                                       SRAN 250
C        REMARKS                                                        SRAN 260
C           T IS SET TO ZERO IF N IS LESS THAN TEN                      SRAN 270
C                                                                       SRAN 280
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  SRAN 290
C           RANK                                                        SRAN 300
C           TIE                                                         SRAN 310
C                                                                       SRAN 320
C        METHOD                                                         SRAN 330
C           DESCRIBED IN S. SIEGEL, 'NONPARAMETRIC STATISTICS FOR THE   SRAN 340
C           BEHAVIORAL SCIENCES', MCGRAW-HILL, NEW YORK, 1956,          SRAN 350
C           CHAPTER 9                                                   SRAN 360
C                                                                       SRAN 370
C     ..................................................................SRAN 380
C                                                                       SRAN 390
      SUBROUTINE SRANK(A,B,R,N,RS,T,NDF,NR)                             SRAN 400
      DIMENSION A(1),B(1),R(1)                                          SRAN 410
C                                                                       SRAN 420
      FNNN=N*N*N-N                                                      SRAN 430
C                                                                       SRAN 440
C        DETERMINE WHETHER DATA IS RANKED                               SRAN 450
C                                                                       SRAN 460
      IF(NR-1) 5, 10, 5                                                 SRAN 470
C                                                                       SRAN 480
C        RANK DATA IN A AND B VECTORS AND ASSIGN TIED OBSERVATIONS      SRAN 490
C        AVERAGE OF TIED RANKS                                          SRAN 500
C                                                                       SRAN 510
    5 CALL RANK (A,R,N)                                                 SRAN 520
      CALL RANK (B,R(N+1),N)                                            SRAN 530
      GO TO 40                                                          SRAN 540
C                                                                       SRAN 550
C        MOVE RANKED DATA TO R VECTOR                                   SRAN 560
C                                                                       SRAN 570
   10 DO 20 I=1,N                                                       SRAN 580
   20 R(I)=A(I)                                                         SRAN 590
      DO 30 I=1,N                                                       SRAN 600
      J=I+N                                                             SRAN 610
   30 R(J)=B(I)                                                         SRAN 620
C                                                                       SRAN 630
C        COMPUTE SUM OF SQUARES OF RANK DIFFERENCES                     SRAN 640
C                                                                       SRAN 650
   40 D=0.0                                                             SRAN 660
      DO 50 I=1,N                                                       SRAN 670
      J=I+N                                                             SRAN 680
   50 D=D+(R(I)-R(J))*(R(I)-R(J))                                       SRAN 690
C                                                                       SRAN 700
C        COMPUTE TIED SCORE INDEX                                       SRAN 710
C                                                                       SRAN 720
      KT=1                                                              SRAN 730
      CALL TIE (R,N,KT,TSA)                                             SRAN 740
      CALL TIE (R(N+1),N,KT,TSB)                                        SRAN 750
C                                                                       SRAN 760
C        COMPUTE SPEARMAN RANK CORRELATION COEFFICIENT                  SRAN 770
C                                                                       SRAN 780
      IF(TSA) 60,55,60                                                  SRAN 790
   55 IF(TSB) 60,57,60                                                  SRAN 800
   57 RS=1.0-6.0*D/FNNN                                                 SRAN 810
      GO TO 70                                                          SRAN 820
   60 X=FNNN/12.0-TSA                                                   SRAN 830
      Y=X+TSA-TSB                                                       SRAN 840
      RS=(X+Y-D)/(2.0*(SQRT(X*Y)))                                      SRAN 850
C                                                                       SRAN 860
C        COMPUTE T AND DEGREES OF FREEDOM IF N IS 10 OR LARGER          SRAN 870
C                                                                       SRAN 880
      T=0.0                                                             SRAN 890
   70 IF(N-10) 80,75,75                                                 SRAN 900
   75 T=RS*SQRT(FLOAT(N-2)/(1.0-RS*RS))                                 SRAN 910
   80 NDF=N-2                                                           SRAN 920
      RETURN                                                            SRAN 930
      END                                                               SRAN 940
