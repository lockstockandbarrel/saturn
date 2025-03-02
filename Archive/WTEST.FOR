C                                                                       WTES  10
C     ..................................................................WTES  20
C                                                                       WTES  30
C        SUBROUTINE WTEST                                               WTES  40
C                                                                       WTES  50
C        PURPOSE                                                        WTES  60
C           TEST DEGREE OF ASSOCIATION AMONG A NUMBER OF VARIABLES BY   WTES  70
C           THE KENDALL COEFFICIENT OF CONCORDANCE                      WTES  80
C                                                                       WTES  90
C        USAGE                                                          WTES 100
C           CALL WTEST(A,R,N,M,WA,W,CS,NDF,NR)                          WTES 110
C                                                                       WTES 120
C        DESCRIPTION OF PARAMETERS                                      WTES 130
C           A   - INPUT MATRIX, N BY M, OF ORIGINAL DATA                WTES 140
C           R   - OUTPUT MATRIX, N BY M, OF RANKED DATA.SMALLEST VALUE  WTES 150
C                 IS RANKED 1, LARGEST IS RANKED N. TIES ARE ASSIGNED   WTES 160
C                 AVERAGE OF TIED RANKS                                 WTES 170
C           N   - NUMBER OF VARIABLES                                   WTES 180
C           M   - NUMBER OF CASES                                       WTES 190
C           WA  - WORK AREA VECTOR OF LENGTH 2*M                        WTES 200
C           W   - KENDALL COEFFICIENT OF CONCORDANCE(OUTPUT)            WTES 210
C           CS  - CHI-SQUARE (OUTPUT)                                   WTES 220
C           NDF - NUMBER OF DEGREES OF FREEDOM (OUTPUT)                 WTES 230
C           NR  - CODE, 0 FOR UNRANKED DATA IN A, 1 FOR RANKED DATA     WTES 240
C                 IN A (INPUT)                                          WTES 250
C                                                                       WTES 260
C        REMARKS                                                        WTES 270
C           CHI-SQUARE IS SET TO ZERO IF M IS 7 OR SMALLER              WTES 280
C                                                                       WTES 290
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  WTES 300
C           RANK                                                        WTES 310
C           TIE                                                         WTES 320
C                                                                       WTES 330
C        METHOD                                                         WTES 340
C           DESCRIBED IN S. SIEGEL, 'NONPARAMETRIC STATISTICS FOR THE   WTES 350
C           BEHAVIORAL SCIENCES', MCGRAW-HILL, NEW YORK, 1956,          WTES 360
C           CHAPTER 9                                                   WTES 370
C     ..................................................................WTES 380
C                                                                       WTES 390
C                                                                       WTES 400
      SUBROUTINE WTEST (A,R,N,M,WA,W,CS,NDF,NR)                         WTES 410
      DIMENSION A(1),R(1),WA(1)                                         WTES 420
C                                                                       WTES 430
      FM=M                                                              WTES 440
      FN=N                                                              WTES 450
C                                                                       WTES 460
C        DETERMINE WHETHER DATA IS RANKED                               WTES 470
C        RANK DATA FOR ALL VARIABLES ASSIGNING TIED OBSERVATIONS AVERAGEWTES 480
C        OF TIED RANKS AND COMPUTE CORRECTION FOR TIED SCORES           WTES 490
C                                                                       WTES 500
      T=0.0                                                             WTES 510
      KT=1                                                              WTES 520
      DO 20 I=1,N                                                       WTES 530
      IJ=I-N                                                            WTES 540
      IK=IJ                                                             WTES 550
      IF(NR-1) 5,2,5                                                    WTES 560
    2 DO 3 J=1,M                                                        WTES 570
      IJ=IJ+N                                                           WTES 580
      K=M+J                                                             WTES 590
    3 WA(K)=A(IJ)                                                       WTES 600
      GO TO 15                                                          WTES 610
    5 DO 10 J=1,M                                                       WTES 620
      IJ=IJ+N                                                           WTES 630
   10 WA(J)=A(IJ)                                                       WTES 640
      CALL RANK(WA,WA(M+1),M)                                           WTES 650
   15 CALL TIE(WA(M+1),M,KT,TI)                                         WTES 660
      T=T+TI                                                            WTES 670
      DO 20 J=1,M                                                       WTES 680
      IK=IK+N                                                           WTES 690
      IW=M+J                                                            WTES 700
   20 R(IK)=WA(IW)                                                      WTES 710
C                                                                       WTES 720
C        CALCULATE VECTOR OF SUMS OF RANKS                              WTES 730
C                                                                       WTES 740
      IR=0                                                              WTES 750
      DO 40 J=1,M                                                       WTES 760
      WA(J)=0.0                                                         WTES 770
      DO 40 I=1,N                                                       WTES 780
      IR=IR+1                                                           WTES 790
   40 WA(J)=WA(J)+R(IR)                                                 WTES 800
C                                                                       WTES 810
C        COMPUTE MEAN OF SUMS OF RANKS                                  WTES 820
C                                                                       WTES 830
      SM=0.0                                                            WTES 840
      DO 50 J=1,M                                                       WTES 850
   50 SM=SM+WA(J)                                                       WTES 860
      SM=SM/FM                                                          WTES 870
C                                                                       WTES 880
C        COMPUTE SUM OF SQUARES OF DEVIATIONS                           WTES 890
C                                                                       WTES 900
      S=0.0                                                             WTES 910
      DO 60 J=1,M                                                       WTES 920
   60 S=S+(WA(J)-SM)*(WA(J)-SM)                                         WTES 930
C                                                                       WTES 940
C        COMPUTE W                                                      WTES 950
C                                                                       WTES 960
      W=S/(((FN*FN)*(FM*FM*FM-FM)/12.0)-FN*T)                           WTES 970
C                                                                       WTES 980
C        COMPUTE DEGREES OF FREEDOM AND CHI-SQUARE IF M IS OVER 7       WTES 990
C                                                                       WTES1000
      CS=0.0                                                            WTES1010
      NDF=0                                                             WTES1020
      IF(M-7) 70,70,65                                                  WTES1030
   65 CS=FN*(FM-1.0)*W                                                  WTES1040
      NDF=M-1                                                           WTES1050
   70 RETURN                                                            WTES1060
      END                                                               WTES1070
