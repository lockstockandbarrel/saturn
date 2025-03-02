C                                                                       TWOA  10
C     ..................................................................TWOA  20
C                                                                       TWOA  30
C        SUBROUTINE TWOAV                                               TWOA  40
C                                                                       TWOA  50
C        PURPOSE                                                        TWOA  60
C           TEST WHETHER A NUMBER OF SAMPLES ARE FROM THE SAME          TWOA  70
C           POPULATION BY THE FRIEDMAN TWO-WAY ANALYSIS OF VARIANCE TESTTWOA  80
C                                                                       TWOA  90
C        USAGE                                                          TWOA 100
C           CALL TWOAV(A,R,N,M,W,XR,NDF,NR)                             TWOA 110
C                                                                       TWOA 120
C        DESCRIPTION OF PARAMETERS                                      TWOA 130
C           A   - INPUT MATRIX, N BY M, OF ORIGINAL DATA                TWOA 140
C           R   - OUTPUT MATRIX, N BY M, OF RANKED DATA                 TWOA 150
C           N   - NUMBER OF GROUPS                                      TWOA 160
C           M   - NUMBER OF CASES IN EACH GROUP                         TWOA 170
C           W   - WORK AREA OF LENGTH 2*M                               TWOA 180
C           XR  - FRIEDMAN STATISTIC (OUTPUT)                           TWOA 190
C           NDF - NUMBER OF DEGREES OF FREEDOM (OUTPUT)                 TWOA 200
C           NR  - CODE, 0 FOR UNRANKED DATA IN A, 1 FOR RANKED DATA     TWOA 210
C                 IN A (INPUT)                                          TWOA 220
C                                                                       TWOA 230
C        REMARKS                                                        TWOA 240
C           NONE                                                        TWOA 250
C                                                                       TWOA 260
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  TWOA 270
C           RANK                                                        TWOA 280
C                                                                       TWOA 290
C        METHOD                                                         TWOA 300
C           DESCRIBED IN S. SIEGEL, 'NONPARAMETRIC STATISTICS FOR THE   TWOA 310
C           BEHAVIORAL SCIENCES', MCGRAW-HILL, NEW YORK, 1956,          TWOA 320
C           CHAPTER 7                                                   TWOA 330
C                                                                       TWOA 340
C     ..................................................................TWOA 350
C                                                                       TWOA 360
      SUBROUTINE TWOAV (A,R,N,M,W,XR,NDF,NR)                            TWOA 370
      DIMENSION A(1),R(1),W(1)                                          TWOA 380
C                                                                       TWOA 390
C        DETERMINE WHETHER DATA IS RANKED                               TWOA 400
C                                                                       TWOA 410
      IF(NR-1) 10, 30, 10                                               TWOA 420
C                                                                       TWOA 430
C        RANK DATA IN EACH GROUP AND ASSIGN TIED OBSERVATIONS AVERAGE   TWOA 440
C        OF TIED RANK                                                   TWOA 450
C                                                                       TWOA 460
   10 DO 20 I=1,N                                                       TWOA 470
      IJ=I-N                                                            TWOA 480
      IK=IJ                                                             TWOA 490
      DO 15 J=1,M                                                       TWOA 500
      IJ=IJ+N                                                           TWOA 510
   15 W(J)=A(IJ)                                                        TWOA 520
      CALL RANK (W,W(M+1),M)                                            TWOA 530
      DO 20 J=1,M                                                       TWOA 540
      IK=IK+N                                                           TWOA 550
      IW=M+J                                                            TWOA 560
   20 R(IK)=W(IW)                                                       TWOA 570
      GO TO 35                                                          TWOA 580
   30 NM=N*M                                                            TWOA 590
      DO 32 I=1,NM                                                      TWOA 600
   32 R(I)=A(I)                                                         TWOA 610
C                                                                       TWOA 620
C        CALCULATE SUM OF SQUARES OF SUMS OF RANKS                      TWOA 630
C                                                                       TWOA 640
   35 RTSQ=0.0                                                          TWOA 650
      IR=0                                                              TWOA 660
      DO 50 J=1,M                                                       TWOA 670
      RT=0.0                                                            TWOA 680
      DO 40 I=1,N                                                       TWOA 690
      IR=IR+1                                                           TWOA 700
   40 RT=RT+R(IR)                                                       TWOA 710
   50 RTSQ=RTSQ+RT*RT                                                   TWOA 720
C                                                                       TWOA 730
C        CALCULATE FRIEDMAN TEST VALUE, XR                              TWOA 740
C                                                                       TWOA 750
      FNM=N*(M+1)                                                       TWOA 760
      FM=M                                                              TWOA 770
      XR=(12.0/(FM*FNM))*RTSQ-3.0*FNM                                   TWOA 780
C                                                                       TWOA 790
C        FIND DEGREES OF FREEDOM                                        TWOA 800
C                                                                       TWOA 810
      NDF=M-1                                                           TWOA 820
      RETURN                                                            TWOA 830
      END                                                               TWOA 840
