C                                                                       TALL  10
C     ..................................................................TALL  20
C                                                                       TALL  30
C        SUBROUTINE TALLY                                               TALL  40
C                                                                       TALL  50
C        PURPOSE                                                        TALL  60
C           CALCULATE TOTAL, MEAN, STANDARD DEVIATION, MINIMUM, MAXIMUM TALL  70
C           FOR EACH VARIABLE IN A SET (OR A SUBSET) OF OBSERVATIONS    TALL  80
C                                                                       TALL  90
C        USAGE                                                          TALL 100
C           CALL TALLY(A,S,TOTAL,AVER,SD,VMIN,VMAX,NO,NV,IER)           TALL 110
C                                                                       TALL 120
C        DESCRIPTION OF PARAMETERS                                      TALL 130
C           A     - OBSERVATION MATRIX, NO BY NV                        TALL 140
C           S     - INPUT VECTOR INDICATING SUBSET OF A. ONLY THOSE     TALL 150
C                   OBSERVATIONS WITH A NON-ZERO S(J) ARE CONSIDERED.   TALL 160
C                   VECTOR LENGTH IS NO.                                TALL 170
C           TOTAL - OUTPUT VECTOR OF TOTALS OF EACH VARIABLE. VECTOR    TALL 180
C                   LENGTH IS NV.                                       TALL 190
C           AVER  - OUTPUT VECTOR OF AVERAGES OF EACH VARIABLE. VECTOR  TALL 200
C                   LENGTH IS NV.                                       TALL 210
C           SD    - OUTPUT VECTOR OF STANDARD DEVIATIONS OF EACH        TALL 220
C                   VARIABLE. VECTOR LENGTH IS NV.                      TALL 230
C           VMIN  - OUTPUT VECTOR OF MINIMA OF EACH VARIABLE. VECTOR    TALL 240
C                   LENGTH IS NV.                                       TALL 250
C           VMAX  - OUTPUT VECTOR OF MAXIMA OF EACH VARIABLE. VECTOR    TALL 260
C                   LENGTH IS NV.                                       TALL 270
C           NO    - NUMBER OF OBSERVATIONS                              TALL 280
C           NV    - NUMBER OF VARIABLES FOR EACH OBSERVATION            TALL 290
C           IER   - ZERO, IF NO ERROR.                                  TALL 291
C                 - 1, IF S IS NULL.  VMIN=-1.E75, VMAX=SD=AVER=1.E75.  TALL 292
C                 - 2, IF S HAS ONLY ONE NON-ZERO ELEMENT. VMIN=VMAX.   TALL 293
C                   SD=0.0                                              TALL 294
C                                                                       TALL 300
C        REMARKS                                                        TALL 310
C           NONE                                                        TALL 320
C                                                                       TALL 330
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  TALL 340
C           NONE                                                        TALL 350
C                                                                       TALL 360
C        METHOD                                                         TALL 370
C           ALL OBSERVATIONS CORRESPONDING TO A NON-ZERO ELEMENT IN S   TALL 380
C           VECTOR ARE ANALYZED FOR EACH VARIABLE IN MATRIX A.          TALL 390
C           TOTALS ARE ACCUMULATED AND MINIMUM AND MAXIMUM VALUES ARE   TALL 400
C           FOUND. FOLLOWING THIS, MEANS AND STANDARD DEVIATIONS ARE    TALL 410
C           CALCULATED.  THE DIVISOR FOR STANDARD DEVIATION IS ONE LESS TALL 420
C           THAN THE NUMBER OF OBSERVATIONS USED.                       TALL 430
C                                                                       TALL 440
C     ..................................................................TALL 450
C                                                                       TALL 460
      SUBROUTINE TALLY(A,S,TOTAL,AVER,SD,VMIN,VMAX,NO,NV)               TALL 470
      DIMENSION A(1),S(1),TOTAL(1),AVER(1),SD(1),VMIN(1),VMAX(1)        TALL 480
C                                                                       TALL 490
C        CLEAR OUTPUT VECTORS AND INITIALIZE VMIN,VMAX                  TALL 500
C                                                                       TALL 510
      IER=0                                                             TALL 511
      DO 1 K=1,NV                                                       TALL 520
      TOTAL(K)=0.0                                                      TALL 530
      AVER(K)=1.0E75                                                    TALL 540
      SD(K)=1.0E75                                                      TALL 550
      VMIN(K)=-1.0E75                                                   TALL 560
    1 VMAX(K)=1.0E75                                                    TALL 570
C                                                                       TALL 580
C        TEST SUBSET VECTOR                                             TALL 590
C                                                                       TALL 600
      SCNT=0.0                                                          TALL 610
      DO 7 J=1,NO                                                       TALL 620
      IJ=J-NO                                                           TALL 630
      IF(S(J)) 2,7,2                                                    TALL 640
    2 SCNT=SCNT+1.0                                                     TALL 650
C                                                                       TALL 660
C        CALCULATE TOTAL, MINIMA, MAXIMA                                TALL 670
C                                                                       TALL 680
      DO 6 I=1,NV                                                       TALL 690
      IJ=IJ+NO                                                          TALL 700
      TOTAL(I)=TOTAL(I)+A(IJ)                                           TALL 710
      IF(A(IJ)-VMIN(I)) 3,4,4                                           TALL 720
    3 VMIN(I)=A(IJ)                                                     TALL 730
    4 IF(A(IJ)-VMAX(I)) 6,6,5                                           TALL 740
    5 VMAX(I)=A(IJ)                                                     TALL 750
    6 SD(I)=SD(I)+A(IJ)*A(IJ)                                           TALL 760
    7 CONTINUE                                                          TALL 770
C                                                                       TALL 780
C        CALCULATE MEANS AND STANDARD DEVIATIONS                        TALL 790
C                                                                       TALL 800
      IF (SCNT)8,8,9                                                    TALL 801
    8 IER=1                                                             TALL 802
      GO TO 15                                                          TALL 803
    9 DO 10 I=1,NV                                                      TALL 810
   10 AVER(I)=TOTAL(I)/SCNT                                             TALL 820
      IF (SCNT-1.0) 13,11,13                                            TALL 821
   11 IER=2                                                             TALL 822
      DO 12 I=1,NV                                                      TALL 823
   12 SD(I)=0.0                                                         TALL 824
      GO TO 15                                                          TALL 825
   13 DO 14 I=1,NV                                                      TALL 826
   14 SD(I)=SQRT(ABS((SD(I)-TOTAL(I)*TOTAL(I)/SCNT)/(SCNT-1.0)))        TALL 830
   15 RETURN                                                            TALL 840
      END                                                               TALL 850
