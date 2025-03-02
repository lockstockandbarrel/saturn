C                                                                       CROS  10
C     ..................................................................CROS  20
C                                                                       CROS  30
C        SUBROUTINE CROSS                                               CROS  40
C                                                                       CROS  50
C        PURPOSE                                                        CROS  60
C           TO FIND THE CROSSCOVARIANCES OF SERIES A WITH SERIES B      CROS  70
C           (WHICH LEADS AND LAGS A).                                   CROS  80
C                                                                       CROS  90
C        USAGE                                                          CROS 100
C           CALL CROSS (A,B,N,L,R,S)                                    CROS 110
C                                                                       CROS 120
C        DESCRIPTION OF PARAMETERS                                      CROS 130
C           A    - INPUT VECTOR OF LENGTH N CONTAINING FIRST TIME       CROS 140
C                  SERIES.                                              CROS 150
C           B    - INPUT VECTOR OF LENGTH N CONTAINING SECOND TIME      CROS 160
C                  SERIES.                                              CROS 170
C           N    - LENGTH OF SERIES A AND B.                            CROS 180
C           L    - CROSSCOVARIANCE IS CALCULATED FOR LAGS AND LEADS OF  CROS 190
C                  0, 1, 2,..., L-1.                                    CROS 200
C           R    - OUTPUT VECTOR OF LENGTH L CONTAINING CROSSCOVARI-    CROS 210
C                  ANCES OF A WITH B, WHERE B LAGS A.                   CROS 220
C           S    - OUTPUT VECTOR OF LENGTH L CONTAINING CROSSCOVARI-    CROS 230
C                  ANCES OF A WITH B, WHERE B LEADS A.                  CROS 240
C                                                                       CROS 250
C        REMARKS                                                        CROS 260
C           N MUST BE GREATER THAN L.  IF NOT, R(1) AND S(1) ARE SET TO CROS 270
C           ZERO AND RETURN IS MADE TO THE CALLING PROGRAM.             CROS 280
C                                                                       CROS 290
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  CROS 300
C           NONE                                                        CROS 310
C                                                                       CROS 320
C        METHOD                                                         CROS 330
C           DESCRIBED IN R.B. BLACKMAN AND J.W. TUKEY, 'THE MEASURMENT  CROS 340
C        OF POWER SPECTRA', DOVER PUBLICATIONS INC., NEW YORK, 1959.    CROS 350
C                                                                       CROS 360
C     ..................................................................CROS 370
C                                                                       CROS 380
      SUBROUTINE CROSS (A,B,N,L,R,S)                                    CROS 390
      DIMENSION A(1),B(1),R(1),S(1)                                     CROS 400
C                                                                       CROS 410
C     CALCULATE AVERAGES OF SERIES A AND B                              CROS 420
C                                                                       CROS 430
      FN=N                                                              CROS 440
      AVERA=0.0                                                         CROS 450
      AVERB=0.0                                                         CROS 460
      IF(N-L)50,50,100                                                  CROS 470
   50 R(1)=0.0                                                          CROS 480
      S(1)=0.0                                                          CROS 490
      RETURN                                                            CROS 500
  100 DO 110 I=1,N                                                      CROS 510
      AVERA=AVERA+A(I)                                                  CROS 520
  110 AVERB=AVERB+B(I)                                                  CROS 530
      AVERA=AVERA/FN                                                    CROS 540
      AVERB=AVERB/FN                                                    CROS 550
C                                                                       CROS 560
C     CALCULATE CROSSCOVARIANCES OF SERIES A AND B                      CROS 570
C                                                                       CROS 580
      DO 130 J=1,L                                                      CROS 590
      NJ=N-J+1                                                          CROS 600
      SUMR=0.0                                                          CROS 610
      SUMS=0.0                                                          CROS 620
      DO 120 I=1,NJ                                                     CROS 630
      IJ=I+J-1                                                          CROS 640
      SUMR=SUMR+(A(I)-AVERA)*(B(IJ)-AVERB)                              CROS 650
  120 SUMS=SUMS+(A(IJ)-AVERA)*(B(I)-AVERB)                              CROS 660
      FNJ=NJ                                                            CROS 670
      R(J)=SUMR/FNJ                                                     CROS 680
  130 S(J)=SUMS/FNJ                                                     CROS 690
      RETURN                                                            CROS 700
      END                                                               CROS 710
