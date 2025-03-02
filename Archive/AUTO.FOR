C                                                                       AUTO  10
C     ..................................................................AUTO  20
C                                                                       AUTO  30
C        SUBROUTINE AUTO                                                AUTO  40
C                                                                       AUTO  50
C        PURPOSE                                                        AUTO  60
C           TO FIND AUTOCOVARIANCES OF SERIES A FOR LAGS 0 TO L-1.      AUTO  70
C                                                                       AUTO  80
C        USAGE                                                          AUTO  90
C           CALL AUTO (A,N,L,R)                                         AUTO 100
C                                                                       AUTO 110
C        DESCRIPTION OF PARAMETERS                                      AUTO 120
C           A    - INPUT VECTOR OF LENGTH N CONTAINING THE TIME SERIES  AUTO 130
C                  WHOSE AUTOCOVARIANCE IS DESIRED.                     AUTO 140
C           N    - LENGTH OF THE VECTOR A.                              AUTO 150
C           L    - AUTOCOVARIANCE IS CALCULATED FOR LAGS OF 0, 1, 2,...,AUTO 160
C                  L-1.                                                 AUTO 170
C           R    - OUTPUT VECTOR OF LENGTH L CONTAINING AUTOCOVARIANCES AUTO 180
C                  OF SERIES A.                                         AUTO 190
C                                                                       AUTO 200
C        REMARKS                                                        AUTO 210
C           THE LENGTH OF R IS DIFFERENT FROM THE LENGTH OF A.  N MUST  AUTO 220
C           BE GREATER THAN L.  IF NOT, R(1) IS SET TO ZERO AND RETURN  AUTO 230
C           IS MADE TO THE CALLING PROGRAM.                             AUTO 240
C                                                                       AUTO 250
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  AUTO 260
C           NONE                                                        AUTO 270
C                                                                       AUTO 280
C        METHOD                                                         AUTO 290
C           DESCRIBED IN R.B. BLACKMAN AND J.W. TUKEY, 'THE MEASURMENT  AUTO 300
C        OF POWER SPECTRA', DOVER PUBLICATIONS INC., NEW YORK, 1959.    AUTO 310
C                                                                       AUTO 320
C     ..................................................................AUTO 330
C                                                                       AUTO 340
      SUBROUTINE AUTO (A,N,L,R)                                         AUTO 350
      DIMENSION A(1),R(1)                                               AUTO 360
C                                                                       AUTO 370
C     CALCULATE AVERAGE OF TIME SERIES A                                AUTO 380
C                                                                       AUTO 390
      AVER=0.0                                                          AUTO 400
      IF(N-L) 50,50,100                                                 AUTO 410
   50 R(1)=0.0                                                          AUTO 420
      RETURN                                                            AUTO 430
  100 DO 110 I=1,N                                                      AUTO 440
  110 AVER=AVER+A(I)                                                    AUTO 450
      FN=N                                                              AUTO 460
      AVER=AVER/FN                                                      AUTO 470
C                                                                       AUTO 480
C     CALCULATE AUTOCOVARIANCES                                         AUTO 490
C                                                                       AUTO 500
      DO 130 J=1,L                                                      AUTO 510
      NJ=N-J+1                                                          AUTO 520
      SUM=0.0                                                           AUTO 530
      DO 120 I=1,NJ                                                     AUTO 540
      IJ=I+J-1                                                          AUTO 550
  120 SUM=SUM+(A(I)-AVER)*(A(IJ)-AVER)                                  AUTO 560
      FNJ=NJ                                                            AUTO 570
  130 R(J)=SUM/FNJ                                                      AUTO 580
      RETURN                                                            AUTO 590
      END                                                               AUTO 600
