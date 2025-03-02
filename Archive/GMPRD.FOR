C                                                                       GMPR  10
C     ..................................................................GMPR  20
C                                                                       GMPR  30
C        SUBROUTINE GMPRD                                               GMPR  40
C                                                                       GMPR  50
C        PURPOSE                                                        GMPR  60
C           MULTIPLY TWO GENERAL MATRICES TO FORM A RESULTANT GENERAL   GMPR  70
C           MATRIX                                                      GMPR  80
C                                                                       GMPR  90
C        USAGE                                                          GMPR 100
C           CALL GMPRD(A,B,R,N,M,L)                                     GMPR 110
C                                                                       GMPR 120
C        DESCRIPTION OF PARAMETERS                                      GMPR 130
C           A - NAME OF FIRST INPUT MATRIX                              GMPR 140
C           B - NAME OF SECOND INPUT MATRIX                             GMPR 150
C           R - NAME OF OUTPUT MATRIX                                   GMPR 160
C           N - NUMBER OF ROWS IN A                                     GMPR 170
C           M - NUMBER OF COLUMNS IN A AND ROWS IN B                    GMPR 180
C           L - NUMBER OF COLUMNS IN B                                  GMPR 190
C                                                                       GMPR 200
C        REMARKS                                                        GMPR 210
C           ALL MATRICES MUST BE STORED AS GENERAL MATRICES             GMPR 220
C           MATRIX R CANNOT BE IN THE SAME LOCATION AS MATRIX A         GMPR 230
C           MATRIX R CANNOT BE IN THE SAME LOCATION AS MATRIX B         GMPR 240
C           NUMBER OF COLUMNS OF MATRIX A MUST BE EQUAL TO NUMBER OF ROWGMPR 250
C           OF MATRIX B                                                 GMPR 260
C                                                                       GMPR 270
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  GMPR 280
C           NONE                                                        GMPR 290
C                                                                       GMPR 300
C        METHOD                                                         GMPR 310
C           THE M BY L MATRIX B IS PREMULTIPLIED BY THE N BY M MATRIX A GMPR 320
C           AND THE RESULT IS STORED IN THE N BY L MATRIX R.            GMPR 330
C                                                                       GMPR 340
C     ..................................................................GMPR 350
C                                                                       GMPR 360
      SUBROUTINE GMPRD(A,B,R,N,M,L)                                     GMPR 370
      DIMENSION A(1),B(1),R(1)                                          GMPR 380
C                                                                       GMPR 390
      IR=0                                                              GMPR 400
      IK=-M                                                             GMPR 410
      DO 10 K=1,L                                                       GMPR 420
      IK=IK+M                                                           GMPR 430
      DO 10 J=1,N                                                       GMPR 440
      IR=IR+1                                                           GMPR 450
      JI=J-N                                                            GMPR 460
      IB=IK                                                             GMPR 470
      R(IR)=0                                                           GMPR 480
      DO 10 I=1,M                                                       GMPR 490
      JI=JI+N                                                           GMPR 500
      IB=IB+1                                                           GMPR 510
   10 R(IR)=R(IR)+A(JI)*B(IB)                                           GMPR 520
      RETURN                                                            GMPR 530
      END                                                               GMPR 540
