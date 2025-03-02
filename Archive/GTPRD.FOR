C                                                                       GTPR  10
C     ..................................................................GTPR  20
C                                                                       GTPR  30
C        SUBROUTINE GTPRD                                               GTPR  40
C                                                                       GTPR  50
C        PURPOSE                                                        GTPR  60
C           PREMULTIPLY A GENERAL MATRIX BY THE TRANSPOSE OF ANOTHER    GTPR  70
C           GENERAL MATRIX                                              GTPR  80
C                                                                       GTPR  90
C        USAGE                                                          GTPR 100
C           CALL GTPRD(A,B,R,N,M,L)                                     GTPR 110
C                                                                       GTPR 120
C        DESCRIPTION OF PARAMETERS                                      GTPR 130
C           A - NAME OF FIRST INPUT MATRIX                              GTPR 140
C           B - NAME OF SECOND INPUT MATRIX                             GTPR 150
C           R - NAME OF OUTPUT MATRIX                                   GTPR 160
C           N - NUMBER OF ROWS IN A AND B                               GTPR 170
C           M - NUMBER OF COLUMNS IN A AND ROWS IN R                    GTPR 180
C           L - NUMBER OF COLUMNS IN B AND R                            GTPR 190
C                                                                       GTPR 200
C        REMARKS                                                        GTPR 210
C           MATRIX R CANNOT BE IN THE SAME LOCATION AS MATRIX A         GTPR 220
C           MATRIX R CANNOT BE IN THE SAME LOCATION AS MATRIX B         GTPR 230
C           ALL MATRICES MUST BE STORED AS GENERAL MATRICES             GTPR 240
C                                                                       GTPR 250
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  GTPR 260
C           NONE                                                        GTPR 270
C                                                                       GTPR 280
C        METHOD                                                         GTPR 290
C           MATRIX TRANSPOSE OF A IS NOT ACTUALLY CALCULATED. INSTEAD,  GTPR 300
C           ELEMENTS OF MATRIX A ARE TAKEN COLUMNWISE RATHER THAN       GTPR 310
C           ROWWISE FOR POSTMULTIPLICATION BY MATRIX B.                 GTPR 320
C                                                                       GTPR 330
C     ..................................................................GTPR 340
C                                                                       GTPR 350
      SUBROUTINE GTPRD(A,B,R,N,M,L)                                     GTPR 360
      DIMENSION A(1),B(1),R(1)                                          GTPR 370
C                                                                       GTPR 380
      IR=0                                                              GTPR 390
      IK=-N                                                             GTPR 400
      DO 10 K=1,L                                                       GTPR 410
      IJ=0                                                              GTPR 420
      IK=IK+N                                                           GTPR 430
      DO 10 J=1,M                                                       GTPR 440
      IB=IK                                                             GTPR 450
      IR=IR+1                                                           GTPR 460
      R(IR)=0                                                           GTPR 470
      DO 10 I=1,N                                                       GTPR 480
      IJ=IJ+1                                                           GTPR 490
      IB=IB+1                                                           GTPR 500
   10 R(IR)=R(IR)+A(IJ)*B(IB)                                           GTPR 510
      RETURN                                                            GTPR 520
      END                                                               GTPR 530
