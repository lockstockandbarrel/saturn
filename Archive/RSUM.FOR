C                                                                       RSUM  10
C     ..................................................................RSUM  20
C                                                                       RSUM  30
C        SUBROUTINE RSUM                                                RSUM  40
C                                                                       RSUM  50
C        PURPOSE                                                        RSUM  60
C           SUM ELEMENTS OF EACH ROW TO FORM COLUMN VECTOR              RSUM  70
C                                                                       RSUM  80
C        USAGE                                                          RSUM  90
C           CALL RSUM (A,R,N,M,MS)                                      RSUM 100
C                                                                       RSUM 110
C        DESCRIPTION OF PARAMETERS                                      RSUM 120
C           A - NAME OF INPUT MATRIX                                    RSUM 130
C           R - NAME OF VECTOR OF LENGTH N                              RSUM 140
C           N - NUMBER OF ROWS IN A                                     RSUM 150
C           M - NUMBER OF COLUMNS IN A                                  RSUM 160
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         RSUM 170
C                  0 - GENERAL                                          RSUM 180
C                  1 - SYMMETRIC                                        RSUM 190
C                  2 - DIAGONAL                                         RSUM 200
C                                                                       RSUM 210
C        REMARKS                                                        RSUM 220
C           VECTOR R CANNOT BE IN THE SAME LOCATION AS MATRIX A         RSUM 230
C           UNLESS A IS GENERAL                                         RSUM 240
C                                                                       RSUM 250
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  RSUM 260
C           LOC                                                         RSUM 270
C                                                                       RSUM 280
C        METHOD                                                         RSUM 290
C           ELEMENTS ARE SUMMED ACROSS EACH ROW INTO A CORRESPONDING    RSUM 300
C           ELEMENT OF OUTPUT COLUMN VECTOR R                           RSUM 310
C                                                                       RSUM 320
C     ..................................................................RSUM 330
C                                                                       RSUM 340
      SUBROUTINE RSUM(A,R,N,M,MS)                                       RSUM 350
      DIMENSION A(1),R(1)                                               RSUM 360
C                                                                       RSUM 370
      DO 3 I=1,N                                                        RSUM 380
C                                                                       RSUM 390
C        CLEAR OUTPUT LOCATION                                          RSUM 400
C                                                                       RSUM 410
      R(I)=0.0                                                          RSUM 420
C                                                                       RSUM 430
      DO 3 J=1,M                                                        RSUM 440
C                                                                       RSUM 450
C        LOCATE ELEMENT FOR ANY MATRIX STORAGE MODE                     RSUM 460
C                                                                       RSUM 470
      CALL LOC(I,J,IJ,N,M,MS)                                           RSUM 480
C                                                                       RSUM 490
C        TEST FOR ZERO ELEMENT IN DIAGONAL MATRIX                       RSUM 500
C                                                                       RSUM 510
      IF(IJ) 2,3,2                                                      RSUM 520
C                                                                       RSUM 530
C        ACCUMULATE IN OUTPUT VECTOR                                    RSUM 540
C                                                                       RSUM 550
    2 R(I)=R(I)+A(IJ)                                                   RSUM 560
    3 CONTINUE                                                          RSUM 570
      RETURN                                                            RSUM 580
      END                                                               RSUM 590
