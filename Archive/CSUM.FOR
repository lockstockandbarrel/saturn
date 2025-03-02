C                                                                       CSUM  10
C     ..................................................................CSUM  20
C                                                                       CSUM  30
C        SUBROUTINE CSUM                                                CSUM  40
C                                                                       CSUM  50
C        PURPOSE                                                        CSUM  60
C           SUM ELEMENTS OF EACH COLUMN TO FORM ROW VECTOR              CSUM  70
C                                                                       CSUM  80
C        USAGE                                                          CSUM  90
C           CALL CSUM(A,R,N,M,MS)                                       CSUM 100
C                                                                       CSUM 110
C        DESCRIPTION OF PARAMETERS                                      CSUM 120
C           A - NAME OF INPUT MATRIX                                    CSUM 130
C           R - NAME OF VECTOR OF LENGTH M                              CSUM 140
C           N - NUMBER OF ROWS IN A                                     CSUM 150
C           M - NUMBER OF COLUMNS IN A                                  CSUM 160
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         CSUM 170
C                  0 - GENERAL                                          CSUM 180
C                  1 - SYMMETRIC                                        CSUM 190
C                  2 - DIAGONAL                                         CSUM 200
C                                                                       CSUM 210
C        REMARKS                                                        CSUM 220
C           VECTOR R CANNOT BE IN THE SAME LOCATION AS MATRIX A         CSUM 230
C           UNLESS A IS GENERAL                                         CSUM 240
C                                                                       CSUM 250
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  CSUM 260
C           LOC                                                         CSUM 270
C                                                                       CSUM 280
C        METHOD                                                         CSUM 290
C           ELEMENTS ARE SUMMED DOWN EACH COLUMN INTO A CORRESPONDING   CSUM 300
C           ELEMENT OF OUTPUT ROW VECTOR R                              CSUM 310
C                                                                       CSUM 320
C     ..................................................................CSUM 330
C                                                                       CSUM 340
      SUBROUTINE CSUM(A,R,N,M,MS)                                       CSUM 350
      DIMENSION A(1),R(1)                                               CSUM 360
C                                                                       CSUM 370
      DO 3 J=1,M                                                        CSUM 380
C                                                                       CSUM 390
C        CLEAR OUTPUT LOCATION                                          CSUM 400
C                                                                       CSUM 410
      R(J)=0.0                                                          CSUM 420
C                                                                       CSUM 430
      DO 3 I=1,N                                                        CSUM 440
C                                                                       CSUM 450
C        LOCATE ELEMENT FOR ANY MATRIX STORAGE MODE                     CSUM 460
C                                                                       CSUM 470
      CALL LOC(I,J,IJ,N,M,MS)                                           CSUM 480
C                                                                       CSUM 490
C        TEST FOR ZERO ELEMENT IN DIAGONAL MATRIX                       CSUM 500
C                                                                       CSUM 510
      IF(IJ) 2,3,2                                                      CSUM 520
C                                                                       CSUM 530
C        ACCUMULATE IN OUTPUT VECTOR                                    CSUM 540
C                                                                       CSUM 550
    2 R(J)=R(J)+A(IJ)                                                   CSUM 560
    3 CONTINUE                                                          CSUM 570
      RETURN                                                            CSUM 580
      END                                                               CSUM 590
