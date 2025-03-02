C                                                                       SADD  10
C     ..................................................................SADD  20
C                                                                       SADD  30
C        SUBROUTINE SADD                                                SADD  40
C                                                                       SADD  50
C        PURPOSE                                                        SADD  60
C           ADD A SCALAR TO EACH ELEMENT OF A MATRIX TO FORM A RESULTANTSADD  70
C           MATRIX                                                      SADD  80
C                                                                       SADD  90
C        USAGE                                                          SADD 100
C           CALL SADD(A,C,R,N,M,MS)                                     SADD 110
C                                                                       SADD 120
C        DESCRIPTION OF PARAMETERS                                      SADD 130
C           A - NAME OF INPUT MATRIX                                    SADD 140
C           C - SCALAR                                                  SADD 150
C           R - NAME OF OUTPUT MATRIX                                   SADD 160
C           N - NUMBER OF ROWS IN MATRIX A AND R                        SADD 170
C           M - NUMBER OF COLUMNS IN MATRIX A AND R                     SADD 180
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A (AND R) SADD 190
C                  0 - GENERAL                                          SADD 200
C                  1 - SYMMETRIC                                        SADD 210
C                  2 - DIAGONAL                                         SADD 220
C                                                                       SADD 230
C        REMARKS                                                        SADD 240
C           NONE                                                        SADD 250
C                                                                       SADD 260
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  SADD 270
C           LOC                                                         SADD 280
C                                                                       SADD 290
C        METHOD                                                         SADD 300
C           SCALAR IS ADDED TO EACH ELEMENT OF MATRIX                   SADD 310
C                                                                       SADD 320
C     ..................................................................SADD 330
C                                                                       SADD 340
      SUBROUTINE SADD(A,C,R,N,M,MS)                                     SADD 350
      DIMENSION A(1),R(1)                                               SADD 360
C                                                                       SADD 370
C        COMPUTE VECTOR LENGTH, IT                                      SADD 380
C                                                                       SADD 390
      CALL LOC(N,M,IT,N,M,MS)                                           SADD 400
C                                                                       SADD 410
C        ADD SCALAR                                                     SADD 420
C                                                                       SADD 430
      DO 1 I=1,IT                                                       SADD 440
    1 R(I)=A(I)+C                                                       SADD 450
      RETURN                                                            SADD 460
      END                                                               SADD 470
