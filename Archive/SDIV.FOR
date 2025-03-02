C                                                                       SDIV  10
C     ..................................................................SDIV  20
C                                                                       SDIV  30
C        SUBROUTINE SDIV                                                SDIV  40
C                                                                       SDIV  50
C        PURPOSE                                                        SDIV  60
C           DIVIDE EACH ELEMENT OF A MATRIX BY A SCALAR TO FORM A       SDIV  70
C           RESULTANT MATRIX                                            SDIV  80
C                                                                       SDIV  90
C        USAGE                                                          SDIV 100
C           CALL SDIV(A,C,R,N,M,MS)                                     SDIV 110
C                                                                       SDIV 120
C        DESCRIPTION OF PARAMETERS                                      SDIV 130
C           A - NAME OF INPUT MATRIX                                    SDIV 140
C           C - SCALAR                                                  SDIV 150
C           R - NAME OF OUTPUT MATRIX                                   SDIV 160
C           N - NUMBER OF ROWS IN MATRIX A AND R                        SDIV 170
C           M - NUMBER OF COLUMNS IN MATRIX A AND R                     SDIV 180
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A (AND R) SDIV 190
C                  0 - GENERAL                                          SDIV 200
C                  1 - SYMMETRIC                                        SDIV 210
C                  2 - DIAGONAL                                         SDIV 220
C                                                                       SDIV 230
C        REMARKS                                                        SDIV 240
C           IF SCALAR IS ZERO, DIVISION IS PERFORMED ONLY ONCE TO CAUSE SDIV 250
C           FLOATING POINT OVERFLOW CONDITION                           SDIV 260
C                                                                       SDIV 270
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  SDIV 280
C           LOC                                                         SDIV 290
C                                                                       SDIV 300
C        METHOD                                                         SDIV 310
C           EACH ELEMENT OF MATRIX IS DIVIDED BY SCALAR                 SDIV 320
C                                                                       SDIV 330
C     ..................................................................SDIV 340
C                                                                       SDIV 350
      SUBROUTINE SDIV(A,C,R,N,M,MS)                                     SDIV 360
      DIMENSION A(1),R(1)                                               SDIV 370
C                                                                       SDIV 380
C        COMPUTE VECTOR LENGTH, IT                                      SDIV 390
C                                                                       SDIV 400
      CALL LOC(N,M,IT,N,M,MS)                                           SDIV 410
C                                                                       SDIV 420
C        DIVIDE BY SCALAR (IF SCALAR IS ZERO, DIVIDE ONLY ONCE)         SDIV 430
C                                                                       SDIV 440
      IF(C) 2,1,2                                                       SDIV 450
    1 IT=1                                                              SDIV 460
    2 DO 3 I=1,IT                                                       SDIV 470
    3 R(I)=A(I)/C                                                       SDIV 480
      RETURN                                                            SDIV 490
      END                                                               SDIV 500
