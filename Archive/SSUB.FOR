C                                                                       SSUB  10
C     ..................................................................SSUB  20
C                                                                       SSUB  30
C        SUBROUTINE SSUB                                                SSUB  40
C                                                                       SSUB  50
C        PURPOSE                                                        SSUB  60
C           SUBTRACT A SCALAR FROM EACH ELEMENT OF A MATRIX TO FORM A   SSUB  70
C           RESULTANT MATRIX                                            SSUB  80
C                                                                       SSUB  90
C        USAGE                                                          SSUB 100
C           CALL SSUB(A,C,R,N,M,MS)                                     SSUB 110
C                                                                       SSUB 120
C        DESCRIPTION OF PARAMETERS                                      SSUB 130
C           A - NAME OF INPUT MATRIX                                    SSUB 140
C           C - SCALAR                                                  SSUB 150
C           R - NAME OF OUTPUT MATRIX                                   SSUB 160
C           N - NUMBER OF ROWS IN MATRIX A AND R                        SSUB 170
C           M - NUMBER OF COLUMNS IN MATRIX A AND R                     SSUB 180
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A (AND R) SSUB 190
C                  0 - GENERAL                                          SSUB 200
C                  1 - SYMMETRIC                                        SSUB 210
C                  2 - DIAGONAL                                         SSUB 220
C                                                                       SSUB 230
C        REMARKS                                                        SSUB 240
C           NONE                                                        SSUB 250
C                                                                       SSUB 260
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  SSUB 270
C           LOC                                                         SSUB 280
C                                                                       SSUB 290
C        METHOD                                                         SSUB 300
C           SCALAR IS SUBTRACTED FROM EACH EACH ELEMENT OF MATRIX       SSUB 310
C                                                                       SSUB 320
C     ..................................................................SSUB 330
C                                                                       SSUB 340
      SUBROUTINE SSUB(A,C,R,N,M,MS)                                     SSUB 350
      DIMENSION A(1),R(1)                                               SSUB 360
C                                                                       SSUB 370
C        COMPUTE VECTOR LENGTH, IT                                      SSUB 380
C                                                                       SSUB 390
      CALL LOC(N,M,IT,N,M,MS)                                           SSUB 400
C                                                                       SSUB 410
C        SUBTRACT SCALAR                                                SSUB 420
C                                                                       SSUB 430
      DO 1 I=1,IT                                                       SSUB 440
    1 R(I)=A(I)-C                                                       SSUB 450
      RETURN                                                            SSUB 460
      END                                                               SSUB 470
