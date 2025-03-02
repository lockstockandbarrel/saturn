C                                                                       SMPY  10
C     ..................................................................SMPY  20
C                                                                       SMPY  30
C        SUBROUTINE SMPY                                                SMPY  40
C                                                                       SMPY  50
C        PURPOSE                                                        SMPY  60
C           MULTIPLY EACH ELEMENT OF A MATRIX BY A SCALAR TO FORM A     SMPY  70
C           RESULTANT MATRIX                                            SMPY  80
C                                                                       SMPY  90
C        USAGE                                                          SMPY 100
C           CALL SMPY(A,C,R,N,M,MS)                                     SMPY 110
C                                                                       SMPY 120
C        DESCRIPTION OF PARAMETERS                                      SMPY 130
C           A - NAME OF INPUT MATRIX                                    SMPY 140
C           C - SCALAR                                                  SMPY 150
C           R - NAME OF OUTPUT MATRIX                                   SMPY 160
C           N - NUMBER OF ROWS IN MATRIX A AND R                        SMPY 170
C           M - NUMBER OF COLUMNS IN MATRIX A AND R                     SMPY 180
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A (AND R) SMPY 190
C                  0 - GENERAL                                          SMPY 200
C                  1 - SYMMETRIC                                        SMPY 210
C                  2 - DIAGONAL                                         SMPY 220
C                                                                       SMPY 230
C        REMARKS                                                        SMPY 240
C           NONE                                                        SMPY 250
C                                                                       SMPY 260
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  SMPY 270
C           LOC                                                         SMPY 280
C                                                                       SMPY 290
C        METHOD                                                         SMPY 300
C           SCALAR IS MULTIPLIED BY EACH ELEMENT OF MATRIX              SMPY 310
C                                                                       SMPY 320
C     ..................................................................SMPY 330
C                                                                       SMPY 340
      SUBROUTINE SMPY(A,C,R,N,M,MS)                                     SMPY 350
      DIMENSION A(1),R(1)                                               SMPY 360
C                                                                       SMPY 370
C        COMPUTE VECTOR LENGTH, IT                                      SMPY 380
C                                                                       SMPY 390
      CALL LOC(N,M,IT,N,M,MS)                                           SMPY 400
C                                                                       SMPY 410
C        MULTIPLY BY SCALAR                                             SMPY 420
C                                                                       SMPY 430
      DO 1 I=1,IT                                                       SMPY 440
    1 R(I)=A(I)*C                                                       SMPY 450
      RETURN                                                            SMPY 460
      END                                                               SMPY 470
