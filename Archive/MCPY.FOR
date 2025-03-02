C                                                                       MCPY  10
C        ...............................................................MCPY  20
C                                                                       MCPY  30
C        SUBROUTINE MCPY                                                MCPY  40
C                                                                       MCPY  50
C        PURPOSE                                                        MCPY  60
C           COPY ENTIRE MATRIX                                          MCPY  70
C                                                                       MCPY  80
C        USAGE                                                          MCPY  90
C           CALL MCPY (A,R,N,M,MS)                                      MCPY 100
C                                                                       MCPY 110
C        DESCRIPTION OF PARAMETERS                                      MCPY 120
C           A - NAME OF INPUT MATRIX                                    MCPY 130
C           R - NAME OF OUTPUT MATRIX                                   MCPY 140
C           N - NUMBER OF ROWS IN A OR R                                MCPY 150
C           M - NUMBER OF COLUMNS IN A OR R                             MCPY 160
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A (AND R) MCPY 170
C                  0 - GENERAL                                          MCPY 180
C                  1 - SYMMETRIC                                        MCPY 190
C                  2 - DIAGONAL                                         MCPY 200
C                                                                       MCPY 210
C        REMARKS                                                        MCPY 220
C           NONE                                                        MCPY 230
C                                                                       MCPY 240
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  MCPY 250
C           LOC                                                         MCPY 260
C                                                                       MCPY 270
C        METHOD                                                         MCPY 280
C           EACH ELEMENT OF MATRIX A IS MOVED TO THE CORRESPONDING      MCPY 290
C           ELEMENT OF MATRIX R                                         MCPY 300
C                                                                       MCPY 310
C     ..................................................................MCPY 320
C                                                                       MCPY 330
      SUBROUTINE MCPY(A,R,N,M,MS)                                       MCPY 340
      DIMENSION A(1),R(1)                                               MCPY 350
C                                                                       MCPY 360
C        COMPUTE VECTOR LENGTH, IT                                      MCPY 370
C                                                                       MCPY 380
      CALL LOC(N,M,IT,N,M,MS)                                           MCPY 390
C                                                                       MCPY 400
C        COPY MATRIX                                                    MCPY 410
C                                                                       MCPY 420
      DO 1 I=1,IT                                                       MCPY 430
    1 R(I)=A(I)                                                         MCPY 440
      RETURN                                                            MCPY 450
      END                                                               MCPY 460
