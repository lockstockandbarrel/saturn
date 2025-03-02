C                                                                       SCLA  10
C     ..................................................................SCLA  20
C                                                                       SCLA  30
C        SUBROUTINE SCLA                                                SCLA  40
C                                                                       SCLA  50
C        PURPOSE                                                        SCLA  60
C           SET EACH ELEMENT OF A MATRIX EQUAL TO A GIVEN SCALAR        SCLA  70
C                                                                       SCLA  80
C        USAGE                                                          SCLA  90
C           CALL SCLA (A,C,N,M,MS)                                      SCLA 100
C                                                                       SCLA 110
C        DESCRIPTION OF PARAMETERS                                      SCLA 120
C           A - NAME OF INPUT MATRIX                                    SCLA 130
C           C - SCALAR                                                  SCLA 140
C           N - NUMBER OF ROWS IN MATRIX A                              SCLA 150
C           M - NUMBER OF COLUMNS IN MATRIX A                           SCLA 160
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         SCLA 170
C                  0 - GENERAL                                          SCLA 180
C                  1 - SYMMETRIC                                        SCLA 190
C                  2 - DIAGONAL                                         SCLA 200
C                                                                       SCLA 210
C        REMARKS                                                        SCLA 220
C           NONE                                                        SCLA 230
C                                                                       SCLA 240
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  SCLA 250
C           LOC                                                         SCLA 260
C                                                                       SCLA 270
C        METHOD                                                         SCLA 280
C           EACH ELEMENT OF MATRIX A IS REPLACED BY SCALAR C            SCLA 290
C                                                                       SCLA 300
C     ..................................................................SCLA 310
C                                                                       SCLA 320
      SUBROUTINE SCLA(A,C,N,M,MS)                                       SCLA 330
      DIMENSION A(1)                                                    SCLA 340
C                                                                       SCLA 350
C        COMPUTE VECTOR LENGTH, IT                                      SCLA 360
C                                                                       SCLA 370
      CALL LOC(N,M,IT,N,M,MS)                                           SCLA 380
C                                                                       SCLA 390
C        REPLACE BY SCALAR                                              SCLA 400
C                                                                       SCLA 410
      DO 1 I=1,IT                                                       SCLA 420
    1 A(I)=C                                                            SCLA 430
      RETURN                                                            SCLA 440
      END                                                               SCLA 450
