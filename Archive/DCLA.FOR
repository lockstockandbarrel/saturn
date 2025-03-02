C                                                                       DCLA  10
C     ..................................................................DCLA  20
C                                                                       DCLA  30
C        SUBROUTINE DCLA                                                DCLA  40
C                                                                       DCLA  50
C        PURPOSE                                                        DCLA  60
C           SET EACH DIAGONAL ELEMENT OF A MATRIX EQUAL TO A SCALAR     DCLA  70
C                                                                       DCLA  80
C        USAGE                                                          DCLA  90
C           CALL DCLA (A,C,N,MS)                                        DCLA 100
C                                                                       DCLA 110
C        DESCRIPTION OF PARAMETERS                                      DCLA 120
C           A - NAME OF INPUT MATRIX                                    DCLA 130
C           C - SCALAR                                                  DCLA 140
C           N - NUMBER OF ROWS AND COLUMNS IN MATRIX A                  DCLA 150
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         DCLA 160
C                  0 - GENERAL                                          DCLA 170
C                  1 - SYMMETRIC                                        DCLA 180
C                  2 - DIAGONAL                                         DCLA 190
C                                                                       DCLA 200
C        REMARKS                                                        DCLA 210
C           INPUT MATRIX MUST BE A SQUARE MATRIX                        DCLA 220
C                                                                       DCLA 230
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DCLA 240
C           LOC                                                         DCLA 250
C                                                                       DCLA 260
C        METHOD                                                         DCLA 270
C           EACH ELEMENT ON DIAGONAL OF MATRIX IS REPLACED BY SCALAR C  DCLA 280
C                                                                       DCLA 290
C     ..................................................................DCLA 300
C                                                                       DCLA 310
      SUBROUTINE DCLA(A,C,N,MS)                                         DCLA 320
      DIMENSION A(1)                                                    DCLA 330
C                                                                       DCLA 340
      DO 3 I=1,N                                                        DCLA 350
C                                                                       DCLA 360
C        LOCATE DIAGONAL ELEMENT FOR ANY MATRIX STORAGE MODE            DCLA 370
C                                                                       DCLA 380
      CALL LOC(I,I,ID,N,N,MS)                                           DCLA 390
C                                                                       DCLA 400
C        REPLACE DIAGONAL ELEMENTS                                      DCLA 410
C                                                                       DCLA 420
    3 A(ID)=C                                                           DCLA 430
      RETURN                                                            DCLA 440
      END                                                               DCLA 450
