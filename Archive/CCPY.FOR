C                                                                       CCPY  10
C     ..................................................................CCPY  20
C                                                                       CCPY  30
C        SUBROUTINE CCPY                                                CCPY  40
C                                                                       CCPY  50
C        PURPOSE                                                        CCPY  60
C           COPY SPECIFIED COLUMN OF A MATRIX INTO A VECTOR             CCPY  70
C                                                                       CCPY  80
C        USAGE                                                          CCPY  90
C           CALL CCPY(A,L,R,N,M,MS)                                     CCPY 100
C                                                                       CCPY 110
C        DESCRIPTION OF PARAMETERS                                      CCPY 120
C           A - NAME OF INPUT MATRIX                                    CCPY 130
C           L - COLUMN OF A TO BE MOVED TO R                            CCPY 140
C           R - NAME OF OUTPUT VECTOR OF LENGTH N                       CCPY 150
C           N - NUMBER OR ROWS IN A                                     CCPY 160
C           M - NUMBER OF COLUMNS IN A                                  CCPY 170
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         CCPY 180
C                  0 - GENERAL                                          CCPY 190
C                  1 - SYMMETRIC                                        CCPY 200
C                  2 - DIAGONAL                                         CCPY 210
C                                                                       CCPY 220
C        REMARKS                                                        CCPY 230
C           NONE                                                        CCPY 240
C                                                                       CCPY 250
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  CCPY 260
C           LOC                                                         CCPY 270
C                                                                       CCPY 280
C        METHOD                                                         CCPY 290
C           ELEMENTS OF COLUMN L ARE MOVED TO CORRESPONDING POSITIONS   CCPY 300
C           OF VECTOR R                                                 CCPY 310
C                                                                       CCPY 320
C     ..................................................................CCPY 330
C                                                                       CCPY 340
      SUBROUTINE CCPY(A,L,R,N,M,MS)                                     CCPY 350
      DIMENSION A(1),R(1)                                               CCPY 360
C                                                                       CCPY 370
      DO 3 I=1,N                                                        CCPY 380
C                                                                       CCPY 390
C        LOCATE ELEMENT FOR ANY MATRIX STORAGE MODE                     CCPY 400
C                                                                       CCPY 410
      CALL LOC(I,L,IL,N,M,MS)                                           CCPY 420
C                                                                       CCPY 430
C        TEST FOR ZERO ELEMENT IN DIAGONAL MATRIX                       CCPY 440
C                                                                       CCPY 450
      IF(IL) 1,2,1                                                      CCPY 460
C                                                                       CCPY 470
C        MOVE ELEMENT TO R                                              CCPY 480
C                                                                       CCPY 490
    1 R(I)=A(IL)                                                        CCPY 500
      GO TO 3                                                           CCPY 510
    2 R(I)=0.0                                                          CCPY 520
    3 CONTINUE                                                          CCPY 530
      RETURN                                                            CCPY 540
      END                                                               CCPY 550
