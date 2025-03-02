C                                                                       RCPY  10
C     ..................................................................RCPY  20
C                                                                       RCPY  30
C        SUBROUTINE RCPY                                                RCPY  40
C                                                                       RCPY  50
C        PURPOSE                                                        RCPY  60
C           COPY SPECIFIED ROW OF A MATRIX INTO A VECTOR                RCPY  70
C                                                                       RCPY  80
C        USAGE                                                          RCPY  90
C           CALL RCPY (A,L,R,N,M,MS)                                    RCPY 100
C                                                                       RCPY 110
C        DESCRIPTION OF PARAMETERS                                      RCPY 120
C           A - NAME OF INPUT MATRIX                                    RCPY 130
C           L - ROW OF A TO BE MOVED TO R                               RCPY 140
C           R - NAME OF OUTPUT VECTOR OF LENGTH M                       RCPY 150
C           N - NUMBER OR ROWS IN A                                     RCPY 160
C           M - NUMBER OF COLUMNS IN A                                  RCPY 170
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         RCPY 180
C                  0 - GENERAL                                          RCPY 190
C                  1 - SYMMETRIC                                        RCPY 200
C                  2 - DIAGONAL                                         RCPY 210
C                                                                       RCPY 220
C        REMARKS                                                        RCPY 230
C           NONE                                                        RCPY 240
C                                                                       RCPY 250
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  RCPY 260
C           LOC                                                         RCPY 270
C                                                                       RCPY 280
C        METHOD                                                         RCPY 290
C           ELEMENTS OF ROW L ARE MOVED TO CORRESPONDING POSITIONS      RCPY 300
C           OF VECTOR R                                                 RCPY 310
C                                                                       RCPY 320
C     ..................................................................RCPY 330
C                                                                       RCPY 340
      SUBROUTINE RCPY(A,L,R,N,M,MS)                                     RCPY 350
      DIMENSION A(1),R(1)                                               RCPY 360
C                                                                       RCPY 370
      DO 3 J=1,M                                                        RCPY 380
C                                                                       RCPY 390
C        LOCATE ELEMENT FOR ANY MATRIX STORAGE MODE                     RCPY 400
C                                                                       RCPY 410
      CALL LOC(L,J,LJ,N,M,MS)                                           RCPY 420
C                                                                       RCPY 430
C        TEST FOR ZERO ELEMENT IN DIAGONAL MATRIX                       RCPY 440
C                                                                       RCPY 450
      IF(LJ) 1,2,1                                                      RCPY 460
C                                                                       RCPY 470
C        MOVE ELEMENT TO R                                              RCPY 480
C                                                                       RCPY 490
    1 R(J)=A(LJ)                                                        RCPY 500
      GO TO 3                                                           RCPY 510
    2 R(J)=0.0                                                          RCPY 520
    3 CONTINUE                                                          RCPY 530
      RETURN                                                            RCPY 540
      END                                                               RCPY 550
