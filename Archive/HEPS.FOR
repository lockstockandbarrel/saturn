C                                                                       HEPS  10
C     ..................................................................HEPS  20
C                                                                       HEPS  30
C        SUBROUTINE HEPS                                                HEPS  40
C                                                                       HEPS  50
C        PURPOSE                                                        HEPS  60
C           COMPUTES THE VALUE OF AN N-TERM EXPANSION IN HERMITE        HEPS  70
C           POLYNOMIALS WITH COEFFICIENT VECTOR C FOR ARGUMENT VALUE X. HEPS  80
C                                                                       HEPS  90
C        USAGE                                                          HEPS 100
C           CALL HEPS(Y,X,C,N)                                          HEPS 110
C                                                                       HEPS 120
C        DESCRIPTION OF PARAMETERS                                      HEPS 130
C           Y     - RESULT VALUE                                        HEPS 140
C           X     - ARGUMENT VALUE                                      HEPS 150
C           C     - COEFFICIENT VECTOR OF GIVEN EXPANSION               HEPS 160
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           HEPS 170
C           N     - DIMENSION OF COEFFICIENT VECTOR C                   HEPS 180
C                                                                       HEPS 190
C        REMARKS                                                        HEPS 200
C           OPERATION IS BYPASSED IN CASE N LESS THAN 1                 HEPS 210
C                                                                       HEPS 220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  HEPS 230
C           NONE                                                        HEPS 240
C                                                                       HEPS 250
C        METHOD                                                         HEPS 260
C           DEFINITION                                                  HEPS 270
C           Y=SUM(C(I)*H(I-1,X), SUMMED OVER I FROM 1 TO N).            HEPS 280
C           EVALUATION IS DONE BY MEANS OF UPWARD RECURSION             HEPS 290
C           USING THE RECURRENCE EQUATION FOR HERMITE POLYNOMIALS       HEPS 300
C           H(N+1,X)=2*(X*H(N,X)-N*H(N-1,X)).                           HEPS 310
C                                                                       HEPS 320
C     ..................................................................HEPS 330
C                                                                       HEPS 340
      SUBROUTINE HEPS(Y,X,C,N)                                          HEPS 350
C                                                                       HEPS 360
      DIMENSION C(1)                                                    HEPS 370
C                                                                       HEPS 380
C        TEST OF DIMENSION                                              HEPS 390
      IF(N)1,1,2                                                        HEPS 400
    1 RETURN                                                            HEPS 410
C                                                                       HEPS 420
    2 Y=C(1)                                                            HEPS 430
      IF(N-2)1,3,3                                                      HEPS 440
C                                                                       HEPS 450
C        INITIALIZATION                                                 HEPS 460
    3 H0=1.                                                             HEPS 470
      H1=X+X                                                            HEPS 480
C                                                                       HEPS 490
      DO 4 I=2,N                                                        HEPS 500
      H2=X*H1-FLOAT(I-1)*H0                                             HEPS 510
      H0=H1                                                             HEPS 520
      H1=H2+H2                                                          HEPS 530
    4 Y=Y+C(I)*H0                                                       HEPS 540
      RETURN                                                            HEPS 550
      END                                                               HEPS 560
