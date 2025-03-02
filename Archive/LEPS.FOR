C                                                                       LEPS  10
C     ..................................................................LEPS  20
C                                                                       LEPS  30
C        SUBROUTINE LEPS                                                LEPS  40
C                                                                       LEPS  50
C        PURPOSE                                                        LEPS  60
C           COMPUTES THE VALUE OF AN N-TERM EXPANSION IN LEGENDRE       LEPS  70
C           POLYNOMIALS WITH COEFFICIENT VECTOR C FOR ARGUMENT VALUE X. LEPS  80
C                                                                       LEPS  90
C        USAGE                                                          LEPS 100
C           CALL LEPS(Y,X,C,N)                                          LEPS 110
C                                                                       LEPS 120
C        DESCRIPTION OF PARAMETERS                                      LEPS 130
C           Y     - RESULT VALUE                                        LEPS 140
C           X     - ARGUMENT VALUE                                      LEPS 150
C           C     - COEFFICIENT VECTOR OF GIVEN EXPANSION               LEPS 160
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           LEPS 170
C           N     - DIMENSION OF COEFFICIENT VECTOR C                   LEPS 180
C                                                                       LEPS 190
C        REMARKS                                                        LEPS 200
C           OPERATION IS BYPASSED IN CASE N LESS THAN 1                 LEPS 210
C                                                                       LEPS 220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  LEPS 230
C           NONE                                                        LEPS 240
C                                                                       LEPS 250
C        METHOD                                                         LEPS 260
C           DEFINITION                                                  LEPS 270
C           Y=SUM(C(I)*P(I-1,X), SUMMED OVER I FROM 1 TO N).            LEPS 280
C           EVALUATION IS DONE BY MEANS OF UPWARD RECURSION             LEPS 290
C           USING THE RECURRENCE EQUATION FOR LEGENDRE POLYNOMIALS      LEPS 300
C           P(N+1,X)=2*X*P(N,X)-P(N-1,X)-(X*P(N,X)-P(N-1,X))/(N+1).     LEPS 310
C                                                                       LEPS 320
C     ..................................................................LEPS 330
C                                                                       LEPS 340
      SUBROUTINE LEPS(Y,X,C,N)                                          LEPS 350
C                                                                       LEPS 360
      DIMENSION C(1)                                                    LEPS 370
C                                                                       LEPS 380
C        TEST OF DIMENSION                                              LEPS 390
      IF(N)1,1,2                                                        LEPS 400
    1 RETURN                                                            LEPS 410
C                                                                       LEPS 420
    2 Y=C(1)                                                            LEPS 430
      IF(N-2)1,3,3                                                      LEPS 440
C                                                                       LEPS 450
C        INITIALIZATION                                                 LEPS 460
    3 H0=1.                                                             LEPS 470
      H1=X                                                              LEPS 480
C                                                                       LEPS 490
      DO 4 I=2,N                                                        LEPS 500
      H2=X*H1                                                           LEPS 510
      H2=H2-H0+H2-(H2-H0)/FLOAT(I)                                      LEPS 520
      H0=H1                                                             LEPS 530
      H1=H2                                                             LEPS 540
    4 Y=Y+C(I)*H0                                                       LEPS 550
      RETURN                                                            LEPS 560
      END                                                               LEPS 570
