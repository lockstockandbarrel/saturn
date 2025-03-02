C                                                                       DEPS  10
C     ..................................................................DEPS  20
C                                                                       DEPS  30
C        SUBROUTINE DLEPS                                               DEPS  40
C                                                                       DEPS  50
C        PURPOSE                                                        DEPS  60
C           COMPUTES THE VALUE OF AN N-TERM EXPANSION IN LEGENDRE       DEPS  70
C           POLYNOMIALS WITH COEFFICIENT VECTOR C FOR ARGUMENT VALUE X. DEPS  80
C                                                                       DEPS  90
C        USAGE                                                          DEPS 100
C           CALL DLEPS(Y,X,C,N)                                         DEPS 110
C                                                                       DEPS 120
C        DESCRIPTION OF PARAMETERS                                      DEPS 130
C           Y     - RESULT VALUE                                        DEPS 140
C                   DOUBLE PRECISION VARIABLE                           DEPS 150
C           X     - ARGUMENT VALUE                                      DEPS 160
C                   DOUBLE PRECISION VARIABLE                           DEPS 170
C           C     - COEFFICIENT VECTOR OF GIVEN EXPANSION               DEPS 180
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           DEPS 190
C                   DOUBLE PRECISION VECTOR                             DEPS 200
C           N     - DIMENSION OF COEFFICIENT VECTOR C                   DEPS 210
C                                                                       DEPS 220
C        REMARKS                                                        DEPS 230
C           OPERATION IS BYPASSED IN CASE N LESS THAN 1                 DEPS 240
C                                                                       DEPS 250
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DEPS 260
C           NONE                                                        DEPS 270
C                                                                       DEPS 280
C        METHOD                                                         DEPS 290
C           DEFINITION                                                  DEPS 300
C           Y=SUM(C(I)*P(I-1,X), SUMMED OVER I FROM 1 TO N).            DEPS 310
C           EVALUATION IS DONE BY MEANS OF UPWARD RECURSION             DEPS 320
C           USING THE RECURRENCE EQUATION FOR LEGENDRE POLYNOMIALS      DEPS 330
C           P(N+1,X)=2*X*P(N,X)-P(N-1,X)-(X*P(N,X)-P(N-1,X))/(N+1).     DEPS 340
C                                                                       DEPS 350
C     ..................................................................DEPS 360
C                                                                       DEPS 370
      SUBROUTINE DLEPS(Y,X,C,N)                                         DEPS 380
C                                                                       DEPS 390
      DIMENSION C(1)                                                    DEPS 400
      DOUBLE PRECISION C,Y,X,H0,H1,H2                                   DEPS 410
C                                                                       DEPS 420
C        TEST OF DIMENSION                                              DEPS 430
      IF(N)1,1,2                                                        DEPS 440
    1 RETURN                                                            DEPS 450
C                                                                       DEPS 460
    2 Y=C(1)                                                            DEPS 470
      IF(N-2)1,3,3                                                      DEPS 480
C                                                                       DEPS 490
C        INITIALIZATION                                                 DEPS 500
    3 H0=1.D0                                                           DEPS 510
      H1=X                                                              DEPS 520
C                                                                       DEPS 530
      DO 4 I=2,N                                                        DEPS 540
      H2=X*H1                                                           DEPS 550
      H2=H2-H0+H2-(H2-H0)/DFLOAT(I)                                     DEPS 560
      H0=H1                                                             DEPS 570
      H1=H2                                                             DEPS 580
    4 Y=Y+C(I)*H0                                                       DEPS 590
      RETURN                                                            DEPS 600
      END                                                               DEPS 610
