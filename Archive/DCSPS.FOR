C                                                                       DSPS  10
C     ..................................................................DSPS  20
C                                                                       DSPS  30
C        SUBROUTINE DCSPS                                               DSPS  40
C                                                                       DSPS  50
C        PURPOSE                                                        DSPS  60
C           COMPUTES THE VALUE OF AN N-TERM EXPANSION IN SHIFTED        DSPS  70
C           CHEBYSHEV POLYNOMIALS WITH COEFFICIENT VECTOR C             DSPS  80
C           FOR ARGUMENT VALUE X.                                       DSPS  90
C                                                                       DSPS 100
C        USAGE                                                          DSPS 110
C           CALL DCSPS(Y,X,C,N)                                         DSPS 120
C                                                                       DSPS 130
C        DESCRIPTION OF PARAMETERS                                      DSPS 140
C           Y     - RESULT VALUE                                        DSPS 150
C                   DOUBLE PRECISION VARIABLE                           DSPS 160
C           X     - ARGUMENT VALUE                                      DSPS 170
C                   DOUBLE PRECISION VARIABLE                           DSPS 180
C           C     - COEFFICIENT VECTOR OF GIVEN EXPANSION               DSPS 190
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           DSPS 200
C                   DOUBLE PRECISION VECTOR                             DSPS 210
C           N     - DIMENSION OF COEFFICIENT VECTOR C                   DSPS 220
C                                                                       DSPS 230
C        REMARKS                                                        DSPS 240
C           OPERATION IS BYPASSED IN CASE N LESS THAN 1                 DSPS 250
C                                                                       DSPS 260
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DSPS 270
C           NONE                                                        DSPS 280
C                                                                       DSPS 290
C        METHOD                                                         DSPS 300
C           DEFINITION                                                  DSPS 310
C           Y=SUM(C(I)*TS(I-1,X), SUMMED OVER I FROM 1 TO N).           DSPS 320
C           EVALUATION IS DONE BY MEANS OF BACKWARD RECURSION           DSPS 330
C           USING THE RECURRENCE EQUATION FOR SHIFTED                   DSPS 340
C           CHEBYSHEV POLYNOMIALS                                       DSPS 350
C           TS(N+1,X)=(4*X-2)*TS(N,X)-TS(N-1,X).                        DSPS 360
C                                                                       DSPS 370
C     ..................................................................DSPS 380
C                                                                       DSPS 390
      SUBROUTINE DCSPS(Y,X,C,N)                                         DSPS 400
C                                                                       DSPS 410
      DIMENSION C(1)                                                    DSPS 420
      DOUBLE PRECISION C,Y,X,H0,H1,H2,ARG                               DSPS 430
C                                                                       DSPS 440
C        TEST OF DIMENSION                                              DSPS 450
      IF(N)1,1,2                                                        DSPS 460
    1 RETURN                                                            DSPS 470
C                                                                       DSPS 480
    2 IF(N-2)3,4,4                                                      DSPS 490
    3 Y=C(1)                                                            DSPS 500
      RETURN                                                            DSPS 510
C                                                                       DSPS 520
C        INITIALIZATION                                                 DSPS 530
    4 ARG=X+X-1.D0                                                      DSPS 540
      ARG=ARG+ARG                                                       DSPS 550
      H1=0.D0                                                           DSPS 560
      H0=0.D0                                                           DSPS 570
      DO 5 I=1,N                                                        DSPS 580
      K=N-I                                                             DSPS 590
      H2=H1                                                             DSPS 600
      H1=H0                                                             DSPS 610
    5 H0=ARG*H1-H2+C(K+1)                                               DSPS 620
      Y=0.5D0*(C(1)-H2+H0)                                              DSPS 630
      RETURN                                                            DSPS 640
      END                                                               DSPS 650
