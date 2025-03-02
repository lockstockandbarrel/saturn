C                                                                       DNPS  10
C     ..................................................................DNPS  20
C                                                                       DNPS  30
C        SUBROUTINE DCNPS                                               DNPS  40
C                                                                       DNPS  50
C        PURPOSE                                                        DNPS  60
C           COMPUTES THE VALUE OF AN N-TERM EXPANSION IN CHEBYSHEV      DNPS  70
C           POLYNOMIALS WITH COEFFICIENT VECTOR C FOR ARGUMENT VALUE X. DNPS  80
C                                                                       DNPS  90
C        USAGE                                                          DNPS 100
C           CALL DCNPS(Y,X,C,N)                                         DNPS 110
C                                                                       DNPS 120
C        DESCRIPTION OF PARAMETERS                                      DNPS 130
C           Y     - RESULT VALUE                                        DNPS 140
C                   DOUBLE PRECISION VARIABLE                           DNPS 150
C           X     - ARGUMENT VALUE                                      DNPS 160
C                   DOUBLE PRECISION VARIABLE                           DNPS 170
C           C     - COEFFICIENT VECTOR OF GIVEN EXPANSION               DNPS 180
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           DNPS 190
C                   DOUBLE PRECISION VECTOR                             DNPS 200
C           N     - DIMENSION OF COEFFICIENT VECTOR C                   DNPS 210
C                                                                       DNPS 220
C        REMARKS                                                        DNPS 230
C           OPERATION IS BYPASSED IN CASE N LESS THAN 1                 DNPS 240
C                                                                       DNPS 250
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DNPS 260
C           NONE                                                        DNPS 270
C                                                                       DNPS 280
C        METHOD                                                         DNPS 290
C           DEFINITION                                                  DNPS 300
C           Y=SUM(C(I)*T(I-1,X), SUMMED OVER I FROM 1 TO N).            DNPS 310
C           EVALUATION IS DONE BY MEANS OF BACKWARD RECURSION           DNPS 320
C           USING THE RECURRENCE EQUATION FOR CHEBYSHEV POLYNOMIALS     DNPS 330
C           T(N+1,X)=2*X*T(N,X)-T(N-1,X).                               DNPS 340
C                                                                       DNPS 350
C     ..................................................................DNPS 360
C                                                                       DNPS 370
      SUBROUTINE DCNPS(Y,X,C,N)                                         DNPS 380
C                                                                       DNPS 390
      DIMENSION C(1)                                                    DNPS 400
      DOUBLE PRECISION C,Y,X,H0,H1,H2,ARG                               DNPS 410
C                                                                       DNPS 420
C        TEST OF DIMENSION                                              DNPS 430
      IF(N)1,1,2                                                        DNPS 440
    1 RETURN                                                            DNPS 450
C                                                                       DNPS 460
    2 IF(N-2)3,4,4                                                      DNPS 470
    3 Y=C(1)                                                            DNPS 480
      RETURN                                                            DNPS 490
C                                                                       DNPS 500
C        INITIALIZATION                                                 DNPS 510
    4 ARG=X+X                                                           DNPS 520
      H1=0.D0                                                           DNPS 530
      H0=0.D0                                                           DNPS 540
C                                                                       DNPS 550
      DO 5 I=1,N                                                        DNPS 560
      K=N-I                                                             DNPS 570
      H2=H1                                                             DNPS 580
      H1=H0                                                             DNPS 590
    5 H0=ARG*H1-H2+C(K+1)                                               DNPS 600
      Y=0.5D0*(C(1)-H2+H0)                                              DNPS 610
      RETURN                                                            DNPS 620
      END                                                               DNPS 630
