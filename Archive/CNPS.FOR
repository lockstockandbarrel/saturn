C                                                                       CNPS  10
C     ..................................................................CNPS  20
C                                                                       CNPS  30
C        SUBROUTINE CNPS                                                CNPS  40
C                                                                       CNPS  50
C        PURPOSE                                                        CNPS  60
C           COMPUTES THE VALUE OF AN N-TERM EXPANSION IN CHEBYSHEV      CNPS  70
C           POLYNOMIALS WITH COEFFICIENT VECTOR C FOR ARGUMENT VALUE X. CNPS  80
C                                                                       CNPS  90
C        USAGE                                                          CNPS 100
C           CALL CNPS(Y,X,C,N)                                          CNPS 110
C                                                                       CNPS 120
C        DESCRIPTION OF PARAMETERS                                      CNPS 130
C           Y     - RESULT VALUE                                        CNPS 140
C           X     - ARGUMENT VALUE                                      CNPS 150
C           C     - COEFFICIENT VECTOR OF GIVEN EXPANSION               CNPS 160
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           CNPS 170
C           N     - DIMENSION OF COEFFICIENT VECTOR C                   CNPS 180
C                                                                       CNPS 190
C        REMARKS                                                        CNPS 200
C           OPERATION IS BYPASSED IN CASE N LESS THAN 1                 CNPS 210
C                                                                       CNPS 220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  CNPS 230
C           NONE                                                        CNPS 240
C                                                                       CNPS 250
C        METHOD                                                         CNPS 260
C           DEFINITION                                                  CNPS 270
C           Y=SUM(C(I)*T(I-1,X), SUMMED OVER I FROM 1 TO N).            CNPS 280
C           EVALUATION IS DONE BY MEANS OF BACKWARD RECURSION           CNPS 290
C           USING THE RECURRENCE EQUATION FOR CHEBYSHEV POLYNOMIALS     CNPS 300
C           T(N+1,X)=2*X*T(N,X)-T(N-1,X).                               CNPS 310
C                                                                       CNPS 320
C     ..................................................................CNPS 330
C                                                                       CNPS 340
      SUBROUTINE CNPS(Y,X,C,N)                                          CNPS 350
C                                                                       CNPS 360
      DIMENSION C(1)                                                    CNPS 370
C                                                                       CNPS 380
C        TEST OF DIMENSION                                              CNPS 390
      IF(N)1,1,2                                                        CNPS 400
    1 RETURN                                                            CNPS 410
C                                                                       CNPS 420
    2 IF(N-2)3,4,4                                                      CNPS 430
    3 Y=C(1)                                                            CNPS 440
      RETURN                                                            CNPS 450
C                                                                       CNPS 460
C        INITIALIZATION                                                 CNPS 470
    4 ARG=X+X                                                           CNPS 480
      H1=0.                                                             CNPS 490
      H0=0.                                                             CNPS 500
C                                                                       CNPS 510
      DO 5 I=1,N                                                        CNPS 520
      K=N-I                                                             CNPS 530
      H2=H1                                                             CNPS 540
      H1=H0                                                             CNPS 550
    5 H0=ARG*H1-H2+C(K+1)                                               CNPS 560
      Y=0.5*(C(1)-H2+H0)                                                CNPS 570
      RETURN                                                            CNPS 580
      END                                                               CNPS 590
