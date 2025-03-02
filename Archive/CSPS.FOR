C                                                                       CSPS  10
C     ..................................................................CSPS  20
C                                                                       CSPS  30
C        SUBROUTINE CSPS                                                CSPS  40
C                                                                       CSPS  50
C        PURPOSE                                                        CSPS  60
C           COMPUTES THE VALUE OF AN N-TERM EXPANSION IN SHIFTED        CSPS  70
C           CHEBYSHEV POLYNOMIALS WITH COEFFICIENT VECTOR C             CSPS  80
C           FOR ARGUMENT VALUE X.                                       CSPS  90
C                                                                       CSPS 100
C        USAGE                                                          CSPS 110
C           CALL CSPS(Y,X,C,N)                                          CSPS 120
C                                                                       CSPS 130
C        DESCRIPTION OF PARAMETERS                                      CSPS 140
C           Y     - RESULT VALUE                                        CSPS 150
C           X     - ARGUMENT VALUE                                      CSPS 160
C           C     - COEFFICIENT VECTOR OF GIVEN EXPANSION               CSPS 170
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           CSPS 180
C           N     - DIMENSION OF COEFFICIENT VECTOR C                   CSPS 190
C                                                                       CSPS 200
C        REMARKS                                                        CSPS 210
C           OPERATION IS BYPASSED IN CASE N LESS THAN 1                 CSPS 220
C                                                                       CSPS 230
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  CSPS 240
C           NONE                                                        CSPS 250
C                                                                       CSPS 260
C        METHOD                                                         CSPS 270
C           DEFINITION                                                  CSPS 280
C           Y=SUM(C(I)*TS(I-1,X), SUMMED OVER I FROM 1 TO N).           CSPS 290
C           EVALUATION IS DONE BY MEANS OF BACKWARD RECURSION           CSPS 300
C           USING THE RECURRENCE EQUATION FOR SHIFTED                   CSPS 310
C           CHEBYSHEV POLYNOMIALS                                       CSPS 320
C           TS(N+1,X)=(4*X-2)*TS(N,X)-TS(N-1,X).                        CSPS 330
C                                                                       CSPS 340
C     ..................................................................CSPS 350
C                                                                       CSPS 360
      SUBROUTINE CSPS(Y,X,C,N)                                          CSPS 370
C                                                                       CSPS 380
      DIMENSION C(1)                                                    CSPS 390
C                                                                       CSPS 400
C        TEST OF DIMENSION                                              CSPS 410
      IF(N)1,1,2                                                        CSPS 420
    1 RETURN                                                            CSPS 430
C                                                                       CSPS 440
    2 IF(N-2)3,4,4                                                      CSPS 450
    3 Y=C(1)                                                            CSPS 460
      RETURN                                                            CSPS 470
C                                                                       CSPS 480
C        INITIALIZATION                                                 CSPS 490
    4 ARG=X+X-1.                                                        CSPS 500
      ARG=ARG+ARG                                                       CSPS 510
      H1=0.                                                             CSPS 520
      H0=0.                                                             CSPS 530
C                                                                       CSPS 540
      DO 5 I=1,N                                                        CSPS 550
      K=N-I                                                             CSPS 560
      H2=H1                                                             CSPS 570
      H1=H0                                                             CSPS 580
    5 H0=ARG*H1-H2+C(K+1)                                               CSPS 590
      Y=0.5*(C(1)-H2+H0)                                                CSPS 600
      RETURN                                                            CSPS 610
      END                                                               CSPS 620
