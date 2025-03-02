C                                                                       DHPS  10
C     ..................................................................DHPS  20
C                                                                       DHPS  30
C        SUBROUTINE DHEPS                                               DHPS  40
C                                                                       DHPS  50
C        PURPOSE                                                        DHPS  60
C           COMPUTES THE VALUE OF AN N-TERM EXPANSION IN HERMITE        DHPS  70
C           POLYNOMIALS WITH COEFFICIENT VECTOR C FOR ARGUMENT VALUE X. DHPS  80
C                                                                       DHPS  90
C        USAGE                                                          DHPS 100
C           CALL DHEPS(Y,X,C,N)                                         DHPS 110
C                                                                       DHPS 120
C        DESCRIPTION OF PARAMETERS                                      DHPS 130
C           Y     - RESULT VALUE                                        DHPS 140
C                   DOUBLE PRECISION VARIABLE                           DHPS 150
C           X     - ARGUMENT VALUE                                      DHPS 160
C                   DOUBLE PRECISION VARIABLE                           DHPS 170
C           C     - COEFFICIENT VECTOR OF GIVEN EXPANSION               DHPS 180
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           DHPS 190
C                   DOUBLE PRECISION VECTOR                             DHPS 200
C           N     - DIMENSION OF COEFFICIENT VECTOR C                   DHPS 210
C                                                                       DHPS 220
C        REMARKS                                                        DHPS 230
C           OPERATION IS BYPASSED IN CASE N LESS THAN 1                 DHPS 240
C                                                                       DHPS 250
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DHPS 260
C           NONE                                                        DHPS 270
C                                                                       DHPS 280
C        METHOD                                                         DHPS 290
C           DEFINITION                                                  DHPS 300
C           Y=SUM(C(I)*H(I-1,X), SUMMED OVER I FROM 1 TO N).            DHPS 310
C           EVALUATION IS DONE BY MEANS OF UPWARD RECURSION             DHPS 320
C           USING THE RECURRENCE EQUATION FOR HERMITE POLYNOMIALS       DHPS 330
C           H(N+1,X)=2*(X*H(N,X)-N*H(N-1,X)).                           DHPS 340
C                                                                       DHPS 350
C     ..................................................................DHPS 360
C                                                                       DHPS 370
      SUBROUTINE DHEPS(Y,X,C,N)                                         DHPS 380
C                                                                       DHPS 390
      DIMENSION C(1)                                                    DHPS 400
      DOUBLE PRECISION C,Y,X,H0,H1,H2                                   DHPS 410
C                                                                       DHPS 420
C        TEST OF DIMENSION                                              DHPS 430
      IF(N)1,1,2                                                        DHPS 440
    1 RETURN                                                            DHPS 450
C                                                                       DHPS 460
    2 Y=C(1)                                                            DHPS 470
      IF(N-2)1,3,3                                                      DHPS 480
C                                                                       DHPS 490
C        INITIALIZATION                                                 DHPS 500
    3 H0=1.D0                                                           DHPS 510
      H1=X+X                                                            DHPS 520
C                                                                       DHPS 530
      DO 4 I=2,N                                                        DHPS 540
      H2=X*H1-DFLOAT(I-1)*H0                                            DHPS 550
      H0=H1                                                             DHPS 560
      H1=H2+H2                                                          DHPS 570
    4 Y=Y+C(I)*H0                                                       DHPS 580
      RETURN                                                            DHPS 590
      END                                                               DHPS 600
