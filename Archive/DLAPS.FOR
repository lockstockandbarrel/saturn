C                                                                       DAPS  10
C     ..................................................................DAPS  20
C                                                                       DAPS  30
C        SUBROUTINE DLAPS                                               DAPS  40
C                                                                       DAPS  50
C        PURPOSE                                                        DAPS  60
C           COMPUTES THE VALUE OF AN N-TERM EXPANSION IN LAGUERRE       DAPS  70
C           POLYNOMIALS WITH COEFFICIENT VECTOR C FOR ARGUMENT VALUE X. DAPS  80
C                                                                       DAPS  90
C        USAGE                                                          DAPS 100
C           CALL DLAPS(Y,X,C,N)                                         DAPS 110
C                                                                       DAPS 120
C        DESCRIPTION OF PARAMETERS                                      DAPS 130
C           Y     - RESULT VALUE                                        DAPS 140
C           X     - ARGUMENT VALUE                                      DAPS 150
C           C     - COEFFICIENT VECTOR OF GIVEN EXPANSION               DAPS 160
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           DAPS 170
C           N     - DIMENSION OF COEFFICIENT VECTOR C                   DAPS 180
C        NOTE: C,Y,X,H0,H1,H2 AND T ARE DOUBLE PRECISION                DAPS 190
C        REMARKS                                                        DAPS 200
C           OPERATION IS BYPASSED IN CASE N LESS THAN 1                 DAPS 210
C                                                                       DAPS 220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DAPS 230
C           NONE                                                        DAPS 240
C                                                                       DAPS 250
C        METHOD                                                         DAPS 260
C           DEFINITION                                                  DAPS 270
C           Y=SUM(C(I)*L(I-1,X), SUMMED OVER I FROM 1 TO N).            DAPS 280
C           EVALUATION IS DONE BY MEANS OF UPWARD RECURSION             DAPS 290
C           USING THE RECURRENCE EQUATION FOR LAGUERRE POLYNOMIALS      DAPS 300
C           L(N+1,X)=2*L(N,X)-L(N-1,X)-((1+X)*L(N,X)-L(N-1,X))/(N+1).   DAPS 310
C                                                                       DAPS 320
C     ..................................................................DAPS 330
C                                                                       DAPS 340
      SUBROUTINE LAPS(Y,X,C,N)                                          DAPS 350
C                                                                       DAPS 360
      DIMENSION C(1)                                                    DAPS 370
      DOUBLE PRECISION C,Y,X,H0,H1,H2,T                                 DAPS 380
C        TEST OF DIMENSION                                              DAPS 390
      IF(N)1,1,2                                                        DAPS 400
    1 RETURN                                                            DAPS 410
C                                                                       DAPS 420
    2 Y=C(1)                                                            DAPS 430
      IF(N-2)1,3,3                                                      DAPS 440
C                                                                       DAPS 450
C        INITIALIZATION                                                 DAPS 460
    3 H0=1.                                                             DAPS 470
      H1=1.-X                                                           DAPS 480
      T=1.+X                                                            DAPS 490
C                                                                       DAPS 500
      DO 4 I=2,N                                                        DAPS 510
      H2=H1-H0+H1-(T*H1-H0)/DFLOAT(I)                                   DAPS 520
      H0=H1                                                             DAPS 530
      H1=H2                                                             DAPS 540
    4 Y=Y+C(I)*H0                                                       DAPS 550
      RETURN                                                            DAPS 560
      END                                                               DAPS 570
