C                                                                       LAPS  10
C     ..................................................................LAPS  20
C                                                                       LAPS  30
C        SUBROUTINE LAPS                                                LAPS  40
C                                                                       LAPS  50
C        PURPOSE                                                        LAPS  60
C           COMPUTES THE VALUE OF AN N-TERM EXPANSION IN LAGUERRE       LAPS  70
C           POLYNOMIALS WITH COEFFICIENT VECTOR C FOR ARGUMENT VALUE X. LAPS  80
C                                                                       LAPS  90
C        USAGE                                                          LAPS 100
C           CALL LAPS(Y,X,C,N)                                          LAPS 110
C                                                                       LAPS 120
C        DESCRIPTION OF PARAMETERS                                      LAPS 130
C           Y     - RESULT VALUE                                        LAPS 140
C           X     - ARGUMENT VALUE                                      LAPS 150
C           C     - COEFFICIENT VECTOR OF GIVEN EXPANSION               LAPS 160
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           LAPS 170
C           N     - DIMENSION OF COEFFICIENT VECTOR C                   LAPS 180
C                                                                       LAPS 190
C        REMARKS                                                        LAPS 200
C           OPERATION IS BYPASSED IN CASE N LESS THAN 1                 LAPS 210
C                                                                       LAPS 220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  LAPS 230
C           NONE                                                        LAPS 240
C                                                                       LAPS 250
C        METHOD                                                         LAPS 260
C           DEFINITION                                                  LAPS 270
C           Y=SUM(C(I)*L(I-1,X), SUMMED OVER I FROM 1 TO N).            LAPS 280
C           EVALUATION IS DONE BY MEANS OF UPWARD RECURSION             LAPS 290
C           USING THE RECURRENCE EQUATION FOR LAGUERRE POLYNOMIALS      LAPS 300
C           L(N+1,X)=2*L(N,X)-L(N-1,X)-((1+X)*L(N,X)-L(N-1,X))/(N+1).   LAPS 310
C                                                                       LAPS 320
C     ..................................................................LAPS 330
C                                                                       LAPS 340
      SUBROUTINE LAPS(Y,X,C,N)                                          LAPS 350
C                                                                       LAPS 360
      DIMENSION C(1)                                                    LAPS 370
C                                                                       LAPS 380
C        TEST OF DIMENSION                                              LAPS 390
      IF(N)1,1,2                                                        LAPS 400
    1 RETURN                                                            LAPS 410
C                                                                       LAPS 420
    2 Y=C(1)                                                            LAPS 430
      IF(N-2)1,3,3                                                      LAPS 440
C                                                                       LAPS 450
C        INITIALIZATION                                                 LAPS 460
    3 H0=1.                                                             LAPS 470
      H1=1.-X                                                           LAPS 480
      T=1.+X                                                            LAPS 490
C                                                                       LAPS 500
      DO 4 I=2,N                                                        LAPS 510
      H2=H1-H0+H1-(T*H1-H0)/FLOAT(I)                                    LAPS 520
      H0=H1                                                             LAPS 530
      H1=H2                                                             LAPS 540
    4 Y=Y+C(I)*H0                                                       LAPS 550
      RETURN                                                            LAPS 560
      END                                                               LAPS 570
