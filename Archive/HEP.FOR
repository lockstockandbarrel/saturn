C                                                                       HEP   10
C     ..................................................................HEP   20
C                                                                       HEP   30
C        SUBROUTINE HEP                                                 HEP   40
C                                                                       HEP   50
C        PURPOSE                                                        HEP   60
C           COMPUTE THE VALUES OF THE HERMITE POLYNOMIALS H(N,X)        HEP   70
C           FOR ARGUMENT VALUE X AND ORDERS 0 UP TO N.                  HEP   80
C                                                                       HEP   90
C        USAGE                                                          HEP  100
C           CALL HEP(Y,X,N)                                             HEP  110
C                                                                       HEP  120
C        DESCRIPTION OF PARAMETERS                                      HEP  130
C           Y     - RESULT VECTOR OF DIMENSION N+1 CONTAINING THE VALUESHEP  140
C                   OF HERMITE POLYNOMIALS OF ORDER 0 UP TO N           HEP  150
C                   FOR GIVEN ARGUMENT X.                               HEP  160
C                   VALUES ARE ORDERED FROM LOW TO HIGH ORDER           HEP  170
C           X     - ARGUMENT OF HERMITE POLYNOMIAL                      HEP  180
C           N     - ORDER OF HERMITE POLYNOMIAL                         HEP  190
C                                                                       HEP  200
C        REMARKS                                                        HEP  210
C           N LESS THAN 0 IS TREATED AS IF N WERE 0                     HEP  220
C                                                                       HEP  230
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  HEP  240
C           NONE                                                        HEP  250
C                                                                       HEP  260
C        METHOD                                                         HEP  270
C           EVALUATION IS BASED ON THE RECURRENCE EQUATION FOR          HEP  280
C           HERMITE POLYNOMIALS H(N,X)                                  HEP  290
C           H(N+1,X)=2*(X*H(N,X)-N*H(N-1,X))                            HEP  300
C           WHERE THE FIRST TERM IN BRACKETS IS THE INDEX,              HEP  310
C           THE SECOND IS THE ARGUMENT.                                 HEP  320
C           STARTING VALUES ARE H(0,X)=1, H(1,X)=2*X.                   HEP  330
C                                                                       HEP  340
C     ..................................................................HEP  350
C                                                                       HEP  360
      SUBROUTINE HEP(Y,X,N)                                             HEP  370
C                                                                       HEP  380
      DIMENSION Y(1)                                                    HEP  390
C                                                                       HEP  400
C        TEST OF ORDER                                                  HEP  410
      Y(1)=1.                                                           HEP  420
      IF(N)1,1,2                                                        HEP  430
    1 RETURN                                                            HEP  440
C                                                                       HEP  450
    2 Y(2)=X+X                                                          HEP  460
      IF(N-1)1,1,3                                                      HEP  470
C                                                                       HEP  480
    3 DO 4 I=2,N                                                        HEP  490
      F=X*Y(I)-FLOAT(I-1)*Y(I-1)                                        HEP  500
    4 Y(I+1)=F+F                                                        HEP  510
      RETURN                                                            HEP  520
      END                                                               HEP  530
