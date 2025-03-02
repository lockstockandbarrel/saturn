C                                                                       LEP   10
C     ..................................................................LEP   20
C                                                                       LEP   30
C        SUBROUTINE LEP                                                 LEP   40
C                                                                       LEP   50
C        PURPOSE                                                        LEP   60
C           COMPUTE THE VALUES OF THE LEGENDRE POLYNOMIALS P(N,X)       LEP   70
C           FOR ARGUMENT VALUE X AND ORDERS 0 UP TO N.                  LEP   80
C                                                                       LEP   90
C        USAGE                                                          LEP  100
C           CALL LEP(Y,X,N)                                             LEP  110
C                                                                       LEP  120
C        DESCRIPTION OF PARAMETERS                                      LEP  130
C           Y     - RESULT VECTOR OF DIMENSION N+1 CONTAINING THE VALUESLEP  140
C                   OF LEGENDRE POLYNOMIALS OF ORDER 0 UP TO N          LEP  150
C                   FOR GIVEN ARGUMENT X.                               LEP  160
C                   VALUES ARE ORDERED FROM LOW TO HIGH ORDER           LEP  170
C           X     - ARGUMENT OF LEGENDRE POLYNOMIAL                     LEP  180
C           N     - ORDER OF LEGENDRE POLYNOMIAL                        LEP  190
C                                                                       LEP  200
C        REMARKS                                                        LEP  210
C           N LESS THAN 0 IS TREATED AS IF N WERE 0                     LEP  220
C                                                                       LEP  230
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  LEP  240
C           NONE                                                        LEP  250
C                                                                       LEP  260
C        METHOD                                                         LEP  270
C           EVALUATION IS BASED ON THE RECURRENCE EQUATION FOR          LEP  280
C           LEGENDRE POLYNOMIALS P(N,X)                                 LEP  290
C           P(N+1,X)=2*X*P(N,X)-P(N-1,X)-(X*P(N,X)-P(N-1,X))/(N+1),     LEP  300
C           WHERE THE FIRST TERM IN BRACKETS IS THE ORDER,              LEP  310
C           THE SECOND IS THE ARGUMENT.                                 LEP  320
C           STARTING VALUES ARE P(0,X)=1, P(1,X)=X.                     LEP  330
C                                                                       LEP  340
C     ..................................................................LEP  350
C                                                                       LEP  360
      SUBROUTINE LEP(Y,X,N)                                             LEP  370
C                                                                       LEP  380
      DIMENSION Y(1)                                                    LEP  390
C                                                                       LEP  400
C        TEST OF ORDER                                                  LEP  410
      Y(1)=1.                                                           LEP  420
      IF(N)1,1,2                                                        LEP  430
    1 RETURN                                                            LEP  440
C                                                                       LEP  450
    2 Y(2)=X                                                            LEP  460
      IF(N-1)1,1,3                                                      LEP  470
C                                                                       LEP  480
    3 DO 4 I=2,N                                                        LEP  490
      G=X*Y(I)                                                          LEP  500
    4 Y(I+1)=G-Y(I-1)+G-(G-Y(I-1))/FLOAT(I)                             LEP  510
      RETURN                                                            LEP  520
      END                                                               LEP  530
