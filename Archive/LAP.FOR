C                                                                       LAP   10
C     ..................................................................LAP   20
C                                                                       LAP   30
C        SUBROUTINE LAP                                                 LAP   40
C                                                                       LAP   50
C        PURPOSE                                                        LAP   60
C           COMPUTE THE VALUES OF THE LAGUERRE POLYNOMIALS L(N,X)       LAP   70
C           FOR ARGUMENT VALUE X AND ORDERS 0 UP TO N.                  LAP   80
C                                                                       LAP   90
C        USAGE                                                          LAP  100
C           CALL LAP(Y,X,N)                                             LAP  110
C                                                                       LAP  120
C        DESCRIPTION OF PARAMETERS                                      LAP  130
C           Y     - RESULT VECTOR OF DIMENSION N+1 CONTAINING THE VALUESLAP  140
C                   OF LAGUERRE POLYNOMIALS OF ORDER 0 UP TO N          LAP  150
C                   FOR GIVEN ARGUMENT X.                               LAP  160
C                   VALUES ARE ORDERED FROM LOW TO HIGH ORDER           LAP  170
C           X     - ARGUMENT OF LAGUERRE POLYNOMIAL                     LAP  180
C           N     - ORDER OF LAGUERRE POLYNOMIAL                        LAP  190
C                                                                       LAP  200
C        REMARKS                                                        LAP  210
C           N LESS THAN 0 IS TREATED AS IF N WERE 0                     LAP  220
C                                                                       LAP  230
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  LAP  240
C           NONE                                                        LAP  250
C                                                                       LAP  260
C        METHOD                                                         LAP  270
C           EVALUATION IS BASED ON THE RECURRENCE EQUATION FOR          LAP  280
C           LAGUERRE POLYNOMIALS L(N,X)                                 LAP  290
C           L(N+1,X)=2*L(N,X)-L(N-1,X)-((1+X)*L(N,X)-L(N-1,X))/(N+1),   LAP  300
C           WHERE THE FIRST TERM IN BRACKETS IS THE ORDER,              LAP  310
C           THE SECOND IS THE ARGUMENT.                                 LAP  320
C           STARTING VALUES ARE L(0,X)=1, L(1,X)=1.-X.                  LAP  330
C                                                                       LAP  340
C     ..................................................................LAP  350
C                                                                       LAP  360
      SUBROUTINE LAP(Y,X,N)                                             LAP  370
C                                                                       LAP  380
      DIMENSION Y(1)                                                    LAP  390
C                                                                       LAP  400
C        TEST OF ORDER                                                  LAP  410
      Y(1)=1.                                                           LAP  420
      IF(N)1,1,2                                                        LAP  430
    1 RETURN                                                            LAP  440
C                                                                       LAP  450
    2 Y(2)=1.-X                                                         LAP  460
      IF(N-1)1,1,3                                                      LAP  470
C                                                                       LAP  480
C        INITIALIZATION                                                 LAP  490
    3 T=1.+X                                                            LAP  500
C                                                                       LAP  510
      DO 4 I=2,N                                                        LAP  520
    4 Y(I+1)=Y(I)-Y(I-1)+Y(I)-(T*Y(I)-Y(I-1))/FLOAT(I)                  LAP  530
      RETURN                                                            LAP  540
      END                                                               LAP  550
