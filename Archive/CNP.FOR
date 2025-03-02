C                                                                       CNP   10
C     ..................................................................CNP   20
C                                                                       CNP   30
C        SUBROUTINE CNP                                                 CNP   40
C                                                                       CNP   50
C        PURPOSE                                                        CNP   60
C           COMPUTE THE VALUES OF THE CHEBYSHEV POLYNOMIALS T(N,X)      CNP   70
C           FOR ARGUMENT VALUE X AND ORDERS 0 UP TO N.                  CNP   80
C                                                                       CNP   90
C        USAGE                                                          CNP  100
C           CALL CNP(Y,X,N)                                             CNP  110
C                                                                       CNP  120
C        DESCRIPTION OF PARAMETERS                                      CNP  130
C           Y     - RESULT VECTOR OF DIMENSION N+1 CONTAINING THE VALUESCNP  140
C                   OF CHEBYSHEV POLYNOMIALS OF ORDER 0 UP TO N         CNP  150
C                   FOR GIVEN ARGUMENT X.                               CNP  160
C           Y     - RESULT VALUE                                        CNP  170
C                   VALUES ARE ORDERED FROM LOW TO HIGH ORDER           CNP  180
C           X     - ARGUMENT OF CHEBYSHEV POLYNOMIAL                    CNP  190
C           N     - ORDER OF CHEBYSHEV POLYNOMIAL                       CNP  200
C                                                                       CNP  210
C        REMARKS                                                        CNP  220
C           N LESS THAN 0 IS TREATED AS IF N WERE 0                     CNP  230
C                                                                       CNP  240
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  CNP  250
C           NONE                                                        CNP  260
C                                                                       CNP  270
C        METHOD                                                         CNP  280
C           EVALUATION IS BASED ON THE RECURRENCE EQUATION FOR          CNP  290
C           CHEBYSHEV POLYNOMIALS T(N,X)                                CNP  300
C           T(N+1,X)=2*X*T(N,X)-T(N-1,X),                               CNP  310
C           WHERE THE FIRST TERM IN BRACKETS IS THE ORDER,              CNP  320
C           THE SECOND IS THE ARGUMENT.                                 CNP  330
C           STARTING VALUES ARE T(0,X)=1, T(1,X)=X.                     CNP  340
C                                                                       CNP  350
C     ..................................................................CNP  360
C                                                                       CNP  370
      SUBROUTINE CNP(Y,X,N)                                             CNP  380
C                                                                       CNP  390
      DIMENSION Y(1)                                                    CNP  400
      Y(1)=1.                                                           CNP  410
      IF(N)1,1,2                                                        CNP  420
    1 RETURN                                                            CNP  430
C                                                                       CNP  440
    2 Y(2)=X                                                            CNP  450
      IF(N-1)1,1,3                                                      CNP  460
C                                                                       CNP  470
C        INITIALIZATION                                                 CNP  480
    3 F=X+X                                                             CNP  490
C                                                                       CNP  500
      DO 4 I=2,N                                                        CNP  510
    4 Y(I+1)=F*Y(I)-Y(I-1)                                              CNP  520
      RETURN                                                            CNP  530
      END                                                               CNP  540
