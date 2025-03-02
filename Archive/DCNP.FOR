C                                                                       DNP   10
C     ..................................................................DNP   20
C                                                                       DNP   30
C        SUBROUTINE DCNP                                                DNP   40
C                                                                       DNP   50
C        PURPOSE                                                        DNP   60
C           COMPUTE THE VALUES OF THE CHEBYSHEV POLYNOMIALS T(N,X)      DNP   70
C           FOR ARGUMENT VALUE X AND ORDERS 0 UP TO N.                  DNP   80
C                                                                       DNP   90
C        USAGE                                                          DNP  100
C           CALL DCNP,Y,X,N)                                            DNP  110
C                                                                       DNP  120
C        DESCRIPTION OF PARAMETERS                                      DNP  130
C           Y     - RESULT VECTOR OF DIMENSION N+1 CONTAINING THE VALUESDNP  140
C                   OF CHEBYSHEV POLYNOMIALS OF ORDER 0 UP TO N         DNP  150
C                   FOR GIVEN ARGUMENT X.                               DNP  160
C                   DOUBLE PRECISION VECTOR.                            DNP  170
C                   VALUES ARE ORDERED FROM LOW TO HIGH ORDER           DNP  180
C           Y     - RESULT VALUE                                        DNP  190
C                   DOUBLE PRECISION VARIABLE.                          DNP  200
C           X     - ARGUMENT OF CHEBYSHEV POLYNOMIAL                    DNP  210
C           N     - ORDER OF CHEBYSHEV POLYNOMIAL                       DNP  220
C                                                                       DNP  230
C        REMARKS                                                        DNP  240
C           N LESS THAN 0 IS TREATED AS IF N WERE 0                     DNP  250
C                                                                       DNP  260
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DNP  270
C           NONE                                                        DNP  280
C                                                                       DNP  290
C        METHOD                                                         DNP  300
C           EVALUATION IS BASED ON THE RECURRENCE EQUATION FOR          DNP  310
C           CHEBYSHEV POLYNOMIALS T(N,X)                                DNP  320
C           T(N+1,X)=2*X*T(N,X)-T(N-1,X),                               DNP  330
C           WHERE THE FIRST TERM IN BRACKETS IS THE ORDER,              DNP  340
C           THE SECOND IS THE ARGUMENT.                                 DNP  350
C           STARTING VALUES ARE T(0,X)=1, T(1,X)=X.                     DNP  360
C                                                                       DNP  370
C     ..................................................................DNP  380
C                                                                       DNP  390
      SUBROUTINE DCNP(Y,X,N)                                            DNP  400
C                                                                       DNP  410
      DIMENSION Y(1)                                                    DNP  420
      DOUBLE PRECISION Y,X,F                                            DNP  430
C                                                                       DNP  440
      Y(1)=1.D0                                                         DNP  450
      IF(N)1,1,2                                                        DNP  460
    1 RETURN                                                            DNP  470
C                                                                       DNP  480
    2 Y(2)=X                                                            DNP  490
      IF(N-1)1,1,3                                                      DNP  500
C                                                                       DNP  510
C        INITIALIZATION                                                 DNP  520
    3 F=X+X                                                             DNP  530
C                                                                       DNP  540
      DO 4 I=2,N                                                        DNP  550
    4 Y(I+1)=F*Y(I)-Y(I-1)                                              DNP  560
      RETURN                                                            DNP  570
      END                                                               DNP  580
