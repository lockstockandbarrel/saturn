C                                                                       CSP   10
C     ..................................................................CSP   20
C                                                                       CSP   30
C        SUBROUTINE CSP                                                 CSP   40
C                                                                       CSP   50
C        PURPOSE                                                        CSP   60
C           COMPUTE THE VALUES OF THE SHIFTED CHEBYSHEV POLYNOMIALS     CSP   70
C           TS(N,X) FOR ARGUMENT X AND ORDERS 0 UP TO N.                CSP   80
C                                                                       CSP   90
C        USAGE                                                          CSP  100
C           CALL CSP(Y,X,N)                                             CSP  110
C                                                                       CSP  120
C        DESCRIPTION OF PARAMETERS                                      CSP  130
C           Y     - RESULT VECTOR OF DIMENSION N+1 CONTAINING THE VALUESCSP  140
C                   OF SHIFTED CHEBYSHEV POLYNOMIALS OF ORDER 0 UP TO N CSP  150
C                   FOR GIVEN ARGUMENT X.                               CSP  160
C                   VALUES ARE ORDERED FROM LOW TO HIGH ORDER           CSP  170
C           X     - ARGUMENT OF SHIFTED CHEBYSHEV POLYNOMIAL            CSP  180
C           N     - ORDER OF SHIFTED CHEBYSHEV POLYNOMIAL               CSP  190
C                                                                       CSP  200
C        REMARKS                                                        CSP  210
C           N LESS THAN 0 IS TREATED AS IF N WERE 0                     CSP  220
C                                                                       CSP  230
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  CSP  240
C           NONE                                                        CSP  250
C                                                                       CSP  260
C        METHOD                                                         CSP  270
C           EVALUATION IS BASED ON THE RECURRENCE EQUATION FOR          CSP  280
C           SHIFTED CHEBYSHEV POLYNOMIALS TS(N,X)                       CSP  290
C           TS(N+1,X)=(4*X-2)*TS(N,X)-TS(N-1,X),                        CSP  300
C           WHERE THE FIRST TERM IN BRACKETS IS THE ORDER,              CSP  310
C           THE SECOND IS THE ARGUMENT.                                 CSP  320
C           STARTING VALUES ARE TS(0,X)=1, TS(1,X)=2*X-1.               CSP  330
C                                                                       CSP  340
C     ..................................................................CSP  350
C                                                                       CSP  360
      SUBROUTINE CSP(Y,X,N)                                             CSP  370
C                                                                       CSP  380
      DIMENSION Y(1)                                                    CSP  390
C                                                                       CSP  400
C        TEST OF ORDER                                                  CSP  410
      Y(1)=1.                                                           CSP  420
      IF(N)1,1,2                                                        CSP  430
    1 RETURN                                                            CSP  440
C                                                                       CSP  450
    2 Y(2)=X+X-1.                                                       CSP  460
      IF(N-1)1,1,3                                                      CSP  470
C                                                                       CSP  480
C        INITIALIZATION                                                 CSP  490
    3 F=Y(2)+Y(2)                                                       CSP  500
C                                                                       CSP  510
      DO 4 I=2,N                                                        CSP  520
    4 Y(I+1)=F*Y(I)-Y(I-1)                                              CSP  530
      RETURN                                                            CSP  540
      END                                                               CSP  550
