C                                                                       DSP   10
C     ..................................................................DSP   20
C                                                                       DSP   30
C        SUBROUTINE DCSP                                                DSP   40
C                                                                       DSP   50
C        PURPOSE                                                        DSP   60
C           COMPUTE THE VALUES OF THE SHIFTED CHEBYSHEV POLYNOMIALS     DSP   70
C           TS(N,X) FOR ARGUMENT X AND ORDERS 0 UP TO N.                DSP   80
C                                                                       DSP   90
C        USAGE                                                          DSP  100
C           CALL DCSP(Y,X,N)                                            DSP  110
C                                                                       DSP  120
C        DESCRIPTION OF PARAMETERS                                      DSP  130
C           Y     - RESULT VECTOR OF DIMENSION N+1 CONTAINING THE VALUESDSP  140
C                   OF SHIFTED CHEBYSHEV POLYNOMIALS OF ORDER 0 UP TO N DSP  150
C                   FOR GIVEN ARGUMENT X.                               DSP  160
C                   DOUBLE PRECISION VECTOR.                            DSP  170
C                   VALUES ARE ORDERED FROM LOW TO HIGH ORDER           DSP  180
C           X     - ARGUMENT OF SHIFTED CHEBYSHEV POLYNOMIAL            DSP  190
C                   DOUBLE PRECISION VARIABLE.                          DSP  200
C           N     - ORDER OF SHIFTED CHEBYSHEV POLYNOMIAL               DSP  210
C                                                                       DSP  220
C        REMARKS                                                        DSP  230
C           N LESS THAN 0 IS TREATED AS IF N WERE 0                     DSP  240
C                                                                       DSP  250
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DSP  260
C           NONE                                                        DSP  270
C                                                                       DSP  280
C        METHOD                                                         DSP  290
C           EVALUATION IS BASED ON THE RECURRENCE EQUATION FOR          DSP  300
C           SHIFTED CHEBYSHEV POLYNOMIALS TS(N,X)                       DSP  310
C           TS(N+1,X)=(4*X-2)*TS(N,X)-TS(N-1,X),                        DSP  320
C           WHERE THE FIRST TERM IN BRACKETS IS THE ORDER,              DSP  330
C           THE SECOND IS THE ARGUMENT.                                 DSP  340
C           STARTING VALUES ARE TS(0,X)=1, TS(1,X)=2*X-1.               DSP  350
C                                                                       DSP  360
C     ..................................................................DSP  370
C                                                                       DSP  380
      SUBROUTINE DCSP(Y,X,N)                                            DSP  390
C                                                                       DSP  400
      DIMENSION Y(1)                                                    DSP  410
      DOUBLE PRECISION Y,X,F                                            DSP  420
C                                                                       DSP  430
C        TEST OF ORDER                                                  DSP  440
      Y(1)=1.D0                                                         DSP  450
      IF(N)1,1,2                                                        DSP  460
    1 RETURN                                                            DSP  470
C                                                                       DSP  480
    2 Y(2)=X+X-1.D0                                                     DSP  490
      IF(N-1)1,1,3                                                      DSP  500
C                                                                       DSP  510
C        INITIALIZATION                                                 DSP  520
    3 F=Y(2)+Y(2)                                                       DSP  530
C                                                                       DSP  540
      DO 4 I=2,N                                                        DSP  550
    4 Y(I+1)=F*Y(I)-Y(I-1)                                              DSP  560
      RETURN                                                            DSP  570
      END                                                               DSP  580
