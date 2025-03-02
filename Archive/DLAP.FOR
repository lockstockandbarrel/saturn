C                                                                       DAP   10
C     ..................................................................DAP   20
C                                                                       DAP   30
C        SUBROUTINE DLAP                                                DAP   40
C                                                                       DAP   50
C        PURPOSE                                                        DAP   60
C           COMPUTE THE VALUES OF THE LAGUERRE POLYNOMIALS L(N,X)       DAP   70
C           FOR ARGUMENT VALUE X AND ORDERS 0 UP TO N.                  DAP   80
C                                                                       DAP   90
C        USAGE                                                          DAP  100
C           CALL LAP(Y,X,N)                                             DAP  110
C                                                                       DAP  120
C        DESCRIPTION OF PARAMETERS                                      DAP  130
C           Y     - RESULT VECTOR OF DIMENSION N+1 CONTAINING THE VALUESDAP  140
C                   OF LAGUERRE POLYNOMIALS OF ORDER 0 UP TO N          DAP  150
C                   FOR GIVEN ARGUMENT X.                               DAP  160
C                   VALUES ARE ORDERED FROM LOW TO HIGH ORDER           DAP  170
C           X     - ARGUMENT OF LAGUERRE POLYNOMIAL                     DAP  180
C           N     - ORDER OF LAGUERRE POLYNOMIAL                        DAP  190
C        NOTE: Y,X AND T ARE DOUBLE PRECISION                           DAP  200
C        REMARKS                                                        DAP  210
C           N LESS THAN 0 IS TREATED AS IF N WERE 0                     DAP  220
C                                                                       DAP  230
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DAP  240
C           NONE                                                        DAP  250
C                                                                       DAP  260
C        METHOD                                                         DAP  270
C           EVALUATION IS BASED ON THE RECURRENCE EQUATION FOR          DAP  280
C           LAGUERRE POLYNOMIALS L(N,X)                                 DAP  290
C           L(N+1,X)=2*L(N,X)-L(N-1,X)-((1+X)*L(N,X)-L(N-1,X))/(N+1),   DAP  300
C           WHERE THE FIRST TERM IN BRACKETS IS THE ORDER,              DAP  310
C           THE SECOND IS THE ARGUMENT.                                 DAP  320
C           STARTING VALUES ARE L(0,X)=1, L(1,X)=1.-X.                  DAP  330
C                                                                       DAP  340
C     ..................................................................DAP  350
C                                                                       DAP  360
      SUBROUTINE DLAP(Y,X,N)                                            DAP  370
      DOUBLE PRECISION Y,X,T                                            DAP  380
      DIMENSION Y(1)                                                    DAP  390
C                                                                       DAP  400
C        TEST OF ORDER                                                  DAP  410
      Y(1)=1.                                                           DAP  420
      IF(N)1,1,2                                                        DAP  430
    1 RETURN                                                            DAP  440
C                                                                       DAP  450
    2 Y(2)=1.-X                                                         DAP  460
      IF(N-1)1,1,3                                                      DAP  470
C                                                                       DAP  480
C        INITIALIZATION                                                 DAP  490
    3 T=1.+X                                                            DAP  500
C                                                                       DAP  510
      DO 4 I=2,N                                                        DAP  520
    4 Y(I+1)=Y(I)-Y(I-1)+Y(I)-(T*Y(I)-Y(I-1))/DFLOAT(I)                 DAP  530
      RETURN                                                            DAP  540
      END                                                               DAP  550
