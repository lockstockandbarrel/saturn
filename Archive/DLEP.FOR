C                                                                       DEP   10
C     ..................................................................DEP   20
C                                                                       DEP   30
C        SUBROUTINE DLEP                                                DEP   40
C                                                                       DEP   50
C        PURPOSE                                                        DEP   60
C           COMPUTE THE VALUES OF THE LEGENDRE POLYNOMIALS P(N,X)       DEP   70
C           FOR ARGUMENT VALUE X AND ORDERS 0 UP TO N.                  DEP   80
C                                                                       DEP   90
C        USAGE                                                          DEP  100
C           CALL DLEP(Y,X,N)                                            DEP  110
C                                                                       DEP  120
C        DESCRIPTION OF PARAMETERS                                      DEP  130
C           Y     - RESULT VECTOR OF DIMENSION N+1 CONTAINING THE VALUESDEP  140
C                   OF LEGENDRE POLYNOMIALS OF ORDER 0 UP TO N          DEP  150
C                   FOR GIVEN ARGUMENT X.                               DEP  160
C                   DOUBLE PRECISION VECTOR.                            DEP  170
C                   VALUES ARE ORDERED FROM LOW TO HIGH ORDER           DEP  180
C           X     - ARGUMENT OF LEGENDRE POLYNOMIAL                     DEP  190
C                   DOUBLE PRECISION VARIABLE.                          DEP  200
C           N     - ORDER OF LEGENDRE POLYNOMIAL                        DEP  210
C                                                                       DEP  220
C        REMARKS                                                        DEP  230
C           N LESS THAN 0 IS TREATED AS IF N WERE 0                     DEP  240
C                                                                       DEP  250
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DEP  260
C           NONE                                                        DEP  270
C                                                                       DEP  280
C        METHOD                                                         DEP  290
C           EVALUATION IS BASED ON THE RECURRENCE EQUATION FOR          DEP  300
C           LEGENDRE POLYNOMIALS P(N,X)                                 DEP  310
C           P(N+1,X)=2*X*P(N,X)-P(N-1,X)-(X*P(N,X)-P(N-1,X))/(N+1),     DEP  320
C           WHERE THE FIRST TERM IN BRACKETS IS THE ORDER,              DEP  330
C           THE SECOND IS THE ARGUMENT.                                 DEP  340
C           STARTING VALUES ARE P(0,X)=1, P(1,X)=X.                     DEP  350
C                                                                       DEP  360
C     ..................................................................DEP  370
C                                                                       DEP  380
      SUBROUTINE DLEP(Y,X,N)                                            DEP  390
C                                                                       DEP  400
      DIMENSION Y(1)                                                    DEP  410
      DOUBLE PRECISION Y,X,G                                            DEP  420
C                                                                       DEP  430
C        TEST OF ORDER                                                  DEP  440
      Y(1)=1.D0                                                         DEP  450
      IF(N)1,1,2                                                        DEP  460
    1 RETURN                                                            DEP  470
C                                                                       DEP  480
    2 Y(2)=X                                                            DEP  490
      IF(N-1)1,1,3                                                      DEP  500
C                                                                       DEP  510
    3 DO 4 I=2,N                                                        DEP  520
      G=X*Y(I)                                                          DEP  530
    4 Y(I+1)=G-Y(I-1)+G-(G-Y(I-1))/DFLOAT(I)                            DEP  540
      RETURN                                                            DEP  550
      END                                                               DEP  560
