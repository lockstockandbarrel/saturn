C                                                                       DHP   10
C     ..................................................................DHP   20
C                                                                       DHP   30
C        SUBROUTINE DHEP                                                DHP   40
C                                                                       DHP   50
C        PURPOSE                                                        DHP   60
C           COMPUTE THE VALUES OF THE HERMITE POLYNOMIALS H(N,X)        DHP   70
C           FOR ARGUMENT VALUE X AND ORDERS 0 UP TO N.                  DHP   80
C                                                                       DHP   90
C        USAGE                                                          DHP  100
C           CALL DHEP(Y,X,N)                                            DHP  110
C                                                                       DHP  120
C        DESCRIPTION OF PARAMETERS                                      DHP  130
C           Y     - RESULT VECTOR OF DIMENSION N+1 CONTAINING THE VALUESDHP  140
C                   OF HERMITE POLYNOMIALS OF ORDER 0 UP TO N           DHP  150
C                   FOR GIVEN ARGUMENT X.                               DHP  160
C                   DOUBLE PRECISION VECTOR.                            DHP  170
C                   VALUES ARE ORDERED FROM LOW TO HIGH ORDER           DHP  180
C           X     - ARGUMENT OF HERMITE POLYNOMIAL                      DHP  190
C                   DOUBLE PRECISION VARIABLE.                          DHP  200
C           N     - ORDER OF HERMITE POLYNOMIAL                         DHP  210
C                                                                       DHP  220
C        REMARKS                                                        DHP  230
C           N LESS THAN 0 IS TREATED AS IF N WERE 0                     DHP  240
C                                                                       DHP  250
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DHP  260
C           NONE                                                        DHP  270
C                                                                       DHP  280
C        METHOD                                                         DHP  290
C           EVALUATION IS BASED ON THE RECURRENCE EQUATION FOR          DHP  300
C           HERMITE POLYNOMIALS H(N,X)                                  DHP  310
C           H(N+1,X)=2*(X*H(N,X)-N*H(N-1,X))                            DHP  320
C           WHERE THE FIRST TERM IN BRACKETS IS THE INDEX,              DHP  330
C           THE SECOND IS THE ARGUMENT.                                 DHP  340
C           STARTING VALUES ARE H(0,X)=1, H(1,X)=2*X.                   DHP  350
C                                                                       DHP  360
C     ..................................................................DHP  370
C                                                                       DHP  380
      SUBROUTINE DHEP(Y,X,N)                                            DHP  390
C                                                                       DHP  400
      DIMENSION Y(1)                                                    DHP  410
      DOUBLE PRECISION Y,X,F                                            DHP  420
C                                                                       DHP  430
C        TEST OF ORDER                                                  DHP  440
      Y(1)=1.D0                                                         DHP  450
      IF(N)1,1,2                                                        DHP  460
    1 RETURN                                                            DHP  470
C                                                                       DHP  480
    2 Y(2)=X+X                                                          DHP  490
      IF(N-1)1,1,3                                                      DHP  500
C                                                                       DHP  510
    3 DO 4 I=2,N                                                        DHP  520
      F=X*Y(I)-DFLOAT(I-1)*Y(I-1)                                       DHP  530
    4 Y(I+1)=F+F                                                        DHP  540
      RETURN                                                            DHP  550
      END                                                               DHP  560
