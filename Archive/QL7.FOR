C                                                                       QL7   10
C     ..................................................................QL7   20
C                                                                       QL7   30
C        SUBROUTINE QL7                                                 QL7   40
C                                                                       QL7   50
C        PURPOSE                                                        QL7   60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X), SUMMED OVER X FROM 0    QL7   70
C                               TO INFINITY).                           QL7   80
C                                                                       QL7   90
C        USAGE                                                          QL7  100
C           CALL QL7 (FCT,Y)                                            QL7  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QL7  120
C                                                                       QL7  130
C        DESCRIPTION OF PARAMETERS                                      QL7  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QL7  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QL7  160
C                                                                       QL7  170
C        REMARKS                                                        QL7  180
C           NONE                                                        QL7  190
C                                                                       QL7  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QL7  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QL7  220
C           BY THE USER.                                                QL7  230
C                                                                       QL7  240
C        METHOD                                                         QL7  250
C           EVALUATION IS DONE BY MEANS OF 7-POINT GAUSSIAN-LAGUERRE    QL7  260
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       QL7  270
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 13.                     QL7  280
C           FOR REFERENCE, SEE                                          QL7  290
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QL7  300
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.130-132 AND 347-352.   QL7  310
C                                                                       QL7  320
C     ..................................................................QL7  330
C                                                                       QL7  340
      SUBROUTINE QL7(FCT,Y)                                             QL7  350
C                                                                       QL7  360
C                                                                       QL7  370
      X=19.39573                                                        QL7  380
      Y=.3170315E-7*FCT(X)                                              QL7  390
      X=12.73418                                                        QL7  400
      Y=Y+.1586546E-4*FCT(X)                                            QL7  410
      X=8.182153                                                        QL7  420
      Y=Y+.1074010E-2*FCT(X)                                            QL7  430
      X=4.900353                                                        QL7  440
      Y=Y+.02063351*FCT(X)                                              QL7  450
      X=2.567877                                                        QL7  460
      Y=Y+.1471263*FCT(X)                                               QL7  470
      X=1.026665                                                        QL7  480
      Y=Y+.4218313*FCT(X)                                               QL7  490
      X=.1930437                                                        QL7  500
      Y=Y+.4093190*FCT(X)                                               QL7  510
      RETURN                                                            QL7  520
      END                                                               QL7  530
