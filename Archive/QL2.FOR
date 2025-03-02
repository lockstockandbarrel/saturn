C                                                                       QL2   10
C     ..................................................................QL2   20
C                                                                       QL2   30
C        SUBROUTINE QL2                                                 QL2   40
C                                                                       QL2   50
C        PURPOSE                                                        QL2   60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X), SUMMED OVER X FROM 0    QL2   70
C                               TO INFINITY).                           QL2   80
C                                                                       QL2   90
C        USAGE                                                          QL2  100
C           CALL QL2 (FCT,Y)                                            QL2  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QL2  120
C                                                                       QL2  130
C        DESCRIPTION OF PARAMETERS                                      QL2  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QL2  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QL2  160
C                                                                       QL2  170
C        REMARKS                                                        QL2  180
C           NONE                                                        QL2  190
C                                                                       QL2  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QL2  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QL2  220
C           BY THE USER.                                                QL2  230
C                                                                       QL2  240
C        METHOD                                                         QL2  250
C           EVALUATION IS DONE BY MEANS OF 2-POINT GAUSSIAN-LAGUERRE    QL2  260
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       QL2  270
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 3.                      QL2  280
C           FOR REFERENCE, SEE                                          QL2  290
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QL2  300
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.130-132 AND 347-352.   QL2  310
C                                                                       QL2  320
C     ..................................................................QL2  330
C                                                                       QL2  340
      SUBROUTINE QL2(FCT,Y)                                             QL2  350
C                                                                       QL2  360
C                                                                       QL2  370
      X=3.414214                                                        QL2  380
      Y=.1464466*FCT(X)                                                 QL2  390
      X=.5857864                                                        QL2  400
      Y=Y+.8535534*FCT(X)                                               QL2  410
      RETURN                                                            QL2  420
      END                                                               QL2  430
