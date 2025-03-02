C                                                                       QL6   10
C     ..................................................................QL6   20
C                                                                       QL6   30
C        SUBROUTINE QL6                                                 QL6   40
C                                                                       QL6   50
C        PURPOSE                                                        QL6   60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X), SUMMED OVER X FROM 0    QL6   70
C                               TO INFINITY).                           QL6   80
C                                                                       QL6   90
C        USAGE                                                          QL6  100
C           CALL QL6 (FCT,Y)                                            QL6  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QL6  120
C                                                                       QL6  130
C        DESCRIPTION OF PARAMETERS                                      QL6  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QL6  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QL6  160
C                                                                       QL6  170
C        REMARKS                                                        QL6  180
C           NONE                                                        QL6  190
C                                                                       QL6  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QL6  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QL6  220
C           BY THE USER.                                                QL6  230
C                                                                       QL6  240
C        METHOD                                                         QL6  250
C           EVALUATION IS DONE BY MEANS OF 6-POINT GAUSSIAN-LAGUERRE    QL6  260
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       QL6  270
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 11.                     QL6  280
C           FOR REFERENCE, SEE                                          QL6  290
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QL6  300
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.130-132 AND 347-352.   QL6  310
C                                                                       QL6  320
C     ..................................................................QL6  330
C                                                                       QL6  340
      SUBROUTINE QL6(FCT,Y)                                             QL6  350
C                                                                       QL6  360
C                                                                       QL6  370
      X=15.98287                                                        QL6  380
      Y=.8985479E-6*FCT(X)                                              QL6  390
      X=9.837467                                                        QL6  400
      Y=Y+.2610172E-3*FCT(X)                                            QL6  410
      X=5.775144                                                        QL6  420
      Y=Y+.01039920*FCT(X)                                              QL6  430
      X=2.992736                                                        QL6  440
      Y=Y+.1133734*FCT(X)                                               QL6  450
      X=1.188932                                                        QL6  460
      Y=Y+.4170008*FCT(X)                                               QL6  470
      X=.2228466                                                        QL6  480
      Y=Y+.4589647*FCT(X)                                               QL6  490
      RETURN                                                            QL6  500
      END                                                               QL6  510
