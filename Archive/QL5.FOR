C                                                                       QL5   10
C     ..................................................................QL5   20
C                                                                       QL5   30
C        SUBROUTINE QL5                                                 QL5   40
C                                                                       QL5   50
C        PURPOSE                                                        QL5   60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X), SUMMED OVER X FROM 0    QL5   70
C                               TO INFINITY).                           QL5   80
C                                                                       QL5   90
C        USAGE                                                          QL5  100
C           CALL QL5 (FCT,Y)                                            QL5  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QL5  120
C                                                                       QL5  130
C        DESCRIPTION OF PARAMETERS                                      QL5  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QL5  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QL5  160
C                                                                       QL5  170
C        REMARKS                                                        QL5  180
C           NONE                                                        QL5  190
C                                                                       QL5  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QL5  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QL5  220
C           BY THE USER.                                                QL5  230
C                                                                       QL5  240
C        METHOD                                                         QL5  250
C           EVALUATION IS DONE BY MEANS OF 5-POINT GAUSSIAN-LAGUERRE    QL5  260
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       QL5  270
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 9.                      QL5  280
C           FOR REFERENCE, SEE                                          QL5  290
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QL5  300
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.130-132 AND 347-352.   QL5  310
C                                                                       QL5  320
C     ..................................................................QL5  330
C                                                                       QL5  340
      SUBROUTINE QL5(FCT,Y)                                             QL5  350
C                                                                       QL5  360
C                                                                       QL5  370
      X=12.64080                                                        QL5  380
      Y=.2336997E-4*FCT(X)                                              QL5  390
      X=7.085810                                                        QL5  400
      Y=Y+.3611759E-2*FCT(X)                                            QL5  410
      X=3.596426                                                        QL5  420
      Y=Y+.07594245*FCT(X)                                              QL5  430
      X=1.413403                                                        QL5  440
      Y=Y+.3986668*FCT(X)                                               QL5  450
      X=.2635603                                                        QL5  460
      Y=Y+.5217556*FCT(X)                                               QL5  470
      RETURN                                                            QL5  480
      END                                                               QL5  490
