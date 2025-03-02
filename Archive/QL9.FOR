C                                                                       QL9   10
C     ..................................................................QL9   20
C                                                                       QL9   30
C        SUBROUTINE QL9                                                 QL9   40
C                                                                       QL9   50
C        PURPOSE                                                        QL9   60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X), SUMMED OVER X FROM 0    QL9   70
C                               TO INFINITY).                           QL9   80
C                                                                       QL9   90
C        USAGE                                                          QL9  100
C           CALL QL9 (FCT,Y)                                            QL9  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QL9  120
C                                                                       QL9  130
C        DESCRIPTION OF PARAMETERS                                      QL9  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QL9  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QL9  160
C                                                                       QL9  170
C        REMARKS                                                        QL9  180
C           NONE                                                        QL9  190
C                                                                       QL9  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QL9  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QL9  220
C           BY THE USER.                                                QL9  230
C                                                                       QL9  240
C        METHOD                                                         QL9  250
C           EVALUATION IS DONE BY MEANS OF 9-POINT GAUSSIAN-LAGUERRE    QL9  260
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       QL9  270
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 17.                     QL9  280
C           FOR REFERENCE, SEE                                          QL9  290
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QL9  300
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.130-132 AND 347-352.   QL9  310
C                                                                       QL9  320
C     ..................................................................QL9  330
C                                                                       QL9  340
      SUBROUTINE QL9(FCT,Y)                                             QL9  350
C                                                                       QL9  360
C                                                                       QL9  370
      X=26.37407                                                        QL9  380
      Y=.3290874E-10*FCT(X)                                             QL9  390
      X=18.83360                                                        QL9  400
      Y=Y+.4110769E-7*FCT(X)                                            QL9  410
      X=13.46624                                                        QL9  420
      Y=Y+.6592123E-5*FCT(X)                                            QL9  430
      X=9.372985                                                        QL9  440
      Y=Y+.3052498E-3*FCT(X)                                            QL9  450
      X=6.204957                                                        QL9  460
      Y=Y+.005599627*FCT(X)                                             QL9  470
      X=3.783474                                                        QL9  480
      Y=Y+.04746056*FCT(X)                                              QL9  490
      X=2.005135                                                        QL9  500
      Y=Y+.1992875*FCT(X)                                               QL9  510
      X=.8072200                                                        QL9  520
      Y=Y+.4112140*FCT(X)                                               QL9  530
      X=.1523222                                                        QL9  540
      Y=Y+.3361264*FCT(X)                                               QL9  550
      RETURN                                                            QL9  560
      END                                                               QL9  570
