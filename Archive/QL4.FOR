C                                                                       QL4   10
C     ..................................................................QL4   20
C                                                                       QL4   30
C        SUBROUTINE QL4                                                 QL4   40
C                                                                       QL4   50
C        PURPOSE                                                        QL4   60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X), SUMMED OVER X FROM 0    QL4   70
C                               TO INFINITY).                           QL4   80
C                                                                       QL4   90
C        USAGE                                                          QL4  100
C           CALL QL4 (FCT,Y)                                            QL4  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QL4  120
C                                                                       QL4  130
C        DESCRIPTION OF PARAMETERS                                      QL4  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QL4  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QL4  160
C                                                                       QL4  170
C        REMARKS                                                        QL4  180
C           NONE                                                        QL4  190
C                                                                       QL4  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QL4  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QL4  220
C           BY THE USER.                                                QL4  230
C                                                                       QL4  240
C        METHOD                                                         QL4  250
C           EVALUATION IS DONE BY MEANS OF 4-POINT GAUSSIAN-LAGUERRE    QL4  260
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       QL4  270
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 7.                      QL4  280
C           FOR REFERENCE, SEE                                          QL4  290
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QL4  300
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.130-132 AND 347-352.   QL4  310
C                                                                       QL4  320
C     ..................................................................QL4  330
C                                                                       QL4  340
      SUBROUTINE QL4(FCT,Y)                                             QL4  350
C                                                                       QL4  360
C                                                                       QL4  370
      X=9.395071                                                        QL4  380
      Y=.5392947E-3*FCT(X)                                              QL4  390
      X=4.536620                                                        QL4  400
      Y=Y+.03888791*FCT(X)                                              QL4  410
      X=1.745761                                                        QL4  420
      Y=Y+.3574187*FCT(X)                                               QL4  430
      X=.3225477                                                        QL4  440
      Y=Y+.6031541*FCT(X)                                               QL4  450
      RETURN                                                            QL4  460
      END                                                               QL4  470
