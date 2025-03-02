C                                                                       QL8   10
C     ..................................................................QL8   20
C                                                                       QL8   30
C        SUBROUTINE QL8                                                 QL8   40
C                                                                       QL8   50
C        PURPOSE                                                        QL8   60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X), SUMMED OVER X FROM 0    QL8   70
C                               TO INFINITY).                           QL8   80
C                                                                       QL8   90
C        USAGE                                                          QL8  100
C           CALL QL8 (FCT,Y)                                            QL8  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QL8  120
C                                                                       QL8  130
C        DESCRIPTION OF PARAMETERS                                      QL8  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QL8  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QL8  160
C                                                                       QL8  170
C        REMARKS                                                        QL8  180
C           NONE                                                        QL8  190
C                                                                       QL8  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QL8  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QL8  220
C           BY THE USER.                                                QL8  230
C                                                                       QL8  240
C        METHOD                                                         QL8  250
C           EVALUATION IS DONE BY MEANS OF 8-POINT GAUSSIAN-LAGUERRE    QL8  260
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       QL8  270
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 15.                     QL8  280
C           FOR REFERENCE, SEE                                          QL8  290
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QL8  300
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.130-132 AND 347-352.   QL8  310
C                                                                       QL8  320
C     ..................................................................QL8  330
C                                                                       QL8  340
      SUBROUTINE QL8(FCT,Y)                                             QL8  350
C                                                                       QL8  360
C                                                                       QL8  370
      X=22.86313                                                        QL8  380
      Y=.1048001E-8*FCT(X)                                              QL8  390
      X=15.74068                                                        QL8  400
      Y=Y+.8485747E-6*FCT(X)                                            QL8  410
      X=10.75852                                                        QL8  420
      Y=Y+.9076509E-4*FCT(X)                                            QL8  430
      X=7.045905                                                        QL8  440
      Y=Y+.2794536E-2*FCT(X)                                            QL8  450
      X=4.266700                                                        QL8  460
      Y=Y+.03334349*FCT(X)                                              QL8  470
      X=2.251087                                                        QL8  480
      Y=Y+.1757950*FCT(X)                                               QL8  490
      X=.9037018                                                        QL8  500
      Y=Y+.4187868*FCT(X)                                               QL8  510
      X=.1702796                                                        QL8  520
      Y=Y+.3691886*FCT(X)                                               QL8  530
      RETURN                                                            QL8  540
      END                                                               QL8  550
