C                                                                       QL3   10
C     ..................................................................QL3   20
C                                                                       QL3   30
C        SUBROUTINE QL3                                                 QL3   40
C                                                                       QL3   50
C        PURPOSE                                                        QL3   60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X), SUMMED OVER X FROM 0    QL3   70
C                               TO INFINITY).                           QL3   80
C                                                                       QL3   90
C        USAGE                                                          QL3  100
C           CALL QL3 (FCT,Y)                                            QL3  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QL3  120
C                                                                       QL3  130
C        DESCRIPTION OF PARAMETERS                                      QL3  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QL3  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QL3  160
C                                                                       QL3  170
C        REMARKS                                                        QL3  180
C           NONE                                                        QL3  190
C                                                                       QL3  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QL3  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QL3  220
C           BY THE USER.                                                QL3  230
C                                                                       QL3  240
C        METHOD                                                         QL3  250
C           EVALUATION IS DONE BY MEANS OF 3-POINT GAUSSIAN-LAGUERRE    QL3  260
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       QL3  270
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 5.                      QL3  280
C           FOR REFERENCE, SEE                                          QL3  290
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QL3  300
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.130-132 AND 347-352.   QL3  310
C                                                                       QL3  320
C     ..................................................................QL3  330
C                                                                       QL3  340
      SUBROUTINE QL3(FCT,Y)                                             QL3  350
C                                                                       QL3  360
C                                                                       QL3  370
      X=6.289945                                                        QL3  380
      Y=.01038926*FCT(X)                                                QL3  390
      X=2.294280                                                        QL3  400
      Y=Y+.2785177*FCT(X)                                               QL3  410
      X=.4157746                                                        QL3  420
      Y=Y+.7110930*FCT(X)                                               QL3  430
      RETURN                                                            QL3  440
      END                                                               QL3  450
