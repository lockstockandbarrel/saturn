C                                                                       DQH8  10
C     ..................................................................DQH8  20
C                                                                       DQH8  30
C        SUBROUTINE DQH8                                                DQH8  40
C                                                                       DQH8  50
C        PURPOSE                                                        DQH8  60
C           TO COMPUTE INTEGRAL(EXP(-X*X)*FCT(X), SUMMED OVER X FROM    DQH8  70
C                               -INFINITY TO +INFINITY).                DQH8  80
C                                                                       DQH8  90
C        USAGE                                                          DQH8 100
C           CALL DQH8 (FCT,Y)                                           DQH8 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DQH8 120
C                                                                       DQH8 130
C        DESCRIPTION OF PARAMETERS                                      DQH8 140
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DQH8 150
C                    SUBPROGRAM USED.                                   DQH8 160
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DQH8 170
C                                                                       DQH8 180
C        REMARKS                                                        DQH8 190
C           NONE                                                        DQH8 200
C                                                                       DQH8 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DQH8 220
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DQH8 230
C           MUST BE FURNISHED BY THE USER.                              DQH8 240
C                                                                       DQH8 250
C        METHOD                                                         DQH8 260
C           EVALUATION IS DONE BY MEANS OF 8-POINT GAUSSIAN-HERMITE     DQH8 270
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       DQH8 280
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 15.                     DQH8 290
C           FOR REFERENCE, SEE                                          DQH8 300
C           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF    DQH8 310
C           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED     DQH8 320
C           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT       DQH8 330
C           TR00.1100 (MARCH 1964), PP.213-214.                         DQH8 340
C                                                                       DQH8 350
C     ..................................................................DQH8 360
C                                                                       DQH8 370
      SUBROUTINE DQH8(FCT,Y)                                            DQH8 380
C                                                                       DQH8 390
C                                                                       DQH8 400
      DOUBLE PRECISION X,Y,Z,FCT                                        DQH8 410
C                                                                       DQH8 420
      X=.29306374202572440D1                                            DQH8 430
      Z=-X                                                              DQH8 440
      Y=.19960407221136762D-3*(FCT(X)+FCT(Z))                           DQH8 450
      X=.19816567566958429D1                                            DQH8 460
      Z=-X                                                              DQH8 470
      Y=Y+.17077983007413475D-1*(FCT(X)+FCT(Z))                         DQH8 480
      X=.11571937124467802D1                                            DQH8 490
      Z=-X                                                              DQH8 500
      Y=Y+.20780232581489188D0*(FCT(X)+FCT(Z))                          DQH8 510
      X=.38118699020732212D0                                            DQH8 520
      Z=-X                                                              DQH8 530
      Y=Y+.66114701255824129D0*(FCT(X)+FCT(Z))                          DQH8 540
      RETURN                                                            DQH8 550
      END                                                               DQH8 560
