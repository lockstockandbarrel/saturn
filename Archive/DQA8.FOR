C                                                                       DQA8  10
C     ..................................................................DQA8  20
C                                                                       DQA8  30
C        SUBROUTINE DQA8                                                DQA8  40
C                                                                       DQA8  50
C        PURPOSE                                                        DQA8  60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X)/SQRT(X), SUMMED OVER X   DQA8  70
C                               FROM 0 TO INFINITY).                    DQA8  80
C                                                                       DQA8  90
C        USAGE                                                          DQA8 100
C           CALL DQA8 (FCT,Y)                                           DQA8 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DQA8 120
C                                                                       DQA8 130
C        DESCRIPTION OF PARAMETERS                                      DQA8 140
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DQA8 150
C                    SUBPROGRAM USED.                                   DQA8 160
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DQA8 170
C                                                                       DQA8 180
C        REMARKS                                                        DQA8 190
C           NONE                                                        DQA8 200
C                                                                       DQA8 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DQA8 220
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DQA8 230
C           MUST BE FURNISHED BY THE USER.                              DQA8 240
C                                                                       DQA8 250
C        METHOD                                                         DQA8 260
C           EVALUATION IS DONE BY MEANS OF 8-POINT GENERALIZED GAUSS-   DQA8 270
C           LAGUERRE QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY       DQA8 280
C           WHENEVER FCT(X) IS A POLYNOMIAL UP TO DEGREE 15.            DQA8 290
C           FOR REFERENCE, SEE                                          DQA8 300
C           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF    DQA8 310
C           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED     DQA8 320
C           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT       DQA8 330
C           TR00.1100 (MARCH 1964), PP.15-16.                           DQA8 340
C                                                                       DQA8 350
C     ..................................................................DQA8 360
C                                                                       DQA8 370
      SUBROUTINE DQA8(FCT,Y)                                            DQA8 380
C                                                                       DQA8 390
C                                                                       DQA8 400
      DOUBLE PRECISION X,Y,FCT                                          DQA8 410
C                                                                       DQA8 420
      X=.21984272840962651D2                                            DQA8 430
      Y=.53096149480223645D-9*FCT(X)                                    DQA8 440
      X=.14972627088426393D2                                            DQA8 450
      Y=Y+.46419616897304213D-6*FCT(X)                                  DQA8 460
      X=.10093323675221343D2                                            DQA8 470
      Y=Y+.54237201850757630D-4*FCT(X)                                  DQA8 480
      X=.64831454286271704D1                                            DQA8 490
      Y=Y+.18645680172483611D-2*FCT(X)                                  DQA8 500
      X=.38094763614849071D1                                            DQA8 510
      Y=Y+.25760623071019947D-1*FCT(X)                                  DQA8 520
      X=.19051136350314284D1                                            DQA8 530
      Y=Y+.16762008279797166D0*FCT(X)                                   DQA8 540
      X=.67724908764928915D0                                            DQA8 550
      Y=Y+.56129491705706735D0*FCT(X)                                   DQA8 560
      X=.7479188259681827D-1                                            DQA8 570
      Y=Y+.10158589580332275D1*FCT(X)                                   DQA8 580
      RETURN                                                            DQA8 590
      END                                                               DQA8 600
