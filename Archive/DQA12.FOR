C                                                                       DA12  10
C     ..................................................................DA12  20
C                                                                       DA12  30
C        SUBROUTINE DQA12                                               DA12  40
C                                                                       DA12  50
C        PURPOSE                                                        DA12  60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X)/SQRT(X), SUMMED OVER X   DA12  70
C                               FROM 0 TO INFINITY).                    DA12  80
C                                                                       DA12  90
C        USAGE                                                          DA12 100
C           CALL DQA12 (FCT,Y)                                          DA12 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DA12 120
C                                                                       DA12 130
C        DESCRIPTION OF PARAMETERS                                      DA12 140
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DA12 150
C                    SUBPROGRAM USED.                                   DA12 160
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DA12 170
C                                                                       DA12 180
C        REMARKS                                                        DA12 190
C           NONE                                                        DA12 200
C                                                                       DA12 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DA12 220
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DA12 230
C           MUST BE FURNISHED BY THE USER.                              DA12 240
C                                                                       DA12 250
C        METHOD                                                         DA12 260
C           EVALUATION IS DONE BY MEANS OF 12-POINT GENERALIZED GAUSS-  DA12 270
C           LAGUERRE QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY       DA12 280
C           WHENEVER FCT(X) IS A POLYNOMIAL UP TO DEGREE 23.            DA12 290
C           FOR REFERENCE, SEE                                          DA12 300
C           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF    DA12 310
C           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED     DA12 320
C           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT       DA12 330
C           TR00.1100 (MARCH 1964), PP.15-16.                           DA12 340
C                                                                       DA12 350
C     ..................................................................DA12 360
C                                                                       DA12 370
      SUBROUTINE DQA12(FCT,Y)                                           DA12 380
C                                                                       DA12 390
C                                                                       DA12 400
      DOUBLE PRECISION X,Y,FCT                                          DA12 410
C                                                                       DA12 420
      X=.36191360360615602D2                                            DA12 430
      Y=.33287369929782177D-15*FCT(X)                                   DA12 440
      X=.27661108779846090D2                                            DA12 450
      Y=Y+.13169240486156340D-11*FCT(X)                                 DA12 460
      X=.21396755936166109D2                                            DA12 470
      Y=Y+.60925085399751278D-9*FCT(X)                                  DA12 480
      X=.16432195087675313D2                                            DA12 490
      Y=Y+.8037942349882859D-7*FCT(X)                                   DA12 500
      X=.12390447963809471D2                                            DA12 510
      Y=Y+.43164914098046673D-5*FCT(X)                                  DA12 520
      X=.9075434230961203D1                                             DA12 530
      Y=Y+.11377383272808760D-3*FCT(X)                                  DA12 540
      X=.63699753880306349D1                                            DA12 550
      Y=Y+.16473849653768349D-2*FCT(X)                                  DA12 560
      X=.41984156448784132D1                                            DA12 570
      Y=Y+.14096711620145342D-1*FCT(X)                                  DA12 580
      X=.25098480972321280D1                                            DA12 590
      Y=Y+.7489094100646149D-1*FCT(X)                                   DA12 600
      X=.12695899401039615D1                                            DA12 610
      Y=Y+.25547924356911832D0*FCT(X)                                   DA12 620
      X=.45450668156378028D0                                            DA12 630
      Y=Y+.57235907069288604D0*FCT(X)                                   DA12 640
      X=.50361889117293951D-1                                           DA12 650
      Y=Y+.8538623277373985D0*FCT(X)                                    DA12 660
      RETURN                                                            DA12 590
      END                                                               DA12 600
