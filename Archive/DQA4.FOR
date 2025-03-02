C                                                                       DQA4  10
C     ..................................................................DQA4  20
C                                                                       DQA4  30
C        SUBROUTINE DQA4                                                DQA4  40
C                                                                       DQA4  50
C        PURPOSE                                                        DQA4  60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X)/SQRT(X), SUMMED OVER X   DQA4  70
C                               FROM 0 TO INFINITY).                    DQA4  80
C                                                                       DQA4  90
C        USAGE                                                          DQA4 100
C           CALL DQA4 (FCT,Y)                                           DQA4 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DQA4 120
C                                                                       DQA4 130
C        DESCRIPTION OF PARAMETERS                                      DQA4 140
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DQA4 150
C                    SUBPROGRAM USED.                                   DQA4 160
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DQA4 170
C                                                                       DQA4 180
C        REMARKS                                                        DQA4 190
C           NONE                                                        DQA4 200
C                                                                       DQA4 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DQA4 220
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DQA4 230
C           MUST BE FURNISHED BY THE USER.                              DQA4 240
C                                                                       DQA4 250
C        METHOD                                                         DQA4 260
C           EVALUATION IS DONE BY MEANS OF 4-POINT GENERALIZED GAUSS-   DQA4 270
C           LAGUERRE QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY       DQA4 280
C           WHENEVER FCT(X) IS A POLYNOMIAL UP TO DEGREE 7.             DQA4 290
C           FOR REFERENCE, SEE                                          DQA4 300
C           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF    DQA4 310
C           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED     DQA4 320
C           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT       DQA4 330
C           TR00.1100 (MARCH 1964), PP.15-16.                           DQA4 340
C                                                                       DQA4 350
C     ..................................................................DQA4 360
C                                                                       DQA4 370
      SUBROUTINE DQA4(FCT,Y)                                            DQA4 380
C                                                                       DQA4 390
C                                                                       DQA4 400
      DOUBLE PRECISION X,Y,FCT                                          DQA4 410
C                                                                       DQA4 420
      X=.8588635689012034D1                                             DQA4 430
      Y=.39920814442273524D-3*FCT(X)                                    DQA4 440
      X=.39269635013582872D1                                            DQA4 450
      Y=Y+.34155966014826951D-1*FCT(X)                                  DQA4 460
      X=.13390972881263614D1                                            DQA4 470
      Y=Y+.41560465162978376D0*FCT(X)                                   DQA4 480
      X=.14530352150331709D0                                            DQA4 490
      Y=Y+.13222940251164826D1*FCT(X)                                   DQA4 500
      RETURN                                                            DQA4 510
      END                                                               DQA4 520
