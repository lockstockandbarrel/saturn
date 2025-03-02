C                                                                       DH32  10
C     ..................................................................DH32  20
C                                                                       DH32  30
C        SUBROUTINE DQH32                                               DH32  40
C                                                                       DH32  50
C        PURPOSE                                                        DH32  60
C           TO COMPUTE INTEGRAL(EXP(-X*X)*FCT(X), SUMMED OVER X FROM    DH32  70
C                               -INFINITY TO +INFINITY).                DH32  80
C                                                                       DH32  90
C        USAGE                                                          DH32 100
C           CALL DQH32 (FCT,Y)                                          DH32 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DH32 120
C                                                                       DH32 130
C        DESCRIPTION OF PARAMETERS                                      DH32 140
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DH32 150
C                    SUBPROGRAM USED.                                   DH32 160
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DH32 170
C                                                                       DH32 180
C        REMARKS                                                        DH32 190
C           NONE                                                        DH32 200
C                                                                       DH32 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DH32 220
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DH32 230
C           MUST BE FURNISHED BY THE USER.                              DH32 240
C                                                                       DH32 250
C        METHOD                                                         DH32 260
C           EVALUATION IS DONE BY MEANS OF 32-POINT GAUSSIAN-HERMITE    DH32 270
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       DH32 280
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 63.                     DH32 290
C           FOR REFERENCE, SEE                                          DH32 300
C           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF    DH32 310
C           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED     DH32 320
C           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT       DH32 330
C           TR00.1100 (MARCH 1964), PP.213-214.                         DH32 340
C                                                                       DH32 350
C     ..................................................................DH32 360
C                                                                       DH32 370
      SUBROUTINE DQH32(FCT,Y)                                           DH32 380
C                                                                       DH32 390
C                                                                       DH32 400
      DOUBLE PRECISION X,Y,Z,FCT                                        DH32 410
C                                                                       DH32 420
      X=.71258139098307276D1                                            DH32 430
      Z=-X                                                              DH32 440
      Y=.7310676427384162D-22*(FCT(X)+FCT(Z))                           DH32 450
      X=.64094981492696604D1                                            DH32 460
      Z=-X                                                              DH32 470
      Y=Y+.9231736536518292D-18*(FCT(X)+FCT(Z))                         DH32 480
      X=.58122259495159138D1                                            DH32 490
      Z=-X                                                              DH32 500
      Y=Y+.11973440170928487D-14*(FCT(X)+FCT(Z))                        DH32 510
      X=.52755509865158801D1                                            DH32 520
      Z=-X                                                              DH32 530
      Y=Y+.42150102113264476D-12*(FCT(X)+FCT(Z))                        DH32 540
      X=.47771645035025964D1                                            DH32 550
      Z=-X                                                              DH32 560
      Y=Y+.59332914633966386D-10*(FCT(X)+FCT(Z))                        DH32 570
      X=.43055479533511984D1                                            DH32 580
      Z=-X                                                              DH32 590
      Y=Y+.40988321647708966D-8*(FCT(X)+FCT(Z))                         DH32 600
      X=.38537554854714446D1                                            DH32 610
      Z=-X                                                              DH32 620
      Y=Y+.15741677925455940D-6*(FCT(X)+FCT(Z))                         DH32 630
      X=.34171674928185707D1                                            DH32 640
      Z=-X                                                              DH32 650
      Y=Y+.36505851295623761D-5*(FCT(X)+FCT(Z))                         DH32 660
      X=.29924908250023742D1                                            DH32 670
      Z=-X                                                              DH32 680
      Y=Y+.54165840618199826D-4*(FCT(X)+FCT(Z))                         DH32 690
      X=.25772495377323175D1                                            DH32 700
      Z=-X                                                              DH32 710
      Y=Y+.53626836552797205D-3*(FCT(X)+FCT(Z))                         DH32 720
      X=.21694991836061122D1                                            DH32 730
      Z=-X                                                              DH32 740
      Y=Y+.36548903266544281D-2*(FCT(X)+FCT(Z))                         DH32 750
      X=.17676541094632016D1                                            DH32 760
      Z=-X                                                              DH32 770
      Y=Y+.17553428831573430D-1*(FCT(X)+FCT(Z))                         DH32 780
      X=.13703764109528718D1                                            DH32 790
      Z=-X                                                              DH32 800
      Y=Y+.60458130955912614D-1*(FCT(X)+FCT(Z))                         DH32 810
      X=.9765004635896828D0                                             DH32 820
      Z=-X                                                              DH32 830
      Y=Y+.15126973407664248D0*(FCT(X)+FCT(Z))                          DH32 840
      X=.58497876543593245D0                                            DH32 850
      Z=-X                                                              DH32 860
      Y=Y+.27745814230252990D0*(FCT(X)+FCT(Z))                          DH32 870
      X=.19484074156939933D0                                            DH32 880
      Z=-X                                                              DH32 890
      Y=Y+.37523835259280239D0*(FCT(X)+FCT(Z))                          DH32 900
      RETURN                                                            DH32 910
      END                                                               DH32 920
