C                                                                       DH24  10
C     ..................................................................DH24  20
C                                                                       DH24  30
C        SUBROUTINE DQH24                                               DH24  40
C                                                                       DH24  50
C        PURPOSE                                                        DH24  60
C           TO COMPUTE INTEGRAL(EXP(-X*X)*FCT(X), SUMMED OVER X FROM    DH24  70
C                               -INFINITY TO +INFINITY).                DH24  80
C                                                                       DH24  90
C        USAGE                                                          DH24 100
C           CALL DQH24 (FCT,Y)                                          DH24 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DH24 120
C                                                                       DH24 130
C        DESCRIPTION OF PARAMETERS                                      DH24 140
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DH24 150
C                    SUBPROGRAM USED.                                   DH24 160
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DH24 170
C                                                                       DH24 180
C        REMARKS                                                        DH24 190
C           NONE                                                        DH24 200
C                                                                       DH24 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DH24 220
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DH24 230
C           MUST BE FURNISHED BY THE USER.                              DH24 240
C                                                                       DH24 250
C        METHOD                                                         DH24 260
C           EVALUATION IS DONE BY MEANS OF 24-POINT GAUSSIAN-HERMITE    DH24 270
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       DH24 280
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 47.                     DH24 290
C           FOR REFERENCE, SEE                                          DH24 300
C           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF    DH24 310
C           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED     DH24 320
C           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT       DH24 330
C           TR00.1100 (MARCH 1964), PP.213-214.                         DH24 340
C                                                                       DH24 350
C     ..................................................................DH24 360
C                                                                       DH24 370
      SUBROUTINE DQH24(FCT,Y)                                           DH24 380
C                                                                       DH24 390
C                                                                       DH24 400
      DOUBLE PRECISION X,Y,Z,FCT                                        DH24 410
C                                                                       DH24 420
      X=.60159255614257397D1                                            DH24 430
      Z=-X                                                              DH24 440
      Y=.16643684964891089D-15*(FCT(X)+FCT(Z))                          DH24 450
      X=.52593829276680444D1                                            DH24 460
      Z=-X                                                              DH24 470
      Y=Y+.65846202430781701D-12*(FCT(X)+FCT(Z))                        DH24 480
      X=.46256627564237873D1                                            DH24 490
      Z=-X                                                              DH24 500
      Y=Y+.30462542699875639D-9*(FCT(X)+FCT(Z))                         DH24 510
      X=.40536644024481495D1                                            DH24 520
      Z=-X                                                              DH24 530
      Y=Y+.40189711749414297D-7*(FCT(X)+FCT(Z))                         DH24 540
      X=.35200068130345247D1                                            DH24 550
      Z=-X                                                              DH24 560
      Y=Y+.21582457049023336D-5*(FCT(X)+FCT(Z))                         DH24 570
      X=.30125461375655648D1                                            DH24 580
      Z=-X                                                              DH24 590
      Y=Y+.56886916364043798D-4*(FCT(X)+FCT(Z))                         DH24 600
      X=.25238810170114270D1                                            DH24 610
      Z=-X                                                              DH24 620
      Y=Y+.8236924826884175D-3*(FCT(X)+FCT(Z))                          DH24 630
      X=.20490035736616989D1                                            DH24 640
      Z=-X                                                              DH24 650
      Y=Y+.70483558100726710D-2*(FCT(X)+FCT(Z))                         DH24 660
      X=.15842500109616941D1                                            DH24 670
      Z=-X                                                              DH24 680
      Y=Y+.37445470503230746D-1*(FCT(X)+FCT(Z))                         DH24 690
      X=.11267608176112451D1                                            DH24 700
      Z=-X                                                              DH24 710
      Y=Y+.12773962178455916D0*(FCT(X)+FCT(Z))                          DH24 720
      X=.67417110703721224D0                                            DH24 730
      Z=-X                                                              DH24 740
      Y=Y+.28617953534644302D0*(FCT(X)+FCT(Z))                          DH24 750
      X=.22441454747251559D0                                            DH24 760
      Z=-X                                                              DH24 770
      Y=Y+.42693116386869925D0*(FCT(X)+FCT(Z))                          DH24 780
      RETURN                                                            DH24 790
      END                                                               DH24 800
