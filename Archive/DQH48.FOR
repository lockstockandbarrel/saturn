C                                                                       DH48  10
C     ..................................................................DH48  20
C                                                                       DH48  30
C        SUBROUTINE DQH48                                               DH48  40
C                                                                       DH48  50
C        PURPOSE                                                        DH48  60
C           TO COMPUTE INTEGRAL(EXP(-X*X)*FCT(X), SUMMED OVER X FROM    DH48  70
C                               -INFINITY TO +INFINITY).                DH48  80
C                                                                       DH48  90
C        USAGE                                                          DH48 100
C           CALL DQH48 (FCT,Y)                                          DH48 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DH48 120
C                                                                       DH48 130
C        DESCRIPTION OF PARAMETERS                                      DH48 140
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DH48 150
C                    SUBPROGRAM USED.                                   DH48 160
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DH48 170
C                                                                       DH48 180
C        REMARKS                                                        DH48 190
C           NONE                                                        DH48 200
C                                                                       DH48 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DH48 220
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DH48 230
C           MUST BE FURNISHED BY THE USER.                              DH48 240
C                                                                       DH48 250
C        METHOD                                                         DH48 260
C           EVALUATION IS DONE BY MEANS OF 48-POINT GAUSSIAN-HERMITE    DH48 270
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       DH48 280
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 95.                     DH48 290
C           FOR REFERENCE, SEE                                          DH48 300
C           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF    DH48 310
C           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED     DH48 320
C           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT       DH48 330
C           TR00.1100 (MARCH 1964), PP.213-214.                         DH48 340
C                                                                       DH48 350
C     ..................................................................DH48 360
C                                                                       DH48 370
      SUBROUTINE DQH48(FCT,Y)                                           DH48 380
C                                                                       DH48 390
C                                                                       DH48 400
      DOUBLE PRECISION X,Y,Z,FCT                                        DH48 410
C                                                                       DH48 420
      X=.8975315081931687D1                                             DH48 430
      Z=-X                                                              DH48 440
      Y=.7935551460773997D-35*(FCT(X)+FCT(Z))                           DH48 450
      X=.8310752190704784D1                                             DH48 460
      Z=-X                                                              DH48 470
      Y=Y+.59846126933138784D-30*(FCT(X)+FCT(Z))                        DH48 480
      X=.7759295519765775D1                                             DH48 490
      Z=-X                                                              DH48 500
      Y=Y+.36850360801506699D-26*(FCT(X)+FCT(Z))                        DH48 510
      X=.7266046554164350D1                                             DH48 520
      Z=-X                                                              DH48 530
      Y=Y+.55645774689022848D-23*(FCT(X)+FCT(Z))                        DH48 540
      X=.68100645780741414D1                                            DH48 550
      Z=-X                                                              DH48 560
      Y=Y+.31883873235051384D-20*(FCT(X)+FCT(Z))                        DH48 570
      X=.63805640961864106D1                                            DH48 580
      Z=-X                                                              DH48 590
      Y=Y+.8730159601186677D-18*(FCT(X)+FCT(Z))                         DH48 600
      X=.59710722250135454D1                                            DH48 610
      Z=-X                                                              DH48 620
      Y=Y+.13151596226584085D-15*(FCT(X)+FCT(Z))                        DH48 630
      X=.55773169812237286D1                                            DH48 640
      Z=-X                                                              DH48 650
      Y=Y+.11975898654791794D-13*(FCT(X)+FCT(Z))                        DH48 660
      X=.51962877187923645D1                                            DH48 670
      Z=-X                                                              DH48 680
      Y=Y+.70469325815458891D-12*(FCT(X)+FCT(Z))                        DH48 690
      X=.48257572281332095D1                                            DH48 700
      Z=-X                                                              DH48 710
      Y=Y+.28152965378381691D-10*(FCT(X)+FCT(Z))                        DH48 720
      X=.44640145469344589D1                                            DH48 730
      Z=-X                                                              DH48 740
      Y=Y+.7930467495165382D-9*(FCT(X)+FCT(Z))                          DH48 750
      X=.41097046035605902D1                                            DH48 760
      Z=-X                                                              DH48 770
      Y=Y+.16225141358957698D-7*(FCT(X)+FCT(Z))                         DH48 780
      X=.37617264902283578D1                                            DH48 790
      Z=-X                                                              DH48 800
      Y=Y+.24686589936697505D-6*(FCT(X)+FCT(Z))                         DH48 810
      X=.34191659693638846D1                                            DH48 820
      Z=-X                                                              DH48 830
      Y=Y+.28472586917348481D-5*(FCT(X)+FCT(Z))                         DH48 840
      X=.30812489886451058D1                                            DH48 850
      Z=-X                                                              DH48 860
      Y=Y+.25285990277484889D-4*(FCT(X)+FCT(Z))                         DH48 870
      X=.27473086248223832D1                                            DH48 880
      Z=-X                                                              DH48 890
      Y=Y+.17515043180117283D-3*(FCT(X)+FCT(Z))                         DH48 900
      X=.24167609048732165D1                                            DH48 910
      Z=-X                                                              DH48 920
      Y=Y+.9563923198194153D-3*(FCT(X)+FCT(Z))                          DH48 930
      X=.20890866609442764D1                                            DH48 940
      Z=-X                                                              DH48 950
      Y=Y+.41530049119775525D-2*(FCT(X)+FCT(Z))                         DH48 960
      X=.17638175798953000D1                                            DH48 970
      Z=-X                                                              DH48 980
      Y=Y+.14444961574981099D-1*(FCT(X)+FCT(Z))                         DH48 990
      X=.14405252201375652D1                                            DH481000
      Z=-X                                                              DH481010
      Y=Y+.40479676984603849D-1*(FCT(X)+FCT(Z))                         DH481020
      X=.11188121524021566D1                                            DH481030
      Z=-X                                                              DH481040
      Y=Y+.9182229707928518D-1*(FCT(X)+FCT(Z))                          DH481050
      X=.7983046277785622D0                                             DH481060
      Z=-X                                                              DH481070
      Y=Y+.16920447194564111D0*(FCT(X)+FCT(Z))                          DH481080
      X=.47864633759449610D0                                            DH481090
      Z=-X                                                              DH481100
      Y=Y+.25396154266475910D0*(FCT(X)+FCT(Z))                          DH481110
      X=.15949293584886247D0                                            DH481120
      Z=-X                                                              DH481130
      Y=Y+.31100103037796308D0*(FCT(X)+FCT(Z))                          DH481140
      RETURN                                                            DH481150
      END                                                               DH481160
