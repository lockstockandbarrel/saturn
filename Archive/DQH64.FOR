C                                                                       DH64  10
C     ..................................................................DH64  20
C                                                                       DH64  30
C        SUBROUTINE DQH64                                               DH64  40
C                                                                       DH64  50
C        PURPOSE                                                        DH64  60
C           TO COMPUTE INTEGRAL(EXP(-X*X)*FCT(X), SUMMED OVER X FROM    DH64  70
C                               -INFINITY TO +INFINITY).                DH64  80
C                                                                       DH64  90
C        USAGE                                                          DH64 100
C           CALL DQH64 (FCT,Y)                                          DH64 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DH64 120
C                                                                       DH64 130
C        DESCRIPTION OF PARAMETERS                                      DH64 140
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DH64 150
C                    SUBPROGRAM USED.                                   DH64 160
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DH64 170
C                                                                       DH64 180
C        REMARKS                                                        DH64 190
C           NONE                                                        DH64 200
C                                                                       DH64 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DH64 220
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DH64 230
C           MUST BE FURNISHED BY THE USER.                              DH64 240
C                                                                       DH64 250
C        METHOD                                                         DH64 260
C           EVALUATION IS DONE BY MEANS OF 64-POINT GAUSSIAN-HERMITE    DH64 270
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       DH64 280
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 127.                    DH64 290
C           FOR REFERENCE, SEE                                          DH64 300
C           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF    DH64 310
C           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED     DH64 320
C           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT       DH64 330
C           TR00.1100 (MARCH 1964), PP.213-214.                         DH64 340
C                                                                       DH64 350
C     ..................................................................DH64 360
C                                                                       DH64 370
      SUBROUTINE DQH64(FCT,Y)                                           DH64 380
C                                                                       DH64 390
C                                                                       DH64 400
      DOUBLE PRECISION X,Y,Z,FCT                                        DH64 410
C                                                                       DH64 420
      X=.10526123167960546D2                                            DH64 430
      Z=-X                                                              DH64 440
      Y=.55357065358569428D-48*(FCT(X)+FCT(Z))                          DH64 450
      X=.9895287586829539D1                                             DH64 460
      Z=-X                                                              DH64 470
      Y=Y+.16797479901081592D-42*(FCT(X)+FCT(Z))                        DH64 480
      X=.9373159549646721D1                                             DH64 490
      Z=-X                                                              DH64 500
      Y=Y+.34211380112557405D-38*(FCT(X)+FCT(Z))                        DH64 510
      X=.8907249099964770D1                                             DH64 520
      Z=-X                                                              DH64 530
      Y=Y+.15573906246297638D-34*(FCT(X)+FCT(Z))                        DH64 540
      X=.8477529083379863D1                                             DH64 550
      Z=-X                                                              DH64 560
      Y=Y+.25496608991129993D-31*(FCT(X)+FCT(Z))                        DH64 570
      X=.8073687285010225D1                                             DH64 580
      Z=-X                                                              DH64 590
      Y=Y+.19291035954649669D-28*(FCT(X)+FCT(Z))                        DH64 600
      X=.7689540164040497D1                                             DH64 610
      Z=-X                                                              DH64 620
      Y=Y+.7861797788925910D-26*(FCT(X)+FCT(Z))                         DH64 630
      X=.7321013032780949D1                                             DH64 640
      Z=-X                                                              DH64 650
      Y=Y+.19117068833006428D-23*(FCT(X)+FCT(Z))                        DH64 660
      X=.69652411205511075D1                                            DH64 670
      Z=-X                                                              DH64 680
      Y=Y+.29828627842798512D-21*(FCT(X)+FCT(Z))                        DH64 690
      X=.66201122626360274D1                                            DH64 700
      Z=-X                                                              DH64 710
      Y=Y+.31522545665037814D-19*(FCT(X)+FCT(Z))                        DH64 720
      X=.62840112287748282D1                                            DH64 730
      Z=-X                                                              DH64 740
      Y=Y+.23518847106758191D-17*(FCT(X)+FCT(Z))                        DH64 750
      X=.59556663267994860D1                                            DH64 760
      Z=-X                                                              DH64 770
      Y=Y+.12800933913224380D-15*(FCT(X)+FCT(Z))                        DH64 780
      X=.56340521643499721D1                                            DH64 790
      Z=-X                                                              DH64 800
      Y=Y+.52186237265908475D-14*(FCT(X)+FCT(Z))                        DH64 810
      X=.53183252246332709D1                                            DH64 820
      Z=-X                                                              DH64 830
      Y=Y+.16283407307097204D-12*(FCT(X)+FCT(Z))                        DH64 840
      X=.50077796021987682D1                                            DH64 850
      Z=-X                                                              DH64 860
      Y=Y+.39591777669477239D-11*(FCT(X)+FCT(Z))                        DH64 870
      X=.47018156474074998D1                                            DH64 880
      Z=-X                                                              DH64 890
      Y=Y+.7615217250145451D-10*(FCT(X)+FCT(Z))                         DH64 900
      X=.43999171682281376D1                                            DH64 910
      Z=-X                                                              DH64 920
      Y=Y+.11736167423215493D-8*(FCT(X)+FCT(Z))                         DH64 930
      X=.41016344745666567D1                                            DH64 940
      Z=-X                                                              DH64 950
      Y=Y+.14651253164761094D-7*(FCT(X)+FCT(Z))                         DH64 960
      X=.38065715139453605D1                                            DH64 970
      Z=-X                                                              DH64 980
      Y=Y+.14955329367272471D-6*(FCT(X)+FCT(Z))                         DH64 990
      X=.35143759357409062D1                                            DH641000
      Z=-X                                                              DH641010
      Y=Y+.12583402510311846D-5*(FCT(X)+FCT(Z))                         DH641020
      X=.32247312919920357D1                                            DH641030
      Z=-X                                                              DH641040
      Y=Y+.8788499230850359D-5*(FCT(X)+FCT(Z))                          DH641050
      X=.29373508230046218D1                                            DH641060
      Z=-X                                                              DH641070
      Y=Y+.51259291357862747D-4*(FCT(X)+FCT(Z))                         DH641080
      X=.26519724354306350D1                                            DH641090
      Z=-X                                                              DH641100
      Y=Y+.25098369851306249D-3*(FCT(X)+FCT(Z))                         DH641110
      X=.23683545886324014D1                                            DH641120
      Z=-X                                                              DH641130
      Y=Y+.10363290995075777D-2*(FCT(X)+FCT(Z))                         DH641140
      X=.20862728798817620D1                                            DH641150
      Z=-X                                                              DH641160
      Y=Y+.36225869785344588D-2*(FCT(X)+FCT(Z))                         DH641170
      X=.18055171714655449D1                                            DH641180
      Z=-X                                                              DH641190
      Y=Y+.10756040509879137D-1*(FCT(X)+FCT(Z))                         DH641200
      X=.15258891402098637D1                                            DH641210
      Z=-X                                                              DH641220
      Y=Y+.27203128953688918D-1*(FCT(X)+FCT(Z))                         DH641230
      X=.12472001569431179D1                                            DH641240
      Z=-X                                                              DH641250
      Y=Y+.58739981964099435D-1*(FCT(X)+FCT(Z))                         DH641260
      X=.9692694230711780D0                                             DH641270
      Z=-X                                                              DH641280
      Y=Y+.10849834930618684D0*(FCT(X)+FCT(Z))                          DH641290
      X=.69192230581004458D0                                            DH641300
      Z=-X                                                              DH641310
      Y=Y+.17168584234908370D0*(FCT(X)+FCT(Z))                          DH641320
      X=.41498882412107868D0                                            DH641330
      Z=-X                                                              DH641340
      Y=Y+.23299478606267805D0*(FCT(X)+FCT(Z))                          DH641350
      X=.13830224498700972D0                                            DH641360
      Z=-X                                                              DH641370
      Y=Y+.27137742494130398D0*(FCT(X)+FCT(Z))                          DH641380
      RETURN                                                            DH641390
      END                                                               DH641400
