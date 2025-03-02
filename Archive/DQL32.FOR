C                                                                       DL32  10
C     ..................................................................DL32  20
C                                                                       DL32  30
C        SUBROUTINE DQL32                                               DL32  40
C                                                                       DL32  50
C        PURPOSE                                                        DL32  60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X), SUMMED OVER X           DL32  70
C                               FROM 0 TO INFINITY).                    DL32  80
C                                                                       DL32  90
C        USAGE                                                          DL32 100
C           CALL DQL32 (FCT,Y)                                          DL32 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DL32 120
C                                                                       DL32 130
C        DESCRIPTION OF PARAMETERS                                      DL32 140
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DL32 150
C                    SUBPROGRAM USED.                                   DL32 160
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DL32 170
C                                                                       DL32 180
C        REMARKS                                                        DL32 190
C           NONE                                                        DL32 200
C                                                                       DL32 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DL32 220
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DL32 230
C           MUST BE FURNISHED BY THE USER.                              DL32 240
C                                                                       DL32 250
C        METHOD                                                         DL32 260
C           EVALUATION IS DONE BY MEANS OF 32-POINT GAUSSIAN-LAGUERRE   DL32 270
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY,               DL32 280
C           WHENEVER FCT(X) IS A POLYNOMIAL UP TO DEGREE 63.            DL32 290
C           FOR REFERENCE, SEE                                          DL32 300
C           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF    DL32 310
C           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED     DL32 320
C           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT       DL32 330
C           TR00.1100 (MARCH 1964), PP.24-25.                           DL32 340
C                                                                       DL32 350
C     ..................................................................DL32 360
C                                                                       DL32 370
      SUBROUTINE DQL32(FCT,Y)                                           DL32 380
C                                                                       DL32 390
C                                                                       DL32 400
      DOUBLE PRECISION X,Y,FCT                                          DL32 410
C                                                                       DL32 420
      X=.11175139809793770D3                                            DL32 430
      Y=.45105361938989742D-47*FCT(X)                                   DL32 440
      X=.9882954286828397D2                                             DL32 450
      Y=Y+.13386169421062563D-41*FCT(X)                                 DL32 460
      X=.8873534041789240D2                                             DL32 470
      Y=Y+.26715112192401370D-37*FCT(X)                                 DL32 480
      X=.8018744697791352D2                                             DL32 490
      Y=Y+.11922487600982224D-33*FCT(X)                                 DL32 500
      X=.7268762809066271D2                                             DL32 510
      Y=Y+.19133754944542243D-30*FCT(X)                                 DL32 520
      X=.65975377287935053D2                                            DL32 530
      Y=Y+.14185605454630369D-27*FCT(X)                                 DL32 540
      X=.59892509162134018D2                                            DL32 550
      Y=Y+.56612941303973594D-25*FCT(X)                                 DL32 560
      X=.54333721333396907D2                                            DL32 570
      Y=Y+.13469825866373952D-22*FCT(X)                                 DL32 580
      X=.49224394987308639D2                                            DL32 590
      Y=Y+.20544296737880454D-20*FCT(X)                                 DL32 600
      X=.44509207995754938D2                                            DL32 610
      Y=Y+.21197922901636186D-18*FCT(X)                                 DL32 620
      X=.40145719771539442D2                                            DL32 630
      Y=Y+.15421338333938234D-16*FCT(X)                                 DL32 640
      X=.36100494805751974D2                                            DL32 650
      Y=Y+.8171823443420719D-15*FCT(X)                                  DL32 660
      X=.32346629153964737D2                                            DL32 670
      Y=Y+.32378016577292665D-13*FCT(X)                                 DL32 680
      X=.28862101816323475D2                                            DL32 690
      Y=Y+.9799379288727094D-12*FCT(X)                                  DL32 700
      X=.25628636022459248D2                                            DL32 710
      Y=Y+.23058994918913361D-10*FCT(X)                                 DL32 720
      X=.22630889013196774D2                                            DL32 730
      Y=Y+.42813829710409289D-9*FCT(X)                                  DL32 740
      X=.19855860940336055D2                                            DL32 750
      Y=Y+.63506022266258067D-8*FCT(X)                                  DL32 760
      X=.17292454336715315D2                                            DL32 770
      Y=Y+.7604567879120781D-7*FCT(X)                                   DL32 780
      X=.14931139755522557D2                                            DL32 790
      Y=Y+.7416404578667552D-6*FCT(X)                                   DL32 800
      X=.12763697986742725D2                                            DL32 810
      Y=Y+.59345416128686329D-5*FCT(X)                                  DL32 820
      X=.10783018632539972D2                                            DL32 830
      Y=Y+.39203419679879472D-4*FCT(X)                                  DL32 840
      X=.8982940924212596D1                                             DL32 850
      Y=Y+.21486491880136419D-3*FCT(X)                                  DL32 860
      X=.7358126733186241D1                                             DL32 870
      Y=Y+.9808033066149551D-3*FCT(X)                                   DL32 880
      X=.59039585041742439D1                                            DL32 890
      Y=Y+.37388162946115248D-2*FCT(X)                                  DL32 900
      X=.46164567697497674D1                                            DL32 910
      Y=Y+.11918214834838557D-1*FCT(X)                                  DL32 920
      X=.34922132730219945D1                                            DL32 930
      Y=Y+.31760912509175070D-1*FCT(X)                                  DL32 940
      X=.25283367064257949D1                                            DL32 950
      Y=Y+.70578623865717442D-1*FCT(X)                                  DL32 960
      X=.17224087764446454D1                                            DL32 970
      Y=Y+.12998378628607176D0*FCT(X)                                   DL32 980
      X=.10724487538178176D1                                            DL32 990
      Y=Y+.19590333597288104D0*FCT(X)                                   DL321000
      X=.57688462930188643D0                                            DL321010
      Y=Y+.23521322966984801D0*FCT(X)                                   DL321020
      X=.23452610951961854D0                                            DL321030
      Y=Y+.21044310793881323D0*FCT(X)                                   DL321040
      X=.44489365833267018D-1                                           DL321050
      Y=Y+.10921834195238497D0*FCT(X)                                   DL321060
      RETURN                                                            DL321070
      END                                                               DL321080
