C                                                                       DA16  10
C     ..................................................................DA16  20
C                                                                       DA16  30
C        SUBROUTINE DQA16                                               DA16  40
C                                                                       DA16  50
C        PURPOSE                                                        DA16  60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X)/SQRT(X), SUMMED OVER X   DA16  70
C                               FROM 0 TO INFINITY).                    DA16  80
C                                                                       DA16  90
C        USAGE                                                          DA16 100
C           CALL DQA16 (FCT,Y)                                          DA16 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DA16 120
C                                                                       DA16 130
C        DESCRIPTION OF PARAMETERS                                      DA16 140
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DA16 150
C                    SUBPROGRAM USED.                                   DA16 160
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DA16 170
C                                                                       DA16 180
C        REMARKS                                                        DA16 190
C           NONE                                                        DA16 200
C                                                                       DA16 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DA16 220
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DA16 230
C           MUST BE FURNISHED BY THE USER.                              DA16 240
C                                                                       DA16 250
C        METHOD                                                         DA16 260
C           EVALUATION IS DONE BY MEANS OF 16-POINT GENERALIZED GAUSS-  DA16 270
C           LAGUERRE QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY       DA16 280
C           WHENEVER FCT(X) IS A POLYNOMIAL UP TO DEGREE 31.            DA16 290
C           FOR REFERENCE, SEE                                          DA16 300
C           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF    DA16 310
C           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED     DA16 320
C           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT       DA16 330
C           TR00.1100 (MARCH 1964), PP.15-16.                           DA16 340
C                                                                       DA16 350
C     ..................................................................DA16 360
C                                                                       DA16 370
      SUBROUTINE DQA16(FCT,Y)                                           DA16 380
C                                                                       DA16 390
C                                                                       DA16 400
      DOUBLE PRECISION X,Y,FCT                                          DA16 410
C                                                                       DA16 420
      X=.50777223877537080D2                                            DA16 430
      Y=.14621352854768325D-21*FCT(X)                                   DA16 440
      X=.41081666525491202D2                                            DA16 450
      Y=Y+.18463473073036584D-17*FCT(X)                                 DA16 460
      X=.33781970488226166D2                                            DA16 470
      Y=Y+.23946880341856973D-14*FCT(X)                                 DA16 480
      X=.27831438211328676D2                                            DA16 490
      Y=Y+.8430020422652895D-12*FCT(X)                                  DA16 500
      X=.22821300693525208D2                                            DA16 510
      Y=Y+.11866582926793277D-9*FCT(X)                                  DA16 520
      X=.18537743178606694D2                                            DA16 530
      Y=Y+.8197664329541793D-8*FCT(X)                                   DA16 540
      X=.14851431341801250D2                                            DA16 550
      Y=Y+.31483355850911881D-6*FCT(X)                                  DA16 560
      X=.11677033673975957D2                                            DA16 570
      Y=Y+.7301170259124752D-5*FCT(X)                                   DA16 580
      X=.8955001337723390D1                                             DA16 590
      Y=Y+.10833168123639965D-3*FCT(X)                                  DA16 600
      X=.66422151797414440D1                                            DA16 610
      Y=Y+.10725367310559441D-2*FCT(X)                                  DA16 620
      X=.47067267076675872D1                                            DA16 630
      Y=Y+.7309780653308856D-2*FCT(X)                                   DA16 640
      X=.31246010507021443D1                                            DA16 650
      Y=Y+.35106857663146861D-1*FCT(X)                                  DA16 660
      X=.18779315076960743D1                                            DA16 670
      Y=Y+.12091626191182523D0*FCT(X)                                   DA16 680
      X=.9535531553908655D0                                             DA16 690
      Y=Y+.30253946815328497D0*FCT(X)                                   DA16 700
      X=.34220015601094768D0                                            DA16 710
      Y=Y+.55491628460505980D0*FCT(X)                                   DA16 720
      X=.37962914575313455D-1                                           DA16 730
      Y=Y+.7504767051856048D0*FCT(X)                                    DA16 740
      RETURN                                                            DA16 750
      END                                                               DA16 760
