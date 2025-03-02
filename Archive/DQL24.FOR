C                                                                       DL24  10
C     ..................................................................DL24  20
C                                                                       DL24  30
C        SUBROUTINE DQL24                                               DL24  40
C                                                                       DL24  50
C        PURPOSE                                                        DL24  60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X), SUMMED OVER X           DL24  70
C                               FROM 0 TO INFINITY).                    DL24  80
C                                                                       DL24  90
C        USAGE                                                          DL24 100
C           CALL DQL24 (FCT,Y)                                          DL24 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DL24 120
C                                                                       DL24 130
C        DESCRIPTION OF PARAMETERS                                      DL24 140
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DL24 150
C                    SUBPROGRAM USED.                                   DL24 160
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DL24 170
C                                                                       DL24 180
C        REMARKS                                                        DL24 190
C           NONE                                                        DL24 200
C                                                                       DL24 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DL24 220
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DL24 230
C           MUST BE FURNISHED BY THE USER.                              DL24 240
C                                                                       DL24 250
C        METHOD                                                         DL24 260
C           EVALUATION IS DONE BY MEANS OF 24-POINT GAUSSIAN-LAGUERRE   DL24 270
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY,               DL24 280
C           WHENEVER FCT(X) IS A POLYNOMIAL UP TO DEGREE 47.            DL24 290
C           FOR REFERENCE, SEE                                          DL24 300
C           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF    DL24 310
C           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED     DL24 320
C           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT       DL24 330
C           TR00.1100 (MARCH 1964), PP.24-25.                           DL24 340
C                                                                       DL24 350
C     ..................................................................DL24 360
C                                                                       DL24 370
      SUBROUTINE DQL24(FCT,Y)                                           DL24 380
C                                                                       DL24 390
C                                                                       DL24 400
      DOUBLE PRECISION X,Y,FCT                                          DL24 410
C                                                                       DL24 420
      X=.8149827923394889D2                                             DL24 430
      Y=.55753457883283568D-34*FCT(X)                                   DL24 440
      X=.69962240035105030D2                                            DL24 450
      Y=Y+.40883015936806578D-29*FCT(X)                                 DL24 460
      X=.61058531447218762D2                                            DL24 470
      Y=Y+.24518188458784027D-25*FCT(X)                                 DL24 480
      X=.53608574544695070D2                                            DL24 490
      Y=Y+.36057658645529590D-22*FCT(X)                                 DL24 500
      X=.47153106445156323D2                                            DL24 510
      Y=Y+.20105174645555035D-19*FCT(X)                                 DL24 520
      X=.41451720484870767D2                                            DL24 530
      Y=Y+.53501888130100376D-17*FCT(X)                                 DL24 540
      X=.36358405801651622D2                                            DL24 550
      Y=Y+.7819800382459448D-15*FCT(X)                                  DL24 560
      X=.31776041352374723D2                                            DL24 570
      Y=Y+.68941810529580857D-13*FCT(X)                                 DL24 580
      X=.27635937174332717D2                                            DL24 590
      Y=Y+.39177365150584514D-11*FCT(X)                                 DL24 600
      X=.23887329848169733D2                                            DL24 610
      Y=Y+.15070082262925849D-9*FCT(X)                                  DL24 620
      X=.20491460082616425D2                                            DL24 630
      Y=Y+.40728589875499997D-8*FCT(X)                                  DL24 640
      X=.17417992646508979D2                                            DL24 650
      Y=Y+.7960812959133630D-7*FCT(X)                                   DL24 660
      X=.14642732289596674D2                                            DL24 670
      Y=Y+.11513158127372799D-5*FCT(X)                                  DL24 680
      X=.12146102711729766D2                                            DL24 690
      Y=Y+.12544721977993333D-4*FCT(X)                                  DL24 700
      X=.9912098015077706D1                                             DL24 710
      Y=Y+.10446121465927518D-3*FCT(X)                                  DL24 720
      X=.7927539247172152D1                                             DL24 730
      Y=Y+.67216256409354789D-3*FCT(X)                                  DL24 740
      X=.61815351187367654D1                                            DL24 750
      Y=Y+.33693490584783036D-2*FCT(X)                                  DL24 760
      X=.46650837034671708D1                                            DL24 770
      Y=Y+.13226019405120157D-1*FCT(X)                                  DL24 780
      X=.33707742642089977D1                                            DL24 790
      Y=Y+.40732478151408646D-1*FCT(X)                                  DL24 800
      X=.22925620586321903D1                                            DL24 810
      Y=Y+.9816627262991889D-1*FCT(X)                                   DL24 820
      X=.14255975908036131D1                                            DL24 830
      Y=Y+.18332268897777802D0*FCT(X)                                   DL24 840
      X=.7660969055459366D0                                             DL24 850
      Y=Y+.25880670727286980D0*FCT(X)                                   DL24 860
      X=.31123914619848373D0                                            DL24 870
      Y=Y+.25877410751742390D0*FCT(X)                                   DL24 880
      X=.59019852181507977D-1                                           DL24 890
      Y=Y+.14281197333478185D0*FCT(X)                                   DL24 900
      RETURN                                                            DL24 910
      END                                                               DL24 920
