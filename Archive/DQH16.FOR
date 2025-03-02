C                                                                       DH16  10
C     ..................................................................DH16  20
C                                                                       DH16  30
C        SUBROUTINE DQH16                                               DH16  40
C                                                                       DH16  50
C        PURPOSE                                                        DH16  60
C           TO COMPUTE INTEGRAL(EXP(-X*X)*FCT(X), SUMMED OVER X FROM    DH16  70
C                               -INFINITY TO +INFINITY).                DH16  80
C                                                                       DH16  90
C        USAGE                                                          DH16 100
C           CALL DQH16 (FCT,Y)                                          DH16 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DH16 120
C                                                                       DH16 130
C        DESCRIPTION OF PARAMETERS                                      DH16 140
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DH16 150
C                    SUBPROGRAM USED.                                   DH16 160
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DH16 170
C                                                                       DH16 180
C        REMARKS                                                        DH16 190
C           NONE                                                        DH16 200
C                                                                       DH16 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DH16 220
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DH16 230
C           MUST BE FURNISHED BY THE USER.                              DH16 240
C                                                                       DH16 250
C        METHOD                                                         DH16 260
C           EVALUATION IS DONE BY MEANS OF 16-POINT GAUSSIAN-HERMITE    DH16 270
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       DH16 280
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 31.                     DH16 290
C           FOR REFERENCE, SEE                                          DH16 300
C           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF    DH16 310
C           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED     DH16 320
C           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT       DH16 330
C           TR00.1100 (MARCH 1964), PP.213-214.                         DH16 340
C                                                                       DH16 350
C     ..................................................................DH16 360
C                                                                       DH16 370
      SUBROUTINE DQH16(FCT,Y)                                           DH16 380
C                                                                       DH16 390
C                                                                       DH16 400
      DOUBLE PRECISION X,Y,Z,FCT                                        DH16 410
C                                                                       DH16 420
      X=.46887389393058184D1                                            DH16 430
      Z=-X                                                              DH16 440
      Y=.26548074740111822D-9*(FCT(X)+FCT(Z))                           DH16 450
      X=.38694479048601227D1                                            DH16 460
      Z=-X                                                              DH16 470
      Y=Y+.23209808448652107D-6*(FCT(X)+FCT(Z))                         DH16 480
      X=.31769991619799560D1                                            DH16 490
      Z=-X                                                              DH16 500
      Y=Y+.27118600925378815D-4*(FCT(X)+FCT(Z))                         DH16 510
      X=.25462021578474814D1                                            DH16 520
      Z=-X                                                              DH16 530
      Y=Y+.9322840086241805D-3*(FCT(X)+FCT(Z))                          DH16 540
      X=.19517879909162540D1                                            DH16 550
      Z=-X                                                              DH16 560
      Y=Y+.12880311535509974D-1*(FCT(X)+FCT(Z))                         DH16 570
      X=.13802585391988808D1                                            DH16 580
      Z=-X                                                              DH16 590
      Y=Y+.8381004139898583D-1*(FCT(X)+FCT(Z))                          DH16 600
      X=.8229514491446559D0                                             DH16 610
      Z=-X                                                              DH16 620
      Y=Y+.28064745852853368D0*(FCT(X)+FCT(Z))                          DH16 630
      X=.27348104613815245D0                                            DH16 640
      Z=-X                                                              DH16 650
      Y=Y+.50792947901661374D0*(FCT(X)+FCT(Z))                          DH16 660
      RETURN                                                            DH16 670
      END                                                               DH16 680
