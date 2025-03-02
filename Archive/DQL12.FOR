C                                                                       DL12  10
C     ..................................................................DL12  20
C                                                                       DL12  30
C        SUBROUTINE DQL12                                               DL12  40
C                                                                       DL12  50
C        PURPOSE                                                        DL12  60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X), SUMMED OVER X           DL12  70
C                               FROM 0 TO INFINITY).                    DL12  80
C                                                                       DL12  90
C        USAGE                                                          DL12 100
C           CALL DQL12 (FCT,Y)                                          DL12 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DL12 120
C                                                                       DL12 130
C        DESCRIPTION OF PARAMETERS                                      DL12 140
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DL12 150
C                    SUBPROGRAM USED.                                   DL12 160
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DL12 170
C                                                                       DL12 180
C        REMARKS                                                        DL12 190
C           NONE                                                        DL12 200
C                                                                       DL12 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DL12 220
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DL12 230
C           MUST BE FURNISHED BY THE USER.                              DL12 240
C                                                                       DL12 250
C        METHOD                                                         DL12 260
C           EVALUATION IS DONE BY MEANS OF 12-POINT GAUSSIAN-LAGUERRE   DL12 270
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY,               DL12 280
C           WHENEVER FCT(X) IS A POLYNOMIAL UP TO DEGREE 23.            DL12 290
C           FOR REFERENCE, SEE                                          DL12 300
C           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF    DL12 310
C           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED     DL12 320
C           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT       DL12 330
C           TR00.1100 (MARCH 1964), PP.24-25.                           DL12 340
C                                                                       DL12 350
C     ..................................................................DL12 360
C                                                                       DL12 370
      SUBROUTINE DQL12(FCT,Y)                                           DL12 380
C                                                                       DL12 390
C                                                                       DL12 400
      DOUBLE PRECISION X,Y,FCT                                          DL12 410
C                                                                       DL12 420
      X=.37099121044466920D2                                            DL12 430
      Y=.8148077467426242D-15*FCT(X)                                    DL12 440
      X=.28487967250984000D2                                            DL12 450
      Y=Y+.30616016350350208D-11*FCT(X)                                 DL12 460
      X=.22151090379397006D2                                            DL12 470
      Y=Y+.13423910305150041D-8*FCT(X)                                  DL12 480
      X=.17116855187462256D2                                            DL12 490
      Y=Y+.16684938765409103D-6*FCT(X)                                  DL12 500
      X=.13006054993306348D2                                            DL12 510
      Y=Y+.8365055856819799D-5*FCT(X)                                   DL12 520
      X=.9621316842456867D1                                             DL12 530
      Y=Y+.20323159266299939D-3*FCT(X)                                  DL12 540
      X=.68445254531151773D1                                            DL12 550
      Y=Y+.26639735418653159D-2*FCT(X)                                  DL12 560
      X=.45992276394183485D1                                            DL12 570
      Y=Y+.20102381154634097D-1*FCT(X)                                  DL12 580
      X=.28337513377435072D1                                            DL12 590
      Y=Y+.9044922221168093D-1*FCT(X)                                   DL12 600
      X=.15126102697764188D1                                            DL12 610
      Y=Y+.24408201131987756D0*FCT(X)                                   DL12 620
      X=.61175748451513067D0                                            DL12 630
      Y=Y+.37775927587313798D0*FCT(X)                                   DL12 640
      X=.11572211735802068D0                                            DL12 650
      Y=Y+.26473137105544319D0*FCT(X)                                   DL12 660
      RETURN                                                            DL12 670
      END                                                               DL12 680
