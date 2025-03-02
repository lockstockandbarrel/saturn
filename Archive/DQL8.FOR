C                                                                       DQL8  10
C     ..................................................................DQL8  20
C                                                                       DQL8  30
C        SUBROUTINE DQL8                                                DQL8  40
C                                                                       DQL8  50
C        PURPOSE                                                        DQL8  60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X), SUMMED OVER X           DQL8  70
C                               FROM 0 TO INFINITY).                    DQL8  80
C                                                                       DQL8  90
C        USAGE                                                          DQL8 100
C           CALL DQL8 (FCT,Y)                                           DQL8 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DQL8 120
C                                                                       DQL8 130
C        DESCRIPTION OF PARAMETERS                                      DQL8 140
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DQL8 150
C                    SUBPROGRAM USED.                                   DQL8 160
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DQL8 170
C                                                                       DQL8 180
C        REMARKS                                                        DQL8 190
C           NONE                                                        DQL8 200
C                                                                       DQL8 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DQL8 220
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DQL8 230
C           MUST BE FURNISHED BY THE USER.                              DQL8 240
C                                                                       DQL8 250
C        METHOD                                                         DQL8 260
C           EVALUATION IS DONE BY MEANS OF 8-POINT GAUSSIAN-LAGUERRE    DQL8 270
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY,               DQL8 280
C           WHENEVER FCT(X) IS A POLYNOMIAL UP TO DEGREE 15.            DQL8 290
C           FOR REFERENCE, SEE                                          DQL8 300
C           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF    DQL8 310
C           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED     DQL8 320
C           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT       DQL8 330
C           TR00.1100 (MARCH 1964), PP.24-25.                           DQL8 340
C                                                                       DQL8 350
C     ..................................................................DQL8 360
C                                                                       DQL8 370
      SUBROUTINE DQL8(FCT,Y)                                            DQL8 380
C                                                                       DQL8 390
C                                                                       DQL8 400
      DOUBLE PRECISION X,Y,FCT                                          DQL8 410
C                                                                       DQL8 420
      X=.22863131736889264D2                                            DQL8 430
      Y=.10480011748715104D-8*FCT(X)                                    DQL8 440
      X=.15740678641278005D2                                            DQL8 450
      Y=Y+.8485746716272532D-6*FCT(X)                                   DQL8 460
      X=.10758516010180995D2                                            DQL8 470
      Y=Y+.9076508773358213D-4*FCT(X)                                   DQL8 480
      X=.70459054023934657D1                                            DQL8 490
      Y=Y+.27945362352256725D-2*FCT(X)                                  DQL8 500
      X=.42667001702876588D1                                            DQL8 510
      Y=Y+.33343492261215652D-1*FCT(X)                                  DQL8 520
      X=.22510866298661307D1                                            DQL8 530
      Y=Y+.17579498663717181D0*FCT(X)                                   DQL8 540
      X=.9037017767993799D0                                             DQL8 550
      Y=Y+.41878678081434296D0*FCT(X)                                   DQL8 560
      X=.17027963230510100D0                                            DQL8 570
      Y=Y+.36918858934163753D0*FCT(X)                                   DQL8 580
      RETURN                                                            DQL8 590
      END                                                               DQL8 600
