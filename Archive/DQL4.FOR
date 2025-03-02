C                                                                       DQL4  10
C     ..................................................................DQL4  20
C                                                                       DQL4  30
C        SUBROUTINE DQL4                                                DQL4  40
C                                                                       DQL4  50
C        PURPOSE                                                        DQL4  60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X), SUMMED OVER X           DQL4  70
C                               FROM 0 TO INFINITY).                    DQL4  80
C                                                                       DQL4  90
C        USAGE                                                          DQL4 100
C           CALL DQL4 (FCT,Y)                                           DQL4 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DQL4 120
C                                                                       DQL4 130
C        DESCRIPTION OF PARAMETERS                                      DQL4 140
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DQL4 150
C                    SUBPROGRAM USED.                                   DQL4 160
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DQL4 170
C                                                                       DQL4 180
C        REMARKS                                                        DQL4 190
C           NONE                                                        DQL4 200
C                                                                       DQL4 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DQL4 220
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DQL4 230
C           MUST BE FURNISHED BY THE USER.                              DQL4 240
C                                                                       DQL4 250
C        METHOD                                                         DQL4 260
C           EVALUATION IS DONE BY MEANS OF 4-POINT GAUSSIAN-LAGUERRE    DQL4 270
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY,               DQL4 280
C           WHENEVER FCT(X) IS A POLYNOMIAL UP TO DEGREE 7.             DQL4 290
C           FOR REFERENCE, SEE                                          DQL4 300
C           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF    DQL4 310
C           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED     DQL4 320
C           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT       DQL4 330
C           TR00.1100 (MARCH 1964), PP.24-25.                           DQL4 340
C                                                                       DQL4 350
C     ..................................................................DQL4 360
C                                                                       DQL4 370
      SUBROUTINE DQL4(FCT,Y)                                            DQL4 380
C                                                                       DQL4 390
C                                                                       DQL4 400
      DOUBLE PRECISION X,Y,FCT                                          DQL4 410
C                                                                       DQL4 420
      X=.9395070912301133D1                                             DQL4 430
      Y=.53929470556132745D-3*FCT(X)                                    DQL4 440
      X=.45366202969211280D1                                            DQL4 450
      Y=Y+.38887908515005384D-1*FCT(X)                                  DQL4 460
      X=.17457611011583466D1                                            DQL4 470
      Y=Y+.35741869243779969D0*FCT(X)                                   DQL4 480
      X=.32254768961939231D0                                            DQL4 490
      Y=Y+.60315410434163360D0*FCT(X)                                   DQL4 500
      RETURN                                                            DQL4 510
      END                                                               DQL4 520
