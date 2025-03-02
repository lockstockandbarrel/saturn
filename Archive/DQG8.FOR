C                                                                       DQG8  10
C     ..................................................................DQG8  20
C                                                                       DQG8  30
C        SUBROUTINE DQG8                                                DQG8  40
C                                                                       DQG8  50
C        PURPOSE                                                        DQG8  60
C           TO COMPUTE INTEGRAL(FCT(X), SUMMED OVER X FROM XL TO XU)    DQG8  70
C                                                                       DQG8  80
C        USAGE                                                          DQG8  90
C           CALL DQG8 (XL,XU,FCT,Y)                                     DQG8 100
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DQG8 110
C                                                                       DQG8 120
C        DESCRIPTION OF PARAMETERS                                      DQG8 130
C           XL     - DOUBLE PRECISION LOWER BOUND OF THE INTERVAL.      DQG8 140
C           XU     - DOUBLE PRECISION UPPER BOUND OF THE INTERVAL.      DQG8 150
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DQG8 160
C                    SUBPROGRAM USED.                                   DQG8 170
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DQG8 180
C                                                                       DQG8 190
C        REMARKS                                                        DQG8 200
C           NONE                                                        DQG8 210
C                                                                       DQG8 220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DQG8 230
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DQG8 240
C           MUST BE FURNISHED BY THE USER.                              DQG8 250
C                                                                       DQG8 260
C        METHOD                                                         DQG8 270
C           EVALUATION IS DONE BY MEANS OF 8-POINT GAUSS QUADRATURE     DQG8 280
C           FORMULA, WHICH INTEGRATES POLYNOMIALS UP TO DEGREE 15       DQG8 290
C           EXACTLY. FOR REFERENCE, SEE                                 DQG8 300
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           DQG8 310
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.100-111 AND 337-340.   DQG8 320
C                                                                       DQG8 330
C     ..................................................................DQG8 340
C                                                                       DQG8 350
      SUBROUTINE DQG8(XL,XU,FCT,Y)                                      DQG8 360
C                                                                       DQG8 370
C                                                                       DQG8 380
      DOUBLE PRECISION XL,XU,Y,A,B,C,FCT                                DQG8 390
C                                                                       DQG8 400
      A=.5D0*(XU+XL)                                                    DQG8 410
      B=XU-XL                                                           DQG8 420
      C=.48014492824876812D0*B                                          DQG8 430
      Y=.50614268145188130D-1*(FCT(A+C)+FCT(A-C))                       DQG8 440
      C=.39833323870681337D0*B                                          DQG8 450
      Y=Y+.11119051722668724D0*(FCT(A+C)+FCT(A-C))                      DQG8 460
      C=.26276620495816449D0*B                                          DQG8 470
      Y=Y+.15685332293894364D0*(FCT(A+C)+FCT(A-C))                      DQG8 480
      C=.9171732124782490D-1*B                                          DQG8 490
      Y=B*(Y+.18134189168918099D0*(FCT(A+C)+FCT(A-C)))                  DQG8 500
      RETURN                                                            DQG8 510
      END                                                               DQG8 520
