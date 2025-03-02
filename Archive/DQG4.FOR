C                                                                       DQG4  10
C     ..................................................................DQG4  20
C                                                                       DQG4  30
C        SUBROUTINE DQG4                                                DQG4  40
C                                                                       DQG4  50
C        PURPOSE                                                        DQG4  60
C           TO COMPUTE INTEGRAL(FCT(X), SUMMED OVER X FROM XL TO XU)    DQG4  70
C                                                                       DQG4  80
C        USAGE                                                          DQG4  90
C           CALL DQG4 (XL,XU,FCT,Y)                                     DQG4 100
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DQG4 110
C                                                                       DQG4 120
C        DESCRIPTION OF PARAMETERS                                      DQG4 130
C           XL     - DOUBLE PRECISION LOWER BOUND OF THE INTERVAL.      DQG4 140
C           XU     - DOUBLE PRECISION UPPER BOUND OF THE INTERVAL.      DQG4 150
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DQG4 160
C                    SUBPROGRAM USED.                                   DQG4 170
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DQG4 180
C                                                                       DQG4 190
C        REMARKS                                                        DQG4 200
C           NONE                                                        DQG4 210
C                                                                       DQG4 220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DQG4 230
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DQG4 240
C           MUST BE FURNISHED BY THE USER.                              DQG4 250
C                                                                       DQG4 260
C        METHOD                                                         DQG4 270
C           EVALUATION IS DONE BY MEANS OF 4-POINT GAUSS QUADRATURE     DQG4 280
C           FORMULA, WHICH INTEGRATES POLYNOMIALS UP TO DEGREE 7        DQG4 290
C           EXACTLY. FOR REFERENCE, SEE                                 DQG4 300
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           DQG4 310
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.100-111 AND 337-340.   DQG4 320
C                                                                       DQG4 330
C     ..................................................................DQG4 340
C                                                                       DQG4 350
      SUBROUTINE DQG4(XL,XU,FCT,Y)                                      DQG4 360
C                                                                       DQG4 370
C                                                                       DQG4 380
      DOUBLE PRECISION XL,XU,Y,A,B,C,FCT                                DQG4 390
C                                                                       DQG4 400
      A=.5D0*(XU+XL)                                                    DQG4 410
      B=XU-XL                                                           DQG4 420
      C=.43056815579702629D0*B                                          DQG4 430
      Y=.17392742256872693D0*(FCT(A+C)+FCT(A-C))                        DQG4 440
      C=.16999052179242813D0*B                                          DQG4 450
      Y=B*(Y+.32607257743127307D0*(FCT(A+C)+FCT(A-C)))                  DQG4 460
      RETURN                                                            DQG4 470
      END                                                               DQG4 480
