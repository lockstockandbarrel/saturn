C                                                                       DG12  10
C     ..................................................................DG12  20
C                                                                       DG12  30
C        SUBROUTINE DQG12                                               DG12  40
C                                                                       DG12  50
C        PURPOSE                                                        DG12  60
C           TO COMPUTE INTEGRAL(FCT(X), SUMMED OVER X FROM XL TO XU)    DG12  70
C                                                                       DG12  80
C        USAGE                                                          DG12  90
C           CALL DQG12 (XL,XU,FCT,Y)                                    DG12 100
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DG12 110
C                                                                       DG12 120
C        DESCRIPTION OF PARAMETERS                                      DG12 130
C           XL     - DOUBLE PRECISION LOWER BOUND OF THE INTERVAL.      DG12 140
C           XU     - DOUBLE PRECISION UPPER BOUND OF THE INTERVAL.      DG12 150
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DG12 160
C                    SUBPROGRAM USED.                                   DG12 170
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DG12 180
C                                                                       DG12 190
C        REMARKS                                                        DG12 200
C           NONE                                                        DG12 210
C                                                                       DG12 220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DG12 230
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DG12 240
C           MUST BE FURNISHED BY THE USER.                              DG12 250
C                                                                       DG12 260
C        METHOD                                                         DG12 270
C           EVALUATION IS DONE BY MEANS OF 12-POINT GAUSS QUADRATURE    DG12 280
C           FORMULA, WHICH INTEGRATES POLYNOMIALS UP TO DEGREE 23       DG12 290
C           EXACTLY. FOR REFERENCE, SEE                                 DG12 300
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           DG12 310
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.100-111 AND 337-340.   DG12 320
C                                                                       DG12 330
C     ..................................................................DG12 340
C                                                                       DG12 350
      SUBROUTINE DQG12(XL,XU,FCT,Y)                                     DG12 360
C                                                                       DG12 370
C                                                                       DG12 380
      DOUBLE PRECISION XL,XU,Y,A,B,C,FCT                                DG12 390
C                                                                       DG12 400
      A=.5D0*(XU+XL)                                                    DG12 410
      B=XU-XL                                                           DG12 420
      C=.49078031712335963D0*B                                          DG12 430
      Y=.23587668193255914D-1*(FCT(A+C)+FCT(A-C))                       DG12 440
      C=.45205862818523743D0*B                                          DG12 450
      Y=Y+.53469662997659215D-1*(FCT(A+C)+FCT(A-C))                     DG12 460
      C=.38495133709715234D0*B                                          DG12 470
      Y=Y+.8003916427167311D-1*(FCT(A+C)+FCT(A-C))                      DG12 480
      C=.29365897714330872D0*B                                          DG12 490
      Y=Y+.10158371336153296D0*(FCT(A+C)+FCT(A-C))                      DG12 500
      C=.18391574949909010D0*B                                          DG12 510
      Y=Y+.11674626826917740D0*(FCT(A+C)+FCT(A-C))                      DG12 520
      C=.62616704255734458D-1*B                                         DG12 530
      Y=B*(Y+.12457352290670139D0*(FCT(A+C)+FCT(A-C)))                  DG12 540
      RETURN                                                            DG12 550
      END                                                               DG12 560
