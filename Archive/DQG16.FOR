C                                                                       DG16  10
C     ..................................................................DG16  20
C                                                                       DG16  30
C        SUBROUTINE DQG16                                               DG16  40
C                                                                       DG16  50
C        PURPOSE                                                        DG16  60
C           TO COMPUTE INTEGRAL(FCT(X), SUMMED OVER X FROM XL TO XU)    DG16  70
C                                                                       DG16  80
C        USAGE                                                          DG16  90
C           CALL DQG16 (XL,XU,FCT,Y)                                    DG16 100
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DG16 110
C                                                                       DG16 120
C        DESCRIPTION OF PARAMETERS                                      DG16 130
C           XL     - DOUBLE PRECISION LOWER BOUND OF THE INTERVAL.      DG16 140
C           XU     - DOUBLE PRECISION UPPER BOUND OF THE INTERVAL.      DG16 150
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DG16 160
C                    SUBPROGRAM USED.                                   DG16 170
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DG16 180
C                                                                       DG16 190
C        REMARKS                                                        DG16 200
C           NONE                                                        DG16 210
C                                                                       DG16 220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DG16 230
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DG16 240
C           MUST BE FURNISHED BY THE USER.                              DG16 250
C                                                                       DG16 260
C        METHOD                                                         DG16 270
C           EVALUATION IS DONE BY MEANS OF 16-POINT GAUSS QUADRATURE    DG16 280
C           FORMULA, WHICH INTEGRATES POLYNOMIALS UP TO DEGREE 31       DG16 290
C           EXACTLY. FOR REFERENCE, SEE                                 DG16 300
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           DG16 310
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.100-111 AND 337-340.   DG16 320
C                                                                       DG16 330
C     ..................................................................DG16 340
C                                                                       DG16 350
      SUBROUTINE DQG16(XL,XU,FCT,Y)                                     DG16 360
C                                                                       DG16 370
C                                                                       DG16 380
      DOUBLE PRECISION XL,XU,Y,A,B,C,FCT                                DG16 390
C                                                                       DG16 400
      A=.5D0*(XU+XL)                                                    DG16 410
      B=XU-XL                                                           DG16 420
      C=.49470046749582497D0*B                                          DG16 430
      Y=.13576229705877047D-1*(FCT(A+C)+FCT(A-C))                       DG16 440
      C=.47228751153661629D0*B                                          DG16 450
      Y=Y+.31126761969323946D-1*(FCT(A+C)+FCT(A-C))                     DG16 460
      C=.43281560119391587D0*B                                          DG16 470
      Y=Y+.47579255841246392D-1*(FCT(A+C)+FCT(A-C))                     DG16 480
      C=.37770220417750152D0*B                                          DG16 490
      Y=Y+.62314485627766936D-1*(FCT(A+C)+FCT(A-C))                     DG16 500
      C=.30893812220132187D0*B                                          DG16 510
      Y=Y+.7479799440828837D-1*(FCT(A+C)+FCT(A-C))                      DG16 520
      C=.22900838882861369D0*B                                          DG16 530
      Y=Y+.8457825969750127D-1*(FCT(A+C)+FCT(A-C))                      DG16 540
      C=.14080177538962946D0*B                                          DG16 550
      Y=Y+.9130170752246179D-1*(FCT(A+C)+FCT(A-C))                      DG16 560
      C=.47506254918818720D-1*B                                         DG16 570
      Y=B*(Y+.9472530522753425D-1*(FCT(A+C)+FCT(A-C)))                  DG16 580
      RETURN                                                            DG16 590
      END                                                               DG16 600
