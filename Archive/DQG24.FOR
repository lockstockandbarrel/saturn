C                                                                       DG24  10
C     ..................................................................DG24  20
C                                                                       DG24  30
C        SUBROUTINE DQG24                                               DG24  40
C                                                                       DG24  50
C        PURPOSE                                                        DG24  60
C           TO COMPUTE INTEGRAL(FCT(X), SUMMED OVER X FROM XL TO XU)    DG24  70
C                                                                       DG24  80
C        USAGE                                                          DG24  90
C           CALL DQG24 (XL,XU,FCT,Y)                                    DG24 100
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DG24 110
C                                                                       DG24 120
C        DESCRIPTION OF PARAMETERS                                      DG24 130
C           XL     - DOUBLE PRECISION LOWER BOUND OF THE INTERVAL.      DG24 140
C           XU     - DOUBLE PRECISION UPPER BOUND OF THE INTERVAL.      DG24 150
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DG24 160
C                    SUBPROGRAM USED.                                   DG24 170
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DG24 180
C                                                                       DG24 190
C        REMARKS                                                        DG24 200
C           NONE                                                        DG24 210
C                                                                       DG24 220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DG24 230
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DG24 240
C           MUST BE FURNISHED BY THE USER.                              DG24 250
C                                                                       DG24 260
C        METHOD                                                         DG24 270
C           EVALUATION IS DONE BY MEANS OF 24-POINT GAUSS QUADRATURE    DG24 280
C           FORMULA, WHICH INTEGRATES POLYNOMIALS UP TO DEGREE 47       DG24 290
C           EXACTLY. FOR REFERENCE, SEE                                 DG24 300
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           DG24 310
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.100-111 AND 337-340.   DG24 320
C                                                                       DG24 330
C     ..................................................................DG24 340
C                                                                       DG24 350
      SUBROUTINE DQG24(XL,XU,FCT,Y)                                     DG24 360
C                                                                       DG24 370
C                                                                       DG24 380
      DOUBLE PRECISION XL,XU,Y,A,B,C,FCT                                DG24 390
C                                                                       DG24 400
      A=.5D0*(XU+XL)                                                    DG24 410
      B=XU-XL                                                           DG24 420
      C=.49759360999851068D0*B                                          DG24 430
      Y=.61706148999935998D-2*(FCT(A+C)+FCT(A-C))                       DG24 440
      C=.48736427798565475D0*B                                          DG24 450
      Y=Y+.14265694314466832D-1*(FCT(A+C)+FCT(A-C))                     DG24 460
      C=.46913727600136638D0*B                                          DG24 470
      Y=Y+.22138719408709903D-1*(FCT(A+C)+FCT(A-C))                     DG24 480
      C=.44320776350220052D0*B                                          DG24 490
      Y=Y+.29649292457718390D-1*(FCT(A+C)+FCT(A-C))                     DG24 500
      C=.41000099298695146D0*B                                          DG24 510
      Y=Y+.36673240705540153D-1*(FCT(A+C)+FCT(A-C))                     DG24 520
      C=.37006209578927718D0*B                                          DG24 530
      Y=Y+.43095080765976638D-1*(FCT(A+C)+FCT(A-C))                     DG24 540
      C=.32404682596848778D0*B                                          DG24 550
      Y=Y+.48809326052056944D-1*(FCT(A+C)+FCT(A-C))                     DG24 560
      C=.27271073569441977D0*B                                          DG24 570
      Y=Y+.53722135057982817D-1*(FCT(A+C)+FCT(A-C))                     DG24 580
      C=.21689675381302257D0*B                                          DG24 590
      Y=Y+.57752834026862801D-1*(FCT(A+C)+FCT(A-C))                     DG24 600
      C=.15752133984808169D0*B                                          DG24 610
      Y=Y+.60835236463901696D-1*(FCT(A+C)+FCT(A-C))                     DG24 620
      C=.9555943373680815D-1*B                                          DG24 630
      Y=Y+.62918728173414148D-1*(FCT(A+C)+FCT(A-C))                     DG24 640
      C=.32028446431302813D-1*B                                         DG24 650
      Y=B*(Y+.63969097673376078D-1*(FCT(A+C)+FCT(A-C)))                 DG24 660
      RETURN                                                            DG24 670
      END                                                               DG24 680
