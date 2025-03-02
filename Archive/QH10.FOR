C                                                                       QH10  10
C     ..................................................................QH10  20
C                                                                       QH10  30
C        SUBROUTINE QH10                                                QH10  40
C                                                                       QH10  50
C        PURPOSE                                                        QH10  60
C           TO COMPUTE INTEGRAL(EXP(-X*X)*FCT(X), SUMMED OVER X FROM    QH10  70
C                               -INFINITY TO +INFINITY).                QH10  80
C                                                                       QH10  90
C        USAGE                                                          QH10 100
C           CALL QH10(FCT,Y)                                            QH10 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QH10 120
C                                                                       QH10 130
C        DESCRIPTION OF PARAMETERS                                      QH10 140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QH10 150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QH10 160
C                                                                       QH10 170
C        REMARKS                                                        QH10 180
C           NONE                                                        QH10 190
C                                                                       QH10 200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QH10 210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QH10 220
C           BY THE USER.                                                QH10 230
C                                                                       QH10 240
C        METHOD                                                         QH10 250
C           EVALUATION IS DONE BY MEANS OF 10-POINT GAUSSIAN-HERMITE    QH10 260
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       QH10 270
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 19.                     QH10 280
C           FOR REFERENCE, SEE                                          QH10 290
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QH10 300
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.129-130 AND 343-346.   QH10 310
C                                                                       QH10 320
C     ..................................................................QH10 330
C                                                                       QH10 340
      SUBROUTINE QH10(FCT,Y)                                            QH10 350
C                                                                       QH10 360
C                                                                       QH10 370
      X=3.436159                                                        QH10 380
      Z=-X                                                              QH10 390
      Y=.7640433E-5*(FCT(X)+FCT(Z))                                     QH10 400
      X=2.532732                                                        QH10 410
      Z=-X                                                              QH10 420
      Y=Y+.001343646*(FCT(X)+FCT(Z))                                    QH10 430
      X=1.756684                                                        QH10 440
      Z=-X                                                              QH10 450
      Y=Y+.03387439*(FCT(X)+FCT(Z))                                     QH10 460
      X=1.036611                                                        QH10 470
      Z=-X                                                              QH10 480
      Y=Y+.2401386*(FCT(X)+FCT(Z))                                      QH10 490
      X=.3429013                                                        QH10 500
      Z=-X                                                              QH10 510
      Y=Y+.6108626*(FCT(X)+FCT(Z))                                      QH10 520
      RETURN                                                            QH10 530
      END                                                               QH10 540
