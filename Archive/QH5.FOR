C                                                                       QH5   10
C     ..................................................................QH5   20
C                                                                       QH5   30
C        SUBROUTINE QH5                                                 QH5   40
C                                                                       QH5   50
C        PURPOSE                                                        QH5   60
C           TO COMPUTE INTEGRAL(EXP(-X*X)*FCT(X), SUMMED OVER X FROM    QH5   70
C                               -INFINITY TO +INFINITY).                QH5   80
C                                                                       QH5   90
C        USAGE                                                          QH5  100
C           CALL QH5 (FCT,Y)                                            QH5  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QH5  120
C                                                                       QH5  130
C        DESCRIPTION OF PARAMETERS                                      QH5  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QH5  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QH5  160
C                                                                       QH5  170
C        REMARKS                                                        QH5  180
C           NONE                                                        QH5  190
C                                                                       QH5  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QH5  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QH5  220
C           BY THE USER.                                                QH5  230
C                                                                       QH5  240
C        METHOD                                                         QH5  250
C           EVALUATION IS DONE BY MEANS OF 5-POINT GAUSSIAN-HERMITE     QH5  260
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       QH5  270
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 9.                      QH5  280
C           FOR REFERENCE, SEE                                          QH5  290
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QH5  300
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.129-130 AND 343-346.   QH5  310
C                                                                       QH5  320
C     ..................................................................QH5  330
C                                                                       QH5  340
      SUBROUTINE QH5(FCT,Y)                                             QH5  350
C                                                                       QH5  360
C                                                                       QH5  370
      X=2.020183                                                        QH5  380
      Z=-X                                                              QH5  390
      Y=.01995324*(FCT(X)+FCT(Z))                                       QH5  400
      X=.9585725                                                        QH5  410
      Z=-X                                                              QH5  420
      Y=Y+.3936193*(FCT(X)+FCT(Z))                                      QH5  430
      X=0.                                                              QH5  440
      Y=Y+.9453087*FCT(X)                                               QH5  450
      RETURN                                                            QH5  460
      END                                                               QH5  470
