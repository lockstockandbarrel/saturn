C                                                                       QH8   10
C     ..................................................................QH8   20
C                                                                       QH8   30
C        SUBROUTINE QH8                                                 QH8   40
C                                                                       QH8   50
C        PURPOSE                                                        QH8   60
C           TO COMPUTE INTEGRAL(EXP(-X*X)*FCT(X), SUMMED OVER X FROM    QH8   70
C                               -INFINITY TO +INFINITY).                QH8   80
C                                                                       QH8   90
C        USAGE                                                          QH8  100
C           CALL QH8 (FCT,Y)                                            QH8  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QH8  120
C                                                                       QH8  130
C        DESCRIPTION OF PARAMETERS                                      QH8  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QH8  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QH8  160
C                                                                       QH8  170
C        REMARKS                                                        QH8  180
C           NONE                                                        QH8  190
C                                                                       QH8  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QH8  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QH8  220
C           BY THE USER.                                                QH8  230
C                                                                       QH8  240
C        METHOD                                                         QH8  250
C           EVALUATION IS DONE BY MEANS OF 8-POINT GAUSSIAN-HERMITE     QH8  260
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       QH8  270
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 15.                     QH8  280
C           FOR REFERENCE, SEE                                          QH8  290
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QH8  300
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.129-130 AND 343-346.   QH8  310
C                                                                       QH8  320
C     ..................................................................QH8  330
C                                                                       QH8  340
      SUBROUTINE QH8(FCT,Y)                                             QH8  350
C                                                                       QH8  360
C                                                                       QH8  370
      X=2.930637                                                        QH8  380
      Z=-X                                                              QH8  390
      Y=.0001996041*(FCT(X)+FCT(Z))                                     QH8  400
      X=1.981657                                                        QH8  410
      Z=-X                                                              QH8  420
      Y=Y+.01707798*(FCT(X)+FCT(Z))                                     QH8  430
      X=1.157194                                                        QH8  440
      Z=-X                                                              QH8  450
      Y=Y+.2078023*(FCT(X)+FCT(Z))                                      QH8  460
      X=.3811870                                                        QH8  470
      Z=-X                                                              QH8  480
      Y=Y+.6611470*(FCT(X)+FCT(Z))                                      QH8  490
      RETURN                                                            QH8  500
      END                                                               QH8  510
