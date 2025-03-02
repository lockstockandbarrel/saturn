C                                                                       QH7   10
C     ..................................................................QH7   20
C                                                                       QH7   30
C        SUBROUTINE QH7                                                 QH7   40
C                                                                       QH7   50
C        PURPOSE                                                        QH7   60
C           TO COMPUTE INTEGRAL(EXP(-X*X)*FCT(X), SUMMED OVER X FROM    QH7   70
C                               -INFINITY TO +INFINITY).                QH7   80
C                                                                       QH7   90
C        USAGE                                                          QH7  100
C           CALL QH7 (FCT,Y)                                            QH7  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QH7  120
C                                                                       QH7  130
C        DESCRIPTION OF PARAMETERS                                      QH7  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QH7  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QH7  160
C                                                                       QH7  170
C        REMARKS                                                        QH7  180
C           NONE                                                        QH7  190
C                                                                       QH7  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QH7  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QH7  220
C           BY THE USER.                                                QH7  230
C                                                                       QH7  240
C        METHOD                                                         QH7  250
C           EVALUATION IS DONE BY MEANS OF 7-POINT GAUSSIAN-HERMITE     QH7  260
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       QH7  270
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 13.                     QH7  280
C           FOR REFERENCE, SEE                                          QH7  290
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QH7  300
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.129-130 AND 343-346.   QH7  310
C                                                                       QH7  320
C     ..................................................................QH7  330
C                                                                       QH7  340
      SUBROUTINE QH7(FCT,Y)                                             QH7  350
C                                                                       QH7  360
C                                                                       QH7  370
      X=2.651961                                                        QH7  380
      Z=-X                                                              QH7  390
      Y=.0009717812*(FCT(X)+FCT(Z))                                     QH7  400
      X=1.673552                                                        QH7  410
      Z=-X                                                              QH7  420
      Y=Y+.05451558*(FCT(X)+FCT(Z))                                     QH7  430
      X=.8162879                                                        QH7  440
      Z=-X                                                              QH7  450
      Y=Y+.4256073*(FCT(X)+FCT(Z))                                      QH7  460
      X=0.                                                              QH7  470
      Y=Y+.8102646*FCT(X)                                               QH7  480
      RETURN                                                            QH7  490
      END                                                               QH7  500
