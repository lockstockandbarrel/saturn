C                                                                       QH4   10
C     ..................................................................QH4   20
C                                                                       QH4   30
C        SUBROUTINE QH4                                                 QH4   40
C                                                                       QH4   50
C        PURPOSE                                                        QH4   60
C           TO COMPUTE INTEGRAL(EXP(-X*X)*FCT(X), SUMMED OVER X FROM    QH4   70
C                               -INFINITY TO +INFINITY).                QH4   80
C                                                                       QH4   90
C        USAGE                                                          QH4  100
C           CALL QH4 (FCT,Y)                                            QH4  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QH4  120
C                                                                       QH4  130
C        DESCRIPTION OF PARAMETERS                                      QH4  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QH4  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QH4  160
C                                                                       QH4  170
C        REMARKS                                                        QH4  180
C           NONE                                                        QH4  190
C                                                                       QH4  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QH4  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QH4  220
C           BY THE USER.                                                QH4  230
C                                                                       QH4  240
C        METHOD                                                         QH4  250
C           EVALUATION IS DONE BY MEANS OF 4-POINT GAUSSIAN-HERMITE     QH4  260
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       QH4  270
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 7.                      QH4  280
C           FOR REFERENCE, SEE                                          QH4  290
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QH4  300
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.129-130 AND 343-346.   QH4  310
C                                                                       QH4  320
C     ..................................................................QH4  330
C                                                                       QH4  340
      SUBROUTINE QH4(FCT,Y)                                             QH4  350
C                                                                       QH4  360
C                                                                       QH4  370
      X=1.650680                                                        QH4  380
      Z=-X                                                              QH4  390
      Y=.08131284*(FCT(X)+FCT(Z))                                       QH4  400
      X=.5246476                                                        QH4  410
      Z=-X                                                              QH4  420
      Y=Y+.8049141*(FCT(X)+FCT(Z))                                      QH4  430
      RETURN                                                            QH4  440
      END                                                               QH4  450
