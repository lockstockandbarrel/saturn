C                                                                       QH9   10
C     ..................................................................QH9   20
C                                                                       QH9   30
C        SUBROUTINE QH9                                                 QH9   40
C                                                                       QH9   50
C        PURPOSE                                                        QH9   60
C           TO COMPUTE INTEGRAL(EXP(-X*X)*FCT(X), SUMMED OVER X FROM    QH9   70
C                               -INFINITY TO +INFINITY).                QH9   80
C                                                                       QH9   90
C        USAGE                                                          QH9  100
C           CALL QH9 (FCT,Y)                                            QH9  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QH9  120
C                                                                       QH9  130
C        DESCRIPTION OF PARAMETERS                                      QH9  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QH9  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QH9  160
C                                                                       QH9  170
C        REMARKS                                                        QH9  180
C           NONE                                                        QH9  190
C                                                                       QH9  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QH9  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QH9  220
C           BY THE USER.                                                QH9  230
C                                                                       QH9  240
C        METHOD                                                         QH9  250
C           EVALUATION IS DONE BY MEANS OF 9-POINT GAUSSIAN-HERMITE     QH9  260
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       QH9  270
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 17.                     QH9  280
C           FOR REFERENCE, SEE                                          QH9  290
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QH9  300
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.129-130 AND 343-346.   QH9  310
C                                                                       QH9  320
C     ..................................................................QH9  330
C                                                                       QH9  340
      SUBROUTINE QH9(FCT,Y)                                             QH9  350
C                                                                       QH9  360
C                                                                       QH9  370
      X=3.190993                                                        QH9  380
      Z=-X                                                              QH9  390
      Y=.3960698E-4*(FCT(X)+FCT(Z))                                     QH9  400
      X=2.266581                                                        QH9  410
      Z=-X                                                              QH9  420
      Y=Y+.004943624*(FCT(X)+FCT(Z))                                    QH9  430
      X=1.468553                                                        QH9  440
      Z=-X                                                              QH9  450
      Y=Y+.08847453*(FCT(X)+FCT(Z))                                     QH9  460
      X=.7235510                                                        QH9  470
      Z=-X                                                              QH9  480
      Y=Y+.4326516*(FCT(X)+FCT(Z))                                      QH9  490
      X=0.                                                              QH9  500
      Y=Y+.7202352*FCT(X)                                               QH9  510
      RETURN                                                            QH9  520
      END                                                               QH9  530
