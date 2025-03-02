C                                                                       QH3   10
C     ..................................................................QH3   20
C                                                                       QH3   30
C        SUBROUTINE QH3                                                 QH3   40
C                                                                       QH3   50
C        PURPOSE                                                        QH3   60
C           TO COMPUTE INTEGRAL(EXP(-X*X)*FCT(X), SUMMED OVER X FROM    QH3   70
C                               -INFINITY TO +INFINITY).                QH3   80
C                                                                       QH3   90
C        USAGE                                                          QH3  100
C           CALL QH3 (FCT,Y)                                            QH3  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QH3  120
C                                                                       QH3  130
C        DESCRIPTION OF PARAMETERS                                      QH3  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QH3  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QH3  160
C                                                                       QH3  170
C        REMARKS                                                        QH3  180
C           NONE                                                        QH3  190
C                                                                       QH3  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QH3  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QH3  220
C           BY THE USER.                                                QH3  230
C                                                                       QH3  240
C        METHOD                                                         QH3  250
C           EVALUATION IS DONE BY MEANS OF 3-POINT GAUSSIAN-HERMITE     QH3  260
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       QH3  270
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 5.                      QH3  280
C           FOR REFERENCE, SEE                                          QH3  290
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QH3  300
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.129-130 AND 343-346.   QH3  310
C                                                                       QH3  320
C     ..................................................................QH3  330
C                                                                       QH3  340
      SUBROUTINE QH3(FCT,Y)                                             QH3  350
C                                                                       QH3  360
C                                                                       QH3  370
      X=1.224745                                                        QH3  380
      Z=-X                                                              QH3  390
      Y=.2954090*(FCT(X)+FCT(Z))                                        QH3  400
      X=0.                                                              QH3  410
      Y=Y+1.181636*FCT(X)                                               QH3  420
      RETURN                                                            QH3  430
      END                                                               QH3  440
