C                                                                       QH2   10
C     ..................................................................QH2   20
C                                                                       QH2   30
C        SUBROUTINE QH2                                                 QH2   40
C                                                                       QH2   50
C        PURPOSE                                                        QH2   60
C           TO COMPUTE INTEGRAL(EXP(-X*X)*FCT(X), SUMMED OVER X FROM    QH2   70
C                               -INFINITY TO +INFINITY).                QH2   80
C                                                                       QH2   90
C        USAGE                                                          QH2  100
C           CALL QH2 (FCT,Y)                                            QH2  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QH2  120
C                                                                       QH2  130
C        DESCRIPTION OF PARAMETERS                                      QH2  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QH2  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QH2  160
C                                                                       QH2  170
C        REMARKS                                                        QH2  180
C           NONE                                                        QH2  190
C                                                                       QH2  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QH2  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QH2  220
C           BY THE USER.                                                QH2  230
C                                                                       QH2  240
C        METHOD                                                         QH2  250
C           EVALUATION IS DONE BY MEANS OF 2-POINT GAUSSIAN-HERMITE     QH2  260
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       QH2  270
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 3.                      QH2  280
C           FOR REFERENCE, SEE                                          QH2  290
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QH2  300
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.129-130 AND 343-346.   QH2  310
C                                                                       QH2  320
C     ..................................................................QH2  330
C                                                                       QH2  340
      SUBROUTINE QH2(FCT,Y)                                             QH2  350
C                                                                       QH2  360
C                                                                       QH2  370
      X=.7071068                                                        QH2  380
      Z=-X                                                              QH2  390
      Y=.8862269*(FCT(X)+FCT(Z))                                        QH2  400
      RETURN                                                            QH2  410
      END                                                               QH2  420
