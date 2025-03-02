C                                                                       QH6   10
C     ..................................................................QH6   20
C                                                                       QH6   30
C        SUBROUTINE QH6                                                 QH6   40
C                                                                       QH6   50
C        PURPOSE                                                        QH6   60
C           TO COMPUTE INTEGRAL(EXP(-X*X)*FCT(X), SUMMED OVER X FROM    QH6   70
C                               -INFINITY TO +INFINITY).                QH6   80
C                                                                       QH6   90
C        USAGE                                                          QH6  100
C           CALL QH6 (FCT,Y)                                            QH6  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QH6  120
C                                                                       QH6  130
C        DESCRIPTION OF PARAMETERS                                      QH6  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QH6  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QH6  160
C                                                                       QH6  170
C        REMARKS                                                        QH6  180
C           NONE                                                        QH6  190
C                                                                       QH6  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QH6  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QH6  220
C           BY THE USER.                                                QH6  230
C                                                                       QH6  240
C        METHOD                                                         QH6  250
C           EVALUATION IS DONE BY MEANS OF 6-POINT GAUSSIAN-HERMITE     QH6  260
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       QH6  270
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 11.                     QH6  280
C           FOR REFERENCE, SEE                                          QH6  290
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QH6  300
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.129-130 AND 343-346.   QH6  310
C                                                                       QH6  320
C     ..................................................................QH6  330
C                                                                       QH6  340
      SUBROUTINE QH6(FCT,Y)                                             QH6  350
C                                                                       QH6  360
C                                                                       QH6  370
      X=2.350605                                                        QH6  380
      Z=-X                                                              QH6  390
      Y=.004530010*(FCT(X)+FCT(Z))                                      QH6  400
      X=1.335849                                                        QH6  410
      Z=-X                                                              QH6  420
      Y=Y+.1570673*(FCT(X)+FCT(Z))                                      QH6  430
      X=.4360774                                                        QH6  440
      Z=-X                                                              QH6  450
      Y=Y+.7246296*(FCT(X)+FCT(Z))                                      QH6  460
      RETURN                                                            QH6  470
      END                                                               QH6  480
