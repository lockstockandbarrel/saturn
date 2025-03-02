C                                                                       QA3   10
C     ..................................................................QA3   20
C                                                                       QA3   30
C        SUBROUTINE QA3                                                 QA3   40
C                                                                       QA3   50
C        PURPOSE                                                        QA3   60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X)/SQRT(X), SUMMED OVER X   QA3   70
C                               FROM 0 TO INFINITY).                    QA3   80
C                                                                       QA3   90
C        USAGE                                                          QA3  100
C           CALL QA3 (FCT,Y)                                            QA3  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QA3  120
C                                                                       QA3  130
C        DESCRIPTION OF PARAMETERS                                      QA3  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QA3  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QA3  160
C                                                                       QA3  170
C        REMARKS                                                        QA3  180
C           NONE                                                        QA3  190
C                                                                       QA3  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QA3  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QA3  220
C           BY THE USER.                                                QA3  230
C                                                                       QA3  240
C        METHOD                                                         QA3  250
C           EVALUATION IS DONE BY MEANS OF 3-POINT GENERALIZED GAUSSIAN-QA3  260
C           LAGUERRE QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY,      QA3  270
C           WHENEVER FCT(X) IS A POLYNOMIAL UP TO DEGREE 5.             QA3  280
C           FOR REFERENCE, SEE                                          QA3  290
C           CONCUS/CASSATT/JAEHNIG/MELBY, TABLES FOR THE EVALUATION OF  QA3  300
C           INTEGRAL(X**BETA*EXP(-X)*F(X), SUMMED OVER X FROM 0 TO      QA3  310
C           INFINITY) BY GAUSS-LAGUERRE QUADRATURE, MTAC, VOL.17,       QA3  320
C           ISS.83 (1963), PP.245-256.                                  QA3  330
C                                                                       QA3  340
C     ..................................................................QA3  350
C                                                                       QA3  360
      SUBROUTINE QA3(FCT,Y)                                             QA3  370
C                                                                       QA3  380
C                                                                       QA3  390
      X=5.525344                                                        QA3  400
      Y=.009060020*FCT(X)                                               QA3  410
      X=1.784493                                                        QA3  420
      Y=Y+.3141346*FCT(X)                                               QA3  430
      X=.1901635                                                        QA3  440
      Y=Y+1.449259*FCT(X)                                               QA3  450
      RETURN                                                            QA3  460
      END                                                               QA3  470
