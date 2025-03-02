C                                                                       QA6   10
C     ..................................................................QA6   20
C                                                                       QA6   30
C        SUBROUTINE QA6                                                 QA6   40
C                                                                       QA6   50
C        PURPOSE                                                        QA6   60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X)/SQRT(X), SUMMED OVER X   QA6   70
C                               FROM 0 TO INFINITY).                    QA6   80
C                                                                       QA6   90
C        USAGE                                                          QA6  100
C           CALL QA6 (FCT,Y)                                            QA6  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QA6  120
C                                                                       QA6  130
C        DESCRIPTION OF PARAMETERS                                      QA6  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QA6  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QA6  160
C                                                                       QA6  170
C        REMARKS                                                        QA6  180
C           NONE                                                        QA6  190
C                                                                       QA6  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QA6  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QA6  220
C           BY THE USER.                                                QA6  230
C                                                                       QA6  240
C        METHOD                                                         QA6  250
C           EVALUATION IS DONE BY MEANS OF 6-POINT GENERALIZED GAUSSIAN-QA6  260
C           LAGUERRE QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY,      QA6  270
C           WHENEVER FCT(X) IS A POLYNOMIAL UP TO DEGREE 11.            QA6  280
C           FOR REFERENCE, SEE                                          QA6  290
C           CONCUS/CASSATT/JAEHNIG/MELBY, TABLES FOR THE EVALUATION OF  QA6  300
C           INTEGRAL(X**BETA*EXP(-X)*F(X), SUMMED OVER X FROM 0 TO      QA6  310
C           INFINITY) BY GAUSS-LAGUERRE QUADRATURE, MTAC, VOL.17,       QA6  320
C           ISS.83 (1963), PP.245-256.                                  QA6  330
C                                                                       QA6  340
C     ..................................................................QA6  350
C                                                                       QA6  360
      SUBROUTINE QA6(FCT,Y)                                             QA6  370
C                                                                       QA6  380
C                                                                       QA6  390
      X=15.12996                                                        QA6  400
      Y=.5317103E-6*FCT(X)                                              QA6  410
      X=9.124248                                                        QA6  420
      Y=Y+.0001714737*FCT(X)                                            QA6  430
      X=5.196153                                                        QA6  440
      Y=Y+.007810781*FCT(X)                                             QA6  450
      X=2.552590                                                        QA6  460
      Y=Y+.1032160*FCT(X)                                               QA6  470
      X=.8983028                                                        QA6  480
      Y=Y+.5209846*FCT(X)                                               QA6  490
      X=.09874701                                                       QA6  500
      Y=Y+1.140270*FCT(X)                                               QA6  510
      RETURN                                                            QA6  520
      END                                                               QA6  530
