C                                                                       QA5   10
C     ..................................................................QA5   20
C                                                                       QA5   30
C        SUBROUTINE QA5                                                 QA5   40
C                                                                       QA5   50
C        PURPOSE                                                        QA5   60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X)/SQRT(X), SUMMED OVER X   QA5   70
C                               FROM 0 TO INFINITY).                    QA5   80
C                                                                       QA5   90
C        USAGE                                                          QA5  100
C           CALL QA5 (FCT,Y)                                            QA5  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QA5  120
C                                                                       QA5  130
C        DESCRIPTION OF PARAMETERS                                      QA5  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QA5  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QA5  160
C                                                                       QA5  170
C        REMARKS                                                        QA5  180
C           NONE                                                        QA5  190
C                                                                       QA5  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QA5  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QA5  220
C           BY THE USER.                                                QA5  230
C                                                                       QA5  240
C        METHOD                                                         QA5  250
C           EVALUATION IS DONE BY MEANS OF 5-POINT GENERALIZED GAUSSIAN-QA5  260
C           LAGUERRE QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY,      QA5  270
C           WHENEVER FCT(X) IS A POLYNOMIAL UP TO DEGREE 9.             QA5  280
C           FOR REFERENCE, SEE                                          QA5  290
C           CONCUS/CASSATT/JAEHNIG/MELBY, TABLES FOR THE EVALUATION OF  QA5  300
C           INTEGRAL(X**BETA*EXP(-X)*F(X), SUMMED OVER X FROM 0 TO      QA5  310
C           INFINITY) BY GAUSS-LAGUERRE QUADRATURE, MTAC, VOL.17,       QA5  320
C           ISS.83 (1963), PP.245-256.                                  QA5  330
C                                                                       QA5  340
C     ..................................................................QA5  350
C                                                                       QA5  360
      SUBROUTINE QA5(FCT,Y)                                             QA5  370
C                                                                       QA5  380
C                                                                       QA5  390
      X=11.80719                                                        QA5  400
      Y=.1528087E-4*FCT(X)                                              QA5  410
      X=6.414730                                                        QA5  420
      Y=Y+.002687291*FCT(X)                                             QA5  430
      X=3.085937                                                        QA5  440
      Y=Y+.06774879*FCT(X)                                              QA5  450
      X=1.074562                                                        QA5  460
      Y=Y+.4802772*FCT(X)                                               QA5  470
      X=.1175813                                                        QA5  480
      Y=Y+1.221725*FCT(X)                                               QA5  490
      RETURN                                                            QA5  500
      END                                                               QA5  510
