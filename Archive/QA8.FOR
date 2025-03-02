C                                                                       QA8   10
C     ..................................................................QA8   20
C                                                                       QA8   30
C        SUBROUTINE QA8                                                 QA8   40
C                                                                       QA8   50
C        PURPOSE                                                        QA8   60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X)/SQRT(X), SUMMED OVER X   QA8   70
C                               FROM 0 TO INFINITY).                    QA8   80
C                                                                       QA8   90
C        USAGE                                                          QA8  100
C           CALL QA8 (FCT,Y)                                            QA8  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QA8  120
C                                                                       QA8  130
C        DESCRIPTION OF PARAMETERS                                      QA8  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QA8  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QA8  160
C                                                                       QA8  170
C        REMARKS                                                        QA8  180
C           NONE                                                        QA8  190
C                                                                       QA8  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QA8  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QA8  220
C           BY THE USER.                                                QA8  230
C                                                                       QA8  240
C        METHOD                                                         QA8  250
C           EVALUATION IS DONE BY MEANS OF 8-POINT GENERALIZED GAUSSIAN-QA8  260
C           LAGUERRE QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY,      QA8  270
C           WHENEVER FCT(X) IS A POLYNOMIAL UP TO DEGREE 15.            QA8  280
C           FOR REFERENCE, SEE                                          QA8  290
C           CONCUS/CASSATT/JAEHNIG/MELBY, TABLES FOR THE EVALUATION OF  QA8  300
C           INTEGRAL(X**BETA*EXP(-X)*F(X), SUMMED OVER X FROM 0 TO      QA8  310
C           INFINITY) BY GAUSS-LAGUERRE QUADRATURE, MTAC, VOL.17,       QA8  320
C           ISS.83 (1963), PP.245-256.                                  QA8  330
C                                                                       QA8  340
C     ..................................................................QA8  350
C                                                                       QA8  360
      SUBROUTINE QA8(FCT,Y)                                             QA8  370
C                                                                       QA8  380
C                                                                       QA8  390
      X=21.98427                                                        QA8  400
      Y=.5309615E-9*FCT(X)                                              QA8  410
      X=14.97262                                                        QA8  420
      Y=Y+.4641962E-6*FCT(X)                                            QA8  430
      X=10.09332                                                        QA8  440
      Y=Y+.5423720E-4*FCT(X)                                            QA8  450
      X=6.483145                                                        QA8  460
      Y=Y+.001864568*FCT(X)                                             QA8  470
      X=3.809476                                                        QA8  480
      Y=Y+.02576062*FCT(X)                                              QA8  490
      X=1.905114                                                        QA8  500
      Y=Y+.1676201*FCT(X)                                               QA8  510
      X=.6772491                                                        QA8  520
      Y=Y+.5612949*FCT(X)                                               QA8  530
      X=.07479188                                                       QA8  540
      Y=Y+1.015859*FCT(X)                                               QA8  550
      RETURN                                                            QA8  560
      END                                                               QA8  570
