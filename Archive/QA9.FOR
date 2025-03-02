C                                                                       QA9   10
C     ..................................................................QA9   20
C                                                                       QA9   30
C        SUBROUTINE QA9                                                 QA9   40
C                                                                       QA9   50
C        PURPOSE                                                        QA9   60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X)/SQRT(X), SUMMED OVER X   QA9   70
C                               FROM 0 TO INFINITY).                    QA9   80
C                                                                       QA9   90
C        USAGE                                                          QA9  100
C           CALL QA9 (FCT,Y)                                            QA9  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QA9  120
C                                                                       QA9  130
C        DESCRIPTION OF PARAMETERS                                      QA9  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QA9  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QA9  160
C                                                                       QA9  170
C        REMARKS                                                        QA9  180
C           NONE                                                        QA9  190
C                                                                       QA9  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QA9  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QA9  220
C           BY THE USER.                                                QA9  230
C                                                                       QA9  240
C        METHOD                                                         QA9  250
C           EVALUATION IS DONE BY MEANS OF 9-POINT GENERALIZED GAUSSIAN-QA9  260
C           LAGUERRE QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY,      QA9  270
C           WHENEVER FCT(X) IS A POLYNOMIAL UP TO DEGREE 17.            QA9  280
C           FOR REFERENCE, SEE                                          QA9  290
C           CONCUS/CASSATT/JAEHNIG/MELBY, TABLES FOR THE EVALUATION OF  QA9  300
C           INTEGRAL(X**BETA*EXP(-X)*F(X), SUMMED OVER X FROM 0 TO      QA9  310
C           INFINITY) BY GAUSS-LAGUERRE QUADRATURE, MTAC, VOL.17,       QA9  320
C           ISS.83 (1963), PP.245-256.                                  QA9  330
C                                                                       QA9  340
C     ..................................................................QA9  350
C                                                                       QA9  360
      SUBROUTINE QA9(FCT,Y)                                             QA9  370
C                                                                       QA9  380
C                                                                       QA9  390
      X=25.48598                                                        QA9  400
      Y=.1565640E-10*FCT(X)                                             QA9  410
      X=18.04651                                                        QA9  420
      Y=Y+.2093441E-7*FCT(X)                                            QA9  430
      X=12.77183                                                        QA9  440
      Y=Y+.3621309E-5*FCT(X)                                            QA9  450
      X=8.769757                                                        QA9  460
      Y=Y+.0001836225*FCT(X)                                            QA9  470
      X=5.694423                                                        QA9  480
      Y=Y+.003777045*FCT(X)                                             QA9  490
      X=3.369176                                                        QA9  500
      Y=Y+.03728008*FCT(X)                                              QA9  510
      X=1.692395                                                        QA9  520
      Y=Y+.1946035*FCT(X)                                               QA9  530
      X=.6032364                                                        QA9  540
      Y=Y+.5696146*FCT(X)                                               QA9  550
      X=.06670223                                                       QA9  560
      Y=Y+.9669914*FCT(X)                                               QA9  570
      RETURN                                                            QA9  580
      END                                                               QA9  590
