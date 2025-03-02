C                                                                       QA7   10
C     ..................................................................QA7   20
C                                                                       QA7   30
C        SUBROUTINE QA7                                                 QA7   40
C                                                                       QA7   50
C        PURPOSE                                                        QA7   60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X)/SQRT(X), SUMMED OVER X   QA7   70
C                               FROM 0 TO INFINITY).                    QA7   80
C                                                                       QA7   90
C        USAGE                                                          QA7  100
C           CALL QA7 (FCT,Y)                                            QA7  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QA7  120
C                                                                       QA7  130
C        DESCRIPTION OF PARAMETERS                                      QA7  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QA7  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QA7  160
C                                                                       QA7  170
C        REMARKS                                                        QA7  180
C           NONE                                                        QA7  190
C                                                                       QA7  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QA7  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QA7  220
C           BY THE USER.                                                QA7  230
C                                                                       QA7  240
C        METHOD                                                         QA7  250
C           EVALUATION IS DONE BY MEANS OF 7-POINT GENERALIZED GAUSSIAN-QA7  260
C           LAGUERRE QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY,      QA7  270
C           WHENEVER FCT(X) IS A POLYNOMIAL UP TO DEGREE 13.            QA7  280
C           FOR REFERENCE, SEE                                          QA7  290
C           CONCUS/CASSATT/JAEHNIG/MELBY, TABLES FOR THE EVALUATION OF  QA7  300
C           INTEGRAL(X**BETA*EXP(-X)*F(X), SUMMED OVER X FROM 0 TO      QA7  310
C           INFINITY) BY GAUSS-LAGUERRE QUADRATURE, MTAC, VOL.17,       QA7  320
C           ISS.83 (1963), PP.245-256.                                  QA7  330
C                                                                       QA7  340
C     ..................................................................QA7  350
C                                                                       QA7  360
      SUBROUTINE QA7(FCT,Y)                                             QA7  370
C                                                                       QA7  380
C                                                                       QA7  390
      X=18.52828                                                        QA7  400
      Y=.1725718E-7*FCT(X)                                              QA7  410
      X=11.98999                                                        QA7  420
      Y=Y+.9432969E-5*FCT(X)                                            QA7  430
      X=7.554091                                                        QA7  440
      Y=Y+.0007101852*FCT(X)                                            QA7  450
      X=4.389793                                                        QA7  460
      Y=Y+.01570011*FCT(X)                                              QA7  470
      X=2.180592                                                        QA7  480
      Y=Y+.1370111*FCT(X)                                               QA7  490
      X=.7721379                                                        QA7  500
      Y=Y+.5462112*FCT(X)                                               QA7  510
      X=.08511544                                                       QA7  520
      Y=Y+1.072812*FCT(X)                                               QA7  530
      RETURN                                                            QA7  540
      END                                                               QA7  550
