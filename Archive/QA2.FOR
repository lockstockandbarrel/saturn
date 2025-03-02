C                                                                       QA2   10
C     ..................................................................QA2   20
C                                                                       QA2   30
C        SUBROUTINE QA2                                                 QA2   40
C                                                                       QA2   50
C        PURPOSE                                                        QA2   60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X)/SQRT(X), SUMMED OVER X   QA2   70
C                               FROM 0 TO INFINITY).                    QA2   80
C                                                                       QA2   90
C        USAGE                                                          QA2  100
C           CALL QA2 (FCT,Y)                                            QA2  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QA2  120
C                                                                       QA2  130
C        DESCRIPTION OF PARAMETERS                                      QA2  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QA2  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QA2  160
C                                                                       QA2  170
C        REMARKS                                                        QA2  180
C           NONE                                                        QA2  190
C                                                                       QA2  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QA2  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QA2  220
C           BY THE USER.                                                QA2  230
C                                                                       QA2  240
C        METHOD                                                         QA2  250
C           EVALUATION IS DONE BY MEANS OF 2-POINT GENERALIZED GAUSSIAN-QA2  260
C           LAGUERRE QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY,      QA2  270
C           WHENEVER FCT(X) IS A POLYNOMIAL UP TO DEGREE 3.             QA2  280
C           FOR REFERENCE, SEE                                          QA2  290
C           CONCUS/CASSATT/JAEHNIG/MELBY, TABLES FOR THE EVALUATION OF  QA2  300
C           INTEGRAL(X**BETA*EXP(-X)*F(X), SUMMED OVER X FROM 0 TO      QA2  310
C           INFINITY) BY GAUSS-LAGUERRE QUADRATURE, MTAC, VOL.17,       QA2  320
C           ISS.83 (1963), PP.245-256.                                  QA2  330
C                                                                       QA2  340
C     ..................................................................QA2  350
C                                                                       QA2  360
      SUBROUTINE QA2(FCT,Y)                                             QA2  370
C                                                                       QA2  380
C                                                                       QA2  390
      X=2.724745                                                        QA2  400
      Y=.1626257*FCT(X)                                                 QA2  410
      X=.2752551                                                        QA2  420
      Y=Y+1.609828*FCT(X)                                               QA2  430
      RETURN                                                            QA2  440
      END                                                               QA2  450
