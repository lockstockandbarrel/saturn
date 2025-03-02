C                                                                       QA4   10
C     ..................................................................QA4   20
C                                                                       QA4   30
C        SUBROUTINE QA4                                                 QA4   40
C                                                                       QA4   50
C        PURPOSE                                                        QA4   60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X)/SQRT(X), SUMMED OVER X   QA4   70
C                               FROM 0 TO INFINITY).                    QA4   80
C                                                                       QA4   90
C        USAGE                                                          QA4  100
C           CALL QA4 (FCT,Y)                                            QA4  110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QA4  120
C                                                                       QA4  130
C        DESCRIPTION OF PARAMETERS                                      QA4  140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QA4  150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QA4  160
C                                                                       QA4  170
C        REMARKS                                                        QA4  180
C           NONE                                                        QA4  190
C                                                                       QA4  200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QA4  210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QA4  220
C           BY THE USER.                                                QA4  230
C                                                                       QA4  240
C        METHOD                                                         QA4  250
C           EVALUATION IS DONE BY MEANS OF 4-POINT GENERALIZED GAUSSIAN-QA4  260
C           LAGUERRE QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY,      QA4  270
C           WHENEVER FCT(X) IS A POLYNOMIAL UP TO DEGREE 7.             QA4  280
C           FOR REFERENCE, SEE                                          QA4  290
C           CONCUS/CASSATT/JAEHNIG/MELBY, TABLES FOR THE EVALUATION OF  QA4  300
C           INTEGRAL(X**BETA*EXP(-X)*F(X), SUMMED OVER X FROM 0 TO      QA4  310
C           INFINITY) BY GAUSS-LAGUERRE QUADRATURE, MTAC, VOL.17,       QA4  320
C           ISS.83 (1963), PP.245-256.                                  QA4  330
C                                                                       QA4  340
C     ..................................................................QA4  350
C                                                                       QA4  360
      SUBROUTINE QA4(FCT,Y)                                             QA4  370
C                                                                       QA4  380
C                                                                       QA4  390
      X=8.588636                                                        QA4  400
      Y=.0003992081*FCT(X)                                              QA4  410
      X=3.926964                                                        QA4  420
      Y=Y+.03415597*FCT(X)                                              QA4  430
      X=1.339097                                                        QA4  440
      Y=Y+.4156047*FCT(X)                                               QA4  450
      X=.1453035                                                        QA4  460
      Y=Y+1.322294*FCT(X)                                               QA4  470
      RETURN                                                            QA4  480
      END                                                               QA4  490
