C                                                                       QA10  10
C     ..................................................................QA10  20
C                                                                       QA10  30
C        SUBROUTINE QA10                                                QA10  40
C                                                                       QA10  50
C        PURPOSE                                                        QA10  60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X)/SQRT(X), SUMMED OVER X   QA10  70
C                               FROM 0 TO INFINITY).                    QA10  80
C                                                                       QA10  90
C        USAGE                                                          QA10 100
C           CALL QA10 (FCT,Y)                                           QA10 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QA10 120
C                                                                       QA10 130
C        DESCRIPTION OF PARAMETERS                                      QA10 140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QA10 150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QA10 160
C                                                                       QA10 170
C        REMARKS                                                        QA10 180
C           NONE                                                        QA10 190
C                                                                       QA10 200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QA10 210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QA10 220
C           BY THE USER.                                                QA10 230
C                                                                       QA10 240
C        METHOD                                                         QA10 250
C           EVALUATION IS DONE BY MEANS OF 10-POINT GENERALIZED GAUSS-  QA10 260
C           LAGUERRE QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY,      QA10 270
C           WHENEVER FCT(X) IS A POLYNOMIAL UP TO DEGREE 19.            QA10 280
C           FOR REFERENCE, SEE                                          QA10 290
C           CONCUS/CASSATT/JAEHNIG/MELBY, TABLES FOR THE EVALUATION OF  QA10 300
C           INTEGRAL(X**BETA*EXP(-X)*F(X), SUMMED OVER X FROM 0 TO      QA10 310
C           INFINITY) BY GAUSS-LAGUERRE QUADRATURE, MTAC, VOL.17,       QA10 320
C           ISS.83 (1963), PP.245-256.                                  QA10 330
C                                                                       QA10 340
C     ..................................................................QA10 350
C                                                                       QA10 360
      SUBROUTINE QA10(FCT,Y)                                            QA10 370
C                                                                       QA10 380
C                                                                       QA10 390
      X=29.02495                                                        QA10 400
      Y=.4458787E-12*FCT(X)                                             QA10 410
      X=21.19389                                                        QA10 420
      Y=Y+.8798682E-9*FCT(X)                                            QA10 430
      X=15.56116                                                        QA10 440
      Y=Y+.2172139E-6*FCT(X)                                            QA10 450
      X=11.20813                                                        QA10 460
      Y=Y+.1560511E-4*FCT(X)                                            QA10 470
      X=7.777439                                                        QA10 480
      Y=Y+.0004566773*FCT(X)                                            QA10 490
      X=5.084908                                                        QA10 500
      Y=Y+.006487547*FCT(X)                                             QA10 510
      X=3.022513                                                        QA10 520
      Y=Y+.04962104*FCT(X)                                              QA10 530
      X=1.522944                                                        QA10 540
      Y=Y+.2180344*FCT(X)                                               QA10 550
      X=.5438675                                                        QA10 560
      Y=Y+.5733510*FCT(X)                                               QA10 570
      X=.06019206                                                       QA10 580
      Y=Y+.9244873*FCT(X)                                               QA10 590
      RETURN                                                            QA10 600
      END                                                               QA10 610
