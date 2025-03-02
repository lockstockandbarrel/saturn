C                                                                       QG9   10
C     ..................................................................QG9   20
C                                                                       QG9   30
C        SUBROUTINE QG9                                                 QG9   40
C                                                                       QG9   50
C        PURPOSE                                                        QG9   60
C           TO COMPUTE INTEGRAL(FCT(X), SUMMED OVER X FROM XL TO XU)    QG9   70
C                                                                       QG9   80
C        USAGE                                                          QG9   90
C           CALL QG9 (XL,XU,FCT,Y)                                      QG9  100
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QG9  110
C                                                                       QG9  120
C        DESCRIPTION OF PARAMETERS                                      QG9  130
C           XL     - THE LOWER BOUND OF THE INTERVAL.                   QG9  140
C           XU     - THE UPPER BOUND OF THE INTERVAL.                   QG9  150
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QG9  160
C           Y      - THE RESULTING INTEGRAL VALUE.                      QG9  170
C                                                                       QG9  180
C        REMARKS                                                        QG9  190
C           NONE                                                        QG9  200
C                                                                       QG9  210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QG9  220
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QG9  230
C           BY THE USER.                                                QG9  240
C                                                                       QG9  250
C        METHOD                                                         QG9  260
C           EVALUATION IS DONE BY MEANS OF 9-POINT GAUSS QUADRATURE     QG9  270
C           FORMULA, WHICH INTEGRATES POLYNOMIALS UP TO DEGREE 17       QG9  280
C           EXACTLY.                                                    QG9  290
C           FOR REFERENCE, SEE                                          QG9  300
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QG9  310
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.100-111 AND 337-338.   QG9  320
C                                                                       QG9  330
C     ..................................................................QG9  340
C                                                                       QG9  350
      SUBROUTINE QG9(XL,XU,FCT,Y)                                       QG9  360
C                                                                       QG9  370
C                                                                       QG9  380
      A=.5*(XU+XL)                                                      QG9  390
      B=XU-XL                                                           QG9  400
      C=.4840801*B                                                      QG9  410
      Y=.04063719*(FCT(A+C)+FCT(A-C))                                   QG9  420
      C=.4180156*B                                                      QG9  430
      Y=Y+.09032408*(FCT(A+C)+FCT(A-C))                                 QG9  440
      C=.3066857*B                                                      QG9  450
      Y=Y+.1303053*(FCT(A+C)+FCT(A-C))                                  QG9  460
      C=.1621267*B                                                      QG9  470
      Y=Y+.1561735*(FCT(A+C)+FCT(A-C))                                  QG9  480
      Y=B*(Y+.1651197*FCT(A))                                           QG9  490
      RETURN                                                            QG9  500
      END                                                               QG9  510
