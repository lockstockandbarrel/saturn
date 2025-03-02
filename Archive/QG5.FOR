C                                                                       QG5   10
C     ..................................................................QG5   20
C                                                                       QG5   30
C        SUBROUTINE QG5                                                 QG5   40
C                                                                       QG5   50
C        PURPOSE                                                        QG5   60
C           TO COMPUTE INTEGRAL(FCT(X), SUMMED OVER X FROM XL TO XU)    QG5   70
C                                                                       QG5   80
C        USAGE                                                          QG5   90
C           CALL QG5 (XL,XU,FCT,Y)                                      QG5  100
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QG5  110
C                                                                       QG5  120
C        DESCRIPTION OF PARAMETERS                                      QG5  130
C           XL     - THE LOWER BOUND OF THE INTERVAL.                   QG5  140
C           XU     - THE UPPER BOUND OF THE INTERVAL.                   QG5  150
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QG5  160
C           Y      - THE RESULTING INTEGRAL VALUE.                      QG5  170
C                                                                       QG5  180
C        REMARKS                                                        QG5  190
C           NONE                                                        QG5  200
C                                                                       QG5  210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QG5  220
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QG5  230
C           BY THE USER.                                                QG5  240
C                                                                       QG5  250
C        METHOD                                                         QG5  260
C           EVALUATION IS DONE BY MEANS OF 5-POINT GAUSS QUADRATURE     QG5  270
C           FORMULA, WHICH INTEGRATES POLYNOMIALS UP TO DEGREE 9        QG5  280
C           EXACTLY.                                                    QG5  290
C           FOR REFERENCE, SEE                                          QG5  300
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QG5  310
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.100-111 AND 337-338.   QG5  320
C                                                                       QG5  330
C     ..................................................................QG5  340
C                                                                       QG5  350
      SUBROUTINE QG5(XL,XU,FCT,Y)                                       QG5  360
C                                                                       QG5  370
C                                                                       QG5  380
      A=.5*(XU+XL)                                                      QG5  390
      B=XU-XL                                                           QG5  400
      C=.4530899*B                                                      QG5  410
      Y=.1184634*(FCT(A+C)+FCT(A-C))                                    QG5  420
      C=.2692347*B                                                      QG5  430
      Y=Y+.2393143*(FCT(A+C)+FCT(A-C))                                  QG5  440
      Y=B*(Y+.2844444*FCT(A))                                           QG5  450
      RETURN                                                            QG5  460
      END                                                               QG5  470
