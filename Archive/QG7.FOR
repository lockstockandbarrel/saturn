C                                                                       QG7   10
C     ..................................................................QG7   20
C                                                                       QG7   30
C        SUBROUTINE QG7                                                 QG7   40
C                                                                       QG7   50
C        PURPOSE                                                        QG7   60
C           TO COMPUTE INTEGRAL(FCT(X), SUMMED OVER X FROM XL TO XU)    QG7   70
C                                                                       QG7   80
C        USAGE                                                          QG7   90
C           CALL QG7 (XL,XU,FCT,Y)                                      QG7  100
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QG7  110
C                                                                       QG7  120
C        DESCRIPTION OF PARAMETERS                                      QG7  130
C           XL     - THE LOWER BOUND OF THE INTERVAL.                   QG7  140
C           XU     - THE UPPER BOUND OF THE INTERVAL.                   QG7  150
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QG7  160
C           Y      - THE RESULTING INTEGRAL VALUE.                      QG7  170
C                                                                       QG7  180
C        REMARKS                                                        QG7  190
C           NONE                                                        QG7  200
C                                                                       QG7  210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QG7  220
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QG7  230
C           BY THE USER.                                                QG7  240
C                                                                       QG7  250
C        METHOD                                                         QG7  260
C           EVALUATION IS DONE BY MEANS OF 7-POINT GAUSS QUADRATURE     QG7  270
C           FORMULA, WHICH INTEGRATES POLYNOMIALS UP TO DEGREE 13       QG7  280
C           EXACTLY.                                                    QG7  290
C           FOR REFERENCE, SEE                                          QG7  300
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QG7  310
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.100-111 AND 337-338.   QG7  320
C                                                                       QG7  330
C     ..................................................................QG7  340
C                                                                       QG7  350
      SUBROUTINE QG7(XL,XU,FCT,Y)                                       QG7  360
C                                                                       QG7  370
C                                                                       QG7  380
      A=.5*(XU+XL)                                                      QG7  390
      B=XU-XL                                                           QG7  400
      C=.4745540*B                                                      QG7  410
      Y=.06474248*(FCT(A+C)+FCT(A-C))                                   QG7  420
      C=.3707656*B                                                      QG7  430
      Y=Y+.1398527*(FCT(A+C)+FCT(A-C))                                  QG7  440
      C=.2029226*B                                                      QG7  450
      Y=Y+.1909150*(FCT(A+C)+FCT(A-C))                                  QG7  460
      Y=B*(Y+.2089796*FCT(A))                                           QG7  470
      RETURN                                                            QG7  480
      END                                                               QG7  490
