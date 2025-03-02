C                                                                       QG2   10
C     ..................................................................QG2   20
C                                                                       QG2   30
C        SUBROUTINE QG2                                                 QG2   40
C                                                                       QG2   50
C        PURPOSE                                                        QG2   60
C           TO COMPUTE INTEGRAL(FCT(X), SUMMED OVER X FROM XL TO XU)    QG2   70
C                                                                       QG2   80
C        USAGE                                                          QG2   90
C           CALL QG2 (XL,XU,FCT,Y)                                      QG2  100
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QG2  110
C                                                                       QG2  120
C        DESCRIPTION OF PARAMETERS                                      QG2  130
C           XL     - THE LOWER BOUND OF THE INTERVAL.                   QG2  140
C           XU     - THE UPPER BOUND OF THE INTERVAL.                   QG2  150
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QG2  160
C           Y      - THE RESULTING INTEGRAL VALUE.                      QG2  170
C                                                                       QG2  180
C        REMARKS                                                        QG2  190
C           NONE                                                        QG2  200
C                                                                       QG2  210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QG2  220
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QG2  230
C           BY THE USER.                                                QG2  240
C                                                                       QG2  250
C        METHOD                                                         QG2  260
C           EVALUATION IS DONE BY MEANS OF 2-POINT GAUSS QUADRATURE     QG2  270
C           FORMULA, WHICH INTEGRATES POLYNOMIALS UP TO DEGREE 3        QG2  280
C           EXACTLY.                                                    QG2  290
C           FOR REFERENCE, SEE                                          QG2  300
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QG2  310
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.100-111 AND 337-338.   QG2  320
C                                                                       QG2  330
C     ..................................................................QG2  340
C                                                                       QG2  350
      SUBROUTINE QG2(XL,XU,FCT,Y)                                       QG2  360
C                                                                       QG2  370
C                                                                       QG2  380
      A=.5*(XU+XL)                                                      QG2  390
      B=XU-XL                                                           QG2  400
      Y=.2886751*B                                                      QG2  410
      Y=.5*B*(FCT(A+Y)+FCT(A-Y))                                        QG2  420
      RETURN                                                            QG2  430
      END                                                               QG2  440
