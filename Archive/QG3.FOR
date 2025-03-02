C                                                                       QG3   10
C     ..................................................................QG3   20
C                                                                       QG3   30
C        SUBROUTINE QG3                                                 QG3   40
C                                                                       QG3   50
C        PURPOSE                                                        QG3   60
C           TO COMPUTE INTEGRAL(FCT(X), SUMMED OVER X FROM XL TO XU)    QG3   70
C                                                                       QG3   80
C        USAGE                                                          QG3   90
C           CALL QG3 (XL,XU,FCT,Y)                                      QG3  100
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QG3  110
C                                                                       QG3  120
C        DESCRIPTION OF PARAMETERS                                      QG3  130
C           XL     - THE LOWER BOUND OF THE INTERVAL.                   QG3  140
C           XU     - THE UPPER BOUND OF THE INTERVAL.                   QG3  150
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QG3  160
C           Y      - THE RESULTING INTEGRAL VALUE.                      QG3  170
C                                                                       QG3  180
C        REMARKS                                                        QG3  190
C           NONE                                                        QG3  200
C                                                                       QG3  210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QG3  220
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QG3  230
C           BY THE USER.                                                QG3  240
C                                                                       QG3  250
C        METHOD                                                         QG3  260
C           EVALUATION IS DONE BY MEANS OF 3-POINT GAUSS QUADRATURE     QG3  270
C           FORMULA, WHICH INTEGRATES POLYNOMIALS UP TO DEGREE 5        QG3  280
C           EXACTLY.                                                    QG3  290
C           FOR REFERENCE, SEE                                          QG3  300
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QG3  310
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.100-111 AND 337-338.   QG3  320
C                                                                       QG3  330
C     ..................................................................QG3  340
C                                                                       QG3  350
      SUBROUTINE QG3(XL,XU,FCT,Y)                                       QG3  360
C                                                                       QG3  370
C                                                                       QG3  380
      A=.5*(XU+XL)                                                      QG3  390
      B=XU-XL                                                           QG3  400
      Y=.3872983*B                                                      QG3  410
      Y=.2777778*(FCT(A+Y)+FCT(A-Y))                                    QG3  420
      Y=B*(Y+.4444444*FCT(A))                                           QG3  430
      RETURN                                                            QG3  440
      END                                                               QG3  450
