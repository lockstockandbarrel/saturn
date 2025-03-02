C                                                                       QG6   10
C     ..................................................................QG6   20
C                                                                       QG6   30
C        SUBROUTINE QG6                                                 QG6   40
C                                                                       QG6   50
C        PURPOSE                                                        QG6   60
C           TO COMPUTE INTEGRAL(FCT(X), SUMMED OVER X FROM XL TO XU)    QG6   70
C                                                                       QG6   80
C        USAGE                                                          QG6   90
C           CALL QG6 (XL,XU,FCT,Y)                                      QG6  100
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QG6  110
C                                                                       QG6  120
C        DESCRIPTION OF PARAMETERS                                      QG6  130
C           XL     - THE LOWER BOUND OF THE INTERVAL.                   QG6  140
C           XU     - THE UPPER BOUND OF THE INTERVAL.                   QG6  150
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QG6  160
C           Y      - THE RESULTING INTEGRAL VALUE.                      QG6  170
C                                                                       QG6  180
C        REMARKS                                                        QG6  190
C           NONE                                                        QG6  200
C                                                                       QG6  210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QG6  220
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QG6  230
C           BY THE USER.                                                QG6  240
C                                                                       QG6  250
C        METHOD                                                         QG6  260
C           EVALUATION IS DONE BY MEANS OF 6-POINT GAUSS QUADRATURE     QG6  270
C           FORMULA, WHICH INTEGRATES POLYNOMIALS UP TO DEGREE 11       QG6  280
C           EXACTLY.                                                    QG6  290
C           FOR REFERENCE, SEE                                          QG6  300
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QG6  310
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.100-111 AND 337-338.   QG6  320
C                                                                       QG6  330
C     ..................................................................QG6  340
C                                                                       QG6  350
      SUBROUTINE QG6(XL,XU,FCT,Y)                                       QG6  360
C                                                                       QG6  370
C                                                                       QG6  380
      A=.5*(XU+XL)                                                      QG6  390
      B=XU-XL                                                           QG6  400
      C=.4662348*B                                                      QG6  410
      Y=.08566225*(FCT(A+C)+FCT(A-C))                                   QG6  420
      C=.3306047*B                                                      QG6  430
      Y=Y+.1803808*(FCT(A+C)+FCT(A-C))                                  QG6  440
      C=.1193096*B                                                      QG6  450
      Y=B*(Y+.2339570*(FCT(A+C)+FCT(A-C)))                              QG6  460
      RETURN                                                            QG6  470
      END                                                               QG6  480
