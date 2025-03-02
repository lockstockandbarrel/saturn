C                                                                       QG4   10
C     ..................................................................QG4   20
C                                                                       QG4   30
C        SUBROUTINE QG4                                                 QG4   40
C                                                                       QG4   50
C        PURPOSE                                                        QG4   60
C           TO COMPUTE INTEGRAL(FCT(X), SUMMED OVER X FROM XL TO XU)    QG4   70
C                                                                       QG4   80
C        USAGE                                                          QG4   90
C           CALL QG4 (XL,XU,FCT,Y)                                      QG4  100
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QG4  110
C                                                                       QG4  120
C        DESCRIPTION OF PARAMETERS                                      QG4  130
C           XL     - THE LOWER BOUND OF THE INTERVAL.                   QG4  140
C           XU     - THE UPPER BOUND OF THE INTERVAL.                   QG4  150
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QG4  160
C           Y      - THE RESULTING INTEGRAL VALUE.                      QG4  170
C                                                                       QG4  180
C        REMARKS                                                        QG4  190
C           NONE                                                        QG4  200
C                                                                       QG4  210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QG4  220
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QG4  230
C           BY THE USER.                                                QG4  240
C                                                                       QG4  250
C        METHOD                                                         QG4  260
C           EVALUATION IS DONE BY MEANS OF 4-POINT GAUSS QUADRATURE     QG4  270
C           FORMULA, WHICH INTEGRATES POLYNOMIALS UP TO DEGREE 7        QG4  280
C           EXACTLY.                                                    QG4  290
C           FOR REFERENCE, SEE                                          QG4  300
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QG4  310
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.100-111 AND 337-338.   QG4  320
C                                                                       QG4  330
C     ..................................................................QG4  340
C                                                                       QG4  350
      SUBROUTINE QG4(XL,XU,FCT,Y)                                       QG4  360
C                                                                       QG4  370
C                                                                       QG4  380
      A=.5*(XU+XL)                                                      QG4  390
      B=XU-XL                                                           QG4  400
      C=.4305682*B                                                      QG4  410
      Y=.1739274*(FCT(A+C)+FCT(A-C))                                    QG4  420
      C=.1699905*B                                                      QG4  430
      Y=B*(Y+.3260726*(FCT(A+C)+FCT(A-C)))                              QG4  440
      RETURN                                                            QG4  450
      END                                                               QG4  460
