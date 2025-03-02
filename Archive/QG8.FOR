C                                                                       QG8   10
C     ..................................................................QG8   20
C                                                                       QG8   30
C        SUBROUTINE QG8                                                 QG8   40
C                                                                       QG8   50
C        PURPOSE                                                        QG8   60
C           TO COMPUTE INTEGRAL(FCT(X), SUMMED OVER X FROM XL TO XU)    QG8   70
C                                                                       QG8   80
C        USAGE                                                          QG8   90
C           CALL QG8 (XL,XU,FCT,Y)                                      QG8  100
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QG8  110
C                                                                       QG8  120
C        DESCRIPTION OF PARAMETERS                                      QG8  130
C           XL     - THE LOWER BOUND OF THE INTERVAL.                   QG8  140
C           XU     - THE UPPER BOUND OF THE INTERVAL.                   QG8  150
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QG8  160
C           Y      - THE RESULTING INTEGRAL VALUE.                      QG8  170
C                                                                       QG8  180
C        REMARKS                                                        QG8  190
C           NONE                                                        QG8  200
C                                                                       QG8  210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QG8  220
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QG8  230
C           BY THE USER.                                                QG8  240
C                                                                       QG8  250
C        METHOD                                                         QG8  260
C           EVALUATION IS DONE BY MEANS OF 8-POINT GAUSS QUADRATURE     QG8  270
C           FORMULA, WHICH INTEGRATES POLYNOMIALS UP TO DEGREE 15       QG8  280
C           EXACTLY.                                                    QG8  290
C           FOR REFERENCE, SEE                                          QG8  300
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QG8  310
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.100-111 AND 337-338.   QG8  320
C                                                                       QG8  330
C     ..................................................................QG8  340
C                                                                       QG8  350
      SUBROUTINE QG8(XL,XU,FCT,Y)                                       QG8  360
C                                                                       QG8  370
C                                                                       QG8  380
      A=.5*(XU+XL)                                                      QG8  390
      B=XU-XL                                                           QG8  400
      C=.4801449*B                                                      QG8  410
      Y=.05061427*(FCT(A+C)+FCT(A-C))                                   QG8  420
      C=.3983332*B                                                      QG8  430
      Y=Y+.1111905*(FCT(A+C)+FCT(A-C))                                  QG8  440
      C=.2627662*B                                                      QG8  450
      Y=Y+.1568533*(FCT(A+C)+FCT(A-C))                                  QG8  460
      C=.09171732*B                                                     QG8  470
      Y=B*(Y+.1813419*(FCT(A+C)+FCT(A-C)))                              QG8  480
      RETURN                                                            QG8  490
      END                                                               QG8  500
