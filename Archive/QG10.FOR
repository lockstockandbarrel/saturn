C                                                                       QG10  10
C     ..................................................................QG10  20
C                                                                       QG10  30
C        SUBROUTINE QG10                                                QG10  40
C                                                                       QG10  50
C        PURPOSE                                                        QG10  60
C           TO COMPUTE INTEGRAL(FCT(X), SUMMED OVER X FROM XL TO XU)    QG10  70
C                                                                       QG10  80
C        USAGE                                                          QG10  90
C           CALL QG10(XL,XU,FCT,Y)                                      QG10 100
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QG10 110
C                                                                       QG10 120
C        DESCRIPTION OF PARAMETERS                                      QG10 130
C           XL     - THE LOWER BOUND OF THE INTERVAL.                   QG10 140
C           XU     - THE UPPER BOUND OF THE INTERVAL.                   QG10 150
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QG10 160
C           Y      - THE RESULTING INTEGRAL VALUE.                      QG10 170
C                                                                       QG10 180
C        REMARKS                                                        QG10 190
C           NONE                                                        QG10 200
C                                                                       QG10 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QG10 220
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QG10 230
C           BY THE USER.                                                QG10 240
C                                                                       QG10 250
C        METHOD                                                         QG10 260
C           EVALUATION IS DONE BY MEANS OF 10-POINT GAUSS QUADRATURE    QG10 270
C           FORMULA, WHICH INTEGRATES POLYNOMIALS UP TO DEGREE 19       QG10 280
C           EXACTLY.                                                    QG10 290
C           FOR REFERENCE, SEE                                          QG10 300
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QG10 310
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.100-111 AND 337-338.   QG10 320
C                                                                       QG10 330
C     ..................................................................QG10 340
C                                                                       QG10 350
      SUBROUTINE QG10(XL,XU,FCT,Y)                                      QG10 360
C                                                                       QG10 370
C                                                                       QG10 380
      A=.5*(XU+XL)                                                      QG10 390
      B=XU-XL                                                           QG10 400
      C=.4869533*B                                                      QG10 410
      Y=.03333567*(FCT(A+C)+FCT(A-C))                                   QG10 420
      C=.4325317*B                                                      QG10 430
      Y=Y+.07472567*(FCT(A+C)+FCT(A-C))                                 QG10 440
      C=.3397048*B                                                      QG10 450
      Y=Y+.1095432*(FCT(A+C)+FCT(A-C))                                  QG10 460
      C=.2166977*B                                                      QG10 470
      Y=Y+.1346334*(FCT(A+C)+FCT(A-C))                                  QG10 480
      C=.07443717*B                                                     QG10 490
      Y=B*(Y+.1477621*(FCT(A+C)+FCT(A-C)))                              QG10 500
      RETURN                                                            QG10 510
      END                                                               QG10 520
