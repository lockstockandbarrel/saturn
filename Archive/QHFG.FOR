C                                                                       QHFG  10
C     ..................................................................QHFG  20
C                                                                       QHFG  30
C        SUBROUTINE QHFG                                                QHFG  40
C                                                                       QHFG  50
C        PURPOSE                                                        QHFG  60
C           TO COMPUTE THE VECTOR OF INTEGRAL VALUES FOR A GIVEN        QHFG  70
C           GENERAL TABLE OF ARGUMENT, FUNCTION, AND DERIVATIVE VALUES. QHFG  80
C                                                                       QHFG  90
C        USAGE                                                          QHFG 100
C           CALL QHFG (X,Y,DERY,Z,NDIM)                                 QHFG 110
C                                                                       QHFG 120
C        DESCRIPTION OF PARAMETERS                                      QHFG 130
C           X      - THE INPUT VECTOR OF ARGUMENT VALUES.               QHFG 140
C           Y      - THE INPUT VECTOR OF FUNCTION VALUES.               QHFG 150
C           DERY   - THE INPUT VECTOR OF DERIVATIVE VALUES.             QHFG 160
C           Z      - THE RESULTING VECTOR OF INTEGRAL VALUES. Z MAY BE  QHFG 170
C                    IDENTICAL WITH X,Y OR DERY.                        QHFG 180
C           NDIM   - THE DIMENSION OF VECTORS X,Y,DERY,Z.               QHFG 190
C                                                                       QHFG 200
C        REMARKS                                                        QHFG 210
C           NO ACTION IN CASE NDIM LESS THAN 1.                         QHFG 220
C                                                                       QHFG 230
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QHFG 240
C           NONE                                                        QHFG 250
C                                                                       QHFG 260
C        METHOD                                                         QHFG 270
C           BEGINNING WITH Z(1)=0, EVALUATION OF VECTOR Z IS DONE BY    QHFG 280
C           MEANS OF HERMITEAN FOURTH ORDER INTEGRATION FORMULA.        QHFG 290
C           FOR REFERENCE, SEE                                          QHFG 300
C           (1) F.B.HILDEBRAND, INTRODUCTION TO NUMERICAL ANALYSIS,     QHFG 310
C               MCGRAW-HILL, NEW YORK/TORONTO/LONDON, 1956, PP.314-319. QHFG 320
C           (2) R.ZURMUEHL, PRAKTISCHE MATHEMATIK FUER INGENIEURE UND   QHFG 330
C               PHYSIKER, SPRINGER, BERLIN/GOETTINGEN/HEIDELBERG, 1963, QHFG 340
C               PP.227-230.                                             QHFG 350
C                                                                       QHFG 360
C     ..................................................................QHFG 370
C                                                                       QHFG 380
      SUBROUTINE QHFG(X,Y,DERY,Z,NDIM)                                  QHFG 390
C                                                                       QHFG 400
C                                                                       QHFG 410
      DIMENSION X(1),Y(1),DERY(1),Z(1)                                  QHFG 420
C                                                                       QHFG 430
      SUM2=0.                                                           QHFG 440
      IF(NDIM-1)4,3,1                                                   QHFG 450
C                                                                       QHFG 460
C     INTEGRATION LOOP                                                  QHFG 470
    1 DO 2 I=2,NDIM                                                     QHFG 480
      SUM1=SUM2                                                         QHFG 490
      SUM2=.5*(X(I)-X(I-1))                                             QHFG 500
      SUM2=SUM1+SUM2*((Y(I)+Y(I-1))+.3333333*SUM2*(DERY(I-1)-DERY(I)))  QHFG 510
    2 Z(I-1)=SUM1                                                       QHFG 520
    3 Z(NDIM)=SUM2                                                      QHFG 530
    4 RETURN                                                            QHFG 540
      END                                                               QHFG 550
