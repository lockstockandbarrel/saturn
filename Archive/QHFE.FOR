C                                                                       QHFE  10
C     ..................................................................QHFE  20
C                                                                       QHFE  30
C        SUBROUTINE QHFE                                                QHFE  40
C                                                                       QHFE  50
C        PURPOSE                                                        QHFE  60
C           TO COMPUTE THE VECTOR OF INTEGRAL VALUES FOR A GIVEN        QHFE  70
C           EQUIDISTANT TABLE OF FUNCTION AND DERIVATIVE VALUES.        QHFE  80
C                                                                       QHFE  90
C        USAGE                                                          QHFE 100
C           CALL QHFE (H,Y,DERY,Z,NDIM)                                 QHFE 110
C                                                                       QHFE 120
C        DESCRIPTION OF PARAMETERS                                      QHFE 130
C           H      - THE INCREMENT OF ARGUMENT VALUES.                  QHFE 140
C           Y      - THE INPUT VECTOR OF FUNCTION VALUES.               QHFE 150
C           DERY   - THE INPUT VECTOR OF DERIVATIVE VALUES.             QHFE 160
C           Z      - THE RESULTING VECTOR OF INTEGRAL VALUES. Z MAY BE  QHFE 170
C                    IDENTICAL WITH Y OR DERY.                          QHFE 180
C           NDIM   - THE DIMENSION OF VECTORS Y,DERY,Z.                 QHFE 190
C                                                                       QHFE 200
C        REMARKS                                                        QHFE 210
C           NO ACTION IN CASE NDIM LESS THAN 1.                         QHFE 220
C                                                                       QHFE 230
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QHFE 240
C           NONE                                                        QHFE 250
C                                                                       QHFE 260
C        METHOD                                                         QHFE 270
C           BEGINNING WITH Z(1)=0, EVALUATION OF VECTOR Z IS DONE BY    QHFE 280
C           MEANS OF HERMITEAN FOURTH ORDER INTEGRATION FORMULA.        QHFE 290
C           FOR REFERENCE, SEE                                          QHFE 300
C           (1) F.B.HILDEBRAND, INTRODUCTION TO NUMERICAL ANALYSIS,     QHFE 310
C               MCGRAW-HILL, NEW YORK/TORONTO/LONDON, 1956, PP.314-319. QHFE 320
C           (2) R.ZURMUEHL, PRAKTISCHE MATHEMATIK FUER INGENIEURE UND   QHFE 330
C               PHYSIKER, SPRINGER, BERLIN/GOETTINGEN/HEIDELBERG, 1963, QHFE 340
C               PP.227-230.                                             QHFE 350
C                                                                       QHFE 360
C     ..................................................................QHFE 370
C                                                                       QHFE 380
      SUBROUTINE QHFE(H,Y,DERY,Z,NDIM)                                  QHFE 390
C                                                                       QHFE 400
C                                                                       QHFE 410
      DIMENSION Y(1),DERY(1),Z(1)                                       QHFE 420
C                                                                       QHFE 430
      SUM2=0.                                                           QHFE 440
      IF(NDIM-1)4,3,1                                                   QHFE 450
    1 HH=.5*H                                                           QHFE 460
      HS=.1666667*H                                                     QHFE 470
C                                                                       QHFE 480
C     INTEGRATION LOOP                                                  QHFE 490
      DO 2 I=2,NDIM                                                     QHFE 500
      SUM1=SUM2                                                         QHFE 510
      SUM2=SUM2+HH*((Y(I)+Y(I-1))+HS*(DERY(I-1)-DERY(I)))               QHFE 520
    2 Z(I-1)=SUM1                                                       QHFE 530
    3 Z(NDIM)=SUM2                                                      QHFE 540
    4 RETURN                                                            QHFE 550
      END                                                               QHFE 560
