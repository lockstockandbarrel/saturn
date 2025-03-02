C                                                                       QTFE  10
C     ..................................................................QTFE  20
C                                                                       QTFE  30
C        SUBROUTINE QTFE                                                QTFE  40
C                                                                       QTFE  50
C        PURPOSE                                                        QTFE  60
C           TO COMPUTE THE VECTOR OF INTEGRAL VALUES FOR A GIVEN        QTFE  70
C           EQUIDISTANT TABLE OF FUNCTION VALUES.                       QTFE  80
C                                                                       QTFE  90
C        USAGE                                                          QTFE 100
C           CALL QTFE (H,Y,Z,NDIM)                                      QTFE 110
C                                                                       QTFE 120
C        DESCRIPTION OF PARAMETERS                                      QTFE 130
C           H      - THE INCREMENT OF ARGUMENT VALUES.                  QTFE 140
C           Y      - THE INPUT VECTOR OF FUNCTION VALUES.               QTFE 150
C           Z      - THE RESULTING VECTOR OF INTEGRAL VALUES. Z MAY BE  QTFE 160
C                    IDENTICAL WITH Y.                                  QTFE 170
C           NDIM   - THE DIMENSION OF VECTORS Y AND Z.                  QTFE 180
C                                                                       QTFE 190
C        REMARKS                                                        QTFE 200
C           NO ACTION IN CASE NDIM LESS THAN 1.                         QTFE 210
C                                                                       QTFE 220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QTFE 230
C           NONE                                                        QTFE 240
C                                                                       QTFE 250
C        METHOD                                                         QTFE 260
C           BEGINNING WITH Z(1)=0, EVALUATION OF VECTOR Z IS DONE BY    QTFE 270
C           MEANS OF TRAPEZOIDAL RULE (SECOND ORDER FORMULA).           QTFE 280
C           FOR REFERENCE, SEE                                          QTFE 290
C           F.B.HILDEBRAND, INTRODUCTION TO NUMERICAL ANALYSIS,         QTFE 300
C           MCGRAW-HILL, NEW YORK/TORONTO/LONDON, 1956, PP.75.          QTFE 310
C                                                                       QTFE 320
C     ..................................................................QTFE 330
C                                                                       QTFE 340
      SUBROUTINE QTFE(H,Y,Z,NDIM)                                       QTFE 350
C                                                                       QTFE 360
C                                                                       QTFE 370
      DIMENSION Y(1),Z(1)                                               QTFE 380
C                                                                       QTFE 390
      SUM2=0.                                                           QTFE 400
      IF(NDIM-1)4,3,1                                                   QTFE 410
    1 HH=.5*H                                                           QTFE 420
C                                                                       QTFE 430
C     INTEGRATION LOOP                                                  QTFE 440
      DO 2 I=2,NDIM                                                     QTFE 450
      SUM1=SUM2                                                         QTFE 460
      SUM2=SUM2+HH*(Y(I)+Y(I-1))                                        QTFE 470
    2 Z(I-1)=SUM1                                                       QTFE 480
    3 Z(NDIM)=SUM2                                                      QTFE 490
    4 RETURN                                                            QTFE 500
      END                                                               QTFE 510
