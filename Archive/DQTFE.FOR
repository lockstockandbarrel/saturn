C                                                                       DTFE  10
C     ..................................................................DTFE  20
C                                                                       DTFE  30
C        SUBROUTINE DQTFE                                               DTFE  40
C                                                                       DTFE  50
C        PURPOSE                                                        DTFE  60
C           TO COMPUTE THE VECTOR OF INTEGRAL VALUES FOR A GIVEN        DTFE  70
C           EQUIDISTANT TABLE OF FUNCTION VALUES.                       DTFE  80
C                                                                       DTFE  90
C        USAGE                                                          DTFE 100
C           CALL DQTFE (H,Y,Z,NDIM)                                     DTFE 110
C                                                                       DTFE 120
C        DESCRIPTION OF PARAMETERS                                      DTFE 130
C           H      - DOUBLE PRECISION INCREMENT OF ARGUMENT VALUES.     DTFE 140
C           Y      - DOUBLE PRECISION INPUT VECTOR OF FUNCTION VALUES.  DTFE 150
C           Z      - RESULTING DOUBLE PRECISION VECTOR OF INTEGRAL      DTFE 160
C                    VALUES. Z MAY BE IDENTICAL WITH Y.                 DTFE 170
C           NDIM   - THE DIMENSION OF VECTORS Y AND Z.                  DTFE 180
C                                                                       DTFE 190
C        REMARKS                                                        DTFE 200
C           NO ACTION IN CASE NDIM LESS THAN 1.                         DTFE 210
C                                                                       DTFE 220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DTFE 230
C           NONE                                                        DTFE 240
C                                                                       DTFE 250
C        METHOD                                                         DTFE 260
C           BEGINNING WITH Z(1)=0, EVALUATION OF VECTOR Z IS DONE BY    DTFE 270
C           MEANS OF TRAPEZOIDAL RULE (SECOND ORDER FORMULA).           DTFE 280
C           FOR REFERENCE, SEE                                          DTFE 290
C           F.B.HILDEBRAND, INTRODUCTION TO NUMERICAL ANALYSIS,         DTFE 300
C           MCGRAW-HILL, NEW YORK/TORONTO/LONDON, 1956, PP.75.          DTFE 310
C                                                                       DTFE 320
C     ..................................................................DTFE 330
C                                                                       DTFE 340
      SUBROUTINE DQTFE(H,Y,Z,NDIM)                                      DTFE 350
C                                                                       DTFE 360
C                                                                       DTFE 370
      DIMENSION Y(1),Z(1)                                               DTFE 380
      DOUBLE PRECISION Y,Z,H,HH,SUM1,SUM2                               DTFE 390
C                                                                       DTFE 400
      SUM2=0.D0                                                         DTFE 410
      IF(NDIM-1)4,3,1                                                   DTFE 420
    1 HH=.5D0*H                                                         DTFE 430
C                                                                       DTFE 440
C     INTEGRATION LOOP                                                  DTFE 450
      DO 2 I=2,NDIM                                                     DTFE 460
      SUM1=SUM2                                                         DTFE 470
      SUM2=SUM2+HH*(Y(I)+Y(I-1))                                        DTFE 480
    2 Z(I-1)=SUM1                                                       DTFE 490
    3 Z(NDIM)=SUM2                                                      DTFE 500
    4 RETURN                                                            DTFE 510
      END                                                               DTFE 520
