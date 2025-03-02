C                                                                       DQHE  10
C     ..................................................................DQHE  20
C                                                                       DQHE  30
C        SUBROUTINE DQHFE                                               DQHE  40
C                                                                       DQHE  50
C        PURPOSE                                                        DQHE  60
C           TO COMPUTE THE VECTOR OF INTEGRAL VALUES FOR A GIVEN        DQHE  70
C           EQUIDISTANT TABLE OF FUNCTION AND DERIVATIVE VALUES.        DQHE  80
C                                                                       DQHE  90
C        USAGE                                                          DQHE 100
C           CALL DQHFE (H,Y,DERY,Z,NDIM)                                DQHE 110
C                                                                       DQHE 120
C        DESCRIPTION OF PARAMETERS                                      DQHE 130
C           H      - THE INCREMENT OF ARGUMENT VALUES.                  DQHE 140
C           Y      - THE INPUT VECTOR OF FUNCTION VALUES.               DQHE 150
C           DERY   - THE INPUT VECTOR OF DERIVATIVE VALUES.             DQHE 160
C           Z      - THE RESULTING VECTOR OF INTEGRAL VALUES. Z MAY BE  DQHE 170
C                    IDENTICAL WITH Y OR DERY.                          DQHE 180
C           NDIM   - THE DIMENSION OF VECTORS Y,DERY,Z.                 DQHE 190
C        NOTE: Y,DERY,Z,H,HH,HS,SUM1,SUM2 DOUBLE PRECISION              DQHE 200
C        REMARKS                                                        DQHE 210
C           NO ACTION IN CASE NDIM LESS THAN 1.                         DQHE 220
C                                                                       DQHE 230
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DQHE 240
C           NONE                                                        DQHE 250
C                                                                       DQHE 260
C        METHOD                                                         DQHE 270
C           BEGINNING WITH Z(1)=0, EVALUATION OF VECTOR Z IS DONE BY    DQHE 280
C           MEANS OF HERMITEAN FOURTH ORDER INTEGRATION FORMULA.        DQHE 290
C           FOR REFERENCE, SEE                                          DQHE 300
C           (1) F.B.HILDEBRAND, INTRODUCTION TO NUMERICAL ANALYSIS,     DQHE 310
C               MCGRAW-HILL, NEW YORK/TORONTO/LONDON, 1956, PP.314-319. DQHE 320
C           (2) R.ZURMUEHL, PRAKTISCHE MATHEMATIK FUER INGENIEURE UND   DQHE 330
C               PHYSIKER, SPRINGER, BERLIN/GOETTINGEN/HEIDELBERG, 1963, DQHE 340
C               PP.227-230.                                             DQHE 350
C                                                                       DQHE 360
C     ..................................................................DQHE 370
C                                                                       DQHE 380
      SUBROUTINE DQHFE(H,Y,DERY,Z,NDIM)                                 DQHE 390
C                                                                       DQHE 400
C                                                                       DQHE 410
      DIMENSION Y(1),DERY(1),Z(1)                                       DQHE 420
C     DOUBLE PRECISION Y,DERY,Z,H,HH,HS,SUM1,SUM2                       DQHE 430
      SUM2=0.D0                                                         DQHE 440
      IF(NDIM-1)4,3,1                                                   DQHE 450
    1 HH=.5D0*H                                                         DQHE 460
      HS=.16666666666666667D0*H                                         DQHE 470
C                                                                       DQHE 480
C     INTEGRATION LOOP                                                  DQHE 490
      DO 2 I=2,NDIM                                                     DQHE 500
      SUM1=SUM2                                                         DQHE 510
      SUM2=SUM2+HH*((Y(I)+Y(I-1))+HS*(DERY(I-1)-DERY(I)))               DQHE 520
    2 Z(I-1)=SUM1                                                       DQHE 530
    3 Z(NDIM)=SUM2                                                      DQHE 540
    4 RETURN                                                            DQHE 550
      END                                                               DQHE 560
