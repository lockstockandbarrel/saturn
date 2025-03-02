C                                                                       QTFG  10
C     ..................................................................QTFG  20
C                                                                       QTFG  30
C        SUBROUTINE QTFG                                                QTFG  40
C                                                                       QTFG  50
C        PURPOSE                                                        QTFG  60
C           TO COMPUTE THE VECTOR OF INTEGRAL VALUES FOR A GIVEN        QTFG  70
C           GENERAL TABLE OF ARGUMENT AND FUNCTION VALUES.              QTFG  80
C                                                                       QTFG  90
C        USAGE                                                          QTFG 100
C           CALL QTFG (X,Y,Z,NDIM)                                      QTFG 110
C                                                                       QTFG 120
C        DESCRIPTION OF PARAMETERS                                      QTFG 130
C           X      - THE INPUT VECTOR OF ARGUMENT VALUES.               QTFG 140
C           Y      - THE INPUT VECTOR OF FUNCTION VALUES.               QTFG 150
C           Z      - THE RESULTING VECTOR OF INTEGRAL VALUES. Z MAY BE  QTFG 160
C                    IDENTICAL WITH X OR Y.                             QTFG 170
C           NDIM   - THE DIMENSION OF VECTORS X,Y,Z.                    QTFG 180
C                                                                       QTFG 190
C        REMARKS                                                        QTFG 200
C           NO ACTION IN CASE NDIM LESS THAN 1.                         QTFG 210
C                                                                       QTFG 220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QTFG 230
C           NONE                                                        QTFG 240
C                                                                       QTFG 250
C        METHOD                                                         QTFG 260
C           BEGINNING WITH Z(1)=0, EVALUATION OF VECTOR Z IS DONE BY    QTFG 270
C           MEANS OF TRAPEZOIDAL RULE (SECOND ORDER FORMULA).           QTFG 280
C           FOR REFERENCE, SEE                                          QTFG 290
C           F.B.HILDEBRAND, INTRODUCTION TO NUMERICAL ANALYSIS,         QTFG 300
C           MCGRAW-HILL, NEW YORK/TORONTO/LONDON, 1956, PP.75.          QTFG 310
C                                                                       QTFG 320
C     ..................................................................QTFG 330
C                                                                       QTFG 340
      SUBROUTINE QTFG(X,Y,Z,NDIM)                                       QTFG 350
C                                                                       QTFG 360
C                                                                       QTFG 370
      DIMENSION X(1),Y(1),Z(1)                                          QTFG 380
C                                                                       QTFG 390
      SUM2=0.                                                           QTFG 400
      IF(NDIM-1)4,3,1                                                   QTFG 410
C                                                                       QTFG 420
C     INTEGRATION LOOP                                                  QTFG 430
    1 DO 2 I=2,NDIM                                                     QTFG 440
      SUM1=SUM2                                                         QTFG 450
      SUM2=SUM2+.5*(X(I)-X(I-1))*(Y(I)+Y(I-1))                          QTFG 460
    2 Z(I-1)=SUM1                                                       QTFG 470
    3 Z(NDIM)=SUM2                                                      QTFG 480
    4 RETURN                                                            QTFG 490
      END                                                               QTFG 500
