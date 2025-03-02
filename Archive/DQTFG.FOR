C                                                                       DTFG  10
C     ..................................................................DTFG  20
C                                                                       DTFG  30
C        SUBROUTINE DQTFG                                               DTFG  40
C                                                                       DTFG  50
C        PURPOSE                                                        DTFG  60
C           TO COMPUTE THE VECTOR OF INTEGRAL VALUES FOR A GIVEN        DTFG  70
C           GENERAL TABLE OF ARGUMENT AND FUNCTION VALUES.              DTFG  80
C                                                                       DTFG  90
C        USAGE                                                          DTFG 100
C           CALL DQTFG (X,Y,Z,NDIM)                                     DTFG 110
C                                                                       DTFG 120
C        DESCRIPTION OF PARAMETERS                                      DTFG 130
C           X      - DOUBLE PRECISION INPUT VECTOR OF ARGUMENT VALUES.  DTFG 140
C           Y      - DOUBLE PRECISION INPUT VECTOR OF FUNCTION VALUES.  DTFG 150
C           Z      - RESULTING DOUBLE PRECISION VECTOR OF INTEGRAL      DTFG 160
C                    VALUES. Z MAY BE IDENTICAL WITH X OR Y.            DTFG 170
C           NDIM   - THE DIMENSION OF VECTORS X,Y,Z.                    DTFG 180
C                                                                       DTFG 190
C        REMARKS                                                        DTFG 200
C           NO ACTION IN CASE NDIM LESS THAN 1.                         DTFG 210
C                                                                       DTFG 220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DTFG 230
C           NONE                                                        DTFG 240
C                                                                       DTFG 250
C        METHOD                                                         DTFG 260
C           BEGINNING WITH Z(1)=0, EVALUATION OF VECTOR Z IS DONE BY    DTFG 270
C           MEANS OF TRAPEZOIDAL RULE (SECOND ORDER FORMULA).           DTFG 280
C           FOR REFERENCE, SEE                                          DTFG 290
C           F.B.HILDEBRAND, INTRODUCTION TO NUMERICAL ANALYSIS,         DTFG 300
C           MCGRAW-HILL, NEW YORK/TORONTO/LONDON, 1956, PP.75.          DTFG 310
C                                                                       DTFG 320
C     ..................................................................DTFG 330
C                                                                       DTFG 340
      SUBROUTINE DQTFG(X,Y,Z,NDIM)                                      DTFG 350
C                                                                       DTFG 360
C                                                                       DTFG 370
      DIMENSION X(1),Y(1),Z(1)                                          DTFG 380
      DOUBLE PRECISION X,Y,Z,SUM1,SUM2                                  DTFG 390
C                                                                       DTFG 400
      SUM2=0.D0                                                         DTFG 410
      IF(NDIM-1)4,3,1                                                   DTFG 420
C                                                                       DTFG 430
C     INTEGRATION LOOP                                                  DTFG 440
    1 DO 2 I=2,NDIM                                                     DTFG 450
      SUM1=SUM2                                                         DTFG 460
      SUM2=SUM2+.5D0*(X(I)-X(I-1))*(Y(I)+Y(I-1))                        DTFG 470
    2 Z(I-1)=SUM1                                                       DTFG 480
    3 Z(NDIM)=SUM2                                                      DTFG 490
    4 RETURN                                                            DTFG 500
      END                                                               DTFG 510
