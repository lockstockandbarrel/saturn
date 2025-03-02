C                                                                       DQHG  10
C     ..................................................................DQHG  20
C                                                                       DQHG  30
C        SUBROUTINE DQHFG                                               DQHG  40
C                                                                       DQHG  50
C        PURPOSE                                                        DQHG  60
C           TO COMPUTE THE VECTOR OF INTEGRAL VALUES FOR A GIVEN        DQHG  70
C           GENERAL TABLE OF ARGUMENT, FUNCTION, AND DERIVATIVE VALUES. DQHG  80
C                                                                       DQHG  90
C        USAGE                                                          DQHG 100
C           CALL DQHFG (X,Y,DERY,Z,NDIM)                                DQHG 110
C                                                                       DQHG 120
C        DESCRIPTION OF PARAMETERS                                      DQHG 130
C           X      - THE INPUT VECTOR OF ARGUMENT VALUES.               DQHG 140
C           Y      - THE INPUT VECTOR OF FUNCTION VALUES.               DQHG 150
C           DERY   - THE INPUT VECTOR OF DERIVATIVE VALUES.             DQHG 160
C           Z      - THE RESULTING VECTOR OF INTEGRAL VALUES. Z MAY BE  DQHG 170
C                    IDENTICAL WITH X,Y OR DERY.                        DQHG 180
C           NDIM   - THE DIMENSION OF VECTORS X,Y,DERY,Z.               DQHG 190
C        NOTE:X,Y,DERY,Z,SUM1,SUM2 ARE DOUBLE PRECISION                 DQHG 200
C        REMARKS                                                        DQHG 210
C           NO ACTION IN CASE NDIM LESS THAN 1.                         DQHG 220
C                                                                       DQHG 230
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DQHG 240
C           NONE                                                        DQHG 250
C                                                                       DQHG 260
C        METHOD                                                         DQHG 270
C           BEGINNING WITH Z(1)=0, EVALUATION OF VECTOR Z IS DONE BY    DQHG 280
C           MEANS OF HERMITEAN FOURTH ORDER INTEGRATION FORMULA.        DQHG 290
C           FOR REFERENCE, SEE                                          DQHG 300
C           (1) F.B.HILDEBRAND, INTRODUCTION TO NUMERICAL ANALYSIS,     DQHG 310
C               MCGRAW-HILL, NEW YORK/TORONTO/LONDON, 1956, PP.314-319. DQHG 320
C           (2) R.ZURMUEHL, PRAKTISCHE MATHEMATIK FUER INGENIEURE UND   DQHG 330
C               PHYSIKER, SPRINGER, BERLIN/GOETTINGEN/HEIDELBERG, 1963, DQHG 340
C               PP.227-230.                                             DQHG 350
C                                                                       DQHG 360
C     ..................................................................DQHG 370
C                                                                       DQHG 380
      SUBROUTINE DQHFG(X,Y,DERY,Z,NDIM)                                 DQHG 390
C                                                                       DQHG 400
C                                                                       DQHG 410
      DIMENSION X(1),Y(1),DERY(1),Z(1)                                  DQHG 420
      DOUBLE PRECISION X,Y,DERY,Z,SUM1,SUM2                             DQHG 430
      SUM2=0.D0                                                         DQHG 440
      IF(NDIM-1)4,3,1                                                   DQHG 450
C                                                                       DQHG 460
C     INTEGRATION LOOP                                                  DQHG 470
    1 DO 2 I=2,NDIM                                                     DQHG 480
      SUM1=SUM2                                                         DQHG 490
      SUM2=.5D0*(X(I)-X(I-1))                                           DQHG 500
      SUM2=SUM1+SUM2*((Y(I)+Y(I-1))+.33333333333333333D0*SUM2*          DQHG 510
    1 (DERY(I-1)-DERY(I)))                                              DQHG 520
    2 Z(I-1)=SUM1                                                       DQHG 530
    3 Z(NDIM)=SUM2                                                      DQHG 540
    4 RETURN                                                            DQHG 550
      END                                                               DQHG 560
