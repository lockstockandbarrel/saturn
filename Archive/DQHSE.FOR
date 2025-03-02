C                                                                       DQHS  10
C     ..................................................................DQHS  20
C                                                                       DQHS  30
C        SUBROUTINE DQHSE                                               DQHS  40
C                                                                       DQHS  50
C        PURPOSE                                                        DQHS  60
C           TO COMPUTE THE VECTOR OF INTEGRAL VALUES FOR A GIVEN        DQHS  70
C           EQUIDISTANT TABLE OF FUNCTION, FIRST DERIVATIVE,            DQHS  80
C           AND SECOND DERIVATIVE VALUES.                               DQHS  90
C                                                                       DQHS 100
C        USAGE                                                          DQHS 110
C           CALL DQHSE (H,Y,FDY,SDY,Z,NDIM)                             DQHS 120
C                                                                       DQHS 130
C        DESCRIPTION OF PARAMETERS                                      DQHS 140
C           H      - DOUBLE PRECISION INCREMENT OF ARGUMENT VALUES.     DQHS 150
C           Y      - DOUBLE PRECISION INPUT VECTOR OF FUNCTION VALUES.  DQHS 160
C           FDY    - DOUBLE PRECISION INPUT VECTOR OF FIRST DERIVATIVE. DQHS 170
C           SDY    - DOUBLE PRECISION INPUT VECTOR OF SECOND DERIVATIVE.DQHS 180
C           Z      - RESULTING DOUBLE PRECISION VECTOR OF INTEGRAL      DQHS 190
C                    VALUES. Z MAY BE IDENTICAL WITH Y, FDY OR SDY.     DQHS 200
C           NDIM   - THE DIMENSION OF VECTORS Y,FDY,SDY,Z.              DQHS 210
C                                                                       DQHS 220
C        REMARKS                                                        DQHS 230
C           NO ACTION IN CASE NDIM LESS THAN 1.                         DQHS 240
C                                                                       DQHS 250
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DQHS 260
C           NONE                                                        DQHS 270
C                                                                       DQHS 280
C        METHOD                                                         DQHS 290
C           BEGINNING WITH Z(1)=0, EVALUATION OF VECTOR Z IS DONE BY    DQHS 300
C           MEANS OF HERMITEAN SIXTH ORDER INTEGRATION FORMULA.         DQHS 310
C           FOR REFERENCE, SEE                                          DQHS 320
C           R.ZURMUEHL, PRAKTISCHE MATHEMATIK FUER INGENIEURE UND       DQHS 330
C           PHYSIKER, SPRINGER, BERLIN/GOETTINGEN/HEIDELBERG, 1963,     DQHS 340
C           PP.227-230.                                                 DQHS 350
C                                                                       DQHS 360
C     ..................................................................DQHS 370
C                                                                       DQHS 380
      SUBROUTINE DQHSE(H,Y,FDY,SDY,Z,NDIM)                              DQHS 390
C                                                                       DQHS 400
C                                                                       DQHS 410
      DIMENSION Y(1),FDY(1),SDY(1),Z(1)                                 DQHS 420
      DOUBLE PRECISION Y,FDY,SDY,Z,H,HH,HF,HT,SUM1,SUM2                 DQHS 430
C                                                                       DQHS 440
      SUM2=0.D0                                                         DQHS 450
      IF(NDIM-1)4,3,1                                                   DQHS 460
    1 HH=.5D0*H                                                         DQHS 470
      HF=.2D0*H                                                         DQHS 480
      HT=.08333333333333333D0*H                                         DQHS 490
C                                                                       DQHS 500
C     INTEGRATION LOOP                                                  DQHS 510
      DO 2 I=2,NDIM                                                     DQHS 520
      SUM1=SUM2                                                         DQHS 530
     0SUM2=SUM2+HH*((Y(I-1)+Y(I))+HF*((FDY(I-1)-FDY(I))+                DQHS 540
     1              HT*(SDY(I-1)+SDY(I))))                              DQHS 550
    2 Z(I-1)=SUM1                                                       DQHS 560
    3 Z(NDIM)=SUM2                                                      DQHS 570
    4 RETURN                                                            DQHS 580
      END                                                               DQHS 590
