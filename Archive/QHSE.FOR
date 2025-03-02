C                                                                       QHSE  10
C     ..................................................................QHSE  20
C                                                                       QHSE  30
C        SUBROUTINE QHSE                                                QHSE  40
C                                                                       QHSE  50
C        PURPOSE                                                        QHSE  60
C           TO COMPUTE THE VECTOR OF INTEGRAL VALUES FOR A GIVEN        QHSE  70
C           EQUIDISTANT TABLE OF FUNCTION, FIRST DERIVATIVE,            QHSE  80
C           AND SECOND DERIVATIVE VALUES.                               QHSE  90
C                                                                       QHSE 100
C        USAGE                                                          QHSE 110
C           CALL QHSE (H,Y,FDY,SDY,Z,NDIM)                              QHSE 120
C                                                                       QHSE 130
C        DESCRIPTION OF PARAMETERS                                      QHSE 140
C           H      - THE INCREMENT OF ARGUMENT VALUES.                  QHSE 150
C           Y      - THE INPUT VECTOR OF FUNCTION VALUES.               QHSE 160
C           FDY    - THE INPUT VECTOR OF FIRST DERIVATIVE.              QHSE 170
C           SDY    - THE INPUT VECTOR OF SECOND DERIVATIVE.             QHSE 180
C           Z      - THE RESULTING VECTOR OF INTEGRAL VALUES. Z MAY BE  QHSE 190
C                    IDENTICAL WITH Y,FDY OR SDY.                       QHSE 200
C           NDIM   - THE DIMENSION OF VECTORS Y,FDY,SDY,Z.              QHSE 210
C                                                                       QHSE 220
C        REMARKS                                                        QHSE 230
C           NO ACTION IN CASE NDIM LESS THAN 1.                         QHSE 240
C                                                                       QHSE 250
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QHSE 260
C           NONE                                                        QHSE 270
C                                                                       QHSE 280
C        METHOD                                                         QHSE 290
C           BEGINNING WITH Z(1)=0, EVALUATION OF VECTOR Z IS DONE BY    QHSE 300
C           MEANS OF HERMITEAN SIXTH ORDER INTEGRATION FORMULA.         QHSE 310
C           FOR REFERENCE, SEE                                          QHSE 320
C           R.ZURMUEHL, PRAKTISCHE MATHEMATIK FUER INGENIEURE UND       QHSE 330
C           PHYSIKER, SPRINGER, BERLIN/GOETTINGEN/HEIDELBERG, 1963,     QHSE 340
C           PP.227-230.                                                 QHSE 350
C                                                                       QHSE 360
C     ..................................................................QHSE 370
C                                                                       QHSE 380
      SUBROUTINE QHSE(H,Y,FDY,SDY,Z,NDIM)                               QHSE 390
C                                                                       QHSE 400
C                                                                       QHSE 410
      DIMENSION Y(1),FDY(1),SDY(1),Z(1)                                 QHSE 420
C                                                                       QHSE 430
      SUM2=0.                                                           QHSE 440
      IF(NDIM-1)4,3,1                                                   QHSE 450
    1 HH=.5*H                                                           QHSE 460
      HF=.2*H                                                           QHSE 470
      HT=.08333333*H                                                    QHSE 480
C                                                                       QHSE 490
C     INTEGRATION LOOP                                                  QHSE 500
      DO 2 I=2,NDIM                                                     QHSE 510
      SUM1=SUM2                                                         QHSE 520
     0SUM2=SUM2+HH*((Y(I-1)+Y(I))+HF*((FDY(I-1)-FDY(I))+                QHSE 530
     1              HT*(SDY(I-1)+SDY(I))))                              QHSE 540
    2 Z(I-1)=SUM1                                                       QHSE 550
    3 Z(NDIM)=SUM2                                                      QHSE 560
    4 RETURN                                                            QHSE 570
      END                                                               QHSE 580
