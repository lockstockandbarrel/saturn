C                                                                       QHSG  10
C     ..................................................................QHSG  20
C                                                                       QHSG  30
C        SUBROUTINE QHSG                                                QHSG  40
C                                                                       QHSG  50
C        PURPOSE                                                        QHSG  60
C           TO COMPUTE THE VECTOR OF INTEGRAL VALUES FOR A GIVEN        QHSG  70
C           GENERAL TABLE OF ARGUMENT, FUNCTION, FIRST DERIVATIVE,      QHSG  80
C           AND SECOND DERIVATIVE VALUES.                               QHSG  90
C                                                                       QHSG 100
C        USAGE                                                          QHSG 110
C           CALL QHSG (X,Y,FDY,SDY,Z,NDIM)                              QHSG 120
C                                                                       QHSG 130
C        DESCRIPTION OF PARAMETERS                                      QHSG 140
C           X      - THE INPUT VECTOR OF ARGUMENT VALUES.               QHSG 150
C           Y      - THE INPUT VECTOR OF FUNCTION VALUES.               QHSG 160
C           FDY    - THE INPUT VECTOR OF FIRST DERIVATIVE.              QHSG 170
C           SDY    - THE INPUT VECTOR OF SECOND DERIVATIVE.             QHSG 180
C           Z      - THE RESULTING VECTOR OF INTEGRAL VALUES. Z MAY BE  QHSG 190
C                    IDENTICAL WITH X,Y,FDY OR SDY.                     QHSG 200
C           NDIM   - THE DIMENSION OF VECTORS X,Y,FDY,SDY,Z.            QHSG 210
C                                                                       QHSG 220
C        REMARKS                                                        QHSG 230
C           NO ACTION IN CASE NDIM LESS THAN 1.                         QHSG 240
C                                                                       QHSG 250
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QHSG 260
C           NONE                                                        QHSG 270
C                                                                       QHSG 280
C        METHOD                                                         QHSG 290
C           BEGINNING WITH Z(1)=0, EVALUATION OF VECTOR Z IS DONE BY    QHSG 300
C           MEANS OF HERMITEAN SIXTH ORDER INTEGRATION FORMULA.         QHSG 310
C           FOR REFERENCE, SEE                                          QHSG 320
C           R.ZURMUEHL, PRAKTISCHE MATHEMATIK FUER INGENIEURE UND       QHSG 330
C           PHYSIKER, SPRINGER, BERLIN/GOETTINGEN/HEIDELBERG, 1963,     QHSG 340
C           PP.227-230.                                                 QHSG 350
C                                                                       QHSG 360
C     ..................................................................QHSG 370
C                                                                       QHSG 380
      SUBROUTINE QHSG(X,Y,FDY,SDY,Z,NDIM)                               QHSG 390
C                                                                       QHSG 400
C                                                                       QHSG 410
      DIMENSION X(1),Y(1),FDY(1),SDY(1),Z(1)                            QHSG 420
C                                                                       QHSG 430
      SUM2=0.                                                           QHSG 440
      IF(NDIM-1)4,3,1                                                   QHSG 450
C                                                                       QHSG 460
C     INTEGRATION LOOP                                                  QHSG 470
    1 DO 2 I=2,NDIM                                                     QHSG 480
      SUM1=SUM2                                                         QHSG 490
      SUM2=.5*(X(I)-X(I-1))                                             QHSG 500
     0SUM2=SUM1+SUM2*((Y(I-1)+Y(I))+.4*SUM2*((FDY(I-1)-FDY(I))+         QHSG 510
     1     .1666667*SUM2*(SDY(I-1)+SDY(I))))                            QHSG 520
    2 Z(I-1)=SUM1                                                       QHSG 530
    3 Z(NDIM)=SUM2                                                      QHSG 540
    4 RETURN                                                            QHSG 550
      END                                                               QHSG 560
