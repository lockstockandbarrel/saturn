C                                                                       DHSG  10
C     ..................................................................DHSG  20
C                                                                       DHSG  30
C        SUBROUTINE DQHSG                                               DHSG  40
C                                                                       DHSG  50
C        PURPOSE                                                        DHSG  60
C           TO COMPUTE THE VECTOR OF INTEGRAL VALUES FOR A GIVEN        DHSG  70
C           GENERAL TABLE OF ARGUMENT, FUNCTION, FIRST DERIVATIVE,      DHSG  80
C           AND SECOND DERIVATIVE VALUES.                               DHSG  90
C                                                                       DHSG 100
C        USAGE                                                          DHSG 110
C           CALL DQHSG (X,Y,FDY,SDY,Z,NDIM)                             DHSG 120
C                                                                       DHSG 130
C        DESCRIPTION OF PARAMETERS                                      DHSG 140
C           X      - DOUBLE PRECISION INPUT VECTOR OF ARGUMENT VALUES.  DHSG 150
C           Y      - DOUBLE PRECISION INPUT VECTOR OF FUNCTION VALUES.  DHSG 160
C           FDY    - DOUBLE PRECISION INPUT VECTOR OF FIRST DERIVATIVE. DHSG 170
C           SDY    - DOUBLE PRECISION INPUT VECTOR OF SECOND DERIVATIVE.DHSG 180
C           Z      - RESULTING DOUBLE PRECISION VECTOR OF INTEGRAL      DHSG 190
C                    VALUES. Z MAY BE IDENTICAL WITH X, Y, FDY OR SDY.  DHSG 200
C           NDIM   - THE DIMENSION OF VECTORS X,Y,FDY,SDY,Z.            DHSG 210
C                                                                       DHSG 220
C        REMARKS                                                        DHSG 230
C           NO ACTION IN CASE NDIM LESS THAN 1.                         DHSG 240
C                                                                       DHSG 250
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DHSG 260
C           NONE                                                        DHSG 270
C                                                                       DHSG 280
C        METHOD                                                         DHSG 290
C           BEGINNING WITH Z(1)=0, EVALUATION OF VECTOR Z IS DONE BY    DHSG 300
C           MEANS OF HERMITEAN SIXTH ORDER INTEGRATION FORMULA.         DHSG 310
C           FOR REFERENCE, SEE                                          DHSG 320
C           R.ZURMUEHL, PRAKTISCHE MATHEMATIK FUER INGENIEURE UND       DHSG 330
C           PHYSIKER, SPRINGER, BERLIN/GOETTINGEN/HEIDELBERG, 1963,     DHSG 340
C           PP.227-230.                                                 DHSG 350
C                                                                       DHSG 360
C     ..................................................................DHSG 370
C                                                                       DHSG 380
      SUBROUTINE DQHSG(X,Y,FDY,SDY,Z,NDIM)                              DHSG 390
C                                                                       DHSG 400
C                                                                       DHSG 410
      DIMENSION X(1),Y(1),FDY(1),SDY(1),Z(1)                            DHSG 420
      DOUBLE PRECISION X,Y,FDY,SDY,Z,SUM1,SUM2                          DHSG 430
C                                                                       DHSG 440
      SUM2=0.D0                                                         DHSG 450
      IF(NDIM-1)4,3,1                                                   DHSG 460
C                                                                       DHSG 470
C     INTEGRATION LOOP                                                  DHSG 480
    1 DO 2 I=2,NDIM                                                     DHSG 490
      SUM1=SUM2                                                         DHSG 500
      SUM2=.5D0*(X(I)-X(I-1))                                           DHSG 510
     0SUM2=SUM1+SUM2*((Y(I-1)+Y(I))+.4D0*SUM2*((FDY(I-1)-FDY(I))+       DHSG 520
     1     .16666666666666667D0*SUM2*(SDY(I-1)+SDY(I))))                DHSG 530
    2 Z(I-1)=SUM1                                                       DHSG 540
    3 Z(NDIM)=SUM2                                                      DHSG 550
    4 RETURN                                                            DHSG 560
      END                                                               DHSG 570
