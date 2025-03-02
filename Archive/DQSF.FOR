C                                                                       DQSF  10
C     ..................................................................DQSF  20
C                                                                       DQSF  30
C        SUBROUTINE DQSF                                                DQSF  40
C                                                                       DQSF  50
C        PURPOSE                                                        DQSF  60
C           TO COMPUTE THE VECTOR OF INTEGRAL VALUES FOR A GIVEN        DQSF  70
C           EQUIDISTANT TABLE OF FUNCTION VALUES.                       DQSF  80
C                                                                       DQSF  90
C        USAGE                                                          DQSF 100
C           CALL DQSF (H,Y,Z,NDIM)                                      DQSF 110
C                                                                       DQSF 120
C        DESCRIPTION OF PARAMETERS                                      DQSF 130
C           H      - DOUBLE PRECISION INCREMENT OF ARGUMENT VALUES.     DQSF 140
C           Y      - DOUBLE PRECISION INPUT VECTOR OF FUNCTION VALUES.  DQSF 150
C           Z      - RESULTING DOUBLE PRECISION VECTOR OF INTEGRAL      DQSF 160
C                    VALUES. Z MAY BE IDENTICAL WITH Y.                 DQSF 170
C           NDIM   - THE DIMENSION OF VECTORS Y AND Z.                  DQSF 180
C                                                                       DQSF 190
C        REMARKS                                                        DQSF 200
C           NO ACTION IN CASE NDIM LESS THAN 3.                         DQSF 210
C                                                                       DQSF 220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DQSF 230
C           NONE                                                        DQSF 240
C                                                                       DQSF 250
C        METHOD                                                         DQSF 260
C           BEGINNING WITH Z(1)=0, EVALUATION OF VECTOR Z IS DONE BY    DQSF 270
C           MEANS OF SIMPSONS RULE TOGETHER WITH NEWTONS 3/8 RULE OR A  DQSF 280
C           COMBINATION OF THESE TWO RULES. TRUNCATION ERROR IS OF      DQSF 290
C           ORDER H**5 (I.E. FOURTH ORDER METHOD). ONLY IN CASE NDIM=3  DQSF 300
C           TRUNCATION ERROR OF Z(2) IS OF ORDER H**4.                  DQSF 310
C           FOR REFERENCE, SEE                                          DQSF 320
C           (1) F.B.HILDEBRAND, INTRODUCTION TO NUMERICAL ANALYSIS,     DQSF 330
C               MCGRAW-HILL, NEW YORK/TORONTO/LONDON, 1956, PP.71-76.   DQSF 340
C           (2) R.ZURMUEHL, PRAKTISCHE MATHEMATIK FUER INGENIEURE UND   DQSF 350
C               PHYSIKER, SPRINGER, BERLIN/GOETTINGEN/HEIDELBERG, 1963, DQSF 360
C               PP.214-221.                                             DQSF 370
C                                                                       DQSF 380
C     ..................................................................DQSF 390
C                                                                       DQSF 400
      SUBROUTINE DQSF(H,Y,Z,NDIM)                                       DQSF 410
C                                                                       DQSF 420
C                                                                       DQSF 430
      DIMENSION Y(1),Z(1)                                               DQSF 440
      DOUBLE PRECISION Y,Z,H,HT,SUM1,SUM2,AUX,AUX1,AUX2                 DQSF 450
C                                                                       DQSF 460
      HT=.33333333333333333D0*H                                         DQSF 470
      IF(NDIM-5)7,8,1                                                   DQSF 480
C                                                                       DQSF 490
C     NDIM IS GREATER THAN 5. PREPARATIONS OF INTEGRATION LOOP          DQSF 500
    1 SUM1=Y(2)+Y(2)                                                    DQSF 510
      SUM1=SUM1+SUM1                                                    DQSF 520
      SUM1=HT*(Y(1)+SUM1+Y(3))                                          DQSF 530
      AUX1=Y(4)+Y(4)                                                    DQSF 540
      AUX1=AUX1+AUX1                                                    DQSF 550
      AUX1=SUM1+HT*(Y(3)+AUX1+Y(5))                                     DQSF 560
      AUX2=HT*(Y(1)+3.875D0*(Y(2)+Y(5))+2.625D0*(Y(3)+Y(4))+Y(6))       DQSF 570
      SUM2=Y(5)+Y(5)                                                    DQSF 580
      SUM2=SUM2+SUM2                                                    DQSF 590
      SUM2=AUX2-HT*(Y(4)+SUM2+Y(6))                                     DQSF 600
      Z(1)=0.D0                                                         DQSF 610
      AUX=Y(3)+Y(3)                                                     DQSF 620
      AUX=AUX+AUX                                                       DQSF 630
      Z(2)=SUM2-HT*(Y(2)+AUX+Y(4))                                      DQSF 640
      Z(3)=SUM1                                                         DQSF 650
      Z(4)=SUM2                                                         DQSF 660
      IF(NDIM-6)5,5,2                                                   DQSF 670
C                                                                       DQSF 680
C     INTEGRATION LOOP                                                  DQSF 690
    2 DO 4 I=7,NDIM,2                                                   DQSF 700
      SUM1=AUX1                                                         DQSF 710
      SUM2=AUX2                                                         DQSF 720
      AUX1=Y(I-1)+Y(I-1)                                                DQSF 730
      AUX1=AUX1+AUX1                                                    DQSF 740
      AUX1=SUM1+HT*(Y(I-2)+AUX1+Y(I))                                   DQSF 750
      Z(I-2)=SUM1                                                       DQSF 760
      IF(I-NDIM)3,6,6                                                   DQSF 770
    3 AUX2=Y(I)+Y(I)                                                    DQSF 780
      AUX2=AUX2+AUX2                                                    DQSF 790
      AUX2=SUM2+HT*(Y(I-1)+AUX2+Y(I+1))                                 DQSF 800
    4 Z(I-1)=SUM2                                                       DQSF 810
    5 Z(NDIM-1)=AUX1                                                    DQSF 820
      Z(NDIM)=AUX2                                                      DQSF 830
      RETURN                                                            DQSF 840
    6 Z(NDIM-1)=SUM2                                                    DQSF 850
      Z(NDIM)=AUX1                                                      DQSF 860
      RETURN                                                            DQSF 870
C     END OF INTEGRATION LOOP                                           DQSF 880
C                                                                       DQSF 890
    7 IF(NDIM-3)12,11,8                                                 DQSF 900
C                                                                       DQSF 910
C     NDIM IS EQUAL TO 4 OR 5                                           DQSF 920
    8 SUM2=1.125D0*HT*(Y(1)+Y(2)+Y(2)+Y(2)+Y(3)+Y(3)+Y(3)+Y(4))         DQSF 930
      SUM1=Y(2)+Y(2)                                                    DQSF 940
      SUM1=SUM1+SUM1                                                    DQSF 950
      SUM1=HT*(Y(1)+SUM1+Y(3))                                          DQSF 960
      Z(1)=0.D0                                                         DQSF 970
      AUX1=Y(3)+Y(3)                                                    DQSF 980
      AUX1=AUX1+AUX1                                                    DQSF 990
      Z(2)=SUM2-HT*(Y(2)+AUX1+Y(4))                                     DQSF1000
      IF(NDIM-5)10,9,9                                                  DQSF1010
    9 AUX1=Y(4)+Y(4)                                                    DQSF1020
      AUX1=AUX1+AUX1                                                    DQSF1030
      Z(5)=SUM1+HT*(Y(3)+AUX1+Y(5))                                     DQSF1040
   10 Z(3)=SUM1                                                         DQSF1050
      Z(4)=SUM2                                                         DQSF1060
      RETURN                                                            DQSF1070
C                                                                       DQSF1080
C     NDIM IS EQUAL TO 3                                                DQSF1090
   11 SUM1=HT*(1.25D0*Y(1)+Y(2)+Y(2)-.25D0*Y(3))                        DQSF1100
      SUM2=Y(2)+Y(2)                                                    DQSF1110
      SUM2=SUM2+SUM2                                                    DQSF1120
      Z(3)=HT*(Y(1)+SUM2+Y(3))                                          DQSF1130
      Z(1)=0.D0                                                         DQSF1140
      Z(2)=SUM1                                                         DQSF1150
   12 RETURN                                                            DQSF1160
      END                                                               DQSF1170
