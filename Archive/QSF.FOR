C                                                                       QSF   10
C     ..................................................................QSF   20
C                                                                       QSF   30
C        SUBROUTINE QSF                                                 QSF   40
C                                                                       QSF   50
C        PURPOSE                                                        QSF   60
C           TO COMPUTE THE VECTOR OF INTEGRAL VALUES FOR A GIVEN        QSF   70
C           EQUIDISTANT TABLE OF FUNCTION VALUES.                       QSF   80
C                                                                       QSF   90
C        USAGE                                                          QSF  100
C           CALL QSF (H,Y,Z,NDIM)                                       QSF  110
C                                                                       QSF  120
C        DESCRIPTION OF PARAMETERS                                      QSF  130
C           H      - THE INCREMENT OF ARGUMENT VALUES.                  QSF  140
C           Y      - THE INPUT VECTOR OF FUNCTION VALUES.               QSF  150
C           Z      - THE RESULTING VECTOR OF INTEGRAL VALUES. Z MAY BE  QSF  160
C                    IDENTICAL WITH Y.                                  QSF  170
C           NDIM   - THE DIMENSION OF VECTORS Y AND Z.                  QSF  180
C                                                                       QSF  190
C        REMARKS                                                        QSF  200
C           NO ACTION IN CASE NDIM LESS THAN 3.                         QSF  210
C                                                                       QSF  220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QSF  230
C           NONE                                                        QSF  240
C                                                                       QSF  250
C        METHOD                                                         QSF  260
C           BEGINNING WITH Z(1)=0, EVALUATION OF VECTOR Z IS DONE BY    QSF  270
C           MEANS OF SIMPSONS RULE TOGETHER WITH NEWTONS 3/8 RULE OR A  QSF  280
C           COMBINATION OF THESE TWO RULES. TRUNCATION ERROR IS OF      QSF  290
C           ORDER H**5 (I.E. FOURTH ORDER METHOD). ONLY IN CASE NDIM=3  QSF  300
C           TRUNCATION ERROR OF Z(2) IS OF ORDER H**4.                  QSF  310
C           FOR REFERENCE, SEE                                          QSF  320
C           (1) F.B.HILDEBRAND, INTRODUCTION TO NUMERICAL ANALYSIS,     QSF  330
C               MCGRAW-HILL, NEW YORK/TORONTO/LONDON, 1956, PP.71-76.   QSF  340
C           (2) R.ZURMUEHL, PRAKTISCHE MATHEMATIK FUER INGENIEURE UND   QSF  350
C               PHYSIKER, SPRINGER, BERLIN/GOETTINGEN/HEIDELBERG, 1963, QSF  360
C               PP.214-221.                                             QSF  370
C                                                                       QSF  380
C     ..................................................................QSF  390
C                                                                       QSF  400
      SUBROUTINE QSF(H,Y,Z,NDIM)                                        QSF  410
C                                                                       QSF  420
C                                                                       QSF  430
      DIMENSION Y(1),Z(1)                                               QSF  440
C                                                                       QSF  450
      HT=.3333333*H                                                     QSF  460
      IF(NDIM-5)7,8,1                                                   QSF  470
C                                                                       QSF  480
C     NDIM IS GREATER THAN 5. PREPARATIONS OF INTEGRATION LOOP          QSF  490
    1 SUM1=Y(2)+Y(2)                                                    QSF  500
      SUM1=SUM1+SUM1                                                    QSF  510
      SUM1=HT*(Y(1)+SUM1+Y(3))                                          QSF  520
      AUX1=Y(4)+Y(4)                                                    QSF  530
      AUX1=AUX1+AUX1                                                    QSF  540
      AUX1=SUM1+HT*(Y(3)+AUX1+Y(5))                                     QSF  550
      AUX2=HT*(Y(1)+3.875*(Y(2)+Y(5))+2.625*(Y(3)+Y(4))+Y(6))           QSF  560
      SUM2=Y(5)+Y(5)                                                    QSF  570
      SUM2=SUM2+SUM2                                                    QSF  580
      SUM2=AUX2-HT*(Y(4)+SUM2+Y(6))                                     QSF  590
      Z(1)=0.                                                           QSF  600
      AUX=Y(3)+Y(3)                                                     QSF  610
      AUX=AUX+AUX                                                       QSF  620
      Z(2)=SUM2-HT*(Y(2)+AUX+Y(4))                                      QSF  630
      Z(3)=SUM1                                                         QSF  640
      Z(4)=SUM2                                                         QSF  650
      IF(NDIM-6)5,5,2                                                   QSF  660
C                                                                       QSF  670
C     INTEGRATION LOOP                                                  QSF  680
    2 DO 4 I=7,NDIM,2                                                   QSF  690
      SUM1=AUX1                                                         QSF  700
      SUM2=AUX2                                                         QSF  710
      AUX1=Y(I-1)+Y(I-1)                                                QSF  720
      AUX1=AUX1+AUX1                                                    QSF  730
      AUX1=SUM1+HT*(Y(I-2)+AUX1+Y(I))                                   QSF  740
      Z(I-2)=SUM1                                                       QSF  750
      IF(I-NDIM)3,6,6                                                   QSF  760
    3 AUX2=Y(I)+Y(I)                                                    QSF  770
      AUX2=AUX2+AUX2                                                    QSF  780
      AUX2=SUM2+HT*(Y(I-1)+AUX2+Y(I+1))                                 QSF  790
    4 Z(I-1)=SUM2                                                       QSF  800
    5 Z(NDIM-1)=AUX1                                                    QSF  810
      Z(NDIM)=AUX2                                                      QSF  820
      RETURN                                                            QSF  830
    6 Z(NDIM-1)=SUM2                                                    QSF  840
      Z(NDIM)=AUX1                                                      QSF  850
      RETURN                                                            QSF  860
C     END OF INTEGRATION LOOP                                           QSF  870
C                                                                       QSF  880
    7 IF(NDIM-3)12,11,8                                                 QSF  890
C                                                                       QSF  900
C     NDIM IS EQUAL TO 4 OR 5                                           QSF  910
    8 SUM2=1.125*HT*(Y(1)+Y(2)+Y(2)+Y(2)+Y(3)+Y(3)+Y(3)+Y(4))           QSF  920
      SUM1=Y(2)+Y(2)                                                    QSF  930
      SUM1=SUM1+SUM1                                                    QSF  940
      SUM1=HT*(Y(1)+SUM1+Y(3))                                          QSF  950
      Z(1)=0.                                                           QSF  960
      AUX1=Y(3)+Y(3)                                                    QSF  970
      AUX1=AUX1+AUX1                                                    QSF  980
      Z(2)=SUM2-HT*(Y(2)+AUX1+Y(4))                                     QSF  990
      IF(NDIM-5)10,9,9                                                  QSF 1000
    9 AUX1=Y(4)+Y(4)                                                    QSF 1010
      AUX1=AUX1+AUX1                                                    QSF 1020
      Z(5)=SUM1+HT*(Y(3)+AUX1+Y(5))                                     QSF 1030
   10 Z(3)=SUM1                                                         QSF 1040
      Z(4)=SUM2                                                         QSF 1050
      RETURN                                                            QSF 1060
C                                                                       QSF 1070
C     NDIM IS EQUAL TO 3                                                QSF 1080
   11 SUM1=HT*(1.25*Y(1)+Y(2)+Y(2)-.25*Y(3))                            QSF 1090
      SUM2=Y(2)+Y(2)                                                    QSF 1100
      SUM2=SUM2+SUM2                                                    QSF 1110
      Z(3)=HT*(Y(1)+SUM2+Y(3))                                          QSF 1120
      Z(1)=0.                                                           QSF 1130
      Z(2)=SUM1                                                         QSF 1140
   12 RETURN                                                            QSF 1150
      END                                                               QSF 1160
