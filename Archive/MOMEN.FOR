C                                                                       MOME  10
C     ..................................................................MOME  20
C                                                                       MOME  30
C        SUBROUTINE MOMEN                                               MOME  40
C                                                                       MOME  50
C        PURPOSE                                                        MOME  60
C           TO FIND THE THE FIRST FOUR MOMENTS FOR GROUPED DATA ON      MOME  70
C           EQUAL CLASS INTERVALS.                                      MOME  80
C                                                                       MOME  90
C        USAGE                                                          MOME 100
C           CALL MOMEN (F,UBO,NOP,ANS)                                  MOME 110
C                                                                       MOME 120
C        DESCRIPTION OF PARAMETERS                                      MOME 130
C           F   - GROUPED DATA (FREQUENCIES).  GIVEN AS A VECTOR OF     MOME 140
C                 LENGTH (UBO(3)-UBO(1))/UBO(2)                         MOME 150
C           UBO - 3 CELL VECTOR, UBO(1) IS LOWER BOUND AND UBO(3) UPPER MOME 160
C                 BOUND ON DATA.  UBO(2) IS CLASS INTERVAL.  NOTE THAT  MOME 170
C                 UBO(3) MUST BE GREATER THAN UBO(1).                   MOME 180
C           NOP - OPTION PARAMETER.  IF NOP = 1, ANS(1) = MEAN.  IF     MOME 190
C                 NOP = 2, ANS(2) = SECOND MOMENT.  IF NOP = 3, ANS(3) =MOME 200
C                 THIRD MOMENT.  IF NOP = 4, ANS(4) = FOURTH MOMENT.    MOME 210
C                 IF NOP = 5, ALL FOUR MOMENTS ARE FILLED IN.           MOME 220
C           ANS - OUTPUT VECTOR OF LENGTH 4 INTO WHICH MOMENTS ARE PUT. MOME 230
C                                                                       MOME 240
C        REMARKS                                                        MOME 250
C           NOTE THAT THE FIRST MOMENT IS NOT CENTRAL BUT THE VALUE OF  MOME 260
C           THE MEAN ITSELF.  THE MEAN IS ALWAYS CALCULATED.  MOMENTS   MOME 270
C           ARE BIASED AND NOT CORRECTED FOR GROUPING.                  MOME 280
C                                                                       MOME 290
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  MOME 300
C           NONE                                                        MOME 310
C                                                                       MOME 320
C        METHOD                                                         MOME 330
C           REFER TO M. G. KENDALL, 'THE ADVANCED THEORY OF STATISTICS',MOME 340
C           V.1, HAFNER PUBLISHING COMPANY, 1958, CHAPTER 3.            MOME 350
C                                                                       MOME 360
C     ..................................................................MOME 370
C                                                                       MOME 380
      SUBROUTINE MOMEN (F,UBO,NOP,ANS)                                  MOME 390
      DIMENSION F(1),UBO(1),ANS(1)                                      MOME 400
C                                                                       MOME 410
      DO 100 I=1,4                                                      MOME 420
  100 ANS(I)=0.0                                                        MOME 430
C                                                                       MOME 440
C     CALCULATE THE NUMBER OF CLASS INTERVALS                           MOME 450
C                                                                       MOME 460
      N=(UBO(3)-UBO(1))/UBO(2)+0.5                                      MOME 470
C                                                                       MOME 480
C     CALCULATE TOTAL FREQUENCY                                         MOME 490
C                                                                       MOME 500
      T=0.0                                                             MOME 510
      DO 110 I=1,N                                                      MOME 520
  110 T=T+F(I)                                                          MOME 530
C                                                                       MOME 540
      IF(NOP-5) 130, 120, 115                                           MOME 550
  115 NOP=5                                                             MOME 560
  120 JUMP=1                                                            MOME 570
      GO TO 150                                                         MOME 580
  130 JUMP=2                                                            MOME 590
C                                                                       MOME 600
C        FIRST MOMENT                                                   MOME 610
C                                                                       MOME 620
  150 DO 160 I=1,N                                                      MOME 630
      FI=I                                                              MOME 640
  160 ANS(1)=ANS(1)+F(I)*(UBO(1)+(FI-0.5)*UBO(2))                       MOME 650
      ANS(1)=ANS(1)/T                                                   MOME 660
C                                                                       MOME 670
      GO TO (350,200,250,300,200), NOP                                  MOME 680
C                                                                       MOME 690
C        SECOND MOMENT                                                  MOME 700
C                                                                       MOME 710
  200 DO 210 I=1,N                                                      MOME 720
      FI=I                                                              MOME 730
  210 ANS(2)=ANS(2)+F(I)*(UBO(1)+(FI-0.5)*UBO(2)-ANS(1))**2             MOME 740
      ANS(2)=ANS(2)/T                                                   MOME 750
      GO TO (250,350), JUMP                                             MOME 760
C                                                                       MOME 770
C        THIRD MOMENT                                                   MOME 780
C                                                                       MOME 790
  250 DO 260 I=1,N                                                      MOME 800
      FI=I                                                              MOME 810
  260 ANS(3)=ANS(3)+F(I)*(UBO(1)+(FI-0.5)*UBO(2)-ANS(1))**3             MOME 820
      ANS(3)=ANS(3)/T                                                   MOME 830
      GO TO (300,350), JUMP                                             MOME 840
C                                                                       MOME 850
C        FOURTH MOMENT                                                  MOME 860
C                                                                       MOME 870
  300 DO 310 I=1,N                                                      MOME 880
      FI=I                                                              MOME 890
  310 ANS(4)=ANS(4)+F(I)*(UBO(1)+(FI-0.5)*UBO(2)-ANS(1))**4             MOME 900
      ANS(4)=ANS(4)/T                                                   MOME 910
  350 RETURN                                                            MOME 920
      END                                                               MOME 930
