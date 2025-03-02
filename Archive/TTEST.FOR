C                                                                       TTES  10
C     ..................................................................TTES  20
C                                                                       TTES  30
C        SUBROUTINE TTEST                                               TTES  40
C                                                                       TTES  50
C        PURPOSE                                                        TTES  60
C           TO FIND CERTAIN T-STATISTICS ON THE MEANS OF POPULATIONS.   TTES  70
C                                                                       TTES  80
C        USAGE                                                          TTES  90
C           CALL TTEST (A,NA,B,NB,NOP,NDF,ANS)                          TTES 100
C                                                                       TTES 110
C        DESCRIPTION OF PARAMETERS                                      TTES 120
C           A   - INPUT VECTOR OF LENGTH NA CONTAINING DATA.            TTES 130
C           NA  - NUMBER OF OBSERVATIONS IN A.                          TTES 140
C           B   - INPUT VECTOR OF LENGTH NB CONTAINING DATA.            TTES 150
C           NB  - NUMBER OF OBSERVATIONS IN B.                          TTES 160
C           NOP - OPTIONS FOR VARIOUS HYPOTHESES..                      TTES 170
C                 NOP=1--- THAT POPULATION MEAN OF B = GIVEN VALUE A.   TTES 180
C                          (SET NA=1)                                   TTES 190
C                 NOP=2--- THAT POPULATION MEAN OF B = POPULATION MEAN  TTES 200
C                          OF A, GIVEN THAT THE VARIANCE OF B = THE     TTES 210
C                          VARIANCE OF A.                               TTES 220
C                 NOP=3--- THAT POPULATION MEAN OF B = POPULATION MEAN  TTES 230
C                          OF A, GIVEN THAT THE VARIANCE OF B IS NOT    TTES 240
C                          EQUAL TO THE VARIANCE OF A.                  TTES 250
C                 NOP=4--- THAT POPULATION MEAN OF B = POPULATION MEAN  TTES 260
C                          OF A, GIVEN NO INFORMATION ABOUT VARIANCES OFTTES 270
C                          A AND B.  (SET NA=NB)                        TTES 280
C           NDF - OUTPUT VARIABLE CONTAINING DEGREES OF FREEDOM ASSOCI- TTES 290
C                 ATED WITH T-STATISTIC CALCULATED.                     TTES 300
C           ANS - T-STATISTIC FOR GIVEN HYPOTHESIS.                     TTES 310
C                                                                       TTES 320
C        REMARKS                                                        TTES 330
C           NA AND NB MUST BE GREATER THAN 1, EXCEPT THAT NA=1 IN       TTES 340
C           OPTION 1. NA AND NB MUST BE THE SAME IN OPTION 4.           TTES 350
C           IF NOP IS OTHER THAN 1, 2, 3 OR 4, DEGREES OF FREEDOM AND   TTES 360
C           T-STATISTIC WILL NOT BE CALCULATED.  NDF AND ANS WILL BE    TTES 370
C           SET TO ZERO.                                                TTES 380
C                                                                       TTES 390
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  TTES 400
C           NONE                                                        TTES 410
C                                                                       TTES 420
C        METHOD                                                         TTES 430
C           REFER TO OSTLE, BERNARD, 'STATISTICS IN RESEARCH', IOWA     TTES 440
C           STATE COLLEGE PRESS, 1954, CHAPTER 5.                       TTES 450
C                                                                       TTES 460
C     ..................................................................TTES 470
C                                                                       TTES 480
      SUBROUTINE TTEST (A,NA,B,NB,NOP,NDF,ANS)                          TTES 490
      DIMENSION A(1),B(1)                                               TTES 500
C                                                                       TTES 510
C     INITIALIZATION                                                    TTES 520
C                                                                       TTES 530
      NDF=0                                                             TTES 540
      ANS=0.0                                                           TTES 550
C                                                                       TTES 560
C     CALCULATE THE MEAN OF A                                           TTES 570
C                                                                       TTES 580
      AMEAN=0.0                                                         TTES 590
      DO 110 I=1,NA                                                     TTES 600
  110 AMEAN=AMEAN+A(I)                                                  TTES 610
      FNA=NA                                                            TTES 620
      AMEAN=AMEAN/FNA                                                   TTES 630
C                                                                       TTES 640
C     CALCULATE THE MEAN OF B                                           TTES 650
C                                                                       TTES 660
  115 BMEAN=0.0                                                         TTES 670
      DO 120 I=1,NB                                                     TTES 680
  120 BMEAN=BMEAN+B(I)                                                  TTES 690
      FNB=NB                                                            TTES 700
      BMEAN=BMEAN/FNB                                                   TTES 710
C                                                                       TTES 720
      IF(NOP-4) 122, 180, 200                                           TTES 730
  122 IF(NOP-1) 200, 135, 125                                           TTES 740
C                                                                       TTES 750
C     CALCULATE THE VARIANCE OF A                                       TTES 760
C                                                                       TTES 770
  125 SA2=0.0                                                           TTES 780
      DO 130 I=1,NA                                                     TTES 790
  130 SA2=SA2+(A(I)-AMEAN)**2                                           TTES 800
      SA2=SA2/(FNA-1.0)                                                 TTES 810
C                                                                       TTES 820
C     CALCULATE THE VARIANCE OF B                                       TTES 830
C                                                                       TTES 840
  135 SB2=0.0                                                           TTES 850
      DO 140 I=1,NB                                                     TTES 860
  140 SB2=SB2+(B(I)-BMEAN)**2                                           TTES 870
      SB2=SB2/(FNB-1.0)                                                 TTES 880
C                                                                       TTES 890
      GO TO (150,160,170), NOP                                          TTES 900
C                                                                       TTES 910
C        OPTION 1                                                       TTES 920
C                                                                       TTES 930
  150 ANS=((BMEAN-AMEAN)/SQRT(SB2))*SQRT(FNB)                           TTES 940
      NDF=NB-1                                                          TTES 950
      GO TO 200                                                         TTES 960
C                                                                       TTES 970
C        OPTION 2                                                       TTES 980
C                                                                       TTES 990
  160 NDF=NA+NB-2                                                       TTES1000
      FNDF=NDF                                                          TTES1010
      S=SQRT(((FNA-1.0)*SA2+(FNB-1.0)*SB2)/FNDF)                        TTES1020
      ANS=((BMEAN-AMEAN)/S)*(1.0/SQRT(1.0/FNA+1.0/FNB))                 TTES1030
      GO TO 200                                                         TTES1040
C                                                                       TTES1050
C        OPTION 3                                                       TTES1060
C                                                                       TTES1070
  170 ANS=(BMEAN-AMEAN)/SQRT(SA2/FNA+SB2/FNB)                           TTES1080
      A1=(SA2/FNA+SB2/FNB)**2                                           TTES1090
      A2=(SA2/FNA)**2/(FNA+1.0)+(SB2/FNB)**2/(FNB+1.0)                  TTES1100
      NDF=A1/A2-2.0+0.5                                                 TTES1110
      GO TO 200                                                         TTES1120
C                                                                       TTES1130
C        OPTION 4                                                       TTES1140
C                                                                       TTES1150
  180 SD=0.0                                                            TTES1160
      D=BMEAN-AMEAN                                                     TTES1170
      DO 190 I=1,NB                                                     TTES1180
  190 SD=SD+(B(I)-A(I)-D)**2                                            TTES1190
      SD=SQRT(SD/(FNB-1.0))                                             TTES1200
      ANS=(D/SD)*SQRT(FNB)                                              TTES1210
      NDF=NB-1                                                          TTES1220
C                                                                       TTES1230
  200 RETURN                                                            TTES1240
      END                                                               TTES1250
