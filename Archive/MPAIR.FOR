C                                                                       MPAI  10
C     ..................................................................MPAI  20
C                                                                       MPAI  30
C        SUBROUTINE MPAIR                                               MPAI  40
C                                                                       MPAI  50
C        PURPOSE                                                        MPAI  60
C           PERFORM THE WILCOXON MATCHED-PAIRS SIGNED-RANKS TEST, GIVEN MPAI  70
C           TWO VECTORS OF N OBSERVATIONS OF THE MATCHED SAMPLES.       MPAI  80
C                                                                       MPAI  90
C        USAGE                                                          MPAI 100
C           CALL MPAIR (N,A,B,K,T,Z,P,D,E,L,IE)                         MPAI 110
C                                                                       MPAI 120
C        DESCRIPTION OF PARAMETERS                                      MPAI 130
C           N - NUMBER OF OBSERVATIONS IN THE VECTORS A AND B           MPAI 140
C           A - INPUT VECTOR OF LENGTH N CONTAINING DATA FROM THE FIRST MPAI 150
C               SAMPLE                                                  MPAI 160
C           B - INPUT VECTOR OF LENGTH N CONTAINING DATA FROM THE SECONDMPAI 170
C               SAMPLE                                                  MPAI 180
C           K - OUTPUT VARIABLE CONTAINING THE NUMBER OF PAIRS OF THE   MPAI 190
C               MATCHED SAMPLES WHOSE DIFFERENCES ARE NON ZERO (0)      MPAI 200
C           T - OUTPUT VARIABLE CONTAINING THE SUM OF THE RANKS OF PLUS MPAI 210
C               OR MINUS DIFFERENCES, WHICHEVER IS SMALLER              MPAI 220
C           Z - VALUE OF THE STANDARDIZED NORMAL SCORE COMPUTED FOR THE MPAI 230
C               WILCOXON MATCHED-PAIRS SIGNED-RANKS TEST                MPAI 240
C           P - COMPUTED PROBABILITY OF OBTAINING A VALUE OF Z AS       MPAI 250
C               EXTREME AS THE ONE FOUND BY THE TEST                    MPAI 260
C           D - WORKING VECTOR OF LENGTH N                              MPAI 270
C           E - WORKING VECTOR OF LENGTH N                              MPAI 280
C           L - WORKING VECTOR OF LENGTH N                              MPAI 290
C           IE- 1, IF SAMPLES A AND B ARE IDENTICAL.                    MPAI 300
C               0 OTHERWISE.  IF IE=1, THEN T=P=0, AND Z=-10**75        MPAI 310
C                                                                       MPAI 320
C        REMARKS                                                        MPAI 330
C           THE COMPUTED PROBABILTY IS FOR A ONE-TAILED TEST.           MPAI 340
C           MULTIPLYING P BY 2 WILL GIVE THE VALUE FOR A TWO-TAILED     MPAI 350
C           TEST.                                                       MPAI 360
C                                                                       MPAI 370
C        SUBROUTINES AND FUNCTIONS SUBPROGRAMS REQUIRED                 MPAI 380
C           RANK                                                        MPAI 390
C           NDTR                                                        MPAI 400
C                                                                       MPAI 410
C        METHOD                                                         MPAI 420
C           REFER TO DIXON AND MASSEY, AN INTRODUCTION TO STATISTICAL   MPAI 430
C           ANALYSIS (MC GRAW-HILL, 1957)                               MPAI 440
C                                                                       MPAI 450
C     ..................................................................MPAI 460
C                                                                       MPAI 470
      SUBROUTINE MPAIR (N,A,B,K,T,Z,P,D,E,L,IE)                         MPAI 480
C                                                                       MPAI 490
      DIMENSION A(1),B(1),D(1),E(1),L(1)                                MPAI 500
C                                                                       MPAI 510
      IE=0                                                              MPAI 520
      K=N                                                               MPAI 530
C                                                                       MPAI 540
C        FIND DIFFERENCES OF MATCHED-PAIRS                              MPAI 550
C                                                                       MPAI 560
      BIG=0.0                                                           MPAI 570
      DO 55 I=1,N                                                       MPAI 580
      DIF=A(I)-B(I)                                                     MPAI 590
      IF(DIF) 10, 20, 30                                                MPAI 600
C                                                                       MPAI 610
C        DIFFERENCE HAS A NEGATIVE SIGN (-)                             MPAI 620
C                                                                       MPAI 630
   10 L(I)=1                                                            MPAI 640
      GO TO 40                                                          MPAI 650
C                                                                       MPAI 660
C        DIFFERENCE IS ZERO (0)                                         MPAI 670
C                                                                       MPAI 680
   20 L(I)=2                                                            MPAI 690
      K=K-1                                                             MPAI 700
      GO TO 40                                                          MPAI 710
C                                                                       MPAI 720
C        DIFFERENCE HAS A POSITIVE SIGN (+)                             MPAI 730
C                                                                       MPAI 740
   30 L(I)=3                                                            MPAI 750
C                                                                       MPAI 760
   40 DIF= ABS(DIF)                                                     MPAI 770
      IF(BIG-DIF) 45, 50, 50                                            MPAI 780
   45 BIG=DIF                                                           MPAI 790
   50 D(I)=DIF                                                          MPAI 800
C                                                                       MPAI 810
   55 CONTINUE                                                          MPAI 820
      IF(K) 57,57,59                                                    MPAI 830
   57 IE=1                                                              MPAI 840
      T=0.0                                                             MPAI 850
      Z=-1.0E75                                                         MPAI 860
      P=0                                                               MPAI 870
      GO TO 100                                                         MPAI 880
C                                                                       MPAI 890
C        STORE A LARGE VALUE IN PLACE OF 0 DIFFERENCE IN ORDER TO       MPAI 900
C        ASSIGN A LARGE RANK (LARGER THAN K), SO THAT ABSOLUTE VALUES   MPAI 910
C        OF SIGNED DIFFERENCES WILL BE PROPERLY RANKED                  MPAI 920
C                                                                       MPAI 930
   59 BIG=BIG*2.0                                                       MPAI 940
      DO 65 I=1,N                                                       MPAI 950
      IF(L(I)-2) 65, 60, 65                                             MPAI 960
   60 D(I)=BIG                                                          MPAI 970
   65 CONTINUE                                                          MPAI 980
C                                                                       MPAI 990
      CALL RANK (D,E,N)                                                 MPAI1000
C                                                                       MPAI1010
C        FIND SUMS OF RANKS OF (+) DIFFERENCES AND (-) DIFFERENCES      MPAI1020
C                                                                       MPAI1030
      SUMP=0.0                                                          MPAI1040
      SUMM=0.0                                                          MPAI1050
      DO 80 I=1,N                                                       MPAI1060
      IF(L(I)-2) 70, 80, 75                                             MPAI1070
   70 SUMM=SUMM+E(I)                                                    MPAI1080
      GO TO 80                                                          MPAI1090
   75 SUMP=SUMP+E(I)                                                    MPAI1100
   80 CONTINUE                                                          MPAI1110
C                                                                       MPAI1120
C        SET T = SMALLER SUM                                            MPAI1130
C                                                                       MPAI1140
      IF(SUMP-SUMM) 85, 85, 90                                          MPAI1150
   85 T=SUMP                                                            MPAI1160
      GO TO 95                                                          MPAI1170
   90 T=SUMM                                                            MPAI1180
C                                                                       MPAI1190
C        COMPUTE MEAN, STANDARD DEVIATION, AND Z                        MPAI1200
C                                                                       MPAI1210
   95 FK=K                                                              MPAI1220
      U=FK*(FK+1.0)/4.0                                                 MPAI1230
      S= SQRT((FK*(FK+1.0)*(2.0*FK+1.0))/24.0)                          MPAI1240
      Z=(T-U)/S                                                         MPAI1250
C                                                                       MPAI1260
C        COMPUTE THE PROBABILITY OF A VALUE AS EXTREME AS Z             MPAI1270
C                                                                       MPAI1280
      CALL NDTR (Z,P,BIG)                                               MPAI1290
C                                                                       MPAI1300
  100 RETURN                                                            MPAI1310
      END                                                               MPAI1320
