C                                                                       SIGN  10
C     ..................................................................SIGN  20
C                                                                       SIGN  30
C        SUBROUTINE SIGNT                                               SIGN  40
C                                                                       SIGN  50
C        PURPOSE                                                        SIGN  60
C           TO PERFORM A NON-PARAMETRIC SIGN TEST, GIVEN TWO SETS OF    SIGN  70
C           MATCHED OBSERVATIONS.  IT TESTS THE NULL HYPOTHESIS THAT THESIGN  80
C           DIFFERENCES BETWEEN EACH PAIR OF MATCHED OBSERVATIONS HAS A SIGN  90
C           MEDIAN EQUAL TO ZERO.                                       SIGN 100
C                                                                       SIGN 110
C        USAGE                                                          SIGN 120
C           CALL SIGNT (N,A,B,K,M,P,IE)                                 SIGN 130
C                                                                       SIGN 140
C        DESCRIPTION OF PARAMETERS                                      SIGN 150
C           N - NUMBER OF OBSERVATIONS IN SETS A AND B                  SIGN 160
C           A - INPUT VECTOR OF LENGTH N CONTAINING DATA FROM THE FIRST SIGN 170
C               SAMPLE, A                                               SIGN 180
C           B - INPUT VECTOR OF LENGTH N CONTAINING DATA FROM THE SECONDSIGN 190
C               SAMPLE, B                                               SIGN 200
C           K - OUTPUT VARIABLE CONTAINING THE NUMBER OF PAIRS OF       SIGN 210
C               OBSERVATIONS FROM THE TWO SAMPLES WHOSE DIFFERENCES ARE SIGN 220
C               NON-ZERO                                                SIGN 230
C           M - OUTPUT VARIABLE CONTAINING THE NUMBER OF PLUS OR MINUS  SIGN 240
C               DIFFERENCES, WHICHEVER IS FEWER.                        SIGN 250
C           P - COMPUTED PROBABILITY OF AS FEW AS M NUMBER OF PAIRS     SIGN 260
C               HAVING THE SAME SIGN, ASSUMING THAT THE SAMPLES CAME    SIGN 270
C               FROM THE SAME POPULATION.                               SIGN 280
C           IE- 0, IF THERE IS NO ERROR.                                SIGN 290
C               1, IF K IS ZERO.  IN THIS CASE, P IS SET TO 1.0 AND     SIGN 300
C               M TO 0.                                                 SIGN 310
C                                                                       SIGN 320
C        REMARKS                                                        SIGN 330
C           IF K IS LESS THAN OR EQUAL TO 25, THE PROBABILITY WILL BE   SIGN 340
C           COMPUTED USING THE BINOMIAL DISTRIBUTION.  IF K IS GREATER  SIGN 350
C           THAN 25, THE PROBABILITY WILL BE COMPUTED USING THE NORMAL  SIGN 360
C           APPROXIMATION TO THE BINOMIAL DISTRIBUTION.                 SIGN 370
C           P COMPUTED IS THE PROBABILITY FOR A ONE-TAILED TEST.  THUS, SIGN 380
C           FOR A TWO TAILED TEST, DOUBLE THE VALUE FOR P.              SIGN 390
C                                                                       SIGN 400
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  SIGN 410
C           NDTR                                                        SIGN 420
C                                                                       SIGN 430
C        METHOD                                                         SIGN 440
C           REFER TO DIXON AND MASSEY, INTRODUCTION TO STATISTICAL      SIGN 450
C           ANALYSIS (MCGRAW-HILL, 1957).                               SIGN 460
C                                                                       SIGN 470
C     ..................................................................SIGN 480
C                                                                       SIGN 490
      SUBROUTINE SIGNT (N,A,B,K,M,P,IE)                                 SIGN 500
C                                                                       SIGN 510
      DIMENSION A(1),B(1)                                               SIGN 520
      DOUBLE PRECISION FN,FD                                            SIGN 530
C                                                                       SIGN 540
C        INITIALIZATION                                                 SIGN 550
C                                                                       SIGN 560
      IE=0                                                              SIGN 570
      K=0                                                               SIGN 580
      MPLUS=0                                                           SIGN 590
      MMINS=0                                                           SIGN 600
C                                                                       SIGN 610
C        FIND (+) OR (-) DIFFERENCE                                     SIGN 620
C                                                                       SIGN 630
      DO 40 I=1,N                                                       SIGN 640
      D=A(I)-B(I)                                                       SIGN 650
      IF(D) 20, 40, 30                                                  SIGN 660
C                                                                       SIGN 670
C        (-) DIFFERENCE                                                 SIGN 680
C                                                                       SIGN 690
   20 K=K+1                                                             SIGN 700
      MMINS=MMINS+1                                                     SIGN 710
      GO TO 40                                                          SIGN 720
C                                                                       SIGN 730
C        (+) DIFFERENCE                                                 SIGN 740
C                                                                       SIGN 750
   30 K=K+1                                                             SIGN 760
      MPLUS=MPLUS+1                                                     SIGN 770
C                                                                       SIGN 780
   40 CONTINUE                                                          SIGN 790
      IF(K) 41,41,42                                                    SIGN 800
   41 IE=1                                                              SIGN 810
      P=1.0                                                             SIGN 820
      M=0                                                               SIGN 830
      GO TO 95                                                          SIGN 840
   42 FK=K                                                              SIGN 850
C                                                                       SIGN 860
C        FIND THE NUMBER OF FEWER SIGNS                                 SIGN 870
C                                                                       SIGN 880
      IF(MPLUS-MMINS) 45, 45, 50                                        SIGN 890
   45 M=MPLUS                                                           SIGN 900
      GO TO 55                                                          SIGN 910
   50 M=MMINS                                                           SIGN 920
C                                                                       SIGN 930
C        TEST WHETHER K IS GREATER THAN 25                              SIGN 940
C                                                                       SIGN 950
   55 IF(K-25) 60, 60, 77                                               SIGN 960
C                                                                       SIGN 970
C        K IS LESS THAN OR EQUAL TO 25                                  SIGN 980
C                                                                       SIGN 990
   60 P=1.0                                                             SIGN1000
      IF(M) 75, 75, 65                                                  SIGN1010
   65 FN=1.0                                                            SIGN1020
      FD=1.0                                                            SIGN1030
      DO 70 I=1,M                                                       SIGN1040
      FI=I                                                              SIGN1050
      FN=FN*(FK-(FI-1.0))                                               SIGN1060
      FD=FD*FI                                                          SIGN1070
   70 P=P+FN/FD                                                         SIGN1080
C                                                                       SIGN1090
   75 P=P/(2.0**K)                                                      SIGN1100
      GO TO 95                                                          SIGN1110
C                                                                       SIGN1120
C        K IS GREATER THAN 25.  COMPUTE MEAN, STANDARD DEVIATION, AND Z SIGN1130
C                                                                       SIGN1140
   77 U=0.5*FK                                                          SIGN1150
      S=0.5*SQRT(FK)                                                    SIGN1160
      FM=M                                                              SIGN1170
      IF(FM-U) 80, 85, 85                                               SIGN1180
   80 CON=0.5                                                           SIGN1190
      GO TO 90                                                          SIGN1200
   85 CON=0.0                                                           SIGN1210
   90 Z=(FM+CON-U)/S                                                    SIGN1220
C                                                                       SIGN1230
C        COMPUTE P ASSOCIATED WITH THE VALUE AS EXTREME AS Z            SIGN1240
C                                                                       SIGN1250
      CALL NDTR (Z,P,D)                                                 SIGN1260
C                                                                       SIGN1270
   95 RETURN                                                            SIGN1280
      END                                                               SIGN1290
