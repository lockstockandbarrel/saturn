C                                                                       TAB1  10
C     ..................................................................TAB1  20
C                                                                       TAB1  30
C        SUBROUTINE TAB1                                                TAB1  40
C                                                                       TAB1  50
C        PURPOSE                                                        TAB1  60
C           TABULATE FOR ONE VARIABLE IN AN OBSERVATION MATRIX (OR A    TAB1  70
C           MATRIX SUBSET), THE FREQUENCY AND PERCENT FREQUENCY OVER    TAB1  80
C           GIVEN CLASS INTERVALS. IN ADDITION, CALCULATE FOR THE SAME  TAB1  90
C           VARIABLE THE TOTAL, AVERAGE, STANDARD DEVIATION, MINIMUM,   TAB1 100
C           AND MAXIMUM.                                                TAB1 110
C                                                                       TAB1 120
C        USAGE                                                          TAB1 130
C           CALL TAB1(A,S,NOVAR,UBO,FREQ,PCT,STATS,NO,NV)               TAB1 140
C                                                                       TAB1 150
C        DESCRIPTION OF PARAMETERS                                      TAB1 160
C           A     - OBSERVATION MATRIX, NO BY NV                        TAB1 170
C           S     - INPUT VECTOR GIVING SUBSET OF A. ONLY THOSE         TAB1 180
C                   OBSERVATIONS WITH A CORRESPONDING NON-ZERO S(J) ARE TAB1 190
C                   CONSIDERED. VECTOR LENGTH IS NO.                    TAB1 200
C           NOVAR - THE VARIABLE TO BE TABULATED. NOVAR MUST BE GREATER TAB1 210
C                   THAN OR EQUAL TO 1 AND LESS THAN OR EQUAL TO NV.    TAB1 220
C                   AND UPPER LIMIT OF VARIABLE TO BE TABULATED         TAB1 230
C                   IN UBO(1), UBO(2) AND UBO(3) RESPECTIVELY. IF       TAB1 240
C                   LOWER LIMIT IS EQUAL TO UPPER LIMIT, THE PROGRAM    TAB1 250
C                   USES THE MINIMUM AND MAXIMUM VALUES OF THE VARIABLE.TAB1 260
C                   NUMBER OF INTERVALS, UBO(2), MUST INCLUDE TWO CELLS TAB1 270
C                   FOR VALUES UNDER AND ABOVE LIMITS. VECTOR LENGTH    TAB1 280
C                   IS 3.                                               TAB1 290
C           FREQ  - OUTPUT VECTOR OF FREQUENCIES. VECTOR LENGTH IS      TAB1 300
C                   UBO(2).                                             TAB1 310
C           PCT   - OUTPUT VECTOR OF RELATIVE FREQUENCIES. VECTOR       TAB1 320
C                   LENGTH IS UBO(2).                                   TAB1 330
C           STATS - OUTPUT VECTOR OF SUMMARY STATISTICS, I.E., TOTAL,   TAB1 340
C                   AVERAGE, STANDARD DEVIATION, MINIMUM AND MAXIMUM.   TAB1 350
C                   VECTOR LENGTH IS 5. IF S IS NULL, THEN TOTAL,AVERAGETAB1 360
C                   AND STANDARD DEVIATION = 0, MIN=1.E75 AND MAX=-1.E75TAB1 361
C           NO    - NUMBER OF OBSERVATIONS. NO MUST BE > OR = TO 1      TAB1 370
C           NV    - NUMBER OF VARIABLES FOR EACH OBSERVATION. NV MUST   TAB1 380
C                   BE GREATER THAN OR EQUAL TO 1.                      TAB1 381
C                                                                       TAB1 390
C        REMARKS                                                        TAB1 400
C           NONE                                                        TAB1 410
C                                                                       TAB1 420
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  TAB1 430
C           NONE                                                        TAB1 440
C                                                                       TAB1 450
C        METHOD                                                         TAB1 460
C           THE INTERVAL SIZE IS CALCULATED FROM THE GIVEN INFORMATION  TAB1 470
C           OR OPTIONALLY FROM THE MINIMUM AND MAXIMUM VALUES FOR       TAB1 480
C           VARIABLE NOVAR. THE FREQUENCIES AND PERCENT FREQUENCIES ARE TAB1 490
C           THEN CALCULATED ALONG WITH SUMMARY STATISTICS.              TAB1 500
C           THE DIVISOR FOR STANDARD DEVIATION IS ONE LESS THAN THE     TAB1 510
C           NUMBER OF OBSERVATIONS USED.                                TAB1 520
C                                                                       TAB1 530
C     ..................................................................TAB1 540
C                                                                       TAB1 550
      SUBROUTINE TAB1(A,S,NOVAR,UBO,FREQ,PCT,STATS,NO,NV)               TAB1 560
      DIMENSION A(1),S(1),UBO(1),FREQ(1),PCT(1),STATS(1)                TAB1 570
      DIMENSION WBO(3)                                                  TAB1 580
      DO 5 I=1,3                                                        TAB1 590
    5 WBO(I)=UBO(I)                                                     TAB1 600
C                                                                       TAB1 610
C        CALCULATE MIN AND MAX                                          TAB1 620
C                                                                       TAB1 630
      VMIN=1.0E75                                                       TAB1 640
      VMAX=-1.0E75                                                      TAB1 650
      IJ=NO*(NOVAR-1)                                                   TAB1 660
      DO 30 J=1,NO                                                      TAB1 670
      IJ=IJ+1                                                           TAB1 680
      IF(S(J)) 10,30,10                                                 TAB1 690
   10 IF(A(IJ)-VMIN) 15,20,20                                           TAB1 700
   15 VMIN=A(IJ)                                                        TAB1 710
   20 IF(A(IJ)-VMAX) 30,30,25                                           TAB1 720
   25 VMAX=A(IJ)                                                        TAB1 730
   30 CONTINUE                                                          TAB1 740
      STATS(4)=VMIN                                                     TAB1 750
      STATS(5)=VMAX                                                     TAB1 760
C                                                                       TAB1 770
C        DETERMINE LIMITS                                               TAB1 780
C                                                                       TAB1 790
      IF(UBO(1)-UBO(3)) 40,35,40                                        TAB1 800
   35 UBO(1)=VMIN                                                       TAB1 810
      UBO(3)=VMAX                                                       TAB1 820
   40 INN=UBO(2)                                                        TAB1 830
C                                                                       TAB1 840
C        CLEAR OUTPUT AREAS                                             TAB1 850
C                                                                       TAB1 860
      DO 45 I=1,INN                                                     TAB1 870
      FREQ(I)=0.0                                                       TAB1 880
   45 PCT(I)=0.0                                                        TAB1 890
      DO 50 I=1,3                                                       TAB1 900
   50 STATS(I)=0.0                                                      TAB1 910
C                                                                       TAB1 920
C        CALCULATE INTERVAL SIZE                                        TAB1 930
C                                                                       TAB1 940
      SINT=ABS((UBO(3)-UBO(1))/(UBO(2)-2.0))                            TAB1 950
C                                                                       TAB1 960
C        TEST SUBSET VECTOR                                             TAB1 970
C                                                                       TAB1 980
      SCNT=0.0                                                          TAB1 990
      IJ=NO*(NOVAR-1)                                                   TAB11000
      DO 75 J=1,NO                                                      TAB11010
      IJ=IJ+1                                                           TAB11020
      IF(S(J)) 55,75,55                                                 TAB11030
   55 SCNT=SCNT+1.0                                                     TAB11040
C                                                                       TAB11050
C        DEVELOP TOTAL AND FREQUENCIES                                  TAB11060
C                                                                       TAB11070
      STATS(1)=STATS(1)+A(IJ)                                           TAB11080
      STATS(3)=STATS(3)+A(IJ)*A(IJ)                                     TAB11090
      TEMP=UBO(1)-SINT                                                  TAB11100
      INTX=INN-1                                                        TAB11110
      DO 60 I=1,INTX                                                    TAB11120
      TEMP=TEMP+SINT                                                    TAB11130
      IF(A(IJ)-TEMP) 70,60,60                                           TAB11140
   60 CONTINUE                                                          TAB11150
      IF(A(IJ)-TEMP) 75,65,65                                           TAB11160
   65 FREQ(INN)=FREQ(INN)+1.0                                           TAB11170
      GO TO 75                                                          TAB11180
   70 FREQ(I)=FREQ(I)+1.0                                               TAB11190
   75 CONTINUE                                                          TAB11200
      IF (SCNT)79,105,79                                                TAB11201
C                                                                       TAB11210
C        CALCULATE RELATIVE FREQUENCIES                                 TAB11220
C                                                                       TAB11230
   79 DO 80 I=1,INN                                                     TAB11240
   80 PCT(I)=FREQ(I)*100.0/SCNT                                         TAB11250
C                                                                       TAB11260
C        CALCULATE MEAN AND STANDARD DEVIATION                          TAB11270
C                                                                       TAB11280
      IF(SCNT-1.0) 85,85,90                                             TAB11290
   85 STATS(2)=STATS(1)                                                 TAB11300
      STATS(3)=0.0                                                      TAB11310
      GO TO 95                                                          TAB11320
   90 STATS(2)=STATS(1)/SCNT                                            TAB11330
      STATS(3)=SQRT(ABS((STATS(3)-STATS(1)*STATS(1)/SCNT)/(SCNT-1.0)))  TAB11340
   95 DO 100 I=1,3                                                      TAB11350
  100 UBO(I)=WBO(I)                                                     TAB11360
  105 RETURN                                                            TAB11370
      END                                                               TAB11380
