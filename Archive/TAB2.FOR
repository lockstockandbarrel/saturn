C                                                                       TAB2  10
C     ..................................................................TAB2  20
C                                                                       TAB2  30
C        SUBROUTINE TAB2                                                TAB2  40
C                                                                       TAB2  50
C        PURPOSE                                                        TAB2  60
C           PERFORM A TWO-WAY CLASSIFICATION FOR TWO VARIABLES IN AN    TAB2  70
C           OBSERVATION MATRIX (OR A MATRIX SUBSET) OF THE FREQUENCY,   TAB2  80
C           PERCENT FREQUENCY, AND OTHER STATISTICS OVER GIVEN CLASS    TAB2  90
C           INTERVALS.                                                  TAB2 100
C                                                                       TAB2 110
C        USAGE                                                          TAB2 120
C           CALL TAB2(A,S,NOV,UBO,FREQ,PCT,STAT1,STAT2,NO,NV)           TAB2 130
C                                                                       TAB2 140
C        DESCRIPTION OF PARAMETERS                                      TAB2 150
C           A     - OBSERVATION MATRIX, NO BY NV                        TAB2 160
C           S     - INPUT VECTOR GIVING SUBSET OF A. ONLY THOSE         TAB2 170
C                   OBSERVATIONS WITH A CORRESPONDING NON-ZERO S(J) ARE TAB2 180
C                   CONSIDERED. VECTOR LENGTH IS NO.                    TAB2 190
C           NOV   - VARIABLES TO BE CROSS-TABULATED. NOV(1) IS VARIABLE TAB2 200
C                       1, NOV(2) IS VARIABLE 2. VECTOR LENGTH IS 2. NOVTAB2 210
C                       MUST BE GREATER THAN OR EQUAL TO 1 AND LESS THANTAB2 211
C                       OR EQUAL TO NV.                                 TAB2 212
C           UBO   - 3 BY 2 MATRIX GIVING LOWER LIMIT, NUMBER OF         TAB2 220
C                   INTERVALS, AND UPPER LIMIT OF BOTH VARIABLES TO BE  TAB2 230
C                   TABULATED (FIRST COLUMN FOR VARIABLE 1, SECOND      TAB2 240
C                   COLUMN FOR VARIABLE 2). IF LOWER LIMIT IS EQUAL TO  TAB2 250
C                   UPPER LIMIT FOR VARIABLE 1, THE PROGRAM USES THE    TAB2 260
C                   MINIMUM AND MAXIMUM VALUES ON EACH VARIABLE. NUMBER TAB2 270
C                   OF INTERVALS MUST INCLUDE TWO CELLS FOR UNDER AND   TAB2 280
C                   ABOVE LIMITS.                                       TAB2 290
C           FREQ  - OUTPUT MATRIX OF FREQUENCIES IN THE TWO-WAY         TAB2 300
C                   CLASSIFICATION. ORDER OF MATRIX IS INT1 BY INT2,    TAB2 310
C                   WHERE INT1 IS THE NUMBER OF INTERVALS OF VARIABLE 1 TAB2 320
C                   AND INT2 IS THE NUMBER OF INTERVALS OF VARIABLE 2.  TAB2 330
C                   INT1 AND INT2 MUST BE SPECIFIED IN THE SECOND       TAB2 340
C                   POSITION OF RESPECTIVE COLUMN OF UBO MATRIX.        TAB2 350
C           PCT   - OUTPUT MATRIX OF PERCENT FREQUENCIES, SAME ORDER    TAB2 360
C                   AS FREQ.                                            TAB2 370
C           STAT1 - OUTPUT MATRIX SUMMARIZING TOTALS, MEANS, AND        TAB2 380
C                   STANDARD DEVIATIONS FOR EACH CLASS INTERVAL OF      TAB2 390
C                   VARIABLE 1. ORDER OF MATRIX IS 3 BY INT1.           TAB2 400
C           STAT2 - SAME AS STAT1 BUT OVER VARIABLE 2. ORDER OF MATRIX  TAB2 410
C                   IS 3 BY INT2.                                       TAB2 420
C           NO    - NUMBER OF OBSERVATIONS. NO MUST BE GREATER THAN     TAB2 430
C                   OR EQUAL TO 1.                                      TAB2 431
C           NV    - NUMBER OF VARIABLES FOR EACH OBSERVATION. NV        TAB2 440
C                   MUST BE GREATER THAN OR EQUAL TO 1.                 TAB2 441
C                                                                       TAB2 450
C        REMARKS                                                        TAB2 460
C           IF S IS NULL, OUTPUT AREAS ARE SET TO ZERO                  TAB2 470
C                                                                       TAB2 480
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  TAB2 490
C           NONE                                                        TAB2 500
C                                                                       TAB2 510
C        METHOD                                                         TAB2 520
C           INTERVAL SIZES FOR BOTH VARIABLES ARE CALCULATED FROM THE   TAB2 530
C           GIVEN INFORMATION OR OPTIONALLY FROM THE MINIMUM AND MAXIMUMTAB2 540
C           VALUES. THE FREQUENCY AND PERCENT FREQUENCY MATRICES ARE    TAB2 550
C           DEVELOPED. MATRICES STAT1 AND STAT2 SUMMARIZING TOTALS,     TAB2 560
C           MEANS, AND STANDARD DEVIATIONS ARE THEN CALCULATED.         TAB2 570
C           THE DIVISOR FOR STANDARD DEVIATION IS ONE LESS THAN THE     TAB2 580
C           NUMBER OF OBSERVATIONS USED IN EACH CLASS INTERVAL.         TAB2 590
C                                                                       TAB2 600
C     ..................................................................TAB2 610
C                                                                       TAB2 620
      SUBROUTINE TAB2(A,S,NOV,UBO,FREQ,PCT,STAT1,STAT2,NO,NV)           TAB2 630
      DIMENSION A(1),S(1),NOV(2),UBO(3,2),FREQ(1),PCT(1),STAT1(1),      TAB2 640
     1STAT2(2),SINT(2)                                                  TAB2 650
      DIMENSION WBO(3,2)                                                TAB2 660
      DO 5 I=1,3                                                        TAB2 670
      DO 5 J=1,2                                                        TAB2 680
    5 WBO(I,J)=UBO(I,J)                                                 TAB2 690
C                                                                       TAB2 700
C        DETERMINE LIMITS                                               TAB2 710
C                                                                       TAB2 720
      DO 40 I=1,2                                                       TAB2 730
      IF(UBO(1,I)-UBO(3,I)) 40, 10, 40                                  TAB2 740
   10 VMIN=1.0E75                                                       TAB2 750
      VMAX=-1.0E75                                                      TAB2 760
      IJ=NO*(NOV(I)-1)                                                  TAB2 770
      DO 35 J=1,NO                                                      TAB2 780
      IJ=IJ+1                                                           TAB2 790
      IF(S(J)) 15,35,15                                                 TAB2 800
   15 IF(A(IJ)-VMIN) 20,25,25                                           TAB2 810
   20 VMIN=A(IJ)                                                        TAB2 820
   25 IF(A(IJ)-VMAX) 35,35,30                                           TAB2 830
   30 VMAX=A(IJ)                                                        TAB2 840
   35 CONTINUE                                                          TAB2 850
      UBO(1,I)=VMIN                                                     TAB2 860
      UBO(3,I)=VMAX                                                     TAB2 870
   40 CONTINUE                                                          TAB2 880
C                                                                       TAB2 890
C        CALCULATE INTERVAL SIZE                                        TAB2 900
C                                                                       TAB2 910
   45 DO 50 I=1,2                                                       TAB2 920
   50 SINT(I)=ABS((UBO(3,I)-UBO(1,I))/(UBO(2,I)-2.0))                   TAB2 930
C                                                                       TAB2 940
C        CLEAR OUTPUT AREAS                                             TAB2 950
C                                                                       TAB2 960
      INT1=UBO(2,1)                                                     TAB2 970
      INT2=UBO(2,2)                                                     TAB2 980
      INTT=INT1*INT2                                                    TAB2 990
      DO 55 I=1,INTT                                                    TAB21000
      FREQ(I)=0.0                                                       TAB21010
   55 PCT(I)=0.0                                                        TAB21020
      INTY=3*INT1                                                       TAB21030
      DO 60 I=1,INTY                                                    TAB21040
   60 STAT1(I)=0.0                                                      TAB21050
      INTZ=3*INT2                                                       TAB21060
      DO 65 I=1,INTZ                                                    TAB21070
   65 STAT2(I)=0.0                                                      TAB21080
C                                                                       TAB21090
C        TEST SUBSET VECTOR                                             TAB21100
C                                                                       TAB21110
      SCNT=0.0                                                          TAB21120
      INTY=INT1-1                                                       TAB21130
      INTX=INT2-1                                                       TAB21140
      IJ=NO*(NOV(1)-1)                                                  TAB21150
      IJX=NO*(NOV(2)-1)                                                 TAB21160
      DO 95 J=1,NO                                                      TAB21170
      IJ=IJ+1                                                           TAB21180
      IJX=IJX+1                                                         TAB21190
      IF(S(J)) 70,95,70                                                 TAB21200
   70 SCNT=SCNT+1.0                                                     TAB21210
C                                                                       TAB21220
C        CALCULATE FREQUENCIES                                          TAB21230
C                                                                       TAB21240
      TEMP1=UBO(1,1)-SINT(1)                                            TAB21250
      DO 75 IY=1,INTY                                                   TAB21260
      TEMP1=TEMP1+SINT(1)                                               TAB21270
      IF(A(IJ)-TEMP1) 80,75,75                                          TAB21280
   75 CONTINUE                                                          TAB21290
      IY=INT1                                                           TAB21300
   80 IYY=3*(IY-1)+1                                                    TAB21310
      STAT1(IYY)=STAT1(IYY)+A(IJ)                                       TAB21320
      IYY=IYY+1                                                         TAB21330
      STAT1(IYY)=STAT1(IYY)+1.0                                         TAB21340
      IYY=IYY+1                                                         TAB21350
      STAT1(IYY)=STAT1(IYY)+A(IJ)*A(IJ)                                 TAB21360
      TEMP2=UBO(1,2)-SINT(2)                                            TAB21370
      DO 85 IX=1,INTX                                                   TAB21380
      TEMP2=TEMP2+SINT(2)                                               TAB21390
      IF(A(IJX)-TEMP2) 90,85,85                                         TAB21400
   85 CONTINUE                                                          TAB21410
      IX=INT2                                                           TAB21420
   90 IJF=INT1*(IX-1)+IY                                                TAB21430
      FREQ(IJF)=FREQ(IJF)+1.0                                           TAB21440
      IX=3*(IX-1)+1                                                     TAB21450
      STAT2(IX)=STAT2(IX)+A(IJX)                                        TAB21460
      IX=IX+1                                                           TAB21470
      STAT2(IX)=STAT2(IX)+1.0                                           TAB21480
      IX=IX+1                                                           TAB21490
      STAT2(IX)=STAT2(IX)+A(IJX)*A(IJX)                                 TAB21500
   95 CONTINUE                                                          TAB21510
      IF (SCNT)98,151,98                                                TAB21511
C                                                                       TAB21520
C        CALCULATE PERCENT FREQUENCIES                                  TAB21530
C                                                                       TAB21540
   98 DO 100 I=1,INTT                                                   TAB21550
  100 PCT(I)=FREQ(I)*100.0/SCNT                                         TAB21560
C                                                                       TAB21570
C        CALCULATE TOTALS, MEANS, STANDARD DEVIATIONS                   TAB21580
C                                                                       TAB21590
      IXY=-1                                                            TAB21600
      DO 120 I=1,INT1                                                   TAB21610
      IXY=IXY+3                                                         TAB21620
      ISD=IXY+1                                                         TAB21630
      TEMP1=STAT1(IXY)                                                  TAB21640
      SUM=STAT1(IXY-1)                                                  TAB21650
      IF(TEMP1-1.0) 120,105,110                                         TAB21660
  105 STAT1(ISD)=0.0                                                    TAB21670
      GO TO 115                                                         TAB21680
  110 STAT1(ISD)=SQRT(ABS((STAT1(ISD)-SUM*SUM/TEMP1)/(TEMP1-1.0)))      TAB21690
  115 STAT1(IXY)=SUM/TEMP1                                              TAB21700
  120 CONTINUE                                                          TAB21710
      IXX=-1                                                            TAB21720
      DO 140 I=1,INT2                                                   TAB21730
      IXX=IXX+3                                                         TAB21740
      ISD=IXX+1                                                         TAB21750
      TEMP2=STAT2(IXX)                                                  TAB21760
      SUM=STAT2(IXX-1)                                                  TAB21770
      IF(TEMP2-1.0) 140,125,130                                         TAB21780
  125 STAT2(ISD)=0.0                                                    TAB21790
      GO TO 135                                                         TAB21800
  130 STAT2(ISD)=SQRT(ABS((STAT2(ISD)-SUM*SUM/TEMP2)/(TEMP2-1.0)))      TAB21810
  135 STAT2(IXX)=SUM/TEMP2                                              TAB21820
  140 CONTINUE                                                          TAB21830
      DO 150 I=1,3                                                      TAB21840
      DO 150 J=1,2                                                      TAB21850
  150 UBO(I,J)=WBO(I,J)                                                 TAB21860
  151 RETURN                                                            TAB21870
      END                                                               TAB21880
