C                                                                       STPR  10
C     ..................................................................STPR  20
C                                                                       STPR  30
C        SUBROUTINE STPRG                                               STPR  40
C                                                                       STPR  50
C        PURPOSE                                                        STPR  60
C           TO PERFORM A STEPWISE MULTIPLE REGRESSION ANALYSIS FOR A    STPR  70
C           DEPENDENT VARIABLE AND A SET OF INDEPENDENT VARIABLES.  AT  STPR  80
C           EACH STEP, THE VARIABLE ENTERED INTO THE REGRESSION EQUATIONSTPR  90
C           IS THAT WHICH EXPLAINS THE GREATEST AMOUNT OF VARIANCE      STPR 100
C           BETWEEN IT AND THE DEPENDENT VARIABLE (I.E. THE VARIABLE    STPR 110
C           WITH THE HIGHEST PARTIAL CORRELATION WITH THE DEPENDENT     STPR 120
C           VARIABLE).  ANY VARIABLE CAN BE DESIGNATED AS THE DEPENDENT STPR 130
C           VARIABLE.  ANY INDEPENDENT VARIABLE CAN BE FORCED INTO OR   STPR 140
C           DELETED FROM THE REGRESSION EQUATION, IRRESPECTIVE OF ITS   STPR 150
C           CONTRIBUTION TO THE EQUATION.                               STPR 160
C                                                                       STPR 170
C        USAGE                                                          STPR 180
C           CALL STPRG (M,N,D,XBAR,IDX,PCT,NSTEP,ANS,L,B,S,T,LL,IER)    STPR 190
C                                                                       STPR 200
C        DESCRIPTION OF PARAMETERS                                      STPR 210
C           M    - TOTAL NUMBER OF VARIABLES IN DATA MATRIX             STPR 220
C           N    - NUMBER OF OBSERVATIONS                               STPR 230
C           D    - INPUT MATRIX (M X M) OF SUMS OF CROSS-PRODUCTS OF    STPR 240
C                  DEVIATIONS FROM MEAN.  THIS MATRIX WILL BE DESTROYED.STPR 250
C           XBAR - INPUT VECTOR OF LENGTH M OF MEANS                    STPR 260
C           IDX  - INPUT VECTOR OF LENGTH M HAVING ONE OF THE FOLLOWING STPR 270
C                  CODES FOR EACH VARIABLE.                             STPR 280
C                    0 - INDEPENDENT VARIABLE AVAILABLE FOR SELECTION   STPR 290
C                    1 - INDEPENDENT VARIABLE TO BE FORCED INTO THE     STPR 300
C                        REGRESSION EQUATION                            STPR 310
C                    2 - VARIABLE NOT TO BE CONSIDERED IN THE EQUATION  STPR 320
C                    3 - DEPENDENT VARIABLE                             STPR 330
C                  THIS VECTOR WILL BE DESTROYED                        STPR 340
C           PCT  - A CONSTANT VALUE INDICATING THE PROPORTION OF THE    STPR 350
C                  TOTAL VARIANCE TO BE EXPLAINED BY ANY INDEPENDENT    STPR 360
C                  VARIABLE.  THOSE INDEPENDENT VARIABLES WHICH FALL    STPR 370
C                  BELOW THIS PROPORTION WILL NOT ENTER THE REGRESSION  STPR 380
C                  EQUATION.  TO ENSURE THAT ALL VARIABLES ENTER THE    STPR 390
C                  EQUATION, SET PCT = 0.0.                             STPR 400
C           NSTEP- OUTPUT VECTOR OF LENGTH 5 CONTAINING THE FOLLOWING   STPR 410
C                  INFORMATION                                          STPR 420
C                     NSTEP(1)- THE NUMBER OF THE DEPENDENT VARIABLE    STPR 430
C                     NSTEP(2)- NUMBER OF VARIABLES FORCED INTO THE     STPR 440
C                               REGRESSION EQUATION                     STPR 450
C                     NSTEP(3)- NUMBER OF VARIABLE DELETED FROM THE     STPR 460
C                               EQUATION                                STPR 470
C                     NSTEP(4)- THE NUMBER OF THE LAST STEP             STPR 480
C                     NSTEP(5)- THE NUMBER OF THE LAST VARIABLE ENTERED STPR 490
C           ANS  - OUTPUT VECTOR OF LENGTH 11 CONTAINING THE FOLLOWING  STPR 500
C                  INFORMATION FOR THE LAST STEP                        STPR 510
C                     ANS(1)- SUM OF SQUARES REDUCED BY THIS STEP       STPR 520
C                     ANS(2)- PROPORTION OF TOTAL SUM OF SQUARES REDUCEDSTPR 530
C                     ANS(3)- CUMULATIVE SUM OF SQUARES REDUCED UP TO   STPR 540
C                             THIS STEP                                 STPR 550
C                     ANS(4)- CUMULATIVE PROPORTION OF TOTAL SUM OF     STPR 560
C                             SQUARES REDUCED                           STPR 570
C                     ANS(5)- SUM OF SQUARES OF THE DEPENDENT VARIABLE  STPR 580
C                     ANS(6)- MULTIPLE CORRELATION COEFFICIENT          STPR 590
C                     ANS(7)- F RATIO FOR SUM OF SQUARES DUE TO         STPR 600
C                             REGRESSION                                STPR 610
C                     ANS(8)- STANDARD ERROR OF THE ESTIMATE (RESIDUAL  STPR 620
C                             MEAN SQUARE)                              STPR 630
C                     ANS(9)- INTERCEPT CONSTANT                        STPR 640
C                     ANS(10)-MULTIPLE CORRELATION COEFFICIENT ADJUSTED STPR 650
C                             FOR DEGREES OF FREEDOM.                   STPR 660
C                     ANS(11)-STANDARD ERROR OF THE ESTIMATE ADJUSTED   STPR 670
C                             FOR DEGREES OF FREEDOM.                   STPR 680
C           L    - OUTPUT VECTOR OF LENGTH K, WHERE K IS THE NUMBER OF  STPR 690
C                  INDEPENDENT VARIABLES IN THE REGRESSION EQUATION.    STPR 700
C                  THIS VECTOR CONTAINS THE NUMBERS OF THE INDEPENDENT  STPR 710
C                  VARIABLES IN THE EQUATION.                           STPR 720
C           B    - OUTPUT VECTOR OF LENGTH K, CONTAINING THE PARTIAL    STPR 730
C                  REGRESSION COEFFICIENTS CORRESPONDING TO THE         STPR 740
C                  VARIABLES IN VECTOR L.                               STPR 750
C           S    - OUTPUT VECTOR OF LENGTH K, CONTAINING THE STANDARD   STPR 760
C                  ERRORS OF THE PARTIAL REGRESSION COEFFICIENTS,       STPR 770
C                  CORRESPONDING TO THE VARIABLES IN VECTOR L.          STPR 780
C           T    - OUTPUT VECTOR OF LENGTH K, CONTAINING THE COMPUTED   STPR 790
C                  T-VALUES CORRESPONDING TO THE VARIABLES IN VECTOR L. STPR 800
C           LL   - WORKING VECTOR OF LENGTH M                           STPR 810
C           IER  - 0, IF THERE IS NO ERROR.                             STPR 820
C                  1, IF RESIDUAL SUM OF SQUARES IS NEGATIVE OR IF THE  STPR 830
C                  PIVOTAL ELEMENT IN THE STEPWISE INVERSION PROCESS IS STPR 840
C                  ZERO.  IN THIS CASE, THE VARIABLE WHICH CAUSES THIS  STPR 850
C                  ERROR IS NOT ENTERED IN THE REGRESSION, THE RESULT   STPR 860
C                  PRIOR TO THIS STEP IS RETAINED, AND THE CURRENT      STPR 870
C                  SELECTION IS TERMINATED.                             STPR 880
C                                                                       STPR 890
C        REMARKS                                                        STPR 900
C           THE NUMBER OF DATA POINTS MUST BE AT LEAST GREATER THAN THE STPR 910
C           NUMBER OF INDEPENDENT VARIABLES PLUS ONE.  FORCED VARIABLES STPR 920
C           ARE ENTERED INTO THE REGRESSION EQUATION BEFORE ALL OTHER   STPR 930
C           INDEPENDENT VARIABLES.  WITHIN THE SET OF FORCED VARIABLES, STPR 940
C           THE ONE TO BE CHOSEN FIRST WILL BE THAT ONE WHICH EXPLAINS  STPR 950
C           THE GREATEST AMOUNT OF VARIANCE.                            STPR 960
C           INSTEAD OF USING, AS A STOPPING CRITERION, A PROPORTION OF  STPR 970
C           THE TOTAL VARIANCE, SOME OTHER CRITERION MAY BE ADDED TO    STPR 980
C           SUBROUTINE STOUT.                                           STPR 990
C                                                                       STPR1000
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  STPR1010
C           STOUT(NSTEP,ANS,L,B,S,T,NSTOP)                              STPR1020
C           THIS SUBROUTINE MUST BE PROVIDED BY THE USER.  IT IS AN     STPR1030
C           OUTPUT ROUTINE WHICH WILL PRINT THE RESULTS OF EACH STEP OF STPR1040
C           THE REGRESSION ANALYSIS.  NSTOP IS AN OPTION CODE WHICH IS  STPR1050
C           ONE IF THE STEPWISE REGRESSION IS TO BE TERMINATED, AND IS  STPR1060
C           ZERO IF IT IS TO CONTINUE.  THE USER MUST CONSIDER THIS IF  STPR1070
C           SOME OTHER STOPPING CRITERION THAN VARIANCE PROPORTION IS TOSTPR1080
C           BE USED.                                                    STPR1090
C                                                                       STPR1100
C        METHOD                                                         STPR1110
C           THE ABBREVIATED DOOLITTLE METHOD IS USED TO (1) DECIDE VARI-STPR1120
C           ABLES ENTERING IN THE REGRESSION AND (2) COMPUTE REGRESSION STPR1130
C           COEFFICIENTS.  REFER TO C. A. BENNETT AND N. L. FRANKLIN,   STPR1140
C           'STATISTICAL ANALYSIS IN CHEMISTRY AND THE CHEMICAL INDUS-  STPR1150
C           TRY', JOHN WILEY AND SONS, 1954, APPENDIX 6A.               STPR1160
C                                                                       STPR1170
C     ..................................................................STPR1180
C                                                                       STPR1190
      SUBROUTINE STPRG (M,N,D,XBAR,IDX,PCT,NSTEP,ANS,L,B,S,T,LL,IER)    STPR1200
C                                                                       STPR1210
      DIMENSION D(1),XBAR(1),IDX(1),NSTEP(1),ANS(1),L(1),B(1),S(1),T(1),STPR1220
     1LL(1)                                                             STPR1230
C                                                                       STPR1240
C     ..................................................................STPR1250
C                                                                       STPR1260
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE  STPR1270
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION      STPR1280
C        STATEMENT WHICH FOLLOWS.                                       STPR1290
C                                                                       STPR1300
C     DOUBLE PRECISION D,XBAR,ANS,B,S,T,RD,RE                           STPR1310
C                                                                       STPR1320
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS    STPR1330
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS      STPR1340
C        ROUTINE.                                                       STPR1350
C                                                                       STPR1360
C        THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO      STPR1370
C        CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  SQRT IN STATEMENTSSTPR1380
C        85,90,114,132,AND 134, MUST BE CHANGED TO DSQRT.               STPR1390
C                                                                       STPR1400
C     ..................................................................STPR1410
C                                                                       STPR1420
C        INITIALIZATION                                                 STPR1430
C                                                                       STPR1440
      IER=0                                                             STPR1450
      ONM=N-1                                                           STPR1460
      NFO=0                                                             STPR1470
      NSTEP(3)=0                                                        STPR1480
      ANS(3)=0.0                                                        STPR1490
      ANS(4)=0.0                                                        STPR1500
      NSTOP=0                                                           STPR1510
C                                                                       STPR1520
C        FIND DEPENDENT VARIABLE, NUMBER OF VARIABLES TO BE FORCED TO   STPR1530
C        ENTER IN THE REGRESSION, AND NUMBER OF VARIABLES TO BE DELETED STPR1540
C                                                                       STPR1550
      DO 30 I=1,M                                                       STPR1560
      LL(I)=1                                                           STPR1570
      IF(IDX(I)) 30, 30, 10                                             STPR1580
   10 IF(IDX(I)-2) 15, 20, 25                                           STPR1590
   15 NFO=NFO+1                                                         STPR1600
      IDX(NFO)=I                                                        STPR1610
      GO TO 30                                                          STPR1620
   20 NSTEP(3)=NSTEP(3)+1                                               STPR1630
      LL(I)=-1                                                          STPR1640
      GO TO 30                                                          STPR1650
   25 MY=I                                                              STPR1660
      NSTEP(1)=MY                                                       STPR1670
      LY=M*(MY-1)                                                       STPR1680
      LYP=LY+MY                                                         STPR1690
      ANS(5)=D(LYP)                                                     STPR1700
   30 CONTINUE                                                          STPR1710
      NSTEP(2)=NFO                                                      STPR1720
C                                                                       STPR1730
C        FIND THE MAXIMUM NUMBER OF STEPS                               STPR1740
C                                                                       STPR1750
      MX=M-NSTEP(3)-1                                                   STPR1760
C                                                                       STPR1770
C        START SELECTION OF VARIABLES                                   STPR1780
C                                                                       STPR1790
      DO 140 NL=1,MX                                                    STPR1800
      RD=0                                                              STPR1810
      IF(NL-NFO) 35, 35, 55                                             STPR1820
C                                                                       STPR1830
C        SELECT NEXT VARIABLE TO ENTER AMONG FORCED VARIABLES           STPR1840
C                                                                       STPR1850
   35 DO 50 I=1,NFO                                                     STPR1860
      K=IDX(I)                                                          STPR1870
      IF(LL(K)) 50, 50, 40                                              STPR1880
   40 LYP=LY+K                                                          STPR1890
      IP=M*(K-1)+K                                                      STPR1900
      RE=D(LYP)*D(LYP)/D(IP)                                            STPR1910
      IF(RD-RE) 45, 50, 50                                              STPR1920
   45 RD=RE                                                             STPR1930
      NEW=K                                                             STPR1940
   50 CONTINUE                                                          STPR1950
      GO TO 75                                                          STPR1960
C                                                                       STPR1970
C        SELECT NEXT VARIABLE TO ENTER AMONG NON-FORCED VARIABLES       STPR1980
C                                                                       STPR1990
   55 DO 70 I=1,M                                                       STPR2000
      IF(I-MY) 60, 70, 60                                               STPR2010
   60 IF(LL(I)) 70, 70, 62                                              STPR2020
   62 LYP=LY+I                                                          STPR2030
      IP=M*(I-1)+I                                                      STPR2040
      RE=D(LYP)*D(LYP)/D(IP)                                            STPR2050
      IF(RD-RE) 64, 70, 70                                              STPR2060
   64 RD=RE                                                             STPR2070
      NEW=I                                                             STPR2080
   70 CONTINUE                                                          STPR2090
C                                                                       STPR2100
C        TEST WHETHER THE PROPORTION OF THE SUM OF SQUARES REDUCED BY   STPR2110
C        THE LAST VARIABLE ENTERED IS GREATER THAN OR EQUAL TO THE      STPR2120
C        SPECIFIED PROPORTION                                           STPR2130
C                                                                       STPR2140
   75 IF(RD) 77,77,76                                                   STPR2150
   76 IF(ANS(5)-(ANS(3)+RD))77,77,78                                    STPR2160
   77 IER=1                                                             STPR2170
      GO TO 150                                                         STPR2180
   78 RE=RD/ANS(5)                                                      STPR2190
      IF(RE-PCT) 150, 80, 80                                            STPR2200
C                                                                       STPR2210
C        IT IS GREATER THAN OR EQUAL                                    STPR2220
C                                                                       STPR2230
   80 LL(NEW)=0                                                         STPR2240
      L(NL)=NEW                                                         STPR2250
      ANS(1)=RD                                                         STPR2260
      ANS(2)=RE                                                         STPR2270
      ANS(3)=ANS(3)+RD                                                  STPR2280
      ANS(4)=ANS(4)+RE                                                  STPR2290
      NSTEP(4)=NL                                                       STPR2300
      NSTEP(5)=NEW                                                      STPR2310
C                                                                       STPR2320
C        COMPUTE MULTIPLE CORRELATION, F-VALUE FOR ANALYSIS OF          STPR2330
C        VARIANCE, AND STANDARD ERROR OF ESTIMATE                       STPR2340
C                                                                       STPR2350
   85 ANS(6)= SQRT(ANS(4))                                              STPR2360
      RD=NL                                                             STPR2370
      RE=ONM-RD                                                         STPR2380
      RE=(ANS(5)-ANS(3))/RE                                             STPR2390
      ANS(7)=(ANS(3)/RD)/RE                                             STPR2400
   90 ANS(8)= SQRT(RE)                                                  STPR2410
C                                                                       STPR2420
C        DIVIDE BY THE PIVOTAL ELEMENT                                  STPR2430
C                                                                       STPR2440
      IP=M*(NEW-1)+NEW                                                  STPR2450
      RD=D(IP)                                                          STPR2460
      LYP=NEW-M                                                         STPR2470
      DO 100 J=1,M                                                      STPR2480
      LYP=LYP+M                                                         STPR2490
      IF(LL(J)) 100, 94, 97                                             STPR2500
   94 IF(J-NEW) 96, 98, 96                                              STPR2510
   96 IJ=M*(J-1)+J                                                      STPR2520
      D(IJ)=D(IJ)+D(LYP)*D(LYP)/RD                                      STPR2530
   97 D(LYP)=D(LYP)/RD                                                  STPR2540
      GO TO 100                                                         STPR2550
   98 D(IP)=1.0/RD                                                      STPR2560
  100 CONTINUE                                                          STPR2570
C                                                                       STPR2580
C        COMPUTE REGRESSION COEFFICIENTS                                STPR2590
C                                                                       STPR2600
      LYP=LY+NEW                                                        STPR2610
      B(NL)=D(LYP)                                                      STPR2620
      IF(NL-1) 112, 112, 105                                            STPR2630
  105 ID=NL-1                                                           STPR2640
      DO 110 J=1,ID                                                     STPR2650
      IJ=NL-J                                                           STPR2660
      KK=L(IJ)                                                          STPR2670
      LYP=LY+KK                                                         STPR2680
      B(IJ)=D(LYP)                                                      STPR2690
      DO 110 K=1,J                                                      STPR2700
      IK=NL-K+1                                                         STPR2710
      MK=L(IK)                                                          STPR2720
      LYP=M*(MK-1)+KK                                                   STPR2730
  110 B(IJ)=B(IJ)-D(LYP)*B(IK)                                          STPR2740
C                                                                       STPR2750
C        COMPUTE INTERCEPT                                              STPR2760
C                                                                       STPR2770
  112 ANS(9)=XBAR(MY)                                                   STPR2780
      DO 115 I=1,NL                                                     STPR2790
      KK=L(I)                                                           STPR2800
      ANS(9)=ANS(9)-B(I)*XBAR(KK)                                       STPR2810
      IJ=M*(KK-1)+KK                                                    STPR2820
  114 S(I)=ANS(8)* SQRT(D(IJ))                                          STPR2830
  115 T(I)=B(I)/S(I)                                                    STPR2840
C                                                                       STPR2850
C        PERFORM A REDUCTION TO ELIMINATE THE LAST VARIABLE ENTERED     STPR2860
C                                                                       STPR2870
      IP=M*(NEW-1)                                                      STPR2880
      DO 130 I=1,M                                                      STPR2890
      IJ=I-M                                                            STPR2900
      IK=NEW-M                                                          STPR2910
      IP=IP+1                                                           STPR2920
      IF(LL(I)) 130, 130, 120                                           STPR2930
  120 DO 126 J=1,M                                                      STPR2940
      IJ=IJ+M                                                           STPR2950
      IK=IK+M                                                           STPR2960
      IF(LL(J)) 126, 122, 122                                           STPR2970
  122 IF(J-NEW) 124, 126, 124                                           STPR2980
  124 D(IJ)=D(IJ)-D(IP)*D(IK)                                           STPR2990
  126 CONTINUE                                                          STPR3000
      D(IP)=D(IP)/(-RD)                                                 STPR3010
  130 CONTINUE                                                          STPR3020
C                                                                       STPR3030
C        ADJUST STANDARD ERROR OF THE ESTIMATE AND MULTIPLE CORRELATION STPR3040
C        COEFFICIENT                                                    STPR3050
C                                                                       STPR3060
      RD=N-NSTEP(4)                                                     STPR3070
      RD=ONM/RD                                                         STPR3080
  132 ANS(10)=SQRT(1.0-(1.0-ANS(6)*ANS(6))*RD)                          STPR3090
  134 ANS(11)=ANS(8)*SQRT(RD)                                           STPR3100
C                                                                       STPR3110
C        CALL THE OUTPUT SUBROUTINE                                     STPR3120
      CALL STOUT (NSTEP,ANS,L,B,S,T,NSTOP)                              STPR3130
C                                                                       STPR3140
C        TEST WHETHER THE STEP-WISE REGRESSION WAS TERMINATED IN        STPR3150
C        SUBROUTINE STOUT                                               STPR3160
C                                                                       STPR3170
      IF(NSTOP) 140, 140, 150                                           STPR3180
C                                                                       STPR3190
  140 CONTINUE                                                          STPR3200
C                                                                       STPR3210
  150 RETURN                                                            STPR3220
      END                                                               STPR3230
