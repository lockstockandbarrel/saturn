C                                                                       DHCG  10
C     ..................................................................DHCG  20
C                                                                       DHCG  30
C        SUBROUTINE DHPCG                                               DHCG  40
C                                                                       DHCG  50
C        PURPOSE                                                        DHCG  60
C           TO SOLVE A SYSTEM OF FIRST ORDER ORDINARY GENERAL           DHCG  70
C           DIFFERENTIAL EQUATIONS WITH GIVEN INITIAL VALUES.           DHCG  80
C                                                                       DHCG  90
C        USAGE                                                          DHCG 100
C           CALL DHPCG (PRMT,Y,DERY,NDIM,IHLF,FCT,OUTP,AUX)             DHCG 110
C           PARAMETERS FCT AND OUTP REQUIRE AN EXTERNAL STATEMENT.      DHCG 120
C                                                                       DHCG 130
C        DESCRIPTION OF PARAMETERS                                      DHCG 140
C           PRMT   - DOUBLE PRECISION INPUT AND OUTPUT VECTOR WITH      DHCG 150
C                    DIMENSION GREATER THAN OR EQUAL TO 5, WHICH        DHCG 160
C                    SPECIFIES THE PARAMETERS OF THE INTERVAL AND OF    DHCG 170
C                    ACCURACY AND WHICH SERVES FOR COMMUNICATION BETWEENDHCG 180
C                    OUTPUT SUBROUTINE (FURNISHED BY THE USER) AND      DHCG 190
C                    SUBROUTINE DHPCG. EXCEPT PRMT(5) THE COMPONENTS    DHCG 200
C                    ARE NOT DESTROYED BY SUBROUTINE DHPCG AND THEY ARE DHCG 210
C           PRMT(1)- LOWER BOUND OF THE INTERVAL (INPUT),               DHCG 220
C           PRMT(2)- UPPER BOUND OF THE INTERVAL (INPUT),               DHCG 230
C           PRMT(3)- INITIAL INCREMENT OF THE INDEPENDENT VARIABLE      DHCG 240
C                    (INPUT),                                           DHCG 250
C           PRMT(4)- UPPER ERROR BOUND (INPUT). IF ABSOLUTE ERROR IS    DHCG 260
C                    GREATER THAN PRMT(4), INCREMENT GETS HALVED.       DHCG 270
C                    IF INCREMENT IS LESS THAN PRMT(3) AND ABSOLUTE     DHCG 280
C                    ERROR LESS THAN PRMT(4)/50, INCREMENT GETS DOUBLED.DHCG 290
C                    THE USER MAY CHANGE PRMT(4) BY MEANS OF HIS        DHCG 300
C                    OUTPUT SUBROUTINE.                                 DHCG 310
C           PRMT(5)- NO INPUT PARAMETER. SUBROUTINE DHPCG INITIALIZES   DHCG 320
C                    PRMT(5)=0. IF THE USER WANTS TO TERMINATE          DHCG 330
C                    SUBROUTINE DHPCG AT ANY OUTPUT POINT, HE HAS TO    DHCG 340
C                    CHANGE PRMT(5) TO NON-ZERO BY MEANS OF SUBROUTINE  DHCG 350
C                    OUTP. FURTHER COMPONENTS OF VECTOR PRMT ARE        DHCG 360
C                    FEASIBLE IF ITS DIMENSION IS DEFINED GREATER       DHCG 370
C                    THAN 5. HOWEVER SUBROUTINE DHPCG DOES NOT REQUIRE  DHCG 380
C                    AND CHANGE THEM. NEVERTHELESS THEY MAY BE USEFUL   DHCG 390
C                    FOR HANDING RESULT VALUES TO THE MAIN PROGRAM      DHCG 400
C                    (CALLING DHPCG) WHICH ARE OBTAINED BY SPECIAL      DHCG 410
C                    MANIPULATIONS WITH OUTPUT DATA IN SUBROUTINE OUTP. DHCG 420
C           Y      - DOUBLE PRECISION INPUT VECTOR OF INITIAL VALUES    DHCG 430
C                    (DESTROYED). LATERON Y IS THE RESULTING VECTOR OF  DHCG 440
C                    DEPENDENT VARIABLES COMPUTED AT INTERMEDIATE       DHCG 450
C                    POINTS X.                                          DHCG 460
C           DERY   - DOUBLE PRECISION INPUT VECTOR OF ERROR WEIGHTS     DHCG 470
C                    (DESTROYED). THE SUM OF ITS COMPONENTS MUST BE     DHCG 480
C                    EQUAL TO 1. LATERON DERY IS THE VECTOR OF          DHCG 490
C                    DERIVATIVES, WHICH BELONG TO FUNCTION VALUES Y AT  DHCG 500
C                    INTERMEDIATE POINTS X.                             DHCG 510
C           NDIM   - AN INPUT VALUE, WHICH SPECIFIES THE NUMBER OF      DHCG 520
C                    EQUATIONS IN THE SYSTEM.                           DHCG 530
C           IHLF   - AN OUTPUT VALUE, WHICH SPECIFIES THE NUMBER OF     DHCG 540
C                    BISECTIONS OF THE INITIAL INCREMENT. IF IHLF GETS  DHCG 550
C                    GREATER THAN 10, SUBROUTINE DHPCG RETURNS WITH     DHCG 560
C                    ERROR MESSAGE IHLF=11 INTO MAIN PROGRAM.           DHCG 570
C                    ERROR MESSAGE IHLF=12 OR IHLF=13 APPEARS IN CASE   DHCG 580
C                    PRMT(3)=0 OR IN CASE SIGN(PRMT(3)).NE.SIGN(PRMT(2)-DHCG 590
C                    PRMT(1)) RESPECTIVELY.                             DHCG 600
C           FCT    - THE NAME OF AN EXTERNAL SUBROUTINE USED. IT        DHCG 610
C                    COMPUTES THE RIGHT HAND SIDES DERY OF THE SYSTEM   DHCG 620
C                    TO GIVEN VALUES OF X AND Y. ITS PARAMETER LIST     DHCG 630
C                    MUST BE X,Y,DERY. THE SUBROUTINE SHOULD NOT        DHCG 640
C                    DESTROY X AND Y.                                   DHCG 650
C           OUTP   - THE NAME OF AN EXTERNAL OUTPUT SUBROUTINE USED.    DHCG 660
C                    ITS PARAMETER LIST MUST BE X,Y,DERY,IHLF,NDIM,PRMT.DHCG 670
C                    NONE OF THESE PARAMETERS (EXCEPT, IF NECESSARY,    DHCG 680
C                    PRMT(4),PRMT(5),...) SHOULD BE CHANGED BY          DHCG 690
C                    SUBROUTINE OUTP. IF PRMT(5) IS CHANGED TO NON-ZERO,DHCG 700
C                    SUBROUTINE DHPCG IS TERMINATED.                    DHCG 710
C           AUX    - DOUBLE PRECISION AUXILIARY STORAGE ARRAY WITH 16   DHCG 720
C                    ROWS AND NDIM COLUMNS.                             DHCG 730
C                                                                       DHCG 740
C        REMARKS                                                        DHCG 750
C           THE PROCEDURE TERMINATES AND RETURNS TO CALLING PROGRAM, IF DHCG 760
C           (1) MORE THAN 10 BISECTIONS OF THE INITIAL INCREMENT ARE    DHCG 770
C               NECESSARY TO GET SATISFACTORY ACCURACY (ERROR MESSAGE   DHCG 780
C               IHLF=11),                                               DHCG 790
C           (2) INITIAL INCREMENT IS EQUAL TO 0 OR HAS WRONG SIGN       DHCG 800
C               (ERROR MESSAGES IHLF=12 OR IHLF=13),                    DHCG 810
C           (3) THE WHOLE INTEGRATION INTERVAL IS WORKED THROUGH,       DHCG 820
C           (4) SUBROUTINE OUTP HAS CHANGED PRMT(5) TO NON-ZERO.        DHCG 830
C                                                                       DHCG 840
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DHCG 850
C           THE EXTERNAL SUBROUTINES FCT(X,Y,DERY) AND                  DHCG 860
C           OUTP(X,Y,DERY,IHLF,NDIM,PRMT) MUST BE FURNISHED BY THE USER.DHCG 870
C                                                                       DHCG 880
C        METHOD                                                         DHCG 890
C           EVALUATION IS DONE BY MEANS OF HAMMINGS MODIFIED PREDICTOR- DHCG 900
C           CORRECTOR METHOD. IT IS A FOURTH ORDER METHOD, USING 4      DHCG 910
C           PRECEEDING POINTS FOR COMPUTATION OF A NEW VECTOR Y OF THE  DHCG 920
C           DEPENDENT VARIABLES.                                        DHCG 930
C           FOURTH ORDER RUNGE-KUTTA METHOD SUGGESTED BY RALSTON IS     DHCG 940
C           USED FOR ADJUSTMENT OF THE INITIAL INCREMENT AND FOR        DHCG 950
C           COMPUTATION OF STARTING VALUES.                             DHCG 960
C           SUBROUTINE DHPCG AUTOMATICALLY ADJUSTS THE INCREMENT DURING DHCG 970
C           THE WHOLE COMPUTATION BY HALVING OR DOUBLING.               DHCG 980
C           TO GET FULL FLEXIBILITY IN OUTPUT, AN OUTPUT SUBROUTINE     DHCG 990
C           MUST BE CODED BY THE USER.                                  DHCG1000
C           FOR REFERENCE, SEE                                          DHCG1010
C           (1)  RALSTON/WILF, MATHEMATICAL METHODS FOR DIGITAL         DHCG1020
C                COMPUTERS, WILEY, NEW YORK/LONDON, 1960, PP.95-109.    DHCG1030
C           (2)  RALSTON, RUNGE-KUTTA METHODS WITH MINIMUM ERROR BOUNDS,DHCG1040
C                MTAC, VOL.16, ISS.80 (1962), PP.431-437.               DHCG1050
C                                                                       DHCG1060
C     ..................................................................DHCG1070
C                                                                       DHCG1080
      SUBROUTINE DHPCG(PRMT,Y,DERY,NDIM,IHLF,FCT,OUTP,AUX)              DHCG1090
C                                                                       DHCG1100
C                                                                       DHCG1110
      DIMENSION PRMT(1),Y(1),DERY(1),AUX(16,1)                          DHCG1120
      DOUBLE PRECISION Y,DERY,AUX,PRMT,X,H,Z,DELT                       DHCG1130
      N=1                                                               DHCG1140
      IHLF=0                                                            DHCG1150
      X=PRMT(1)                                                         DHCG1160
      H=PRMT(3)                                                         DHCG1170
      PRMT(5)=0.D0                                                      DHCG1180
      DO 1 I=1,NDIM                                                     DHCG1190
      AUX(16,I)=0.D0                                                    DHCG1200
      AUX(15,I)=DERY(I)                                                 DHCG1210
    1 AUX(1,I)=Y(I)                                                     DHCG1220
      IF(H*(PRMT(2)-X))3,2,4                                            DHCG1230
C                                                                       DHCG1240
C     ERROR RETURNS                                                     DHCG1250
    2 IHLF=12                                                           DHCG1260
      GOTO 4                                                            DHCG1270
    3 IHLF=13                                                           DHCG1280
C                                                                       DHCG1290
C     COMPUTATION OF DERY FOR STARTING VALUES                           DHCG1300
    4 CALL FCT(X,Y,DERY)                                                DHCG1310
C                                                                       DHCG1320
C     RECORDING OF STARTING VALUES                                      DHCG1330
      CALL OUTP(X,Y,DERY,IHLF,NDIM,PRMT)                                DHCG1340
      IF(PRMT(5))6,5,6                                                  DHCG1350
    5 IF(IHLF)7,7,6                                                     DHCG1360
    6 RETURN                                                            DHCG1370
    7 DO 8 I=1,NDIM                                                     DHCG1380
    8 AUX(8,I)=DERY(I)                                                  DHCG1390
C                                                                       DHCG1400
C     COMPUTATION OF AUX(2,I)                                           DHCG1410
      ISW=1                                                             DHCG1420
      GOTO 100                                                          DHCG1430
C                                                                       DHCG1440
    9 X=X+H                                                             DHCG1450
      DO 10 I=1,NDIM                                                    DHCG1460
   10 AUX(2,I)=Y(I)                                                     DHCG1470
C                                                                       DHCG1480
C     INCREMENT H IS TESTED BY MEANS OF BISECTION                       DHCG1490
   11 IHLF=IHLF+1                                                       DHCG1500
      X=X-H                                                             DHCG1510
      DO 12 I=1,NDIM                                                    DHCG1520
   12 AUX(4,I)=AUX(2,I)                                                 DHCG1530
      H=.5D0*H                                                          DHCG1540
      N=1                                                               DHCG1550
      ISW=2                                                             DHCG1560
      GOTO 100                                                          DHCG1570
C                                                                       DHCG1580
   13 X=X+H                                                             DHCG1590
      CALL FCT(X,Y,DERY)                                                DHCG1600
      N=2                                                               DHCG1610
      DO 14 I=1,NDIM                                                    DHCG1620
      AUX(2,I)=Y(I)                                                     DHCG1630
   14 AUX(9,I)=DERY(I)                                                  DHCG1640
      ISW=3                                                             DHCG1650
      GOTO 100                                                          DHCG1660
C                                                                       DHCG1670
C     COMPUTATION OF TEST VALUE DELT                                    DHCG1680
   15 DELT=0.D0                                                         DHCG1690
      DO 16 I=1,NDIM                                                    DHCG1700
   16 DELT=DELT+AUX(15,I)*DABS(Y(I)-AUX(4,I))                           DHCG1710
      DELT=.066666666666666667D0*DELT                                   DHCG1720
      IF(DELT-PRMT(4))19,19,17                                          DHCG1730
   17 IF(IHLF-10)11,18,18                                               DHCG1740
C                                                                       DHCG1750
C     NO SATISFACTORY ACCURACY AFTER 10 BISECTIONS. ERROR MESSAGE.      DHCG1760
   18 IHLF=11                                                           DHCG1770
      X=X+H                                                             DHCG1780
      GOTO 4                                                            DHCG1790
C                                                                       DHCG1800
C     THERE IS SATISFACTORY ACCURACY AFTER LESS THAN 11 BISECTIONS.     DHCG1810
   19 X=X+H                                                             DHCG1820
      CALL FCT(X,Y,DERY)                                                DHCG1830
      DO 20 I=1,NDIM                                                    DHCG1840
      AUX(3,I)=Y(I)                                                     DHCG1850
   20 AUX(10,I)=DERY(I)                                                 DHCG1860
      N=3                                                               DHCG1870
      ISW=4                                                             DHCG1880
      GOTO 100                                                          DHCG1890
C                                                                       DHCG1900
   21 N=1                                                               DHCG1910
      X=X+H                                                             DHCG1920
      CALL FCT(X,Y,DERY)                                                DHCG1930
      X=PRMT(1)                                                         DHCG1940
      DO 22 I=1,NDIM                                                    DHCG1950
      AUX(11,I)=DERY(I)                                                 DHCG1960
   220Y(I)=AUX(1,I)+H*(.375D0*AUX(8,I)+.7916666666666667D0*AUX(9,I)     DHCG1970
     1-.20833333333333333D0*AUX(10,I)+.041666666666666667D0*DERY(I))    DHCG1980
   23 X=X+H                                                             DHCG1990
      N=N+1                                                             DHCG2000
      CALL FCT(X,Y,DERY)                                                DHCG2010
      CALL OUTP(X,Y,DERY,IHLF,NDIM,PRMT)                                DHCG2020
      IF(PRMT(5))6,24,6                                                 DHCG2030
   24 IF(N-4)25,200,200                                                 DHCG2040
   25 DO 26 I=1,NDIM                                                    DHCG2050
      AUX(N,I)=Y(I)                                                     DHCG2060
   26 AUX(N+7,I)=DERY(I)                                                DHCG2070
      IF(N-3)27,29,200                                                  DHCG2080
C                                                                       DHCG2090
   27 DO 28 I=1,NDIM                                                    DHCG2100
      DELT=AUX(9,I)+AUX(9,I)                                            DHCG2110
      DELT=DELT+DELT                                                    DHCG2120
   28 Y(I)=AUX(1,I)+.33333333333333333D0*H*(AUX(8,I)+DELT+AUX(10,I))    DHCG2130
      GOTO 23                                                           DHCG2140
C                                                                       DHCG2150
   29 DO 30 I=1,NDIM                                                    DHCG2160
      DELT=AUX(9,I)+AUX(10,I)                                           DHCG2170
      DELT=DELT+DELT+DELT                                               DHCG2180
   30 Y(I)=AUX(1,I)+.375D0*H*(AUX(8,I)+DELT+AUX(11,I))                  DHCG2190
      GOTO 23                                                           DHCG2200
C                                                                       DHCG2210
C     THE FOLLOWING PART OF SUBROUTINE DHPCG COMPUTES BY MEANS OF       DHCG2220
C     RUNGE-KUTTA METHOD STARTING VALUES FOR THE NOT SELF-STARTING      DHCG2230
C     PREDICTOR-CORRECTOR METHOD.                                       DHCG2240
  100 DO 101 I=1,NDIM                                                   DHCG2250
      Z=H*AUX(N+7,I)                                                    DHCG2260
      AUX(5,I)=Z                                                        DHCG2270
  101 Y(I)=AUX(N,I)+.4D0*Z                                              DHCG2280
C     Z IS AN AUXILIARY STORAGE LOCATION                                DHCG2290
C                                                                       DHCG2300
      Z=X+.4D0*H                                                        DHCG2310
      CALL FCT(Z,Y,DERY)                                                DHCG2320
      DO 102 I=1,NDIM                                                   DHCG2330
      Z=H*DERY(I)                                                       DHCG2340
      AUX(6,I)=Z                                                        DHCG2350
  102 Y(I)=AUX(N,I)+.29697760924775360D0*AUX(5,I)+.15875964497103583D0*ZDHCG2360
C                                                                       DHCG2370
      Z=X+.45573725421878943D0*H                                        DHCG2380
      CALL FCT(Z,Y,DERY)                                                DHCG2390
      DO 103 I=1,NDIM                                                   DHCG2400
      Z=H*DERY(I)                                                       DHCG2410
      AUX(7,I)=Z                                                        DHCG2420
  103 Y(I)=AUX(N,I)+.21810038822592047D0*AUX(5,I)-3.0509651486929308D0* DHCG2430
     1AUX(6,I)+3.8328647604670103D0*Z                                   DHCG2440
C                                                                       DHCG2450
      Z=X+H                                                             DHCG2460
      CALL FCT(Z,Y,DERY)                                                DHCG2470
      DO 104 I=1,NDIM                                                   DHCG2480
  1040Y(I)=AUX(N,I)+.17476028226269037D0*AUX(5,I)-.55148066287873294D0* DHCG2490
     1AUX(6,I)+1.2055355993965235D0*AUX(7,I)+.17118478121951903D0*      DHCG2500
     2H*DERY(I)                                                         DHCG2510
      GOTO(9,13,15,21),ISW                                              DHCG2520
C                                                                       DHCG2530
C     POSSIBLE BREAK-POINT FOR LINKAGE                                  DHCG2540
C                                                                       DHCG2550
C     STARTING VALUES ARE COMPUTED.                                     DHCG2560
C     NOW START HAMMINGS MODIFIED PREDICTOR-CORRECTOR METHOD.           DHCG2570
  200 ISTEP=3                                                           DHCG2580
  201 IF(N-8)204,202,204                                                DHCG2590
C                                                                       DHCG2600
C     N=8 CAUSES THE ROWS OF AUX TO CHANGE THEIR STORAGE LOCATIONS      DHCG2610
  202 DO 203 N=2,7                                                      DHCG2620
      DO 203 I=1,NDIM                                                   DHCG2630
      AUX(N-1,I)=AUX(N,I)                                               DHCG2640
  203 AUX(N+6,I)=AUX(N+7,I)                                             DHCG2650
      N=7                                                               DHCG2660
C                                                                       DHCG2670
C     N LESS THAN 8 CAUSES N+1 TO GET N                                 DHCG2680
  204 N=N+1                                                             DHCG2690
C                                                                       DHCG2700
C     COMPUTATION OF NEXT VECTOR Y                                      DHCG2710
      DO 205 I=1,NDIM                                                   DHCG2720
      AUX(N-1,I)=Y(I)                                                   DHCG2730
  205 AUX(N+6,I)=DERY(I)                                                DHCG2740
      X=X+H                                                             DHCG2750
  206 ISTEP=ISTEP+1                                                     DHCG2760
      DO 207 I=1,NDIM                                                   DHCG2770
     0DELT=AUX(N-4,I)+1.3333333333333333D0*H*(AUX(N+6,I)+AUX(N+6,I)-    DHCG2780
     1AUX(N+5,I)+AUX(N+4,I)+AUX(N+4,I))                                 DHCG2790
      Y(I)=DELT-.9256198347107438D0*AUX(16,I)                           DHCG2800
  207 AUX(16,I)=DELT                                                    DHCG2810
C     PREDICTOR IS NOW GENERATED IN ROW 16 OF AUX, MODIFIED PREDICTOR   DHCG2820
C     IS GENERATED IN Y. DELT MEANS AN AUXILIARY STORAGE.               DHCG2830
C                                                                       DHCG2840
      CALL FCT(X,Y,DERY)                                                DHCG2850
C     DERIVATIVE OF MODIFIED PREDICTOR IS GENERATED IN DERY             DHCG2860
C                                                                       DHCG2870
      DO 208 I=1,NDIM                                                   DHCG2880
     0DELT=.125D0*(9.D0*AUX(N-1,I)-AUX(N-3,I)+3.D0*H*(DERY(I)+AUX(N+6,I)DHCG2890
     1+AUX(N+6,I)-AUX(N+5,I)))                                          DHCG2900
      AUX(16,I)=AUX(16,I)-DELT                                          DHCG2910
  208 Y(I)=DELT+.07438016528925620D0*AUX(16,I)                          DHCG2920
C                                                                       DHCG2930
C     TEST WHETHER H MUST BE HALVED OR DOUBLED                          DHCG2940
      DELT=0.D0                                                         DHCG2950
      DO 209 I=1,NDIM                                                   DHCG2960
  209 DELT=DELT+AUX(15,I)*DABS(AUX(16,I))                               DHCG2970
      IF(DELT-PRMT(4))210,222,222                                       DHCG2980
C                                                                       DHCG2990
C     H MUST NOT BE HALVED. THAT MEANS Y(I) ARE GOOD.                   DHCG3000
  210 CALL FCT(X,Y,DERY)                                                DHCG3010
      CALL OUTP(X,Y,DERY,IHLF,NDIM,PRMT)                                DHCG3020
      IF(PRMT(5))212,211,212                                            DHCG3030
  211 IF(IHLF-11)213,212,212                                            DHCG3040
  212 RETURN                                                            DHCG3050
  213 IF(H*(X-PRMT(2)))214,212,212                                      DHCG3060
  214 IF(DABS(X-PRMT(2))-.1D0*DABS(H))212,215,215                       DHCG3070
  215 IF(DELT-.02D0*PRMT(4))216,216,201                                 DHCG3080
C                                                                       DHCG3090
C                                                                       DHCG3100
C     H COULD BE DOUBLED IF ALL NECESSARY PRECEEDING VALUES ARE         DHCG3110
C     AVAILABLE                                                         DHCG3120
  216 IF(IHLF)201,201,217                                               DHCG3130
  217 IF(N-7)201,218,218                                                DHCG3140
  218 IF(ISTEP-4)201,219,219                                            DHCG3150
  219 IMOD=ISTEP/2                                                      DHCG3160
      IF(ISTEP-IMOD-IMOD)201,220,201                                    DHCG3170
  220 H=H+H                                                             DHCG3180
      IHLF=IHLF-1                                                       DHCG3190
      ISTEP=0                                                           DHCG3200
      DO 221 I=1,NDIM                                                   DHCG3210
      AUX(N-1,I)=AUX(N-2,I)                                             DHCG3220
      AUX(N-2,I)=AUX(N-4,I)                                             DHCG3230
      AUX(N-3,I)=AUX(N-6,I)                                             DHCG3240
      AUX(N+6,I)=AUX(N+5,I)                                             DHCG3250
      AUX(N+5,I)=AUX(N+3,I)                                             DHCG3260
      AUX(N+4,I)=AUX(N+1,I)                                             DHCG3270
      DELT=AUX(N+6,I)+AUX(N+5,I)                                        DHCG3280
      DELT=DELT+DELT+DELT                                               DHCG3290
  2210AUX(16,I)=8.962962962962963D0*(Y(I)-AUX(N-3,I))                   DHCG3300
     1-3.3611111111111111D0*H*(DERY(I)+DELT+AUX(N+4,I))                 DHCG3310
      GOTO 201                                                          DHCG3320
C                                                                       DHCG3330
C                                                                       DHCG3340
C     H MUST BE HALVED                                                  DHCG3350
  222 IHLF=IHLF+1                                                       DHCG3360
      IF(IHLF-10)223,223,210                                            DHCG3370
  223 H=.5D0*H                                                          DHCG3380
      ISTEP=0                                                           DHCG3390
      DO 224 I=1,NDIM                                                   DHCG3400
     0Y(I)=.390625D-2*(8.D1*AUX(N-1,I)+135.D0*AUX(N-2,I)+4.D1*AUX(N-3,I)DHCG3410
     1+AUX(N-4,I))-.1171875D0*(AUX(N+6,I)-6.D0*AUX(N+5,I)-AUX(N+4,I))*H DHCG3420
     0AUX(N-4,I)=.390625D-2*(12.D0*AUX(N-1,I)+135.D0*AUX(N-2,I)+        DHCG3430
     1108.D0*AUX(N-3,I)+AUX(N-4,I))-.0234375D0*(AUX(N+6,I)+             DHCG3440
     218.D0*AUX(N+5,I)-9.D0*AUX(N+4,I))*H                               DHCG3450
      AUX(N-3,I)=AUX(N-2,I)                                             DHCG3460
  224 AUX(N+4,I)=AUX(N+5,I)                                             DHCG3470
      X=X-H                                                             DHCG3480
      DELT=X-(H+H)                                                      DHCG3490
      CALL FCT(DELT,Y,DERY)                                             DHCG3500
      DO 225 I=1,NDIM                                                   DHCG3510
      AUX(N-2,I)=Y(I)                                                   DHCG3520
      AUX(N+5,I)=DERY(I)                                                DHCG3530
  225 Y(I)=AUX(N-4,I)                                                   DHCG3540
      DELT=DELT-(H+H)                                                   DHCG3550
      CALL FCT(DELT,Y,DERY)                                             DHCG3560
      DO 226 I=1,NDIM                                                   DHCG3570
      DELT=AUX(N+5,I)+AUX(N+4,I)                                        DHCG3580
      DELT=DELT+DELT+DELT                                               DHCG3590
     0AUX(16,I)=8.962962962962963D0*(AUX(N-1,I)-Y(I))                   DHCG3600
     1-3.3611111111111111D0*H*(AUX(N+6,I)+DELT+DERY(I))                 DHCG3610
  226 AUX(N+3,I)=DERY(I)                                                DHCG3620
      GOTO 206                                                          DHCG3630
      END                                                               DHCG3640
