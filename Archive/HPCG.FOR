C                                                                       HPCG  10
C                                                                       HPCG  20
C     ..................................................................HPCG  30
C                                                                       HPCG  40
C        SUBROUTINE HPCG                                                HPCG  50
C                                                                       HPCG  60
C        PURPOSE                                                        HPCG  70
C           TO SOLVE A SYSTEM OF FIRST ORDER ORDINARY GENERAL           HPCG  80
C           DIFFERENTIAL EQUATIONS WITH GIVEN INITIAL VALUES.           HPCG  90
C                                                                       HPCG 100
C        USAGE                                                          HPCG 110
C           CALL HPCG (PRMT,Y,DERY,NDIM,IHLF,FCT,OUTP,AUX)              HPCG 120
C           PARAMETERS FCT AND OUTP REQUIRE AN EXTERNAL STATEMENT.      HPCG 130
C                                                                       HPCG 140
C        DESCRIPTION OF PARAMETERS                                      HPCG 150
C           PRMT   - AN INPUT AND OUTPUT VECTOR WITH DIMENSION GREATER  HPCG 160
C                    OR EQUAL TO 5, WHICH SPECIFIES THE PARAMETERS OF   HPCG 170
C                    THE INTERVAL AND OF ACCURACY AND WHICH SERVES FOR  HPCG 180
C                    COMMUNICATION BETWEEN OUTPUT SUBROUTINE (FURNISHED HPCG 190
C                    BY THE USER) AND SUBROUTINE HPCG. EXCEPT PRMT(5)   HPCG 200
C                    THE COMPONENTS ARE NOT DESTROYED BY SUBROUTINE     HPCG 210
C                    HPCG AND THEY ARE                                  HPCG 220
C           PRMT(1)- LOWER BOUND OF THE INTERVAL (INPUT),               HPCG 230
C           PRMT(2)- UPPER BOUND OF THE INTERVAL (INPUT),               HPCG 240
C           PRMT(3)- INITIAL INCREMENT OF THE INDEPENDENT VARIABLE      HPCG 250
C                    (INPUT),                                           HPCG 260
C           PRMT(4)- UPPER ERROR BOUND (INPUT). IF ABSOLUTE ERROR IS    HPCG 270
C                    GREATER THAN PRMT(4), INCREMENT GETS HALVED.       HPCG 280
C                    IF INCREMENT IS LESS THAN PRMT(3) AND ABSOLUTE     HPCG 290
C                    ERROR LESS THAN PRMT(4)/50, INCREMENT GETS DOUBLED.HPCG 300
C                    THE USER MAY CHANGE PRMT(4) BY MEANS OF HIS        HPCG 310
C                    OUTPUT SUBROUTINE.                                 HPCG 320
C           PRMT(5)- NO INPUT PARAMETER. SUBROUTINE HPCG INITIALIZES    HPCG 330
C                    PRMT(5)=0. IF THE USER WANTS TO TERMINATE          HPCG 340
C                    SUBROUTINE HPCG AT ANY OUTPUT POINT, HE HAS TO     HPCG 350
C                    CHANGE PRMT(5) TO NON-ZERO BY MEANS OF SUBROUTINE  HPCG 360
C                    OUTP. FURTHER COMPONENTS OF VECTOR PRMT ARE        HPCG 370
C                    FEASIBLE IF ITS DIMENSION IS DEFINED GREATER       HPCG 380
C                    THAN 5. HOWEVER SUBROUTINE HPCG DOES NOT REQUIRE   HPCG 390
C                    AND CHANGE THEM. NEVERTHELESS THEY MAY BE USEFUL   HPCG 400
C                    FOR HANDING RESULT VALUES TO THE MAIN PROGRAM      HPCG 410
C                    (CALLING HPCG) WHICH ARE OBTAINED BY SPECIAL       HPCG 420
C                    MANIPULATIONS WITH OUTPUT DATA IN SUBROUTINE OUTP. HPCG 430
C           Y      - INPUT VECTOR OF INITIAL VALUES.  (DESTROYED)       HPCG 440
C                    LATERON Y IS THE RESULTING VECTOR OF DEPENDENT     HPCG 450
C                    VARIABLES COMPUTED AT INTERMEDIATE POINTS X.       HPCG 460
C           DERY   - INPUT VECTOR OF ERROR WEIGHTS.  (DESTROYED)        HPCG 470
C                    THE SUM OF ITS COMPONENTS MUST BE EQUAL TO 1.      HPCG 480
C                    LATERON DERY IS THE VECTOR OF DERIVATIVES, WHICH   HPCG 490
C                    BELONG TO FUNCTION VALUES Y AT A POINT X.          HPCG 500
C           NDIM   - AN INPUT VALUE, WHICH SPECIFIES THE NUMBER OF      HPCG 510
C                    EQUATIONS IN THE SYSTEM.                           HPCG 520
C           IHLF   - AN OUTPUT VALUE, WHICH SPECIFIES THE NUMBER OF     HPCG 530
C                    BISECTIONS OF THE INITIAL INCREMENT. IF IHLF GETS  HPCG 540
C                    GREATER THAN 10, SUBROUTINE HPCG RETURNS WITH      HPCG 550
C                    ERROR MESSAGE IHLF=11 INTO MAIN PROGRAM.           HPCG 560
C                    ERROR MESSAGE IHLF=12 OR IHLF=13 APPEARS IN CASE   HPCG 570
C                    PRMT(3)=0 OR IN CASE SIGN(PRMT(3)).NE.SIGN(PRMT(2)-HPCG 580
C                    PRMT(1)) RESPECTIVELY.                             HPCG 590
C           FCT    - THE NAME OF AN EXTERNAL SUBROUTINE USED. IT        HPCG 600
C                    COMPUTES THE RIGHT HAND SIDES DERY OF THE SYSTEM   HPCG 610
C                    TO GIVEN VALUES OF X AND Y. ITS PARAMETER LIST     HPCG 620
C                    MUST BE X,Y,DERY. THE SUBROUTINE SHOULD NOT        HPCG 630
C                    DESTROY X AND Y.                                   HPCG 640
C           OUTP   - THE NAME OF AN EXTERNAL OUTPUT SUBROUTINE USED.    HPCG 650
C                    ITS PARAMETER LIST MUST BE X,Y,DERY,IHLF,NDIM,PRMT.HPCG 660
C                    NONE OF THESE PARAMETERS (EXCEPT, IF NECESSARY,    HPCG 670
C                    PRMT(4),PRMT(5),...) SHOULD BE CHANGED BY          HPCG 680
C                    SUBROUTINE OUTP. IF PRMT(5) IS CHANGED TO NON-ZERO,HPCG 690
C                    SUBROUTINE HPCG IS TERMINATED.                     HPCG 700
C           AUX    - AN AUXILIARY STORAGE ARRAY WITH 16 ROWS AND NDIM   HPCG 710
C                    COLUMNS.                                           HPCG 720
C                                                                       HPCG 730
C        REMARKS                                                        HPCG 740
C           THE PROCEDURE TERMINATES AND RETURNS TO CALLING PROGRAM, IF HPCG 750
C           (1) MORE THAN 10 BISECTIONS OF THE INITIAL INCREMENT ARE    HPCG 760
C               NECESSARY TO GET SATISFACTORY ACCURACY (ERROR MESSAGE   HPCG 770
C               IHLF=11),                                               HPCG 780
C           (2) INITIAL INCREMENT IS EQUAL TO 0 OR HAS WRONG SIGN       HPCG 790
C               (ERROR MESSAGES IHLF=12 OR IHLF=13),                    HPCG 800
C           (3) THE WHOLE INTEGRATION INTERVAL IS WORKED THROUGH,       HPCG 810
C           (4) SUBROUTINE OUTP HAS CHANGED PRMT(5) TO NON-ZERO.        HPCG 820
C                                                                       HPCG 830
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  HPCG 840
C           THE EXTERNAL SUBROUTINES FCT(X,Y,DERY) AND                  HPCG 850
C           OUTP(X,Y,DERY,IHLF,NDIM,PRMT) MUST BE FURNISHED BY THE USER.HPCG 860
C                                                                       HPCG 870
C        METHOD                                                         HPCG 880
C           EVALUATION IS DONE BY MEANS OF HAMMINGS MODIFIED PREDICTOR- HPCG 890
C           CORRECTOR METHOD. IT IS A FOURTH ORDER METHOD, USING 4      HPCG 900
C           PRECEEDING POINTS FOR COMPUTATION OF A NEW VECTOR Y OF THE  HPCG 910
C           DEPENDENT VARIABLES.                                        HPCG 920
C           FOURTH ORDER RUNGE-KUTTA METHOD SUGGESTED BY RALSTON IS     HPCG 930
C           USED FOR ADJUSTMENT OF THE INITIAL INCREMENT AND FOR        HPCG 940
C           COMPUTATION OF STARTING VALUES.                             HPCG 950
C           SUBROUTINE HPCG AUTOMATICALLY ADJUSTS THE INCREMENT DURING  HPCG 960
C           THE WHOLE COMPUTATION BY HALVING OR DOUBLING.               HPCG 970
C           TO GET FULL FLEXIBILITY IN OUTPUT, AN OUTPUT SUBROUTINE     HPCG 980
C           MUST BE CODED BY THE USER.                                  HPCG 990
C           FOR REFERENCE, SEE                                          HPCG1000
C           (1)  RALSTON/WILF, MATHEMATICAL METHODS FOR DIGITAL         HPCG1010
C                COMPUTERS, WILEY, NEW YORK/LONDON, 1960, PP.95-109.    HPCG1020
C           (2)  RALSTON, RUNGE-KUTTA METHODS WITH MINIMUM ERROR BOUNDS,HPCG1030
C                MTAC, VOL.16, ISS.80 (1962), PP.431-437.               HPCG1040
C                                                                       HPCG1050
C     ..................................................................HPCG1060
C                                                                       HPCG1070
      SUBROUTINE HPCG(PRMT,Y,DERY,NDIM,IHLF,FCT,OUTP,AUX)               HPCG1080
C                                                                       HPCG1090
C                                                                       HPCG1100
      DIMENSION PRMT(1),Y(1),DERY(1),AUX(16,1)                          HPCG1110
      N=1                                                               HPCG1120
      IHLF=0                                                            HPCG1130
      X=PRMT(1)                                                         HPCG1140
      H=PRMT(3)                                                         HPCG1150
      PRMT(5)=0.                                                        HPCG1160
      DO 1 I=1,NDIM                                                     HPCG1170
      AUX(16,I)=0.                                                      HPCG1180
      AUX(15,I)=DERY(I)                                                 HPCG1190
    1 AUX(1,I)=Y(I)                                                     HPCG1200
      IF(H*(PRMT(2)-X))3,2,4                                            HPCG1210
C                                                                       HPCG1220
C     ERROR RETURNS                                                     HPCG1230
    2 IHLF=12                                                           HPCG1240
      GOTO 4                                                            HPCG1250
    3 IHLF=13                                                           HPCG1260
C                                                                       HPCG1270
C     COMPUTATION OF DERY FOR STARTING VALUES                           HPCG1280
    4 CALL FCT(X,Y,DERY)                                                HPCG1290
C                                                                       HPCG1300
C     RECORDING OF STARTING VALUES                                      HPCG1310
      CALL OUTP(X,Y,DERY,IHLF,NDIM,PRMT)                                HPCG1320
      IF(PRMT(5))6,5,6                                                  HPCG1330
    5 IF(IHLF)7,7,6                                                     HPCG1340
    6 RETURN                                                            HPCG1350
    7 DO 8 I=1,NDIM                                                     HPCG1360
    8 AUX(8,I)=DERY(I)                                                  HPCG1370
C                                                                       HPCG1380
C     COMPUTATION OF AUX(2,I)                                           HPCG1390
      ISW=1                                                             HPCG1400
      GOTO 100                                                          HPCG1410
C                                                                       HPCG1420
    9 X=X+H                                                             HPCG1430
      DO 10 I=1,NDIM                                                    HPCG1440
   10 AUX(2,I)=Y(I)                                                     HPCG1450
C                                                                       HPCG1460
C     INCREMENT H IS TESTED BY MEANS OF BISECTION                       HPCG1470
   11 IHLF=IHLF+1                                                       HPCG1480
      X=X-H                                                             HPCG1490
      DO 12 I=1,NDIM                                                    HPCG1500
   12 AUX(4,I)=AUX(2,I)                                                 HPCG1510
      H=.5*H                                                            HPCG1520
      N=1                                                               HPCG1530
      ISW=2                                                             HPCG1540
      GOTO 100                                                          HPCG1550
C                                                                       HPCG1560
   13 X=X+H                                                             HPCG1570
      CALL FCT(X,Y,DERY)                                                HPCG1580
      N=2                                                               HPCG1590
      DO 14 I=1,NDIM                                                    HPCG1600
      AUX(2,I)=Y(I)                                                     HPCG1610
   14 AUX(9,I)=DERY(I)                                                  HPCG1620
      ISW=3                                                             HPCG1630
      GOTO 100                                                          HPCG1640
C                                                                       HPCG1650
C     COMPUTATION OF TEST VALUE DELT                                    HPCG1660
   15 DELT=0.                                                           HPCG1670
      DO 16 I=1,NDIM                                                    HPCG1680
   16 DELT=DELT+AUX(15,I)*ABS(Y(I)-AUX(4,I))                            HPCG1690
      DELT=.06666667*DELT                                               HPCG1700
      IF(DELT-PRMT(4))19,19,17                                          HPCG1710
   17 IF(IHLF-10)11,18,18                                               HPCG1720
C                                                                       HPCG1730
C     NO SATISFACTORY ACCURACY AFTER 10 BISECTIONS. ERROR MESSAGE.      HPCG1740
   18 IHLF=11                                                           HPCG1750
      X=X+H                                                             HPCG1760
      GOTO 4                                                            HPCG1770
C                                                                       HPCG1780
C     THERE IS SATISFACTORY ACCURACY AFTER LESS THAN 11 BISECTIONS.     HPCG1790
   19 X=X+H                                                             HPCG1800
      CALL FCT(X,Y,DERY)                                                HPCG1810
      DO 20 I=1,NDIM                                                    HPCG1820
      AUX(3,I)=Y(I)                                                     HPCG1830
   20 AUX(10,I)=DERY(I)                                                 HPCG1840
      N=3                                                               HPCG1850
      ISW=4                                                             HPCG1860
      GOTO 100                                                          HPCG1870
C                                                                       HPCG1880
   21 N=1                                                               HPCG1890
      X=X+H                                                             HPCG1900
      CALL FCT(X,Y,DERY)                                                HPCG1910
      X=PRMT(1)                                                         HPCG1920
      DO 22 I=1,NDIM                                                    HPCG1930
      AUX(11,I)=DERY(I)                                                 HPCG1940
   220Y(I)=AUX(1,I)+H*(.375*AUX(8,I)+.7916667*AUX(9,I)                  HPCG1950
     1-.2083333*AUX(10,I)+.04166667*DERY(I))                            HPCG1960
   23 X=X+H                                                             HPCG1970
      N=N+1                                                             HPCG1980
      CALL FCT(X,Y,DERY)                                                HPCG1990
      CALL OUTP(X,Y,DERY,IHLF,NDIM,PRMT)                                HPCG2000
      IF(PRMT(5))6,24,6                                                 HPCG2010
   24 IF(N-4)25,200,200                                                 HPCG2020
   25 DO 26 I=1,NDIM                                                    HPCG2030
      AUX(N,I)=Y(I)                                                     HPCG2040
   26 AUX(N+7,I)=DERY(I)                                                HPCG2050
      IF(N-3)27,29,200                                                  HPCG2060
C                                                                       HPCG2070
   27 DO 28 I=1,NDIM                                                    HPCG2080
      DELT=AUX(9,I)+AUX(9,I)                                            HPCG2090
      DELT=DELT+DELT                                                    HPCG2100
   28 Y(I)=AUX(1,I)+.3333333*H*(AUX(8,I)+DELT+AUX(10,I))                HPCG2110
      GOTO 23                                                           HPCG2120
C                                                                       HPCG2130
   29 DO 30 I=1,NDIM                                                    HPCG2140
      DELT=AUX(9,I)+AUX(10,I)                                           HPCG2150
      DELT=DELT+DELT+DELT                                               HPCG2160
   30 Y(I)=AUX(1,I)+.375*H*(AUX(8,I)+DELT+AUX(11,I))                    HPCG2170
      GOTO 23                                                           HPCG2180
C                                                                       HPCG2190
C     THE FOLLOWING PART OF SUBROUTINE HPCG COMPUTES BY MEANS OF        HPCG2200
C     RUNGE-KUTTA METHOD STARTING VALUES FOR THE NOT SELF-STARTING      HPCG2210
C     PREDICTOR-CORRECTOR METHOD.                                       HPCG2220
  100 DO 101 I=1,NDIM                                                   HPCG2230
      Z=H*AUX(N+7,I)                                                    HPCG2240
      AUX(5,I)=Z                                                        HPCG2250
  101 Y(I)=AUX(N,I)+.4*Z                                                HPCG2260
C     Z IS AN AUXILIARY STORAGE LOCATION                                HPCG2270
C                                                                       HPCG2280
      Z=X+.4*H                                                          HPCG2290
      CALL FCT(Z,Y,DERY)                                                HPCG2300
      DO 102 I=1,NDIM                                                   HPCG2310
      Z=H*DERY(I)                                                       HPCG2320
      AUX(6,I)=Z                                                        HPCG2330
  102 Y(I)=AUX(N,I)+.2969776*AUX(5,I)+.1587596*Z                        HPCG2340
C                                                                       HPCG2350
      Z=X+.4557372*H                                                    HPCG2360
      CALL FCT(Z,Y,DERY)                                                HPCG2370
      DO 103 I=1,NDIM                                                   HPCG2380
      Z=H*DERY(I)                                                       HPCG2390
      AUX(7,I)=Z                                                        HPCG2400
  103 Y(I)=AUX(N,I)+.2181004*AUX(5,I)-3.050965*AUX(6,I)+3.832865*Z      HPCG2410
C                                                                       HPCG2420
      Z=X+H                                                             HPCG2430
      CALL FCT(Z,Y,DERY)                                                HPCG2440
      DO 104 I=1,NDIM                                                   HPCG2450
  1040Y(I)=AUX(N,I)+.1747603*AUX(5,I)-.5514807*AUX(6,I)                 HPCG2460
     1+1.205536*AUX(7,I)+.1711848*H*DERY(I)                             HPCG2470
      GOTO(9,13,15,21),ISW                                              HPCG2480
C                                                                       HPCG2490
C     POSSIBLE BREAK-POINT FOR LINKAGE                                  HPCG2500
C                                                                       HPCG2510
C     STARTING VALUES ARE COMPUTED.                                     HPCG2520
C     NOW START HAMMINGS MODIFIED PREDICTOR-CORRECTOR METHOD.           HPCG2530
  200 ISTEP=3                                                           HPCG2540
  201 IF(N-8)204,202,204                                                HPCG2550
C                                                                       HPCG2560
C     N=8 CAUSES THE ROWS OF AUX TO CHANGE THEIR STORAGE LOCATIONS      HPCG2570
  202 DO 203 N=2,7                                                      HPCG2580
      DO 203 I=1,NDIM                                                   HPCG2590
      AUX(N-1,I)=AUX(N,I)                                               HPCG2600
  203 AUX(N+6,I)=AUX(N+7,I)                                             HPCG2610
      N=7                                                               HPCG2620
C                                                                       HPCG2630
C     N LESS THAN 8 CAUSES N+1 TO GET N                                 HPCG2640
  204 N=N+1                                                             HPCG2650
C                                                                       HPCG2660
C     COMPUTATION OF NEXT VECTOR Y                                      HPCG2670
      DO 205 I=1,NDIM                                                   HPCG2680
      AUX(N-1,I)=Y(I)                                                   HPCG2690
  205 AUX(N+6,I)=DERY(I)                                                HPCG2700
      X=X+H                                                             HPCG2710
  206 ISTEP=ISTEP+1                                                     HPCG2720
      DO 207 I=1,NDIM                                                   HPCG2730
     0DELT=AUX(N-4,I)+1.333333*H*(AUX(N+6,I)+AUX(N+6,I)-AUX(N+5,I)+     HPCG2740
     1AUX(N+4,I)+AUX(N+4,I))                                            HPCG2750
      Y(I)=DELT-.9256198*AUX(16,I)                                      HPCG2760
  207 AUX(16,I)=DELT                                                    HPCG2770
C     PREDICTOR IS NOW GENERATED IN ROW 16 OF AUX, MODIFIED PREDICTOR   HPCG2780
C     IS GENERATED IN Y. DELT MEANS AN AUXILIARY STORAGE.               HPCG2790
C                                                                       HPCG2800
      CALL FCT(X,Y,DERY)                                                HPCG2810
C     DERIVATIVE OF MODIFIED PREDICTOR IS GENERATED IN DERY             HPCG2820
C                                                                       HPCG2830
      DO 208 I=1,NDIM                                                   HPCG2840
     0DELT=.125*(9.*AUX(N-1,I)-AUX(N-3,I)+3.*H*(DERY(I)+AUX(N+6,I)+     HPCG2850
     1AUX(N+6,I)-AUX(N+5,I)))                                           HPCG2860
      AUX(16,I)=AUX(16,I)-DELT                                          HPCG2870
  208 Y(I)=DELT+.07438017*AUX(16,I)                                     HPCG2880
C                                                                       HPCG2890
C     TEST WHETHER H MUST BE HALVED OR DOUBLED                          HPCG2900
      DELT=0.                                                           HPCG2910
      DO 209 I=1,NDIM                                                   HPCG2920
  209 DELT=DELT+AUX(15,I)*ABS(AUX(16,I))                                HPCG2930
      IF(DELT-PRMT(4))210,222,222                                       HPCG2940
C                                                                       HPCG2950
C     H MUST NOT BE HALVED. THAT MEANS Y(I) ARE GOOD.                   HPCG2960
  210 CALL FCT(X,Y,DERY)                                                HPCG2970
      CALL OUTP(X,Y,DERY,IHLF,NDIM,PRMT)                                HPCG2980
      IF(PRMT(5))212,211,212                                            HPCG2990
  211 IF(IHLF-11)213,212,212                                            HPCG3000
  212 RETURN                                                            HPCG3010
  213 IF(H*(X-PRMT(2)))214,212,212                                      HPCG3020
  214 IF(ABS(X-PRMT(2))-.1*ABS(H))212,215,215                           HPCG3030
  215 IF(DELT-.02*PRMT(4))216,216,201                                   HPCG3040
C                                                                       HPCG3050
C                                                                       HPCG3060
C     H COULD BE DOUBLED IF ALL NECESSARY PRECEEDING VALUES ARE         HPCG3070
C     AVAILABLE                                                         HPCG3080
  216 IF(IHLF)201,201,217                                               HPCG3090
  217 IF(N-7)201,218,218                                                HPCG3100
  218 IF(ISTEP-4)201,219,219                                            HPCG3110
  219 IMOD=ISTEP/2                                                      HPCG3120
      IF(ISTEP-IMOD-IMOD)201,220,201                                    HPCG3130
  220 H=H+H                                                             HPCG3140
      IHLF=IHLF-1                                                       HPCG3150
      ISTEP=0                                                           HPCG3160
      DO 221 I=1,NDIM                                                   HPCG3170
      AUX(N-1,I)=AUX(N-2,I)                                             HPCG3180
      AUX(N-2,I)=AUX(N-4,I)                                             HPCG3190
      AUX(N-3,I)=AUX(N-6,I)                                             HPCG3200
      AUX(N+6,I)=AUX(N+5,I)                                             HPCG3210
      AUX(N+5,I)=AUX(N+3,I)                                             HPCG3220
      AUX(N+4,I)=AUX(N+1,I)                                             HPCG3230
      DELT=AUX(N+6,I)+AUX(N+5,I)                                        HPCG3240
      DELT=DELT+DELT+DELT                                               HPCG3250
  2210AUX(16,I)=8.962963*(Y(I)-AUX(N-3,I))-3.361111*H*(DERY(I)+DELT     HPCG3260
     1+AUX(N+4,I))                                                      HPCG3270
      GOTO 201                                                          HPCG3280
C                                                                       HPCG3290
C                                                                       HPCG3300
C     H MUST BE HALVED                                                  HPCG3310
  222 IHLF=IHLF+1                                                       HPCG3320
      IF(IHLF-10)223,223,210                                            HPCG3330
  223 H=.5*H                                                            HPCG3340
      ISTEP=0                                                           HPCG3350
      DO 224 I=1,NDIM                                                   HPCG3360
     0Y(I)=.00390625*(80.*AUX(N-1,I)+135.*AUX(N-2,I)+40.*AUX(N-3,I)+    HPCG3370
     1AUX(N-4,I))-.1171875*(AUX(N+6,I)-6.*AUX(N+5,I)-AUX(N+4,I))*H      HPCG3380
     0AUX(N-4,I)=.00390625*(12.*AUX(N-1,I)+135.*AUX(N-2,I)+             HPCG3390
     1108.*AUX(N-3,I)+AUX(N-4,I))-.0234375*(AUX(N+6,I)+18.*AUX(N+5,I)-  HPCG3400
     29.*AUX(N+4,I))*H                                                  HPCG3410
      AUX(N-3,I)=AUX(N-2,I)                                             HPCG3420
  224 AUX(N+4,I)=AUX(N+5,I)                                             HPCG3430
      X=X-H                                                             HPCG3440
      DELT=X-(H+H)                                                      HPCG3450
      CALL FCT(DELT,Y,DERY)                                             HPCG3460
      DO 225 I=1,NDIM                                                   HPCG3470
      AUX(N-2,I)=Y(I)                                                   HPCG3480
      AUX(N+5,I)=DERY(I)                                                HPCG3490
  225 Y(I)=AUX(N-4,I)                                                   HPCG3500
      DELT=DELT-(H+H)                                                   HPCG3510
      CALL FCT(DELT,Y,DERY)                                             HPCG3520
      DO 226 I=1,NDIM                                                   HPCG3530
      DELT=AUX(N+5,I)+AUX(N+4,I)                                        HPCG3540
      DELT=DELT+DELT+DELT                                               HPCG3550
     0AUX(16,I)=8.962963*(AUX(N-1,I)-Y(I))-3.361111*H*(AUX(N+6,I)+DELT  HPCG3560
     1+DERY(I))                                                         HPCG3570
  226 AUX(N+3,I)=DERY(I)                                                HPCG3580
      GOTO 206                                                          HPCG3590
      END                                                               HPCG3600
