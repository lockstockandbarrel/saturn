C                                                                       HPCL  10
C     ..................................................................HPCL  20
C                                                                       HPCL  30
C        SUBROUTINE HPCL                                                HPCL  40
C                                                                       HPCL  50
C        PURPOSE                                                        HPCL  60
C           TO SOLVE A SYSTEM OF FIRST ORDER ORDINARY LINEAR            HPCL  70
C           DIFFERENTIAL EQUATIONS WITH GIVEN INITIAL VALUES.           HPCL  80
C                                                                       HPCL  90
C        USAGE                                                          HPCL 100
C           CALL HPCL (PRMT,Y,DERY,NDIM,IHLF,AFCT,FCT,OUTP,AUX,A)       HPCL 110
C           PARAMETERS AFCT,FCT AND OUTP REQUIRE AN EXTERNAL STATEMENT. HPCL 120
C                                                                       HPCL 130
C        DESCRIPTION OF PARAMETERS                                      HPCL 140
C           PRMT   - AN INPUT AND OUTPUT VECTOR WITH DIMENSION GREATER  HPCL 150
C                    OR EQUAL TO 5, WHICH SPECIFIES THE PARAMETERS OF   HPCL 160
C                    THE INTERVAL AND OF ACCURACY AND WHICH SERVES FOR  HPCL 170
C                    COMMUNICATION BETWEEN OUTPUT SUBROUTINE (FURNISHED HPCL 180
C                    BY THE USER) AND SUBROUTINE HPCL. EXCEPT PRMT(5)   HPCL 190
C                    THE COMPONENTS ARE NOT DESTROYED BY SUBROUTINE     HPCL 200
C                    HPCL AND THEY ARE                                  HPCL 210
C           PRMT(1)- LOWER BOUND OF THE INTERVAL (INPUT),               HPCL 220
C           PRMT(2)- UPPER BOUND OF THE INTERVAL (INPUT),               HPCL 230
C           PRMT(3)- INITIAL INCREMENT OF THE INDEPENDENT VARIABLE      HPCL 240
C                    (INPUT),                                           HPCL 250
C           PRMT(4)- UPPER ERROR BOUND (INPUT). IF ABSOLUTE ERROR IS    HPCL 260
C                    GREATER THAN PRMT(4), INCREMENT GETS HALVED.       HPCL 270
C                    IF INCREMENT IS LESS THAN PRMT(3) AND ABSOLUTE     HPCL 280
C                    ERROR LESS THAN PRMT(4)/50, INCREMENT GETS DOUBLED.HPCL 290
C                    THE USER MAY CHANGE PRMT(4) BY MEANS OF HIS        HPCL 300
C                    OUTPUT SUBROUTINE.                                 HPCL 310
C           PRMT(5)- NO INPUT PARAMETER. SUBROUTINE HPCL INITIALIZES    HPCL 320
C                    PRMT(5)=0. IF THE USER WANTS TO TERMINATE          HPCL 330
C                    SUBROUTINE HPCL AT ANY OUTPUT POINT, HE HAS TO     HPCL 340
C                    CHANGE PRMT(5) TO NON-ZERO BY MEANS OF SUBROUTINE  HPCL 350
C                    OUTP. FURTHER COMPONENTS OF VECTOR PRMT ARE        HPCL 360
C                    FEASIBLE IF ITS DIMENSION IS DEFINED GREATER       HPCL 370
C                    THAN 5. HOWEVER SUBROUTINE HPCL DOES NOT REQUIRE   HPCL 380
C                    AND CHANGE THEM. NEVERTHELESS THEY MAY BE USEFUL   HPCL 390
C                    FOR HANDING RESULT VALUES TO THE MAIN PROGRAM      HPCL 400
C                    (CALLING HPCL) WHICH ARE OBTAINED BY SPECIAL       HPCL 410
C                    MANIPULATIONS WITH OUTPUT DATA IN SUBROUTINE OUTP. HPCL 420
C           Y      - INPUT VECTOR OF INITIAL VALUES.  (DESTROYED)       HPCL 430
C                    LATERON Y IS THE RESULTING VECTOR OF DEPENDENT     HPCL 440
C                    VARIABLES COMPUTED AT INTERMEDIATE POINTS X.       HPCL 450
C           DERY   - INPUT VECTOR OF ERROR WEIGHTS.  (DESTROYED)        HPCL 460
C                    THE SUM OF ITS COMPONENTS MUST BE EQUAL TO 1.      HPCL 470
C                    LATERON DERY IS THE VECTOR OF DERIVATIVES, WHICH   HPCL 480
C                    BELONG TO FUNCTION VALUES Y AT A POINT X.          HPCL 490
C           NDIM   - AN INPUT VALUE, WHICH SPECIFIES THE NUMBER OF      HPCL 500
C                    EQUATIONS IN THE SYSTEM.                           HPCL 510
C           IHLF   - AN OUTPUT VALUE, WHICH SPECIFIES THE NUMBER OF     HPCL 520
C                    BISECTIONS OF THE INITIAL INCREMENT. IF IHLF GETS  HPCL 530
C                    GREATER THAN 10, SUBROUTINE HPCL RETURNS WITH      HPCL 540
C                    ERROR MESSAGE IHLF=11 INTO MAIN PROGRAM.           HPCL 550
C                    ERROR MESSAGE IHLF=12 OR IHLF=13 APPEARS IN CASE   HPCL 560
C                    PRMT(3)=0 OR IN CASE SIGN(PRMT(3)).NE.SIGN(PRMT(2)-HPCL 570
C                    PRMT(1)) RESPECTIVELY.                             HPCL 580
C           AFCT   - THE NAME OF AN EXTERNAL SUBROUTINE USED. IT        HPCL 590
C                    COMPUTES MATRIX A (FACTOR OF VECTOR Y ON THE       HPCL 600
C                    RIGHT HAND SIDE OF THE SYSTEM) FOR A GIVEN X-VALUE.HPCL 610
C                    ITS PARAMETER LIST MUST BE X,A. THE SUBROUTINE     HPCL 620
C                    SHOULD NOT DESTROY X.                              HPCL 630
C           FCT    - THE NAME OF AN EXTERNAL SUBROUTINE USED. IT        HPCL 640
C                    COMPUTES VECTOR F (INHOMOGENEOUS PART OF THE       HPCL 650
C                    RIGHT HAND SIDE OF THE SYSTEM) FOR A GIVEN X-VALUE.HPCL 660
C                    ITS PARAMETER LIST MUST BE X,F. THE SUBROUTINE     HPCL 670
C                    SHOULD NOT DESTROY X.                              HPCL 680
C           OUTP   - THE NAME OF AN EXTERNAL OUTPUT SUBROUTINE USED.    HPCL 690
C                    ITS PARAMETER LIST MUST BE X,Y,DERY,IHLF,NDIM,PRMT.HPCL 700
C                    NONE OF THESE PARAMETERS (EXCEPT, IF NECESSARY,    HPCL 710
C                    PRMT(4),PRMT(5),...) SHOULD BE CHANGED BY          HPCL 720
C                    SUBROUTINE OUTP. IF PRMT(5) IS CHANGED TO NON-ZERO,HPCL 730
C                    SUBROUTINE HPCL IS TERMINATED.                     HPCL 740
C           AUX    - AN AUXILIARY STORAGE ARRAY WITH 16 ROWS AND NDIM   HPCL 750
C                    COLUMNS.                                           HPCL 760
C           A      - AN NDIM BY NDIM MATRIX, WHICH IS USED AS AUXILIARY HPCL 770
C                    STORAGE ARRAY.                                     HPCL 780
C                                                                       HPCL 790
C        REMARKS                                                        HPCL 800
C           THE PROCEDURE TERMINATES AND RETURNS TO CALLING PROGRAM, IF HPCL 810
C           (1) MORE THAN 10 BISECTIONS OF THE INITIAL INCREMENT ARE    HPCL 820
C               NECESSARY TO GET SATISFACTORY ACCURACY (ERROR MESSAGE   HPCL 830
C               IHLF=11),                                               HPCL 840
C           (2) INITIAL INCREMENT IS EQUAL TO 0 OR HAS WRONG SIGN       HPCL 850
C               (ERROR MESSAGES IHLF=12 OR IHLF=13),                    HPCL 860
C           (3) THE WHOLE INTEGRATION INTERVAL IS WORKED THROUGH,       HPCL 870
C           (4) SUBROUTINE OUTP HAS CHANGED PRMT(5) TO NON-ZERO.        HPCL 880
C                                                                       HPCL 890
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  HPCL 900
C           THE EXTERNAL SUBROUTINES AFCT(X,A), FCT(X,F) AND            HPCL 910
C           OUTP(X,Y,DERY,IHLF,NDIM,PRMT) MUST BE FURNISHED BY THE USER.HPCL 920
C                                                                       HPCL 930
C        METHOD                                                         HPCL 940
C           EVALUATION IS DONE BY MEANS OF HAMMINGS MODIFIED PREDICTOR- HPCL 950
C           CORRECTOR METHOD. IT IS A FOURTH ORDER METHOD, USING 4      HPCL 960
C           PRECEEDING POINTS FOR COMPUTATION OF A NEW VECTOR Y OF THE  HPCL 970
C           DEPENDENT VARIABLES.                                        HPCL 980
C           FOURTH ORDER RUNGE-KUTTA METHOD SUGGESTED BY RALSTON IS     HPCL 990
C           USED FOR ADJUSTMENT OF THE INITIAL INCREMENT AND FOR        HPCL1000
C           COMPUTATION OF STARTING VALUES.                             HPCL1010
C           SUBROUTINE HPCL AUTOMATICALLY ADJUSTS THE INCREMENT DURING  HPCL1020
C           THE WHOLE COMPUTATION BY HALVING OR DOUBLING.               HPCL1030
C           TO GET FULL FLEXIBILITY IN OUTPUT, AN OUTPUT SUBROUTINE     HPCL1040
C           MUST BE CODED BY THE USER.                                  HPCL1050
C           FOR REFERENCE, SEE                                          HPCL1060
C           (1)  RALSTON/WILF, MATHEMATICAL METHODS FOR DIGITAL         HPCL1070
C                COMPUTERS, WILEY, NEW YORK/LONDON, 1960, PP.95-109.    HPCL1080
C           (2)  RALSTON, RUNGE-KUTTA METHODS WITH MINIMUM ERROR BOUNDS,HPCL1090
C                MTAC, VOL.16, ISS.80 (1962), PP.431-437.               HPCL1100
C                                                                       HPCL1110
C     ..................................................................HPCL1120
C                                                                       HPCL1130
      SUBROUTINE HPCL(PRMT,Y,DERY,NDIM,IHLF,AFCT,FCT,OUTP,AUX,A)        HPCL1140
C                                                                       HPCL1150
C                                                                       HPCL1160
C     THE FOLLOWING FIRST PART OF SUBROUTINE HPCL (UNTIL FIRST BREAK-   HPCL1170
C     POINT FOR LINKAGE) HAS TO STAY IN CORE DURING THE WHOLE           HPCL1180
C     COMPUTATION                                                       HPCL1190
C                                                                       HPCL1200
      DIMENSION PRMT(1),Y(1),DERY(1),AUX(16,1),A(1)                     HPCL1210
      GOTO 100                                                          HPCL1220
C                                                                       HPCL1230
C     THIS PART OF SUBROUTINE HPCL COMPUTES THE RIGHT HAND SIDE DERY OF HPCL1240
C     THE GIVEN SYSTEM OF LINEAR DIFFERENTIAL EQUATIONS.                HPCL1250
    1 CALL AFCT(X,A)                                                    HPCL1260
      CALL FCT(X,DERY)                                                  HPCL1270
      DO 3 M=1,NDIM                                                     HPCL1280
      LL=M-NDIM                                                         HPCL1290
      HS=0.                                                             HPCL1300
      DO 2 L=1,NDIM                                                     HPCL1310
      LL=LL+NDIM                                                        HPCL1320
    2 HS=HS+A(LL)*Y(L)                                                  HPCL1330
    3 DERY(M)=HS+DERY(M)                                                HPCL1340
      GOTO(105,202,204,206,115,122,125,308,312,327,329,128),ISW2        HPCL1350
C                                                                       HPCL1360
C     POSSIBLE BREAK-POINT FOR LINKAGE                                  HPCL1370
C                                                                       HPCL1380
  100 N=1                                                               HPCL1390
      IHLF=0                                                            HPCL1400
      X=PRMT(1)                                                         HPCL1410
      H=PRMT(3)                                                         HPCL1420
      PRMT(5)=0.                                                        HPCL1430
      DO 101 I=1,NDIM                                                   HPCL1440
      AUX(16,I)=0.                                                      HPCL1450
      AUX(15,I)=DERY(I)                                                 HPCL1460
  101 AUX(1,I)=Y(I)                                                     HPCL1470
      IF(H*(PRMT(2)-X))103,102,104                                      HPCL1480
C                                                                       HPCL1490
C     ERROR RETURNS                                                     HPCL1500
  102 IHLF=12                                                           HPCL1510
      GOTO 104                                                          HPCL1520
  103 IHLF=13                                                           HPCL1530
C                                                                       HPCL1540
C     COMPUTATION OF DERY FOR STARTING VALUES                           HPCL1550
  104 ISW2=1                                                            HPCL1560
      GOTO 1                                                            HPCL1570
C                                                                       HPCL1580
C     RECORDING OF STARTING VALUES                                      HPCL1590
  105 CALL OUTP(X,Y,DERY,IHLF,NDIM,PRMT)                                HPCL1600
      IF(PRMT(5))107,106,107                                            HPCL1610
  106 IF(IHLF)108,108,107                                               HPCL1620
  107 RETURN                                                            HPCL1630
  108 DO 109 I=1,NDIM                                                   HPCL1640
  109 AUX(8,I)=DERY(I)                                                  HPCL1650
C                                                                       HPCL1660
C     COMPUTATION OF AUX(2,I)                                           HPCL1670
      ISW1=1                                                            HPCL1680
      GOTO 200                                                          HPCL1690
C                                                                       HPCL1700
  110 X=X+H                                                             HPCL1710
      DO 111 I=1,NDIM                                                   HPCL1720
  111 AUX(2,I)=Y(I)                                                     HPCL1730
C                                                                       HPCL1740
C     INCREMENT H IS TESTED BY MEANS OF BISECTION                       HPCL1750
  112 IHLF=IHLF+1                                                       HPCL1760
      X=X-H                                                             HPCL1770
      DO 113 I=1,NDIM                                                   HPCL1780
  113 AUX(4,I)=AUX(2,I)                                                 HPCL1790
      H=.5*H                                                            HPCL1800
      N=1                                                               HPCL1810
      ISW1=2                                                            HPCL1820
      GOTO 200                                                          HPCL1830
C                                                                       HPCL1840
  114 X=X+H                                                             HPCL1850
      ISW2=5                                                            HPCL1860
      GOTO 1                                                            HPCL1870
  115 N=2                                                               HPCL1880
      DO 116 I=1,NDIM                                                   HPCL1890
      AUX(2,I)=Y(I)                                                     HPCL1900
  116 AUX(9,I)=DERY(I)                                                  HPCL1910
      ISW1=3                                                            HPCL1920
      GOTO 200                                                          HPCL1930
C                                                                       HPCL1940
C     COMPUTATION OF TEST VALUE DELT                                    HPCL1950
  117 DELT=0.                                                           HPCL1960
      DO 118 I=1,NDIM                                                   HPCL1970
  118 DELT=DELT+AUX(15,I)*ABS(Y(I)-AUX(4,I))                            HPCL1980
      DELT=.06666667*DELT                                               HPCL1990
      IF(DELT-PRMT(4))121,121,119                                       HPCL2000
  119 IF(IHLF-10)112,120,120                                            HPCL2010
C                                                                       HPCL2020
C     NO SATISFACTORY ACCURACY AFTER 10 BISECTIONS. ERROR MESSAGE.      HPCL2030
  120 IHLF=11                                                           HPCL2040
      X=X+H                                                             HPCL2050
      GOTO 104                                                          HPCL2060
C                                                                       HPCL2070
C     SATISFACTORY ACCURACY AFTER LESS THAN 11 BISECTIONS               HPCL2080
  121 X=X+H                                                             HPCL2090
      ISW2=6                                                            HPCL2100
      GOTO 1                                                            HPCL2110
  122 DO 123 I=1,NDIM                                                   HPCL2120
      AUX(3,I)=Y(I)                                                     HPCL2130
  123 AUX(10,I)=DERY(I)                                                 HPCL2140
      N=3                                                               HPCL2150
      ISW1=4                                                            HPCL2160
      GOTO 200                                                          HPCL2170
C                                                                       HPCL2180
  124 N=1                                                               HPCL2190
      X=X+H                                                             HPCL2200
      ISW2=7                                                            HPCL2210
      GOTO 1                                                            HPCL2220
  125 X=PRMT(1)                                                         HPCL2230
      DO 126 I=1,NDIM                                                   HPCL2240
      AUX(11,I)=DERY(I)                                                 HPCL2250
  1260Y(I)=AUX(1,I)+H*(.375*AUX(8,I)+.7916667*AUX(9,I)                  HPCL2260
     1-.2083333*AUX(10,I)+.04166667*DERY(I))                            HPCL2270
  127 X=X+H                                                             HPCL2280
      N=N+1                                                             HPCL2290
      ISW2=12                                                           HPCL2300
      GOTO 1                                                            HPCL2310
  128 CALL OUTP(X,Y,DERY,IHLF,NDIM,PRMT)                                HPCL2320
      IF(PRMT(5))107,129,107                                            HPCL2330
  129 IF(N-4)130,300,300                                                HPCL2340
  130 DO 131 I=1,NDIM                                                   HPCL2350
      AUX(N,I)=Y(I)                                                     HPCL2360
  131 AUX(N+7,I)=DERY(I)                                                HPCL2370
      IF(N-3)132,134,300                                                HPCL2380
C                                                                       HPCL2390
  132 DO 133 I=1,NDIM                                                   HPCL2400
      DELT=AUX(9,I)+AUX(9,I)                                            HPCL2410
      DELT=DELT+DELT                                                    HPCL2420
  133 Y(I)=AUX(1,I)+.3333333*H*(AUX(8,I)+DELT+AUX(10,I))                HPCL2430
      GOTO 127                                                          HPCL2440
C                                                                       HPCL2450
  134 DO 135 I=1,NDIM                                                   HPCL2460
      DELT=AUX(9,I)+AUX(10,I)                                           HPCL2470
      DELT=DELT+DELT+DELT                                               HPCL2480
  135 Y(I)=AUX(1,I)+.375*H*(AUX(8,I)+DELT+AUX(11,I))                    HPCL2490
      GOTO 127                                                          HPCL2500
C                                                                       HPCL2510
C     THE FOLLOWING PART OF SUBROUTINE HPCL COMPUTES BY MEANS OF        HPCL2520
C     RUNGE-KUTTA METHOD STARTING VALUES FOR THE NOT SELF-STARTING      HPCL2530
C     PREDICTOR-CORRECTOR METHOD.                                       HPCL2540
  200 Z=X                                                               HPCL2550
      DO 201 I=1,NDIM                                                   HPCL2560
      X=H*AUX(N+7,I)                                                    HPCL2570
      AUX(5,I)=X                                                        HPCL2580
  201 Y(I)=AUX(N,I)+.4*X                                                HPCL2590
C     X IS AN AUXILIARY STORAGE LOCATION                                HPCL2600
C                                                                       HPCL2610
      X=Z+.4*H                                                          HPCL2620
      ISW2=2                                                            HPCL2630
      GOTO 1                                                            HPCL2640
  202 DO 203 I=1,NDIM                                                   HPCL2650
      X=H*DERY(I)                                                       HPCL2660
      AUX(6,I)=X                                                        HPCL2670
  203 Y(I)=AUX(N,I)+.2969776*AUX(5,I)+.1587596*X                        HPCL2680
C                                                                       HPCL2690
      X=Z+.4557372*H                                                    HPCL2700
      ISW2=3                                                            HPCL2710
      GOTO 1                                                            HPCL2720
  204 DO 205 I=1,NDIM                                                   HPCL2730
      X=H*DERY(I)                                                       HPCL2740
      AUX(7,I)=X                                                        HPCL2750
  205 Y(I)=AUX(N,I)+.2181004*AUX(5,I)-3.050965*AUX(6,I)+3.832865*X      HPCL2760
C                                                                       HPCL2770
      X=Z+H                                                             HPCL2780
      ISW2=4                                                            HPCL2790
      GOTO 1                                                            HPCL2800
  206 DO 207 I=1,NDIM                                                   HPCL2810
  2070Y(I)=AUX(N,I)+.1747603*AUX(5,I)-.5514807*AUX(6,I)                 HPCL2820
     1+1.205536*AUX(7,I)+.1711848*H*DERY(I)                             HPCL2830
      X=Z                                                               HPCL2840
      GOTO(110,114,117,124),ISW1                                        HPCL2850
C                                                                       HPCL2860
C     POSSIBLE BREAK-POINT FOR LINKAGE                                  HPCL2870
C                                                                       HPCL2880
C     STARTING VALUES ARE COMPUTED.                                     HPCL2890
C     NOW START HAMMINGS MODIFIED PREDICTOR-CORRECTOR METHOD.           HPCL2900
  300 ISTEP=3                                                           HPCL2910
  301 IF(N-8)304,302,304                                                HPCL2920
C                                                                       HPCL2930
C     N=8 CAUSES THE ROWS OF AUX TO CHANGE THEIR STORAGE LOCATIONS      HPCL2940
  302 DO 303 N=2,7                                                      HPCL2950
      DO 303 I=1,NDIM                                                   HPCL2960
      AUX(N-1,I)=AUX(N,I)                                               HPCL2970
  303 AUX(N+6,I)=AUX(N+7,I)                                             HPCL2980
      N=7                                                               HPCL2990
C                                                                       HPCL3000
C     N LESS THAN 8 CAUSES N+1 TO GET N                                 HPCL3010
  304 N=N+1                                                             HPCL3020
C                                                                       HPCL3030
C     COMPUTATION OF NEXT VECTOR Y                                      HPCL3040
      DO 305 I=1,NDIM                                                   HPCL3050
      AUX(N-1,I)=Y(I)                                                   HPCL3060
  305 AUX(N+6,I)=DERY(I)                                                HPCL3070
      X=X+H                                                             HPCL3080
  306 ISTEP=ISTEP+1                                                     HPCL3090
      DO 307 I=1,NDIM                                                   HPCL3100
     0DELT=AUX(N-4,I)+1.333333*H*(AUX(N+6,I)+AUX(N+6,I)-AUX(N+5,I)+     HPCL3110
     1AUX(N+4,I)+AUX(N+4,I))                                            HPCL3120
      Y(I)=DELT-.9256198*AUX(16,I)                                      HPCL3130
  307 AUX(16,I)=DELT                                                    HPCL3140
C     PREDICTOR IS NOW GENERATED IN ROW 16 OF AUX, MODIFIED PREDICTOR   HPCL3150
C     IS GENERATED IN Y. DELT MEANS AN AUXILIARY STORAGE.               HPCL3160
      ISW2=8                                                            HPCL3170
      GOTO 1                                                            HPCL3180
C     DERIVATIVE OF MODIFIED PREDICTOR IS GENERATED IN DERY             HPCL3190
C                                                                       HPCL3200
  308 DO 309 I=1,NDIM                                                   HPCL3210
     0DELT=.125*(9.*AUX(N-1,I)-AUX(N-3,I)+3.*H*(DERY(I)+AUX(N+6,I)+     HPCL3220
     1AUX(N+6,I)-AUX(N+5,I)))                                           HPCL3230
      AUX(16,I)=AUX(16,I)-DELT                                          HPCL3240
  309 Y(I)=DELT+.07438017*AUX(16,I)                                     HPCL3250
C                                                                       HPCL3260
C     TEST WHETHER H MUST BE HALVED OR DOUBLED                          HPCL3270
      DELT=0.                                                           HPCL3280
      DO 310 I=1,NDIM                                                   HPCL3290
  310 DELT=DELT+AUX(15,I)*ABS(AUX(16,I))                                HPCL3300
      IF(DELT-PRMT(4))311,324,324                                       HPCL3310
C                                                                       HPCL3320
C     H MUST NOT BE HALVED. THAT MEANS Y(I) ARE GOOD.                   HPCL3330
  311 ISW2=9                                                            HPCL3340
      GOTO 1                                                            HPCL3350
  312 CALL OUTP(X,Y,DERY,IHLF,NDIM,PRMT)                                HPCL3360
      IF(PRMT(5))314,313,314                                            HPCL3370
  313 IF(IHLF-11)315,314,314                                            HPCL3380
  314 RETURN                                                            HPCL3390
  315 IF(H*(X-PRMT(2)))316,314,314                                      HPCL3400
  316 IF(ABS(X-PRMT(2))-.1*ABS(H))314,317,317                           HPCL3410
  317 IF(DELT-.02*PRMT(4))318,318,301                                   HPCL3420
C                                                                       HPCL3430
C                                                                       HPCL3440
C     H COULD BE DOUBLED IF ALL NECESSARY PRECEEDING VALUES ARE         HPCL3450
C     AVAILABLE                                                         HPCL3460
  318 IF(IHLF)301,301,319                                               HPCL3470
  319 IF(N-7)301,320,320                                                HPCL3480
  320 IF(ISTEP-4)301,321,321                                            HPCL3490
  321 IMOD=ISTEP/2                                                      HPCL3500
      IF(ISTEP-IMOD-IMOD)301,322,301                                    HPCL3510
  322 H=H+H                                                             HPCL3520
      IHLF=IHLF-1                                                       HPCL3530
      ISTEP=0                                                           HPCL3540
      DO 323 I=1,NDIM                                                   HPCL3550
      AUX(N-1,I)=AUX(N-2,I)                                             HPCL3560
      AUX(N-2,I)=AUX(N-4,I)                                             HPCL3570
      AUX(N-3,I)=AUX(N-6,I)                                             HPCL3580
      AUX(N+6,I)=AUX(N+5,I)                                             HPCL3590
      AUX(N+5,I)=AUX(N+3,I)                                             HPCL3600
      AUX(N+4,I)=AUX(N+1,I)                                             HPCL3610
      DELT=AUX(N+6,I)+AUX(N+5,I)                                        HPCL3620
      DELT=DELT+DELT+DELT                                               HPCL3630
  3230AUX(16,I)=8.962963*(Y(I)-AUX(N-3,I))-3.361111*H*(DERY(I)+DELT     HPCL3640
     1+AUX(N+4,I))                                                      HPCL3650
      GOTO 301                                                          HPCL3660
C                                                                       HPCL3670
C                                                                       HPCL3680
C     H MUST BE HALVED                                                  HPCL3690
  324 IHLF=IHLF+1                                                       HPCL3700
      IF(IHLF-10)325,325,311                                            HPCL3710
  325 H=.5*H                                                            HPCL3720
      ISTEP=0                                                           HPCL3730
      DO 326 I=1,NDIM                                                   HPCL3740
     0Y(I)=.00390625*(80.*AUX(N-1,I)+135.*AUX(N-2,I)+40.*AUX(N-3,I)+    HPCL3750
     1AUX(N-4,I))-.1171875*(AUX(N+6,I)-6.*AUX(N+5,I)-AUX(N+4,I))*H      HPCL3760
     0AUX(N-4,I)=.00390625*(12.*AUX(N-1,I)+135.*AUX(N-2,I)+             HPCL3770
     1108.*AUX(N-3,I)+AUX(N-4,I))-.0234375*(AUX(N+6,I)+18.*AUX(N+5,I)-  HPCL3780
     29.*AUX(N+4,I))*H                                                  HPCL3790
      AUX(N-3,I)=AUX(N-2,I)                                             HPCL3800
  326 AUX(N+4,I)=AUX(N+5,I)                                             HPCL3810
      DELT=X-H                                                          HPCL3820
      X=DELT-(H+H)                                                      HPCL3830
      ISW2=10                                                           HPCL3840
      GOTO 1                                                            HPCL3850
  327 DO 328 I=1,NDIM                                                   HPCL3860
      AUX(N-2,I)=Y(I)                                                   HPCL3870
      AUX(N+5,I)=DERY(I)                                                HPCL3880
  328 Y(I)=AUX(N-4,I)                                                   HPCL3890
      X=X-(H+H)                                                         HPCL3900
      ISW2=11                                                           HPCL3910
      GOTO 1                                                            HPCL3920
  329 X=DELT                                                            HPCL3930
      DO 330 I=1,NDIM                                                   HPCL3940
      DELT=AUX(N+5,I)+AUX(N+4,I)                                        HPCL3950
      DELT=DELT+DELT+DELT                                               HPCL3960
     0AUX(16,I)=8.962963*(AUX(N-1,I)-Y(I))-3.361111*H*(AUX(N+6,I)+DELT  HPCL3970
     1+DERY(I))                                                         HPCL3980
  330 AUX(N+3,I)=DERY(I)                                                HPCL3990
      GOTO 306                                                          HPCL4000
      END                                                               HPCL4010
