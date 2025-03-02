C                                                                       DHCL  10
C     ..................................................................DHCL  20
C                                                                       DHCL  30
C        SUBROUTINE DHPCL                                               DHCL  40
C                                                                       DHCL  50
C        PURPOSE                                                        DHCL  60
C           TO SOLVE A SYSTEM OF FIRST ORDER ORDINARY LINEAR            DHCL  70
C           DIFFERENTIAL EQUATIONS WITH GIVEN INITIAL VALUES.           DHCL  80
C                                                                       DHCL  90
C        USAGE                                                          DHCL 100
C           CALL DHPCL (PRMT,Y,DERY,NDIM,IHLF,AFCT,FCT,OUTP,AUX,A)      DHCL 110
C           PARAMETERS AFCT,FCT AND OUTP REQUIRE AN EXTERNAL STATEMENT. DHCL 120
C                                                                       DHCL 130
C        DESCRIPTION OF PARAMETERS                                      DHCL 140
C           PRMT   - DOUBLE PRECISION INPUT AND OUTPUT VECTOR WITH      DHCL 150
C                    DIMENSION GREATER THAN OR EQUAL TO 5, WHICH        DHCL 160
C                    SPECIFIES THE PARAMETERS OF THE INTERVAL AND OF    DHCL 170
C                    ACCURACY AND WHICH SERVES FOR COMMUNICATION BETWEENDHCL 180
C                    OUTPUT SUBROUTINE (FURNISHED BY THE USER) AND      DHCL 190
C                    SUBROUTINE DHPCL. EXCEPT PRMT(5) THE COMPONENTS    DHCL 200
C                    ARE NOT DESTROYED BY SUBROUTINE DHPCL AND THEY ARE DHCL 210
C           PRMT(1)- LOWER BOUND OF THE INTERVAL (INPUT),               DHCL 220
C           PRMT(2)- UPPER BOUND OF THE INTERVAL (INPUT),               DHCL 230
C           PRMT(3)- INITIAL INCREMENT OF THE INDEPENDENT VARIABLE      DHCL 240
C                    (INPUT),                                           DHCL 250
C           PRMT(4)- UPPER ERROR BOUND (INPUT). IF ABSOLUTE ERROR IS    DHCL 260
C                    GREATER THAN PRMT(4), INCREMENT GETS HALVED.       DHCL 270
C                    IF INCREMENT IS LESS THAN PRMT(3) AND ABSOLUTE     DHCL 280
C                    ERROR LESS THAN PRMT(4)/50, INCREMENT GETS DOUBLED.DHCL 290
C                    THE USER MAY CHANGE PRMT(4) BY MEANS OF HIS        DHCL 300
C                    OUTPUT SUBROUTINE.                                 DHCL 310
C           PRMT(5)- NO INPUT PARAMETER. SUBROUTINE DHPCL INITIALIZES   DHCL 320
C                    PRMT(5)=0. IF THE USER WANTS TO TERMINATE          DHCL 330
C                    SUBROUTINE DHPCL AT ANY OUTPUT POINT, HE HAS TO    DHCL 340
C                    CHANGE PRMT(5) TO NON-ZERO BY MEANS OF SUBROUTINE  DHCL 350
C                    OUTP. FURTHER COMPONENTS OF VECTOR PRMT ARE        DHCL 360
C                    FEASIBLE IF ITS DIMENSION IS DEFINED GREATER       DHCL 370
C                    THAN 5. HOWEVER SUBROUTINE DHPCL DOES NOT REQUIRE  DHCL 380
C                    AND CHANGE THEM. NEVERTHELESS THEY MAY BE USEFUL   DHCL 390
C                    FOR HANDING RESULT VALUES TO THE MAIN PROGRAM      DHCL 400
C                    (CALLING DHPCL) WHICH ARE OBTAINED BY SPECIAL      DHCL 410
C                    MANIPULATIONS WITH OUTPUT DATA IN SUBROUTINE OUTP. DHCL 420
C           Y      - DOUBLE PRECISION INPUT VECTOR OF INITIAL VALUES    DHCL 430
C                    (DESTROYED). LATERON Y IS THE RESULTING VECTOR OF  DHCL 440
C                    DEPENDENT VARIABLES COMPUTED AT INTERMEDIATE       DHCL 450
C                    POINTS X.                                          DHCL 460
C           DERY   - DOUBLE PRECISION INPUT VECTOR OF ERROR WEIGHTS     DHCL 470
C                    (DESTROYED). THE SUM OF ITS COMPONENTS MUST BE     DHCL 480
C                    EQUAL TO 1. LATERON DERY IS THE VECTOR OF          DHCL 490
C                    DERIVATIVES, WHICH BELONG TO FUNCTION VALUES Y AT  DHCL 500
C                    INTERMEDIATE POINTS X.                             DHCL 510
C           NDIM   - AN INPUT VALUE, WHICH SPECIFIES THE NUMBER OF      DHCL 520
C                    EQUATIONS IN THE SYSTEM.                           DHCL 530
C           IHLF   - AN OUTPUT VALUE, WHICH SPECIFIES THE NUMBER OF     DHCL 540
C                    BISECTIONS OF THE INITIAL INCREMENT. IF IHLF GETS  DHCL 550
C                    GREATER THAN 10, SUBROUTINE DHPCL RETURNS WITH     DHCL 560
C                    ERROR MESSAGE IHLF=11 INTO MAIN PROGRAM.           DHCL 570
C                    ERROR MESSAGE IHLF=12 OR IHLF=13 APPEARS IN CASE   DHCL 580
C                    PRMT(3)=0 OR IN CASE SIGN(PRMT(3)).NE.SIGN(PRMT(2)-DHCL 590
C                    PRMT(1)) RESPECTIVELY.                             DHCL 600
C           AFCT   - THE NAME OF AN EXTERNAL SUBROUTINE USED. IT        DHCL 610
C                    COMPUTES MATRIX A (FACTOR OF VECTOR Y ON THE       DHCL 620
C                    RIGHT HAND SIDE OF THE SYSTEM) FOR A GIVEN X-VALUE.DHCL 630
C                    ITS PARAMETER LIST MUST BE X,A. THE SUBROUTINE     DHCL 640
C                    SHOULD NOT DESTROY X.                              DHCL 650
C           FCT    - THE NAME OF AN EXTERNAL SUBROUTINE USED. IT        DHCL 660
C                    COMPUTES VECTOR F (INHOMOGENEOUS PART OF THE       DHCL 670
C                    RIGHT HAND SIDE OF THE SYSTEM) FOR A GIVEN X-VALUE.DHCL 680
C                    ITS PARAMETER LIST MUST BE X,F. THE SUBROUTINE     DHCL 690
C                    SHOULD NOT DESTROY X.                              DHCL 700
C           OUTP   - THE NAME OF AN EXTERNAL OUTPUT SUBROUTINE USED.    DHCL 710
C                    ITS PARAMETER LIST MUST BE X,Y,DERY,IHLF,NDIM,PRMT.DHCL 720
C                    NONE OF THESE PARAMETERS (EXCEPT, IF NECESSARY,    DHCL 730
C                    PRMT(4),PRMT(5),...) SHOULD BE CHANGED BY          DHCL 740
C                    SUBROUTINE OUTP. IF PRMT(5) IS CHANGED TO NON-ZERO,DHCL 750
C                    SUBROUTINE DHPCL IS TERMINATED.                    DHCL 760
C           AUX    - DOUBLE PRECISION AUXILIARY STORAGE ARRAY WITH 16   DHCL 770
C                    ROWS AND NDIM COLUMNS.                             DHCL 780
C           A      - DOUBLE PRECISION NDIM BY NDIM MATRIX, WHICH IS USEDDHCL 790
C                    AS AUXILIARY STORAGE ARRAY.                        DHCL 800
C                                                                       DHCL 810
C        REMARKS                                                        DHCL 820
C           THE PROCEDURE TERMINATES AND RETURNS TO CALLING PROGRAM, IF DHCL 830
C           (1) MORE THAN 10 BISECTIONS OF THE INITIAL INCREMENT ARE    DHCL 840
C               NECESSARY TO GET SATISFACTORY ACCURACY (ERROR MESSAGE   DHCL 850
C               IHLF=11),                                               DHCL 860
C           (2) INITIAL INCREMENT IS EQUAL TO 0 OR HAS WRONG SIGN       DHCL 870
C               (ERROR MESSAGES IHLF=12 OR IHLF=13),                    DHCL 880
C           (3) THE WHOLE INTEGRATION INTERVAL IS WORKED THROUGH,       DHCL 890
C           (4) SUBROUTINE OUTP HAS CHANGED PRMT(5) TO NON-ZERO.        DHCL 900
C                                                                       DHCL 910
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DHCL 920
C           THE EXTERNAL SUBROUTINES AFCT(X,A), FCT(X,F) AND            DHCL 930
C           OUTP(X,Y,DERY,IHLF,NDIM,PRMT) MUST BE FURNISHED BY THE USER.DHCL 940
C                                                                       DHCL 950
C        METHOD                                                         DHCL 960
C           EVALUATION IS DONE BY MEANS OF HAMMINGS MODIFIED PREDICTOR- DHCL 970
C           CORRECTOR METHOD. IT IS A FOURTH ORDER METHOD, USING 4      DHCL 980
C           PRECEEDING POINTS FOR COMPUTATION OF A NEW VECTOR Y OF THE  DHCL 990
C           DEPENDENT VARIABLES.                                        DHCL1000
C           FOURTH ORDER RUNGE-KUTTA METHOD SUGGESTED BY RALSTON IS     DHCL1010
C           USED FOR ADJUSTMENT OF THE INITIAL INCREMENT AND FOR        DHCL1020
C           COMPUTATION OF STARTING VALUES.                             DHCL1030
C           SUBROUTINE DHPCL AUTOMATICALLY ADJUSTS THE INCREMENT DURING DHCL1040
C           THE WHOLE COMPUTATION BY HALVING OR DOUBLING.               DHCL1050
C           TO GET FULL FLEXIBILITY IN OUTPUT, AN OUTPUT SUBROUTINE     DHCL1060
C           MUST BE CODED BY THE USER.                                  DHCL1070
C           FOR REFERENCE, SEE                                          DHCL1080
C           (1)  RALSTON/WILF, MATHEMATICAL METHODS FOR DIGITAL         DHCL1090
C                COMPUTERS, WILEY, NEW YORK/LONDON, 1960, PP.95-109.    DHCL1100
C           (2)  RALSTON, RUNGE-KUTTA METHODS WITH MINIMUM ERROR BOUNDS,DHCL1110
C                MTAC, VOL.16, ISS.80 (1962), PP.431-437.               DHCL1120
C                                                                       DHCL1130
C     ..................................................................DHCL1140
C                                                                       DHCL1150
      SUBROUTINE DHPCL(PRMT,Y,DERY,NDIM,IHLF,AFCT,FCT,OUTP,AUX,A)       DHCL1160
C                                                                       DHCL1170
C                                                                       DHCL1180
C     THE FOLLOWING FIRST PART OF SUBROUTINE DHPCL (UNTIL FIRST BREAK-  DHCL1190
C     POINT FOR LINKAGE) HAS TO STAY IN CORE DURING THE WHOLE           DHCL1200
C     COMPUTATION                                                       DHCL1210
C                                                                       DHCL1220
      DIMENSION PRMT(1),Y(1),DERY(1),AUX(16,1),A(1)                     DHCL1230
      DOUBLE PRECISION PRMT,Y,DERY,AUX,X,H,Z,DELT,A,HS                  DHCL1240
      GOTO 100                                                          DHCL1250
C                                                                       DHCL1260
C     THIS PART OF SUBROUTINE DHPCL COMPUTES THE RIGHT HAND SIDE DERY OFDHCL1270
C     THE GIVEN SYSTEM OF LINEAR DIFFERENTIAL EQUATIONS.                DHCL1280
    1 CALL AFCT(X,A)                                                    DHCL1290
      CALL FCT(X,DERY)                                                  DHCL1300
      DO 3 M=1,NDIM                                                     DHCL1310
      LL=M-NDIM                                                         DHCL1320
      HS=0.D0                                                           DHCL1330
      DO 2 L=1,NDIM                                                     DHCL1340
      LL=LL+NDIM                                                        DHCL1350
    2 HS=HS+A(LL)*Y(L)                                                  DHCL1360
    3 DERY(M)=HS+DERY(M)                                                DHCL1370
      GOTO(105,202,204,206,115,122,125,308,312,327,329,128),ISW2        DHCL1380
C                                                                       DHCL1390
C     POSSIBLE BREAK-POINT FOR LINKAGE                                  DHCL1400
C                                                                       DHCL1410
  100 N=1                                                               DHCL1420
      IHLF=0                                                            DHCL1430
      X=PRMT(1)                                                         DHCL1440
      H=PRMT(3)                                                         DHCL1450
      PRMT(5)=0.D0                                                      DHCL1460
      DO 101 I=1,NDIM                                                   DHCL1470
      AUX(16,I)=0.D0                                                    DHCL1480
      AUX(15,I)=DERY(I)                                                 DHCL1490
  101 AUX(1,I)=Y(I)                                                     DHCL1500
      IF(H*(PRMT(2)-X))103,102,104                                      DHCL1510
C                                                                       DHCL1520
C     ERROR RETURNS                                                     DHCL1530
  102 IHLF=12                                                           DHCL1540
      GOTO 104                                                          DHCL1550
  103 IHLF=13                                                           DHCL1560
C                                                                       DHCL1570
C     COMPUTATION OF DERY FOR STARTING VALUES                           DHCL1580
  104 ISW2=1                                                            DHCL1590
      GOTO 1                                                            DHCL1600
C                                                                       DHCL1610
C     RECORDING OF STARTING VALUES                                      DHCL1620
  105 CALL OUTP(X,Y,DERY,IHLF,NDIM,PRMT)                                DHCL1630
      IF(PRMT(5))107,106,107                                            DHCL1640
  106 IF(IHLF)108,108,107                                               DHCL1650
  107 RETURN                                                            DHCL1660
  108 DO 109 I=1,NDIM                                                   DHCL1670
C                                                                       DHCL1680
  109 AUX(8,I)=DERY(I)                                                  DHCL1690
C     COMPUTATION OF AUX(2,I)                                           DHCL1700
      ISW1=1                                                            DHCL1710
      GOTO 200                                                          DHCL1720
  110 X=X+H                                                             DHCL1730
      DO 111 I=1,NDIM                                                   DHCL1740
  111 AUX(2,I)=Y(I)                                                     DHCL1750
C                                                                       DHCL1760
C     INCREMENT H IS TESTED BY MEANS OF BISECTION                       DHCL1770
  112 IHLF=IHLF+1                                                       DHCL1780
      X=X-H                                                             DHCL1790
      DO 113 I=1,NDIM                                                   DHCL1800
  113 AUX(4,I)=AUX(2,I)                                                 DHCL1810
      H=.5D0*H                                                          DHCL1820
      N=1                                                               DHCL1830
      ISW1=2                                                            DHCL1840
      GOTO 200                                                          DHCL1850
C                                                                       DHCL1860
  114 X=X+H                                                             DHCL1870
      ISW2=5                                                            DHCL1880
      GOTO 1                                                            DHCL1890
  115 N=2                                                               DHCL1900
      DO 116 I=1,NDIM                                                   DHCL1910
      AUX(2,I)=Y(I)                                                     DHCL1920
  116 AUX(9,I)=DERY(I)                                                  DHCL1930
      ISW1=3                                                            DHCL1940
      GOTO 200                                                          DHCL1950
C                                                                       DHCL1960
C     COMPUTATION OF TEST VALUE DELT                                    DHCL1970
  117 DELT=0.D0                                                         DHCL1980
      DO 118 I=1,NDIM                                                   DHCL1990
  118 DELT=DELT+AUX(15,I)*DABS(Y(I)-AUX(4,I))                           DHCL2000
      DELT=.066666666666666667D0*DELT                                   DHCL2010
      IF(DELT-PRMT(4))121,121,119                                       DHCL2020
  119 IF(IHLF-10)112,120,120                                            DHCL2030
C                                                                       DHCL2040
C     NO SATISFACTORY ACCURACY AFTER 10 BISECTIONS. ERROR MESSAGE.      DHCL2050
  120 IHLF=11                                                           DHCL2060
      X=X+H                                                             DHCL2070
      GOTO 104                                                          DHCL2080
C                                                                       DHCL2090
C     SATISFACTORY ACCURACY AFTER LESS THAN 11 BISECTIONS               DHCL2100
  121 X=X+H                                                             DHCL2110
      ISW2=6                                                            DHCL2120
      GOTO 1                                                            DHCL2130
  122 DO 123 I=1,NDIM                                                   DHCL2140
      AUX(3,I)=Y(I)                                                     DHCL2150
  123 AUX(10,I)=DERY(I)                                                 DHCL2160
      N=3                                                               DHCL2170
      ISW1=4                                                            DHCL2180
      GOTO 200                                                          DHCL2190
C                                                                       DHCL2200
  124 N=1                                                               DHCL2210
      X=X+H                                                             DHCL2220
      ISW2=7                                                            DHCL2230
      GOTO 1                                                            DHCL2240
  125 X=PRMT(1)                                                         DHCL2250
      DO 126 I=1,NDIM                                                   DHCL2260
      AUX(11,I)=DERY(I)                                                 DHCL2270
  1260Y(I)=AUX(1,I)+H*(.375D0*AUX(8,I)+.7916666666666667D0*AUX(9,I)     DHCL2280
     1-.20833333333333333D0*AUX(10,I)+.041666666666666667D0*DERY(I))    DHCL2290
  127 X=X+H                                                             DHCL2300
      N=N+1                                                             DHCL2310
      ISW2=12                                                           DHCL2320
      GOTO 1                                                            DHCL2330
  128 CALL OUTP(X,Y,DERY,IHLF,NDIM,PRMT)                                DHCL2340
      IF(PRMT(5))107,129,107                                            DHCL2350
  129 IF(N-4)130,300,300                                                DHCL2360
  130 DO 131 I=1,NDIM                                                   DHCL2370
      AUX(N,I)=Y(I)                                                     DHCL2380
  131 AUX(N+7,I)=DERY(I)                                                DHCL2390
      IF(N-3)132,134,300                                                DHCL2400
C                                                                       DHCL2410
  132 DO 133 I=1,NDIM                                                   DHCL2420
      DELT=AUX(9,I)+AUX(9,I)                                            DHCL2430
      DELT=DELT+DELT                                                    DHCL2440
  133 Y(I)=AUX(1,I)+.33333333333333333D0*H*(AUX(8,I)+DELT+AUX(10,I))    DHCL2450
      GOTO 127                                                          DHCL2460
C                                                                       DHCL2470
  134 DO 135 I=1,NDIM                                                   DHCL2480
      DELT=AUX(9,I)+AUX(10,I)                                           DHCL2490
      DELT=DELT+DELT+DELT                                               DHCL2500
  135 Y(I)=AUX(1,I)+.375D0*H*(AUX(8,I)+DELT+AUX(11,I))                  DHCL2510
      GOTO 127                                                          DHCL2520
C                                                                       DHCL2530
C     THE FOLLOWING PART OF SUBROUTINE DHPCL COMPUTES BY MEANS OF       DHCL2540
C     RUNGE-KUTTA METHOD STARTING VALUES FOR THE NOT SELF-STARTING      DHCL2550
C     PREDICTOR-CORRECTOR METHOD.                                       DHCL2560
  200 Z=X                                                               DHCL2570
      DO 201 I=1,NDIM                                                   DHCL2580
      X=H*AUX(N+7,I)                                                    DHCL2590
      AUX(5,I)=X                                                        DHCL2600
  201 Y(I)=AUX(N,I)+.4D0*X                                              DHCL2610
C     X IS AN AUXILIARY STORAGE LOCATION                                DHCL2620
C                                                                       DHCL2630
      X=Z+.4D0*H                                                        DHCL2640
      ISW2=2                                                            DHCL2650
      GOTO 1                                                            DHCL2660
  202 DO 203 I=1,NDIM                                                   DHCL2670
      X=H*DERY(I)                                                       DHCL2680
      AUX(6,I)=X                                                        DHCL2690
  203 Y(I)=AUX(N,I)+.29697760924775360D0*AUX(5,I)+.15875964497103583D0*XDHCL2700
C                                                                       DHCL2710
      X=Z+.45573725421878943D0*H                                        DHCL2720
      ISW2=3                                                            DHCL2730
      GOTO 1                                                            DHCL2740
  204 DO 205 I=1,NDIM                                                   DHCL2750
      X=H*DERY(I)                                                       DHCL2760
      AUX(7,I)=X                                                        DHCL2770
  205 Y(I)=AUX(N,I)+.21810038822592047D0*AUX(5,I)-3.0509651486929308D0* DHCL2780
     1AUX(6,I)+3.8328647604670103D0*X                                   DHCL2790
C                                                                       DHCL2800
      X=Z+H                                                             DHCL2810
      ISW2=4                                                            DHCL2820
      GOTO 1                                                            DHCL2830
  206 DO 207 I=1,NDIM                                                   DHCL2840
  2070Y(I)=AUX(N,I)+.17476028226269037D0*AUX(5,I)-.55148066287873294D0* DHCL2850
     1AUX(6,I)+1.2055355993965235D0*AUX(7,I)+.17118478121951903D0*      DHCL2860
     2H*DERY(I)                                                         DHCL2870
      X=Z                                                               DHCL2880
      GOTO(110,114,117,124),ISW1                                        DHCL2890
C                                                                       DHCL2900
C     POSSIBLE BREAK-POINT FOR LINKAGE                                  DHCL2910
C                                                                       DHCL2920
C     STARTING VALUES ARE COMPUTED.                                     DHCL2930
C     NOW START HAMMINGS MODIFIED PREDICTOR-CORRECTOR METHOD.           DHCL2940
  300 ISTEP=3                                                           DHCL2950
  301 IF(N-8)304,302,304                                                DHCL2960
C                                                                       DHCL2970
C     N=8 CAUSES THE ROWS OF AUX TO CHANGE THEIR STORAGE LOCATIONS      DHCL2980
  302 DO 303 N=2,7                                                      DHCL2990
      DO 303 I=1,NDIM                                                   DHCL3000
      AUX(N-1,I)=AUX(N,I)                                               DHCL3010
  303 AUX(N+6,I)=AUX(N+7,I)                                             DHCL3020
      N=7                                                               DHCL3030
C                                                                       DHCL3040
C     N LESS THAN 8 CAUSES N+1 TO GET N                                 DHCL3050
  304 N=N+1                                                             DHCL3060
C                                                                       DHCL3070
C     COMPUTATION OF NEXT VECTOR Y                                      DHCL3080
      DO 305 I=1,NDIM                                                   DHCL3090
      AUX(N-1,I)=Y(I)                                                   DHCL3100
  305 AUX(N+6,I)=DERY(I)                                                DHCL3110
      X=X+H                                                             DHCL3120
  306 ISTEP=ISTEP+1                                                     DHCL3130
      DO 307 I=1,NDIM                                                   DHCL3140
     0DELT=AUX(N-4,I)+1.3333333333333333D0*H*(AUX(N+6,I)+AUX(N+6,I)-    DHCL3150
     1AUX(N+5,I)+AUX(N+4,I)+AUX(N+4,I))                                 DHCL3160
      Y(I)=DELT-.9256198347107438D0*AUX(16,I)                           DHCL3170
  307 AUX(16,I)=DELT                                                    DHCL3180
C     PREDICTOR IS NOW GENERATED IN ROW 16 OF AUX, MODIFIED PREDICTOR   DHCL3190
C     IS GENERATED IN Y. DELT MEANS AN AUXILIARY STORAGE.               DHCL3200
      ISW2=8                                                            DHCL3210
      GOTO 1                                                            DHCL3220
C     DERIVATIVE OF MODIFIED PREDICTOR IS GENERATED IN DERY             DHCL3230
C                                                                       DHCL3240
  308 DO 309 I=1,NDIM                                                   DHCL3250
     0DELT=.125D0*(9.D0*AUX(N-1,I)-AUX(N-3,I)+3.D0*H*(DERY(I)+AUX(N+6,I)DHCL3260
     1+AUX(N+6,I)-AUX(N+5,I)))                                          DHCL3270
      AUX(16,I)=AUX(16,I)-DELT                                          DHCL3280
  309 Y(I)=DELT+.07438016528925620D0*AUX(16,I)                          DHCL3290
C                                                                       DHCL3300
C     TEST WHETHER H MUST BE HALVED OR DOUBLED                          DHCL3310
      DELT=0.D0                                                         DHCL3320
      DO 310 I=1,NDIM                                                   DHCL3330
  310 DELT=DELT+AUX(15,I)*DABS(AUX(16,I))                               DHCL3340
      IF(DELT-PRMT(4))311,324,324                                       DHCL3350
C                                                                       DHCL3360
C     H MUST NOT BE HALVED. THAT MEANS Y(I) ARE GOOD.                   DHCL3370
  311 ISW2=9                                                            DHCL3380
      GOTO 1                                                            DHCL3390
  312 CALL OUTP(X,Y,DERY,IHLF,NDIM,PRMT)                                DHCL3400
      IF(PRMT(5))314,313,314                                            DHCL3410
  313 IF(IHLF-11)315,314,314                                            DHCL3420
  314 RETURN                                                            DHCL3430
  315 IF(H*(X-PRMT(2)))316,314,314                                      DHCL3440
  316 IF(DABS(X-PRMT(2))-.1D0*DABS(H))314,317,317                       DHCL3450
  317 IF(DELT-.02D0*PRMT(4))318,318,301                                 DHCL3460
C                                                                       DHCL3470
C     H COULD BE DOUBLED IF ALL NECESSARY PRECEEDING VALUES ARE         DHCL3480
C     AVAILABLE                                                         DHCL3490
  318 IF(IHLF)301,301,319                                               DHCL3500
  319 IF(N-7)301,320,320                                                DHCL3510
  320 IF(ISTEP-4)301,321,321                                            DHCL3520
  321 IMOD=ISTEP/2                                                      DHCL3530
      IF(ISTEP-IMOD-IMOD)301,322,301                                    DHCL3540
  322 H=H+H                                                             DHCL3550
      IHLF=IHLF-1                                                       DHCL3560
      ISTEP=0                                                           DHCL3570
      DO 323 I=1,NDIM                                                   DHCL3580
      AUX(N-1,I)=AUX(N-2,I)                                             DHCL3590
      AUX(N-2,I)=AUX(N-4,I)                                             DHCL3600
      AUX(N-3,I)=AUX(N-6,I)                                             DHCL3610
      AUX(N+6,I)=AUX(N+5,I)                                             DHCL3620
      AUX(N+5,I)=AUX(N+3,I)                                             DHCL3630
      AUX(N+4,I)=AUX(N+1,I)                                             DHCL3640
      DELT=AUX(N+6,I)+AUX(N+5,I)                                        DHCL3650
      DELT=DELT+DELT+DELT                                               DHCL3660
  3230AUX(16,I)=8.962962962962963D0*(Y(I)-AUX(N-3,I))                   DHCL3670
     1-3.3611111111111111D0*H*(DERY(I)+DELT+AUX(N+4,I))                 DHCL3680
      GOTO 301                                                          DHCL3690
C                                                                       DHCL3700
C     H MUST BE HALVED                                                  DHCL3710
  324 IHLF=IHLF+1                                                       DHCL3720
      IF(IHLF-10)325,325,311                                            DHCL3730
  325 H=.5D0*H                                                          DHCL3740
      ISTEP=0                                                           DHCL3750
      DO 326 I=1,NDIM                                                   DHCL3760
     0Y(I)=.390625D-2*(8.D1*AUX(N-1,I)+135.D0*AUX(N-2,I)+4.D1*AUX(N-3,I)DHCL3770
     1+AUX(N-4,I))-.1171875D0*(AUX(N+6,I)-6.D0*AUX(N+5,I)-AUX(N+4,I))*H DHCL3780
     0AUX(N-4,I)=.390625D-2*(12.D0*AUX(N-1,I)+135.D0*AUX(N-2,I)+        DHCL3790
     1108.D0*AUX(N-3,I)+AUX(N-4,I))-.0234375D0*(AUX(N+6,I)+             DHCL3800
     218.D0*AUX(N+5,I)-9.D0*AUX(N+4,I))*H                               DHCL3810
      AUX(N-3,I)=AUX(N-2,I)                                             DHCL3820
  326 AUX(N+4,I)=AUX(N+5,I)                                             DHCL3830
      DELT=X-H                                                          DHCL3840
      X=DELT-(H+H)                                                      DHCL3850
      ISW2=10                                                           DHCL3860
      GOTO 1                                                            DHCL3870
  327 DO 328 I=1,NDIM                                                   DHCL3880
      AUX(N-2,I)=Y(I)                                                   DHCL3890
      AUX(N+5,I)=DERY(I)                                                DHCL3900
  328 Y(I)=AUX(N-4,I)                                                   DHCL3910
      X=X-(H+H)                                                         DHCL3920
      ISW2=11                                                           DHCL3930
      GOTO 1                                                            DHCL3940
  329 X=DELT                                                            DHCL3950
      DO 330 I=1,NDIM                                                   DHCL3960
      DELT=AUX(N+5,I)+AUX(N+4,I)                                        DHCL3970
      DELT=DELT+DELT+DELT                                               DHCL3980
     0AUX(16,I)=8.962962962962963D0*(AUX(N-1,I)-Y(I))                   DHCL3990
     1-3.3611111111111111D0*H*(AUX(N+6,I)+DELT+DERY(I))                 DHCL4000
  330 AUX(N+3,I)=DERY(I)                                                DHCL4010
      GOTO 306                                                          DHCL4020
      END                                                               DHCL4030
