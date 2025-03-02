C                                                                       RKGS  10
C     ..................................................................RKGS  20
C                                                                       RKGS  30
C        SUBROUTINE RKGS                                                RKGS  40
C                                                                      RKGS  50
C        PURPOSE                                                        RKGS  60
C           TO SOLVE A SYSTEM OF FIRST ORDER ORDINARY DIFFERENTIAL      RKGS  70
C           EQUATIONS WITH GIVEN INITIAL VALUES.                        RKGS  80
C                                                                       RKGS  90
C        USAGE                                                          RKGS 100
C           CALL RKGS (PRMT,Y,DERY,NDIM,IHLF,FCT,OUTP,AUX)              RKGS 110
C           PARAMETERS FCT AND OUTP REQUIRE AN EXTERNAL STATEMENT.      RKGS 120
C                                                                       RKGS 130
C        DESCRIPTION OF PARAMETERS                                      RKGS 140
C           PRMT   - AN INPUT AND OUTPUT VECTOR WITH DIMENSION GREATER  RKGS 150
C                    OR EQUAL TO 5, WHICH SPECIFIES THE PARAMETERS OF   RKGS 160
C                    THE INTERVAL AND OF ACCURACY AND WHICH SERVES FOR  RKGS 170
C                    COMMUNICATION BETWEEN OUTPUT SUBROUTINE (FURNISHED RKGS 180
C                    BY THE USER) AND SUBROUTINE RKGS. EXCEPT PRMT(5)   RKGS 190
C                    THE COMPONENTS ARE NOT DESTROYED BY SUBROUTINE     RKGS 200
C                    RKGS AND THEY ARE                                  RKGS 210
C           PRMT(1)- LOWER BOUND OF THE INTERVAL (INPUT),               RKGS 220
C           PRMT(2)- UPPER BOUND OF THE INTERVAL (INPUT),               RKGS 230
C           PRMT(3)- INITIAL INCREMENT OF THE INDEPENDENT VARIABLE      RKGS 240
C                    (INPUT),                                           RKGS 250
C           PRMT(4)- UPPER ERROR BOUND (INPUT). IF ABSOLUTE ERROR IS    RKGS 260
C                    GREATER THAN PRMT(4), INCREMENT GETS HALVED.       RKGS 270
C                    IF INCREMENT IS LESS THAN PRMT(3) AND ABSOLUTE     RKGS 280
C                    ERROR LESS THAN PRMT(4)/50, INCREMENT GETS DOUBLED.RKGS 290
C                    THE USER MAY CHANGE PRMT(4) BY MEANS OF HIS        RKGS 300
C                    OUTPUT SUBROUTINE.                                 RKGS 310
C           PRMT(5)- NO INPUT PARAMETER. SUBROUTINE RKGS INITIALIZES    RKGS 320
C                    PRMT(5)=0. IF THE USER WANTS TO TERMINATE          RKGS 330
C                    SUBROUTINE RKGS AT ANY OUTPUT POINT, HE HAS TO     RKGS 340
C                    CHANGE PRMT(5) TO NON-ZERO BY MEANS OF SUBROUTINE  RKGS 350
C                    OUTP. FURTHER COMPONENTS OF VECTOR PRMT ARE        RKGS 360
C                    FEASIBLE IF ITS DIMENSION IS DEFINED GREATER       RKGS 370
C                    THAN 5. HOWEVER SUBROUTINE RKGS DOES NOT REQUIRE   RKGS 380
C                    AND CHANGE THEM. NEVERTHELESS THEY MAY BE USEFUL   RKGS 390
C                    FOR HANDING RESULT VALUES TO THE MAIN PROGRAM      RKGS 400
C                    (CALLING RKGS) WHICH ARE OBTAINED BY SPECIAL       RKGS 410
C                    MANIPULATIONS WITH OUTPUT DATA IN SUBROUTINE OUTP. RKGS 420
C           Y      - INPUT VECTOR OF INITIAL VALUES.  (DESTROYED)       RKGS 430
C                    LATERON Y IS THE RESULTING VECTOR OF DEPENDENT     RKGS 440
C                    VARIABLES COMPUTED AT INTERMEDIATE POINTS X.       RKGS 450
C           DERY   - INPUT VECTOR OF ERROR WEIGHTS.  (DESTROYED)        RKGS 460
C                    THE SUM OF ITS COMPONENTS MUST BE EQUAL TO 1.      RKGS 470
C                    LATERON DERY IS THE VECTOR OF DERIVATIVES, WHICH   RKGS 480
C                    BELONG TO FUNCTION VALUES Y AT A POINT X.          RKGS 490
C           NDIM   - AN INPUT VALUE, WHICH SPECIFIES THE NUMBER OF      RKGS 500
C                    EQUATIONS IN THE SYSTEM.                           RKGS 510
C           IHLF   - AN OUTPUT VALUE, WHICH SPECIFIES THE NUMBER OF     RKGS 520
C                    BISECTIONS OF THE INITIAL INCREMENT. IF IHLF GETS  RKGS 530
C                    GREATER THAN 10, SUBROUTINE RKGS RETURNS WITH      RKGS 540
C                    ERROR MESSAGE IHLF=11 INTO MAIN PROGRAM. ERROR     RKGS 550
C                    MESSAGE IHLF=12 OR IHLF=13 APPEARS IN CASE         RKGS 560
C                    PRMT(3)=0 OR IN CASE SIGN(PRMT(3)).NE.SIGN(PRMT(2)-RKGS 570
C                    PRMT(1)) RESPECTIVELY.                             RKGS 580
C           FCT    - THE NAME OF AN EXTERNAL SUBROUTINE USED. THIS      RKGS 590
C                    SUBROUTINE COMPUTES THE RIGHT HAND SIDES DERY OF   RKGS 600
C                    THE SYSTEM TO GIVEN VALUES X AND Y. ITS PARAMETER  RKGS 610
C                    LIST MUST BE X,Y,DERY. SUBROUTINE FCT SHOULD       RKGS 620
C                    NOT DESTROY X AND Y.                               RKGS 630
C           OUTP   - THE NAME OF AN EXTERNAL OUTPUT SUBROUTINE USED.    RKGS 640
C                    ITS PARAMETER LIST MUST BE X,Y,DERY,IHLF,NDIM,PRMT.RKGS 650
C                    NONE OF THESE PARAMETERS (EXCEPT, IF NECESSARY,    RKGS 660
C                    PRMT(4),PRMT(5),...) SHOULD BE CHANGED BY          RKGS 670
C                    SUBROUTINE OUTP. IF PRMT(5) IS CHANGED TO NON-ZERO,RKGS 680
C                    SUBROUTINE RKGS IS TERMINATED.                     RKGS 690
C           AUX    - AN AUXILIARY STORAGE ARRAY WITH 8 ROWS AND NDIM    RKGS 700
C                    COLUMNS.                                           RKGS 710
C                                                                       RKGS 720
C        REMARKS                                                        RKGS 730
C           THE PROCEDURE TERMINATES AND RETURNS TO CALLING PROGRAM, IF RKGS 740
C           (1) MORE THAN 10 BISECTIONS OF THE INITIAL INCREMENT ARE    RKGS 750
C               NECESSARY TO GET SATISFACTORY ACCURACY (ERROR MESSAGE   RKGS 760
C               IHLF=11),                                               RKGS 770
C           (2) INITIAL INCREMENT IS EQUAL TO 0 OR HAS WRONG SIGN       RKGS 780
C               (ERROR MESSAGES IHLF=12 OR IHLF=13),                    RKGS 790
C           (3) THE WHOLE INTEGRATION INTERVAL IS WORKED THROUGH,       RKGS 800
C           (4) SUBROUTINE OUTP HAS CHANGED PRMT(5) TO NON-ZERO.        RKGS 810
C                                                                       RKGS 820
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  RKGS 830
C           THE EXTERNAL SUBROUTINES FCT(X,Y,DERY) AND                  RKGS 840
C           OUTP(X,Y,DERY,IHLF,NDIM,PRMT) MUST BE FURNISHED BY THE USER.RKGS 850
C                                                                       RKGS 860
C        METHOD                                                         RKGS 870
C           EVALUATION IS DONE BY MEANS OF FOURTH ORDER RUNGE-KUTTA     RKGS 880
C           FORMULAE IN THE MODIFICATION DUE TO GILL. ACCURACY IS       RKGS 890
C           TESTED COMPARING THE RESULTS OF THE PROCEDURE WITH SINGLE   RKGS 900
C           AND DOUBLE INCREMENT.                                       RKGS 910
C           SUBROUTINE RKGS AUTOMATICALLY ADJUSTS THE INCREMENT DURING  RKGS 920
C           THE WHOLE COMPUTATION BY HALVING OR DOUBLING. IF MORE THAN  RKGS 930
C           10 BISECTIONS OF THE INCREMENT ARE NECESSARY TO GET         RKGS 940
C           SATISFACTORY ACCURACY, THE SUBROUTINE RETURNS WITH          RKGS 950
C           ERROR MESSAGE IHLF=11 INTO MAIN PROGRAM.                    RKGS 960
C           TO GET FULL FLEXIBILITY IN OUTPUT, AN OUTPUT SUBROUTINE     RKGS 970
C           MUST BE FURNISHED BY THE USER.                              RKGS 980
C           FOR REFERENCE, SEE                                          RKGS 990
C           RALSTON/WILF, MATHEMATICAL METHODS FOR DIGITAL COMPUTERS,   RKGS1000
C           WILEY, NEW YORK/LONDON, 1960, PP.110-120.                   RKGS1010
C                                                                       RKGS1020
C     ..................................................................RKGS1030
C                                                                       RKGS1040
      SUBROUTINE RKGS(PRMT,Y,DERY,NDIM,IHLF,FCT,OUTP,AUX)               RKGS1050
C                                                                       RKGS1060
C                                                                       RKGS1070
      DIMENSION Y(1),DERY(1),AUX(8,1),A(4),B(4),C(4),PRMT(1)            RKGS1080
      DO 1 I=1,NDIM                                                     RKGS1090
    1 AUX(8,I)=.06666667*DERY(I)                                        RKGS1100
      X=PRMT(1)                                                         RKGS1110
      XEND=PRMT(2)                                                      RKGS1120
      H=PRMT(3)                                                         RKGS1130
      PRMT(5)=0.                                                        RKGS1140
      CALL FCT(X,Y,DERY)                                                RKGS1150
C                                                                       RKGS1160
C     ERROR TEST                                                        RKGS1170
      IF(H*(XEND-X))38,37,2                                             RKGS1180
C                                                                       RKGS1190
C     PREPARATIONS FOR RUNGE-KUTTA METHOD                               RKGS1200
    2 A(1)=.5                                                           RKGS1210
      A(2)=.2928932                                                     RKGS1220
      A(3)=1.707107                                                     RKGS1230
      A(4)=.1666667                                                     RKGS1240
      B(1)=2.                                                           RKGS1250
      B(2)=1.                                                           RKGS1260
      B(3)=1.                                                           RKGS1270
      B(4)=2.                                                           RKGS1280
      C(1)=.5                                                           RKGS1290
      C(2)=.2928932                                                     RKGS1300
      C(3)=1.707107                                                     RKGS1310
      C(4)=.5                                                           RKGS1320
C                                                                       RKGS1330
C     PREPARATIONS OF FIRST RUNGE-KUTTA STEP                            RKGS1340
      DO 3 I=1,NDIM                                                     RKGS1350
      AUX(1,I)=Y(I)                                                     RKGS1360
      AUX(2,I)=DERY(I)                                                  RKGS1370
      AUX(3,I)=0.                                                       RKGS1380
    3 AUX(6,I)=0.                                                       RKGS1390
      IREC=0                                                            RKGS1400
      H=H+H                                                             RKGS1410
      IHLF=-1                                                           RKGS1420
      ISTEP=0                                                           RKGS1430
      IEND=0                                                            RKGS1440
C                                                                       RKGS1450
C                                                                       RKGS1460
C     START OF A RUNGE-KUTTA STEP                                       RKGS1470
    4 IF((X+H-XEND)*H)7,6,5                                             RKGS1480
    5 H=XEND-X                                                          RKGS1490
    6 IEND=1                                                            RKGS1500
C                                                                       RKGS1510
C     RECORDING OF INITIAL VALUES OF THIS STEP                          RKGS1520
    7 CALL OUTP(X,Y,DERY,IREC,NDIM,PRMT)                                RKGS1530
      IF(PRMT(5))40,8,40                                                RKGS1540
    8 ITEST=0                                                           RKGS1550
    9 ISTEP=ISTEP+1                                                     RKGS1560
C                                                                       RKGS1570
C                                                                       RKGS1580
C     START OF INNERMOST RUNGE-KUTTA LOOP                               RKGS1590
      J=1                                                               RKGS1600
   10 AJ=A(J)                                                           RKGS1610
      BJ=B(J)                                                           RKGS1620
      CJ=C(J)                                                           RKGS1630
      DO 11 I=1,NDIM                                                    RKGS1640
      R1=H*DERY(I)                                                      RKGS1650
      R2=AJ*(R1-BJ*AUX(6,I))                                            RKGS1660
      Y(I)=Y(I)+R2                                                      RKGS1670
      R2=R2+R2+R2                                                       RKGS1680
   11 AUX(6,I)=AUX(6,I)+R2-CJ*R1                                        RKGS1690
      IF(J-4)12,15,15                                                   RKGS1700
   12 J=J+1                                                             RKGS1710
      IF(J-3)13,14,13                                                   RKGS1720
   13 X=X+.5*H                                                          RKGS1730
   14 CALL FCT(X,Y,DERY)                                                RKGS1740
      GOTO 10                                                           RKGS1750
C     END OF INNERMOST RUNGE-KUTTA LOOP                                 RKGS1760
C                                                                       RKGS1770
C                                                                       RKGS1780
C     TEST OF ACCURACY                                                  RKGS1790
   15 IF(ITEST)16,16,20                                                 RKGS1800
C                                                                       RKGS1810
C     IN CASE ITEST=0 THERE IS NO POSSIBILITY FOR TESTING OF ACCURACY   RKGS1820
   16 DO 17 I=1,NDIM                                                    RKGS1830
   17 AUX(4,I)=Y(I)                                                     RKGS1840
      ITEST=1                                                           RKGS1850
      ISTEP=ISTEP+ISTEP-2                                               RKGS1860
   18 IHLF=IHLF+1                                                       RKGS1870
      X=X-H                                                             RKGS1880
      H=.5*H                                                            RKGS1890
      DO 19 I=1,NDIM                                                    RKGS1900
      Y(I)=AUX(1,I)                                                     RKGS1910
      DERY(I)=AUX(2,I)                                                  RKGS1920
   19 AUX(6,I)=AUX(3,I)                                                 RKGS1930
      GOTO 9                                                            RKGS1940
C                                                                       RKGS1950
C     IN CASE ITEST=1 TESTING OF ACCURACY IS POSSIBLE                   RKGS1960
   20 IMOD=ISTEP/2                                                      RKGS1970
      IF(ISTEP-IMOD-IMOD)21,23,21                                       RKGS1980
   21 CALL FCT(X,Y,DERY)                                                RKGS1990
      DO 22 I=1,NDIM                                                    RKGS2000
      AUX(5,I)=Y(I)                                                     RKGS2010
   22 AUX(7,I)=DERY(I)                                                  RKGS2020
      GOTO 9                                                            RKGS2030
C                                                                       RKGS2040
C     COMPUTATION OF TEST VALUE DELT                                    RKGS2050
   23 DELT=0.                                                           RKGS2060
      DO 24 I=1,NDIM                                                    RKGS2070
   24 DELT=DELT+AUX(8,I)*ABS(AUX(4,I)-Y(I))                             RKGS2080
      IF(DELT-PRMT(4))28,28,25                                          RKGS2090
C                                                                       RKGS2100
C     ERROR IS TOO GREAT                                                RKGS2110
   25 IF(IHLF-10)26,36,36                                               RKGS2120
   26 DO 27 I=1,NDIM                                                    RKGS2130
   27 AUX(4,I)=AUX(5,I)                                                 RKGS2140
      ISTEP=ISTEP+ISTEP-4                                               RKGS2150
      X=X-H                                                             RKGS2160
      IEND=0                                                            RKGS2170
      GOTO 18                                                           RKGS2180
C                                                                       RKGS2190
C     RESULT VALUES ARE GOOD                                            RKGS2200
   28 CALL FCT(X,Y,DERY)                                                RKGS2210
      DO 29 I=1,NDIM                                                    RKGS2220
      AUX(1,I)=Y(I)                                                     RKGS2230
      AUX(2,I)=DERY(I)                                                  RKGS2240
      AUX(3,I)=AUX(6,I)                                                 RKGS2250
      Y(I)=AUX(5,I)                                                     RKGS2260
   29 DERY(I)=AUX(7,I)                                                  RKGS2270
      CALL OUTP(X-H,Y,DERY,IHLF,NDIM,PRMT)                              RKGS2280
      IF(PRMT(5))40,30,40                                               RKGS2290
   30 DO 31 I=1,NDIM                                                    RKGS2300
      Y(I)=AUX(1,I)                                                     RKGS2310
   31 DERY(I)=AUX(2,I)                                                  RKGS2320
      IREC=IHLF                                                         RKGS2330
      IF(IEND)32,32,39                                                  RKGS2340
C                                                                       RKGS2350
C     INCREMENT GETS DOUBLED                                            RKGS2360
   32 IHLF=IHLF-1                                                       RKGS2370
      ISTEP=ISTEP/2                                                     RKGS2380
      H=H+H                                                             RKGS2390
      IF(IHLF)4,33,33                                                   RKGS2400
   33 IMOD=ISTEP/2                                                      RKGS2410
      IF(ISTEP-IMOD-IMOD)4,34,4                                         RKGS2420
   34 IF(DELT-.02*PRMT(4))35,35,4                                       RKGS2430
   35 IHLF=IHLF-1                                                       RKGS2440
      ISTEP=ISTEP/2                                                     RKGS2450
      H=H+H                                                             RKGS2460
      GOTO 4                                                            RKGS2470
C                                                                       RKGS2480
C                                                                       RKGS2490
C     RETURNS TO CALLING PROGRAM                                        RKGS2500
   36 IHLF=11                                                           RKGS2510
      CALL FCT(X,Y,DERY)                                                RKGS2520
      GOTO 39                                                           RKGS2530
   37 IHLF=12                                                           RKGS2540
      GOTO 39                                                           RKGS2550
   38 IHLF=13                                                           RKGS2560
   39 CALL OUTP(X,Y,DERY,IHLF,NDIM,PRMT)                                RKGS2570
   40 RETURN                                                            RKGS2580
      END                                                               RKGS2590
