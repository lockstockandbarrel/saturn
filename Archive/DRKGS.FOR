C                                                                       DRKG  10
C     ..................................................................DRKG  20
C                                                                       DRKG  30
C        SUBROUTINE DRKGS                                               DRKG  40
C                                                                       DRKG  50
C        PURPOSE                                                        DRKG  60
C           TO SOLVE A SYSTEM OF FIRST ORDER ORDINARY DIFFERENTIAL      DRKG  70
C           EQUATIONS WITH GIVEN INITIAL VALUES.                        DRKG  80
C                                                                       DRKG  90
C        USAGE                                                          DRKG 100
C           CALL DRKGS (PRMT,Y,DERY,NDIM,IHLF,FCT,OUTP,AUX)             DRKG 110
C           PARAMETERS FCT AND OUTP REQUIRE AN EXTERNAL STATEMENT.      DRKG 120
C                                                                       DRKG 130
C        DESCRIPTION OF PARAMETERS                                      DRKG 140
C           PRMT   - DOUBLE PRECISION INPUT AND OUTPUT VECTOR WITH      DRKG 150
C                    DIMENSION GREATER THAN OR EQUAL TO 5, WHICH        DRKG 160
C                    SPECIFIES THE PARAMETERS OF THE INTERVAL AND OF    DRKG 170
C                    ACCURACY AND WHICH SERVES FOR COMMUNICATION BETWEENDRKG 180
C                    OUTPUT SUBROUTINE (FURNISHED BY THE USER) AND      DRKG 190
C                    SUBROUTINE DRKGS. EXCEPT PRMT(5) THE COMPONENTS    DRKG 200
C                    ARE NOT DESTROYED BY SUBROUTINE DRKGS AND THEY ARE DRKG 210
C           PRMT(1)- LOWER BOUND OF THE INTERVAL (INPUT),               DRKG 220
C           PRMT(2)- UPPER BOUND OF THE INTERVAL (INPUT),               DRKG 230
C           PRMT(3)- INITIAL INCREMENT OF THE INDEPENDENT VARIABLE      DRKG 240
C                    (INPUT),                                           DRKG 250
C           PRMT(4)- UPPER ERROR BOUND (INPUT). IF ABSOLUTE ERROR IS    DRKG 260
C                    GREATER THAN PRMT(4), INCREMENT GETS HALVED.       DRKG 270
C                    IF INCREMENT IS LESS THAN PRMT(3) AND ABSOLUTE     DRKG 280
C                    ERROR LESS THAN PRMT(4)/50, INCREMENT GETS DOUBLED.DRKG 290
C                    THE USER MAY CHANGE PRMT(4) BY MEANS OF HIS        DRKG 300
C                    OUTPUT SUBROUTINE.                                 DRKG 310
C           PRMT(5)- NO INPUT PARAMETER. SUBROUTINE DRKGS INITIALIZES   DRKG 320
C                    PRMT(5)=0. IF THE USER WANTS TO TERMINATE          DRKG 330
C                    SUBROUTINE DRKGS AT ANY OUTPUT POINT, HE HAS TO    DRKG 340
C                    CHANGE PRMT(5) TO NON-ZERO BY MEANS OF SUBROUTINE  DRKG 350
C                    OUTP. FURTHER COMPONENTS OF VECTOR PRMT ARE        DRKG 360
C                    FEASIBLE IF ITS DIMENSION IS DEFINED GREATER       DRKG 370
C                    THAN 5. HOWEVER SUBROUTINE DRKGS DOES NOT REQUIRE  DRKG 380
C                    AND CHANGE THEM. NEVERTHELESS THEY MAY BE USEFUL   DRKG 390
C                    FOR HANDING RESULT VALUES TO THE MAIN PROGRAM      DRKG 400
C                    (CALLING DRKGS) WHICH ARE OBTAINED BY SPECIAL      DRKG 410
C                    MANIPULATIONS WITH OUTPUT DATA IN SUBROUTINE OUTP. DRKG 420
C           Y      - DOUBLE PRECISION INPUT VECTOR OF INITIAL VALUES    DRKG 430
C                    (DESTROYED). LATERON Y IS THE RESULTING VECTOR OF  DRKG 440
C                    DEPENDENT VARIABLES COMPUTED AT INTERMEDIATE       DRKG 450
C                    POINTS X.                                          DRKG 460
C           DERY   - DOUBLE PRECISION INPUT VECTOR OF ERROR WEIGHTS     DRKG 470
C                    (DESTROYED). THE SUM OF ITS COMPONENTS MUST BE     DRKG 480
C                    EQUAL TO 1. LATERON DERY IS THE VECTOR OF          DRKG 490
C                    DERIVATIVES, WHICH BELONG TO FUNCTION VALUES Y AT  DRKG 500
C                    INTERMEDIATE POINTS X.                             DRKG 510
C           NDIM   - AN INPUT VALUE, WHICH SPECIFIES THE NUMBER OF      DRKG 520
C                    EQUATIONS IN THE SYSTEM.                           DRKG 530
C           IHLF   - AN OUTPUT VALUE, WHICH SPECIFIES THE NUMBER OF     DRKG 540
C                    BISECTIONS OF THE INITIAL INCREMENT. IF IHLF GETS  DRKG 550
C                    GREATER THAN 10, SUBROUTINE DRKGS RETURNS WITH     DRKG 560
C                    ERROR MESSAGE IHLF=11 INTO MAIN PROGRAM. ERROR     DRKG 570
C                    MESSAGE IHLF=12 OR IHLF=13 APPEARS IN CASE         DRKG 580
C                    PRMT(3)=0 OR IN CASE SIGN(PRMT(3)).NE.SIGN(PRMT(2)-DRKG 590
C                    PRMT(1)) RESPECTIVELY.                             DRKG 600
C           FCT    - THE NAME OF AN EXTERNAL SUBROUTINE USED. THIS      DRKG 610
C                    SUBROUTINE COMPUTES THE RIGHT HAND SIDES DERY OF   DRKG 620
C                    THE SYSTEM TO GIVEN VALUES X AND Y. ITS PARAMETER  DRKG 630
C                    LIST MUST BE X,Y,DERY. SUBROUTINE FCT SHOULD       DRKG 640
C                    NOT DESTROY X AND Y.                               DRKG 650
C           OUTP   - THE NAME OF AN EXTERNAL OUTPUT SUBROUTINE USED.    DRKG 660
C                    ITS PARAMETER LIST MUST BE X,Y,DERY,IHLF,NDIM,PRMT.DRKG 670
C                    NONE OF THESE PARAMETERS (EXCEPT, IF NECESSARY,    DRKG 680
C                    PRMT(4),PRMT(5),...) SHOULD BE CHANGED BY          DRKG 690
C                    SUBROUTINE OUTP. IF PRMT(5) IS CHANGED TO NON-ZERO,DRKG 700
C                    SUBROUTINE DRKGS IS TERMINATED.                    DRKG 710
C           AUX    - DOUBLE PRECISION AUXILIARY STORAGE ARRAY WITH 8    DRKG 720
C                    ROWS AND NDIM COLUMNS.                             DRKG 730
C                                                                       DRKG 740
C        REMARKS                                                        DRKG 750
C           THE PROCEDURE TERMINATES AND RETURNS TO CALLING PROGRAM, IF DRKG 760
C           (1) MORE THAN 10 BISECTIONS OF THE INITIAL INCREMENT ARE    DRKG 770
C               NECESSARY TO GET SATISFACTORY ACCURACY (ERROR MESSAGE   DRKG 780
C               IHLF=11),                                               DRKG 790
C           (2) INITIAL INCREMENT IS EQUAL TO 0 OR HAS WRONG SIGN       DRKG 800
C               (ERROR MESSAGES IHLF=12 OR IHLF=13),                    DRKG 810
C           (3) THE WHOLE INTEGRATION INTERVAL IS WORKED THROUGH,       DRKG 820
C           (4) SUBROUTINE OUTP HAS CHANGED PRMT(5) TO NON-ZERO.        DRKG 830
C                                                                       DRKG 840
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DRKG 850
C           THE EXTERNAL SUBROUTINES FCT(X,Y,DERY) AND                  DRKG 860
C           OUTP(X,Y,DERY,IHLF,NDIM,PRMT) MUST BE FURNISHED BY THE USER.DRKG 870
C                                                                       DRKG 880
C        METHOD                                                         DRKG 890
C           EVALUATION IS DONE BY MEANS OF FOURTH ORDER RUNGE-KUTTA     DRKG 900
C           FORMULAE IN THE MODIFICATION DUE TO GILL. ACCURACY IS       DRKG 910
C           TESTED COMPARING THE RESULTS OF THE PROCEDURE WITH SINGLE   DRKG 920
C           AND DOUBLE INCREMENT.                                       DRKG 930
C           SUBROUTINE DRKGS AUTOMATICALLY ADJUSTS THE INCREMENT DURING DRKG 940
C           THE WHOLE COMPUTATION BY HALVING OR DOUBLING. IF MORE THAN  DRKG 950
C           10 BISECTIONS OF THE INCREMENT ARE NECESSARY TO GET         DRKG 960
C           SATISFACTORY ACCURACY, THE SUBROUTINE RETURNS WITH          DRKG 970
C           ERROR MESSAGE IHLF=11 INTO MAIN PROGRAM.                    DRKG 980
C           TO GET FULL FLEXIBILITY IN OUTPUT, AN OUTPUT SUBROUTINE     DRKG 990
C           MUST BE FURNISHED BY THE USER.                              DRKG1000
C           FOR REFERENCE, SEE                                          DRKG1010
C           RALSTON/WILF, MATHEMATICAL METHODS FOR DIGITAL COMPUTERS,   DRKG1020
C           WILEY, NEW YORK/LONDON, 1960, PP.110-120.                   DRKG1030
C                                                                       DRKG1040
C     ..................................................................DRKG1050
C                                                                       DRKG1060
      SUBROUTINE DRKGS(PRMT,Y,DERY,NDIM,IHLF,FCT,OUTP,AUX)              DRKG1070
C                                                                       DRKG1080
C                                                                       DRKG1090
      DIMENSION Y(1),DERY(1),AUX(8,1),A(4),B(4),C(4),PRMT(1)            DRKG1100
      DOUBLE PRECISION PRMT,Y,DERY,AUX,A,B,C,X,XEND,H,AJ,BJ,CJ,R1,R2,   DRKG1110
     1DELT                                                              DRKG1120
      DO 1 I=1,NDIM                                                     DRKG1130
    1 AUX(8,I)=.066666666666666667D0*DERY(I)                            DRKG1140
      X=PRMT(1)                                                         DRKG1150
      XEND=PRMT(2)                                                      DRKG1160
      H=PRMT(3)                                                         DRKG1170
      PRMT(5)=0.D0                                                      DRKG1180
      CALL FCT(X,Y,DERY)                                                DRKG1190
C                                                                       DRKG1200
C     ERROR TEST                                                        DRKG1210
      IF(H*(XEND-X))38,37,2                                             DRKG1220
C                                                                       DRKG1230
C     PREPARATIONS FOR RUNGE-KUTTA METHOD                               DRKG1240
    2 A(1)=.5D0                                                         DRKG1250
      A(2)=.29289321881345248D0                                         DRKG1260
      A(3)=1.7071067811865475D0                                         DRKG1270
      A(4)=.16666666666666667D0                                         DRKG1280
      B(1)=2.D0                                                         DRKG1290
      B(2)=1.D0                                                         DRKG1300
      B(3)=1.D0                                                         DRKG1310
      B(4)=2.D0                                                         DRKG1320
      C(1)=.5D0                                                         DRKG1330
      C(2)=.29289321881345248D0                                         DRKG1340
      C(3)=1.7071067811865475D0                                         DRKG1350
      C(4)=.5D0                                                         DRKG1360
C                                                                       DRKG1370
C     PREPARATIONS OF FIRST RUNGE-KUTTA STEP                            DRKG1380
      DO 3 I=1,NDIM                                                     DRKG1390
      AUX(1,I)=Y(I)                                                     DRKG1400
      AUX(2,I)=DERY(I)                                                  DRKG1410
      AUX(3,I)=0.D0                                                     DRKG1420
    3 AUX(6,I)=0.D0                                                     DRKG1430
      IREC=0                                                            DRKG1440
      H=H+H                                                             DRKG1450
      IHLF=-1                                                           DRKG1460
      ISTEP=0                                                           DRKG1470
      IEND=0                                                            DRKG1480
C                                                                       DRKG1490
C                                                                       DRKG1500
C     START OF A RUNGE-KUTTA STEP                                       DRKG1510
    4 IF((X+H-XEND)*H)7,6,5                                             DRKG1520
    5 H=XEND-X                                                          DRKG1530
    6 IEND=1                                                            DRKG1540
C                                                                       DRKG1550
C     RECORDING OF INITIAL VALUES OF THIS STEP                          DRKG1560
    7 CALL OUTP(X,Y,DERY,IREC,NDIM,PRMT)                                DRKG1570
      IF(PRMT(5))40,8,40                                                DRKG1580
    8 ITEST=0                                                           DRKG1590
    9 ISTEP=ISTEP+1                                                     DRKG1600
C                                                                       DRKG1610
C                                                                       DRKG1620
C     START OF INNERMOST RUNGE-KUTTA LOOP                               DRKG1630
      J=1                                                               DRKG1640
   10 AJ=A(J)                                                           DRKG1650
      BJ=B(J)                                                           DRKG1660
      CJ=C(J)                                                           DRKG1670
      DO 11 I=1,NDIM                                                    DRKG1680
      R1=H*DERY(I)                                                      DRKG1690
      R2=AJ*(R1-BJ*AUX(6,I))                                            DRKG1700
      Y(I)=Y(I)+R2                                                      DRKG1710
      R2=R2+R2+R2                                                       DRKG1720
   11 AUX(6,I)=AUX(6,I)+R2-CJ*R1                                        DRKG1730
      IF(J-4)12,15,15                                                   DRKG1740
   12 J=J+1                                                             DRKG1750
      IF(J-3)13,14,13                                                   DRKG1760
   13 X=X+.5D0*H                                                        DRKG1770
   14 CALL FCT(X,Y,DERY)                                                DRKG1780
      GOTO 10                                                           DRKG1790
C     END OF INNERMOST RUNGE-KUTTA LOOP                                 DRKG1800
C                                                                       DRKG1810
C                                                                       DRKG1820
C     TEST OF ACCURACY                                                  DRKG1830
   15 IF(ITEST)16,16,20                                                 DRKG1840
C                                                                       DRKG1850
C     IN CASE ITEST=0 THERE IS NO POSSIBILITY FOR TESTING OF ACCURACY   DRKG1860
   16 DO 17 I=1,NDIM                                                    DRKG1870
   17 AUX(4,I)=Y(I)                                                     DRKG1880
      ITEST=1                                                           DRKG1890
      ISTEP=ISTEP+ISTEP-2                                               DRKG1900
   18 IHLF=IHLF+1                                                       DRKG1910
      X=X-H                                                             DRKG1920
      H=.5D0*H                                                          DRKG1930
      DO 19 I=1,NDIM                                                    DRKG1940
      Y(I)=AUX(1,I)                                                     DRKG1950
      DERY(I)=AUX(2,I)                                                  DRKG1960
   19 AUX(6,I)=AUX(3,I)                                                 DRKG1970
      GOTO 9                                                            DRKG1980
C                                                                       DRKG1990
C     IN CASE ITEST=1 TESTING OF ACCURACY IS POSSIBLE                   DRKG2000
   20 IMOD=ISTEP/2                                                      DRKG2010
      IF(ISTEP-IMOD-IMOD)21,23,21                                       DRKG2020
   21 CALL FCT(X,Y,DERY)                                                DRKG2030
      DO 22 I=1,NDIM                                                    DRKG2040
      AUX(5,I)=Y(I)                                                     DRKG2050
   22 AUX(7,I)=DERY(I)                                                  DRKG2060
      GOTO 9                                                            DRKG2070
C                                                                       DRKG2080
C     COMPUTATION OF TEST VALUE DELT                                    DRKG2090
   23 DELT=0.D0                                                         DRKG2100
      DO 24 I=1,NDIM                                                    DRKG2110
   24 DELT=DELT+AUX(8,I)*DABS(AUX(4,I)-Y(I))                            DRKG2120
      IF(DELT-PRMT(4))28,28,25                                          DRKG2130
C                                                                       DRKG2140
C     ERROR IS TOO GREAT                                                DRKG2150
   25 IF(IHLF-10)26,36,36                                               DRKG2160
   26 DO 27 I=1,NDIM                                                    DRKG2170
   27 AUX(4,I)=AUX(5,I)                                                 DRKG2180
      ISTEP=ISTEP+ISTEP-4                                               DRKG2190
      X=X-H                                                             DRKG2200
      IEND=0                                                            DRKG2210
      GOTO 18                                                           DRKG2220
C                                                                       DRKG2230
C     RESULT VALUES ARE GOOD                                            DRKG2240
   28 CALL FCT(X,Y,DERY)                                                DRKG2250
      DO 29 I=1,NDIM                                                    DRKG2260
      AUX(1,I)=Y(I)                                                     DRKG2270
      AUX(2,I)=DERY(I)                                                  DRKG2280
      AUX(3,I)=AUX(6,I)                                                 DRKG2290
      Y(I)=AUX(5,I)                                                     DRKG2300
   29 DERY(I)=AUX(7,I)                                                  DRKG2310
      CALL OUTP(X-H,Y,DERY,IHLF,NDIM,PRMT)                              DRKG2320
      IF(PRMT(5))40,30,40                                               DRKG2330
   30 DO 31 I=1,NDIM                                                    DRKG2340
      Y(I)=AUX(1,I)                                                     DRKG2350
   31 DERY(I)=AUX(2,I)                                                  DRKG2360
      IREC=IHLF                                                         DRKG2370
      IF(IEND)32,32,39                                                  DRKG2380
C                                                                       DRKG2390
C     INCREMENT GETS DOUBLED                                            DRKG2400
   32 IHLF=IHLF-1                                                       DRKG2410
      ISTEP=ISTEP/2                                                     DRKG2420
      H=H+H                                                             DRKG2430
      IF(IHLF)4,33,33                                                   DRKG2440
   33 IMOD=ISTEP/2                                                      DRKG2450
      IF(ISTEP-IMOD-IMOD)4,34,4                                         DRKG2460
   34 IF(DELT-.02D0*PRMT(4))35,35,4                                     DRKG2470
   35 IHLF=IHLF-1                                                       DRKG2480
      ISTEP=ISTEP/2                                                     DRKG2490
      H=H+H                                                             DRKG2500
      GOTO 4                                                            DRKG2510
C                                                                       DRKG2520
C                                                                       DRKG2530
C     RETURNS TO CALLING PROGRAM                                        DRKG2540
   36 IHLF=11                                                           DRKG2550
      CALL FCT(X,Y,DERY)                                                DRKG2560
      GOTO 39                                                           DRKG2570
   37 IHLF=12                                                           DRKG2580
      GOTO 39                                                           DRKG2590
   38 IHLF=13                                                           DRKG2600
   39 CALL OUTP(X,Y,DERY,IHLF,NDIM,PRMT)                                DRKG2610
   40 RETURN                                                            DRKG2620
      END                                                               DRKG2630
