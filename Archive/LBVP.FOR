C                                                                       LBVP  10
C     ..................................................................LBVP  20
C                                                                       LBVP  30
C        SUBROUTINE LBVP                                                LBVP  40
C                                                                       LBVP  50
C        PURPOSE                                                        LBVP  60
C           TO SOLVE A LINEAR BOUNDARY VALUE PROBLEM, WHICH CONSISTS OF LBVP  70
C           A SYSTEM OF NDIM LINEAR FIRST ORDER DIFFERENTIAL EQUATIONS  LBVP  80
C                  DY/DX=A(X)*Y(X)+F(X)                                 LBVP  90
C           AND NDIM LINEAR BOUNDARY CONDITIONS                         LBVP 100
C                  B*Y(XL)+C*Y(XU)=R.                                   LBVP 110
C                                                                       LBVP 120
C        USAGE                                                          LBVP 130
C           CALL LBVP (PRMT,B,C,R,Y,DERY,NDIM,IHLF,AFCT,FCT,DFCT,OUTP,  LBVP 140
C                      AUX,A)                                           LBVP 150
C           PARAMETERS AFCT,FCT,DFCT,OUTP REQUIRE AN EXTERNAL STATEMENT.LBVP 160
C                                                                       LBVP 170
C        DESCRIPTION OF PARAMETERS                                      LBVP 180
C           PRMT   - AN INPUT AND OUTPUT VECTOR WITH DIMENSION GREATER  LBVP 190
C                    OR EQUAL TO 5, WHICH SPECIFIES THE PARAMETERS OF   LBVP 200
C                    THE INTERVAL AND OF ACCURACY AND WHICH SERVES FOR  LBVP 210
C                    COMMUNICATION BETWEEN OUTPUT SUBROUTINE (FURNISHED LBVP 220
C                    BY THE USER) AND SUBROUTINE LBVP.                  LBVP 230
C                    THE COMPONENTS ARE                                 LBVP 240
C           PRMT(1)- LOWER BOUND XL OF THE INTERVAL (INPUT),            LBVP 250
C           PRMT(1)- UPPER BOUND XU OF THE INTERVAL (INPUT),            LBVP 260
C           PRMT(3)- INITIAL INCREMENT OF THE INDEPENDENT VARIABLE      LBVP 270
C                    (INPUT),                                           LBVP 280
C           PRMT(4)- UPPER ERROR BOUND (INPUT). IF RELATIVE ERROR IS    LBVP 290
C                    GREATER THAN PRMT(4), INCREMENT GETS HALVED.       LBVP 300
C                    IF INCREMENT IS LESS THAN PRMT(3) AND RELATIVE     LBVP 310
C                    ERROR LESS THAN PRMT(4)/50, INCREMENT GETS DOUBLED.LBVP 320
C                    THE USER MAY CHANGE PRMT(4) BY MEANS OF HIS        LBVP 330
C                    OUTPUT SUBROUTINE.                                 LBVP 340
C           PRMT(5)- NO INPUT PARAMETER. SUBROUTINE LBVP INITIALIZES    LBVP 350
C                    PRMT(5)=0. IF THE USER WANTS TO TERMINATE          LBVP 360
C                    SUBROUTINE LBVP AT ANY OUTPUT POINT, HE HAS TO     LBVP 370
C                    CHANGE PRMT(5) TO NON-ZERO BY MEANS OF SUBROUTINE  LBVP 380
C                    OUTP. FURTHER COMPONENTS OF VECTOR PRMT ARE        LBVP 390
C                    FEASIBLE IF ITS DIMENSION IS DEFINED GREATER       LBVP 400
C                    THAN 5. HOWEVER SUBROUTINE LBVP DOES NOT REQUIRE   LBVP 410
C                    AND CHANGE THEM. NEVERTHELESS THEY MAY BE USEFUL   LBVP 420
C                    FOR HANDING RESULT VALUES TO THE MAIN PROGRAM      LBVP 430
C                    (CALLING LBVP) WHICH ARE OBTAINED BY SPECIAL       LBVP 440
C                    MANIPULATIONS WITH OUTPUT DATA IN SUBROUTINE OUTP. LBVP 450
C           B      - AN NDIM BY NDIM INPUT MATRIX.  (DESTROYED)         LBVP 460
C                    IT IS THE COEFFICIENT MATRIX OF Y(XL) IN           LBVP 470
C                    THE BOUNDARY CONDITIONS.                           LBVP 480
C           C      - AN NDIM BY NDIM INPUT MATRIX (POSSIBLY DESTROYED). LBVP 490
C                    IT IS THE COEFFICIENT MATRIX OF Y(XU) IN           LBVP 500
C                    THE BOUNDARY CONDITIONS.                           LBVP 510
C           R      - AN INPUT VECTOR WITH DIMENSION NDIM.  (DESTROYED)  LBVP 520
C                    IT SPECIFIES THE RIGHT HAND SIDE OF THE            LBVP 530
C                    BOUNDARY CONDITIONS.                               LBVP 540
C           Y      - AN AUXILIARY VECTOR WITH DIMENSION NDIM.           LBVP 550
C                    IT IS USED AS STORAGE LOCATION FOR THE RESULTING   LBVP 560
C                    VALUES OF DEPENDENT VARIABLES COMPUTED AT          LBVP 570
C                    INTERMEDIATE POINTS.                               LBVP 580
C           DERY   - INPUT VECTOR OF ERROR WEIGHTS.  (DESTROYED)        LBVP 590
C                    ITS MAXIMAL COMPONENT SHOULD BE EQUAL TO 1.        LBVP 600
C                    LATERON DERY IS THE VECTOR OF DERIVATIVES, WHICH   LBVP 610
C                    BELONG TO FUNCTION VALUES Y AT INTERMEDIATE POINTS.LBVP 620
C           NDIM   - AN INPUT VALUE, WHICH SPECIFIES THE NUMBER OF      LBVP 630
C                    DIFFERENTIAL EQUATIONS IN THE SYSTEM.              LBVP 640
C           IHLF   - AN OUTPUT VALUE, WHICH SPECIFIES THE NUMBER OF     LBVP 650
C                    BISECTIONS OF THE INITIAL INCREMENT. IF IHLF GETS  LBVP 660
C                    GREATER THAN 10, SUBROUTINE LBVP RETURNS WITH      LBVP 670
C                    ERROR MESSAGE IHLF=11 INTO MAIN PROGRAM.           LBVP 680
C                    ERROR MESSAGE IHLF=12 OR IHLF=13 APPEARS IN CASE   LBVP 690
C                    PRMT(3)=0 OR IN CASE SIGN(PRMT(3)).NE.SIGN(PRMT(2)-LBVP 700
C                    PRMT(1)) RESPECTIVELY. FINALLY ERROR MESSAGE       LBVP 710
C                    IHLF=14 INDICATES, THAT THERE IS NO SOLUTION OR    LBVP 720
C                    THAT THERE ARE MORE THAN ONE SOLUTION OF THE       LBVP 730
C                    PROBLEM.                                           LBVP 740
C                    A NEGATIVE VALUE OF IHLF HANDED TO SUBROUTINE OUTP LBVP 750
C                    TOGETHER WITH INITIAL VALUES OF FINALLY GENERATED  LBVP 760
C                    INITIAL VALUE PROBLEM INDICATES, THAT THERE WAS    LBVP 770
C                    POSSIBLE LOSS OF SIGNIFICANCE IN THE SOLUTION OF   LBVP 780
C                    THE SYSTEM OF SIMULTANEOUS LINEAR EQUATIONS FOR    LBVP 790
C                    THESE INITIAL VALUES. THE ABSOLUTE VALUE OF IHLF   LBVP 800
C                    SHOWS, AFTER WHICH ELIMINATION STEP OF GAUSS       LBVP 810
C                    ALGORITHM POSSIBLE LOSS OF SIGNIFICANCE WAS        LBVP 820
C                    DETECTED.                                          LBVP 830
C           AFCT   - THE NAME OF AN EXTERNAL SUBROUTINE USED. IT        LBVP 840
C                    COMPUTES THE COEFFICIENT MATRIX A OF VECTOR Y ON   LBVP 850
C                    THE RIGHT HAND SIDE OF THE SYSTEM OF DIFFERENTIAL  LBVP 860
C                    EQUATIONS FOR A GIVEN X-VALUE. ITS PARAMETER LIST  LBVP 870
C                    MUST BE X,A. SUBROUTINE AFCT SHOULD NOT DESTROY X. LBVP 880
C           FCT    - THE NAME OF AN EXTERNAL SUBROUTINE USED. IT        LBVP 890
C                    COMPUTES VECTOR F (INHOMOGENEOUS PART OF THE       LBVP 900
C                    RIGHT HAND SIDE OF THE SYSTEM OF DIFFERENTIAL      LBVP 910
C                    EQUATIONS) FOR A GIVEN X-VALUE. ITS PARAMETER LIST LBVP 920
C                    MUST BE X,F. SUBROUTINE FCT SHOULD NOT DESTROY X.  LBVP 930
C           DFCT   - THE NAME OF AN EXTERNAL SUBROUTINE USED. IT        LBVP 940
C                    COMPUTES VECTOR DF (DERIVATIVE OF THE INHOMOGENEOUSLBVP 950
C                    PART ON THE RIGHT HAND SIDE OF THE SYSTEM OF       LBVP 960
C                    DIFFERENTIAL EQUATIONS) FOR A GIVEN X-VALUE. ITS   LBVP 970
C                    PARAMETER LIST MUST BE X,DF. SUBROUTINE DFCT       LBVP 980
C                    SHOULD NOT DESTROY X.                              LBVP 990
C           OUTP   - THE NAME OF AN EXTERNAL OUTPUT SUBROUTINE USED.    LBVP1000
C                    ITS PARAMETER LIST MUST BE X,Y,DERY,IHLF,NDIM,PRMT.LBVP1010
C                    NONE OF THESE PARAMETERS (EXCEPT, IF NECESSARY,    LBVP1020
C                    PRMT(4),PRMT(5),...) SHOULD BE CHANGED BY          LBVP1030
C                    SUBROUTINE OUTP. IF PRMT(5) IS CHANGED TO NON-ZERO,LBVP1040
C                    SUBROUTINE LBVP IS TERMINATED.                     LBVP1050
C           AUX    - AN AUXILIARY STORAGE ARRAY WIRH 20 ROWS AND        LBVP1060
C                    NDIM COLUMNS.                                      LBVP1070
C           A      - AN NDIM BY NDIM MATRIX, WHICH IS USED AS AUXILIARY LBVP1080
C                    STORAGE ARRAY.                                     LBVP1090
C                                                                       LBVP1100
C        REMARKS                                                        LBVP1110
C           THE PROCEDURE TERMINATES AND RETURNS TO CALLING PROGRAM, IF LBVP1120
C           (1) MORE THAN 10 BISECTIONS OF THE INITIAL INCREMENT ARE    LBVP1130
C               NECESSARY TO GET SATISFACTORY ACCURACY (ERROR MESSAGE   LBVP1140
C               IHLF=11),                                               LBVP1150
C           (2) INITIAL INCREMENT IS EQUAL TO 0 OR IF IT HAS WRONG SIGN LBVP1160
C               (ERROR MESSAGES IHLF=12 OR IHLF=13),                    LBVP1170
C           (3) THERE IS NO OR MORE THAN ONE SOLUTION OF THE PROBLEM    LBVP1180
C               (ERROR MESSAGE IHLF=14),                                LBVP1190
C           (4) THE WHOLE INTEGRATION INTERVAL IS WORKED THROUGH,       LBVP1200
C           (5) SUBROUTINE OUTP HAS CHANGED PRMT(5) TO NON-ZERO.        LBVP1210
C                                                                       LBVP1220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  LBVP1230
C           SUBROUTINE GELG     SYSTEM OF LINEAR EQUATIONS.             LBVP1240
C           THE EXTERNAL SUBROUTINES AFCT(X,A), FCT(X,F), DFCT(X,DF),   LBVP1250
C           AND OUTP(X,Y,DERY,IHLF,NDIM,PRMT) MUST BE FURNISHED         LBVP1260
C           BY THE USER.                                                LBVP1270
C                                                                       LBVP1280
C        METHOD                                                         LBVP1290
C           EVALUATION IS DONE USING THE METHOD OF ADJOINT EQUATIONS.   LBVP1300
C           HAMMINGS FOURTH ORDER MODIFIED PREDICTOR-CORRECTOR METHOD   LBVP1310
C           IS USED TO SOLVE THE ADJOINT INITIAL VALUE PROBLEMS AND FI- LBVP1320
C           NALLY TO SOLVE THE GENERATED INITIAL VALUE PROBLEM FOR Y(X).LBVP1330
C           THE INITIAL INCREMENT PRMT(3) IS AUTOMATICALLY ADJUSTED.    LBVP1340
C           FOR COMPUTATION OF INTEGRAL SUM, A FOURTH ORDER HERMITEAN   LBVP1350
C           INTEGRATION FORMULA IS USED.                                LBVP1360
C           FOR REFERENCE, SEE                                          LBVP1370
C           (1) LANCE, NUMERICAL METHODS FOR HIGH SPEED COMPUTERS,      LBVP1380
C               ILIFFE, LONDON, 1960, PP.64-67.                         LBVP1390
C           (2) RALSTON/WILF, MATHEMATICAL METHODS FOR DIGITAL          LBVP1400
C               COMPUTERS, WILEY, NEW YORK/LONDON, 1960, PP.95-109.     LBVP1410
C           (3) RALSTON, RUNGE-KUTTA METHODS WITH MINIMUM ERROR BOUNDS, LBVP1420
C               MTAC, VOL.16, ISS.80 (1962), PP.431-437.                LBVP1430
C           (4) ZURMUEHL, PRAKTISCHE MATHEMATIK FUER INGENIEURE UND     LBVP1440
C               PHYSIKER, SPRINGER, BERLIN/GOETTINGEN/HEIDELBERG, 1963, LBVP1450
C               PP.227-232.                                             LBVP1460
C                                                                       LBVP1470
C     ..................................................................LBVP1480
C                                                                       LBVP1490
     0SUBROUTINE LBVP(PRMT,B,C,R,Y,DERY,NDIM,IHLF,AFCT,FCT,DFCT,OUTP,   LBVP1500
     1AUX,A)                                                            LBVP1510
C                                                                       LBVP1520
      DIMENSION PRMT(1),B(1),C(1),R(1),Y(1),DERY(1),AUX(20,1),A(1)      LBVP1530
C                                                                       LBVP1540
C     ERROR TEST                                                        LBVP1550
      IF(PRMT(3)*(PRMT(2)-PRMT(1)))2,1,3                                LBVP1560
    1 IHLF=12                                                           LBVP1570
      RETURN                                                            LBVP1580
    2 IHLF=13                                                           LBVP1590
      RETURN                                                            LBVP1600
C                                                                       LBVP1610
C     SEARCH FOR ZERO-COLUMNS IN MATRICES B AND C                       LBVP1620
    3 KK=-NDIM                                                          LBVP1630
      IB=0                                                              LBVP1640
      IC=0                                                              LBVP1650
      DO 7 K=1,NDIM                                                     LBVP1660
      AUX(15,K)=DERY(K)                                                 LBVP1670
      AUX(1,K)=1.                                                       LBVP1680
      AUX(17,K)=1.                                                      LBVP1690
      KK=KK+NDIM                                                        LBVP1700
      DO 4 I=1,NDIM                                                     LBVP1710
      II=KK+I                                                           LBVP1720
      IF(B(II))5,4,5                                                    LBVP1730
    4 CONTINUE                                                          LBVP1740
      IB=IB+1                                                           LBVP1750
      AUX(1,K)=0.                                                       LBVP1760
    5 DO 6 I=1,NDIM                                                     LBVP1770
      II=KK+I                                                           LBVP1780
      IF(C(II))7,6,7                                                    LBVP1790
    6 CONTINUE                                                          LBVP1800
      IC=IC+1                                                           LBVP1810
      AUX(17,K)=0.                                                      LBVP1820
    7 CONTINUE                                                          LBVP1830
C                                                                       LBVP1840
C     DETERMINATION OF LOWER AND UPPER BOUND                            LBVP1850
      IF(IC-IB)8,11,11                                                  LBVP1860
    8 H=PRMT(2)                                                         LBVP1870
      PRMT(2)=PRMT(1)                                                   LBVP1880
      PRMT(1)=H                                                         LBVP1890
      PRMT(3)=-PRMT(3)                                                  LBVP1900
      DO 9 I=1,NDIM                                                     LBVP1910
    9 AUX(17,I)=AUX(1,I)                                                LBVP1920
      II=NDIM*NDIM                                                      LBVP1930
      DO 10 I=1,II                                                      LBVP1940
      H=B(I)                                                            LBVP1950
      B(I)=C(I)                                                         LBVP1960
   10 C(I)=H                                                            LBVP1970
C                                                                       LBVP1980
C     PREPARATIONS FOR CONSTRUCTION OF ADJOINT INITIAL VALUE PROBLEMS   LBVP1990
   11 X=PRMT(2)                                                         LBVP2000
      CALL FCT(X,Y)                                                     LBVP2010
      CALL DFCT(X,DERY)                                                 LBVP2020
      DO 12 I=1,NDIM                                                    LBVP2030
      AUX(18,I)=Y(I)                                                    LBVP2040
   12 AUX(19,I)=DERY(I)                                                 LBVP2050
C                                                                       LBVP2060
C     POSSIBLE BREAK-POINT FOR LINKAGE                                  LBVP2070
C                                                                       LBVP2080
C     THE FOLLOWING PART OF SUBROUTINE LBVP UNTIL NEXT BREAK-POINT FOR  LBVP2090
C     LINKAGE HAS TO REMAIN IN CORE DURING THE WHOLE REST OF THE        LBVP2100
C     COMPUTATIONS                                                      LBVP2110
C                                                                       LBVP2120
C     START LOOP FOR GENERATING ADJOINT INITIAL VALUE PROBLEMS          LBVP2130
      K=0                                                               LBVP2140
      KK=0                                                              LBVP2150
  100 K=K+1                                                             LBVP2160
      IF(AUX(17,K))108,108,101                                          LBVP2170
C                                                                       LBVP2180
C     INITIALIZATION OF ADJOINT INITIAL VALUE PROBLEM                   LBVP2190
  101 X=PRMT(2)                                                         LBVP2200
      CALL AFCT(X,A)                                                    LBVP2210
      SUM=0.                                                            LBVP2220
      GL=AUX(18,K)                                                      LBVP2230
      DGL=AUX(19,K)                                                     LBVP2240
      II=K                                                              LBVP2250
      DO 104 I=1,NDIM                                                   LBVP2260
      H=-A(II)                                                          LBVP2270
      DERY(I)=H                                                         LBVP2280
      AUX(20,I)=R(I)                                                    LBVP2290
      Y(I)=0.                                                           LBVP2300
      IF(I-K)103,102,103                                                LBVP2310
  102 Y(I)=1.                                                           LBVP2320
  103 DGL=DGL+H*AUX(18,I)                                               LBVP2330
  104 II=II+NDIM                                                        LBVP2340
      XEND=PRMT(1)                                                      LBVP2350
      H=.0625*(XEND-X)                                                  LBVP2360
      ISW=0                                                             LBVP2370
      GOTO 400                                                          LBVP2380
C     THIS IS BRANCH TO ADJOINT LINEAR INITIAL VALUE PROBLEM            LBVP2390
C                                                                       LBVP2400
C     THIS IS RETURN FROM ADJOINT LINEAR INITIAL VALUE PROBLEM          LBVP2410
  105 IF(IHLF-10)106,106,117                                            LBVP2420
C                                                                       LBVP2430
C     UPDATING OF COEFFICIENT MATRIX B AND VECTOR R                     LBVP2440
  106 DO 107 I=1,NDIM                                                   LBVP2450
      KK=KK+1                                                           LBVP2460
      H=C(KK)                                                           LBVP2470
      R(I)=AUX(20,I)+H*SUM                                              LBVP2480
      II=I                                                              LBVP2490
      DO 107 J=1,NDIM                                                   LBVP2500
      B(II)=B(II)+H*Y(J)                                                LBVP2510
  107 II=II+NDIM                                                        LBVP2520
      GOTO 109                                                          LBVP2530
  108 KK=KK+NDIM                                                        LBVP2540
  109 IF(K-NDIM)100,110,110                                             LBVP2550
C                                                                       LBVP2560
C     GENERATION OF LAST INITIAL VALUE PROBLEM                          LBVP2570
  110 X=PRMT(4)                                                         LBVP2580
      CALL GELG(R,B,NDIM,1,X,I)                                         LBVP2590
      IF(I)111,112,112                                                  LBVP2600
  111 IHLF=14                                                           LBVP2610
      RETURN                                                            LBVP2620
C                                                                       LBVP2630
  112 PRMT(5)=0.                                                        LBVP2640
      IHLF=-I                                                           LBVP2650
      X=PRMT(1)                                                         LBVP2660
      XEND=PRMT(2)                                                      LBVP2670
      H=PRMT(3)                                                         LBVP2680
      DO 113 I=1,NDIM                                                   LBVP2690
  113 Y(I)=R(I)                                                         LBVP2700
      ISW=1                                                             LBVP2710
  114 ISW2=12                                                           LBVP2720
      GOTO 200                                                          LBVP2730
  115 ISW3=-1                                                           LBVP2740
      GOTO 300                                                          LBVP2750
  116 IF(IHLF)400,400,117                                               LBVP2760
C     THIS WAS BRANCH INTO INITIAL VALUE PROBLEM                        LBVP2770
C                                                                       LBVP2780
C     THIS IS RETURN FROM INITIAL VALUE PROBLEM                         LBVP2790
  117 RETURN                                                            LBVP2800
C                                                                       LBVP2810
C     THIS PART OF LINEAR BOUNDARY VALUE PROBLEM COMPUTES THE RIGHT     LBVP2820
C     HAND SIDE DERY OF THE SYSTEM OF ADJOINT LINEAR DIFFERENTIAL       LBVP2830
C     EQUATIONS (IN CASE ISW=0) OR OF THE GIVEN SYSTEM (IN CASE ISW=1). LBVP2840
  200 CALL AFCT(X,A)                                                    LBVP2850
      IF(ISW)201,201,205                                                LBVP2860
C                                                                       LBVP2870
C     ADJOINT SYSTEM                                                    LBVP2880
  201 LL=0                                                              LBVP2890
      DO 203 M=1,NDIM                                                   LBVP2900
      HS=0.                                                             LBVP2910
      DO 202 L=1,NDIM                                                   LBVP2920
      LL=LL+1                                                           LBVP2930
  202 HS=HS-A(LL)*Y(L)                                                  LBVP2940
  203 DERY(M)=HS                                                        LBVP2950
  204 GOTO(502,504,506,407,415,418,608,617,632,634,421,115),ISW2        LBVP2960
C                                                                       LBVP2970
C     GIVEN SYSTEM                                                      LBVP2980
  205 CALL FCT(X,DERY)                                                  LBVP2990
      DO 207 M=1,NDIM                                                   LBVP3000
      LL=M-NDIM                                                         LBVP3010
      HS=0.                                                             LBVP3020
      DO 206 L=1,NDIM                                                   LBVP3030
      LL=LL+NDIM                                                        LBVP3040
  206 HS=HS+A(LL)*Y(L)                                                  LBVP3050
  207 DERY(M)=HS+DERY(M)                                                LBVP3060
      GOTO 204                                                          LBVP3070
C                                                                       LBVP3080
C     THIS PART OF LINEAR BOUNDARY VALUE PROBLEM COMPUTES THE VALUE OF  LBVP3090
C     INTEGRAL SUM, WHICH IS A PART OF THE OUTPUT OF ADJOINT INITIAL    LBVP3100
C     VALUE PROBLEM (IN CASE ISW=0) OR RECORDS RESULT VALUES OF THE     LBVP3110
C     FINAL INITIAL VALUE PROBLEM (IN CASE ISW=1).                      LBVP3120
  300 IF(ISW)301,301,305                                                LBVP3130
C                                                                       LBVP3140
C     ADJOINT PROBLEM                                                   LBVP3150
  301 CALL FCT(X,R)                                                     LBVP3160
      GU=0.                                                             LBVP3170
      DGU=0.                                                            LBVP3180
      DO 302 L=1,NDIM                                                   LBVP3190
      GU=GU+Y(L)*R(L)                                                   LBVP3200
  302 DGU=DGU+DERY(L)*R(L)                                              LBVP3210
      CALL DFCT(X,R)                                                    LBVP3220
      DO 303 L=1,NDIM                                                   LBVP3230
  303 DGU=DGU+Y(L)*R(L)                                                 LBVP3240
      SUM=SUM+.5*H*((GL+GU)+.1666667*H*(DGL-DGU))                       LBVP3250
      GL=GU                                                             LBVP3260
      DGL=DGU                                                           LBVP3270
  304 IF(ISW3)116,422,618                                               LBVP3280
C                                                                       LBVP3290
C     GIVEN PROBLEM                                                     LBVP3300
  305 CALL OUTP(X,Y,DERY,IHLF,NDIM,PRMT)                                LBVP3310
      IF(PRMT(5))117,304,117                                            LBVP3320
C                                                                       LBVP3330
C     POSSIBLE BREAK-POINT FOR LINKAGE                                  LBVP3340
C                                                                       LBVP3350
C     THE FOLLOWING PART OF SUBROUTINE LBVP SOLVES IN CASE ISW=0 THE    LBVP3360
C     ADJOINT INITIAL VALUE PROBLEM. IT COMPUTES INTEGRAL SUM AND       LBVP3370
C     THE VECTOR Y OF DEPENDENT VARIABLES AT THE LOWER BOUND PRMT(1).   LBVP3380
C     IN CASE ISW=1 IT SOLVES FINALLY GENERATED INITIAL VALUE PROBLEM.  LBVP3390
  400 N=1                                                               LBVP3400
      XST=X                                                             LBVP3410
      IHLF=0                                                            LBVP3420
      DO 401 I=1,NDIM                                                   LBVP3430
      AUX(16,I)=0.                                                      LBVP3440
      AUX(1,I)=Y(I)                                                     LBVP3450
  401 AUX(8,I)=DERY(I)                                                  LBVP3460
      ISW1=1                                                            LBVP3470
      GOTO 500                                                          LBVP3480
C                                                                       LBVP3490
  402 X=X+H                                                             LBVP3500
      DO 403 I=1,NDIM                                                   LBVP3510
  403 AUX(2,I)=Y(I)                                                     LBVP3520
C                                                                       LBVP3530
C     INCREMENT H IS TESTED BY MEANS OF BISECTION                       LBVP3540
  404 IHLF=IHLF+1                                                       LBVP3550
      X=X-H                                                             LBVP3560
      DO 405 I=1,NDIM                                                   LBVP3570
  405 AUX(4,I)=AUX(2,I)                                                 LBVP3580
      H=.5*H                                                            LBVP3590
      N=1                                                               LBVP3600
      ISW1=2                                                            LBVP3610
      GOTO 500                                                          LBVP3620
C                                                                       LBVP3630
  406 X=X+H                                                             LBVP3640
      ISW2=4                                                            LBVP3650
      GOTO 200                                                          LBVP3660
  407 N=2                                                               LBVP3670
      DO 408 I=1,NDIM                                                   LBVP3680
      AUX(2,I)=Y(I)                                                     LBVP3690
  408 AUX(9,I)=DERY(I)                                                  LBVP3700
      ISW1=3                                                            LBVP3710
      GOTO 500                                                          LBVP3720
C                                                                       LBVP3730
C     TEST ON SATISFACTORY ACCURACY                                     LBVP3740
  409 DO 414 I=1,NDIM                                                   LBVP3750
      Z=ABS(Y(I))                                                       LBVP3760
      IF(Z-1.)410,411,411                                               LBVP3770
  410 Z=1.                                                              LBVP3780
  411 DELT=.06666667*ABS(Y(I)-AUX(4,I))                                 LBVP3790
      IF(ISW)413,413,412                                                LBVP3800
  412 DELT=AUX(15,I)*DELT                                               LBVP3810
  413 IF(DELT-Z*PRMT(4))414,414,429                                     LBVP3820
  414 CONTINUE                                                          LBVP3830
C                                                                       LBVP3840
C     SATISFACTORY ACCURACY AFTER LESS THAN 11 BISECTIONS               LBVP3850
      X=X+H                                                             LBVP3860
      ISW2=5                                                            LBVP3870
      GOTO 200                                                          LBVP3880
  415 DO 416 I=1,NDIM                                                   LBVP3890
      AUX(3,I)=Y(I)                                                     LBVP3900
  416 AUX(10,I)=DERY(I)                                                 LBVP3910
      N=3                                                               LBVP3920
      ISW1=4                                                            LBVP3930
      GOTO 500                                                          LBVP3940
C                                                                       LBVP3950
  417 N=1                                                               LBVP3960
      X=X+H                                                             LBVP3970
      ISW2=6                                                            LBVP3980
      GOTO 200                                                          LBVP3990
  418 X=XST                                                             LBVP4000
      DO 419 I=1,NDIM                                                   LBVP4010
      AUX(11,I)=DERY(I)                                                 LBVP4020
  4190Y(I)=AUX(1,I)+H*(.375*AUX(8,I)+.7916667*AUX(9,I)                  LBVP4030
     1-.2083333*AUX(10,I)+.04166667*DERY(I))                            LBVP4040
  420 X=X+H                                                             LBVP4050
      N=N+1                                                             LBVP4060
      ISW2=11                                                           LBVP4070
      GOTO 200                                                          LBVP4080
  421 ISW3=0                                                            LBVP4090
      GOTO 300                                                          LBVP4100
  422 IF(N-4)423,600,600                                                LBVP4110
  423 DO 424 I=1,NDIM                                                   LBVP4120
      AUX(N,I)=Y(I)                                                     LBVP4130
  424 AUX(N+7,I)=DERY(I)                                                LBVP4140
      IF(N-3)425,427,600                                                LBVP4150
C                                                                       LBVP4160
  425 DO 426 I=1,NDIM                                                   LBVP4170
      DELT=AUX(9,I)+AUX(9,I)                                            LBVP4180
      DELT=DELT+DELT                                                    LBVP4190
  426 Y(I)=AUX(1,I)+.3333333*H*(AUX(8,I)+DELT+AUX(10,I))                LBVP4200
      GOTO 420                                                          LBVP4210
C                                                                       LBVP4220
  427 DO 428 I=1,NDIM                                                   LBVP4230
      DELT=AUX(9,I)+AUX(10,I)                                           LBVP4240
      DELT=DELT+DELT+DELT                                               LBVP4250
  428 Y(I)=AUX(1,I)+.375*H*(AUX(8,I)+DELT+AUX(11,I))                    LBVP4260
      GOTO 420                                                          LBVP4270
C                                                                       LBVP4280
C     NO SATISFACTORY ACCURACY. H MUST BE HALVED.                       LBVP4290
  429 IF(IHLF-10)404,430,430                                            LBVP4300
C                                                                       LBVP4310
C     NO SATISFACTORY ACCURACY AFTER 10 BISECTIONS. ERROR MESSAGE.      LBVP4320
  430 IHLF=11                                                           LBVP4330
      X=X+H                                                             LBVP4340
      IF(ISW)105,105,114                                                LBVP4350
C                                                                       LBVP4360
C     THIS PART OF LINEAR INITIAL VALUE PROBLEM COMPUTES                LBVP4370
C     STARTING VALUES BY MEANS OF RUNGE-KUTTA METHOD.                   LBVP4380
  500 Z=X                                                               LBVP4390
      DO 501 I=1,NDIM                                                   LBVP4400
      X=H*AUX(N+7,I)                                                    LBVP4410
      AUX(5,I)=X                                                        LBVP4420
  501 Y(I)=AUX(N,I)+.4*X                                                LBVP4430
C                                                                       LBVP4440
      X=Z+.4*H                                                          LBVP4450
      ISW2=1                                                            LBVP4460
      GOTO 200                                                          LBVP4470
  502 DO 503 I=1,NDIM                                                   LBVP4480
      X=H*DERY(I)                                                       LBVP4490
      AUX(6,I)=X                                                        LBVP4500
  503 Y(I)=AUX(N,I)+.2969776*AUX(5,I)+.1587596*X                        LBVP4510
C                                                                       LBVP4520
      X=Z+.4557372*H                                                    LBVP4530
      ISW2=2                                                            LBVP4540
      GOTO 200                                                          LBVP4550
  504 DO 505 I=1,NDIM                                                   LBVP4560
      X=H*DERY(I)                                                       LBVP4570
      AUX(7,I)=X                                                        LBVP4580
  505 Y(I)=AUX(N,I)+.2181004*AUX(5,I)-3.050965*AUX(6,I)+3.832865*X      LBVP4590
C                                                                       LBVP4600
      X=Z+H                                                             LBVP4610
      ISW2=3                                                            LBVP4620
      GOTO 200                                                          LBVP4630
  506 DO 507 I=1,NDIM                                                   LBVP4640
  5070Y(I)=AUX(N,I)+.1747603*AUX(5,I)-.5514807*AUX(6,I)                 LBVP4650
     1+1.205536*AUX(7,I)+.1711848*H*DERY(I)                             LBVP4660
      X=Z                                                               LBVP4670
      GOTO(402,406,409,417),ISW1                                        LBVP4680
C                                                                       LBVP4690
C     POSSIBLE BREAK-POINT FOR LINKAGE                                  LBVP4700
C                                                                       LBVP4710
C     STARTING VALUES ARE COMPUTED.                                     LBVP4720
C     NOW START HAMMINGS MODIFIED PREDICTOR-CORRECTOR METHOD.           LBVP4730
  600 ISTEP=3                                                           LBVP4740
  601 IF(N-8)604,602,604                                                LBVP4750
C                                                                       LBVP4760
C     N=8 CAUSES THE ROWS OF AUX TO CHANGE THEIR STORAGE LOCATIONS      LBVP4770
  602 DO 603 N=2,7                                                      LBVP4780
      DO 603 I=1,NDIM                                                   LBVP4790
      AUX(N-1,I)=AUX(N,I)                                               LBVP4800
  603 AUX(N+6,I)=AUX(N+7,I)                                             LBVP4810
      N=7                                                               LBVP4820
C                                                                       LBVP4830
C     N LESS THAN 8 CAUSES N+1 TO GET N                                 LBVP4840
  604 N=N+1                                                             LBVP4850
C                                                                       LBVP4860
C     COMPUTATION OF NEXT VECTOR Y                                      LBVP4870
      DO 605 I=1,NDIM                                                   LBVP4880
      AUX(N-1,I)=Y(I)                                                   LBVP4890
  605 AUX(N+6,I)=DERY(I)                                                LBVP4900
      X=X+H                                                             LBVP4910
  606 ISTEP=ISTEP+1                                                     LBVP4920
      DO 607 I=1,NDIM                                                   LBVP4930
     0DELT=AUX(N-4,I)+1.333333*H*(AUX(N+6,I)+AUX(N+6,I)-AUX(N+5,I)+     LBVP4940
     1AUX(N+4,I)+AUX(N+4,I))                                            LBVP4950
      Y(I)=DELT-.9256198*AUX(16,I)                                      LBVP4960
  607 AUX(16,I)=DELT                                                    LBVP4970
C     PREDICTOR IS NOW GENERATED IN ROW 16 OF AUX, MODIFIED PREDICTOR   LBVP4980
C     IS GENERATED IN Y. DELT MEANS AN AUXILIARY STORAGE.               LBVP4990
C                                                                       LBVP5000
      ISW2=7                                                            LBVP5010
      GOTO 200                                                          LBVP5020
C     DERIVATIVE OF MODIFIED PREDICTOR IS GENERATED IN DERY.            LBVP5030
C                                                                       LBVP5040
  608 DO 609 I=1,NDIM                                                   LBVP5050
     0DELT=.125*(9.*AUX(N-1,I)-AUX(N-3,I)+3.*H*(DERY(I)+AUX(N+6,I)+     LBVP5060
     1AUX(N+6,I)-AUX(N+5,I)))                                           LBVP5070
      AUX(16,I)=AUX(16,I)-DELT                                          LBVP5080
  609 Y(I)=DELT+.07438017*AUX(16,I)                                     LBVP5090
C                                                                       LBVP5100
C     TEST WHETHER H MUST BE HALVED OR DOUBLED                          LBVP5110
      DELT=0.                                                           LBVP5120
      DO 616 I=1,NDIM                                                   LBVP5130
      Z=ABS(Y(I))                                                       LBVP5140
      IF(Z-1.)610,611,611                                               LBVP5150
  610 Z=1.                                                              LBVP5160
  611 Z=ABS(AUX(16,I))/Z                                                LBVP5170
      IF(ISW)613,613,612                                                LBVP5180
  612 Z=AUX(15,I)*Z                                                     LBVP5190
  613 IF(Z-PRMT(4))614,614,628                                          LBVP5200
  614 IF(DELT-Z)615,616,616                                             LBVP5210
  615 DELT=Z                                                            LBVP5220
  616 CONTINUE                                                          LBVP5230
C                                                                       LBVP5240
C     H MUST NOT BE HALVED. THAT MEANS Y(I) ARE GOOD.                   LBVP5250
      ISW2=8                                                            LBVP5260
      GOTO 200                                                          LBVP5270
  617 ISW3=1                                                            LBVP5280
      GOTO 300                                                          LBVP5290
  618 IF(H*(X-XEND))619,621,621                                         LBVP5300
  619 IF(ABS(X-XEND)-.1*ABS(H))621,620,620                              LBVP5310
  620 IF(DELT-.02*PRMT(4))622,622,601                                   LBVP5320
  621 IF(ISW)105,105,117                                                LBVP5330
C                                                                       LBVP5340
C                                                                       LBVP5350
C     H COULD BE DOUBLED IF ALL NECESSARY PRECEEDING VALUES ARE         LBVP5360
C     AVAILABLE.                                                        LBVP5370
  622 IF(IHLF)601,601,623                                               LBVP5380
  623 IF(N-7)601,624,624                                                LBVP5390
  624 IF(ISTEP-4)601,625,625                                            LBVP5400
  625 IMOD=ISTEP/2                                                      LBVP5410
      IF(ISTEP-IMOD-IMOD)601,626,601                                    LBVP5420
  626 H=H+H                                                             LBVP5430
      IHLF=IHLF-1                                                       LBVP5440
      ISTEP=0                                                           LBVP5450
      DO 627 I=1,NDIM                                                   LBVP5460
      AUX(N-1,I)=AUX(N-2,I)                                             LBVP5470
      AUX(N-2,I)=AUX(N-4,I)                                             LBVP5480
      AUX(N-3,I)=AUX(N-6,I)                                             LBVP5490
      AUX(N+6,I)=AUX(N+5,I)                                             LBVP5500
      AUX(N+5,I)=AUX(N+3,I)                                             LBVP5510
      AUX(N+4,I)=AUX(N+1,I)                                             LBVP5520
      DELT=AUX(N+6,I)+AUX(N+5,I)                                        LBVP5530
      DELT=DELT+DELT+DELT                                               LBVP5540
  6270AUX(16,I)=8.962963*(Y(I)-AUX(N-3,I))-3.361111*H*(DERY(I)+DELT     LBVP5550
     1+AUX(N+4,I))                                                      LBVP5560
      GOTO 601                                                          LBVP5570
C                                                                       LBVP5580
C                                                                       LBVP5590
C     H MUST BE HALVED                                                  LBVP5600
  628 IHLF=IHLF+1                                                       LBVP5610
      IF(IHLF-10)630,630,629                                            LBVP5620
  629 IF(ISW)105,105,114                                                LBVP5630
  630 H=.5*H                                                            LBVP5640
      ISTEP=0                                                           LBVP5650
      DO 631 I=1,NDIM                                                   LBVP5660
     0Y(I)=.00390625*(80.*AUX(N-1,I)+135.*AUX(N-2,I)+40.*AUX(N-3,I)+    LBVP5670
     1AUX(N-4,I))-.1171875*(AUX(N+6,I)-6.*AUX(N+5,I)-AUX(N+4,I))*H      LBVP5680
     0AUX(N-4,I)=.00390625*(12.*AUX(N-1,I)+135.*AUX(N-2,I)+             LBVP5690
     1108.*AUX(N-3,I)+AUX(N-4,I))-.0234375*(AUX(N+6,I)+18.*AUX(N+5,I)-  LBVP5700
     29.*AUX(N+4,I))*H                                                  LBVP5710
      AUX(N-3,I)=AUX(N-2,I)                                             LBVP5720
  631 AUX(N+4,I)=AUX(N+5,I)                                             LBVP5730
      DELT=X-H                                                          LBVP5740
      X=DELT-(H+H)                                                      LBVP5750
      ISW2=9                                                            LBVP5760
      GOTO 200                                                          LBVP5770
  632 DO 633 I=1,NDIM                                                   LBVP5780
      AUX(N-2,I)=Y(I)                                                   LBVP5790
      AUX(N+5,I)=DERY(I)                                                LBVP5800
  633 Y(I)=AUX(N-4,I)                                                   LBVP5810
      X=X-(H+H)                                                         LBVP5820
      ISW2=10                                                           LBVP5830
      GOTO 200                                                          LBVP5840
  634 X=DELT                                                            LBVP5850
      DO 635 I=1,NDIM                                                   LBVP5860
      DELT=AUX(N+5,I)+AUX(N+4,I)                                        LBVP5870
      DELT=DELT+DELT+DELT                                               LBVP5880
     0AUX(16,I)=8.962963*(AUX(N-1,I)-Y(I))-3.361111*H*(AUX(N+6,I)+DELT  LBVP5890
     1+DERY(I))                                                         LBVP5900
  635 AUX(N+3,I)=DERY(I)                                                LBVP5910
      GOTO 606                                                          LBVP5920
C                                                                       LBVP5930
C     END OF INITIAL VALUE PROBLEM                                      LBVP5940
      END                                                               LBVP5950
