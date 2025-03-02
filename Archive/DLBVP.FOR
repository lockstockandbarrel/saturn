C                                                                       DLBV  10
C     ..................................................................DLBV  20
C                                                                       DLBV  30
C        SUBROUTINE DLBVP                                               DLBV  40
C                                                                       DLBV  50
C        PURPOSE                                                        DLBV  60
C           TO SOLVE A LINEAR BOUNDARY VALUE PROBLEM, WHICH CONSISTS OF DLBV  70
C           A SYSTEM OF NDIM LINEAR FIRST ORDER DIFFERENTIAL EQUATIONS  DLBV  80
C                  DY/DX=A(X)*Y(X)+F(X)                                 DLBV  90
C           AND NDIM LINEAR BOUNDARY CONDITIONS                         DLBV 100
C                  B*Y(XL)+C*Y(XU)=R.                                   DLBV 110
C                                                                       DLBV 120
C        USAGE                                                          DLBV 130
C           CALL DLBVP (PRMT,B,C,R,Y,DERY,NDIM,IHLF,AFCT,FCT,DFCT,OUTP, DLBV 140
C                      AUX,A)                                           DLBV 150
C           PARAMETERS AFCT,FCT,DFCT,OUTP REQUIRE AN EXTERNAL STATEMENT.DLBV 160
C                                                                       DLBV 170
C        DESCRIPTION OF PARAMETERS                                      DLBV 180
C           PRMT   - DOUBLE PRECISION INPUT AND OUTPUT VECTOR WITH      DLBV 190
C                    DIMENSION GREATER THAN OR EQUAL TO 5, WHICH        DLBV 200
C                    SPECIFIES THE PARAMETERS OF THE INTERVAL AND OF    DLBV 210
C                    ACCURACY AND WHICH SERVES FOR COMMUNICATION BETWEENDLBV 220
C                    OUTPUT SUBROUTINE (FURNISHED BY THE USER) AND      DLBV 230
C                    SUBROUTINE DLBVP. EXCEPT PRMT(5) THE COMPONENTS    DLBV 240
C                    ARE NOT DESTROYED BY SUBROUTINE DLBVP AND THEY ARE DLBV 250
C           PRMT(1)- LOWER BOUND XL OF THE INTERVAL (INPUT),            DLBV 260
C           PRMT(1)- UPPER BOUND XU OF THE INTERVAL (INPUT),            DLBV 270
C           PRMT(3)- INITIAL INCREMENT OF THE INDEPENDENT VARIABLE      DLBV 280
C                    (INPUT),                                           DLBV 290
C           PRMT(4)- UPPER ERROR BOUND (INPUT). IF RELATIVE ERROR IS    DLBV 300
C                    GREATER THAN PRMT(4), INCREMENT GETS HALVED.       DLBV 310
C                    IF INCREMENT IS LESS THAN PRMT(3) AND RELATIVE     DLBV 320
C                    ERROR LESS THAN PRMT(4)/50, INCREMENT GETS DOUBLED.DLBV 330
C                    THE USER MAY CHANGE PRMT(4) BY MEANS OF HIS        DLBV 340
C                    OUTPUT SUBROUTINE.                                 DLBV 350
C           PRMT(5)- NO INPUT PARAMETER. SUBROUTINE DLBVP INITIALIZES   DLBV 360
C                    PRMT(5)=0. IF THE USER WANTS TO TERMINATE          DLBV 370
C                    SUBROUTINE DLBVP AT ANY OUTPUT POINT, HE HAS TO    DLBV 380
C                    CHANGE PRMT(5) TO NON-ZERO BY MEANS OF SUBROUTINE  DLBV 390
C                    OUTP. FURTHER COMPONENTS OF VECTOR PRMT ARE        DLBV 400
C                    FEASIBLE IF ITS DIMENSION IS DEFINED GREATER       DLBV 410
C                    THAN 5. HOWEVER SUBROUTINE DLBVP DOES NOT REQUIRE  DLBV 420
C                    AND CHANGE THEM. NEVERTHELESS THEY MAY BE USEFUL   DLBV 430
C                    FOR HANDING RESULT VALUES TO THE MAIN PROGRAM      DLBV 440
C                    (CALLING DLBVP) WHICH ARE OBTAINED BY SPECIAL      DLBV 450
C                    MANIPULATIONS WITH OUTPUT DATA IN SUBROUTINE OUTP. DLBV 460
C           B      - DOUBLE PRECISION NDIM BY NDIM INPUT MATRIX         DLBV 470
C                    (DESTROYED). IT IS THE COEFFICIENT MATRIX OF Y(XL) DLBV 480
C                    IN THE BOUNDARY CONDITIONS.                        DLBV 490
C           C      - DOUBLE PRECISION NDIM BY NDIM INPUT MATRIX         DLBV 500
C                    (POSSIBLY DESTROYED). IT IS THE COEFFICIENT MATRIX DLBV 510
C                    OF Y(XU) IN THE BOUNDARY CONDITIONS.               DLBV 520
C           R      - DOUBLE PRECISION INPUT VECTOR WITH DIMENSION NDIM  DLBV 530
C                    (DESTROYED). IT SPECIFIES THE RIGHT HAND SIDE OF   DLBV 540
C                    THE BOUNDARY CONDITIONS.                           DLBV 550
C           Y      - DOUBLE PRECISION AUXILIARY VECTOR WITH             DLBV 560
C                    DIMENSION NDIM. IT IS USED AS STORAGE LOCATION     DLBV 570
C                    FOR THE RESULTING VALUES OF DEPENDENT VARIABLES    DLBV 580
C                    COMPUTED AT INTERMEDIATE POINTS X.                 DLBV 590
C           DERY   - DOUBLE PRECISION INPUT VECTOR OF ERROR WEIGHTS     DLBV 600
C                    (DESTROYED). ITS MAXIMAL COMPONENT SHOULD BE       DLBV 610
C                    EQUAL TO 1. LATERON DERY IS THE VECTOR OF          DLBV 620
C                    DERIVATIVES, WHICH BELONG TO FUNCTION VALUES Y AT  DLBV 630
C                    INTERMEDIATE POINTS X.                             DLBV 640
C           NDIM   - AN INPUT VALUE, WHICH SPECIFIES THE NUMBER OF      DLBV 650
C                    DIFFERENTIAL EQUATIONS IN THE SYSTEM.              DLBV 660
C           IHLF   - AN OUTPUT VALUE, WHICH SPECIFIES THE NUMBER OF     DLBV 670
C                    BISECTIONS OF THE INITIAL INCREMENT. IF IHLF GETS  DLBV 680
C                    GREATER THAN 10, SUBROUTINE DLBVP RETURNS WITH     DLBV 690
C                    ERROR MESSAGE IHLF=11 INTO MAIN PROGRAM.           DLBV 700
C                    ERROR MESSAGE IHLF=12 OR IHLF=13 APPEARS IN CASE   DLBV 710
C                    PRMT(3)=0 OR IN CASE SIGN(PRMT(3)).NE.SIGN(PRMT(2)-DLBV 720
C                    PRMT(1)) RESPECTIVELY. FINALLY ERROR MESSAGE       DLBV 730
C                    IHLF=14 INDICATES, THAT THERE IS NO SOLUTION OR    DLBV 740
C                    THAT THERE ARE MORE THAN ONE SOLUTION OF THE       DLBV 750
C                    PROBLEM.                                           DLBV 760
C                    A NEGATIVE VALUE OF IHLF HANDED TO SUBROUTINE OUTP DLBV 770
C                    TOGETHER WITH INITIAL VALUES OF FINALLY GENERATED  DLBV 780
C                    INITIAL VALUE PROBLEM INDICATES, THAT THERE WAS    DLBV 790
C                    POSSIBLE LOSS OF SIGNIFICANCE IN THE SOLUTION OF   DLBV 800
C                    THE SYSTEM OF SIMULTANEOUS LINEAR EQUATIONS FOR    DLBV 810
C                    THESE INITIAL VALUES. THE ABSOLUTE VALUE OF IHLF   DLBV 820
C                    SHOWS, AFTER WHICH ELIMINATION STEP OF GAUSS       DLBV 830
C                    ALGORITHM POSSIBLE LOSS OF SIGNIFICANCE WAS        DLBV 840
C                    DETECTED.                                          DLBV 850
C           AFCT   - THE NAME OF AN EXTERNAL SUBROUTINE USED. IT        DLBV 860
C                    COMPUTES THE COEFFICIENT MATRIX A OF VECTOR Y ON   DLBV 870
C                    THE RIGHT HAND SIDE OF THE SYSTEM OF DIFFERENTIAL  DLBV 880
C                    EQUATIONS FOR A GIVEN X-VALUE. ITS PARAMETER LIST  DLBV 890
C                    MUST BE X,A. SUBROUTINE AFCT SHOULD NOT DESTROY X. DLBV 900
C           FCT    - THE NAME OF AN EXTERNAL SUBROUTINE USED. IT        DLBV 910
C                    COMPUTES VECTOR F (INHOMOGENEOUS PART OF THE       DLBV 920
C                    RIGHT HAND SIDE OF THE SYSTEM OF DIFFERENTIAL      DLBV 930
C                    EQUATIONS) FOR A GIVEN X-VALUE. ITS PARAMETER LIST DLBV 940
C                    MUST BE X,F. SUBROUTINE FCT SHOULD NOT DESTROY X.  DLBV 950
C           DFCT   - THE NAME OF AN EXTERNAL SUBROUTINE USED. IT        DLBV 960
C                    COMPUTES VECTOR DF (DERIVATIVE OF THE INHOMOGENEOUSDLBV 970
C                    PART ON THE RIGHT HAND SIDE OF THE SYSTEM OF       DLBV 980
C                    DIFFERENTIAL EQUATIONS) FOR A GIVEN X-VALUE. ITS   DLBV 990
C                    PARAMETER LIST MUST BE X,DF. SUBROUTINE DFCT       DLBV1000
C                    SHOULD NOT DESTROY X.                              DLBV1010
C           OUTP   - THE NAME OF AN EXTERNAL OUTPUT SUBROUTINE USED.    DLBV1020
C                    ITS PARAMETER LIST MUST BE X,Y,DERY,IHLF,NDIM,PRMT.DLBV1030
C                    NONE OF THESE PARAMETERS (EXCEPT, IF NECESSARY,    DLBV1040
C                    PRMT(4),PRMT(5),...) SHOULD BE CHANGED BY          DLBV1050
C                    SUBROUTINE OUTP. IF PRMT(5) IS CHANGED TO NON-ZERO,DLBV1060
C                    SUBROUTINE DLBVP IS TERMINATED.                    DLBV1070
C           AUX    - DOUBLE PRECISION AUXILIARY STORAGE ARRAY WITH 20   DLBV1080
C                    ROWS AND NDIM COLUMNS.                             DLBV1090
C           A      - DOUBLE PRECISION NDIM BY NDIM MATRIX, WHICH IS USEDDLBV1100
C                    AS AUXILIARY STORAGE ARRAY.                        DLBV1110
C                                                                       DLBV1120
C        REMARKS                                                        DLBV1130
C           THE PROCEDURE TERMINATES AND RETURNS TO CALLING PROGRAM, IF DLBV1140
C           (1) MORE THAN 10 BISECTIONS OF THE INITIAL INCREMENT ARE    DLBV1150
C               NECESSARY TO GET SATISFACTORY ACCURACY (ERROR MESSAGE   DLBV1160
C               IHLF=11),                                               DLBV1170
C           (2) INITIAL INCREMENT IS EQUAL TO 0 OR IF IT HAS WRONG SIGN DLBV1180
C               (ERROR MESSAGES IHLF=12 OR IHLF=13),                    DLBV1190
C           (3) THERE IS NO OR MORE THAN ONE SOLUTION OF THE PROBLEM    DLBV1200
C               (ERROR MESSAGE IHLF=14),                                DLBV1210
C           (4) THE WHOLE INTEGRATION INTERVAL IS WORKED THROUGH,       DLBV1220
C           (5) SUBROUTINE OUTP HAS CHANGED PRMT(5) TO NON-ZERO.        DLBV1230
C                                                                       DLBV1240
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DLBV1250
C           SUBROUTINE DGELG     SYSTEM OF LINEAR EQUATIONS.            DLBV1260
C           THE EXTERNAL SUBROUTINES AFCT(X,A), FCT(X,F), DFCT(X,DF),   DLBV1270
C           AND OUTP(X,Y,DERY,IHLF,NDIM,PRMT) MUST BE FURNISHED         DLBV1280
C           BY THE USER.                                                DLBV1290
C                                                                       DLBV1300
C        METHOD                                                         DLBV1310
C           EVALUATION IS DONE USING THE METHOD OF ADJOINT EQUATIONS.   DLBV1320
C           HAMMINGS FOURTH ORDER MODIFIED PREDICTOR-CORRECTOR METHOD   DLBV1330
C           IS USED TO SOLVE THE ADJOINT INITIAL VALUE PROBLEMS AND FI- DLBV1340
C           NALLY TO SOLVE THE GENERATED INITIAL VALUE PROBLEM FOR Y(X).DLBV1350
C           THE INITIAL INCREMENT PRMT(3) IS AUTOMATICALLY ADJUSTED.    DLBV1360
C           FOR COMPUTATION OF INTEGRAL SUM, A FOURTH ORDER HERMITEAN   DLBV1370
C           INTEGRATION FORMULA IS USED.                                DLBV1380
C           FOR REFERENCE, SEE                                          DLBV1390
C           (1) LANCE, NUMERICAL METHODS FOR HIGH SPEED COMPUTERS,      DLBV1400
C               ILIFFE, LONDON, 1960, PP.64-67.                         DLBV1410
C           (2) RALSTON/WILF, MATHEMATICAL METHODS FOR DIGITAL          DLBV1420
C               COMPUTERS, WILEY, NEW YORK/LONDON, 1960, PP.95-109.     DLBV1430
C           (3) RALSTON, RUNGE-KUTTA METHODS WITH MINIMUM ERROR BOUNDS, DLBV1440
C               MTAC, VOL.16, ISS.80 (1962), PP.431-437.                DLBV1450
C           (4) ZURMUEHL, PRAKTISCHE MATHEMATIK FUER INGENIEURE UND     DLBV1460
C               PHYSIKER, SPRINGER, BERLIN/GOETTINGEN/HEIDELBERG, 1963, DLBV1470
C               PP.227-232.                                             DLBV1480
C                                                                       DLBV1490
C     ..................................................................DLBV1500
C                                                                       DLBV1510
     0SUBROUTINE DLBVP(PRMT,B,C,R,Y,DERY,NDIM,IHLF,AFCT,FCT,DFCT,OUTP,  DLBV1520
     1AUX,A)                                                            DLBV1530
C                                                                       DLBV1540
C                                                                       DLBV1550
      DIMENSION PRMT(1),B(1),C(1),R(1),Y(1),DERY(1),AUX(20,1),A(1)      DLBV1560
      DOUBLE PRECISION PRMT,B,C,R,Y,DERY,AUX,A,H,X,Z,GL,HS,GU,SUM,      DLBV1570
     1DGL,DGU,XST,XEND,DELT                                             DLBV1580
C                                                                       DLBV1590
C     ERROR TEST                                                        DLBV1600
      IF(PRMT(3)*(PRMT(2)-PRMT(1)))2,1,3                                DLBV1610
    1 IHLF=12                                                           DLBV1620
      RETURN                                                            DLBV1630
    2 IHLF=13                                                           DLBV1640
      RETURN                                                            DLBV1650
C                                                                       DLBV1660
C     SEARCH FOR ZERO-COLUMNS IN MATRICES B AND C                       DLBV1670
    3 KK=-NDIM                                                          DLBV1680
      IB=0                                                              DLBV1690
      IC=0                                                              DLBV1700
      DO 7 K=1,NDIM                                                     DLBV1710
      AUX(15,K)=DERY(K)                                                 DLBV1720
      AUX(1,K)=1.D0                                                     DLBV1730
      AUX(17,K)=1.D0                                                    DLBV1740
      KK=KK+NDIM                                                        DLBV1750
      DO 4 I=1,NDIM                                                     DLBV1760
      II=KK+I                                                           DLBV1770
      IF(B(II))5,4,5                                                    DLBV1780
    4 CONTINUE                                                          DLBV1790
      IB=IB+1                                                           DLBV1800
      AUX(1,K)=0.D0                                                     DLBV1810
    5 DO 6 I=1,NDIM                                                     DLBV1820
      II=KK+I                                                           DLBV1830
      IF(C(II))7,6,7                                                    DLBV1840
    6 CONTINUE                                                          DLBV1850
      IC=IC+1                                                           DLBV1860
      AUX(17,K)=0.D0                                                    DLBV1870
    7 CONTINUE                                                          DLBV1880
C                                                                       DLBV1890
C     DETERMINATION OF LOWER AND UPPER BOUND                            DLBV1900
      IF(IC-IB)8,11,11                                                  DLBV1910
    8 H=PRMT(2)                                                         DLBV1920
      PRMT(2)=PRMT(1)                                                   DLBV1930
      PRMT(1)=H                                                         DLBV1940
      PRMT(3)=-PRMT(3)                                                  DLBV1950
      DO 9 I=1,NDIM                                                     DLBV1960
    9 AUX(17,I)=AUX(1,I)                                                DLBV1970
      II=NDIM*NDIM                                                      DLBV1980
      DO 10 I=1,II                                                      DLBV1990
      H=B(I)                                                            DLBV2000
      B(I)=C(I)                                                         DLBV2010
   10 C(I)=H                                                            DLBV2020
C                                                                       DLBV2030
C     PREPARATIONS FOR CONSTRUCTION OF ADJOINT INITIAL VALUE PROBLEMS   DLBV2040
   11 X=PRMT(2)                                                         DLBV2050
      CALL FCT(X,Y)                                                     DLBV2060
      CALL DFCT(X,DERY)                                                 DLBV2070
      DO 12 I=1,NDIM                                                    DLBV2080
      AUX(18,I)=Y(I)                                                    DLBV2090
   12 AUX(19,I)=DERY(I)                                                 DLBV2100
C                                                                       DLBV2110
C     POSSIBLE BREAK-POINT FOR LINKAGE                                  DLBV2120
C                                                                       DLBV2130
C     THE FOLLOWING PART OF SUBROUTINE DLBVP UNTIL NEXT BREAK-POINT FOR DLBV2140
C     LINKAGE HAS TO REMAIN IN CORE DURING THE WHOLE REST OF THE        DLBV2150
C     COMPUTATIONS                                                      DLBV2160
C                                                                       DLBV2170
C     START LOOP FOR GENERATING ADJOINT INITIAL VALUE PROBLEMS          DLBV2180
      K=0                                                               DLBV2190
      KK=0                                                              DLBV2200
  100 K=K+1                                                             DLBV2210
      IF(AUX(17,K))108,108,101                                          DLBV2220
C                                                                       DLBV2230
C     INITIALIZATION OF ADJOINT INITIAL VALUE PROBLEM                   DLBV2240
  101 X=PRMT(2)                                                         DLBV2250
      CALL AFCT(X,A)                                                    DLBV2260
      SUM=0.D0                                                          DLBV2270
      GL=AUX(18,K)                                                      DLBV2280
      DGL=AUX(19,K)                                                     DLBV2290
      II=K                                                              DLBV2300
      DO 104 I=1,NDIM                                                   DLBV2310
      H=-A(II)                                                          DLBV2320
      DERY(I)=H                                                         DLBV2330
      AUX(20,I)=R(I)                                                    DLBV2340
      Y(I)=0.D0                                                         DLBV2350
      IF(I-K)103,102,103                                                DLBV2360
  102 Y(I)=1.D0                                                         DLBV2370
  103 DGL=DGL+H*AUX(18,I)                                               DLBV2380
  104 II=II+NDIM                                                        DLBV2390
      XEND=PRMT(1)                                                      DLBV2400
      H=.0625D0*(XEND-X)                                                DLBV2410
      ISW=0                                                             DLBV2420
      GOTO 400                                                          DLBV2430
C     THIS IS BRANCH TO ADJOINT LINEAR INITIAL VALUE PROBLEM            DLBV2440
C                                                                       DLBV2450
C     THIS IS RETURN FROM ADJOINT LINEAR INITIAL VALUE PROBLEM          DLBV2460
  105 IF(IHLF-10)106,106,117                                            DLBV2470
C                                                                       DLBV2480
C     UPDATING OF COEFFICIENT MATRIX B AND VECTOR R                     DLBV2490
  106 DO 107 I=1,NDIM                                                   DLBV2500
      KK=KK+1                                                           DLBV2510
      H=C(KK)                                                           DLBV2520
      R(I)=AUX(20,I)+H*SUM                                              DLBV2530
      II=I                                                              DLBV2540
      DO 107 J=1,NDIM                                                   DLBV2550
      B(II)=B(II)+H*Y(J)                                                DLBV2560
  107 II=II+NDIM                                                        DLBV2570
      GOTO 109                                                          DLBV2580
  108 KK=KK+NDIM                                                        DLBV2590
  109 IF(K-NDIM)100,110,110                                             DLBV2600
C                                                                       DLBV2610
C                                                                       DLBV2620
C     GENERATION OF LAST INITIAL VALUE PROBLEM                          DLBV2630
  110 EPS=PRMT(4)                                                       DLBV2640
      CALL DGELG(R,B,NDIM,1,EPS,I)                                      DLBV2650
      IF(I)111,112,112                                                  DLBV2660
  111 IHLF=14                                                           DLBV2670
      RETURN                                                            DLBV2680
C                                                                       DLBV2690
  112 PRMT(5)=0.D0                                                      DLBV2700
      IHLF=-I                                                           DLBV2710
      X=PRMT(1)                                                         DLBV2720
      XEND=PRMT(2)                                                      DLBV2730
      H=PRMT(3)                                                         DLBV2740
      DO 113 I=1,NDIM                                                   DLBV2750
  113 Y(I)=R(I)                                                         DLBV2760
      ISW=1                                                             DLBV2770
  114 ISW2=12                                                           DLBV2780
      GOTO 200                                                          DLBV2790
  115 ISW3=-1                                                           DLBV2800
      GOTO 300                                                          DLBV2810
  116 IF(IHLF)400,400,117                                               DLBV2820
C     THIS WAS BRANCH INTO INITIAL VALUE PROBLEM                        DLBV2830
C                                                                       DLBV2840
C     THIS IS RETURN FROM INITIAL VALUE PROBLEM                         DLBV2850
  117 RETURN                                                            DLBV2860
C                                                                       DLBV2870
C     THIS PART OF LINEAR BOUNDARY VALUE PROBLEM COMPUTES THE RIGHT     DLBV2880
C     HAND SIDE DERY OF THE SYSTEM OF ADJOINT LINEAR DIFFERENTIAL       DLBV2890
C     EQUATIONS (IN CASE ISW=0) OR OF THE GIVEN SYSTEM (IN CASE ISW=1). DLBV2900
  200 CALL AFCT(X,A)                                                    DLBV2910
      IF(ISW)201,201,205                                                DLBV2920
C                                                                       DLBV2930
C     ADJOINT SYSTEM                                                    DLBV2940
  201 LL=0                                                              DLBV2950
      DO 203 M=1,NDIM                                                   DLBV2960
      HS=0.D0                                                           DLBV2970
      DO 202 L=1,NDIM                                                   DLBV2980
      LL=LL+1                                                           DLBV2990
  202 HS=HS-A(LL)*Y(L)                                                  DLBV3000
  203 DERY(M)=HS                                                        DLBV3010
  204 GOTO(502,504,506,407,415,418,608,617,632,634,421,115),ISW2        DLBV3020
C                                                                       DLBV3030
C     GIVEN SYSTEM                                                      DLBV3040
  205 CALL FCT(X,DERY)                                                  DLBV3050
      DO 207 M=1,NDIM                                                   DLBV3060
      LL=M-NDIM                                                         DLBV3070
      HS=0.D0                                                           DLBV3080
      DO 206 L=1,NDIM                                                   DLBV3090
      LL=LL+NDIM                                                        DLBV3100
  206 HS=HS+A(LL)*Y(L)                                                  DLBV3110
  207 DERY(M)=HS+DERY(M)                                                DLBV3120
      GOTO 204                                                          DLBV3130
C                                                                       DLBV3140
C     THIS PART OF LINEAR BOUNDARY VALUE PROBLEM COMPUTES THE VALUE OF  DLBV3150
C     INTEGRAL SUM, WHICH IS A PART OF THE OUTPUT OF ADJOINT INITIAL    DLBV3160
C     VALUE PROBLEM (IN CASE ISW=0) OR RECORDS RESULT VALUES OF THE     DLBV3170
C     FINAL INITIAL VALUE PROBLEM (IN CASE ISW=1).                      DLBV3180
  300 IF(ISW)301,301,305                                                DLBV3190
C                                                                       DLBV3200
C     ADJOINT PROBLEM                                                   DLBV3210
  301 CALL FCT(X,R)                                                     DLBV3220
      GU=0.D0                                                           DLBV3230
      DGU=0.D0                                                          DLBV3240
      DO 302 L=1,NDIM                                                   DLBV3250
      GU=GU+Y(L)*R(L)                                                   DLBV3260
  302 DGU=DGU+DERY(L)*R(L)                                              DLBV3270
      CALL DFCT(X,R)                                                    DLBV3280
      DO 303 L=1,NDIM                                                   DLBV3290
  303 DGU=DGU+Y(L)*R(L)                                                 DLBV3300
      SUM=SUM+.5D0*H*((GL+GU)+.16666666666666667D0*H*(DGL-DGU))         DLBV3310
      GL=GU                                                             DLBV3320
      DGL=DGU                                                           DLBV3330
  304 IF(ISW3)116,422,618                                               DLBV3340
C                                                                       DLBV3350
C     GIVEN PROBLEM                                                     DLBV3360
  305 CALL OUTP(X,Y,DERY,IHLF,NDIM,PRMT)                                DLBV3370
      IF(PRMT(5))117,304,117                                            DLBV3380
C                                                                       DLBV3390
C     POSSIBLE BREAK-POINT FOR LINKAGE                                  DLBV3400
C                                                                       DLBV3410
C     THE FOLLOWING PART OF SUBROUTINE DLBVP SOLVES IN CASE ISW=0 THE   DLBV3420
C     ADJOINT INITIAL VALUE PROBLEM. IT COMPUTES INTEGRAL SUM AND       DLBV3430
C     THE VECTOR Y OF DEPENDENT VARIABLES AT THE LOWER BOUND PRMT(1).   DLBV3440
C     IN CASE ISW=1 IT SOLVES FINALLY GENERATED INITIAL VALUE PROBLEM.  DLBV3450
  400 N=1                                                               DLBV3460
      XST=X                                                             DLBV3470
      IHLF=0                                                            DLBV3480
      DO 401 I=1,NDIM                                                   DLBV3490
      AUX(16,I)=0.D0                                                    DLBV3500
      AUX(1,I)=Y(I)                                                     DLBV3510
  401 AUX(8,I)=DERY(I)                                                  DLBV3520
      ISW1=1                                                            DLBV3530
      GOTO 500                                                          DLBV3540
C                                                                       DLBV3550
  402 X=X+H                                                             DLBV3560
      DO 403 I=1,NDIM                                                   DLBV3570
  403 AUX(2,I)=Y(I)                                                     DLBV3580
C                                                                       DLBV3590
C     INCREMENT H IS TESTED BY MEANS OF BISECTION                       DLBV3600
  404 IHLF=IHLF+1                                                       DLBV3610
      X=X-H                                                             DLBV3620
      DO 405 I=1,NDIM                                                   DLBV3630
  405 AUX(4,I)=AUX(2,I)                                                 DLBV3640
      H=.5D0*H                                                          DLBV3650
      N=1                                                               DLBV3660
      ISW1=2                                                            DLBV3670
      GOTO 500                                                          DLBV3680
C                                                                       DLBV3690
  406 X=X+H                                                             DLBV3700
      ISW2=4                                                            DLBV3710
      GOTO 200                                                          DLBV3720
  407 N=2                                                               DLBV3730
      DO 408 I=1,NDIM                                                   DLBV3740
      AUX(2,I)=Y(I)                                                     DLBV3750
  408 AUX(9,I)=DERY(I)                                                  DLBV3760
      ISW1=3                                                            DLBV3770
      GOTO 500                                                          DLBV3780
C                                                                       DLBV3790
C     TEST ON SATISFACTORY ACCURACY                                     DLBV3800
  409 DO 414 I=1,NDIM                                                   DLBV3810
      Z=DABS(Y(I))                                                      DLBV3820
      IF(Z-1.D0)410,411,411                                             DLBV3830
  410 Z=1.D0                                                            DLBV3840
  411 DELT=.066666666666666667D0*DABS(Y(I)-AUX(4,I))                    DLBV3850
      IF(ISW)413,413,412                                                DLBV3860
  412 DELT=AUX(15,I)*DELT                                               DLBV3870
  413 IF(DELT-Z*PRMT(4))414,414,429                                     DLBV3880
  414 CONTINUE                                                          DLBV3890
C                                                                       DLBV3900
C     SATISFACTORY ACCURACY AFTER LESS THAN 11 BISECTIONS               DLBV3910
      X=X+H                                                             DLBV3920
      ISW2=5                                                            DLBV3930
      GOTO 200                                                          DLBV3940
  415 DO 416 I=1,NDIM                                                   DLBV3950
      AUX(3,I)=Y(I)                                                     DLBV3960
  416 AUX(10,I)=DERY(I)                                                 DLBV3970
      N=3                                                               DLBV3980
      ISW1=4                                                            DLBV3990
      GOTO 500                                                          DLBV4000
C                                                                       DLBV4010
  417 N=1                                                               DLBV4020
      X=X+H                                                             DLBV4030
      ISW2=6                                                            DLBV4040
      GOTO 200                                                          DLBV4050
  418 X=XST                                                             DLBV4060
      DO 419 I=1,NDIM                                                   DLBV4070
      AUX(11,I)=DERY(I)                                                 DLBV4080
  4190Y(I)=AUX(1,I)+H*(.375D0*AUX(8,I)+.7916666666666667D0*AUX(9,I)     DLBV4090
     1-.20833333333333333D0*AUX(10,I)+.041666666666666667D0*DERY(I))    DLBV4100
  420 X=X+H                                                             DLBV4110
      N=N+1                                                             DLBV4120
      ISW2=11                                                           DLBV4130
      GOTO 200                                                          DLBV4140
  421 ISW3=0                                                            DLBV4150
      GOTO 300                                                          DLBV4160
  422 IF(N-4)423,600,600                                                DLBV4170
  423 DO 424 I=1,NDIM                                                   DLBV4180
      AUX(N,I)=Y(I)                                                     DLBV4190
  424 AUX(N+7,I)=DERY(I)                                                DLBV4200
      IF(N-3)425,427,600                                                DLBV4210
C                                                                       DLBV4220
  425 DO 426 I=1,NDIM                                                   DLBV4230
      DELT=AUX(9,I)+AUX(9,I)                                            DLBV4240
      DELT=DELT+DELT                                                    DLBV4250
  426 Y(I)=AUX(1,I)+.33333333333333333D0*H*(AUX(8,I)+DELT+AUX(10,I))    DLBV4260
      GOTO 420                                                          DLBV4270
C                                                                       DLBV4280
  427 DO 428 I=1,NDIM                                                   DLBV4290
      DELT=AUX(9,I)+AUX(10,I)                                           DLBV4300
      DELT=DELT+DELT+DELT                                               DLBV4310
  428 Y(I)=AUX(1,I)+.375D0*H*(AUX(8,I)+DELT+AUX(11,I))                  DLBV4320
      GOTO 420                                                          DLBV4330
C                                                                       DLBV4340
C     NO SATISFACTORY ACCURACY. H MUST BE HALVED.                       DLBV4350
  429 IF(IHLF-10)404,430,430                                            DLBV4360
C                                                                       DLBV4370
C     NO SATISFACTORY ACCURACY AFTER 10 BISECTIONS. ERROR MESSAGE.      DLBV4380
  430 IHLF=11                                                           DLBV4390
      X=X+H                                                             DLBV4400
      IF(ISW)105,105,114                                                DLBV4410
C                                                                       DLBV4420
C     THIS PART OF LINEAR INITIAL VALUE PROBLEM COMPUTES                DLBV4430
C     STARTING VALUES BY MEANS OF RUNGE-KUTTA METHOD.                   DLBV4440
  500 Z=X                                                               DLBV4450
      DO 501 I=1,NDIM                                                   DLBV4460
      X=H*AUX(N+7,I)                                                    DLBV4470
      AUX(5,I)=X                                                        DLBV4480
  501 Y(I)=AUX(N,I)+.4D0*X                                              DLBV4490
C                                                                       DLBV4500
      X=Z+.4D0*H                                                        DLBV4510
      ISW2=1                                                            DLBV4520
      GOTO 200                                                          DLBV4530
  502 DO 503 I=1,NDIM                                                   DLBV4540
      X=H*DERY(I)                                                       DLBV4550
      AUX(6,I)=X                                                        DLBV4560
  503 Y(I)=AUX(N,I)+.29697760924775360D0*AUX(5,I)+.15875964497103583D0*XDLBV4570
C                                                                       DLBV4580
      X=Z+.45573725421878943D0*H                                        DLBV4590
      ISW2=2                                                            DLBV4600
      GOTO 200                                                          DLBV4610
  504 DO 505 I=1,NDIM                                                   DLBV4620
      X=H*DERY(I)                                                       DLBV4630
      AUX(7,I)=X                                                        DLBV4640
  505 Y(I)=AUX(N,I)+.21810038822592047D0*AUX(5,I)-3.0509651486929308D0* DLBV4650
     1AUX(6,I)+3.8328647604670103D0*X                                   DLBV4660
C                                                                       DLBV4670
      X=Z+H                                                             DLBV4680
      ISW2=3                                                            DLBV4690
      GOTO 200                                                          DLBV4700
  506 DO 507 I=1,NDIM                                                   DLBV4710
  5070Y(I)=AUX(N,I)+.17476028226269037D0*AUX(5,I)-.55148066287873294D0* DLBV4720
     1AUX(6,I)+1.2055355993965235D0*AUX(7,I)+.17118478121951903D0*      DLBV4730
     2H*DERY(I)                                                         DLBV4740
      X=Z                                                               DLBV4750
      GOTO(402,406,409,417),ISW1                                        DLBV4760
C                                                                       DLBV4770
C     POSSIBLE BREAK-POINT FOR LINKAGE                                  DLBV4780
C                                                                       DLBV4790
C     STARTING VALUES ARE COMPUTED.                                     DLBV4800
C     NOW START HAMMINGS MODIFIED PREDICTOR-CORRECTOR METHOD.           DLBV4810
  600 ISTEP=3                                                           DLBV4820
  601 IF(N-8)604,602,604                                                DLBV4830
C                                                                       DLBV4840
C     N=8 CAUSES THE ROWS OF AUX TO CHANGE THEIR STORAGE LOCATIONS      DLBV4850
  602 DO 603 N=2,7                                                      DLBV4860
      DO 603 I=1,NDIM                                                   DLBV4870
      AUX(N-1,I)=AUX(N,I)                                               DLBV4880
  603 AUX(N+6,I)=AUX(N+7,I)                                             DLBV4890
      N=7                                                               DLBV4900
C                                                                       DLBV4910
C     N LESS THAN 8 CAUSES N+1 TO GET N                                 DLBV4920
  604 N=N+1                                                             DLBV4930
C                                                                       DLBV4940
C     COMPUTATION OF NEXT VECTOR Y                                      DLBV4950
      DO 605 I=1,NDIM                                                   DLBV4960
      AUX(N-1,I)=Y(I)                                                   DLBV4970
  605 AUX(N+6,I)=DERY(I)                                                DLBV4980
      X=X+H                                                             DLBV4990
  606 ISTEP=ISTEP+1                                                     DLBV5000
      DO 607 I=1,NDIM                                                   DLBV5010
     0DELT=AUX(N-4,I)+1.3333333333333333D0*H*(AUX(N+6,I)+AUX(N+6,I)-    DLBV5020
     1AUX(N+5,I)+AUX(N+4,I)+AUX(N+4,I))                                 DLBV5030
      Y(I)=DELT-.9256198347107438D0*AUX(16,I)                           DLBV5040
  607 AUX(16,I)=DELT                                                    DLBV5050
C     PREDICTOR IS NOW GENERATED IN ROW 16 OF AUX, MODIFIED PREDICTOR   DLBV5060
C     IS GENERATED IN Y. DELT MEANS AN AUXILIARY STORAGE.               DLBV5070
C                                                                       DLBV5080
      ISW2=7                                                            DLBV5090
      GOTO 200                                                          DLBV5100
C     DERIVATIVE OF MODIFIED PREDICTOR IS GENERATED IN DERY.            DLBV5110
C                                                                       DLBV5120
  608 DO 609 I=1,NDIM                                                   DLBV5130
     0DELT=.125D0*(9.D0*AUX(N-1,I)-AUX(N-3,I)+3.D0*H*(DERY(I)+AUX(N+6,I)DLBV5140
     1+AUX(N+6,I)-AUX(N+5,I)))                                          DLBV5150
      AUX(16,I)=AUX(16,I)-DELT                                          DLBV5160
  609 Y(I)=DELT+.07438016528925620D0*AUX(16,I)                          DLBV5170
C                                                                       DLBV5180
C     TEST WHETHER H MUST BE HALVED OR DOUBLED                          DLBV5190
      DELT=0.D0                                                         DLBV5200
      DO 616 I=1,NDIM                                                   DLBV5210
      Z=DABS(Y(I))                                                      DLBV5220
      IF(Z-1.D0)610,611,611                                             DLBV5230
  610 Z=1.D0                                                            DLBV5240
  611 Z=DABS(AUX(16,I))/Z                                               DLBV5250
      IF(ISW)613,613,612                                                DLBV5260
  612 Z=AUX(15,I)*Z                                                     DLBV5270
  613 IF(Z-PRMT(4))614,614,628                                          DLBV5280
  614 IF(DELT-Z)615,616,616                                             DLBV5290
  615 DELT=Z                                                            DLBV5300
  616 CONTINUE                                                          DLBV5310
C                                                                       DLBV5320
C     H MUST NOT BE HALVED. THAT MEANS Y(I) ARE GOOD.                   DLBV5330
      ISW2=8                                                            DLBV5340
      GOTO 200                                                          DLBV5350
  617 ISW3=1                                                            DLBV5360
      GOTO 300                                                          DLBV5370
  618 IF(H*(X-XEND))619,621,621                                         DLBV5380
  619 IF(DABS(X-XEND)-.1D0*DABS(H))621,620,620                          DLBV5390
  620 IF(DELT-.02D0*PRMT(4))622,622,601                                 DLBV5400
  621 IF(ISW)105,105,117                                                DLBV5410
C                                                                       DLBV5420
C     H COULD BE DOUBLED IF ALL NECESSARY PRECEEDING VALUES ARE         DLBV5430
C     AVAILABLE.                                                        DLBV5440
  622 IF(IHLF)601,601,623                                               DLBV5450
  623 IF(N-7)601,624,624                                                DLBV5460
  624 IF(ISTEP-4)601,625,625                                            DLBV5470
  625 IMOD=ISTEP/2                                                      DLBV5480
      IF(ISTEP-IMOD-IMOD)601,626,601                                    DLBV5490
  626 H=H+H                                                             DLBV5500
      IHLF=IHLF-1                                                       DLBV5510
      ISTEP=0                                                           DLBV5520
      DO 627 I=1,NDIM                                                   DLBV5530
      AUX(N-1,I)=AUX(N-2,I)                                             DLBV5540
      AUX(N-2,I)=AUX(N-4,I)                                             DLBV5550
      AUX(N-3,I)=AUX(N-6,I)                                             DLBV5560
      AUX(N+6,I)=AUX(N+5,I)                                             DLBV5570
      AUX(N+5,I)=AUX(N+3,I)                                             DLBV5580
      AUX(N+4,I)=AUX(N+1,I)                                             DLBV5590
      DELT=AUX(N+6,I)+AUX(N+5,I)                                        DLBV5600
      DELT=DELT+DELT+DELT                                               DLBV5610
  6270AUX(16,I)=8.962962962962963D0*(Y(I)-AUX(N-3,I))                   DLBV5620
     1-3.3611111111111111D0*H*(DERY(I)+DELT+AUX(N+4,I))                 DLBV5630
      GOTO 601                                                          DLBV5640
C                                                                       DLBV5650
C     H MUST BE HALVED                                                  DLBV5660
  628 IHLF=IHLF+1                                                       DLBV5670
      IF(IHLF-10)630,630,629                                            DLBV5680
  629 IF(ISW)105,105,114                                                DLBV5690
  630 H=.5D0*H                                                          DLBV5700
      ISTEP=0                                                           DLBV5710
      DO 631 I=1,NDIM                                                   DLBV5720
     0Y(I)=.390625D-2*(8.D1*AUX(N-1,I)+135.D0*AUX(N-2,I)+4.D1*AUX(N-3,I)DLBV5730
     1+AUX(N-4,I))-.1171875D0*(AUX(N+6,I)-6.D0*AUX(N+5,I)-AUX(N+4,I))*H DLBV5740
     0AUX(N-4,I)=.390625D-2*(12.D0*AUX(N-1,I)+135.D0*AUX(N-2,I)+        DLBV5750
     1108.D0*AUX(N-3,I)+AUX(N-4,I))-.0234375D0*(AUX(N+6,I)+             DLBV5760
     218.D0*AUX(N+5,I)-9.D0*AUX(N+4,I))*H                               DLBV5770
      AUX(N-3,I)=AUX(N-2,I)                                             DLBV5780
  631 AUX(N+4,I)=AUX(N+5,I)                                             DLBV5790
      DELT=X-H                                                          DLBV5800
      X=DELT-(H+H)                                                      DLBV5810
      ISW2=9                                                            DLBV5820
      GOTO 200                                                          DLBV5830
  632 DO 633 I=1,NDIM                                                   DLBV5840
      AUX(N-2,I)=Y(I)                                                   DLBV5850
      AUX(N+5,I)=DERY(I)                                                DLBV5860
  633 Y(I)=AUX(N-4,I)                                                   DLBV5870
      X=X-(H+H)                                                         DLBV5880
      ISW2=10                                                           DLBV5890
      GOTO 200                                                          DLBV5900
  634 X=DELT                                                            DLBV5910
      DO 635 I=1,NDIM                                                   DLBV5920
      DELT=AUX(N+5,I)+AUX(N+4,I)                                        DLBV5930
      DELT=DELT+DELT+DELT                                               DLBV5940
     0AUX(16,I)=8.962962962962963D0*(AUX(N-1,I)-Y(I))                   DLBV5950
     1-3.3611111111111111D0*H*(AUX(N+6,I)+DELT+DERY(I))                 DLBV5960
  635 AUX(N+3,I)=DERY(I)                                                DLBV5970
      GOTO 606                                                          DLBV5980
C     END OF INITIAL VALUE PROBLEM                                      DLBV5990
      END                                                               DLBV6000
