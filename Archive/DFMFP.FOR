C                                                                       DFMF  10
C     ..................................................................DFMF  20
C                                                                       DFMF  30
C        SUBROUTINE DFMFP                                               DFMF  40
C                                                                       DFMF  50
C        PURPOSE                                                        DFMF  60
C           TO FIND A LOCAL MINIMUM OF A FUNCTION OF SEVERAL VARIABLES  DFMF  70
C           BY THE METHOD OF FLETCHER AND POWELL                        DFMF  80
C                                                                       DFMF  90
C        USAGE                                                          DFMF 100
C           CALL DFMFP(FUNCT,N,X,F,G,EST,EPS,LIMIT,IER,H)               DFMF 110
C                                                                       DFMF 120
C        DESCRIPTION OF PARAMETERS                                      DFMF 130
C           FUNCT  - USER-WRITTEN SUBROUTINE CONCERNING THE FUNCTION TO DFMF 140
C                    BE MINIMIZED. IT MUST BE OF THE FORM               DFMF 150
C                    SUBROUTINE FUNCT(N,ARG,VAL,GRAD)                   DFMF 160
C                    AND MUST SERVE THE FOLLOWING PURPOSE               DFMF 170
C                    FOR EACH N-DIMENSIONAL ARGUMENT VECTOR  ARG,       DFMF 180
C                    FUNCTION VALUE AND GRADIENT VECTOR MUST BE COMPUTEDDFMF 190
C                    AND, ON RETURN, STORED IN VAL AND GRAD RESPECTIVELYDFMF 200
C                    ARG,VAL AND GRAD MUST BE OF DOUBLE PRECISION.      DFMF 210
C           N      - NUMBER OF VARIABLES                                DFMF 220
C           X      - VECTOR OF DIMENSION N CONTAINING THE INITIAL       DFMF 230
C                    ARGUMENT WHERE THE ITERATION STARTS. ON RETURN,    DFMF 240
C                    X HOLDS THE ARGUMENT CORRESPONDING TO THE          DFMF 250
C                    COMPUTED MINIMUM FUNCTION VALUE                    DFMF 260
C                    DOUBLE PRECISION VECTOR.                           DFMF 270
C           F      - SINGLE VARIABLE CONTAINING THE MINIMUM FUNCTION    DFMF 280
C                    VALUE ON RETURN, I.E. F=F(X).                      DFMF 290
C                    DOUBLE PRECISION VARIABLE.                         DFMF 300
C           G      - VECTOR OF DIMENSION N CONTAINING THE GRADIENT      DFMF 310
C                    VECTOR CORRESPONDING TO THE MINIMUM ON RETURN,     DFMF 320
C                    I.E. G=G(X).                                       DFMF 330
C                    DOUBLE PRECISION VECTOR.                           DFMF 340
C           EST    - IS AN ESTIMATE OF THE MINIMUM FUNCTION VALUE.      DFMF 350
C                    SINGLE PRECISION VARIABLE.                         DFMF 360
C           EPS    - TESTVALUE REPRESENTING THE EXPECTED ABSOLUTE ERROR.DFMF 370
C                    A REASONABLE CHOICE IS 10**(-16), I.E.             DFMF 380
C                    SOMEWHAT GREATER THAN 10**(-D), WHERE D IS THE     DFMF 390
C                    NUMBER OF SIGNIFICANT DIGITS IN FLOATING POINT     DFMF 400
C                    REPRESENTATION.                                    DFMF 410
C                    SINGLE PRECISION VARIABLE.                         DFMF 420
C           LIMIT  - MAXIMUM NUMBER OF ITERATIONS.                      DFMF 430
C           IER    - ERROR PARAMETER                                    DFMF 440
C                    IER = 0 MEANS CONVERGENCE WAS OBTAINED             DFMF 450
C                    IER = 1 MEANS NO CONVERGENCE IN LIMIT ITERATIONS   DFMF 460
C                    IER =-1 MEANS ERRORS IN GRADIENT CALCULATION       DFMF 470
C                    IER = 2 MEANS LINEAR SEARCH TECHNIQUE INDICATES    DFMF 480
C                    IT IS LIKELY THAT THERE EXISTS NO MINIMUM.         DFMF 490
C           H      - WORKING STORAGE OF DIMENSION N*(N+7)/2.            DFMF 500
C                    DOUBLE PRECISION ARRAY.                            DFMF 510
C                                                                       DFMF 520
C        REMARKS                                                        DFMF 530
C            I) THE SUBROUTINE NAME REPLACING THE DUMMY ARGUMENT  FUNCT DFMF 540
C               MUST BE DECLARED AS EXTERNAL IN THE CALLING PROGRAM.    DFMF 550
C           II) IER IS SET TO 2 IF , STEPPING IN ONE OF THE COMPUTED    DFMF 560
C               DIRECTIONS, THE FUNCTION WILL NEVER INCREASE WITHIN     DFMF 570
C               A TOLERABLE RANGE OF ARGUMENT.                          DFMF 580
C               IER = 2 MAY OCCUR ALSO IF THE INTERVAL WHERE F          DFMF 590
C               INCREASES IS SMALL AND THE INITIAL ARGUMENT WAS         DFMF 600
C               RELATIVELY FAR AWAY FROM THE MINIMUM SUCH THAT THE      DFMF 610
C               MINIMUM WAS OVERLEAPED. THIS IS DUE TO THE SEARCH       DFMF 620
C               TECHNIQUE WHICH DOUBLES THE STEPSIZE UNTIL A POINT      DFMF 630
C               IS FOUND WHERE THE FUNCTION INCREASES.                  DFMF 640
C                                                                       DFMF 650
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DFMF 660
C           FUNCT                                                       DFMF 670
C                                                                       DFMF 680
C        METHOD                                                         DFMF 690
C           THE METHOD IS DESCRIBED IN THE FOLLOWING ARTICLE            DFMF 700
C           R. FLETCHER AND M.J.D. POWELL, A RAPID DESCENT METHOD FOR   DFMF 710
C           MINIMIZATION,                                               DFMF 720
C           COMPUTER JOURNAL VOL.6, ISS. 2, 1963, PP.163-168.           DFMF 730
C                                                                       DFMF 740
C     ..................................................................DFMF 750
C                                                                       DFMF 760
      SUBROUTINE DFMFP(FUNCT,N,X,F,G,EST,EPS,LIMIT,IER,H)               DFMF 770
C                                                                       DFMF 780
C        DIMENSIONED DUMMY VARIABLES                                    DFMF 790
      DIMENSION H(1),X(1),G(1)                                          DFMF 800
      DOUBLE PRECISION X,F,FX,FY,OLDF,HNRM,GNRM,H,G,DX,DY,ALFA,DALFA,   DFMF 810
     1AMBDA,T,Z,W                                                       DFMF 820
C                                                                       DFMF 830
C        COMPUTE FUNCTION VALUE AND GRADIENT VECTOR FOR INITIAL ARGUMENTDFMF 840
      CALL FUNCT(N,X,F,G)                                               DFMF 850
C                                                                       DFMF 860
C        RESET ITERATION COUNTER AND GENERATE IDENTITY MATRIX           DFMF 870
      IER=0                                                             DFMF 880
      KOUNT=0                                                           DFMF 890
      N2=N+N                                                            DFMF 900
      N3=N2+N                                                           DFMF 910
      N31=N3+1                                                          DFMF 920
    1 K=N31                                                             DFMF 930
      DO 4 J=1,N                                                        DFMF 940
      H(K)=1.D0                                                         DFMF 950
      NJ=N-J                                                            DFMF 960
      IF(NJ)5,5,2                                                       DFMF 970
    2 DO 3 L=1,NJ                                                       DFMF 980
      KL=K+L                                                            DFMF 990
    3 H(KL)=0.D0                                                        DFMF1000
    4 K=KL+1                                                            DFMF1010
C                                                                       DFMF1020
C        START ITERATION LOOP                                           DFMF1030
    5 KOUNT=KOUNT +1                                                    DFMF1040
C                                                                       DFMF1050
C        SAVE FUNCTION VALUE, ARGUMENT VECTOR AND GRADIENT VECTOR       DFMF1060
      OLDF=F                                                            DFMF1070
      DO 9 J=1,N                                                        DFMF1080
      K=N+J                                                             DFMF1090
      H(K)=G(J)                                                         DFMF1100
      K=K+N                                                             DFMF1110
      H(K)=X(J)                                                         DFMF1120
C                                                                       DFMF1130
C        DETERMINE DIRECTION VECTOR H                                   DFMF1140
      K=J+N3                                                            DFMF1150
      T=0.D0                                                            DFMF1160
      DO 8 L=1,N                                                        DFMF1170
      T=T-G(L)*H(K)                                                     DFMF1180
      IF(L-J)6,7,7                                                      DFMF1190
    6 K=K+N-L                                                           DFMF1200
      GO TO 8                                                           DFMF1210
    7 K=K+1                                                             DFMF1220
    8 CONTINUE                                                          DFMF1230
    9 H(J)=T                                                            DFMF1240
C                                                                       DFMF1250
C        CHECK WHETHER FUNCTION WILL DECREASE STEPPING ALONG H.         DFMF1260
      DY=0.D0                                                           DFMF1270
      HNRM=0.D0                                                         DFMF1280
      GNRM=0.D0                                                         DFMF1290
C                                                                       DFMF1300
C        CALCULATE DIRECTIONAL DERIVATIVE AND TESTVALUES FOR DIRECTION  DFMF1310
C        VECTOR H AND GRADIENT VECTOR G.                                DFMF1320
      DO 10 J=1,N                                                       DFMF1330
      HNRM=HNRM+DABS(H(J))                                              DFMF1340
      GNRM=GNRM+DABS(G(J))                                              DFMF1350
   10 DY=DY+H(J)*G(J)                                                   DFMF1360
C                                                                       DFMF1370
C        REPEAT SEARCH IN DIRECTION OF STEEPEST DESCENT IF DIRECTIONAL  DFMF1380
C        DERIVATIVE APPEARS TO BE POSITIVE OR ZERO.                     DFMF1390
      IF(DY)11,51,51                                                    DFMF1400
C                                                                       DFMF1410
C        REPEAT SEARCH IN DIRECTION OF STEEPEST DESCENT IF DIRECTION    DFMF1420
C        VECTOR H IS SMALL COMPARED TO GRADIENT VECTOR G.               DFMF1430
   11 IF(HNRM/GNRM-EPS)51,51,12                                         DFMF1440
C                                                                       DFMF1450
C        SEARCH MINIMUM ALONG DIRECTION H                               DFMF1460
C                                                                       DFMF1470
C        SEARCH ALONG H FOR POSITIVE DIRECTIONAL DERIVATIVE             DFMF1480
   12 FY=F                                                              DFMF1490
      ALFA=2.D0*(EST-F)/DY                                              DFMF1500
      AMBDA=1.D0                                                        DFMF1510
C                                                                       DFMF1520
C        USE ESTIMATE FOR STEPSIZE ONLY IF IT IS POSITIVE AND LESS THAN DFMF1530
C        1. OTHERWISE TAKE 1. AS STEPSIZE                               DFMF1540
      IF(ALFA)15,15,13                                                  DFMF1550
   13 IF(ALFA-AMBDA)14,15,15                                            DFMF1560
   14 AMBDA=ALFA                                                        DFMF1570
   15 ALFA=0.D0                                                         DFMF1580
C                                                                       DFMF1590
C        SAVE FUNCTION AND DERIVATIVE VALUES FOR OLD ARGUMENT           DFMF1600
   16 FX=FY                                                             DFMF1610
      DX=DY                                                             DFMF1620
C                                                                       DFMF1630
C        STEP ARGUMENT ALONG H                                          DFMF1640
      DO 17 I=1,N                                                       DFMF1650
   17 X(I)=X(I)+AMBDA*H(I)                                              DFMF1660
C                                                                       DFMF1670
C        COMPUTE FUNCTION VALUE AND GRADIENT FOR NEW ARGUMENT           DFMF1680
      CALL FUNCT(N,X,F,G)                                               DFMF1690
      FY=F                                                              DFMF1700
C                                                                       DFMF1710
C        COMPUTE DIRECTIONAL DERIVATIVE DY FOR NEW ARGUMENT.  TERMINATE DFMF1720
C        SEARCH, IF DY IS POSITIVE. IF DY IS ZERO THE MINIMUM IS FOUND  DFMF1730
      DY=0.D0                                                           DFMF1740
      DO 18 I=1,N                                                       DFMF1750
   18 DY=DY+G(I)*H(I)                                                   DFMF1760
      IF(DY)19,36,22                                                    DFMF1770
C                                                                       DFMF1780
C        TERMINATE SEARCH ALSO IF THE FUNCTION VALUE INDICATES THAT     DFMF1790
C        A MINIMUM HAS BEEN PASSED                                      DFMF1800
   19 IF(FY-FX)20,22,22                                                 DFMF1810
C                                                                       DFMF1820
C        REPEAT SEARCH AND DOUBLE STEPSIZE FOR FURTHER SEARCHES         DFMF1830
   20 AMBDA=AMBDA+ALFA                                                  DFMF1840
      ALFA=AMBDA                                                        DFMF1850
C        END OF SEARCH LOOP                                             DFMF1860
C                                                                       DFMF1870
C        TERMINATE IF THE CHANGE IN ARGUMENT GETS VERY LARGE            DFMF1880
      IF(HNRM*AMBDA-1.D10)16,16,21                                      DFMF1890
C                                                                       DFMF1900
C        LINEAR SEARCH TECHNIQUE INDICATES THAT NO MINIMUM EXISTS       DFMF1910
   21 IER=2                                                             DFMF1920
      RETURN                                                            DFMF1930
C                                                                       DFMF1940
C        INTERPOLATE CUBICALLY IN THE INTERVAL DEFINED BY THE SEARCH    DFMF1950
C        ABOVE AND COMPUTE THE ARGUMENT X FOR WHICH THE INTERPOLATION   DFMF1960
C        POLYNOMIAL IS MINIMIZED                                        DFMF1970
   22 T=0.D0                                                            DFMF1980
   23 IF(AMBDA)24,36,24                                                 DFMF1990
   24 Z=3.D0*(FX-FY)/AMBDA+DX+DY                                        DFMF2000
      ALFA=DMAX1(DABS(Z),DABS(DX),DABS(DY))                             DFMF2010
      DALFA=Z/ALFA                                                      DFMF2020
      DALFA=DALFA*DALFA-DX/ALFA*DY/ALFA                                 DFMF2030
      IF(DALFA)51,25,25                                                 DFMF2040
   25 W=ALFA*DSQRT(DALFA)                                               DFMF2050
      ALFA=DY-DX+W+W                                                    DFMF2060
      IF(ALFA) 250,251,250                                              DFMF2061
  250 ALFA=(DY-Z+W)/ALFA                                                DFMF2062
      GO TO 252                                                         DFMF2063
  251 ALFA=(Z+DY-W)/(Z+DX+Z+DY)                                         DFMF2064
  252 ALFA=ALFA*AMBDA                                                   DFMF2065
      DO 26 I=1,N                                                       DFMF2070
   26 X(I)=X(I)+(T-ALFA)*H(I)                                           DFMF2080
C                                                                       DFMF2090
C        TERMINATE, IF THE VALUE OF THE ACTUAL FUNCTION AT X IS LESS    DFMF2100
C        THAN THE FUNCTION VALUES AT THE INTERVAL ENDS. OTHERWISE REDUCEDFMF2110
C        THE INTERVAL BY CHOOSING ONE END-POINT EQUAL TO X AND REPEAT   DFMF2120
C        THE INTERPOLATION.  WHICH END-POINT IS CHOOSEN DEPENDS ON THE  DFMF2130
C        VALUE OF THE FUNCTION AND ITS GRADIENT AT X                    DFMF2140
C                                                                       DFMF2150
      CALL FUNCT(N,X,F,G)                                               DFMF2160
      IF(F-FX)27,27,28                                                  DFMF2170
   27 IF(F-FY)36,36,28                                                  DFMF2180
   28 DALFA=0.D0                                                        DFMF2190
      DO 29 I=1,N                                                       DFMF2200
   29 DALFA=DALFA+G(I)*H(I)                                             DFMF2210
      IF(DALFA)30,33,33                                                 DFMF2220
   30 IF(F-FX)32,31,33                                                  DFMF2230
   31 IF(DX-DALFA)32,36,32                                              DFMF2240
   32 FX=F                                                              DFMF2250
      DX=DALFA                                                          DFMF2260
      T=ALFA                                                            DFMF2270
      AMBDA=ALFA                                                        DFMF2280
      GO TO 23                                                          DFMF2290
   33 IF(FY-F)35,34,35                                                  DFMF2300
   34 IF(DY-DALFA)35,36,35                                              DFMF2310
   35 FY=F                                                              DFMF2320
      DY=DALFA                                                          DFMF2330
      AMBDA=AMBDA-ALFA                                                  DFMF2340
      GO TO 22                                                          DFMF2350
C                                                                       DFMF2360
C        TERMINATE, IF FUNCTION HAS NOT DECREASED DURING LAST ITERATION DFMF2370
   36 IF(OLDF-F+EPS)51,38,38                                            DFMF2380
C                                                                       DFMF2390
C        COMPUTE DIFFERENCE VECTORS OF ARGUMENT AND GRADIENT FROM       DFMF2400
C        TWO CONSECUTIVE ITERATIONS                                     DFMF2410
   38 DO 37 J=1,N                                                       DFMF2420
      K=N+J                                                             DFMF2430
      H(K)=G(J)-H(K)                                                    DFMF2440
      K=N+K                                                             DFMF2450
   37 H(K)=X(J)-H(K)                                                    DFMF2460
C                                                                       DFMF2470
C        TEST LENGTH OF ARGUMENT DIFFERENCE VECTOR AND DIRECTION VECTOR DFMF2480
C        IF AT LEAST N ITERATIONS HAVE BEEN EXECUTED. TERMINATE, IF     DFMF2490
C        BOTH ARE LESS THAN  EPS                                        DFMF2500
      IER=0                                                             DFMF2510
      IF(KOUNT-N)42,39,39                                               DFMF2520
   39 T=0.D0                                                            DFMF2530
      Z=0.D0                                                            DFMF2540
      DO 40 J=1,N                                                       DFMF2550
      K=N+J                                                             DFMF2560
      W=H(K)                                                            DFMF2570
      K=K+N                                                             DFMF2580
      T=T+DABS(H(K))                                                    DFMF2590
   40 Z=Z+W*H(K)                                                        DFMF2600
      IF(HNRM-EPS)41,41,42                                              DFMF2610
   41 IF(T-EPS)56,56,42                                                 DFMF2620
C                                                                       DFMF2630
C        TERMINATE, IF NUMBER OF ITERATIONS WOULD EXCEED  LIMIT         DFMF2640
   42 IF(KOUNT-LIMIT)43,50,50                                           DFMF2650
C                                                                       DFMF2660
C        PREPARE UPDATING OF MATRIX H                                   DFMF2670
   43 ALFA=0.D0                                                         DFMF2680
      DO 47 J=1,N                                                       DFMF2690
      K=J+N3                                                            DFMF2700
      W=0.D0                                                            DFMF2710
      DO 46 L=1,N                                                       DFMF2720
      KL=N+L                                                            DFMF2730
      W=W+H(KL)*H(K)                                                    DFMF2740
      IF(L-J)44,45,45                                                   DFMF2750
   44 K=K+N-L                                                           DFMF2760
      GO TO 46                                                          DFMF2770
   45 K=K+1                                                             DFMF2780
   46 CONTINUE                                                          DFMF2790
      K=N+J                                                             DFMF2800
      ALFA=ALFA+W*H(K)                                                  DFMF2810
   47 H(J)=W                                                            DFMF2820
C                                                                       DFMF2830
C        REPEAT SEARCH IN DIRECTION OF STEEPEST DESCENT IF RESULTS      DFMF2840
C        ARE NOT SATISFACTORY                                           DFMF2850
      IF(Z*ALFA)48,1,48                                                 DFMF2860
C                                                                       DFMF2870
C        UPDATE MATRIX H                                                DFMF2880
   48 K=N31                                                             DFMF2890
      DO 49 L=1,N                                                       DFMF2900
      KL=N2+L                                                           DFMF2910
      DO 49 J=L,N                                                       DFMF2920
      NJ=N2+J                                                           DFMF2930
      H(K)=H(K)+H(KL)*H(NJ)/Z-H(L)*H(J)/ALFA                            DFMF2940
   49 K=K+1                                                             DFMF2950
      GO TO 5                                                           DFMF2960
C        END OF ITERATION LOOP                                          DFMF2970
C                                                                       DFMF2980
C        NO CONVERGENCE AFTER  LIMIT  ITERATIONS                        DFMF2990
   50 IER=1                                                             DFMF3000
      RETURN                                                            DFMF3010
C                                                                       DFMF3020
C        RESTORE OLD VALUES OF FUNCTION AND ARGUMENTS                   DFMF3030
   51 DO 52 J=1,N                                                       DFMF3040
      K=N2+J                                                            DFMF3050
   52 X(J)=H(K)                                                         DFMF3060
      CALL FUNCT(N,X,F,G)                                               DFMF3070
C                                                                       DFMF3080
C        REPEAT SEARCH IN DIRECTION OF STEEPEST DESCENT IF DERIVATIVE   DFMF3090
C        FAILS TO BE SUFFICIENTLY SMALL                                 DFMF3100
      IF(GNRM-EPS)55,55,53                                              DFMF3110
C                                                                       DFMF3120
C        TEST FOR REPEATED FAILURE OF ITERATION                         DFMF3130
   53 IF(IER)56,54,54                                                   DFMF3140
   54 IER=-1                                                            DFMF3150
      GOTO 1                                                            DFMF3160
   55 IER=0                                                             DFMF3170
   56 RETURN                                                            DFMF3180
      END                                                               DFMF3190
