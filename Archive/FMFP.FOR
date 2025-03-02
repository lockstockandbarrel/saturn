C                                                                       FMFP  10
C     ..................................................................FMFP  20
C                                                                       FMFP  30
C        SUBROUTINE FMFP                                                FMFP  40
C                                                                       FMFP  50
C        PURPOSE                                                        FMFP  60
C           TO FIND A LOCAL MINIMUM OF A FUNCTION OF SEVERAL VARIABLES  FMFP  70
C           BY THE METHOD OF FLETCHER AND POWELL                        FMFP  80
C                                                                       FMFP  90
C        USAGE                                                          FMFP 100
C           CALL FMFP(FUNCT,N,X,F,G,EST,EPS,LIMIT,IER,H)                FMFP 110
C                                                                       FMFP 120
C        DESCRIPTION OF PARAMETERS                                      FMFP 130
C           FUNCT  - USER-WRITTEN SUBROUTINE CONCERNING THE FUNCTION TO FMFP 140
C                    BE MINIMIZED. IT MUST BE OF THE FORM               FMFP 150
C                    SUBROUTINE FUNCT(N,ARG,VAL,GRAD)                   FMFP 160
C                    AND MUST SERVE THE FOLLOWING PURPOSE               FMFP 170
C                    FOR EACH N-DIMENSIONAL ARGUMENT VECTOR  ARG,       FMFP 180
C                    FUNCTION VALUE AND GRADIENT VECTOR MUST BE COMPUTEDFMFP 190
C                    AND, ON RETURN, STORED IN VAL AND GRAD RESPECTIVELYFMFP 200
C           N      - NUMBER OF VARIABLES                                FMFP 210
C           X      - VECTOR OF DIMENSION N CONTAINING THE INITIAL       FMFP 220
C                    ARGUMENT WHERE THE ITERATION STARTS. ON RETURN,    FMFP 230
C                    X HOLDS THE ARGUMENT CORRESPONDING TO THE          FMFP 240
C                    COMPUTED MINIMUM FUNCTION VALUE                    FMFP 250
C           F      - SINGLE VARIABLE CONTAINING THE MINIMUM FUNCTION    FMFP 260
C                    VALUE ON RETURN, I.E. F=F(X).                      FMFP 270
C           G      - VECTOR OF DIMENSION N CONTAINING THE GRADIENT      FMFP 280
C                    VECTOR CORRESPONDING TO THE MINIMUM ON RETURN,     FMFP 290
C                    I.E. G=G(X).                                       FMFP 300
C           EST    - IS AN ESTIMATE OF THE MINIMUM FUNCTION VALUE.      FMFP 310
C           EPS    - TESTVALUE REPRESENTING THE EXPECTED ABSOLUTE ERROR.FMFP 320
C                    A REASONABLE CHOICE IS 10**(-6), I.E.              FMFP 330
C                    SOMEWHAT GREATER THAN 10**(-D), WHERE D IS THE     FMFP 340
C                    NUMBER OF SIGNIFICANT DIGITS IN FLOATING POINT     FMFP 350
C                    REPRESENTATION.                                    FMFP 360
C           LIMIT  - MAXIMUM NUMBER OF ITERATIONS.                      FMFP 370
C           IER    - ERROR PARAMETER                                    FMFP 380
C                    IER = 0 MEANS CONVERGENCE WAS OBTAINED             FMFP 390
C                    IER = 1 MEANS NO CONVERGENCE IN LIMIT ITERATIONS   FMFP 400
C                    IER =-1 MEANS ERRORS IN GRADIENT CALCULATION       FMFP 410
C                    IER = 2 MEANS LINEAR SEARCH TECHNIQUE INDICATES    FMFP 420
C                    IT IS LIKELY THAT THERE EXISTS NO MINIMUM.         FMFP 430
C           H      - WORKING STORAGE OF DIMENSION N*(N+7)/2.            FMFP 440
C                                                                       FMFP 450
C        REMARKS                                                        FMFP 460
C            I) THE SUBROUTINE NAME REPLACING THE DUMMY ARGUMENT  FUNCT FMFP 470
C               MUST BE DECLARED AS EXTERNAL IN THE CALLING PROGRAM.    FMFP 480
C           II) IER IS SET TO 2 IF , STEPPING IN ONE OF THE COMPUTED    FMFP 490
C               DIRECTIONS, THE FUNCTION WILL NEVER INCREASE WITHIN     FMFP 500
C               A TOLERABLE RANGE OF ARGUMENT.                          FMFP 510
C               IER = 2 MAY OCCUR ALSO IF THE INTERVAL WHERE F          FMFP 520
C               INCREASES IS SMALL AND THE INITIAL ARGUMENT WAS         FMFP 530
C               RELATIVELY FAR AWAY FROM THE MINIMUM SUCH THAT THE      FMFP 540
C               MINIMUM WAS OVERLEAPED. THIS IS DUE TO THE SEARCH       FMFP 550
C               TECHNIQUE WHICH DOUBLES THE STEPSIZE UNTIL A POINT      FMFP 560
C               IS FOUND WHERE THE FUNCTION INCREASES.                  FMFP 570
C                                                                       FMFP 580
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  FMFP 590
C           FUNCT                                                       FMFP 600
C                                                                       FMFP 610
C        METHOD                                                         FMFP 620
C           THE METHOD IS DESCRIBED IN THE FOLLOWING ARTICLE            FMFP 630
C           R. FLETCHER AND M.J.D. POWELL, A RAPID DESCENT METHOD FOR   FMFP 640
C           MINIMIZATION,                                               FMFP 650
C           COMPUTER JOURNAL VOL.6, ISS. 2, 1963, PP.163-168.           FMFP 660
C                                                                       FMFP 670
C     ..................................................................FMFP 680
C                                                                       FMFP 690
      SUBROUTINE FMFP(FUNCT,N,X,F,G,EST,EPS,LIMIT,IER,H)                FMFP 700
C                                                                       FMFP 710
C        DIMENSIONED DUMMY VARIABLES                                    FMFP 720
      DIMENSION H(1),X(1),G(1)                                          FMFP 730
C                                                                       FMFP 740
C        COMPUTE FUNCTION VALUE AND GRADIENT VECTOR FOR INITIAL ARGUMENTFMFP 750
      CALL FUNCT(N,X,F,G)                                               FMFP 760
C                                                                       FMFP 770
C        RESET ITERATION COUNTER AND GENERATE IDENTITY MATRIX           FMFP 780
      IER=0                                                             FMFP 790
      KOUNT=0                                                           FMFP 800
      N2=N+N                                                            FMFP 810
      N3=N2+N                                                           FMFP 820
      N31=N3+1                                                          FMFP 830
    1 K=N31                                                             FMFP 840
      DO 4 J=1,N                                                        FMFP 850
      H(K)=1.                                                           FMFP 860
      NJ=N-J                                                            FMFP 870
      IF(NJ)5,5,2                                                       FMFP 880
    2 DO 3 L=1,NJ                                                       FMFP 890
      KL=K+L                                                            FMFP 900
    3 H(KL)=0.                                                          FMFP 910
    4 K=KL+1                                                            FMFP 920
C                                                                       FMFP 930
C        START ITERATION LOOP                                           FMFP 940
    5 KOUNT=KOUNT +1                                                    FMFP 950
C                                                                       FMFP 960
C        SAVE FUNCTION VALUE, ARGUMENT VECTOR AND GRADIENT VECTOR       FMFP 970
      OLDF=F                                                            FMFP 980
      DO 9 J=1,N                                                        FMFP 990
      K=N+J                                                             FMFP1000
      H(K)=G(J)                                                         FMFP1010
      K=K+N                                                             FMFP1020
      H(K)=X(J)                                                         FMFP1030
C                                                                       FMFP1040
C        DETERMINE DIRECTION VECTOR H                                   FMFP1050
      K=J+N3                                                            FMFP1060
      T=0.                                                              FMFP1070
      DO 8 L=1,N                                                        FMFP1080
      T=T-G(L)*H(K)                                                     FMFP1090
      IF(L-J)6,7,7                                                      FMFP1100
    6 K=K+N-L                                                           FMFP1110
      GO TO 8                                                           FMFP1120
    7 K=K+1                                                             FMFP1130
    8 CONTINUE                                                          FMFP1140
    9 H(J)=T                                                            FMFP1150
C                                                                       FMFP1160
C        CHECK WHETHER FUNCTION WILL DECREASE STEPPING ALONG H.         FMFP1170
      DY=0.                                                             FMFP1180
      HNRM=0.                                                           FMFP1190
      GNRM=0.                                                           FMFP1200
C                                                                       FMFP1210
C        CALCULATE DIRECTIONAL DERIVATIVE AND TESTVALUES FOR DIRECTION  FMFP1220
C        VECTOR H AND GRADIENT VECTOR G.                                FMFP1230
      DO 10 J=1,N                                                       FMFP1240
      HNRM=HNRM+ABS(H(J))                                               FMFP1250
      GNRM=GNRM+ABS(G(J))                                               FMFP1260
   10 DY=DY+H(J)*G(J)                                                   FMFP1270
C                                                                       FMFP1280
C        REPEAT SEARCH IN DIRECTION OF STEEPEST DESCENT IF DIRECTIONAL  FMFP1290
C        DERIVATIVE APPEARS TO BE POSITIVE OR ZERO.                     FMFP1300
      IF(DY)11,51,51                                                    FMFP1310
C                                                                       FMFP1320
C        REPEAT SEARCH IN DIRECTION OF STEEPEST DESCENT IF DIRECTION    FMFP1330
C        VECTOR H IS SMALL COMPARED TO GRADIENT VECTOR G.               FMFP1340
   11 IF(HNRM/GNRM-EPS)51,51,12                                         FMFP1350
C                                                                       FMFP1360
C        SEARCH MINIMUM ALONG DIRECTION H                               FMFP1370
C                                                                       FMFP1380
C        SEARCH ALONG H FOR POSITIVE DIRECTIONAL DERIVATIVE             FMFP1390
   12 FY=F                                                              FMFP1400
      ALFA=2.*(EST-F)/DY                                                FMFP1410
      AMBDA=1.                                                          FMFP1420
C                                                                       FMFP1430
C        USE ESTIMATE FOR STEPSIZE ONLY IF IT IS POSITIVE AND LESS THAN FMFP1440
C        1. OTHERWISE TAKE 1. AS STEPSIZE                               FMFP1450
      IF(ALFA)15,15,13                                                  FMFP1460
   13 IF(ALFA-AMBDA)14,15,15                                            FMFP1470
   14 AMBDA=ALFA                                                        FMFP1480
   15 ALFA=0.                                                           FMFP1490
C                                                                       FMFP1500
C        SAVE FUNCTION AND DERIVATIVE VALUES FOR OLD ARGUMENT           FMFP1510
   16 FX=FY                                                             FMFP1520
      DX=DY                                                             FMFP1530
C                                                                       FMFP1540
C        STEP ARGUMENT ALONG H                                          FMFP1550
      DO 17 I=1,N                                                       FMFP1560
   17 X(I)=X(I)+AMBDA*H(I)                                              FMFP1570
C                                                                       FMFP1580
C        COMPUTE FUNCTION VALUE AND GRADIENT FOR NEW ARGUMENT           FMFP1590
      CALL FUNCT(N,X,F,G)                                               FMFP1600
      FY=F                                                              FMFP1610
C                                                                       FMFP1620
C        COMPUTE DIRECTIONAL DERIVATIVE DY FOR NEW ARGUMENT.  TERMINATE FMFP1630
C        SEARCH, IF DY IS POSITIVE. IF DY IS ZERO THE MINIMUM IS FOUND  FMFP1640
      DY=0.                                                             FMFP1650
      DO 18 I=1,N                                                       FMFP1660
   18 DY=DY+G(I)*H(I)                                                   FMFP1670
      IF(DY)19,36,22                                                    FMFP1680
C                                                                       FMFP1690
C        TERMINATE SEARCH ALSO IF THE FUNCTION VALUE INDICATES THAT     FMFP1700
C        A MINIMUM HAS BEEN PASSED                                      FMFP1710
   19 IF(FY-FX)20,22,22                                                 FMFP1720
C                                                                       FMFP1730
C        REPEAT SEARCH AND DOUBLE STEPSIZE FOR FURTHER SEARCHES         FMFP1740
   20 AMBDA=AMBDA+ALFA                                                  FMFP1750
      ALFA=AMBDA                                                        FMFP1760
C        END OF SEARCH LOOP                                             FMFP1770
C                                                                       FMFP1780
C        TERMINATE IF THE CHANGE IN ARGUMENT GETS VERY LARGE            FMFP1790
      IF(HNRM*AMBDA-1.E10)16,16,21                                      FMFP1800
C                                                                       FMFP1810
C        LINEAR SEARCH TECHNIQUE INDICATES THAT NO MINIMUM EXISTS       FMFP1820
   21 IER=2                                                             FMFP1830
      RETURN                                                            FMFP1840
C                                                                       FMFP1850
C        INTERPOLATE CUBICALLY IN THE INTERVAL DEFINED BY THE SEARCH    FMFP1860
C        ABOVE AND COMPUTE THE ARGUMENT X FOR WHICH THE INTERPOLATION   FMFP1870
C        POLYNOMIAL IS MINIMIZED                                        FMFP1880
   22 T=0.                                                              FMFP1890
   23 IF(AMBDA)24,36,24                                                 FMFP1900
   24 Z=3.*(FX-FY)/AMBDA+DX+DY                                          FMFP1910
      ALFA=AMAX1(ABS(Z),ABS(DX),ABS(DY))                                FMFP1920
      DALFA=Z/ALFA                                                      FMFP1930
      DALFA=DALFA*DALFA-DX/ALFA*DY/ALFA                                 FMFP1940
      IF(DALFA)51,25,25                                                 FMFP1950
   25 W=ALFA*SQRT(DALFA)                                                FMFP1960
      ALFA=DY-DX+W+W                                                    FMFP1970
      IF(ALFA) 250,251,250                                              FMFP1971
  250 ALFA=(DY-Z+W)/ALFA                                                FMFP1972
      GO TO 252                                                         FMFP1973
  251 ALFA=(Z+DY-W)/(Z+DX+Z+DY)                                         FMFP1974
  252 ALFA=ALFA*AMBDA                                                   FMFP1975
      DO 26 I=1,N                                                       FMFP1980
   26 X(I)=X(I)+(T-ALFA)*H(I)                                           FMFP1990
C                                                                       FMFP2000
C        TERMINATE, IF THE VALUE OF THE ACTUAL FUNCTION AT X IS LESS    FMFP2010
C        THAN THE FUNCTION VALUES AT THE INTERVAL ENDS. OTHERWISE REDUCEFMFP2020
C        THE INTERVAL BY CHOOSING ONE END-POINT EQUAL TO X AND REPEAT   FMFP2030
C        THE INTERPOLATION.  WHICH END-POINT IS CHOOSEN DEPENDS ON THE  FMFP2040
C        VALUE OF THE FUNCTION AND ITS GRADIENT AT X                    FMFP2050
C                                                                       FMFP2060
      CALL FUNCT(N,X,F,G)                                               FMFP2070
      IF(F-FX)27,27,28                                                  FMFP2080
   27 IF(F-FY)36,36,28                                                  FMFP2090
   28 DALFA=0.                                                          FMFP2100
      DO 29 I=1,N                                                       FMFP2110
   29 DALFA=DALFA+G(I)*H(I)                                             FMFP2120
      IF(DALFA)30,33,33                                                 FMFP2130
   30 IF(F-FX)32,31,33                                                  FMFP2140
   31 IF(DX-DALFA)32,36,32                                              FMFP2150
   32 FX=F                                                              FMFP2160
      DX=DALFA                                                          FMFP2170
      T=ALFA                                                            FMFP2180
      AMBDA=ALFA                                                        FMFP2190
      GO TO 23                                                          FMFP2200
   33 IF(FY-F)35,34,35                                                  FMFP2210
   34 IF(DY-DALFA)35,36,35                                              FMFP2220
   35 FY=F                                                              FMFP2230
      DY=DALFA                                                          FMFP2240
      AMBDA=AMBDA-ALFA                                                  FMFP2250
      GO TO 22                                                          FMFP2260
C                                                                       FMFP2270
C        TERMINATE, IF FUNCTION HAS NOT DECREASED DURING LAST ITERATION FMFP2280
   36 IF(OLDF-F+EPS)51,38,38                                            FMFP2290
C                                                                       FMFP2300
C        COMPUTE DIFFERENCE VECTORS OF ARGUMENT AND GRADIENT FROM       FMFP2310
C        TWO CONSECUTIVE ITERATIONS                                     FMFP2320
   38 DO 37 J=1,N                                                       FMFP2330
      K=N+J                                                             FMFP2340
      H(K)=G(J)-H(K)                                                    FMFP2350
      K=N+K                                                             FMFP2360
   37 H(K)=X(J)-H(K)                                                    FMFP2370
C                                                                       FMFP2380
C        TEST LENGTH OF ARGUMENT DIFFERENCE VECTOR AND DIRECTION VECTOR FMFP2390
C        IF AT LEAST N ITERATIONS HAVE BEEN EXECUTED. TERMINATE, IF     FMFP2400
C        BOTH ARE LESS THAN  EPS                                        FMFP2410
      IER=0                                                             FMFP2420
      IF(KOUNT-N)42,39,39                                               FMFP2430
   39 T=0.                                                              FMFP2440
      Z=0.                                                              FMFP2450
      DO 40 J=1,N                                                       FMFP2460
      K=N+J                                                             FMFP2470
      W=H(K)                                                            FMFP2480
      K=K+N                                                             FMFP2490
      T=T+ABS(H(K))                                                     FMFP2500
   40 Z=Z+W*H(K)                                                        FMFP2510
      IF(HNRM-EPS)41,41,42                                              FMFP2520
   41 IF(T-EPS)56,56,42                                                 FMFP2530
C                                                                       FMFP2540
C        TERMINATE, IF NUMBER OF ITERATIONS WOULD EXCEED  LIMIT         FMFP2550
   42 IF(KOUNT-LIMIT)43,50,50                                           FMFP2560
C                                                                       FMFP2570
C        PREPARE UPDATING OF MATRIX H                                   FMFP2580
   43 ALFA=0.                                                           FMFP2590
      DO 47 J=1,N                                                       FMFP2600
      K=J+N3                                                            FMFP2610
      W=0.                                                              FMFP2620
      DO 46 L=1,N                                                       FMFP2630
      KL=N+L                                                            FMFP2640
      W=W+H(KL)*H(K)                                                    FMFP2650
      IF(L-J)44,45,45                                                   FMFP2660
   44 K=K+N-L                                                           FMFP2670
      GO TO 46                                                          FMFP2680
   45 K=K+1                                                             FMFP2690
   46 CONTINUE                                                          FMFP2700
      K=N+J                                                             FMFP2710
      ALFA=ALFA+W*H(K)                                                  FMFP2720
   47 H(J)=W                                                            FMFP2730
C                                                                       FMFP2740
C        REPEAT SEARCH IN DIRECTION OF STEEPEST DESCENT IF RESULTS      FMFP2750
C        ARE NOT SATISFACTORY                                           FMFP2760
      IF(Z*ALFA)48,1,48                                                 FMFP2770
C                                                                       FMFP2780
C        UPDATE MATRIX H                                                FMFP2790
   48 K=N31                                                             FMFP2800
      DO 49 L=1,N                                                       FMFP2810
      KL=N2+L                                                           FMFP2820
      DO 49 J=L,N                                                       FMFP2830
      NJ=N2+J                                                           FMFP2840
      H(K)=H(K)+H(KL)*H(NJ)/Z-H(L)*H(J)/ALFA                            FMFP2850
   49 K=K+1                                                             FMFP2860
      GO TO 5                                                           FMFP2870
C        END OF ITERATION LOOP                                          FMFP2880
C                                                                       FMFP2890
C        NO CONVERGENCE AFTER  LIMIT  ITERATIONS                        FMFP2900
   50 IER=1                                                             FMFP2910
      RETURN                                                            FMFP2920
C                                                                       FMFP2930
C        RESTORE OLD VALUES OF FUNCTION AND ARGUMENTS                   FMFP2940
   51 DO 52 J=1,N                                                       FMFP2950
      K=N2+J                                                            FMFP2960
   52 X(J)=H(K)                                                         FMFP2970
      CALL FUNCT(N,X,F,G)                                               FMFP2980
C                                                                       FMFP2990
C        REPEAT SEARCH IN DIRECTION OF STEEPEST DESCENT IF DERIVATIVE   FMFP3000
C        FAILS TO BE SUFFICIENTLY SMALL                                 FMFP3010
      IF(GNRM-EPS)55,55,53                                              FMFP3020
C                                                                       FMFP3030
C        TEST FOR REPEATED FAILURE OF ITERATION                         FMFP3040
   53 IF(IER)56,54,54                                                   FMFP3050
   54 IER=-1                                                            FMFP3060
      GOTO 1                                                            FMFP3070
   55 IER=0                                                             FMFP3080
   56 RETURN                                                            FMFP3090
      END                                                               FMFP3100
