C                                                                       DFMC  10
C     ..................................................................DFMC  20
C                                                                       DFMC  30
C        SUBROUTINE DFMCG                                               DFMC  40
C                                                                       DFMC  50
C        PURPOSE                                                        DFMC  60
C           TO FIND A LOCAL MINIMUM OF A FUNCTION OF SEVERAL VARIABLES  DFMC  70
C           BY THE METHOD OF CONJUGATE GRADIENTS                        DFMC  80
C                                                                       DFMC  90
C        USAGE                                                          DFMC 100
C           CALL DFMCG(FUNCT,N,X,F,G,EST,EPS,LIMIT,IER,H)               DFMC 110
C                                                                       DFMC 120
C        DESCRIPTION OF PARAMETERS                                      DFMC 130
C           FUNCT  - USER-WRITTEN SUBROUTINE CONCERNING THE FUNCTION TO DFMC 140
C                    BE MINIMIZED. IT MUST BE OF THE FORM               DFMC 150
C                    SUBROUTINE FUNCT(N,ARG,VAL,GRAD)                   DFMC 160
C                    AND MUST SERVE THE FOLLOWING PURPOSE               DFMC 170
C                    FOR EACH N-DIMENSIONAL ARGUMENT VECTOR  ARG,       DFMC 180
C                    FUNCTION VALUE AND GRADIENT VECTOR MUST BE COMPUTEDDFMC 190
C                    AND, ON RETURN, STORED IN VAL AND GRAD RESPECTIVELYDFMC 200
C                    ARG,VAL AND GRAD MUST BE OF DOUBLE PRECISION.      DFMC 210
C           N      - NUMBER OF VARIABLES                                DFMC 220
C           X      - VECTOR OF DIMENSION N CONTAINING THE INITIAL       DFMC 230
C                    ARGUMENT WHERE THE ITERATION STARTS. ON RETURN,    DFMC 240
C                    X HOLDS THE ARGUMENT CORRESPONDING TO THE          DFMC 250
C                    COMPUTED MINIMUM FUNCTION VALUE                    DFMC 260
C                    DOUBLE PRECISION VECTOR.                           DFMC 270
C           F      - SINGLE VARIABLE CONTAINING THE MINIMUM FUNCTION    DFMC 280
C                    VALUE ON RETURN, I.E. F=F(X).                      DFMC 290
C                    DOUBLE PRECISION VARIABLE.                         DFMC 300
C           G      - VECTOR OF DIMENSION N CONTAINING THE GRADIENT      DFMC 310
C                    VECTOR CORRESPONDING TO THE MINIMUM ON RETURN,     DFMC 320
C                    I.E. G=G(X).                                       DFMC 330
C                    DOUBLE PRECISION VECTOR.                           DFMC 340
C           EST    - IS AN ESTIMATE OF THE MINIMUM FUNCTION VALUE.      DFMC 350
C                    SINGLE PRECISION VARIABLE.                         DFMC 360
C           EPS    - TESTVALUE REPRESENTING THE EXPECTED ABSOLUTE ERROR.DFMC 370
C                    A REASONABLE CHOICE IS 10**(-16), I.E.             DFMC 380
C                    SOMEWHAT GREATER THAN 10**(-D), WHERE D IS THE     DFMC 390
C                    NUMBER OF SIGNIFICANT DIGITS IN FLOATING POINT     DFMC 400
C                    REPRESENTATION.                                    DFMC 410
C                    SINGLE PRECISION VARIABLE.                         DFMC 420
C           LIMIT  - MAXIMUM NUMBER OF ITERATIONS.                      DFMC 430
C           IER    - ERROR PARAMETER                                    DFMC 440
C                    IER = 0 MEANS CONVERGENCE WAS OBTAINED             DFMC 450
C                    IER = 1 MEANS NO CONVERGENCE IN LIMIT ITERATIONS   DFMC 460
C                    IER =-1 MEANS ERRORS IN GRADIENT CALCULATION       DFMC 470
C                    IER = 2 MEANS LINEAR SEARCH TECHNIQUE INDICATES    DFMC 480
C                    IT IS LIKELY THAT THERE EXISTS NO MINIMUM.         DFMC 490
C           H      - WORKING STORAGE OF DIMENSION 2*N.                  DFMC 500
C                    DOUBLE PRECISION ARRAY.                            DFMC 510
C                                                                       DFMC 520
C        REMARKS                                                        DFMC 530
C            I) THE SUBROUTINE NAME REPLACING THE DUMMY ARGUMENT  FUNCT DFMC 540
C               MUST BE DECLARED AS EXTERNAL IN THE CALLING PROGRAM.    DFMC 550
C           II) IER IS SET TO 2 IF , STEPPING IN ONE OF THE COMPUTED    DFMC 560
C               DIRECTIONS, THE FUNCTION WILL NEVER INCREASE WITHIN     DFMC 570
C               A TOLERABLE RANGE OF ARGUMENT.                          DFMC 580
C               IER = 2 MAY OCCUR ALSO IF THE INTERVAL WHERE F          DFMC 590
C               INCREASES IS SMALL AND THE INITIAL ARGUMENT WAS         DFMC 600
C               RELATIVELY FAR AWAY FROM THE MINIMUM SUCH THAT THE      DFMC 610
C               MINIMUM WAS OVERLEAPED. THIS IS DUE TO THE SEARCH       DFMC 620
C               TECHNIQUE WHICH DOUBLES THE STEPSIZE UNTIL A POINT      DFMC 630
C               IS FOUND WHERE THE FUNCTION INCREASES.                  DFMC 640
C                                                                       DFMC 650
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DFMC 660
C           FUNCT                                                       DFMC 670
C                                                                       DFMC 680
C        METHOD                                                         DFMC 690
C           THE METHOD IS DESCRIBED IN THE FOLLOWING ARTICLE            DFMC 700
C           R.FLETCHER AND C.M.REEVES, FUNCTION MINIMIZATION BY         DFMC 710
C           CONJUGATE GRADIENTS,                                        DFMC 720
C           COMPUTER JOURNAL VOL.7, ISS.2, 1964, PP.149-154.            DFMC 730
C                                                                       DFMC 740
C     ..................................................................DFMC 750
C                                                                       DFMC 760
      SUBROUTINE DFMCG(FUNCT,N,X,F,G,EST,EPS,LIMIT,IER,H)               DFMC 770
C                                                                       DFMC 780
C        DIMENSIONED DUMMY VARIABLES                                    DFMC 790
      DIMENSION X(1),G(1),H(1)                                          DFMC 800
      DOUBLE PRECISION X,G,GNRM,H,HNRM,F,FX,FY,OLDF,OLDG,SNRM,AMBDA,    DFMC 810
     1ALFA,DALFA,T,Z,W,DX,DY                                            DFMC 820
C                                                                       DFMC 830
C        COMPUTE FUNCTION VALUE AND GRADIENT VECTOR FOR INITIAL ARGUMENTDFMC 840
      CALL FUNCT(N,X,F,G)                                               DFMC 850
C                                                                       DFMC 860
C        RESET ITERATION COUNTER                                        DFMC 870
      KOUNT=0                                                           DFMC 880
      IER=0                                                             DFMC 890
      N1=N+1                                                            DFMC 900
C                                                                       DFMC 910
C        START ITERATION CYCLE FOR EVERY N+1 ITERATIONS                 DFMC 920
    1 DO 43 II=1,N1                                                     DFMC 930
C                                                                       DFMC 940
C        STEP ITERATION COUNTER AND SAVE FUNCTION VALUE                 DFMC 950
      KOUNT=KOUNT+1                                                     DFMC 960
      OLDF=F                                                            DFMC 970
C                                                                       DFMC 980
C        COMPUTE SQUARE OF GRADIENT AND TERMINATE IF ZERO               DFMC 990
      GNRM=0.D0                                                         DFMC1000
      DO 2 J=1,N                                                        DFMC1010
    2 GNRM=GNRM+G(J)*G(J)                                               DFMC1020
      IF(GNRM)46,46,3                                                   DFMC1030
C                                                                       DFMC1040
C        EACH TIME THE ITERATION LOOP IS EXECUTED , THE FIRST STEP WILL DFMC1050
C        BE IN DIRECTION OF STEEPEST DESCENT                            DFMC1060
    3 IF(II-1)4,4,6                                                     DFMC1070
    4 DO 5 J=1,N                                                        DFMC1080
    5 H(J)=-G(J)                                                        DFMC1090
      GO TO 8                                                           DFMC1100
C                                                                       DFMC1110
C        FURTHER DIRECTION VECTORS H WILL BE CHOOSEN CORRESPONDING      DFMC1120
C        TO THE CONJUGATE GRADIENT METHOD                               DFMC1130
    6 AMBDA=GNRM/OLDG                                                   DFMC1140
      DO 7 J=1,N                                                        DFMC1150
    7 H(J)=AMBDA*H(J)-G(J)                                              DFMC1160
C                                                                       DFMC1170
C        COMPUTE TESTVALUE FOR DIRECTIONAL VECTOR AND DIRECTIONAL       DFMC1180
C        DERIVATIVE                                                     DFMC1190
    8 DY=0.D0                                                           DFMC1200
      HNRM=0.D0                                                         DFMC1210
      DO 9 J=1,N                                                        DFMC1220
      K=J+N                                                             DFMC1230
C                                                                       DFMC1240
C        SAVE ARGUMENT VECTOR                                           DFMC1250
      H(K)=X(J)                                                         DFMC1260
      HNRM=HNRM+DABS(H(J))                                              DFMC1270
    9 DY=DY+H(J)*G(J)                                                   DFMC1280
C                                                                       DFMC1290
C        CHECK WHETHER FUNCTION WILL DECREASE STEPPING ALONG H AND      DFMC1300
C        SKIP LINEAR SEARCH ROUTINE IF NOT                              DFMC1310
      IF(DY)10,42,42                                                    DFMC1320
C                                                                       DFMC1330
C        COMPUTE SCALE FACTOR USED IN LINEAR SEARCH SUBROUTINE          DFMC1340
   10 SNRM=1.D0/HNRM                                                    DFMC1350
C                                                                       DFMC1360
C        SEARCH MINIMUM ALONG DIRECTION H                               DFMC1370
C                                                                       DFMC1380
C        SEARCH ALONG H FOR POSITIVE DIRECTIONAL DERIVATIVE             DFMC1390
      FY=F                                                              DFMC1400
      ALFA=2.D0*(EST-F)/DY                                              DFMC1410
      AMBDA=SNRM                                                        DFMC1420
C                                                                       DFMC1430
C        USE ESTIMATE FOR STEPSIZE ONLY IF IT IS POSITIVE AND LESS THAN DFMC1440
C        SNRM. OTHERWISE TAKE SNRM AS STEPSIZE.                         DFMC1450
      IF(ALFA)13,13,11                                                  DFMC1460
   11 IF(ALFA-AMBDA)12,13,13                                            DFMC1470
   12 AMBDA=ALFA                                                        DFMC1480
   13 ALFA=0.D0                                                         DFMC1490
C                                                                       DFMC1500
C        SAVE FUNCTION AND DERIVATIVE VALUES FOR OLD ARGUMENT           DFMC1510
   14 FX=FY                                                             DFMC1520
      DX=DY                                                             DFMC1530
C                                                                       DFMC1540
C        STEP ARGUMENT ALONG H                                          DFMC1550
      DO 15 I=1,N                                                       DFMC1560
   15 X(I)=X(I)+AMBDA*H(I)                                              DFMC1570
C                                                                       DFMC1580
C        COMPUTE FUNCTION VALUE AND GRADIENT FOR NEW ARGUMENT           DFMC1590
      CALL FUNCT(N,X,F,G)                                               DFMC1600
      FY=F                                                              DFMC1610
C                                                                       DFMC1620
C        COMPUTE DIRECTIONAL DERIVATIVE DY FOR NEW ARGUMENT.  TERMINATE DFMC1630
C        SEARCH, IF DY POSITIVE. IF DY IS ZERO THE MINIMUM IS FOUND     DFMC1640
      DY=0.D0                                                           DFMC1650
      DO 16 I=1,N                                                       DFMC1660
   16 DY=DY+G(I)*H(I)                                                   DFMC1670
      IF(DY)17,38,20                                                    DFMC1680
C                                                                       DFMC1690
C        TERMINATE SEARCH ALSO IF THE FUNCTION VALUE INDICATES THAT     DFMC1700
C        A MINIMUM HAS BEEN PASSED                                      DFMC1710
   17 IF(FY-FX)18,20,20                                                 DFMC1720
C                                                                       DFMC1730
C        REPEAT SEARCH AND DOUBLE STEPSIZE FOR FURTHER SEARCHES         DFMC1740
   18 AMBDA=AMBDA+ALFA                                                  DFMC1750
      ALFA=AMBDA                                                        DFMC1760
C                                                                       DFMC1770
C        TERMINATE IF THE CHANGE IN ARGUMENT GETS VERY LARGE            DFMC1780
      IF(HNRM*AMBDA-1.D10)14,14,19                                      DFMC1790
C                                                                       DFMC1800
C        LINEAR SEARCH TECHNIQUE INDICATES THAT NO MINIMUM EXISTS       DFMC1810
   19 IER=2                                                             DFMC1820
C                                                                       DFMC1821
C        RESTORE OLD VALUES OF FUNCTION AND ARGUMENTS                   DFMC1822
      F=OLDF                                                            DFMC1823
      DO 100 J=1,N                                                      DFMC1824
      G(J)=H(J)                                                         DFMC1825
      K=N+J                                                             DFMC1826
  100 X(J)=H(K)                                                         DFMC1827
      RETURN                                                            DFMC1830
C        END OF SEARCH LOOP                                             DFMC1840
C                                                                       DFMC1850
C        INTERPOLATE CUBICALLY IN THE INTERVAL DEFINED BY THE SEARCH    DFMC1860
C        ABOVE AND COMPUTE THE ARGUMENT X FOR WHICH THE INTERPOLATION   DFMC1870
C        POLYNOMIAL IS MINIMIZED                                        DFMC1880
C                                                                       DFMC1890
   20 T=0.                                                              DFMC1900
   21 IF(AMBDA)22,38,22                                                 DFMC1910
   22 Z=3.D0*(FX-FY)/AMBDA+DX+DY                                        DFMC1920
      ALFA=DMAX1(DABS(Z),DABS(DX),DABS(DY))                             DFMC1930
      DALFA=Z/ALFA                                                      DFMC1940
      DALFA=DALFA*DALFA-DX/ALFA*DY/ALFA                                 DFMC1950
      IF(DALFA)23,27,27                                                 DFMC1960
C                                                                       DFMC1970
C        RESTORE OLD VALUES OF FUNCTION AND ARGUMENTS                   DFMC1980
   23 DO 24 J=1,N                                                       DFMC1990
      K=N+J                                                             DFMC2000
   24 X(J)=H(K)                                                         DFMC2010
      CALL FUNCT(N,X,F,G)                                               DFMC2020
C                                                                       DFMC2030
C        TEST FOR REPEATED FAILURE OF ITERATION                         DFMC2040
   25 IF(IER)47,26,47                                                   DFMC2050
   26 IER=-1                                                            DFMC2060
      GOTO 1                                                            DFMC2070
   27 W=ALFA*DSQRT(DALFA)                                               DFMC2080
      ALFA=DY-DX+W+W                                                    DFMC2090
      IF(ALFA)270,271,270                                               DFMC2091
  270 ALFA=(DY-Z+W)/ALFA                                                DFMC2092
      GO TO 272                                                         DFMC2093
  271 ALFA=(Z+DY-W)/(Z+DX+Z+DY)                                         DFMC2094
  272 ALFA=ALFA*AMBDA                                                   DFMC2095
      DO 28 I=1,N                                                       DFMC2100
   28 X(I)=X(I)+(T-ALFA)*H(I)                                           DFMC2110
C                                                                       DFMC2120
C        TERMINATE, IF THE VALUE OF THE ACTUAL FUNCTION AT X IS LESS    DFMC2130
C        THAN THE FUNCTION VALUES AT THE INTERVAL ENDS. OTHERWISE REDUCEDFMC2140
C        THE INTERVAL BY CHOOSING ONE END-POINT EQUAL TO X AND REPEAT   DFMC2150
C        THE INTERPOLATION.  WHICH END-POINT IS CHOOSEN DEPENDS ON THE  DFMC2160
C        VALUE OF THE FUNCTION AND ITS GRADIENT AT X                    DFMC2170
C                                                                       DFMC2180
      CALL FUNCT(N,X,F,G)                                               DFMC2190
      IF(F-FX)29,29,30                                                  DFMC2200
   29 IF(F-FY)38,38,30                                                  DFMC2210
C                                                                       DFMC2220
C        COMPUTE DIRECTIONAL DERIVATIVE                                 DFMC2230
   30 DALFA=0.D0                                                        DFMC2240
      DO 31 I=1,N                                                       DFMC2250
   31 DALFA=DALFA+G(I)*H(I)                                             DFMC2260
      IF(DALFA)32,35,35                                                 DFMC2270
   32 IF(F-FX)34,33,35                                                  DFMC2280
   33 IF(DX-DALFA)34,38,34                                              DFMC2290
   34 FX=F                                                              DFMC2300
      DX=DALFA                                                          DFMC2310
      T=ALFA                                                            DFMC2320
      AMBDA=ALFA                                                        DFMC2330
      GO TO 21                                                          DFMC2340
   35 IF(FY-F)37,36,37                                                  DFMC2350
   36 IF(DY-DALFA)37,38,37                                              DFMC2360
   37 FY=F                                                              DFMC2370
      DY=DALFA                                                          DFMC2380
      AMBDA=AMBDA-ALFA                                                  DFMC2390
      GO TO 20                                                          DFMC2400
C                                                                       DFMC2410
C        TERMINATE, IF FUNCTION HAS NOT DECREASED DURING LAST ITERATION DFMC2420
C        OTHERWISE SAVE GRADIENT NORM                                   DFMC2430
   38 IF(OLDF-F+EPS)19,25,39                                            DFMC2440
   39 OLDG=GNRM                                                         DFMC2450
C                                                                       DFMC2460
C        COMPUTE DIFFERENCE OF NEW AND OLD ARGUMENT VECTOR              DFMC2470
      T=0.D0                                                            DFMC2480
      DO 40 J=1,N                                                       DFMC2490
      K=J+N                                                             DFMC2500
      H(K)=X(J)-H(K)                                                    DFMC2510
   40 T=T+DABS(H(K))                                                    DFMC2520
C                                                                       DFMC2530
C        TEST LENGTH OF DIFFERENCE VECTOR IF AT LEAST N+1 ITERATIONS    DFMC2540
C        HAVE BEEN EXECUTED. TERMINATE, IF LENGTH IS LESS THAN EPS      DFMC2550
      IF(KOUNT-N1)42,41,41                                              DFMC2560
   41 IF(T-EPS)45,45,42                                                 DFMC2570
C                                                                       DFMC2580
C        TERMINATE, IF NUMBER OF ITERATIONS WOULD EXCEED  LIMIT         DFMC2590
   42 IF(KOUNT-LIMIT)43,44,44                                           DFMC2600
   43 IER=0                                                             DFMC2610
C        END OF ITERATION CYCLE                                         DFMC2620
C                                                                       DFMC2630
C        START NEXT ITERATION CYCLE                                     DFMC2640
      GO TO 1                                                           DFMC2650
C                                                                       DFMC2660
C        NO CONVERGENCE AFTER  LIMIT  ITERATIONS                        DFMC2670
   44 IER=1                                                             DFMC2680
      IF(GNRM-EPS)46,46,47                                              DFMC2690
C                                                                       DFMC2700
C        TEST FOR SUFFICIENTLY SMALL GRADIENT                           DFMC2710
   45 IF(GNRM-EPS)46,46,25                                              DFMC2720
   46 IER=0                                                             DFMC2730
   47 RETURN                                                            DFMC2740
      END                                                               DFMC2750
