C                                                                       FMCG  10
C     ..................................................................FMCG  20
C                                                                       FMCG  30
C        SUBROUTINE FMCG                                                FMCG  40
C                                                                       FMCG  50
C        PURPOSE                                                        FMCG  60
C           TO FIND A LOCAL MINIMUM OF A FUNCTION OF SEVERAL VARIABLES  FMCG  70
C           BY THE METHOD OF CONJUGATE GRADIENTS                        FMCG  80
C                                                                       FMCG  90
C        USAGE                                                          FMCG 100
C           CALL FMCG(FUNCT,N,X,F,G,EST,EPS,LIMIT,IER,H)                FMCG 110
C                                                                       FMCG 120
C        DESCRIPTION OF PARAMETERS                                      FMCG 130
C           FUNCT  - USER-WRITTEN SUBROUTINE CONCERNING THE FUNCTION TO FMCG 140
C                    BE MINIMIZED. IT MUST BE OF THE FORM               FMCG 150
C                    SUBROUTINE FUNCT(N,ARG,VAL,GRAD)                   FMCG 160
C                    AND MUST SERVE THE FOLLOWING PURPOSE               FMCG 170
C                    FOR EACH N-DIMENSIONAL ARGUMENT VECTOR  ARG,       FMCG 180
C                    FUNCTION VALUE AND GRADIENT VECTOR MUST BE COMPUTEDFMCG 190
C                    AND, ON RETURN, STORED IN VAL AND GRAD RESPECTIVELYFMCG 200
C           N      - NUMBER OF VARIABLES                                FMCG 210
C           X      - VECTOR OF DIMENSION N CONTAINING THE INITIAL       FMCG 220
C                    ARGUMENT WHERE THE ITERATION STARTS. ON RETURN,    FMCG 230
C                    X HOLDS THE ARGUMENT CORRESPONDING TO THE          FMCG 240
C                    COMPUTED MINIMUM FUNCTION VALUE                    FMCG 250
C           F      - SINGLE VARIABLE CONTAINING THE MINIMUM FUNCTION    FMCG 260
C                    VALUE ON RETURN, I.E. F=F(X).                      FMCG 270
C           G      - VECTOR OF DIMENSION N CONTAINING THE GRADIENT      FMCG 280
C                    VECTOR CORRESPONDING TO THE MINIMUM ON RETURN,     FMCG 290
C                    I.E. G=G(X).                                       FMCG 300
C           EST    - IS AN ESTIMATE OF THE MINIMUM FUNCTION VALUE.      FMCG 310
C           EPS    - TESTVALUE REPRESENTING THE EXPECTED ABSOLUTE ERROR.FMCG 320
C                    A REASONABLE CHOICE IS 10**(-6), I.E.              FMCG 330
C                    SOMEWHAT GREATER THAN 10**(-D), WHERE D IS THE     FMCG 340
C                    NUMBER OF SIGNIFICANT DIGITS IN FLOATING POINT     FMCG 350
C                    REPRESENTATION.                                    FMCG 360
C           LIMIT  - MAXIMUM NUMBER OF ITERATIONS.                      FMCG 370
C           IER    - ERROR PARAMETER                                    FMCG 380
C                    IER = 0 MEANS CONVERGENCE WAS OBTAINED             FMCG 390
C                    IER = 1 MEANS NO CONVERGENCE IN LIMIT ITERATIONS   FMCG 400
C                    IER =-1 MEANS ERRORS IN GRADIENT CALCULATION       FMCG 410
C                    IER = 2 MEANS LINEAR SEARCH TECHNIQUE INDICATES    FMCG 420
C                    IT IS LIKELY THAT THERE EXISTS NO MINIMUM.         FMCG 430
C           H      - WORKING STORAGE OF DIMENSION 2*N.                  FMCG 440
C                                                                       FMCG 450
C        REMARKS                                                        FMCG 460
C            I) THE SUBROUTINE NAME REPLACING THE DUMMY ARGUMENT  FUNCT FMCG 470
C               MUST BE DECLARED AS EXTERNAL IN THE CALLING PROGRAM.    FMCG 480
C           II) IER IS SET TO 2 IF , STEPPING IN ONE OF THE COMPUTED    FMCG 490
C               DIRECTIONS, THE FUNCTION WILL NEVER INCREASE WITHIN     FMCG 500
C               A TOLERABLE RANGE OF ARGUMENT.                          FMCG 510
C               IER = 2 MAY OCCUR ALSO IF THE INTERVAL WHERE F          FMCG 520
C               INCREASES IS SMALL AND THE INITIAL ARGUMENT WAS         FMCG 530
C               RELATIVELY FAR AWAY FROM THE MINIMUM SUCH THAT THE      FMCG 540
C               MINIMUM WAS OVERLEAPED. THIS IS DUE TO THE SEARCH       FMCG 550
C               TECHNIQUE WHICH DOUBLES THE STEPSIZE UNTIL A POINT      FMCG 560
C               IS FOUND WHERE THE FUNCTION INCREASES.                  FMCG 570
C                                                                       FMCG 580
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  FMCG 590
C           FUNCT                                                       FMCG 600
C                                                                       FMCG 610
C        METHOD                                                         FMCG 620
C           THE METHOD IS DESCRIBED IN THE FOLLOWING ARTICLE            FMCG 630
C           R.FLETCHER AND C.M.REEVES, FUNCTION MINIMIZATION BY         FMCG 640
C           CONJUGATE GRADIENTS,                                        FMCG 650
C           COMPUTER JOURNAL VOL.7, ISS.2, 1964, PP.149-154.            FMCG 660
C                                                                       FMCG 670
C     ..................................................................FMCG 680
C                                                                       FMCG 690
      SUBROUTINE FMCG(FUNCT,N,X,F,G,EST,EPS,LIMIT,IER,H)                FMCG 700
C                                                                       FMCG 710
C        DIMENSIONED DUMMY VARIABLES                                    FMCG 720
      DIMENSION X(1),G(1),H(1)                                          FMCG 730
C                                                                       FMCG 740
C                                                                       FMCG 750
C        COMPUTE FUNCTION VALUE AND GRADIENT VECTOR FOR INITIAL ARGUMENTFMCG 760
      CALL FUNCT(N,X,F,G)                                               FMCG 770
C                                                                       FMCG 780
C        RESET ITERATION COUNTER                                        FMCG 790
      KOUNT=0                                                           FMCG 800
      IER=0                                                             FMCG 810
      N1=N+1                                                            FMCG 820
C                                                                       FMCG 830
C        START ITERATION CYCLE FOR EVERY N+1 ITERATIONS                 FMCG 840
    1 DO 43 II=1,N1                                                     FMCG 850
C                                                                       FMCG 860
C        STEP ITERATION COUNTER AND SAVE FUNCTION VALUE                 FMCG 870
      KOUNT=KOUNT+1                                                     FMCG 880
      OLDF=F                                                            FMCG 890
C                                                                       FMCG 900
C        COMPUTE SQUARE OF GRADIENT AND TERMINATE IF ZERO               FMCG 910
      GNRM=0.                                                           FMCG 920
      DO 2 J=1,N                                                        FMCG 930
    2 GNRM=GNRM+G(J)*G(J)                                               FMCG 940
      IF(GNRM)46,46,3                                                   FMCG 950
C                                                                       FMCG 960
C        EACH TIME THE ITERATION LOOP IS EXECUTED , THE FIRST STEP WILL FMCG 970
C        BE IN DIRECTION OF STEEPEST DESCENT                            FMCG 980
    3 IF(II-1)4,4,6                                                     FMCG 990
    4 DO 5 J=1,N                                                        FMCG1000
    5 H(J)=-G(J)                                                        FMCG1010
      GO TO 8                                                           FMCG1020
C                                                                       FMCG1030
C        FURTHER DIRECTION VECTORS H WILL BE CHOOSEN CORRESPONDING      FMCG1040
C        TO THE CONJUGATE GRADIENT METHOD                               FMCG1050
    6 AMBDA=GNRM/OLDG                                                   FMCG1060
      DO 7 J=1,N                                                        FMCG1070
    7 H(J)=AMBDA*H(J)-G(J)                                              FMCG1080
C                                                                       FMCG1090
C        COMPUTE TESTVALUE FOR DIRECTIONAL VECTOR AND DIRECTIONAL       FMCG1100
C        DERIVATIVE                                                     FMCG1110
    8 DY=0.                                                             FMCG1120
      HNRM=0.                                                           FMCG1130
      DO 9 J=1,N                                                        FMCG1140
      K=J+N                                                             FMCG1150
C                                                                       FMCG1160
C        SAVE ARGUMENT VECTOR                                           FMCG1170
      H(K)=X(J)                                                         FMCG1180
      HNRM=HNRM+ABS(H(J))                                               FMCG1190
    9 DY=DY+H(J)*G(J)                                                   FMCG1200
C                                                                       FMCG1210
C        CHECK WHETHER FUNCTION WILL DECREASE STEPPING ALONG H AND      FMCG1220
C        SKIP LINEAR SEARCH ROUTINE IF NOT                              FMCG1230
      IF(DY)10,42,42                                                    FMCG1240
C                                                                       FMCG1250
C        COMPUTE SCALE FACTOR USED IN LINEAR SEARCH SUBROUTINE          FMCG1260
   10 SNRM=1./HNRM                                                      FMCG1270
C                                                                       FMCG1280
C        SEARCH MINIMUM ALONG DIRECTION H                               FMCG1290
C                                                                       FMCG1300
C        SEARCH ALONG H FOR POSITIVE DIRECTIONAL DERIVATIVE             FMCG1310
      FY=F                                                              FMCG1320
      ALFA=2.*(EST-F)/DY                                                FMCG1330
      AMBDA=SNRM                                                        FMCG1340
C                                                                       FMCG1350
C        USE ESTIMATE FOR STEPSIZE ONLY IF IT IS POSITIVE AND LESS THAN FMCG1360
C        SNRM. OTHERWISE TAKE SNRM AS STEPSIZE.                         FMCG1370
      IF(ALFA)13,13,11                                                  FMCG1380
   11 IF(ALFA-AMBDA)12,13,13                                            FMCG1390
   12 AMBDA=ALFA                                                        FMCG1400
   13 ALFA=0.                                                           FMCG1410
C                                                                       FMCG1420
C        SAVE FUNCTION AND DERIVATIVE VALUES FOR OLD ARGUMENT           FMCG1430
   14 FX=FY                                                             FMCG1440
      DX=DY                                                             FMCG1450
C                                                                       FMCG1460
C        STEP ARGUMENT ALONG H                                          FMCG1470
      DO 15 I=1,N                                                       FMCG1480
   15 X(I)=X(I)+AMBDA*H(I)                                              FMCG1490
C                                                                       FMCG1500
C        COMPUTE FUNCTION VALUE AND GRADIENT FOR NEW ARGUMENT           FMCG1510
      CALL FUNCT(N,X,F,G)                                               FMCG1520
      FY=F                                                              FMCG1530
C                                                                       FMCG1540
C        COMPUTE DIRECTIONAL DERIVATIVE DY FOR NEW ARGUMENT.  TERMINATE FMCG1550
C        SEARCH, IF DY POSITIVE. IF DY IS ZERO THE MINIMUM IS FOUND     FMCG1560
      DY=0.                                                             FMCG1570
      DO 16 I=1,N                                                       FMCG1580
   16 DY=DY+G(I)*H(I)                                                   FMCG1590
      IF(DY)17,38,20                                                    FMCG1600
C                                                                       FMCG1610
C        TERMINATE SEARCH ALSO IF THE FUNCTION VALUE INDICATES THAT     FMCG1620
C        A MINIMUM HAS BEEN PASSED                                      FMCG1630
   17 IF(FY-FX)18,20,20                                                 FMCG1640
C                                                                       FMCG1650
C        REPEAT SEARCH AND DOUBLE STEPSIZE FOR FURTHER SEARCHES         FMCG1660
   18 AMBDA=AMBDA+ALFA                                                  FMCG1670
      ALFA=AMBDA                                                        FMCG1680
C                                                                       FMCG1690
C        TERMINATE IF THE CHANGE IN ARGUMENT GETS VERY LARGE            FMCG1700
      IF(HNRM*AMBDA-1.E10)14,14,19                                      FMCG1710
C                                                                       FMCG1720
C        LINEAR SEARCH TECHNIQUE INDICATES THAT NO MINIMUM EXISTS       FMCG1730
   19 IER=2                                                             FMCG1740
C                                                                       FMCG1741
C        RESTORE OLD VALUES OF FUNCTION AND ARGUMENTS                   FMCG1742
      F=OLDF                                                            FMCG1743
      DO 100 J=1,N                                                      FMCG1744
      G(J)=H(J)                                                         FMCG1745
      K=N+J                                                             FMCG1746
  100 X(J)=H(K)                                                         FMCG1747
      RETURN                                                            FMCG1750
C        END OF SEARCH LOOP                                             FMCG1760
C                                                                       FMCG1770
C        INTERPOLATE CUBICALLY IN THE INTERVAL DEFINED BY THE SEARCH    FMCG1780
C        ABOVE AND COMPUTE THE ARGUMENT X FOR WHICH THE INTERPOLATION   FMCG1790
C        POLYNOMIAL IS MINIMIZED                                        FMCG1800
C                                                                       FMCG1810
   20 T=0.                                                              FMCG1820
   21 IF(AMBDA)22,38,22                                                 FMCG1830
   22 Z=3.*(FX-FY)/AMBDA+DX+DY                                          FMCG1840
      ALFA=AMAX1(ABS(Z),ABS(DX),ABS(DY))                                FMCG1850
      DALFA=Z/ALFA                                                      FMCG1860
      DALFA=DALFA*DALFA-DX/ALFA*DY/ALFA                                 FMCG1870
      IF(DALFA)23,27,27                                                 FMCG1880
C                                                                       FMCG1890
C        RESTORE OLD VALUES OF FUNCTION AND ARGUMENTS                   FMCG1900
   23 DO 24 J=1,N                                                       FMCG1910
      K=N+J                                                             FMCG1920
   24 X(J)=H(K)                                                         FMCG1930
      CALL FUNCT(N,X,F,G)                                               FMCG1940
C                                                                       FMCG1950
C        TEST FOR REPEATED FAILURE OF ITERATION                         FMCG1960
   25 IF(IER)47,26,47                                                   FMCG1970
   26 IER=-1                                                            FMCG1980
      GOTO 1                                                            FMCG1990
   27 W=ALFA*SQRT(DALFA)                                                FMCG2000
      ALFA=DY-DX+W+W                                                    FMCG2010
      IF(ALFA)270,271,270                                               FMCG2011
  270 ALFA=(DY-Z+W)/ALFA                                                FMCG2012
      GO TO 272                                                         FMCG2013
  271 ALFA=(Z+DY-W)/(Z+DX+Z+DY)                                         FMCG2014
  272 ALFA=ALFA*AMBDA                                                   FMCG2015
      DO 28 I=1,N                                                       FMCG2020
   28 X(I)=X(I)+(T-ALFA)*H(I)                                           FMCG2030
C                                                                       FMCG2040
C        TERMINATE, IF THE VALUE OF THE ACTUAL FUNCTION AT X IS LESS    FMCG2050
C        THAN THE FUNCTION VALUES AT THE INTERVAL ENDS. OTHERWISE REDUCEFMCG2060
C        THE INTERVAL BY CHOOSING ONE END-POINT EQUAL TO X AND REPEAT   FMCG2070
C        THE INTERPOLATION.  WHICH END-POINT IS CHOOSEN DEPENDS ON THE  FMCG2080
C        VALUE OF THE FUNCTION AND ITS GRADIENT AT X                    FMCG2090
C                                                                       FMCG2100
      CALL FUNCT(N,X,F,G)                                               FMCG2110
      IF(F-FX)29,29,30                                                  FMCG2120
   29 IF(F-FY)38,38,30                                                  FMCG2130
C                                                                       FMCG2140
C        COMPUTE DIRECTIONAL DERIVATIVE                                 FMCG2150
   30 DALFA=0.                                                          FMCG2160
      DO 31 I=1,N                                                       FMCG2170
   31 DALFA=DALFA+G(I)*H(I)                                             FMCG2180
      IF(DALFA)32,35,35                                                 FMCG2190
   32 IF(F-FX)34,33,35                                                  FMCG2200
   33 IF(DX-DALFA)34,38,34                                              FMCG2210
   34 FX=F                                                              FMCG2220
      DX=DALFA                                                          FMCG2230
      T=ALFA                                                            FMCG2240
      AMBDA=ALFA                                                        FMCG2250
      GO TO 21                                                          FMCG2260
   35 IF(FY-F)37,36,37                                                  FMCG2270
   36 IF(DY-DALFA)37,38,37                                              FMCG2280
   37 FY=F                                                              FMCG2290
      DY=DALFA                                                          FMCG2300
      AMBDA=AMBDA-ALFA                                                  FMCG2310
      GO TO 20                                                          FMCG2320
C                                                                       FMCG2330
C        TERMINATE, IF FUNCTION HAS NOT DECREASED DURING LAST ITERATION FMCG2340
C        OTHERWISE SAVE GRADIENT NORM                                   FMCG2350
   38 IF(OLDF-F+EPS)19,25,39                                            FMCG2360
   39 OLDG=GNRM                                                         FMCG2370
C                                                                       FMCG2380
C        COMPUTE DIFFERENCE OF NEW AND OLD ARGUMENT VECTOR              FMCG2390
      T=0.                                                              FMCG2400
      DO 40 J=1,N                                                       FMCG2410
      K=J+N                                                             FMCG2420
      H(K)=X(J)-H(K)                                                    FMCG2430
   40 T=T+ABS(H(K))                                                     FMCG2440
C                                                                       FMCG2450
C        TEST LENGTH OF DIFFERENCE VECTOR IF AT LEAST N+1 ITERATIONS    FMCG2460
C        HAVE BEEN EXECUTED. TERMINATE, IF LENGTH IS LESS THAN EPS      FMCG2470
      IF(KOUNT-N1)42,41,41                                              FMCG2480
   41 IF(T-EPS)45,45,42                                                 FMCG2490
C                                                                       FMCG2500
C        TERMINATE, IF NUMBER OF ITERATIONS WOULD EXCEED  LIMIT         FMCG2510
   42 IF(KOUNT-LIMIT)43,44,44                                           FMCG2520
   43 IER=0                                                             FMCG2530
C        END OF ITERATION CYCLE                                         FMCG2540
C                                                                       FMCG2550
C        START NEXT ITERATION CYCLE                                     FMCG2560
      GO TO 1                                                           FMCG2570
C                                                                       FMCG2580
C        NO CONVERGENCE AFTER  LIMIT  ITERATIONS                        FMCG2590
   44 IER=1                                                             FMCG2600
      IF(GNRM-EPS)46,46,47                                              FMCG2610
C                                                                       FMCG2620
C        TEST FOR SUFFICIENTLY SMALL GRADIENT                           FMCG2630
   45 IF(GNRM-EPS)46,46,25                                              FMCG2640
   46 IER=0                                                             FMCG2650
   47 RETURN                                                            FMCG2660
      END                                                               FMCG2670
