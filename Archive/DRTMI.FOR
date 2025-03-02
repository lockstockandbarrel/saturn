C                                                                       DRTM  10
C     ..................................................................DRTM  20
C                                                                       DRTM  30
C        SUBROUTINE DRTMI                                               DRTM  40
C                                                                       DRTM  50
C        PURPOSE                                                        DRTM  60
C           TO SOLVE GENERAL NONLINEAR EQUATIONS OF THE FORM FCT(X)=0   DRTM  70
C           BY MEANS OF MUELLER-S ITERATION METHOD.                     DRTM  80
C                                                                       DRTM  90
C        USAGE                                                          DRTM 100
C           CALL DRTMI (X,F,FCT,XLI,XRI,EPS,IEND,IER)                   DRTM 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT.               DRTM 120
C                                                                       DRTM 130
C        DESCRIPTION OF PARAMETERS                                      DRTM 140
C           X      - DOUBLE PRECISION RESULTANT ROOT OF EQUATION        DRTM 150
C                    FCT(X)=0.                                          DRTM 160
C           F      - DOUBLE PRECISION RESULTANT FUNCTION VALUE          DRTM 170
C                    AT ROOT X.                                         DRTM 180
C           FCT    - NAME OF THE EXTERNAL DOUBLE PRECISION FUNCTION     DRTM 190
C                    SUBPROGRAM USED.                                   DRTM 200
C           XLI    - DOUBLE PRECISION INPUT VALUE WHICH SPECIFIES THE   DRTM 210
C                    INITIAL LEFT BOUND OF THE ROOT X.                  DRTM 220
C           XRI    - DOUBLE PRECISION INPUT VALUE WHICH SPECIFIES THE   DRTM 230
C                    INITIAL RIGHT BOUND OF THE ROOT X.                 DRTM 240
C           EPS    - SINGLE PRECISION INPUT VALUE WHICH SPECIFIES THE   DRTM 250
C                    UPPER BOUND OF THE ERROR OF RESULT X.              DRTM 260
C           IEND   - MAXIMUM NUMBER OF ITERATION STEPS SPECIFIED.       DRTM 270
C           IER    - RESULTANT ERROR PARAMETER CODED AS FOLLOWS         DRTM 280
C                     IER=0 - NO ERROR,                                 DRTM 290
C                     IER=1 - NO CONVERGENCE AFTER IEND ITERATION STEPS DRTM 300
C                             FOLLOWED BY IEND SUCCESSIVE STEPS OF      DRTM 310
C                             BISECTION,                                DRTM 320
C                     IER=2 - BASIC ASSUMPTION FCT(XLI)*FCT(XRI) LESS   DRTM 330
C                             THAN OR EQUAL TO ZERO IS NOT SATISFIED.   DRTM 340
C                                                                       DRTM 350
C        REMARKS                                                        DRTM 360
C           THE PROCEDURE ASSUMES THAT FUNCTION VALUES AT INITIAL       DRTM 370
C           BOUNDS XLI AND XRI HAVE NOT THE SAME SIGN. IF THIS BASIC    DRTM 380
C           ASSUMPTION IS NOT SATISFIED BY INPUT VALUES XLI AND XRI, THEDRTM 390
C           PROCEDURE IS BYPASSED AND GIVES THE ERROR MESSAGE IER=2.    DRTM 400
C                                                                       DRTM 410
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DRTM 420
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DRTM 430
C           MUST BE FURNISHED BY THE USER.                              DRTM 440
C                                                                       DRTM 450
C        METHOD                                                         DRTM 460
C           SOLUTION OF EQUATION FCT(X)=0 IS DONE BY MEANS OF MUELLER-S DRTM 470
C           ITERATION METHOD OF SUCCESSIVE BISECTIONS AND INVERSE       DRTM 480
C           PARABOLIC INTERPOLATION, WHICH STARTS AT THE INITIAL BOUNDS DRTM 490
C           XLI AND XRI. CONVERGENCE IS QUADRATIC IF THE DERIVATIVE OF  DRTM 500
C           FCT(X) AT ROOT X IS NOT EQUAL TO ZERO. ONE ITERATION STEP   DRTM 510
C           REQUIRES TWO EVALUATIONS OF FCT(X). FOR TEST ON SATISFACTORYDRTM 520
C           ACCURACY SEE FORMULAE (3,4) OF MATHEMATICAL DESCRIPTION.    DRTM 530
C           FOR REFERENCE, SEE G. K. KRISTIANSEN, ZERO OF ARBITRARY     DRTM 540
C           FUNCTION, BIT, VOL. 3 (1963), PP.205-206.                   DRTM 550
C                                                                       DRTM 560
C     ..................................................................DRTM 570
C                                                                       DRTM 580
      SUBROUTINE DRTMI(X,F,FCT,XLI,XRI,EPS,IEND,IER)                    DRTM 590
C                                                                       DRTM 600
C                                                                       DRTM 610
      DOUBLE PRECISION X,F,FCT,XLI,XRI,XL,XR,FL,FR,TOL,TOLF,A,DX,XM,FM  DRTM 620
C                                                                       DRTM 630
C     PREPARE ITERATION                                                 DRTM 640
      IER=0                                                             DRTM 650
      XL=XLI                                                            DRTM 660
      XR=XRI                                                            DRTM 670
      X=XL                                                              DRTM 680
      TOL=X                                                             DRTM 690
      F=FCT(TOL)                                                        DRTM 700
      IF(F)1,16,1                                                       DRTM 710
    1 FL=F                                                              DRTM 720
      X=XR                                                              DRTM 730
      TOL=X                                                             DRTM 740
      F=FCT(TOL)                                                        DRTM 750
      IF(F)2,16,2                                                       DRTM 760
    2 FR=F                                                              DRTM 770
      IF(DSIGN(1.D0,FL)+DSIGN(1.D0,FR))25,3,25                          DRTM 780
C                                                                       DRTM 790
C     BASIC ASSUMPTION FL*FR LESS THAN 0 IS SATISFIED.                  DRTM 800
C     GENERATE TOLERANCE FOR FUNCTION VALUES.                           DRTM 810
    3 I=0                                                               DRTM 820
      TOLF=100.*EPS                                                     DRTM 830
C                                                                       DRTM 840
C                                                                       DRTM 850
C     START ITERATION LOOP                                              DRTM 860
    4 I=I+1                                                             DRTM 870
C                                                                       DRTM 880
C     START BISECTION LOOP                                              DRTM 890
      DO 13 K=1,IEND                                                    DRTM 900
      X=.5D0*(XL+XR)                                                    DRTM 910
      TOL=X                                                             DRTM 920
      F=FCT(TOL)                                                        DRTM 930
      IF(F)5,16,5                                                       DRTM 940
    5 IF(DSIGN(1.D0,F)+DSIGN(1.D0,FR))7,6,7                             DRTM 950
C                                                                       DRTM 960
C     INTERCHANGE XL AND XR IN ORDER TO GET THE SAME SIGN IN F AND FR   DRTM 970
    6 TOL=XL                                                            DRTM 980
      XL=XR                                                             DRTM 990
      XR=TOL                                                            DRTM1000
      TOL=FL                                                            DRTM1010
      FL=FR                                                             DRTM1020
      FR=TOL                                                            DRTM1030
    7 TOL=F-FL                                                          DRTM1040
      A=F*TOL                                                           DRTM1050
      A=A+A                                                             DRTM1060
      IF(A-FR*(FR-FL))8,9,9                                             DRTM1070
    8 IF(I-IEND)17,17,9                                                 DRTM1080
    9 XR=X                                                              DRTM1090
      FR=F                                                              DRTM1100
C                                                                       DRTM1110
C     TEST ON SATISFACTORY ACCURACY IN BISECTION LOOP                   DRTM1120
      TOL=EPS                                                           DRTM1130
      A=DABS(XR)                                                        DRTM1140
      IF(A-1.D0)11,11,10                                                DRTM1150
   10 TOL=TOL*A                                                         DRTM1160
   11 IF(DABS(XR-XL)-TOL)12,12,13                                       DRTM1170
   12 IF(DABS(FR-FL)-TOLF)14,14,13                                      DRTM1180
   13 CONTINUE                                                          DRTM1190
C     END OF BISECTION LOOP                                             DRTM1200
C                                                                       DRTM1210
C     NO CONVERGENCE AFTER IEND ITERATION STEPS FOLLOWED BY IEND        DRTM1220
C     SUCCESSIVE STEPS OF BISECTION OR STEADILY INCREASING FUNCTION     DRTM1230
C     VALUES AT RIGHT BOUNDS. ERROR RETURN.                             DRTM1240
      IER=1                                                             DRTM1250
   14 IF(DABS(FR)-DABS(FL))16,16,15                                     DRTM1260
   15 X=XL                                                              DRTM1270
      F=FL                                                              DRTM1280
   16 RETURN                                                            DRTM1290
C                                                                       DRTM1300
C     COMPUTATION OF ITERATED X-VALUE BY INVERSE PARABOLIC INTERPOLATIONDRTM1310
   17 A=FR-F                                                            DRTM1320
      DX=(X-XL)*FL*(1.D0+F*(A-TOL)/(A*(FR-FL)))/TOL                     DRTM1330
      XM=X                                                              DRTM1340
      FM=F                                                              DRTM1350
      X=XL-DX                                                           DRTM1360
      TOL=X                                                             DRTM1370
      F=FCT(TOL)                                                        DRTM1380
      IF(F)18,16,18                                                     DRTM1390
C                                                                       DRTM1400
C     TEST ON SATISFACTORY ACCURACY IN ITERATION LOOP                   DRTM1410
   18 TOL=EPS                                                           DRTM1420
      A=DABS(X)                                                         DRTM1430
      IF(A-1.D0)20,20,19                                                DRTM1440
   19 TOL=TOL*A                                                         DRTM1450
   20 IF(DABS(DX)-TOL)21,21,22                                          DRTM1460
   21 IF(DABS(F)-TOLF)16,16,22                                          DRTM1470
C                                                                       DRTM1480
C     PREPARATION OF NEXT BISECTION LOOP                                DRTM1490
   22 IF(DSIGN(1.D0,F)+DSIGN(1.D0,FL))24,23,24                          DRTM1500
   23 XR=X                                                              DRTM1510
      FR=F                                                              DRTM1520
      GO TO 4                                                           DRTM1530
   24 XL=X                                                              DRTM1540
      FL=F                                                              DRTM1550
      XR=XM                                                             DRTM1560
      FR=FM                                                             DRTM1570
      GO TO 4                                                           DRTM1580
C     END OF ITERATION LOOP                                             DRTM1590
C                                                                       DRTM1600
C                                                                       DRTM1610
C     ERROR RETURN IN CASE OF WRONG INPUT DATA                          DRTM1620
   25 IER=2                                                             DRTM1630
      RETURN                                                            DRTM1640
      END                                                               DRTM1650
