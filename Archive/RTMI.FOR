C                                                                       RTMI  10
C     ..................................................................RTMI  20
C                                                                       RTMI  30
C        SUBROUTINE RTMI                                                RTMI  40
C                                                                       RTMI  50
C        PURPOSE                                                        RTMI  60
C           TO SOLVE GENERAL NONLINEAR EQUATIONS OF THE FORM FCT(X)=0   RTMI  70
C           BY MEANS OF MUELLER-S ITERATION METHOD.                     RTMI  80
C                                                                       RTMI  90
C        USAGE                                                          RTMI 100
C           CALL RTMI (X,F,FCT,XLI,XRI,EPS,IEND,IER)                    RTMI 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT.               RTMI 120
C                                                                       RTMI 130
C        DESCRIPTION OF PARAMETERS                                      RTMI 140
C           X      - RESULTANT ROOT OF EQUATION FCT(X)=0.               RTMI 150
C           F      - RESULTANT FUNCTION VALUE AT ROOT X.                RTMI 160
C           FCT    - NAME OF THE EXTERNAL FUNCTION SUBPROGRAM USED.     RTMI 170
C           XLI    - INPUT VALUE WHICH SPECIFIES THE INITIAL LEFT BOUND RTMI 180
C                    OF THE ROOT X.                                     RTMI 190
C           XRI    - INPUT VALUE WHICH SPECIFIES THE INITIAL RIGHT BOUNDRTMI 200
C                    OF THE ROOT X.                                     RTMI 210
C           EPS    - INPUT VALUE WHICH SPECIFIES THE UPPER BOUND OF THE RTMI 220
C                    ERROR OF RESULT X.                                 RTMI 230
C           IEND   - MAXIMUM NUMBER OF ITERATION STEPS SPECIFIED.       RTMI 240
C           IER    - RESULTANT ERROR PARAMETER CODED AS FOLLOWS         RTMI 250
C                     IER=0 - NO ERROR,                                 RTMI 260
C                     IER=1 - NO CONVERGENCE AFTER IEND ITERATION STEPS RTMI 270
C                             FOLLOWED BY IEND SUCCESSIVE STEPS OF      RTMI 280
C                             BISECTION,                                RTMI 290
C                     IER=2 - BASIC ASSUMPTION FCT(XLI)*FCT(XRI) LESS   RTMI 300
C                             THAN OR EQUAL TO ZERO IS NOT SATISFIED.   RTMI 310
C                                                                       RTMI 320
C        REMARKS                                                        RTMI 330
C           THE PROCEDURE ASSUMES THAT FUNCTION VALUES AT INITIAL       RTMI 340
C           BOUNDS XLI AND XRI HAVE NOT THE SAME SIGN. IF THIS BASIC    RTMI 350
C           ASSUMPTION IS NOT SATISFIED BY INPUT VALUES XLI AND XRI, THERTMI 360
C           PROCEDURE IS BYPASSED AND GIVES THE ERROR MESSAGE IER=2.    RTMI 370
C                                                                       RTMI 380
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  RTMI 390
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   RTMI 400
C           BY THE USER.                                                RTMI 410
C                                                                       RTMI 420
C        METHOD                                                         RTMI 430
C           SOLUTION OF EQUATION FCT(X)=0 IS DONE BY MEANS OF MUELLER-S RTMI 440
C           ITERATION METHOD OF SUCCESSIVE BISECTIONS AND INVERSE       RTMI 450
C           PARABOLIC INTERPOLATION, WHICH STARTS AT THE INITIAL BOUNDS RTMI 460
C           XLI AND XRI. CONVERGENCE IS QUADRATIC IF THE DERIVATIVE OF  RTMI 470
C           FCT(X) AT ROOT X IS NOT EQUAL TO ZERO. ONE ITERATION STEP   RTMI 480
C           REQUIRES TWO EVALUATIONS OF FCT(X). FOR TEST ON SATISFACTORYRTMI 490
C           ACCURACY SEE FORMULAE (3,4) OF MATHEMATICAL DESCRIPTION.    RTMI 500
C           FOR REFERENCE, SEE G. K. KRISTIANSEN, ZERO OF ARBITRARY     RTMI 510
C           FUNCTION, BIT, VOL. 3 (1963), PP.205-206.                   RTMI 520
C                                                                       RTMI 530
C     ..................................................................RTMI 540
C                                                                       RTMI 550
      SUBROUTINE RTMI(X,F,FCT,XLI,XRI,EPS,IEND,IER)                     RTMI 560
C                                                                       RTMI 570
C                                                                       RTMI 580
C     PREPARE ITERATION                                                 RTMI 590
      IER=0                                                             RTMI 600
      XL=XLI                                                            RTMI 610
      XR=XRI                                                            RTMI 620
      X=XL                                                              RTMI 630
      TOL=X                                                             RTMI 640
      F=FCT(TOL)                                                        RTMI 650
      IF(F)1,16,1                                                       RTMI 660
    1 FL=F                                                              RTMI 670
      X=XR                                                              RTMI 680
      TOL=X                                                             RTMI 690
      F=FCT(TOL)                                                        RTMI 700
      IF(F)2,16,2                                                       RTMI 710
    2 FR=F                                                              RTMI 720
      IF(SIGN(1.,FL)+SIGN(1.,FR))25,3,25                                RTMI 730
C                                                                       RTMI 740
C     BASIC ASSUMPTION FL*FR LESS THAN 0 IS SATISFIED.                  RTMI 750
C     GENERATE TOLERANCE FOR FUNCTION VALUES.                           RTMI 760
    3 I=0                                                               RTMI 770
      TOLF=100.*EPS                                                     RTMI 780
C                                                                       RTMI 790
C                                                                       RTMI 800
C     START ITERATION LOOP                                              RTMI 810
    4 I=I+1                                                             RTMI 820
C                                                                       RTMI 830
C     START BISECTION LOOP                                              RTMI 840
      DO 13 K=1,IEND                                                    RTMI 850
      X=.5*(XL+XR)                                                      RTMI 860
      TOL=X                                                             RTMI 870
      F=FCT(TOL)                                                        RTMI 880
      IF(F)5,16,5                                                       RTMI 890
    5 IF(SIGN(1.,F)+SIGN(1.,FR))7,6,7                                   RTMI 900
C                                                                       RTMI 910
C     INTERCHANGE XL AND XR IN ORDER TO GET THE SAME SIGN IN F AND FR   RTMI 920
    6 TOL=XL                                                            RTMI 930
      XL=XR                                                             RTMI 940
      XR=TOL                                                            RTMI 950
      TOL=FL                                                            RTMI 960
      FL=FR                                                             RTMI 970
      FR=TOL                                                            RTMI 980
    7 TOL=F-FL                                                          RTMI 990
      A=F*TOL                                                           RTMI1000
      A=A+A                                                             RTMI1010
      IF(A-FR*(FR-FL))8,9,9                                             RTMI1020
    8 IF(I-IEND)17,17,9                                                 RTMI1030
    9 XR=X                                                              RTMI1040
      FR=F                                                              RTMI1050
C                                                                       RTMI1060
C     TEST ON SATISFACTORY ACCURACY IN BISECTION LOOP                   RTMI1070
      TOL=EPS                                                           RTMI1080
      A=ABS(XR)                                                         RTMI1090
      IF(A-1.)11,11,10                                                  RTMI1100
   10 TOL=TOL*A                                                         RTMI1110
   11 IF(ABS(XR-XL)-TOL)12,12,13                                        RTMI1120
   12 IF(ABS(FR-FL)-TOLF)14,14,13                                       RTMI1130
   13 CONTINUE                                                          RTMI1140
C     END OF BISECTION LOOP                                             RTMI1150
C                                                                       RTMI1160
C     NO CONVERGENCE AFTER IEND ITERATION STEPS FOLLOWED BY IEND        RTMI1170
C     SUCCESSIVE STEPS OF BISECTION OR STEADILY INCREASING FUNCTION     RTMI1180
C     VALUES AT RIGHT BOUNDS. ERROR RETURN.                             RTMI1190
      IER=1                                                             RTMI1200
   14 IF(ABS(FR)-ABS(FL))16,16,15                                       RTMI1210
   15 X=XL                                                              RTMI1220
      F=FL                                                              RTMI1230
   16 RETURN                                                            RTMI1240
C                                                                       RTMI1250
C     COMPUTATION OF ITERATED X-VALUE BY INVERSE PARABOLIC INTERPOLATIONRTMI1260
   17 A=FR-F                                                            RTMI1270
      DX=(X-XL)*FL*(1.+F*(A-TOL)/(A*(FR-FL)))/TOL                       RTMI1280
      XM=X                                                              RTMI1290
      FM=F                                                              RTMI1300
      X=XL-DX                                                           RTMI1310
      TOL=X                                                             RTMI1320
      F=FCT(TOL)                                                        RTMI1330
      IF(F)18,16,18                                                     RTMI1340
C                                                                       RTMI1350
C     TEST ON SATISFACTORY ACCURACY IN ITERATION LOOP                   RTMI1360
   18 TOL=EPS                                                           RTMI1370
      A=ABS(X)                                                          RTMI1380
      IF(A-1.)20,20,19                                                  RTMI1390
   19 TOL=TOL*A                                                         RTMI1400
   20 IF(ABS(DX)-TOL)21,21,22                                           RTMI1410
   21 IF(ABS(F)-TOLF)16,16,22                                           RTMI1420
C                                                                       RTMI1430
C     PREPARATION OF NEXT BISECTION LOOP                                RTMI1440
   22 IF(SIGN(1.,F)+SIGN(1.,FL))24,23,24                                RTMI1450
   23 XR=X                                                              RTMI1460
      FR=F                                                              RTMI1470
      GO TO 4                                                           RTMI1480
   24 XL=X                                                              RTMI1490
      FL=F                                                              RTMI1500
      XR=XM                                                             RTMI1510
      FR=FM                                                             RTMI1520
      GO TO 4                                                           RTMI1530
C     END OF ITERATION LOOP                                             RTMI1540
C                                                                       RTMI1550
C                                                                       RTMI1560
C     ERROR RETURN IN CASE OF WRONG INPUT DATA                          RTMI1570
   25 IER=2                                                             RTMI1580
      RETURN                                                            RTMI1590
      END                                                               RTMI1600
