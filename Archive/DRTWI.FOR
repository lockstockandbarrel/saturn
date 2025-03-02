C                                                                       DRTW  10
C     ..................................................................DRTW  20
C                                                                       DRTW  30
C        SUBROUTINE DRTWI                                               DRTW  40
C                                                                       DRTW  50
C        PURPOSE                                                        DRTW  60
C           TO SOLVE GENERAL NONLINEAR EQUATIONS OF THE FORM X=FCT(X)   DRTW  70
C           BY MEANS OF WEGSTEIN-S ITERATION METHOD.                    DRTW  80
C                                                                       DRTW  90
C        USAGE                                                          DRTW 100
C           CALL DRTWI (X,VAL,FCT,XST,EPS,IEND,IER)                     DRTW 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT.               DRTW 120
C                                                                       DRTW 130
C        DESCRIPTION OF PARAMETERS                                      DRTW 140
C           X      - DOUBLE PRECISION RESULTANT ROOT OF EQUATION        DRTW 150
C                    X=FCT(X).                                          DRTW 160
C           VAL    - DOUBLE PRECISION RESULTANT VALUE OF X-FCT(X)       DRTW 170
C                    AT ROOT X.                                         DRTW 180
C           FCT    - NAME OF THE EXTERNAL DOUBLE PRECISION FUNCTION     DRTW 190
C                    SUBPROGRAM USED.                                   DRTW 200
C           XST    - DOUBLE PRECISION INPUT VALUE WHICH SPECIFIES THE   DRTW 210
C                    INITIAL GUESS OF THE ROOT X.                       DRTW 220
C           EPS    - SINGLE PRECISION INPUT VALUE WHICH SPECIFIES THE   DRTW 230
C                    UPPER BOUND OF THE ERROR OF RESULT X.              DRTW 240
C           IEND   - MAXIMUM NUMBER OF ITERATION STEPS SPECIFIED.       DRTW 250
C           IER    - RESULTANT ERROR PARAMETER CODED AS FOLLOWS         DRTW 260
C                     IER=0 - NO ERROR,                                 DRTW 270
C                     IER=1 - NO CONVERGENCE AFTER IEND ITERATION STEPS,DRTW 280
C                     IER=2 - AT ANY ITERATION STEP THE DENOMINATOR OF  DRTW 290
C                             ITERATION FORMULA WAS EQUAL TO ZERO.      DRTW 300
C                                                                       DRTW 310
C        REMARKS                                                        DRTW 320
C           THE PROCEDURE IS BYPASSED AND GIVES THE ERROR MESSAGE IER=2 DRTW 330
C           IF AT ANY ITERATION STEP THE DENOMINATOR OF ITERATION       DRTW 340
C           FORMULA WAS EQUAL TO ZERO. THAT MEANS THAT THERE IS AT      DRTW 350
C           LEAST ONE POINT IN THE RANGE IN WHICH ITERATION MOVES WITH  DRTW 360
C           DERIVATIVE OF FCT(X) EQUAL TO 1.                            DRTW 370
C                                                                       DRTW 380
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DRTW 390
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DRTW 400
C           MUST BE FURNISHED BY THE USER.                              DRTW 410
C                                                                       DRTW 420
C        METHOD                                                         DRTW 430
C           SOLUTION OF EQUATION X=FCT(X) IS DONE BY MEANS OF           DRTW 440
C           WEGSTEIN-S ITERATION METHOD, WHICH STARTS AT THE INITIAL    DRTW 450
C           GUESS XST OF A ROOT X. ONE ITERATION STEP REQUIRES ONE      DRTW 460
C           EVALUATION OF FCT(X). FOR TEST ON SATISFACTORY ACCURACY SEE DRTW 470
C           FORMULAE (2) OF MATHEMATICAL DESCRIPTION.                   DRTW 480
C           FOR REFERENCE, SEE                                          DRTW 490
C           (1) G. N. LANCE, NUMERICAL METHODS FOR HIGH SPEED COMPUTERS,DRTW 500
C               ILIFFE, LONDON, 1960, PP.134-138,                       DRTW 510
C           (2) J. WEGSTEIN, ALGORITHM 2, CACM, VOL.3, ISS.2 (1960),    DRTW 520
C               PP.74,                                                  DRTW 530
C           (3) H.C. THACHER, ALGORITHM 15, CACM, VOL.3, ISS.8 (1960),  DRTW 540
C               PP.475,                                                 DRTW 550
C           (4) J.G. HERRIOT, ALGORITHM 26, CACM, VOL.3, ISS.11 (1960), DRTW 560
C               PP.603.                                                 DRTW 570
C                                                                       DRTW 580
C     ..................................................................DRTW 590
C                                                                       DRTW 600
      SUBROUTINE DRTWI(X,VAL,FCT,XST,EPS,IEND,IER)                      DRTW 610
C                                                                       DRTW 620
C                                                                       DRTW 630
      DOUBLE PRECISION X,VAL,FCT,XST,A,B,D,TOL                          DRTW 640
C                                                                       DRTW 650
C     PREPARE ITERATION                                                 DRTW 660
      IER=0                                                             DRTW 670
      TOL=XST                                                           DRTW 680
      X=FCT(TOL)                                                        DRTW 690
      A=X-XST                                                           DRTW 700
      B=-A                                                              DRTW 710
      TOL=X                                                             DRTW 720
      VAL=X-FCT(TOL)                                                    DRTW 730
C                                                                       DRTW 740
C                                                                       DRTW 750
C     START ITERATION LOOP                                              DRTW 760
      DO 6 I=1,IEND                                                     DRTW 770
      IF(VAL)1,7,1                                                      DRTW 780
C                                                                       DRTW 790
C     EQUATION IS NOT SATISFIED BY X                                    DRTW 800
    1 B=B/VAL-1.D0                                                      DRTW 810
      IF(B)2,8,2                                                        DRTW 820
C                                                                       DRTW 830
C     ITERATION IS POSSIBLE                                             DRTW 840
    2 A=A/B                                                             DRTW 850
      X=X+A                                                             DRTW 860
      B=VAL                                                             DRTW 870
      TOL=X                                                             DRTW 880
      VAL=X-FCT(TOL)                                                    DRTW 890
C                                                                       DRTW 900
C     TEST ON SATISFACTORY ACCURACY                                     DRTW 910
      TOL=EPS                                                           DRTW 920
      D=DABS(X)                                                         DRTW 930
      IF(D-1.D0)4,4,3                                                   DRTW 940
    3 TOL=TOL*D                                                         DRTW 950
    4 IF(DABS(A)-TOL)5,5,6                                              DRTW 960
    5 IF(DABS(VAL)-1.D1*TOL)7,7,6                                       DRTW 970
    6 CONTINUE                                                          DRTW 980
C     END OF ITERATION LOOP                                             DRTW 990
C                                                                       DRTW1000
C                                                                       DRTW1010
C     NO CONVERGENCE AFTER IEND ITERATION STEPS. ERROR RETURN.          DRTW1020
      IER=1                                                             DRTW1030
    7 RETURN                                                            DRTW1040
C                                                                       DRTW1050
C     ERROR RETURN IN CASE OF ZERO DIVISOR                              DRTW1060
    8 IER=2                                                             DRTW1070
      RETURN                                                            DRTW1080
      END                                                               DRTW1090
