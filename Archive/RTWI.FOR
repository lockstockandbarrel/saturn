C                                                                       RTWI  10
C     ..................................................................RTWI  20
C                                                                       RTWI  30
C        SUBROUTINE RTWI                                                RTWI  40
C                                                                       RTWI  50
C        PURPOSE                                                        RTWI  60
C           TO SOLVE GENERAL NONLINEAR EQUATIONS OF THE FORM X=FCT(X)   RTWI  70
C           BY MEANS OF WEGSTEIN-S ITERATION METHOD.                    RTWI  80
C                                                                       RTWI  90
C        USAGE                                                          RTWI 100
C           CALL RTWI (X,VAL,FCT,XST,EPS,IEND,IER)                      RTWI 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT.               RTWI 120
C                                                                       RTWI 130
C        DESCRIPTION OF PARAMETERS                                      RTWI 140
C           X      - RESULTANT ROOT OF EQUATION X=FCT(X).               RTWI 150
C           VAL    - RESULTANT VALUE OF X-FCT(X) AT ROOT X.             RTWI 160
C           FCT    - NAME OF THE EXTERNAL FUNCTION SUBPROGRAM USED.     RTWI 170
C           XST    - INPUT VALUE WHICH SPECIFIES THE INITIAL GUESS OF   RTWI 180
C                    THE ROOT X.                                        RTWI 190
C           EPS    - INPUT VALUE WHICH SPECIFIES THE UPPER BOUND OF THE RTWI 200
C                    ERROR OF RESULT X.                                 RTWI 210
C           IEND   - MAXIMUM NUMBER OF ITERATION STEPS SPECIFIED.       RTWI 220
C           IER    - RESULTANT ERROR PARAMETER CODED AS FOLLOWS         RTWI 230
C                     IER=0 - NO ERROR,                                 RTWI 240
C                     IER=1 - NO CONVERGENCE AFTER IEND ITERATION STEPS,RTWI 250
C                     IER=2 - AT ANY ITERATION STEP THE DENOMINATOR OF  RTWI 260
C                             ITERATION FORMULA WAS EQUAL TO ZERO.      RTWI 270
C                                                                       RTWI 280
C        REMARKS                                                        RTWI 290
C           THE PROCEDURE IS BYPASSED AND GIVES THE ERROR MESSAGE IER=2 RTWI 300
C           IF AT ANY ITERATION STEP THE DENOMINATOR OF ITERATION       RTWI 310
C           FORMULA WAS EQUAL TO ZERO. THAT MEANS THAT THERE IS AT      RTWI 320
C           LEAST ONE POINT IN THE RANGE IN WHICH ITERATION MOVES WITH  RTWI 330
C           DERIVATIVE OF FCT(X) EQUAL TO 1.                            RTWI 340
C                                                                       RTWI 350
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  RTWI 360
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   RTWI 370
C           BY THE USER.                                                RTWI 380
C                                                                       RTWI 390
C        METHOD                                                         RTWI 400
C           SOLUTION OF EQUATION X=FCT(X) IS DONE BY MEANS OF           RTWI 410
C           WEGSTEIN-S ITERATION METHOD, WHICH STARTS AT THE INITIAL    RTWI 420
C           GUESS XST OF A ROOT X. ONE ITERATION STEP REQUIRES ONE      RTWI 430
C           EVALUATION OF FCT(X). FOR TEST ON SATISFACTORY ACCURACY SEE RTWI 440
C           FORMULAE (2) OF MATHEMATICAL DESCRIPTION.                   RTWI 450
C           FOR REFERENCE, SEE                                          RTWI 460
C           (1) G. N. LANCE, NUMERICAL METHODS FOR HIGH SPEED COMPUTERS,RTWI 470
C               ILIFFE, LONDON, 1960, PP.134-138,                       RTWI 480
C           (2) J. WEGSTEIN, ALGORITHM 2, CACM, VOL.3, ISS.2 (1960),    RTWI 490
C               PP.74,                                                  RTWI 500
C           (3) H.C. THACHER, ALGORITHM 15, CACM, VOL.3, ISS.8 (1960),  RTWI 510
C               PP.475,                                                 RTWI 520
C           (4) J.G. HERRIOT, ALGORITHM 26, CACM, VOL.3, ISS.11 (1960), RTWI 530
C               PP.603.                                                 RTWI 540
C                                                                       RTWI 550
C     ..................................................................RTWI 560
C                                                                       RTWI 570
      SUBROUTINE RTWI(X,VAL,FCT,XST,EPS,IEND,IER)                       RTWI 580
C                                                                       RTWI 590
C                                                                       RTWI 600
C     PREPARE ITERATION                                                 RTWI 610
      IER=0                                                             RTWI 620
      TOL=XST                                                           RTWI 630
      X=FCT(TOL)                                                        RTWI 640
      A=X-XST                                                           RTWI 650
      B=-A                                                              RTWI 660
      TOL=X                                                             RTWI 670
      VAL=X-FCT(TOL)                                                    RTWI 680
C                                                                       RTWI 690
C                                                                       RTWI 700
C     START ITERATION LOOP                                              RTWI 710
      DO 6 I=1,IEND                                                     RTWI 720
      IF(VAL)1,7,1                                                      RTWI 730
C                                                                       RTWI 740
C     EQUATION IS NOT SATISFIED BY X                                    RTWI 750
    1 B=B/VAL-1.                                                        RTWI 760
      IF(B)2,8,2                                                        RTWI 770
C                                                                       RTWI 780
C     ITERATION IS POSSIBLE                                             RTWI 790
    2 A=A/B                                                             RTWI 800
      X=X+A                                                             RTWI 810
      B=VAL                                                             RTWI 820
      TOL=X                                                             RTWI 830
      VAL=X-FCT(TOL)                                                    RTWI 840
C                                                                       RTWI 850
C     TEST ON SATISFACTORY ACCURACY                                     RTWI 860
      TOL=EPS                                                           RTWI 870
      D=ABS(X)                                                          RTWI 880
      IF(D-1.)4,4,3                                                     RTWI 890
    3 TOL=TOL*D                                                         RTWI 900
    4 IF(ABS(A)-TOL)5,5,6                                               RTWI 910
    5 IF(ABS(VAL)-10.*TOL)7,7,6                                         RTWI 920
    6 CONTINUE                                                          RTWI 930
C     END OF ITERATION LOOP                                             RTWI 940
C                                                                       RTWI 950
C                                                                       RTWI 960
C     NO CONVERGENCE AFTER IEND ITERATION STEPS. ERROR RETURN.          RTWI 970
      IER=1                                                             RTWI 980
    7 RETURN                                                            RTWI 990
C                                                                       RTWI1000
C     ERROR RETURN IN CASE OF ZERO DIVISOR                              RTWI1010
    8 IER=2                                                             RTWI1020
      RETURN                                                            RTWI1030
      END                                                               RTWI1040
