C                                                                       DRTN  10
C     ..................................................................DRTN  20
C                                                                       DRTN  30
C        SUBROUTINE DRTNI                                               DRTN  40
C                                                                       DRTN  50
C        PURPOSE                                                        DRTN  60
C           TO SOLVE GENERAL NONLINEAR EQUATIONS OF THE FORM F(X)=0     DRTN  70
C           BY MEANS OF NEWTON-S ITERATION METHOD.                      DRTN  80
C                                                                       DRTN  90
C        USAGE                                                          DRTN 100
C           CALL DRTNI (X,F,DERF,FCT,XST,EPS,IEND,IER)                  DRTN 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT.               DRTN 120
C                                                                       DRTN 130
C        DESCRIPTION OF PARAMETERS                                      DRTN 140
C           X      - DOUBLE PRECISION RESULTANT ROOT OF EQUATION F(X)=0.DRTN 150
C           F      - DOUBLE PRECISION RESULTANT FUNCTION VALUE AT       DRTN 160
C                    ROOT X.                                            DRTN 170
C           DERF   - DOUBLE PRECISION RESULTANT VALUE OF DERIVATIVE     DRTN 180
C                    AT ROOT X.                                         DRTN 190
C           FCT    - NAME OF THE EXTERNAL SUBROUTINE USED. IT COMPUTES  DRTN 200
C                    TO GIVEN ARGUMENT X FUNCTION VALUE F AND DERIVATIVEDRTN 210
C                    DERF. ITS PARAMETER LIST MUST BE X,F,DERF, WHERE   DRTN 220
C                    ALL PARAMETERS ARE DOUBLE PRECISION.               DRTN 230
C           XST    - DOUBLE PRECISION INPUT VALUE WHICH SPECIFIES THE   DRTN 240
C                    INITIAL GUESS OF THE ROOT X.                       DRTN 250
C           EPS    - SINGLE PRECISION INPUT VALUE WHICH SPECIFIES THE   DRTN 260
C                    UPPER BOUND OF THE ERROR OF RESULT X.              DRTN 270
C           IEND   - MAXIMUM NUMBER OF ITERATION STEPS SPECIFIED.       DRTN 280
C           IER    - RESULTANT ERROR PARAMETER CODED AS FOLLOWS         DRTN 290
C                     IER=0 - NO ERROR,                                 DRTN 300
C                     IER=1 - NO CONVERGENCE AFTER IEND ITERATION STEPS,DRTN 310
C                     IER=2 - AT ANY ITERATION STEP DERIVATIVE DERF WAS DRTN 320
C                             EQUAL TO ZERO.                            DRTN 330
C                                                                       DRTN 340
C        REMARKS                                                        DRTN 350
C           THE PROCEDURE IS BYPASSED AND GIVES THE ERROR MESSAGE IER=2 DRTN 360
C           IF AT ANY ITERATION STEP DERIVATIVE OF F(X) IS EQUAL TO 0.  DRTN 370
C           POSSIBLY THE PROCEDURE WOULD BE SUCCESSFUL IF IT IS STARTED DRTN 380
C           ONCE MORE WITH ANOTHER INITIAL GUESS XST.                   DRTN 390
C                                                                       DRTN 400
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DRTN 410
C           THE EXTERNAL SUBROUTINE FCT(X,F,DERF) MUST BE FURNISHED     DRTN 420
C           BY THE USER.                                                DRTN 430
C                                                                       DRTN 440
C        METHOD                                                         DRTN 450
C           SOLUTION OF EQUATION F(X)=0 IS DONE BY MEANS OF NEWTON-S    DRTN 460
C           ITERATION METHOD, WHICH STARTS AT THE INITIAL GUESS XST OF  DRTN 470
C           A ROOT X. CONVERGENCE IS QUADRATIC IF THE DERIVATIVE OF     DRTN 480
C           F(X) AT ROOT X IS NOT EQUAL TO ZERO. ONE ITERATION STEP     DRTN 490
C           REQUIRES ONE EVALUATION OF F(X) AND ONE EVALUATION OF THE   DRTN 500
C           DERIVATIVE OF F(X). FOR TEST ON SATISFACTORY ACCURACY SEE   DRTN 510
C           FORMULAE (2) OF MATHEMATICAL DESCRIPTION.                   DRTN 520
C           FOR REFERENCE, SEE R. ZURMUEHL, PRAKTISCHE MATHEMATIK FUER  DRTN 530
C           INGENIEURE UND PHYSIKER, SPRINGER, BERLIN/GOETTINGEN/       DRTN 540
C           HEIDELBERG, 1963, PP.12-17.                                 DRTN 550
C                                                                       DRTN 560
C     ..................................................................DRTN 570
C                                                                       DRTN 580
      SUBROUTINE DRTNI(X,F,DERF,FCT,XST,EPS,IEND,IER)                   DRTN 590
C                                                                       DRTN 600
C                                                                       DRTN 610
      DOUBLE PRECISION X,F,DERF,XST,TOL,TOLF,DX,A                       DRTN 620
C                                                                       DRTN 630
C     PREPARE ITERATION                                                 DRTN 640
      IER=0                                                             DRTN 650
      X=XST                                                             DRTN 660
      TOL=X                                                             DRTN 670
      CALL FCT(TOL,F,DERF)                                              DRTN 680
      TOLF=100.*EPS                                                     DRTN 690
C                                                                       DRTN 700
C                                                                       DRTN 710
C     START ITERATION LOOP                                              DRTN 720
      DO 6 I=1,IEND                                                     DRTN 730
      IF(F)1,7,1                                                        DRTN 740
C                                                                       DRTN 750
C     EQUATION IS NOT SATISFIED BY X                                    DRTN 760
    1 IF(DERF)2,8,2                                                     DRTN 770
C                                                                       DRTN 780
C     ITERATION IS POSSIBLE                                             DRTN 790
    2 DX=F/DERF                                                         DRTN 800
      X=X-DX                                                            DRTN 810
      TOL=X                                                             DRTN 820
      CALL FCT(TOL,F,DERF)                                              DRTN 830
C                                                                       DRTN 840
C     TEST ON SATISFACTORY ACCURACY                                     DRTN 850
      TOL=EPS                                                           DRTN 860
      A=DABS(X)                                                         DRTN 870
      IF(A-1.D0)4,4,3                                                   DRTN 880
    3 TOL=TOL*A                                                         DRTN 890
    4 IF(DABS(DX)-TOL)5,5,6                                             DRTN 900
    5 IF(DABS(F)-TOLF)7,7,6                                             DRTN 910
    6 CONTINUE                                                          DRTN 920
C     END OF ITERATION LOOP                                             DRTN 930
C                                                                       DRTN 940
C                                                                       DRTN 950
C     NO CONVERGENCE AFTER IEND ITERATION STEPS. ERROR RETURN.          DRTN 960
      IER=1                                                             DRTN 970
    7 RETURN                                                            DRTN 980
C                                                                       DRTN 990
C     ERROR RETURN IN CASE OF ZERO DIVISOR                              DRTN1000
    8 IER=2                                                             DRTN1010
      RETURN                                                            DRTN1020
      END                                                               DRTN1030
