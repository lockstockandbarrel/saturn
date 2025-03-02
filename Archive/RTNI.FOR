C                                                                       RTNI  10
C     ..................................................................RTNI  20
C                                                                       RTNI  30
C        SUBROUTINE RTNI                                                RTNI  40
C                                                                       RTNI  50
C        PURPOSE                                                        RTNI  60
C           TO SOLVE GENERAL NONLINEAR EQUATIONS OF THE FORM F(X)=0     RTNI  70
C           BY MEANS OF NEWTON-S ITERATION METHOD.                      RTNI  80
C                                                                       RTNI  90
C        USAGE                                                          RTNI 100
C           CALL RTNI (X,F,DERF,FCT,XST,EPS,IEND,IER)                   RTNI 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT.               RTNI 120
C                                                                       RTNI 130
C        DESCRIPTION OF PARAMETERS                                      RTNI 140
C           X      - RESULTANT ROOT OF EQUATION F(X)=0.                 RTNI 150
C           F      - RESULTANT FUNCTION VALUE AT ROOT X.                RTNI 160
C           DERF   - RESULTANT VALUE OF DERIVATIVE AT ROOT X.           RTNI 170
C           FCT    - NAME OF THE EXTERNAL SUBROUTINE USED. IT COMPUTES  RTNI 180
C                    TO GIVEN ARGUMENT X FUNCTION VALUE F AND DERIVATIVERTNI 190
C                    DERF. ITS PARAMETER LIST MUST BE X,F,DERF.         RTNI 200
C           XST    - INPUT VALUE WHICH SPECIFIES THE INITIAL GUESS OF   RTNI 210
C                    THE ROOT X.                                        RTNI 220
C           EPS    - INPUT VALUE WHICH SPECIFIES THE UPPER BOUND OF THE RTNI 230
C                    ERROR OF RESULT X.                                 RTNI 240
C           IEND   - MAXIMUM NUMBER OF ITERATION STEPS SPECIFIED.       RTNI 250
C           IER    - RESULTANT ERROR PARAMETER CODED AS FOLLOWS         RTNI 260
C                     IER=0 - NO ERROR,                                 RTNI 270
C                     IER=1 - NO CONVERGENCE AFTER IEND ITERATION STEPS,RTNI 280
C                     IER=2 - AT ANY ITERATION STEP DERIVATIVE DERF WAS RTNI 290
C                             EQUAL TO ZERO.                            RTNI 300
C                                                                       RTNI 310
C        REMARKS                                                        RTNI 320
C           THE PROCEDURE IS BYPASSED AND GIVES THE ERROR MESSAGE IER=2 RTNI 330
C           IF AT ANY ITERATION STEP DERIVATIVE OF F(X) IS EQUAL TO 0.  RTNI 340
C           POSSIBLY THE PROCEDURE WOULD BE SUCCESSFUL IF IT IS STARTED RTNI 350
C           ONCE MORE WITH ANOTHER INITIAL GUESS XST.                   RTNI 360
C                                                                       RTNI 370
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  RTNI 380
C           THE EXTERNAL SUBROUTINE FCT(X,F,DERF) MUST BE FURNISHED     RTNI 390
C           BY THE USER.                                                RTNI 400
C                                                                       RTNI 410
C        METHOD                                                         RTNI 420
C           SOLUTION OF EQUATION F(X)=0 IS DONE BY MEANS OF NEWTON-S    RTNI 430
C           ITERATION METHOD, WHICH STARTS AT THE INITIAL GUESS XST OF  RTNI 440
C           A ROOT X. CONVERGENCE IS QUADRATIC IF THE DERIVATIVE OF     RTNI 450
C           F(X) AT ROOT X IS NOT EQUAL TO ZERO. ONE ITERATION STEP     RTNI 460
C           REQUIRES ONE EVALUATION OF F(X) AND ONE EVALUATION OF THE   RTNI 470
C           DERIVATIVE OF F(X). FOR TEST ON SATISFACTORY ACCURACY SEE   RTNI 480
C           FORMULAE (2) OF MATHEMATICAL DESCRIPTION.                   RTNI 490
C           FOR REFERENCE, SEE R. ZURMUEHL, PRAKTISCHE MATHEMATIK FUER  RTNI 500
C           INGENIEURE UND PHYSIKER, SPRINGER, BERLIN/GOETTINGEN/       RTNI 510
C           HEIDELBERG, 1963, PP.12-17.                                 RTNI 520
C                                                                       RTNI 530
C     ..................................................................RTNI 540
C                                                                       RTNI 550
      SUBROUTINE RTNI(X,F,DERF,FCT,XST,EPS,IEND,IER)                    RTNI 560
C                                                                       RTNI 570
C                                                                       RTNI 580
C     PREPARE ITERATION                                                 RTNI 590
      IER=0                                                             RTNI 600
      X=XST                                                             RTNI 610
      TOL=X                                                             RTNI 620
      CALL FCT(TOL,F,DERF)                                              RTNI 630
      TOLF=100.*EPS                                                     RTNI 640
C                                                                       RTNI 650
C                                                                       RTNI 660
C     START ITERATION LOOP                                              RTNI 670
      DO 6 I=1,IEND                                                     RTNI 680
      IF(F)1,7,1                                                        RTNI 690
C                                                                       RTNI 700
C     EQUATION IS NOT SATISFIED BY X                                    RTNI 710
    1 IF(DERF)2,8,2                                                     RTNI 720
C                                                                       RTNI 730
C     ITERATION IS POSSIBLE                                             RTNI 740
    2 DX=F/DERF                                                         RTNI 750
      X=X-DX                                                            RTNI 760
      TOL=X                                                             RTNI 770
      CALL FCT(TOL,F,DERF)                                              RTNI 780
C                                                                       RTNI 790
C     TEST ON SATISFACTORY ACCURACY                                     RTNI 800
      TOL=EPS                                                           RTNI 810
      A=ABS(X)                                                          RTNI 820
      IF(A-1.)4,4,3                                                     RTNI 830
    3 TOL=TOL*A                                                         RTNI 840
    4 IF(ABS(DX)-TOL)5,5,6                                              RTNI 850
    5 IF(ABS(F)-TOLF)7,7,6                                              RTNI 860
    6 CONTINUE                                                          RTNI 870
C     END OF ITERATION LOOP                                             RTNI 880
C                                                                       RTNI 890
C                                                                       RTNI 900
C     NO CONVERGENCE AFTER IEND ITERATION STEPS. ERROR RETURN.          RTNI 910
      IER=1                                                             RTNI 920
    7 RETURN                                                            RTNI 930
C                                                                       RTNI 940
C     ERROR RETURN IN CASE OF ZERO DIVISOR                              RTNI 950
    8 IER=2                                                             RTNI 960
      RETURN                                                            RTNI 970
      END                                                               RTNI 980
