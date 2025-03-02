C                                                                       DALI  10
C     ................................................................. DALI  20
C                                                                       DALI  30
C        SUBROUTINE DALI                                                DALI  40
C                                                                       DALI  50
C        PURPOSE                                                        DALI  60
C           TO INTERPOLATE FUNCTION VALUE Y FOR A GIVEN ARGUMENT VALUE  DALI  70
C           X USING A GIVEN TABLE (ARG,VAL) OF ARGUMENT AND FUNCTION    DALI  80
C           VALUES.                                                     DALI  90
C                                                                       DALI 100
C        USAGE                                                          DALI 110
C           CALL ALI (X,ARG,VAL,Y,NDIM,EPS,IER)                         DALI 120
C                                                                       DALI 130
C        DESCRIPTION OF PARAMETERS                                      DALI 140
C           X      - THE ARGUMENT VALUE SPECIFIED BY INPUT.             DALI 150
C           ARG    - THE INPUT VECTOR (DIMENSION NDIM) OF ARGUMENT      DALI 160
C                    VALUES OF THE TABLE (NOT DESTROYED).               DALI 170
C           VAL    - THE INPUT VECTOR (DIMENSION NDIM) OF FUNCTION      DALI 180
C                    VALUES OF THE TABLE (DESTROYED).                   DALI 190
C           Y      - THE RESULTING INTERPOLATED FUNCTION VALUE.         DALI 200
C           NDIM   - AN INPUT VALUE WHICH SPECIFIES THE NUMBER OF       DALI 210
C                    POINTS IN TABLE (ARG,VAL).                         DALI 220
C           EPS    - AN INPUT CONSTANT WHICH IS USED AS UPPER BOUND     DALI 230
C                    FOR THE ABSOLUTE ERROR.                            DALI 240
C           IER    - A RESULTING ERROR PARAMETER.                       DALI 250
C           NOTE:   X,ARG,YAL,X,Y,AND H ARE DOUBLE PRECISION            DALI 260
C        REMARKS                                                        DALI 270
C           (1) TABLE (ARG,VAL) SHOULD REPRESENT A SINGLE-VALUED        DALI 280
C               FUNCTION AND SHOULD BE STORED IN SUCH A WAY, THAT THE   DALI 290
C               DISTANCES ABS(ARG(I)-X) INCREASE WITH INCREASING        DALI 300
C               SUBSCRIPT I. TO GENERATE THIS ORDER IN TABLE (ARG,VAL), DALI 310
C               SUBROUTINES ATSG, ATSM OR ATSE COULD BE USED IN A       DALI 320
C               PREVIOUS STAGE.                                         DALI 330
C           (2) NO ACTION BESIDES ERROR MESSAGE IN CASE NDIM LESS       DALI 340
C               THAN 1.                                                 DALI 350
C           (3) INTERPOLATION IS TERMINATED EITHER IF THE DIFFERENCE    DALI 360
C               BETWEEN TWO SUCCESSIVE INTERPOLATED VALUES IS           DALI 370
C               ABSOLUTELY LESS THAN TOLERANCE EPS, OR IF THE ABSOLUTE  DALI 380
C               VALUE OF THIS DIFFERENCE STOPS DIMINISHING, OR AFTER    DALI 390
C               (NDIM-1) STEPS. FURTHER IT IS TERMINATED IF THE         DALI 400
C               PROCEDURE DISCOVERS TWO ARGUMENT VALUES IN VECTOR ARG   DALI 410
C               WHICH ARE IDENTICAL. DEPENDENT ON THESE FOUR CASES,     DALI 420
C               ERROR PARAMETER IER IS CODED IN THE FOLLOWING FORM      DALI 430
C                IER=0 - IT WAS POSSIBLE TO REACH THE REQUIRED          DALI 440
C                        ACCURACY (NO ERROR).                           DALI 450
C                IER=1 - IT WAS IMPOSSIBLE TO REACH THE REQUIRED        DALI 460
C                        ACCURACY BECAUSE OF ROUNDING ERRORS.           DALI 470
C                IER=2 - IT WAS IMPOSSIBLE TO CHECK ACCURACY BECAUSE    DALI 480
C                        NDIM IS LESS THAN 3, OR THE REQUIRED ACCURACY  DALI 490
C                        COULD NOT BE REACHED BY MEANS OF THE GIVEN     DALI 500
C                        TABLE. NDIM SHOULD BE INCREASED.               DALI 510
C                IER=3 - THE PROCEDURE DISCOVERED TWO ARGUMENT VALUES   DALI 520
C                        IN VECTOR ARG WHICH ARE IDENTICAL.             DALI 530
C                                                                       DALI 540
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DALI 550
C           NONE                                                        DALI 560
C                                                                       DALI 570
C        METHOD                                                         DALI 580
C           INTERPOLATION IS DONE BY MEANS OF AITKENS SCHEME OF         DALI 590
C          LAGRANGE INTERPOLATION. ON RETURN Y CONTAINS AN INTERPOLATED DALI 600
C           FUNCTION VALUE AT POINT X, WHICH IS IN THE SENSE OF REMARK  DALI 610
C           (3) OPTIMAL WITH RESPECT TO GIVEN TABLE. FOR REFERENCE, SEE DALI 620
C           F.B.HILDEBRAND, INTRODUCTION TO NUMERICAL ANALYSIS,         DALI 630
C           MCGRAW-HILL, NEW YORK/TORONTO/LONDON, 1956, PP.49-50.       DALI 640
C                                                                       DALI 650
C     ................................................................. DALI 660
C                                                                       DALI 670
      SUBROUTINE DALI(X,ARG,VAL,Y,NDIM,EPS,IER)                         DALI 680
C                                                                       DALI 690
      DOUBLE PRECISION ARG,VAL,X,Y,H                                    DALI 700
      DIMENSION ARG(1),VAL(1)                                           DALI 710
      IER=2                                                             DALI 720
      DELT2=0.                                                          DALI 730
      IF(NDIM-1)9,7,1                                                   DALI 740
C                                                                       DALI 750
C     START OF AITKEN-LOOP                                              DALI 760
    1 DO 6 J=2,NDIM                                                     DALI 770
      DELT1=DELT2                                                       DALI 780
      IEND=J-1                                                          DALI 790
      DO 2 I=1,IEND                                                     DALI 800
      H=ARG(I)-ARG(J)                                                   DALI 810
      IF(H)2,13,2                                                       DALI 820
    2 VAL(J)=(VAL(I)*(X-ARG(J))-VAL(J)*(X-ARG(I)))/H                    DALI 830
      DELT2=DABS(VAL(J)-VAL(IEND))                                      DALI 840
      IF(J-2)6,6,3                                                      DALI 850
    3 IF(DELT2-EPS)10,10,4                                              DALI 860
    4 IF(J-5)6,5,5                                                      DALI 870
    5 IF(DELT2-DELT1)6,11,11                                            DALI 880
    6 CONTINUE                                                          DALI 890
C     END OF AITKEN-LOOP                                                DALI 900
C                                                                       DALI 910
    7 J=NDIM                                                            DALI 920
    8 Y=VAL(J)                                                          DALI 930
    9 RETURN                                                            DALI 940
C                                                                       DALI 950
C     THERE IS SUFFICIENT ACCURACY WITHIN NDIM-1 ITERATION STEPS        DALI 960
   10 IER=0                                                             DALI 970
      GOTO 8                                                            DALI 980
C                                                                       DALI 990
C     TEST VALUE DELT2 STARTS OSCILLATING                               DALI1000
   11 IER=1                                                             DALI1010
   12 J=IEND                                                            DALI1020
      GOTO 8                                                            DALI1030
C                                                                       DALI1040
C     THERE ARE TWO IDENTICAL ARGUMENT VALUES IN VECTOR ARG             DALI1050
   13 IER=3                                                             DALI1060
      GOTO 12                                                           DALI1070
      END                                                               DALI1080
