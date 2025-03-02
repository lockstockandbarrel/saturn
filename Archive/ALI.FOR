C                                                                       ALI   10
C     ..................................................................ALI   20
C                                                                       ALI   30
C        SUBROUTINE ALI                                                 ALI   40
C                                                                       ALI   50
C        PURPOSE                                                        ALI   60
C           TO INTERPOLATE FUNCTION VALUE Y FOR A GIVEN ARGUMENT VALUE  ALI   70
C           X USING A GIVEN TABLE (ARG,VAL) OF ARGUMENT AND FUNCTION    ALI   80
C           VALUES.                                                     ALI   90
C                                                                       ALI  100
C        USAGE                                                          ALI  110
C           CALL ALI (X,ARG,VAL,Y,NDIM,EPS,IER)                         ALI  120
C                                                                       ALI  130
C        DESCRIPTION OF PARAMETERS                                      ALI  140
C           X      - THE ARGUMENT VALUE SPECIFIED BY INPUT.             ALI  150
C           ARG    - THE INPUT VECTOR (DIMENSION NDIM) OF ARGUMENT      ALI  160
C                    VALUES OF THE TABLE (NOT DESTROYED).               ALI  170
C           VAL    - THE INPUT VECTOR (DIMENSION NDIM) OF FUNCTION      ALI  180
C                    VALUES OF THE TABLE (DESTROYED).                   ALI  190
C           Y      - THE RESULTING INTERPOLATED FUNCTION VALUE.         ALI  200
C           NDIM   - AN INPUT VALUE WHICH SPECIFIES THE NUMBER OF       ALI  210
C                    POINTS IN TABLE (ARG,VAL).                         ALI  220
C           EPS    - AN INPUT CONSTANT WHICH IS USED AS UPPER BOUND     ALI  230
C                    FOR THE ABSOLUTE ERROR.                            ALI  240
C           IER    - A RESULTING ERROR PARAMETER.                       ALI  250
C                                                                       ALI  260
C        REMARKS                                                        ALI  270
C           (1) TABLE (ARG,VAL) SHOULD REPRESENT A SINGLE-VALUED        ALI  280
C               FUNCTION AND SHOULD BE STORED IN SUCH A WAY, THAT THE   ALI  290
C               DISTANCES ABS(ARG(I)-X) INCREASE WITH INCREASING        ALI  300
C               SUBSCRIPT I. TO GENERATE THIS ORDER IN TABLE (ARG,VAL), ALI  310
C               SUBROUTINES ATSG, ATSM OR ATSE COULD BE USED IN A       ALI  320
C               PREVIOUS STAGE.                                         ALI  330
C           (2) NO ACTION BESIDES ERROR MESSAGE IN CASE NDIM LESS       ALI  340
C               THAN 1.                                                 ALI  350
C           (3) INTERPOLATION IS TERMINATED EITHER IF THE DIFFERENCE    ALI  360
C               BETWEEN TWO SUCCESSIVE INTERPOLATED VALUES IS           ALI  370
C               ABSOLUTELY LESS THAN TOLERANCE EPS, OR IF THE ABSOLUTE  ALI  380
C               VALUE OF THIS DIFFERENCE STOPS DIMINISHING, OR AFTER    ALI  390
C               (NDIM-1) STEPS. FURTHER IT IS TERMINATED IF THE         ALI  400
C               PROCEDURE DISCOVERS TWO ARGUMENT VALUES IN VECTOR ARG   ALI  410
C               WHICH ARE IDENTICAL. DEPENDENT ON THESE FOUR CASES,     ALI  420
C               ERROR PARAMETER IER IS CODED IN THE FOLLOWING FORM      ALI  430
C                IER=0 - IT WAS POSSIBLE TO REACH THE REQUIRED          ALI  440
C                        ACCURACY (NO ERROR).                           ALI  450
C                IER=1 - IT WAS IMPOSSIBLE TO REACH THE REQUIRED        ALI  460
C                        ACCURACY BECAUSE OF ROUNDING ERRORS.           ALI  470
C                IER=2 - IT WAS IMPOSSIBLE TO CHECK ACCURACY BECAUSE    ALI  480
C                        NDIM IS LESS THAN 3, OR THE REQUIRED ACCURACY  ALI  490
C                        COULD NOT BE REACHED BY MEANS OF THE GIVEN     ALI  500
C                        TABLE. NDIM SHOULD BE INCREASED.               ALI  510
C                IER=3 - THE PROCEDURE DISCOVERED TWO ARGUMENT VALUES   ALI  520
C                        IN VECTOR ARG WHICH ARE IDENTICAL.             ALI  530
C                                                                       ALI  540
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  ALI  550
C           NONE                                                        ALI  560
C                                                                       ALI  570
C        METHOD                                                         ALI  580
C           INTERPOLATION IS DONE BY MEANS OF AITKENS SCHEME OF         ALI  590
C           LAGRANGE INTERPOLATION. ON RETURN Y CONTAINS AN INTERPOLATEDALI  600
C           FUNCTION VALUE AT POINT X, WHICH IS IN THE SENSE OF REMARK  ALI  610
C           (3) OPTIMAL WITH RESPECT TO GIVEN TABLE. FOR REFERENCE, SEE ALI  620
C           F.B.HILDEBRAND, INTRODUCTION TO NUMERICAL ANALYSIS,         ALI  630
C           MCGRAW-HILL, NEW YORK/TORONTO/LONDON, 1956, PP.49-50.       ALI  640
C                                                                       ALI  650
C     ..................................................................ALI  660
C                                                                       ALI  670
      SUBROUTINE ALI(X,ARG,VAL,Y,NDIM,EPS,IER)                          ALI  680
C                                                                       ALI  690
C                                                                       ALI  700
      DIMENSION ARG(1),VAL(1)                                           ALI  710
      IER=2                                                             ALI  720
      DELT2=0.                                                          ALI  730
      IF(NDIM-1)9,7,1                                                   ALI  740
C                                                                       ALI  750
C     START OF AITKEN-LOOP                                              ALI  760
    1 DO 6 J=2,NDIM                                                     ALI  770
      DELT1=DELT2                                                       ALI  780
      IEND=J-1                                                          ALI  790
      DO 2 I=1,IEND                                                     ALI  800
      H=ARG(I)-ARG(J)                                                   ALI  810
      IF(H)2,13,2                                                       ALI  820
    2 VAL(J)=(VAL(I)*(X-ARG(J))-VAL(J)*(X-ARG(I)))/H                    ALI  830
      DELT2=ABS(VAL(J)-VAL(IEND))                                       ALI  840
      IF(J-2)6,6,3                                                      ALI  850
    3 IF(DELT2-EPS)10,10,4                                              ALI  860
    4 IF(J-5)6,5,5                                                      ALI  870
    5 IF(DELT2-DELT1)6,11,11                                            ALI  880
    6 CONTINUE                                                          ALI  890
C     END OF AITKEN-LOOP                                                ALI  900
C                                                                       ALI  910
    7 J=NDIM                                                            ALI  920
    8 Y=VAL(J)                                                          ALI  930
    9 RETURN                                                            ALI  940
C                                                                       ALI  950
C     THERE IS SUFFICIENT ACCURACY WITHIN NDIM-1 ITERATION STEPS        ALI  960
   10 IER=0                                                             ALI  970
      GOTO 8                                                            ALI  980
C                                                                       ALI  990
C     TEST VALUE DELT2 STARTS OSCILLATING                               ALI 1000
   11 IER=1                                                             ALI 1010
   12 J=IEND                                                            ALI 1020
      GOTO 8                                                            ALI 1030
C                                                                       ALI 1040
C     THERE ARE TWO IDENTICAL ARGUMENT VALUES IN VECTOR ARG             ALI 1050
   13 IER=3                                                             ALI 1060
      GOTO 12                                                           ALI 1070
      END                                                               ALI 1080
