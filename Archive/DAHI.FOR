C                                                                       DAHI  10
C     ..................................................................DAHI  20
C                                                                       DAHI  30
C        SUBROUTINE DAHI                                                DAHI  40
C                                                                       DAHI  50
C        PURPOSE                                                        DAHI  60
C           TO INTERPOLATE FUNCTION VALUE Y FOR A GIVEN ARGUMENT VALUE  DAHI  70
C           X USING A GIVEN TABLE (ARG,VAL) OF ARGUMENT, FUNCTION, AND  DAHI  80
C           DERIVATIVE VALUES.                                          DAHI  90
C                                                                       DAHI 100
C        USAGE                                                          DAHI 110
C           CALL DAHI (X,ARG,VAL,Y,NDIM,EPS,IER)                        DAHI 120
C                                                                       DAHI 130
C        DESCRIPTION OF PARAMETERS                                      DAHI 140
C           X      - THE ARGUMENT VALUE SPECIFIED BY INPUT.             DAHI 150
C           ARG    - THE INPUT VECTOR (DIMENSION NDIM) OF ARGUMENT      DAHI 160
C                    VALUES OF THE TABLE (NOT DESTROYED).               DAHI 170
C           VAL    - THE INPUT VECTOR (DIMENSION 2*NDIM) OF FUNCTION    DAHI 180
C                    AND DERIVATIVE VALUES OF THE TABLE (DESTROYED).    DAHI 190
C                    FUNCTION AND DERIVATIVE VALUES MUST BE STORED IN   DAHI 200
C                    PAIRS, THAT MEANS BEGINNING WITH FUNCTION VALUE AT DAHI 210
C                    POINT ARG(1) EVERY FUNCTION VALUE MUST BE FOLLOWED DAHI 220
C                    BY THE VALUE OF DERIVATIVE AT THE SAME POINT.      DAHI 230
C           Y      - THE RESULTING INTERPOLATED FUNCTION VALUE.         DAHI 240
C           NDIM   - AN INPUT VALUE WHICH SPECIFIES THE NUMBER OF       DAHI 250
C                    POINTS IN TABLE (ARG,VAL).                         DAHI 260
C           EPS    - AN INPUT CONSTANT WHICH IS USED AS UPPER BOUND     DAHI 270
C                    FOR THE ABSOLUTE ERROR.                            DAHI 280
C           IER    - A RESULTING ERROR PARAMETER.                       DAHI 290
C           NOTE: ARG,VAL,X,Y,H,H1 AND H2 ARE DOUBLE PRECISION          DAHI 300
C        REMARKS                                                        DAHI 310
C           (1) TABLE (ARG,VAL) SHOULD REPRESENT A SINGLE-VALUED        DAHI 320
C               FUNCTION AND SHOULD BE STORED IN SUCH A WAY, THAT THE   DAHI 330
C               DISTANCES ABS(ARG(I)-X) INCREASE WITH INCREASING        DAHI 340
C               SUBSCRIPT I. TO GENERATE THIS ORDER IN TABLE (ARG,VAL), DAHI 350
C               SUBROUTINES ATSG, ATSM OR ATSE COULD BE USED IN A       DAHI 360
C               PREVIOUS STAGE.                                         DAHI 370
C           (2) NO ACTION BESIDES ERROR MESSAGE IN CASE NDIM LESS       DAHI 380
C               THAN 1.                                                 DAHI 390
C           (3) INTERPOLATION IS TERMINATED EITHER IF THE DIFFERENCE    DAHI 400
C               BETWEEN TWO SUCCESSIVE INTERPOLATED VALUES IS           DAHI 410
C               ABSOLUTELY LESS THAN TOLERANCE EPS, OR IF THE ABSOLUTE  DAHI 420
C               VALUE OF THIS DIFFERENCE STOPS DIMINISHING, OR AFTER    DAHI 430
C               (2*NDIM-2) STEPS. FURTHER IT IS TERMINATED IF THE       DAHI 440
C               PROCEDURE DISCOVERS TWO ARGUMENT VALUES IN VECTOR ARG   DAHI 450
C               WHICH ARE IDENTICAL. DEPENDENT ON THESE FOUR CASES,     DAHI 460
C               ERROR PARAMETER IER IS CODED IN THE FOLLOWING FORM      DAHI 470
C                IER=0 - IT WAS POSSIBLE TO REACH THE REQUIRED          DAHI 480
C                        ACCURACY (NO ERROR).                           DAHI 490
C                IER=1 - IT WAS IMPOSSIBLE TO REACH THE REQUIRED        DAHI 500
C                        ACCURACY BECAUSE OF ROUNDING ERRORS.           DAHI 510
C                IER=2 - IT WAS IMPOSSIBLE TO CHECK ACCURACY BECAUSE    DAHI 520
C                        NDIM IS LESS THAN 2, OR THE REQUIRED ACCURACY  DAHI 530
C                        COULD NOT BE REACHED BY MEANS OF THE GIVEN     DAHI 540
C                        TABLE. NDIM SHOULD BE INCREASED.               DAHI 550
C                IER=3 - THE PROCEDURE DISCOVERED TWO ARGUMENT VALUES   DAHI 560
C                        IN VECTOR ARG WHICH ARE IDENTICAL.             DAHI 570
C                                                                       DAHI 580
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DAHI 590
C           NONE                                                        DAHI 600
C                                                                       DAHI 610
C        METHOD                                                         DAHI 620
C           INTERPOLATION IS DONE BY MEANS OF AITKENS SCHEME OF         DAHI 630
C           HERMITE INTERPOLATION. ON RETURN Y CONTAINS AN INTERPOLATED DAHI 640
C           FUNCTION VALUE AT POINT X, WHICH IS IN THE SENSE OF REMARK  DAHI 650
C           (3) OPTIMAL WITH RESPECT TO GIVEN TABLE. FOR REFERENCE, SEE DAHI 660
C           F.B.HILDEBRAND, INTRODUCTION TO NUMERICAL ANALYSIS,         DAHI 670
C           MCGRAW-HILL, NEW YORK/TORONTO/LONDON, 1956, PP.314-317, AND DAHI 680
C           GERSHINSKY/LEVINE, AITKEN-HERMITE INTERPOLATION,            DAHI 690
C           JACM, VOL.11, ISS.3 (1964), PP.352-356.                     DAHI 700
C                                                                       DAHI 710
C     ..................................................................DAHI 720
C                                                                       DAHI 730
      SUBROUTINE DAHI(X,ARG,VAL,Y,NDIM,EPS,IER)                         DAHI 740
C                                                                       DAHI 750
      DOUBLE PRECISION ARG,VAL,X,Y,H,H1,H2                              DAHI 760
      DIMENSION ARG(1),VAL(1)                                           DAHI 770
      IER=2                                                             DAHI 780
      H2=X-ARG(1)                                                       DAHI 790
      IF(NDIM-1)2,1,3                                                   DAHI 800
    1 Y=VAL(1)+VAL(2)*H2                                                DAHI 810
    2 RETURN                                                            DAHI 820
C                                                                       DAHI 830
C     VECTOR ARG HAS MORE THAN 1 ELEMENT.                               DAHI 840
C     THE FIRST STEP PREPARES VECTOR VAL SUCH THAT AITKEN SCHEME CAN BE DAHI 850
C     USED.                                                             DAHI 860
    3 I=1                                                               DAHI 870
      DO 5 J=2,NDIM                                                     DAHI 880
      H1=H2                                                             DAHI 890
      H2=X-ARG(J)                                                       DAHI 900
      Y=VAL(I)                                                          DAHI 910
      VAL(I)=Y+VAL(I+1)*H1                                              DAHI 920
      H=H1-H2                                                           DAHI 930
      IF(H)4,13,4                                                       DAHI 940
    4 VAL(I+1)=Y+(VAL(I+2)-Y)*H1/H                                      DAHI 950
    5 I=I+2                                                             DAHI 960
      VAL(I)=VAL(I)+VAL(I+1)*H2                                         DAHI 970
C     END OF FIRST STEP                                                 DAHI 980
C                                                                       DAHI 990
C     PREPARE AITKEN SCHEME                                             DAHI1000
      DELT2=0.                                                          DAHI1010
      IEND=I-1                                                          DAHI1020
C                                                                       DAHI1030
C     START AITKEN-LOOP                                                 DAHI1040
      DO 9 I=1,IEND                                                     DAHI1050
      DELT1=DELT2                                                       DAHI1060
      Y=VAL(1)                                                          DAHI1070
      M=(I+3)/2                                                         DAHI1080
      H1=ARG(M)                                                         DAHI1090
      DO 6 J=1,I                                                        DAHI1100
      K=I+1-J                                                           DAHI1110
      L=(K+1)/2                                                         DAHI1120
      H=ARG(L)-H1                                                       DAHI1130
      IF(H)6,14,6                                                       DAHI1140
    6 VAL(K)=(VAL(K)*(X-H1)-VAL(K+1)*(X-ARG(L)))/H                      DAHI1150
      DELT2=DABS(Y-VAL(1))                                              DAHI1160
      IF(DELT2-EPS)11,11,7                                              DAHI1170
    7 IF(I-5)9,8,8                                                      DAHI1180
    8 IF(DELT2-DELT1)9,12,12                                            DAHI1190
    9 CONTINUE                                                          DAHI1200
C     END OF AITKEN-LOOP                                                DAHI1210
C                                                                       DAHI1220
   10 Y=VAL(1)                                                          DAHI1230
      RETURN                                                            DAHI1240
C                                                                       DAHI1250
C     THERE IS SUFFICIENT ACCURACY WITHIN 2*NDIM-2 ITERATION STEPS      DAHI1260
   11 IER=0                                                             DAHI1270
      GOTO 10                                                           DAHI1280
C                                                                       DAHI1290
C     TEST VALUE DELT2 STARTS OSCILLATING                               DAHI1300
   12 IER=1                                                             DAHI1310
      RETURN                                                            DAHI1320
C                                                                       DAHI1330
C     THERE ARE TWO IDENTICAL ARGUMENT VALUES IN VECTOR ARG             DAHI1340
   13 Y=VAL(1)                                                          DAHI1350
   14 IER=3                                                             DAHI1360
      RETURN                                                            DAHI1370
      END                                                               DAHI1380
