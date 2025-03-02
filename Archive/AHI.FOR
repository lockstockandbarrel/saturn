C                                                                       AHI   10
C     ..................................................................AHI   20
C                                                                       AHI   30
C        SUBROUTINE AHI                                                 AHI   40
C                                                                       AHI   50
C        PURPOSE                                                        AHI   60
C           TO INTERPOLATE FUNCTION VALUE Y FOR A GIVEN ARGUMENT VALUE  AHI   70
C           X USING A GIVEN TABLE (ARG,VAL) OF ARGUMENT, FUNCTION, AND  AHI   80
C           DERIVATIVE VALUES.                                          AHI   90
C                                                                       AHI  100
C        USAGE                                                          AHI  110
C           CALL AHI (X,ARG,VAL,Y,NDIM,EPS,IER)                         AHI  120
C                                                                       AHI  130
C        DESCRIPTION OF PARAMETERS                                      AHI  140
C           X      - THE ARGUMENT VALUE SPECIFIED BY INPUT.             AHI  150
C           ARG    - THE INPUT VECTOR (DIMENSION NDIM) OF ARGUMENT      AHI  160
C                    VALUES OF THE TABLE (NOT DESTROYED).               AHI  170
C           VAL    - THE INPUT VECTOR (DIMENSION 2*NDIM) OF FUNCTION    AHI  180
C                    AND DERIVATIVE VALUES OF THE TABLE (DESTROYED).    AHI  190
C                    FUNCTION AND DERIVATIVE VALUES MUST BE STORED IN   AHI  200
C                    PAIRS, THAT MEANS BEGINNING WITH FUNCTION VALUE AT AHI  210
C                    POINT ARG(1) EVERY FUNCTION VALUE MUST BE FOLLOWED AHI  220
C                    BY THE VALUE OF DERIVATIVE AT THE SAME POINT.      AHI  230
C           Y      - THE RESULTING INTERPOLATED FUNCTION VALUE.         AHI  240
C           NDIM   - AN INPUT VALUE WHICH SPECIFIES THE NUMBER OF       AHI  250
C                    POINTS IN TABLE (ARG,VAL).                         AHI  260
C           EPS    - AN INPUT CONSTANT WHICH IS USED AS UPPER BOUND     AHI  270
C                    FOR THE ABSOLUTE ERROR.                            AHI  280
C           IER    - A RESULTING ERROR PARAMETER.                       AHI  290
C                                                                       AHI  300
C        REMARKS                                                        AHI  310
C           (1) TABLE (ARG,VAL) SHOULD REPRESENT A SINGLE-VALUED        AHI  320
C               FUNCTION AND SHOULD BE STORED IN SUCH A WAY, THAT THE   AHI  330
C               DISTANCES ABS(ARG(I)-X) INCREASE WITH INCREASING        AHI  340
C               SUBSCRIPT I. TO GENERATE THIS ORDER IN TABLE (ARG,VAL), AHI  350
C               SUBROUTINES ATSG, ATSM OR ATSE COULD BE USED IN A       AHI  360
C               PREVIOUS STAGE.                                         AHI  370
C           (2) NO ACTION BESIDES ERROR MESSAGE IN CASE NDIM LESS       AHI  380
C               THAN 1.                                                 AHI  390
C           (3) INTERPOLATION IS TERMINATED EITHER IF THE DIFFERENCE    AHI  400
C               BETWEEN TWO SUCCESSIVE INTERPOLATED VALUES IS           AHI  410
C               ABSOLUTELY LESS THAN TOLERANCE EPS, OR IF THE ABSOLUTE  AHI  420
C               VALUE OF THIS DIFFERENCE STOPS DIMINISHING, OR AFTER    AHI  430
C               (2*NDIM-2) STEPS. FURTHER IT IS TERMINATED IF THE       AHI  440
C               PROCEDURE DISCOVERS TWO ARGUMENT VALUES IN VECTOR ARG   AHI  450
C               WHICH ARE IDENTICAL. DEPENDENT ON THESE FOUR CASES,     AHI  460
C               ERROR PARAMETER IER IS CODED IN THE FOLLOWING FORM      AHI  470
C                IER=0 - IT WAS POSSIBLE TO REACH THE REQUIRED          AHI  480
C                        ACCURACY (NO ERROR).                           AHI  490
C                IER=1 - IT WAS IMPOSSIBLE TO REACH THE REQUIRED        AHI  500
C                        ACCURACY BECAUSE OF ROUNDING ERRORS.           AHI  510
C                IER=2 - IT WAS IMPOSSIBLE TO CHECK ACCURACY BECAUSE    AHI  520
C                        NDIM IS LESS THAN 2, OR THE REQUIRED ACCURACY  AHI  530
C                        COULD NOT BE REACHED BY MEANS OF THE GIVEN     AHI  540
C                        TABLE. NDIM SHOULD BE INCREASED.               AHI  550
C                IER=3 - THE PROCEDURE DISCOVERED TWO ARGUMENT VALUES   AHI  560
C                        IN VECTOR ARG WHICH ARE IDENTICAL.             AHI  570
C                                                                       AHI  580
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  AHI  590
C           NONE                                                        AHI  600
C                                                                       AHI  610
C        METHOD                                                         AHI  620
C           INTERPOLATION IS DONE BY MEANS OF AITKENS SCHEME OF         AHI  630
C           HERMITE INTERPOLATION. ON RETURN Y CONTAINS AN INTERPOLATED AHI  640
C           FUNCTION VALUE AT POINT X, WHICH IS IN THE SENSE OF REMARK  AHI  650
C           (3) OPTIMAL WITH RESPECT TO GIVEN TABLE. FOR REFERENCE, SEE AHI  660
C           F.B.HILDEBRAND, INTRODUCTION TO NUMERICAL ANALYSIS,         AHI  670
C           MCGRAW-HILL, NEW YORK/TORONTO/LONDON, 1956, PP.314-317, AND AHI  680
C           GERSHINSKY/LEVINE, AITKEN-HERMITE INTERPOLATION,            AHI  690
C           JACM, VOL.11, ISS.3 (1964), PP.352-356.                     AHI  700
C                                                                       AHI  710
C     ..................................................................AHI  720
C                                                                       AHI  730
      SUBROUTINE AHI(X,ARG,VAL,Y,NDIM,EPS,IER)                          AHI  740
C                                                                       AHI  750
C                                                                       AHI  760
      DIMENSION ARG(1),VAL(1)                                           AHI  770
      IER=2                                                             AHI  780
      H2=X-ARG(1)                                                       AHI  790
      IF(NDIM-1)2,1,3                                                   AHI  800
    1 Y=VAL(1)+VAL(2)*H2                                                AHI  810
    2 RETURN                                                            AHI  820
C                                                                       AHI  830
C     VECTOR ARG HAS MORE THAN 1 ELEMENT.                               AHI  840
C     THE FIRST STEP PREPARES VECTOR VAL SUCH THAT AITKEN SCHEME CAN BE AHI  850
C     USED.                                                             AHI  860
    3 I=1                                                               AHI  870
      DO 5 J=2,NDIM                                                     AHI  880
      H1=H2                                                             AHI  890
      H2=X-ARG(J)                                                       AHI  900
      Y=VAL(I)                                                          AHI  910
      VAL(I)=Y+VAL(I+1)*H1                                              AHI  920
      H=H1-H2                                                           AHI  930
      IF(H)4,13,4                                                       AHI  940
    4 VAL(I+1)=Y+(VAL(I+2)-Y)*H1/H                                      AHI  950
    5 I=I+2                                                             AHI  960
      VAL(I)=VAL(I)+VAL(I+1)*H2                                         AHI  970
C     END OF FIRST STEP                                                 AHI  980
C                                                                       AHI  990
C     PREPARE AITKEN SCHEME                                             AHI 1000
      DELT2=0.                                                          AHI 1010
      IEND=I-1                                                          AHI 1020
C                                                                       AHI 1030
C     START AITKEN-LOOP                                                 AHI 1040
      DO 9 I=1,IEND                                                     AHI 1050
      DELT1=DELT2                                                       AHI 1060
      Y=VAL(1)                                                          AHI 1070
      M=(I+3)/2                                                         AHI 1080
      H1=ARG(M)                                                         AHI 1090
      DO 6 J=1,I                                                        AHI 1100
      K=I+1-J                                                           AHI 1110
      L=(K+1)/2                                                         AHI 1120
      H=ARG(L)-H1                                                       AHI 1130
      IF(H)6,14,6                                                       AHI 1140
    6 VAL(K)=(VAL(K)*(X-H1)-VAL(K+1)*(X-ARG(L)))/H                      AHI 1150
      DELT2=ABS(Y-VAL(1))                                               AHI 1160
      IF(DELT2-EPS)11,11,7                                              AHI 1170
    7 IF(I-5)9,8,8                                                      AHI 1180
    8 IF(DELT2-DELT1)9,12,12                                            AHI 1190
    9 CONTINUE                                                          AHI 1200
C     END OF AITKEN-LOOP                                                AHI 1210
C                                                                       AHI 1220
   10 Y=VAL(1)                                                          AHI 1230
      RETURN                                                            AHI 1240
C                                                                       AHI 1250
C     THERE IS SUFFICIENT ACCURACY WITHIN 2*NDIM-2 ITERATION STEPS      AHI 1260
   11 IER=0                                                             AHI 1270
      GOTO 10                                                           AHI 1280
C                                                                       AHI 1290
C     TEST VALUE DELT2 STARTS OSCILLATING                               AHI 1300
   12 IER=1                                                             AHI 1310
      RETURN                                                            AHI 1320
C                                                                       AHI 1330
C     THERE ARE TWO IDENTICAL ARGUMENT VALUES IN VECTOR ARG             AHI 1340
   13 Y=VAL(1)                                                          AHI 1350
   14 IER=3                                                             AHI 1360
      RETURN                                                            AHI 1370
      END                                                               AHI 1380
