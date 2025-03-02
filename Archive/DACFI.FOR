C                                                                       DACF  10
C     ..................................................................DACF  20
C                                                                       DACF  30
C        SUBROUTINE DACFI                                               DACF  40
C                                                                       DACF  50
C        PURPOSE                                                        DACF  60
C           TO INTERPOLATE FUNCTION VALUE Y FOR A GIVEN ARGUMENT VALUE  DACF  70
C           X USING A GIVEN TABLE (ARG,VAL) OF ARGUMENT AND FUNCTION    DACF  80
C           VALUES.                                                     DACF  90
C                                                                       DACF 100
C        USAGE                                                          DACF 110
C           CALL DACFI (X,ARG,VAL,Y,NDIM,EPS,IER)                       DACF 120
C                                                                       DACF 130
C        DESCRIPTION OF PARAMETERS                                      DACF 140
C           X      - THE ARGUMENT VALUE SPECIFIED BY INPUT.             DACF 150
C           ARG    - THE INPUT VECTOR (DIMENSION NDIM) OF ARGUMENT      DACF 160
C                    VALUES OF THE TABLE (POSSIBLY DESTROYED).          DACF 170
C           VAL    - THE INPUT VECTOR (DIMENSION NDIM) OF FUNCTION      DACF 180
C                    VALUES OF THE TABLE (DESTROYED).                   DACF 190
C           Y      - THE RESULTING INTERPOLATED FUNCTION VALUE.         DACF 200
C           NDIM   - AN INPUT VALUE WHICH SPECIFIES THE NUMBER OF       DACF 210
C                    POINTS IN TABLE (ARG,VAL).                         DACF 220
C           EPS    - AN INPUT CONSTANT WHICH IS USED AS UPPER BOUND     DACF 230
C                    FOR THE ABSOLUTE ERROR.                            DACF 240
C           IER    - A RESULTING ERROR PARAMETER.                       DACF 250
C        NOTE: ARG,VAL,X,Y,Z,P1,P2,P3,Q1,Q2,Q3,AUX,H ARE DOUBLE PREC.   DACF 260
C        REMARKS                                                        DACF 270
C           (1) TABLE (ARG,VAL) SHOULD REPRESENT A SINGLE-VALUED        DACF 280
C               FUNCTION AND SHOULD BE STORED IN SUCH A WAY, THAT THE   DACF 290
C               DISTANCES ABS(ARG(I)-X) INCREASE WITH INCREASING        DACF 300
C               SUBSCRIPT I. TO GENERATE THIS ORDER IN TABLE (ARG,VAL), DACF 310
C               SUBROUTINES ATSG, ATSM OR ATSE COULD BE USED IN A       DACF 320
C               PREVIOUS STAGE.                                         DACF 330
C           (2) NO ACTION BESIDES ERROR MESSAGE IN CASE NDIM LESS       DACF 340
C               THAN 1.                                                 DACF 350
C           (3) INTERPOLATION IS TERMINATED EITHER IF THE DIFFERENCE    DACF 360
C               BETWEEN TWO SUCCESSIVE INTERPOLATED VALUES IS           DACF 370
C               ABSOLUTELY LESS THAN TOLERANCE EPS, OR IF THE ABSOLUTE  DACF 380
C               VALUE OF THIS DIFFERENCE STOPS DIMINISHING, OR AFTER    DACF 390
C               (NDIM-1) STEPS (THE NUMBER OF POSSIBLE STEPS IS         DACF 400
C               DIMINISHED IF AT ANY STAGE INFINITY ELEMENT APPEARS IN  DACF 410
C               THE DOWNWARD DIAGONAL OF INVERTED-DIFFERENCES-SCHEME    DACF 420
C               AND IF IT IS IMPOSSIBLE TO ELIMINATE THIS INFINITY      DACF 430
C               ELEMENT BY INTERCHANGING OF TABLE POINTS).              DACF 440
C               FURTHER IT IS TERMINATED IF THE PROCEDURE DISCOVERS TWO DACF 450
C               ARGUMENT VALUES IN VECTOR ARG WHICH ARE IDENTICAL.      DACF 460
C               DEPENDENT ON THESE FOUR CASES, ERROR PARAMETER IER IS   DACF 470
C               CODED IN THE FOLLOWING FORM                             DACF 480
C                IER=0 - IT WAS POSSIBLE TO REACH THE REQUIRED          DACF 490
C                        ACCURACY (NO ERROR).                           DACF 500
C                IER=1 - IT WAS IMPOSSIBLE TO REACH THE REQUIRED        DACF 510
C                        ACCURACY BECAUSE OF ROUNDING ERRORS.           DACF 520
C                IER=2 - IT WAS IMPOSSIBLE TO CHECK ACCURACY BECAUSE    DACF 530
C                        NDIM IS LESS THAN 2, OR THE REQUIRED ACCURACY  DACF 540
C                        COULD NOT BE REACHED BY MEANS OF THE GIVEN     DACF 550
C                        TABLE. NDIM SHOULD BE INCREASED.               DACF 560
C                IER=3 - THE PROCEDURE DISCOVERED TWO ARGUMENT VALUES   DACF 570
C                        IN VECTOR ARG WHICH ARE IDENTICAL.             DACF 580
C                                                                       DACF 590
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DACF 600
C           NONE                                                        DACF 610
C                                                                       DACF 620
C        METHOD                                                         DACF 630
C           INTERPOLATION IS DONE BY CONTINUED FRACTIONS AND INVERTED-  DACF 640
C           DIFFERENCES-SCHEME. ON RETURN Y CONTAINS AN INTERPOLATED    DACF 650
C           FUNCTION VALUE AT POINT X, WHICH IS IN THE SENSE OF REMARK  DACF 660
C           (3) OPTIMAL WITH RESPECT TO GIVEN TABLE. FOR REFERENCE, SEE DACF 670
C           F.B.HILDEBRAND, INTRODUCTION TO NUMERICAL ANALYSIS,         DACF 680
C           MCGRAW-HILL, NEW YORK/TORONTO/LONDON, 1956, PP.395-406.     DACF 690
C                                                                       DACF 700
C     ..................................................................DACF 710
C                                                                       DACF 720
      SUBROUTINE DACFI(X,ARG,VAL,Y,NDIM,EPS,IER)                        DACF 730
C                                                                       DACF 740
      DOUBLE PRECISION ARG,VAL,X,Y,Z,P1,P2,P3,Q1,Q2,Q3,AUX,H            DACF 750
      DIMENSION ARG(1),VAL(1)                                           DACF 760
      IER=2                                                             DACF 770
      IF(NDIM)20,20,1                                                   DACF 780
    1 Y=VAL(1)                                                          DACF 790
      DELT2=0.                                                          DACF 800
      IF(NDIM-1)20,20,2                                                 DACF 810
C                                                                       DACF 820
C     PREPARATIONS FOR INTERPOLATION LOOP                               DACF 830
    2 P2=1.0D0                                                          DACF 840
      P3=Y                                                              DACF 850
      Q2=0.0D0                                                          DACF 860
      Q3=1.0D0                                                          DACF 870
C                                                                       DACF 880
C                                                                       DACF 890
C     START INTERPOLATION LOOP                                          DACF 900
      DO 16 I=2,NDIM                                                    DACF 910
      II=0                                                              DACF 920
      P1=P2                                                             DACF 930
      P2=P3                                                             DACF 940
      Q1=Q2                                                             DACF 950
      Q2=Q3                                                             DACF 960
      Z=Y                                                               DACF 970
      DELT1=DELT2                                                       DACF 980
      JEND=I-1                                                          DACF 990
C                                                                       DACF1000
C     COMPUTATION OF INVERTED DIFFERENCES                               DACF1010
    3 AUX=VAL(I)                                                        DACF1020
      DO 10 J=1,JEND                                                    DACF1030
      H=VAL(I)-VAL(J)                                                   DACF1040
      IF(DABS(H)-1.E-6*DABS(VAL(I)))4,4,9                               DACF1050
    4 IF(ARG(I)-ARG(J))5,17,5                                           DACF1060
    5 IF(J-JEND)8,6,6                                                   DACF1070
C                                                                       DACF1080
C     INTERCHANGE ROW I WITH ROW I+II                                   DACF1090
    6 II=II+1                                                           DACF1100
      III=I+II                                                          DACF1110
      IF(III-NDIM)7,7,19                                                DACF1120
    7 VAL(I)=VAL(III)                                                   DACF1130
      VAL(III)=AUX                                                      DACF1140
      AUX=ARG(I)                                                        DACF1150
      ARG(I)=ARG(III)                                                   DACF1160
      ARG(III)=AUX                                                      DACF1170
      GOTO 3                                                            DACF1180
C                                                                       DACF1190
C     COMPUTATION OF VAL(I) IN CASE VAL(I)=VAL(J) AND J LESS THAN I-1   DACF1200
    8 VAL(I)=1.D75                                                      DACF1210
      GOTO 10                                                           DACF1220
C                                                                       DACF1230
C     COMPUTATION OF VAL(I) IN CASE VAL(I) NOT EQUAL TO VAL(J)          DACF1240
    9 VAL(I)=(ARG(I)-ARG(J))/H                                          DACF1250
   10 CONTINUE                                                          DACF1260
C     INVERTED DIFFERENCES ARE COMPUTED                                 DACF1270
C                                                                       DACF1280
C     COMPUTATION OF NEW Y                                              DACF1290
      P3=VAL(I)*P2+(X-ARG(I-1))*P1                                      DACF1300
      Q3=VAL(I)*Q2+(X-ARG(I-1))*Q1                                      DACF1310
      IF(Q3)11,12,11                                                    DACF1320
   11 Y=P3/Q3                                                           DACF1330
      GOTO 13                                                           DACF1340
   12 Y=1.D75                                                           DACF1350
   13 DELT2=DABS(Z-Y)                                                   DACFI1360
      IF(DELT2-EPS)19,19,14                                             DACF1370
   14 IF(I-8)16,15,15                                                   DACF1380
   15 IF(DELT2-DELT1)16,18,18                                           DACF1390
   16 CONTINUE                                                          DACF1400
C     END OF INTERPOLATION LOOP                                         DACF1410
C                                                                       DACF1420
C                                                                       DACF1430
      RETURN                                                            DACF1440
C                                                                       DACF1450
C     THERE ARE TWO IDENTICAL ARGUMENT VALUES IN VECTOR ARG             DACF1460
   17 IER=3                                                             DACF1470
      RETURN                                                            DACF1480
C                                                                       DACF1490
C     TEST VALUE DELT2 STARTS OSCILLATING                               DACF1500
   18 Y=Z                                                               DACF1510
      IER=1                                                             DACF1520
      RETURN                                                            DACF1530
C                                                                       DACF1540
C     THERE IS SATISFACTORY ACCURACY WITHIN NDIM-1 STEPS                DACF1550
   19 IER=0                                                             DACF1560
   20 RETURN                                                            DACF1570
      END                                                               DACF1580
