C                                                                       ACFI  10
C     ..................................................................ACFI  20
C                                                                       ACFI  30
C        SUBROUTINE ACFI                                                ACFI  40
C                                                                       ACFI  50
C        PURPOSE                                                        ACFI  60
C           TO INTERPOLATE FUNCTION VALUE Y FOR A GIVEN ARGUMENT VALUE  ACFI  70
C           X USING A GIVEN TABLE (ARG,VAL) OF ARGUMENT AND FUNCTION    ACFI  80
C           VALUES.                                                     ACFI  90
C                                                                       ACFI 100
C        USAGE                                                          ACFI 110
C           CALL ACFI (X,ARG,VAL,Y,NDIM,EPS,IER)                        ACFI 120
C                                                                       ACFI 130
C        DESCRIPTION OF PARAMETERS                                      ACFI 140
C           X      - THE ARGUMENT VALUE SPECIFIED BY INPUT.             ACFI 150
C           ARG    - THE INPUT VECTOR (DIMENSION NDIM) OF ARGUMENT      ACFI 160
C                    VALUES OF THE TABLE (POSSIBLY DESTROYED).          ACFI 170
C           VAL    - THE INPUT VECTOR (DIMENSION NDIM) OF FUNCTION      ACFI 180
C                    VALUES OF THE TABLE (DESTROYED).                   ACFI 190
C           Y      - THE RESULTING INTERPOLATED FUNCTION VALUE.         ACFI 200
C           NDIM   - AN INPUT VALUE WHICH SPECIFIES THE NUMBER OF       ACFI 210
C                    POINTS IN TABLE (ARG,VAL).                         ACFI 220
C           EPS    - AN INPUT CONSTANT WHICH IS USED AS UPPER BOUND     ACFI 230
C                    FOR THE ABSOLUTE ERROR.                            ACFI 240
C           IER    - A RESULTING ERROR PARAMETER.                       ACFI 250
C                                                                       ACFI 260
C        REMARKS                                                        ACFI 270
C           (1) TABLE (ARG,VAL) SHOULD REPRESENT A SINGLE-VALUED        ACFI 280
C               FUNCTION AND SHOULD BE STORED IN SUCH A WAY, THAT THE   ACFI 290
C               DISTANCES ABS(ARG(I)-X) INCREASE WITH INCREASING        ACFI 300
C               SUBSCRIPT I. TO GENERATE THIS ORDER IN TABLE (ARG,VAL), ACFI 310
C               SUBROUTINES ATSG, ATSM OR ATSE COULD BE USED IN A       ACFI 320
C               PREVIOUS STAGE.                                         ACFI 330
C           (2) NO ACTION BESIDES ERROR MESSAGE IN CASE NDIM LESS       ACFI 340
C               THAN 1.                                                 ACFI 350
C           (3) INTERPOLATION IS TERMINATED EITHER IF THE DIFFERENCE    ACFI 360
C               BETWEEN TWO SUCCESSIVE INTERPOLATED VALUES IS           ACFI 370
C               ABSOLUTELY LESS THAN TOLERANCE EPS, OR IF THE ABSOLUTE  ACFI 380
C               VALUE OF THIS DIFFERENCE STOPS DIMINISHING, OR AFTER    ACFI 390
C               (NDIM-1) STEPS (THE NUMBER OF POSSIBLE STEPS IS         ACFI 400
C               DIMINISHED IF AT ANY STAGE INFINITY ELEMENT APPEARS IN  ACFI 410
C               THE DOWNWARD DIAGONAL OF INVERTED-DIFFERENCES-SCHEME    ACFI 420
C               AND IF IT IS IMPOSSIBLE TO ELIMINATE THIS INFINITY      ACFI 430
C               ELEMENT BY INTERCHANGING OF TABLE POINTS).              ACFI 440
C               FURTHER IT IS TERMINATED IF THE PROCEDURE DISCOVERS TWO ACFI 450
C               ARGUMENT VALUES IN VECTOR ARG WHICH ARE IDENTICAL.      ACFI 460
C               DEPENDENT ON THESE FOUR CASES, ERROR PARAMETER IER IS   ACFI 470
C               CODED IN THE FOLLOWING FORM                             ACFI 480
C                IER=0 - IT WAS POSSIBLE TO REACH THE REQUIRED          ACFI 490
C                        ACCURACY (NO ERROR).                           ACFI 500
C                IER=1 - IT WAS IMPOSSIBLE TO REACH THE REQUIRED        ACFI 510
C                        ACCURACY BECAUSE OF ROUNDING ERRORS.           ACFI 520
C                IER=2 - IT WAS IMPOSSIBLE TO CHECK ACCURACY BECAUSE    ACFI 530
C                        NDIM IS LESS THAN 2, OR THE REQUIRED ACCURACY  ACFI 540
C                        COULD NOT BE REACHED BY MEANS OF THE GIVEN     ACFI 550
C                        TABLE. NDIM SHOULD BE INCREASED.               ACFI 560
C                IER=3 - THE PROCEDURE DISCOVERED TWO ARGUMENT VALUES   ACFI 570
C                        IN VECTOR ARG WHICH ARE IDENTICAL.             ACFI 580
C                                                                       ACFI 590
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  ACFI 600
C           NONE                                                        ACFI 610
C                                                                       ACFI 620
C        METHOD                                                         ACFI 630
C           INTERPOLATION IS DONE BY CONTINUED FRACTIONS AND INVERTED-  ACFI 640
C           DIFFERENCES-SCHEME. ON RETURN Y CONTAINS AN INTERPOLATED    ACFI 650
C           FUNCTION VALUE AT POINT X, WHICH IS IN THE SENSE OF REMARK  ACFI 660
C           (3) OPTIMAL WITH RESPECT TO GIVEN TABLE. FOR REFERENCE, SEE ACFI 670
C           F.B.HILDEBRAND, INTRODUCTION TO NUMERICAL ANALYSIS,         ACFI 680
C           MCGRAW-HILL, NEW YORK/TORONTO/LONDON, 1956, PP.395-406.     ACFI 690
C                                                                       ACFI 700
C     ..................................................................ACFI 710
C                                                                       ACFI 720
      SUBROUTINE ACFI(X,ARG,VAL,Y,NDIM,EPS,IER)                         ACFI 730
C                                                                       ACFI 740
C                                                                       ACFI 750
      DIMENSION ARG(1),VAL(1)                                           ACFI 760
      IER=2                                                             ACFI 770
      IF(NDIM)20,20,1                                                   ACFI 780
    1 Y=VAL(1)                                                          ACFI 790
      DELT2=0.                                                          ACFI 800
      IF(NDIM-1)20,20,2                                                 ACFI 810
C                                                                       ACFI 820
C     PREPARATIONS FOR INTERPOLATION LOOP                               ACFI 830
    2 P2=1.                                                             ACFI 840
      P3=Y                                                              ACFI 850
      Q2=0.                                                             ACFI 860
      Q3=1.                                                             ACFI 870
C                                                                       ACFI 880
C                                                                       ACFI 890
C     START INTERPOLATION LOOP                                          ACFI 900
      DO 16 I=2,NDIM                                                    ACFI 910
      II=0                                                              ACFI 920
      P1=P2                                                             ACFI 930
      P2=P3                                                             ACFI 940
      Q1=Q2                                                             ACFI 950
      Q2=Q3                                                             ACFI 960
      Z=Y                                                               ACFI 970
      DELT1=DELT2                                                       ACFI 980
      JEND=I-1                                                          ACFI 990
C                                                                       ACFI1000
C     COMPUTATION OF INVERTED DIFFERENCES                               ACFI1010
    3 AUX=VAL(I)                                                        ACFI1020
      DO 10 J=1,JEND                                                    ACFI1030
      H=VAL(I)-VAL(J)                                                   ACFI1040
      IF(ABS(H)-1.E-6*ABS(VAL(I)))4,4,9                                 ACFI1050
    4 IF(ARG(I)-ARG(J))5,17,5                                           ACFI1060
    5 IF(J-JEND)8,6,6                                                   ACFI1070
C                                                                       ACFI1080
C     INTERCHANGE ROW I WITH ROW I+II                                   ACFI1090
    6 II=II+1                                                           ACFI1100
      III=I+II                                                          ACFI1110
      IF(III-NDIM)7,7,19                                                ACFI1120
    7 VAL(I)=VAL(III)                                                   ACFI1130
      VAL(III)=AUX                                                      ACFI1140
      AUX=ARG(I)                                                        ACFI1150
      ARG(I)=ARG(III)                                                   ACFI1160
      ARG(III)=AUX                                                      ACFI1170
      GOTO 3                                                            ACFI1180
C                                                                       ACFI1190
C     COMPUTATION OF VAL(I) IN CASE VAL(I)=VAL(J) AND J LESS THAN I-1   ACFI1200
    8 VAL(I)=1.E75                                                      ACFI1210
      GOTO 10                                                           ACFI1220
C                                                                       ACFI1230
C     COMPUTATION OF VAL(I) IN CASE VAL(I) NOT EQUAL TO VAL(J)          ACFI1240
    9 VAL(I)=(ARG(I)-ARG(J))/H                                          ACFI1250
   10 CONTINUE                                                          ACFI1260
C     INVERTED DIFFERENCES ARE COMPUTED                                 ACFI1270
C                                                                       ACFI1280
C     COMPUTATION OF NEW Y                                              ACFI1290
      P3=VAL(I)*P2+(X-ARG(I-1))*P1                                      ACFI1300
      Q3=VAL(I)*Q2+(X-ARG(I-1))*Q1                                      ACFI1310
      IF(Q3)11,12,11                                                    ACFI1320
   11 Y=P3/Q3                                                           ACFI1330
      GOTO 13                                                           ACFI1340
   12 Y=1.E75                                                           ACFI1350
   13 DELT2=ABS(Z-Y)                                                    ACFI1360
      IF(DELT2-EPS)19,19,14                                             ACFI1370
   14 IF(I-8)16,15,15                                                   ACFI1380
   15 IF(DELT2-DELT1)16,18,18                                           ACFI1390
   16 CONTINUE                                                          ACFI1400
C     END OF INTERPOLATION LOOP                                         ACFI1410
C                                                                       ACFI1420
C                                                                       ACFI1430
      RETURN                                                            ACFI1440
C                                                                       ACFI1450
C     THERE ARE TWO IDENTICAL ARGUMENT VALUES IN VECTOR ARG             ACFI1460
   17 IER=3                                                             ACFI1470
      RETURN                                                            ACFI1480
C                                                                       ACFI1490
C     TEST VALUE DELT2 STARTS OSCILLATING                               ACFI1500
   18 Y=Z                                                               ACFI1510
      IER=1                                                             ACFI1520
      RETURN                                                            ACFI1530
C                                                                       ACFI1540
C     THERE IS SATISFACTORY ACCURACY WITHIN NDIM-1 STEPS                ACFI1550
   19 IER=0                                                             ACFI1560
   20 RETURN                                                            ACFI1570
      END                                                               ACFI1580
