C                                                                       ATSM  10
C     ..................................................................ATSM  20
C                                                                       ATSM  30
C        SUBROUTINE ATSM                                                ATSM  40
C                                                                       ATSM  50
C        PURPOSE                                                        ATSM  60
C           NDIM POINTS OF A GIVEN TABLE WITH MONOTONIC ARGUMENTS ARE   ATSM  70
C           SELECTED AND ORDERED SUCH THAT                              ATSM  80
C           ABS(ARG(I)-X).GE.ABS(ARG(J)-X) IF I.GT.J.                   ATSM  90
C                                                                       ATSM 100
C        USAGE                                                          ATSM 110
C           CALL ATSM (X,Z,F,IROW,ICOL,ARG,VAL,NDIM)                    ATSM 120
C                                                                       ATSM 130
C        DESCRIPTION OF PARAMETERS                                      ATSM 140
C           X      - THE SEARCH ARGUMENT.                               ATSM 150
C           Z      - THE VECTOR OF ARGUMENT VALUES (DIMENSION IROW).    ATSM 160
C                    THE ARGUMENT VALUES MUST BE STORED IN INCREASING   ATSM 170
C                    OR DECREASING SEQUENCE.                            ATSM 180
C           F      - IN CASE ICOL=1, F IS THE VECTOR OF FUNCTION VALUES ATSM 190
C                    (DIMENSION IROW).                                  ATSM 200
C                    IN CASE ICOL=2, F IS AN IROW BY 2 MATRIX. THE FIRSTATSM 210
C                    COLUMN SPECIFIES THE VECTOR OF FUNCTION VALUES AND ATSM 220
C                    THE SECOND THE VECTOR OF DERIVATIVES.              ATSM 230
C           IROW   - THE DIMENSION OF VECTOR Z AND OF EACH COLUMN       ATSM 240
C                    IN MATRIX F.                                       ATSM 250
C           ICOL   - THE NUMBER OF COLUMNS IN F (I.E. 1 OR 2).          ATSM 260
C           ARG    - THE RESULTING VECTOR OF SELECTED AND ORDERED       ATSM 270
C                    ARGUMENT VALUES (DIMENSION NDIM).                  ATSM 280
C           VAL    - THE RESULTING VECTOR OF SELECTED FUNCTION VALUES   ATSM 290
C                    (DIMENSION NDIM) IN CASE ICOL=1. IN CASE ICOL=2,   ATSM 300
C                    VAL IS THE VECTOR OF FUNCTION AND DERIVATIVE VALUESATSM 310
C                    (DIMENSION 2*NDIM) WHICH ARE STORED IN PAIRS (I.E. ATSM 320
C                    EACH FUNCTION VALUE IS FOLLOWED BY ITS DERIVATIVE  ATSM 330
C                    VALUE).                                            ATSM 340
C           NDIM   - THE NUMBER OF POINTS WHICH MUST BE SELECTED OUT OF ATSM 350
C                    THE GIVEN TABLE (Z,F).                             ATSM 360
C                                                                       ATSM 370
C        REMARKS                                                        ATSM 380
C           NO ACTION IN CASE IROW LESS THAN 1.                         ATSM 390
C           IF INPUT VALUE NDIM IS GREATER THAN IROW, THE PROGRAM       ATSM 400
C           SELECTS ONLY A MAXIMUM TABLE OF IROW POINTS.  THEREFORE THE ATSM 410
C           USER OUGHT TO CHECK CORRESPONDENCE BETWEEN TABLE (ARG,VAL)  ATSM 420
C           AND ITS DIMENSION BY COMPARISON OF NDIM AND IROW, IN ORDER  ATSM 430
C           TO GET CORRECT RESULTS IN FURTHER WORK WITH TABLE (ARG,VAL).ATSM 440
C           THIS TEST MAY BE DONE BEFORE OR AFTER CALLING               ATSM 450
C           SUBROUTINE ATSM.                                            ATSM 460
C           SUBROUTINE ATSM ESPECIALLY CAN BE USED FOR GENERATING THE   ATSM 470
C           TABLE (ARG,VAL) NEEDED IN SUBROUTINES ALI, AHI, AND ACFI.   ATSM 480
C                                                                       ATSM 490
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  ATSM 500
C           NONE                                                        ATSM 510
C                                                                       ATSM 520
C        METHOD                                                         ATSM 530
C           SELECTION IS DONE BY SEARCHING THE SUBSCRIPT J OF THAT      ATSM 540
C           ARGUMENT, WHICH IS NEXT TO X (BINARY SEARCH).               ATSM 550
C           AFTERWARDS NEIGHBOURING ARGUMENT VALUES ARE TESTED AND      ATSM 560
C           SELECTED IN THE ABOVE SENSE.                                ATSM 570
C                                                                       ATSM 580
C     ..................................................................ATSM 590
C                                                                       ATSM 600
      SUBROUTINE ATSM(X,Z,F,IROW,ICOL,ARG,VAL,NDIM)                     ATSM 610
C                                                                       ATSM 620
C                                                                       ATSM 630
      DIMENSION Z(1),F(1),ARG(1),VAL(1)                                 ATSM 640
C                                                                       ATSM 650
C     CASE IROW=1 IS CHECKED OUT                                        ATSM 660
      IF(IROW-1)23,21,1                                                 ATSM 670
    1 N=NDIM                                                            ATSM 680
C                                                                       ATSM 690
C     IF N IS GREATER THAN IROW, N IS SET EQUAL TO IROW.                ATSM 700
      IF(N-IROW)3,3,2                                                   ATSM 710
    2 N=IROW                                                            ATSM 720
C                                                                       ATSM 730
C     CASE IROW.GE.2                                                    ATSM 740
C     SEARCHING FOR SUBSCRIPT J SUCH THAT Z(J) IS NEXT TO X.            ATSM 750
    3 IF(Z(IROW)-Z(1))5,4,4                                             ATSM 760
    4 J=IROW                                                            ATSM 770
      I=1                                                               ATSM 780
      GOTO 6                                                            ATSM 790
    5 I=IROW                                                            ATSM 800
      J=1                                                               ATSM 810
    6 K=(J+I)/2                                                         ATSM 820
      IF(X-Z(K))7,7,8                                                   ATSM 830
    7 J=K                                                               ATSM 840
      GOTO 9                                                            ATSM 850
    8 I=K                                                               ATSM 860
    9 IF(IABS(J-I)-1)10,10,6                                            ATSM 870
   10 IF(ABS(Z(J)-X)-ABS(Z(I)-X))12,12,11                               ATSM 880
   11 J=I                                                               ATSM 890
C                                                                       ATSM 900
C     TABLE SELECTION                                                   ATSM 910
   12 K=J                                                               ATSM 920
      JL=0                                                              ATSM 930
      JR=0                                                              ATSM 940
      DO 20 I=1,N                                                       ATSM 950
      ARG(I)=Z(K)                                                       ATSM 960
      IF(ICOL-1)14,14,13                                                ATSM 970
   13 VAL(2*I-1)=F(K)                                                   ATSM 980
      KK=K+IROW                                                         ATSM 990
      VAL(2*I)=F(KK)                                                    ATSM1000
      GOTO 15                                                           ATSM1010
   14 VAL(I)=F(K)                                                       ATSM1020
   15 JJR=J+JR                                                          ATSM1030
      IF(JJR-IROW)16,18,18                                              ATSM1040
   16 JJL=J-JL                                                          ATSM1050
      IF(JJL-1)19,19,17                                                 ATSM1060
   17 IF(ABS(Z(JJR+1)-X)-ABS(Z(JJL-1)-X))19,19,18                       ATSM1070
   18 JL=JL+1                                                           ATSM1080
      K=J-JL                                                            ATSM1090
      GOTO 20                                                           ATSM1100
   19 JR=JR+1                                                           ATSM1110
      K=J+JR                                                            ATSM1120
   20 CONTINUE                                                          ATSM1130
      RETURN                                                            ATSM1140
C                                                                       ATSM1150
C     CASE IROW=1                                                       ATSM1160
   21 ARG(1)=Z(1)                                                       ATSM1170
      VAL(1)=F(1)                                                       ATSM1180
      IF(ICOL-2)23,22,23                                                ATSM1190
   22 VAL(2)=F(2)                                                       ATSM1200
   23 RETURN                                                            ATSM1210
      END                                                               ATSM1220
