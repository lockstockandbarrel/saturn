C                                                                       DTSM  10
C     ..................................................................DTSM  20
C                                                                       DTSM  30
C        SUBROUTINE DATSM                                               DTSM  40
C                                                                       DTSM  50
C        PURPOSE                                                        DTSM  60
C           NDIM POINTS OF A GIVEN TABLE WITH MONOTONIC ARGUMENTS ARE   DTSM  70
C           SELECTED AND ORDERED SUCH THAT                              DTSM  80
C           ABS(ARG(I)-X).GE.ABS(ARG(J)-X) IF I.GT.J.                   DTSM  90
C                                                                       DTSM 100
C        USAGE                                                          DTSM 110
C           CALL DATSM (X,Z,F,IROW,ICOL,ARG,VAL,NDIM)                   DTSM 120
C                                                                       DTSM 130
C        DESCRIPTION OF PARAMETERS                                      DTSM 140
C           X      - DOUBLE PRECISION SEARCH ARGUMENT.                  DTSM 150
C           Z      - DOUBLE PRECISION VECTOR OF ARGUMENT VALUES (DIMEN- DTSM 160
C                    SION IROW). THE ARGUMENT VALUES MUST BE STORED IN  DTSM 170
C                    INCREASING OR DECREASING SEQUENCE.                 DTSM 180
C           F      - IN CASE ICOL=1, F IS THE DOUBLE PRECISION VECTOR   DTSM 190
C                    OF FUNCTION VALUES (DIMENSION IROW).               DTSM 200
C                    IN CASE ICOL=2, F IS A DOUBLE PRECISION IROW BY 2  DTSM 210
C                    MATRIX. THE FIRST COLUMN SPECIFIES VECTOR OF FUNC- DTSM 220
C                   TION VALUES AND THE SECOND VECTOR OF DERIVATIVES.   DTSM 230
C           IROW   - THE DIMENSION OF VECTOR Z AND OF EACH COLUMN       DTSM 240
C                    IN MATRIX F.                                       DTSM 250
C           ICOL   - THE NUMBER OF COLUMNS IN F (I.E. 1 OR 2).          DTSM 260
C           ARG    - RESULTING DOUBLE PRECISION VECTOR OF SELECTED AND  DTSM 270
C                    ORDERED ARGUMENT VALUES (DIMENSION NDIM).          DTSM 280
C           VAL    - RESULTING DOUBLE PRECISION VECTOR OF SELECTED      DTSM 290
C                    FUNCTION VALUES (DIMENSION NDIM) IN CASE ICOL=1.   DTSM 300
C                    IN CASE ICOL=2, VAL IS THE DOUBLE PRECISION VECTOR DTSM 310
C                    OF FUNCTION AND DERIVATIVE VALUES (DIMENSION       DTSM 320
C                    2*NDIM) WHICH ARE STORED IN PAIRS (I.E. EACH FUNC- DTSM 330
C                    TION VALUE IS FOLLOWED BY ITS DERIVATIVE VALUE).   DTSM 340
C           NDIM   - THE NUMBER OF POINTS WHICH MUST BE SELECTED OUT OF DTSM 350
C                    THE GIVEN TABLE (Z,F).                             DTSM 360
C                                                                       DTSM 370
C        REMARKS                                                        DTSM 380
C           NO ACTION IN CASE IROW LESS THAN 1.                         DTSM 390
C           IF INPUT VALUE NDIM IS GREATER THAN IROW, THE PROGRAM       DTSM 400
C           SELECTS ONLY A MAXIMUM TABLE OF IROW POINTS.  THEREFORE THE DTSM 410
C           USER OUGHT TO CHECK CORRESPONDENCE BETWEEN TABLE (ARG,VAL)  DTSM 420
C           AND ITS DIMENSION BY COMPARISON OF NDIM AND IROW, IN ORDER  DTSM 430
C           TO GET CORRECT RESULTS IN FURTHER WORK WITH TABLE (ARG,VAL).DTSM 440
C           THIS TEST MAY BE DONE BEFORE OR AFTER CALLING               DTSM 450
C           SUBROUTINE DATSM.                                           DTSM 460
C           SUBROUTINE DATSM ESPECIALLY CAN BE USED FOR GENERATING THE  DTSM 470
C           TABLE (ARG,VAL) NEEDED IN SUBROUTINES DALI, DAHI, AND DACFI.DTSM 480
C                                                                       DTSM 490
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DTSM 500
C           NONE                                                        DTSM 510
C                                                                       DTSM 520
C        METHOD                                                         DTSM 530
C           SELECTION IS DONE BY SEARCHING THE SUBSCRIPT J OF THAT      DTSM 540
C           ARGUMENT, WHICH IS NEXT TO X (BINARY SEARCH).               DTSM 550
C           AFTERWARDS NEIGHBOURING ARGUMENT VALUES ARE TESTED AND      DTSM 560
C           SELECTED IN THE ABOVE SENSE.                                DTSM 570
C                                                                       DTSM 580
C     ..................................................................DTSM 590
C                                                                       DTSM 600
      SUBROUTINE DATSM(X,Z,F,IROW,ICOL,ARG,VAL,NDIM)                    DTSM 610
C                                                                       DTSM 620
C                                                                       DTSM 630
      DIMENSION Z(1),F(1),ARG(1),VAL(1)                                 DTSM 640
      DOUBLE PRECISION X,Z,F,ARG,VAL                                    DTSM 650
C                                                                       DTSM 660
C     CASE IROW=1 IS CHECKED OUT                                        DTSM 670
      IF(IROW-1)23,21,1                                                 DTSM 680
    1 N=NDIM                                                            DTSM 690
C                                                                       DTSM 700
C     IF N IS GREATER THAN IROW, N IS SET EQUAL TO IROW.                DTSM 710
      IF(N-IROW)3,3,2                                                   DTSM 720
    2 N=IROW                                                            DTSM 730
C                                                                       DTSM 740
C     CASE IROW.GE.2                                                    DTSM 750
C     SEARCHING FOR SUBSCRIPT J SUCH THAT Z(J) IS NEXT TO X.            DTSM 760
    3 IF(Z(IROW)-Z(1))5,4,4                                             DTSM 770
    4 J=IROW                                                            DTSM 780
      I=1                                                               DTSM 790
      GOTO 6                                                            DTSM 800
    5 I=IROW                                                            DTSM 810
      J=1                                                               DTSM 820
    6 K=(J+I)/2                                                         DTSM 830
      IF(X-Z(K))7,7,8                                                   DTSM 840
    7 J=K                                                               DTSM 850
      GOTO 9                                                            DTSM 860
    8 I=K                                                               DTSM 870
    9 IF(IABS(J-I)-1)10,10,6                                            DTSM 880
   10 IF(DABS(Z(J)-X)-DABS(Z(I)-X))12,12,11                             DTSM 890
   11 J=I                                                               DTSM 900
C                                                                       DTSM 910
C     TABLE SELECTION                                                   DTSM 920
   12 K=J                                                               DTSM 930
      JL=0                                                              DTSM 940
      JR=0                                                              DTSM 950
      DO 20 I=1,N                                                       DTSM 960
      ARG(I)=Z(K)                                                       DTSM 970
      IF(ICOL-1)14,14,13                                                DTSM 980
   13 VAL(2*I-1)=F(K)                                                   DTSM 990
      KK=K+IROW                                                         DTSM1000
      VAL(2*I)=F(KK)                                                    DTSM1010
      GOTO 15                                                           DTSM1020
   14 VAL(I)=F(K)                                                       DTSM1030
   15 JJR=J+JR                                                          DTSM1040
      IF(JJR-IROW)16,18,18                                              DTSM1050
   16 JJL=J-JL                                                          DTSM1060
      IF(JJL-1)19,19,17                                                 DTSM1070
   17 IF(DABS(Z(JJR+1)-X)-DABS(Z(JJL-1)-X))19,19,18                     DTSM1080
   18 JL=JL+1                                                           DTSM1090
      K=J-JL                                                            DTSM1100
      GOTO 20                                                           DTSM1110
   19 JR=JR+1                                                           DTSM1120
      K=J+JR                                                            DTSM1130
   20 CONTINUE                                                          DTSM1140
      RETURN                                                            DTSM1150
C                                                                       DTSM1160
C     CASE IROW=1                                                       DTSM1170
   21 ARG(1)=Z(1)                                                       DTSM1180
      VAL(1)=F(1)                                                       DTSM1190
      IF(ICOL-2)23,22,23                                                DTSM1200
   22 VAL(2)=F(2)                                                       DTSM1210
   23 RETURN                                                            DTSM1220
      END                                                               DTSM1230
