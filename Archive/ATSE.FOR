C                                                                       ATSE  10
C     ..................................................................ATSE  20
C                                                                       ATSE  30
C        SUBROUTINE ATSE                                                ATSE  40
C                                                                       ATSE  50
C        PURPOSE                                                        ATSE  60
C           NDIM POINTS OF A GIVEN TABLE WITH EQUIDISTANT ARGUMENTS ARE ATSE  70
C           SELECTED AND ORDERED SUCH THAT                              ATSE  80
C           ABS(ARG(I)-X).GE.ABS(ARG(J)-X) IF I.GT.J.                   ATSE  90
C                                                                       ATSE 100
C        USAGE                                                          ATSE 110
C           CALL ATSE (X,ZS,DZ,F,IROW,ICOL,ARG,VAL,NDIM)                ATSE 120
C                                                                       ATSE 130
C        DESCRIPTION OF PARAMETERS                                      ATSE 140
C           X      - THE SEARCH ARGUMENT.                               ATSE 150
C           ZS     - THE STARTING VALUE OF ARGUMENTS.                   ATSE 160
C           DZ     - THE INCREMENT OF ARGUMENT VALUES.                  ATSE 170
C           F      - IN CASE ICOL=1, F IS THE VECTOR OF FUNCTION VALUES ATSE 180
C                    (DIMENSION IROW).                                  ATSE 190
C                    IN CASE ICOL=2, F IS AN IROW BY 2 MATRIX. THE FIRSTATSE 200
C                    COLUMN SPECIFIES THE VECTOR OF FUNCTION VALUES AND ATSE 210
C                    THE SECOND THE VECTOR OF DERIVATIVES.              ATSE 220
C           IROW   - THE DIMENSION OF EACH COLUMN IN MATRIX F.          ATSE 230
C           ICOL   - THE NUMBER OF COLUMNS IN F (I.E. 1 OR 2).          ATSE 240
C           ARG    - THE RESULTING VECTOR OF SELECTED AND ORDERED       ATSE 250
C                    ARGUMENT VALUES (DIMENSION NDIM).                  ATSE 260
C           VAL    - THE RESULTING VECTOR OF SELECTED FUNCTION VALUES   ATSE 270
C                    (DIMENSION NDIM) IN CASE ICOL=1. IN CASE ICOL=2,   ATSE 280
C                    VAL IS THE VECTOR OF FUNCTION AND DERIVATIVE VALUESATSE 290
C                    (DIMENSION 2*NDIM) WHICH ARE STORED IN PAIRS (I.E. ATSE 300
C                    EACH FUNCTION VALUE IS FOLLOWED BY ITS DERIVATIVE  ATSE 310
C                    VALUE).                                            ATSE 320
C           NDIM   - THE NUMBER OF POINTS WHICH MUST BE SELECTED OUT OF ATSE 330
C                    THE GIVEN TABLE.                                   ATSE 340
C                                                                       ATSE 350
C        REMARKS                                                        ATSE 360
C           NO ACTION IN CASE IROW LESS THAN 1.                         ATSE 370
C           IF INPUT VALUE NDIM IS GREATER THAN IROW, THE PROGRAM       ATSE 380
C           SELECTS ONLY A MAXIMUM TABLE OF IROW POINTS.  THEREFORE THE ATSE 390
C           USER OUGHT TO CHECK CORRESPONDENCE BETWEEN TABLE (ARG,VAL)  ATSE 400
C           AND ITS DIMENSION BY COMPARISON OF NDIM AND IROW, IN ORDER  ATSE 410
C           TO GET CORRECT RESULTS IN FURTHER WORK WITH TABLE (ARG,VAL).ATSE 420
C           THIS TEST MAY BE DONE BEFORE OR AFTER CALLING               ATSE 430
C           SUBROUTINE ATSE.                                            ATSE 440
C           SUBROUTINE ATSE ESPECIALLY CAN BE USED FOR GENERATING THE   ATSE 450
C           TABLE (ARG,VAL) NEEDED IN SUBROUTINES ALI, AHI, AND ACFI.   ATSE 460
C                                                                       ATSE 470
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  ATSE 480
C           NONE                                                        ATSE 490
C                                                                       ATSE 500
C        METHOD                                                         ATSE 510
C           SELECTION IS DONE BY COMPUTING THE SUBSCRIPT J OF THAT      ATSE 520
C           ARGUMENT, WHICH IS NEXT TO X.                               ATSE 530
C           AFTERWARDS NEIGHBOURING ARGUMENT VALUES ARE TESTED AND      ATSE 540
C           SELECTED IN THE ABOVE SENSE.                                ATSE 550
C                                                                       ATSE 560
C     ..................................................................ATSE 570
C                                                                       ATSE 580
      SUBROUTINE ATSE(X,ZS,DZ,F,IROW,ICOL,ARG,VAL,NDIM)                 ATSE 590
C                                                                       ATSE 600
C                                                                       ATSE 610
      DIMENSION F(1),ARG(1),VAL(1)                                      ATSE 620
      IF(IROW-1)19,17,1                                                 ATSE 630
C                                                                       ATSE 640
C     CASE DZ=0 IS CHECKED OUT                                          ATSE 650
    1 IF(DZ)2,17,2                                                      ATSE 660
    2 N=NDIM                                                            ATSE 670
C                                                                       ATSE 680
C     IF N IS GREATER THAN IROW, N IS SET EQUAL TO IROW.                ATSE 690
      IF(N-IROW)4,4,3                                                   ATSE 700
    3 N=IROW                                                            ATSE 710
C                                                                       ATSE 720
C     COMPUTATION OF STARTING SUBSCRIPT J.                              ATSE 730
    4 J=(X-ZS)/DZ+1.5                                                   ATSE 740
      IF(J)5,5,6                                                        ATSE 750
    5 J=1                                                               ATSE 760
    6 IF(J-IROW)8,8,7                                                   ATSE 770
    7 J=IROW                                                            ATSE 780
C                                                                       ATSE 790
C     GENERATION OF TABLE ARG,VAL IN CASE DZ.NE.0.                      ATSE 800
    8 II=J                                                              ATSE 810
      JL=0                                                              ATSE 820
      JR=0                                                              ATSE 830
      DO 16 I=1,N                                                       ATSE 840
      ARG(I)=ZS+FLOAT(II-1)*DZ                                          ATSE 850
      IF(ICOL-2)9,10,10                                                 ATSE 860
    9 VAL(I)=F(II)                                                      ATSE 870
      GOTO 11                                                           ATSE 880
   10 VAL(2*I-1)=F(II)                                                  ATSE 890
      III=II+IROW                                                       ATSE 900
      VAL(2*I)=F(III)                                                   ATSE 910
   11 IF(J+JR-IROW)12,15,12                                             ATSE 920
   12 IF(J-JL-1)13,14,13                                                ATSE 930
   13 IF((ARG(I)-X)*DZ)14,15,15                                         ATSE 940
   14 JR=JR+1                                                           ATSE 950
      II=J+JR                                                           ATSE 960
      GOTO 16                                                           ATSE 970
   15 JL=JL+1                                                           ATSE 980
      II=J-JL                                                           ATSE 990
   16 CONTINUE                                                          ATSE1000
      RETURN                                                            ATSE1010
C                                                                       ATSE1020
C     CASE DZ=0                                                         ATSE1030
   17 ARG(1)=ZS                                                         ATSE1040
      VAL(1)=F(1)                                                       ATSE1050
      IF(ICOL-2)19,19,18                                                ATSE1060
   18 VAL(2)=F(2)                                                       ATSE1070
   19 RETURN                                                            ATSE1080
      END                                                               ATSE1090
