C                                                                       DTSE  10
C     ..................................................................DTSE  20
C                                                                       DTSE  30
C        SUBROUTINE DATSE                                               DTSE  40
C                                                                       DTSE  50
C        PURPOSE                                                        DTSE  60
C           NDIM POINTS OF A GIVEN TABLE WITH EQUIDISTANT ARGUMENTS ARE DTSE  70
C           SELECTED AND ORDERED SUCH THAT                              DTSE  80
C           ABS(ARG(I)-X).GE.ABS(ARG(J)-X) IF I.GT.J.                   DTSE  90
C                                                                       DTSE 100
C        USAGE                                                          DTSE 110
C           CALL DATSE (X,ZS,DZ,F,IROW,ICOL,ARG,VAL,NDIM)               DTSE 120
C                                                                       DTSE 130
C        DESCRIPTION OF PARAMETERS                                      DTSE 140
C           X      - DOUBLE PRECISION SEARCH ARGUMENT.                  DTSE 150
C           ZS     - DOUBLE PRECISION STARTING VALUE OF ARGUMENTS.      DTSE 160
C           DZ     - DOUBLE PRECISION INCREMENT OF ARGUMENT VALUES.     DTSE 170
C           F      - IN CASE ICOL=1, F IS THE DOUBLE PRECISION VECTOR   DTSE 180
C                    OF FUNCTION VALUES (DIMENSION IROW).               DTSE 190
C                    IN CASE ICOL=2, F IS A DOUBLE PRECISION IROW BY 2  DTSE 200
C                    MATRIX. THE FIRST COLUMN SPECIFIES VECTOR OF FUNC- DTSE 210
C                    TION VALUES AND THE SECOND VECTOR OF DERIVATIVES.  DTSE 220
C           IROW   - THE DIMENSION OF EACH COLUMN IN MATRIX F.          DTSE 230
C           ICOL   - THE NUMBER OF COLUMNS IN F (I.E. 1 OR 2).          DTSE 240
C           ARG    - RESULTING DOUBLE PRECISION VECTOR OF SELECTED AND  DTSE 250
C                    ORDERED ARGUMENT VALUES (DIMENSION NDIM).          DTSE 260
C           VAL    - RESULTING DOUBLE PRECISION VECTOR OF SELECTED      DTSE 270
C                    FUNCTION VALUES (DIMENSION NDIM) IN CASE ICOL=1.   DTSE 280
C                    IN CASE ICOL=2, VAL IS THE DOUBLE PRECISION VECTOR DTSE 290
C                    OF FUNCTION AND DERIVATIVE VALUES (DIMENSION       DTSE 300
C                    2*NDIM) WHICH ARE STORED IN PAIRS (I.E. EACH FUNC- DTSE 310
C                    TION VALUE IS FOLLOWED BY ITS DERIVATIVE VALUE).   DTSE 320
C           NDIM   - THE NUMBER OF POINTS WHICH MUST BE SELECTED OUT OF DTSE 330
C                    THE GIVEN TABLE.                                   DTSE 340
C                                                                       DTSE 350
C        REMARKS                                                        DTSE 360
C           NO ACTION IN CASE IROW LESS THAN 1.                         DTSE 370
C           IF INPUT VALUE NDIM IS GREATER THAN IROW, THE PROGRAM       DTSE 380
C           SELECTS ONLY A MAXIMUM TABLE OF IROW POINTS.  THEREFORE THE DTSE 390
C           USER OUGHT TO CHECK CORRESPONDENCE BETWEEN TABLE (ARG,VAL)  DTSE 400
C           AND ITS DIMENSION BY COMPARISON OF NDIM AND IROW, IN ORDER  DTSE 410
C           TO GET CORRECT RESULTS IN FURTHER WORK WITH TABLE (ARG,VAL).DTSE 420
C           THIS TEST MAY BE DONE BEFORE OR AFTER CALLING               DTSE 430
C           SUBROUTINE DATSE.                                           DTSE 440
C           SUBROUTINE DATSE ESPECIALLY CAN BE USED FOR GENERATING THE  DTSE 450
C           TABLE (ARG,VAL) NEEDED IN SUBROUTINES DALI, DAHI, AND DACFI.DTSE 460
C                                                                       DTSE 470
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DTSE 480
C           NONE                                                        DTSE 490
C                                                                       DTSE 500
C        METHOD                                                         DTSE 510
C           SELECTION IS DONE BY COMPUTING THE SUBSCRIPT J OF THAT      DTSE 520
C           ARGUMENT, WHICH IS NEXT TO X.                               DTSE 530
C           AFTERWARDS NEIGHBOURING ARGUMENT VALUES ARE TESTED AND      DTSE 540
C           SELECTED IN THE ABOVE SENSE.                                DTSE 550
C                                                                       DTSE 560
C     ..................................................................DTSE 570
C                                                                       DTSE 580
      SUBROUTINE DATSE(X,ZS,DZ,F,IROW,ICOL,ARG,VAL,NDIM)                DTSE 590
C                                                                       DTSE 600
C                                                                       DTSE 610
      DIMENSION F(1),ARG(1),VAL(1)                                      DTSE 620
      DOUBLE PRECISION X,ZS,DZ,F,ARG,VAL                                DTSE 630
      IF(IROW-1)19,17,1                                                 DTSE 640
C                                                                       DTSE 650
C     CASE DZ=0 IS CHECKED OUT                                          DTSE 660
    1 IF(DZ)2,17,2                                                      DTSE 670
    2 N=NDIM                                                            DTSE 680
C                                                                       DTSE 690
C     IF N IS GREATER THAN IROW, N IS SET EQUAL TO IROW.                DTSE 700
      IF(N-IROW)4,4,3                                                   DTSE 710
    3 N=IROW                                                            DTSE 720
C                                                                       DTSE 730
C     COMPUTATION OF STARTING SUBSCRIPT J.                              DTSE 740
    4 J=(X-ZS)/DZ+1.5D0                                                 DTSE 750
      IF(J)5,5,6                                                        DTSE 760
    5 J=1                                                               DTSE 770
    6 IF(J-IROW)8,8,7                                                   DTSE 780
    7 J=IROW                                                            DTSE 790
C                                                                       DTSE 800
C     GENERATION OF TABLE ARG,VAL IN CASE DZ.NE.0.                      DTSE 810
    8 II=J                                                              DTSE 820
      JL=0                                                              DTSE 830
      JR=0                                                              DTSE 840
      DO 16 I=1,N                                                       DTSE 850
      ARG(I)=ZS+DFLOAT(II-1)*DZ                                         DTSE 860
      IF(ICOL-2)9,10,10                                                 DTSE 870
    9 VAL(I)=F(II)                                                      DTSE 880
      GOTO 11                                                           DTSE 890
   10 VAL(2*I-1)=F(II)                                                  DTSE 900
      III=II+IROW                                                       DTSE 910
      VAL(2*I)=F(III)                                                   DTSE 920
   11 IF(J+JR-IROW)12,15,12                                             DTSE 930
   12 IF(J-JL-1)13,14,13                                                DTSE 940
   13 IF((ARG(I)-X)*DZ)14,15,15                                         DTSE 950
   14 JR=JR+1                                                           DTSE 960
      II=J+JR                                                           DTSE 970
      GOTO 16                                                           DTSE 980
   15 JL=JL+1                                                           DTSE 990
      II=J-JL                                                           DTSE1000
   16 CONTINUE                                                          DTSE1010
      RETURN                                                            DTSE1020
C                                                                       DTSE1030
C     CASE DZ=0                                                         DTSE1040
   17 ARG(1)=ZS                                                         DTSE1050
      VAL(1)=F(1)                                                       DTSE1060
      IF(ICOL-2)19,19,18                                                DTSE1070
   18 VAL(2)=F(2)                                                       DTSE1080
   19 RETURN                                                            DTSE1090
      END                                                               DTSE1100
