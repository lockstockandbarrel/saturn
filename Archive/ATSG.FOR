C                                                                       ATSG  10
C     ..................................................................ATSG  20
C                                                                       ATSG  30
C        SUBROUTINE ATSG                                                ATSG  40
C                                                                       ATSG  50
C        PURPOSE                                                        ATSG  60
C           NDIM POINTS OF A GIVEN GENERAL TABLE ARE SELECTED AND       ATSG  70
C           ORDERED SUCH THAT ABS(ARG(I)-X).GE.ABS(ARG(J)-X) IF I.GT.J. ATSG  80
C                                                                       ATSG  90
C        USAGE                                                          ATSG 100
C           CALL ATSG (X,Z,F,WORK,IROW,ICOL,ARG,VAL,NDIM)               ATSG 110
C                                                                       ATSG 120
C        DESCRIPTION OF PARAMETERS                                      ATSG 130
C           X      - THE SEARCH ARGUMENT.                               ATSG 140
C           Z      - THE VECTOR OF ARGUMENT VALUES (DIMENSION IROW).    ATSG 150
C           F      - IN CASE ICOL=1, F IS THE VECTOR OF FUNCTION VALUES ATSG 160
C                    (DIMENSION IROW).                                  ATSG 170
C                    IN CASE ICOL=2, F IS AN IROW BY 2 MATRIX. THE FIRSTATSG 180
C                    COLUMN SPECIFIES THE VECTOR OF FUNCTION VALUES AND ATSG 190
C                    THE SECOND THE VECTOR OF DERIVATIVES.              ATSG 200
C           WORK   - A WORKING STORAGE (DIMENSION IROW).                ATSG 210
C           IROW   - THE DIMENSION OF VECTORS Z AND WORK AND OF EACH    ATSG 220
C                    COLUMN IN MATRIX F.                                ATSG 230
C           ICOL   - THE NUMBER OF COLUMNS IN F (I.E. 1 OR 2).          ATSG 240
C           ARG    - THE RESULTING VECTOR OF SELECTED AND ORDERED       ATSG 250
C                    ARGUMENT VALUES (DIMENSION NDIM).                  ATSG 260
C           VAL    - THE RESULTING VECTOR OF SELECTED FUNCTION VALUES   ATSG 270
C                    (DIMENSION NDIM) IN CASE ICOL=1. IN CASE ICOL=2,   ATSG 280
C                    VAL IS THE VECTOR OF FUNCTION AND DERIVATIVE VALUESATSG 290
C                    (DIMENSION 2*NDIM) WHICH ARE STORED IN PAIRS (I.E. ATSG 300
C                    EACH FUNCTION VALUE IS FOLLOWED BY ITS DERIVATIVE  ATSG 310
C                    VALUE).                                            ATSG 320
C           NDIM   - THE NUMBER OF POINTS WHICH MUST BE SELECTED OUT OF ATSG 330
C                    THE GIVEN TABLE (Z,F).                             ATSG 340
C                                                                       ATSG 350
C        REMARKS                                                        ATSG 360
C           NO ACTION IN CASE IROW LESS THAN 1.                         ATSG 370
C           IF INPUT VALUE NDIM IS GREATER THAN IROW, THE PROGRAM       ATSG 380
C           SELECTS ONLY A MAXIMUM TABLE OF IROW POINTS.  THEREFORE THE ATSG 390
C           USER OUGHT TO CHECK CORRESPONDENCE BETWEEN TABLE (ARG,VAL)  ATSG 400
C           AND ITS DIMENSION BY COMPARISON OF NDIM AND IROW, IN ORDER  ATSG 410
C           TO GET CORRECT RESULTS IN FURTHER WORK WITH TABLE (ARG,VAL).ATSG 420
C           THIS TEST MAY BE DONE BEFORE OR AFTER CALLING               ATSG 430
C           SUBROUTINE ATSG.                                            ATSG 440
C           SUBROUTINE ATSG ESPECIALLY CAN BE USED FOR GENERATING THE   ATSG 450
C           TABLE (ARG,VAL) NEEDED IN SUBROUTINES ALI, AHI, AND ACFI.   ATSG 460
C                                                                       ATSG 470
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  ATSG 480
C           NONE                                                        ATSG 490
C                                                                       ATSG 500
C        METHOD                                                         ATSG 510
C           SELECTION IS DONE BY GENERATING THE VECTOR WORK WITH        ATSG 520
C           COMPONENTS WORK(I)=ABS(Z(I)-X) AND AT EACH OF THE NDIM STEPSATSG 530
C           (OR IROW STEPS IF NDIM IS GREATER THAN IROW)                ATSG 540
C           SEARCHING FOR THE SUBSCRIPT OF THE SMALLEST COMPONENT, WHICHATSG 550
C           IS AFTERWARDS REPLACED BY A NUMBER GREATER THAN             ATSG 560
C           MAX(WORK(I)).                                               ATSG 570
C                                                                       ATSG 580
C     ..................................................................ATSG 590
C                                                                       ATSG 600
      SUBROUTINE ATSG(X,Z,F,WORK,IROW,ICOL,ARG,VAL,NDIM)                ATSG 610
C                                                                       ATSG 620
C                                                                       ATSG 630
      DIMENSION Z(1),F(1),WORK(1),ARG(1),VAL(1)                         ATSG 640
      IF(IROW)11,11,1                                                   ATSG 650
    1 N=NDIM                                                            ATSG 660
C     IF N IS GREATER THAN IROW, N IS SET EQUAL TO IROW.                ATSG 670
      IF(N-IROW)3,3,2                                                   ATSG 680
    2 N=IROW                                                            ATSG 690
C                                                                       ATSG 700
C     GENERATION OF VECTOR WORK AND COMPUTATION OF ITS GREATEST ELEMENT.ATSG 710
    3 B=0.                                                              ATSG 720
      DO 5 I=1,IROW                                                     ATSG 730
      DELTA=ABS(Z(I)-X)                                                 ATSG 740
      IF(DELTA-B)5,5,4                                                  ATSG 750
    4 B=DELTA                                                           ATSG 760
    5 WORK(I)=DELTA                                                     ATSG 770
C                                                                       ATSG 780
C     GENERATION OF TABLE (ARG,VAL)                                     ATSG 790
      B=B+1.                                                            ATSG 800
      DO 10 J=1,N                                                       ATSG 810
      DELTA=B                                                           ATSG 820
      DO 7 I=1,IROW                                                     ATSG 830
      IF(WORK(I)-DELTA)6,7,7                                            ATSG 840
    6 II=I                                                              ATSG 850
      DELTA=WORK(I)                                                     ATSG 860
    7 CONTINUE                                                          ATSG 870
      ARG(J)=Z(II)                                                      ATSG 880
      IF(ICOL-1)8,9,8                                                   ATSG 890
    8 VAL(2*J-1)=F(II)                                                  ATSG 900
      III=II+IROW                                                       ATSG 910
      VAL(2*J)=F(III)                                                   ATSG 920
      GOTO 10                                                           ATSG 930
    9 VAL(J)=F(II)                                                      ATSG 940
   10 WORK(II)=B                                                        ATSG 950
   11 RETURN                                                            ATSG 960
      END                                                               ATSG 970
