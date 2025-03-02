C                                                                       DTSG  10
C     ..................................................................DTSG  20
C                                                                       DTSG  30
C        SUBROUTINE DATSG                                               DTSG  40
C                                                                       DTSG  50
C        PURPOSE                                                        DTSG  60
C           NDIM POINTS OF A GIVEN GENERAL TABLE ARE SELECTED AND       DTSG  70
C           ORDERED SUCH THAT ABS(ARG(I)-X).GE.ABS(ARG(J)-X) IF I.GT.J. DTSG  80
C                                                                       DTSG  90
C        USAGE                                                          DTSG 100
C           CALL DATSG (X,Z,F,WORK,IROW,ICOL,ARG,VAL,NDIM)              DTSG 110
C                                                                       DTSG 120
C        DESCRIPTION OF PARAMETERS                                      DTSG 130
C           X      - DOUBLE PRECISION SEARCH ARGUMENT.                  DTSG 140
C           Z      - DOUBLE PRECISION VECTOR OD ARGUMENT VALUES         DTSG 150
C                    (DIMENSION IROW).                                  DTSG 160
C           F      - IN CASE ICOL=1, F IS THE DOUBLE PRECISION VECTOR   DTSG 170
C                    OF FUNCTION VALUES (DIMENSION IROW).               DTSG 180
C                    IN CASE ICOL=2, F IS A DOUBLE PRECISION IROW BY 2  DTSG 190
C                    MATRIX. THE FIRST COLUMN SPECIFIES VECTOR OF FUNC- DTSG 200
C                    TION VALUES AND THE SECOND VECTOR OF DERIVATIVES.  DTSG 210
C           WORK   - DOUBLE PRECISION WORKING STORAGE (DIMENSION IROW). DTSG 220
C           IROW   - THE DIMENSION OF VECTORS Z AND WORK AND OF EACH    DTSG 230
C                    COLUMN IN MATRIX F.                                DTSG 240
C           ICOL   - THE NUMBER OF COLUMNS IN F (I.E. 1 OR 2).          DTSG 250
C           ARG    - RESULTING DOUBLE PRECISION VECTOR OF SELECTED AND  DTSG 260
C                    ORDERED ARGUMENT VALUES (DIMENSION NDIM).          DTSG 270
C           VAL    - RESULTING DOUBLE PRECISION VECTOR OF SELECTED      DTSG 280
C                    FUNCTION VALUES (DIMENSION NDIM) IN CASE ICOL=1.   DTSG 290
C                    IN CASE ICOL=2, VAL IS THE DOUBLE PRECISION VECTOR DTSG 300
C                    OF FUNCTION AND DERIVATIVE VALUES (DIMENSION       DTSG 310
C                    2*NDIM) WHICH ARE STORED IN PAIRS (I.E. EACH FUNC- DTSG 320
C                    TION VALUE IS FOLLOWED BY ITS DERIVATIVE VALUE).   DTSG 330
C           NDIM   - THE NUMBER OF POINTS WHICH MUST BE SELECTED OUT OF DTSG 340
C                    THE GIVEN TABLE (Z,F).                             DTSG 350
C                                                                       DTSG 360
C        REMARKS                                                        DTSG 370
C           NO ACTION IN CASE IROW LESS THAN 1.                         DTSG 380
C           IF INPUT VALUE NDIM IS GREATER THAN IROW, THE PROGRAM       DTSG 390
C           SELECTS ONLY A MAXIMUM TABLE OF IROW POINTS.  THEREFORE THE DTSG 400
C           USER OUGHT TO CHECK CORRESPONDENCE BETWEEN TABLE (ARG,VAL)  DTSG 410
C           AND ITS DIMENSION BY COMPARISON OF NDIM AND IROW, IN ORDER  DTSG 420
C           TO GET CORRECT RESULTS IN FURTHER WORK WITH TABLE (ARG,VAL).DTSG 430
C           THIS TEST MAY BE DONE BEFORE OR AFTER CALLING               DTSG 440
C           SUBROUTINE DATSG.                                           DTSG 450
C           SUBROUTINE DATSG ESPECIALLY CAN BE USED FOR GENERATING THE  DTSG 460
C           TABLE (ARG,VAL) NEEDED IN SUBROUTINES DALI, DAHI, AND DACFI.DTSG 470
C                                                                       DTSG 480
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DTSG 490
C           NONE                                                        DTSG 500
C                                                                       DTSG 510
C        METHOD                                                         DTSG 520
C           SELECTION IS DONE BY GENERATING THE VECTOR WORK WITH        DTSG 530
C           COMPONENTS WORK(I)=ABS(Z(I)-X) AND AT EACH OF THE NDIM STEPSDTSG 540
C           (OR IROW STEPS IF NDIM IS GREATER THAN IROW)                DTSG 550
C           SEARCHING FOR THE SUBSCRIPT OF THE SMALLEST COMPONENT, WHICHDTSG 560
C           IS AFTERWARDS REPLACED BY A NUMBER GREATER THAN             DTSG 570
C           MAX(WORK(I)).                                               DTSG 580
C                                                                       DTSG 590
C     ..................................................................DTSG 600
C                                                                       DTSG 610
      SUBROUTINE DATSG(X,Z,F,WORK,IROW,ICOL,ARG,VAL,NDIM)               DTSG 620
C                                                                       DTSG 630
C                                                                       DTSG 640
      DIMENSION Z(1),F(1),WORK(1),ARG(1),VAL(1)                         DTSG 650
      DOUBLE PRECISION X,Z,F,WORK,ARG,VAL,B,DELTA                       DTSG 660
      IF(IROW)11,11,1                                                   DTSG 670
    1 N=NDIM                                                            DTSG 680
C     IF N IS GREATER THAN IROW, N IS SET EQUAL TO IROW.                DTSG 690
      IF(N-IROW)3,3,2                                                   DTSG 700
    2 N=IROW                                                            DTSG 710
C                                                                       DTSG 720
C     GENERATION OF VECTOR WORK AND COMPUTATION OF ITS GREATEST ELEMENT.DTSG 730
    3 B=0.D0                                                            DTSG 740
      DO 5 I=1,IROW                                                     DTSG 750
      DELTA=DABS(Z(I)-X)                                                DTSG 760
      IF(DELTA-B)5,5,4                                                  DTSG 770
    4 B=DELTA                                                           DTSG 780
    5 WORK(I)=DELTA                                                     DTSG 790
C                                                                       DTSG 800
C     GENERATION OF TABLE (ARG,VAL)                                     DTSG 810
      B=B+1.D0                                                          DTSG 820
      DO 10 J=1,N                                                       DTSG 830
      DELTA=B                                                           DTSG 840
      DO 7 I=1,IROW                                                     DTSG 850
      IF(WORK(I)-DELTA)6,7,7                                            DTSG 860
    6 II=I                                                              DTSG 870
      DELTA=WORK(I)                                                     DTSG 880
    7 CONTINUE                                                          DTSG 890
      ARG(J)=Z(II)                                                      DTSG 900
      IF(ICOL-1)8,9,8                                                   DTSG 910
    8 VAL(2*J-1)=F(II)                                                  DTSG 920
      III=II+IROW                                                       DTSG 930
      VAL(2*J)=F(III)                                                   DTSG 940
      GOTO 10                                                           DTSG 950
    9 VAL(J)=F(II)                                                      DTSG 960
   10 WORK(II)=B                                                        DTSG 970
   11 RETURN                                                            DTSG 980
      END                                                               DTSG 990
