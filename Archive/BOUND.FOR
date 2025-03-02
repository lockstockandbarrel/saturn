C                                                                       BOUN  10
C     ..................................................................BOUN  20
C                                                                       BOUN  30
C        SUBROUTINE BOUND                                               BOUN  40
C                                                                       BOUN  50
C        PURPOSE                                                        BOUN  60
C           SELECT FROM A SET (OR A SUBSET) OF OBSERVATIONS THE NUMBER  BOUN  70
C           OF OBSERVATIONS UNDER, BETWEEN AND OVER TWO GIVEN BOUNDS    BOUN  80
C           FOR EACH VARIABLE                                           BOUN  90
C                                                                       BOUN 100
C        USAGE                                                          BOUN 110
C           CALL BOUND (A,S,BLO,BHI,UNDER,BETW,OVER,NO,NV,IER)          BOUN 120
C                                                                       BOUN 130
C        DESCRIPTION OF PARAMETERS                                      BOUN 140
C           A     - OBSERVATION MATRIX, NO BY NV                        BOUN 150
C           S     - VECTOR INDICATING SUBSET OF A. ONLY THOSE           BOUN 160
C                   OBSERVATIONS WITH A NON-ZERO S(J) ARE CONSIDERED.   BOUN 170
C                   VECTOR LENGTH IS NO.                                BOUN 180
C           BLO   - INPUT VECTOR OF LOWER BOUNDS ON ALL VARIABLES.      BOUN 190
C                   VECTOR LENGTH IS NV.                                BOUN 200
C           BHI   - INPUT VECTOR OF UPPER BOUNDS ON ALL VARIABLES.      BOUN 210
C                   VECTOR LENGTH IS NV.                                BOUN 220
C           UNDER - OUTPUT VECTOR INDICATING, FOR EACH VARIABLE, NUMBER BOUN 230
C                   OF OBSERVATIONS UNDER LOWER BOUNDS. VECTOR LENGTH   BOUN 240
C                   IS NV.                                              BOUN 250
C           BETW  - OUTPUT VECTOR INDICATING, FOR EACH VARIABLE, NUMBER BOUN 260
C                   OF OBSERVATIONS EQUAL TO OR BETWEEN LOWER AND UPPER BOUN 270
C                   BOUNDS. VECTOR LENGTH IS NV.                        BOUN 280
C           OVER  - OUTPUT VECTOR INDICATING, FOR EACH VARIABLE, NUMBER BOUN 290
C                   OF OBSERVATIONS OVER UPPER BOUNDS. VECTOR LENGTH    BOUN 300
C                   IS NV.                                              BOUN 310
C           NO    - NUMBER OF OBSERVATIONS                              BOUN 320
C           NV    - NUMBER OF VARIABLES FOR EACH OBSERVATION            BOUN 330
C           IER   - ZERO, IF NO ERROR.                                  BOUN 331
C                 - 1, IF LOWER BOUND IS GREATER THAN THE UPPER BOUND   BOUN 332
C                   FOR SOME VARIABLE                                   BOUN 333
C                                                                       BOUN 340
C        REMARKS                                                        BOUN 350
C           NONE                                                        BOUN 360
C                                                                       BOUN 370
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  BOUN 380
C           NONE                                                        BOUN 390
C                                                                       BOUN 400
C        METHOD                                                         BOUN 410
C           EACH ROW (OBSERVATION) OF MATRIX A WITH CORRESPONDING       BOUN 420
C           NON-ZERO ELEMENT IN S VECTOR IS TESTED. OBSERVATIONS ARE    BOUN 430
C           COMPARED WITH SPECIFIED LOWER AND UPPER VARIABLE BOUNDS AND BOUN 440
C           A COUNT IS KEPT IN VECTORS UNDER, BETWEEN, AND OVER.        BOUN 450
C                                                                       BOUN 460
C     ..................................................................BOUN 470
C                                                                       BOUN 480
      SUBROUTINE BOUND(A,S,BLO,BHI,UNDER,BETW,OVER,NO,NV)               BOUN 490
      DIMENSION A(1),S(1),BLO(1),BHI(1),UNDER(1),BETW(1),OVER(1)        BOUN 500
C                                                                       BOUN 510
C        CLEAR OUTPUT VECTORS.                                          BOUN 520
C                                                                       BOUN 530
      IER=0                                                             BOUN 531
      DO 10 I=1,NV                                                      BOUN 532
      IF (BLO(I)-BHI(I)) 10,10,11                                       BOUN 533
   11 IER=1                                                             BOUN 534
      GO TO 12                                                          BOUN 535
   10 CONTINUE                                                          BOUN 536
      DO 1 K=1,NV                                                       BOUN 540
      UNDER(K)=0.0                                                      BOUN 550
      BETW(K)=0.0                                                       BOUN 560
    1 OVER(K)=0.0                                                       BOUN 570
C                                                                       BOUN 580
C        TEST SUBSET VECTOR                                             BOUN 590
C                                                                       BOUN 600
      DO 8 J=1,NO                                                       BOUN 610
      IJ=J-NO                                                           BOUN 620
      IF(S(J)) 2,8,2                                                    BOUN 630
C                                                                       BOUN 640
C        COMPARE OBSERVATIONS WITH BOUNDS                               BOUN 650
C                                                                       BOUN 660
    2 DO 7 I=1,NV                                                       BOUN 670
      IJ=IJ+NO                                                          BOUN 680
      IF(A(IJ)-BLO(I)) 5,3,3                                            BOUN 690
    3 IF(A(IJ)-BHI(I)) 4,4,6                                            BOUN 700
C                                                                       BOUN 710
C       COUNT                                                           BOUN 720
C                                                                       BOUN 730
    4 BETW(I)=BETW(I)+1.0                                               BOUN 740
      GO TO 7                                                           BOUN 750
    5 UNDER(I)=UNDER(I)+1.0                                             BOUN 760
      GO TO 7                                                           BOUN 770
    6 OVER(I)=OVER(I)+1.0                                               BOUN 780
    7 CONTINUE                                                          BOUN 790
    8 CONTINUE                                                          BOUN 800
   12 RETURN                                                            BOUN 810
      END                                                               BOUN 820
