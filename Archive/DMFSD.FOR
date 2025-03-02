C                                                                       DMSD  10
C     ..................................................................DMSD  20
C                                                                       DMSD  30
C        SUBROUTINE DMFSD                                               DMSD  40
C                                                                       DMSD  50
C        PURPOSE                                                        DMSD  60
C           FACTOR A GIVEN SYMMETRIC POSITIVE DEFINITE MATRIX           DMSD  70
C                                                                       DMSD  80
C        USAGE                                                          DMSD  90
C           CALL DMFSD(A,N,EPS,IER)                                     DMSD 100
C                                                                       DMSD 110
C        DESCRIPTION OF PARAMETERS                                      DMSD 120
C           A      - DOUBLE PRECISION UPPER TRIANGULAR PART OF GIVEN    DMSD 130
C                    SYMMETRIC POSITIVE DEFINITE N BY N COEFFICIENT     DMSD 140
C                    MATRIX.                                            DMSD 150
C                    ON RETURN A CONTAINS THE RESULTANT UPPER           DMSD 160
C                    TRIANGULAR MATRIX IN DOUBLE PRECISION.             DMSD 170
C           N      - THE NUMBER OF ROWS (COLUMNS) IN GIVEN MATRIX.      DMSD 180
C           EPS    - SINGLE PRECISION INPUT CONSTANT WHICH IS USED      DMSD 190
C                    AS RELATIVE TOLERANCE FOR TEST ON LOSS OF          DMSD 200
C                    SIGNIFICANCE.                                      DMSD 210
C           IER    - RESULTING ERROR PARAMETER CODED AS FOLLOWS         DMSD 220
C                    IER=0  - NO ERROR                                  DMSD 230
C                    IER=-1 - NO RESULT BECAUSE OF WRONG INPUT PARAME-  DMSD 240
C                             TER N OR BECAUSE SOME RADICAND IS NON-    DMSD 250
C                             POSITIVE (MATRIX A IS NOT POSITIVE        DMSD 260
C                             DEFINITE, POSSIBLY DUE TO LOSS OF SIGNI-  DMSD 270
C                             FICANCE)                                  DMSD 280
C                    IER=K  - WARNING WHICH INDICATES LOSS OF SIGNIFI-  DMSD 290
C                             CANCE. THE RADICAND FORMED AT FACTORIZA-  DMSD 300
C                             TION STEP K+1 WAS STILL POSITIVE BUT NO   DMSD 310
C                             LONGER GREATER THAN ABS(EPS*A(K+1,K+1)).  DMSD 320
C                                                                       DMSD 330
C        REMARKS                                                        DMSD 340
C           THE UPPER TRIANGULAR PART OF GIVEN MATRIX IS ASSUMED TO BE  DMSD 350
C           STORED COLUMNWISE IN N*(N+1)/2 SUCCESSIVE STORAGE LOCATIONS.DMSD 360
C           IN THE SAME STORAGE LOCATIONS THE RESULTING UPPER TRIANGU-  DMSD 370
C           LAR MATRIX IS STORED COLUMNWISE TOO.                        DMSD 380
C           THE PROCEDURE GIVES RESULTS IF N IS GREATER THAN 0 AND ALL  DMSD 390
C           CALCULATED RADICANDS ARE POSITIVE.                          DMSD 400
C           THE PRODUCT OF RETURNED DIAGONAL TERMS IS EQUAL TO THE      DMSD 410
C           SQUARE-ROOT OF THE DETERMINANT OF THE GIVEN MATRIX.         DMSD 420
C                                                                       DMSD 430
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DMSD 440
C           NONE                                                        DMSD 450
C                                                                       DMSD 460
C        METHOD                                                         DMSD 470
C           SOLUTION IS DONE USING THE SQUARE-ROOT METHOD OF CHOLESKY.  DMSD 480
C           THE GIVEN MATRIX IS REPRESENTED AS PRODUCT OF TWO TRIANGULARDMSD 490
C           MATRICES, WHERE THE LEFT HAND FACTOR IS THE TRANSPOSE OF    DMSD 500
C           THE RETURNED RIGHT HAND FACTOR.                             DMSD 510
C                                                                       DMSD 520
C     ..................................................................DMSD 530
C                                                                       DMSD 540
      SUBROUTINE DMFSD(A,N,EPS,IER)                                     DMSD 550
C                                                                       DMSD 560
C                                                                       DMSD 570
      DIMENSION A(1)                                                    DMSD 580
      DOUBLE PRECISION DPIV,DSUM,A                                      DMSD 590
C                                                                       DMSD 600
C        TEST ON WRONG INPUT PARAMETER N                                DMSD 610
      IF(N-1) 12,1,1                                                    DMSD 620
    1 IER=0                                                             DMSD 630
C                                                                       DMSD 640
C        INITIALIZE DIAGONAL-LOOP                                       DMSD 650
      KPIV=0                                                            DMSD 660
      DO 11 K=1,N                                                       DMSD 670
      KPIV=KPIV+K                                                       DMSD 680
      IND=KPIV                                                          DMSD 690
      LEND=K-1                                                          DMSD 700
C                                                                       DMSD 710
C        CALCULATE TOLERANCE                                            DMSD 720
      TOL=ABS(EPS*SNGL(A(KPIV)))                                        DMSD 730
C                                                                       DMSD 740
C        START FACTORIZATION-LOOP OVER K-TH ROW                         DMSD 750
      DO 11 I=K,N                                                       DMSD 760
      DSUM=0.D0                                                         DMSD 770
      IF(LEND) 2,4,2                                                    DMSD 780
C                                                                       DMSD 790
C        START INNER LOOP                                               DMSD 800
    2 DO 3 L=1,LEND                                                     DMSD 810
      LANF=KPIV-L                                                       DMSD 820
      LIND=IND-L                                                        DMSD 830
    3 DSUM=DSUM+A(LANF)*A(LIND)                                         DMSD 840
C        END OF INNER LOOP                                              DMSD 850
C                                                                       DMSD 860
C        TRANSFORM ELEMENT A(IND)                                       DMSD 870
    4 DSUM=A(IND)-DSUM                                                  DMSD 880
      IF(I-K) 10,5,10                                                   DMSD 890
C                                                                       DMSD 900
C        TEST FOR NEGATIVE PIVOT ELEMENT AND FOR LOSS OF SIGNIFICANCE   DMSD 910
    5 IF(SNGL(DSUM)-TOL) 6,6,9                                          DMSD 920
    6 IF(DSUM) 12,12,7                                                  DMSD 930
    7 IF(IER) 8,8,9                                                     DMSD 940
    8 IER=K-1                                                           DMSD 950
C                                                                       DMSD 960
C        COMPUTE PIVOT ELEMENT                                          DMSD 970
    9 DPIV=DSQRT(DSUM)                                                  DMSD 980
      A(KPIV)=DPIV                                                      DMSD 990
      DPIV=1.D0/DPIV                                                    DMSD1000
      GO TO 11                                                          DMSD1010
C                                                                       DMSD1020
C        CALCULATE TERMS IN ROW                                         DMSD1030
   10 A(IND)=DSUM*DPIV                                                  DMSD1040
   11 IND=IND+I                                                         DMSD1050
C        END OF DIAGONAL-LOOP                                           DMSD1060
C                                                                       DMSD1070
      RETURN                                                            DMSD1080
   12 IER=-1                                                            DMSD1090
      RETURN                                                            DMSD1100
      END                                                               DMSD1110
