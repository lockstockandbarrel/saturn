C                                                                       MFSD  10
C     ..................................................................MFSD  20
C                                                                       MFSD  30
C        SUBROUTINE MFSD                                                MFSD  40
C                                                                       MFSD  50
C        PURPOSE                                                        MFSD  60
C           FACTOR A GIVEN SYMMETRIC POSITIVE DEFINITE MATRIX           MFSD  70
C                                                                       MFSD  80
C        USAGE                                                          MFSD  90
C           CALL MFSD(A,N,EPS,IER)                                      MFSD 100
C                                                                       MFSD 110
C        DESCRIPTION OF PARAMETERS                                      MFSD 120
C           A      - UPPER TRIANGULAR PART OF THE GIVEN SYMMETRIC       MFSD 130
C                    POSITIVE DEFINITE N BY N COEFFICIENT MATRIX.       MFSD 140
C                    ON RETURN A CONTAINS THE RESULTANT UPPER           MFSD 150
C                    TRIANGULAR MATRIX.                                 MFSD 160
C           N      - THE NUMBER OF ROWS (COLUMNS) IN GIVEN MATRIX.      MFSD 170
C           EPS    - AN INPUT CONSTANT WHICH IS USED AS RELATIVE        MFSD 180
C                    TOLERANCE FOR TEST ON LOSS OF SIGNIFICANCE.        MFSD 190
C           IER    - RESULTING ERROR PARAMETER CODED AS FOLLOWS         MFSD 200
C                    IER=0  - NO ERROR                                  MFSD 210
C                    IER=-1 - NO RESULT BECAUSE OF WRONG INPUT PARAME-  MFSD 220
C                             TER N OR BECAUSE SOME RADICAND IS NON-    MFSD 230
C                             POSITIVE (MATRIX A IS NOT POSITIVE        MFSD 240
C                             DEFINITE, POSSIBLY DUE TO LOSS OF SIGNI-  MFSD 250
C                             FICANCE)                                  MFSD 260
C                    IER=K  - WARNING WHICH INDICATES LOSS OF SIGNIFI-  MFSD 270
C                             CANCE. THE RADICAND FORMED AT FACTORIZA-  MFSD 280
C                             TION STEP K+1 WAS STILL POSITIVE BUT NO   MFSD 290
C                             LONGER GREATER THAN ABS(EPS*A(K+1,K+1)).  MFSD 300
C                                                                       MFSD 310
C        REMARKS                                                        MFSD 320
C           THE UPPER TRIANGULAR PART OF GIVEN MATRIX IS ASSUMED TO BE  MFSD 330
C           STORED COLUMNWISE IN N*(N+1)/2 SUCCESSIVE STORAGE LOCATIONS.MFSD 340
C           IN THE SAME STORAGE LOCATIONS THE RESULTING UPPER TRIANGU-  MFSD 350
C           LAR MATRIX IS STORED COLUMNWISE TOO.                        MFSD 360
C           THE PROCEDURE GIVES RESULTS IF N IS GREATER THAN 0 AND ALL  MFSD 370
C           CALCULATED RADICANDS ARE POSITIVE.                          MFSD 380
C           THE PRODUCT OF RETURNED DIAGONAL TERMS IS EQUAL TO THE      MFSD 390
C           SQUARE-ROOT OF THE DETERMINANT OF THE GIVEN MATRIX.         MFSD 400
C                                                                       MFSD 410
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  MFSD 420
C           NONE                                                        MFSD 430
C                                                                       MFSD 440
C        METHOD                                                         MFSD 450
C           SOLUTION IS DONE USING THE SQUARE-ROOT METHOD OF CHOLESKY.  MFSD 460
C           THE GIVEN MATRIX IS REPRESENTED AS PRODUCT OF TWO TRIANGULARMFSD 470
C           MATRICES, WHERE THE LEFT HAND FACTOR IS THE TRANSPOSE OF    MFSD 480
C           THE RETURNED RIGHT HAND FACTOR.                             MFSD 490
C                                                                       MFSD 500
C     ..................................................................MFSD 510
C                                                                       MFSD 520
      SUBROUTINE MFSD(A,N,EPS,IER)                                      MFSD 530
C                                                                       MFSD 540
C                                                                       MFSD 550
      DIMENSION A(1)                                                    MFSD 560
      DOUBLE PRECISION DPIV,DSUM                                        MFSD 570
C                                                                       MFSD 580
C        TEST ON WRONG INPUT PARAMETER N                                MFSD 590
      IF(N-1) 12,1,1                                                    MFSD 600
    1 IER=0                                                             MFSD 610
C                                                                       MFSD 620
C        INITIALIZE DIAGONAL-LOOP                                       MFSD 630
      KPIV=0                                                            MFSD 640
      DO 11 K=1,N                                                       MFSD 650
      KPIV=KPIV+K                                                       MFSD 660
      IND=KPIV                                                          MFSD 670
      LEND=K-1                                                          MFSD 680
C                                                                       MFSD 690
C        CALCULATE TOLERANCE                                            MFSD 700
      TOL=ABS(EPS*A(KPIV))                                              MFSD 710
C                                                                       MFSD 720
C        START FACTORIZATION-LOOP OVER K-TH ROW                         MFSD 730
      DO 11 I=K,N                                                       MFSD 740
      DSUM=0.D0                                                         MFSD 750
      IF(LEND) 2,4,2                                                    MFSD 760
C                                                                       MFSD 770
C        START INNER LOOP                                               MFSD 780
    2 DO 3 L=1,LEND                                                     MFSD 790
      LANF=KPIV-L                                                       MFSD 800
      LIND=IND-L                                                        MFSD 810
    3 DSUM=DSUM+DBLE(A(LANF)*A(LIND))                                   MFSD 820
C        END OF INNER LOOP                                              MFSD 830
C                                                                       MFSD 840
C        TRANSFORM ELEMENT A(IND)                                       MFSD 850
    4 DSUM=DBLE(A(IND))-DSUM                                            MFSD 860
      IF(I-K) 10,5,10                                                   MFSD 870
C                                                                       MFSD 880
C        TEST FOR NEGATIVE PIVOT ELEMENT AND FOR LOSS OF SIGNIFICANCE   MFSD 890
    5 IF(SNGL(DSUM)-TOL) 6,6,9                                          MFSD 900
    6 IF(DSUM) 12,12,7                                                  MFSD 910
    7 IF(IER) 8,8,9                                                     MFSD 920
    8 IER=K-1                                                           MFSD 930
C                                                                       MFSD 940
C        COMPUTE PIVOT ELEMENT                                          MFSD 950
    9 DPIV=DSQRT(DSUM)                                                  MFSD 960
      A(KPIV)=DPIV                                                      MFSD 970
      DPIV=1.D0/DPIV                                                    MFSD 980
      GO TO 11                                                          MFSD 990
C                                                                       MFSD1000
C        CALCULATE TERMS IN ROW                                         MFSD1010
   10 A(IND)=DSUM*DPIV                                                  MFSD1020
   11 IND=IND+I                                                         MFSD1030
C                                                                       MFSD1040
C        END OF DIAGONAL-LOOP                                           MFSD1050
      RETURN                                                            MFSD1060
   12 IER=-1                                                            MFSD1070
      RETURN                                                            MFSD1080
      END                                                               MFSD1090
