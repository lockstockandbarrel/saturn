C                                                                       SINV  10
C     ..................................................................SINV  20
C                                                                       SINV  30
C        SUBROUTINE SINV                                                SINV  40
C                                                                       SINV  50
C        PURPOSE                                                        SINV  60
C           INVERT A GIVEN SYMMETRIC POSITIVE DEFINITE MATRIX           SINV  70
C                                                                       SINV  80
C        USAGE                                                          SINV  90
C           CALL SINV(A,N,EPS,IER)                                      SINV 100
C                                                                       SINV 110
C        DESCRIPTION OF PARAMETERS                                      SINV 120
C           A      - UPPER TRIANGULAR PART OF THE GIVEN SYMMETRIC       SINV 130
C                    POSITIVE DEFINITE N BY N COEFFICIENT MATRIX.       SINV 140
C                    ON RETURN A CONTAINS THE RESULTANT UPPER           SINV 150
C                    TRIANGULAR MATRIX.                                 SINV 160
C           N      - THE NUMBER OF ROWS (COLUMNS) IN GIVEN MATRIX.      SINV 170
C           EPS    - AN INPUT CONSTANT WHICH IS USED AS RELATIVE        SINV 180
C                    TOLERANCE FOR TEST ON LOSS OF SIGNIFICANCE.        SINV 190
C           IER    - RESULTING ERROR PARAMETER CODED AS FOLLOWS         SINV 200
C                    IER=0  - NO ERROR                                  SINV 210
C                    IER=-1 - NO RESULT BECAUSE OF WRONG INPUT PARAME-  SINV 220
C                             TER N OR BECAUSE SOME RADICAND IS NON-    SINV 230
C                             POSITIVE (MATRIX A IS NOT POSITIVE        SINV 240
C                             DEFINITE, POSSIBLY DUE TO LOSS OF SIGNI-  SINV 250
C                             FICANCE)                                  SINV 260
C                    IER=K  - WARNING WHICH INDICATES LOSS OF SIGNIFI-  SINV 270
C                             CANCE. THE RADICAND FORMED AT FACTORIZA-  SINV 280
C                             TION STEP K+1 WAS STILL POSITIVE BUT NO   SINV 290
C                             LONGER GREATER THAN ABS(EPS*A(K+1,K+1)).  SINV 300
C                                                                       SINV 310
C        REMARKS                                                        SINV 320
C           THE UPPER TRIANGULAR PART OF GIVEN MATRIX IS ASSUMED TO BE  SINV 330
C           STORED COLUMNWISE IN N*(N+1)/2 SUCCESSIVE STORAGE LOCATIONS.SINV 340
C           IN THE SAME STORAGE LOCATIONS THE RESULTING UPPER TRIANGU-  SINV 350
C           LAR MATRIX IS STORED COLUMNWISE TOO.                        SINV 360
C           THE PROCEDURE GIVES RESULTS IF N IS GREATER THAN 0 AND ALL  SINV 370
C           CALCULATED RADICANDS ARE POSITIVE.                          SINV 380
C                                                                       SINV 390
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  SINV 400
C           MFSD                                                        SINV 410
C                                                                       SINV 420
C        METHOD                                                         SINV 430
C           SOLUTION IS DONE USING THE FACTORIZATION BY SUBROUTINE MFSD.SINV 440
C                                                                       SINV 450
C     ..................................................................SINV 460
C                                                                       SINV 470
      SUBROUTINE SINV(A,N,EPS,IER)                                      SINV 480
C                                                                       SINV 490
C                                                                       SINV 500
      DIMENSION A(1)                                                    SINV 510
      DOUBLE PRECISION DIN,WORK                                         SINV 520
C                                                                       SINV 530
C        FACTORIZE GIVEN MATRIX BY MEANS OF SUBROUTINE MFSD             SINV 540
C        A = TRANSPOSE(T) * T                                           SINV 550
      CALL MFSD(A,N,EPS,IER)                                            SINV 560
      IF(IER) 9,1,1                                                     SINV 570
C                                                                       SINV 580
C        INVERT UPPER TRIANGULAR MATRIX T                               SINV 590
C        PREPARE INVERSION-LOOP                                         SINV 600
    1 IPIV=N*(N+1)/2                                                    SINV 610
      IND=IPIV                                                          SINV 620
C                                                                       SINV 630
C        INITIALIZE INVERSION-LOOP                                      SINV 640
      DO 6 I=1,N                                                        SINV 650
      DIN=1.D0/DBLE(A(IPIV))                                            SINV 660
      A(IPIV)=DIN                                                       SINV 670
      MIN=N                                                             SINV 680
      KEND=I-1                                                          SINV 690
      LANF=N-KEND                                                       SINV 700
      IF(KEND) 5,5,2                                                    SINV 710
    2 J=IND                                                             SINV 720
C                                                                       SINV 730
C        INITIALIZE ROW-LOOP                                            SINV 740
      DO 4 K=1,KEND                                                     SINV 750
      WORK=0.D0                                                         SINV 760
      MIN=MIN-1                                                         SINV 770
      LHOR=IPIV                                                         SINV 780
      LVER=J                                                            SINV 790
C                                                                       SINV 800
C        START INNER LOOP                                               SINV 810
      DO 3 L=LANF,MIN                                                   SINV 820
      LVER=LVER+1                                                       SINV 830
      LHOR=LHOR+L                                                       SINV 840
    3 WORK=WORK+DBLE(A(LVER)*A(LHOR))                                   SINV 850
C        END OF INNER LOOP                                              SINV 860
C                                                                       SINV 870
      A(J)=-WORK*DIN                                                    SINV 880
    4 J=J-MIN                                                           SINV 890
C        END OF ROW-LOOP                                                SINV 900
C                                                                       SINV 910
    5 IPIV=IPIV-MIN                                                     SINV 920
    6 IND=IND-1                                                         SINV 930
C        END OF INVERSION-LOOP                                          SINV 940
C                                                                       SINV 950
C        CALCULATE INVERSE(A) BY MEANS OF INVERSE(T)                    SINV 960
C        INVERSE(A) = INVERSE(T) * TRANSPOSE(INVERSE(T))                SINV 970
C        INITIALIZE MULTIPLICATION-LOOP                                 SINV 980
      DO 8 I=1,N                                                        SINV 990
      IPIV=IPIV+I                                                       SINV1000
      J=IPIV                                                            SINV1010
C                                                                       SINV1020
C        INITIALIZE ROW-LOOP                                            SINV1030
      DO 8 K=I,N                                                        SINV1040
      WORK=0.D0                                                         SINV1050
      LHOR=J                                                            SINV1060
C                                                                       SINV1070
C        START INNER LOOP                                               SINV1080
      DO 7 L=K,N                                                        SINV1090
      LVER=LHOR+K-I                                                     SINV1100
      WORK=WORK+DBLE(A(LHOR)*A(LVER))                                   SINV1110
    7 LHOR=LHOR+L                                                       SINV1120
C        END OF INNER LOOP                                              SINV1130
C                                                                       SINV1140
      A(J)=WORK                                                         SINV1150
    8 J=J+K                                                             SINV1160
C        END OF ROW- AND MULTIPLICATION-LOOP                            SINV1170
C                                                                       SINV1180
    9 RETURN                                                            SINV1190
      END                                                               SINV1200
