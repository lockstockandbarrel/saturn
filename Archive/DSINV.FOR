C                                                                       DSIN  10
C     ..................................................................DSIN  20
C                                                                       DSIN  30
C        SUBROUTINE DSINV                                               DSIN  40
C                                                                       DSIN  50
C        PURPOSE                                                        DSIN  60
C           INVERT A GIVEN SYMMETRIC POSITIVE DEFINITE MATRIX           DSIN  70
C                                                                       DSIN  80
C        USAGE                                                          DSIN  90
C           CALL DSINV(A,N,EPS,IER)                                     DSIN 100
C                                                                       DSIN 110
C        DESCRIPTION OF PARAMETERS                                      DSIN 120
C           A      - DOUBLE PRECISION UPPER TRIANGULAR PART OF GIVEN    DSIN 130
C                    SYMMETRIC POSITIVE DEFINITE N BY N COEFFICIENT     DSIN 140
C                    MATRIX.                                            DSIN 150
C                    ON RETURN A CONTAINS THE RESULTANT UPPER           DSIN 160
C                    TRIANGULAR MATRIX IN DOUBLE PRECISION.             DSIN 170
C           N      - THE NUMBER OF ROWS (COLUMNS) IN GIVEN MATRIX.      DSIN 180
C           EPS    - SINGLE PRECISION INPUT CONSTANT WHICH IS USED      DSIN 190
C                    AS RELATIVE TOLERANCE FOR TEST ON LOSS OF          DSIN 200
C                    SIGNIFICANCE.                                      DSIN 210
C           IER    - RESULTING ERROR PARAMETER CODED AS FOLLOWS         DSIN 220
C                    IER=0  - NO ERROR                                  DSIN 230
C                    IER=-1 - NO RESULT BECAUSE OF WRONG INPUT PARAME-  DSIN 240
C                             TER N OR BECAUSE SOME RADICAND IS NON-    DSIN 250
C                             POSITIVE (MATRIX A IS NOT POSITIVE        DSIN 260
C                             DEFINITE, POSSIBLY DUE TO LOSS OF SIGNI-  DSIN 270
C                             FICANCE)                                  DSIN 280
C                    IER=K  - WARNING WHICH INDICATES LOSS OF SIGNIFI-  DSIN 290
C                             CANCE. THE RADICAND FORMED AT FACTORIZA-  DSIN 300
C                             TION STEP K+1 WAS STILL POSITIVE BUT NO   DSIN 310
C                             LONGER GREATER THAN ABS(EPS*A(K+1,K+1)).  DSIN 320
C                                                                       DSIN 330
C        REMARKS                                                        DSIN 340
C           THE UPPER TRIANGULAR PART OF GIVEN MATRIX IS ASSUMED TO BE  DSIN 350
C           STORED COLUMNWISE IN N*(N+1)/2 SUCCESSIVE STORAGE LOCATIONS.DSIN 360
C           IN THE SAME STORAGE LOCATIONS THE RESULTING UPPER TRIANGU-  DSIN 370
C           LAR MATRIX IS STORED COLUMNWISE TOO.                        DSIN 380
C           THE PROCEDURE GIVES RESULTS IF N IS GREATER THAN 0 AND ALL  DSIN 390
C           CALCULATED RADICANDS ARE POSITIVE.                          DSIN 400
C                                                                       DSIN 410
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DSIN 420
C           DMFSD                                                       DSIN 430
C                                                                       DSIN 440
C        METHOD                                                         DSIN 450
C           SOLUTION IS DONE USING FACTORIZATION BY SUBROUTINE DMFSD.   DSIN 460
C                                                                       DSIN 470
C     ..................................................................DSIN 480
C                                                                       DSIN 490
      SUBROUTINE DSINV(A,N,EPS,IER)                                     DSIN 500
C                                                                       DSIN 510
C                                                                       DSIN 520
      DIMENSION A(1)                                                    DSIN 530
      DOUBLE PRECISION A,DIN,WORK                                       DSIN 540
C                                                                       DSIN 550
C        FACTORIZE GIVEN MATRIX BY MEANS OF SUBROUTINE DMFSD            DSIN 560
C        A = TRANSPOSE(T) * T                                           DSIN 570
      CALL DMFSD(A,N,EPS,IER)                                           DSIN 580
      IF(IER) 9,1,1                                                     DSIN 590
C                                                                       DSIN 600
C        INVERT UPPER TRIANGULAR MATRIX T                               DSIN 610
C        PREPARE INVERSION-LOOP                                         DSIN 620
    1 IPIV=N*(N+1)/2                                                    DSIN 630
      IND=IPIV                                                          DSIN 640
C                                                                       DSIN 650
C        INITIALIZE INVERSION-LOOP                                      DSIN 660
      DO 6 I=1,N                                                        DSIN 670
      DIN=1.D0/A(IPIV)                                                  DSIN 680
      A(IPIV)=DIN                                                       DSIN 690
      MIN=N                                                             DSIN 700
      KEND=I-1                                                          DSIN 710
      LANF=N-KEND                                                       DSIN 720
      IF(KEND) 5,5,2                                                    DSIN 730
    2 J=IND                                                             DSIN 740
C                                                                       DSIN 750
C        INITIALIZE ROW-LOOP                                            DSIN 760
      DO 4 K=1,KEND                                                     DSIN 770
      WORK=0.D0                                                         DSIN 780
      MIN=MIN-1                                                         DSIN 790
      LHOR=IPIV                                                         DSIN 800
      LVER=J                                                            DSIN 810
C                                                                       DSIN 820
C        START INNER LOOP                                               DSIN 830
      DO 3 L=LANF,MIN                                                   DSIN 840
      LVER=LVER+1                                                       DSIN 850
      LHOR=LHOR+L                                                       DSIN 860
    3 WORK=WORK+A(LVER)*A(LHOR)                                         DSIN 870
C        END OF INNER LOOP                                              DSIN 880
C                                                                       DSIN 890
      A(J)=-WORK*DIN                                                    DSIN 900
    4 J=J-MIN                                                           DSIN 910
C        END OF ROW-LOOP                                                DSIN 920
C                                                                       DSIN 930
    5 IPIV=IPIV-MIN                                                     DSIN 940
    6 IND=IND-1                                                         DSIN 950
C        END OF INVERSION-LOOP                                          DSIN 960
C                                                                       DSIN 970
C        CALCULATE INVERSE(A) BY MEANS OF INVERSE(T)                    DSIN 980
C        INVERSE(A) = INVERSE(T) * TRANSPOSE(INVERSE(T))                DSIN 990
C        INITIALIZE MULTIPLICATION-LOOP                                 DSIN1000
      DO 8 I=1,N                                                        DSIN1010
      IPIV=IPIV+I                                                       DSIN1020
      J=IPIV                                                            DSIN1030
C                                                                       DSIN1040
C        INITIALIZE ROW-LOOP                                            DSIN1050
      DO 8 K=I,N                                                        DSIN1060
      WORK=0.D0                                                         DSIN1070
      LHOR=J                                                            DSIN1080
C                                                                       DSIN1090
C        START INNER LOOP                                               DSIN1100
      DO 7 L=K,N                                                        DSIN1110
      LVER=LHOR+K-I                                                     DSIN1120
      WORK=WORK+A(LHOR)*A(LVER)                                         DSIN1130
    7 LHOR=LHOR+L                                                       DSIN1140
C        END OF INNER LOOP                                              DSIN1150
C                                                                       DSIN1160
      A(J)=WORK                                                         DSIN1170
    8 J=J+K                                                             DSIN1180
C        END OF ROW- AND MULTIPLICATION-LOOP                            DSIN1190
C                                                                       DSIN1200
    9 RETURN                                                            DSIN1210
      END                                                               DSIN1220
