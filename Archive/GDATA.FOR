C                                                                       GDAT  10
C     ..................................................................GDAT  20
C                                                                       GDAT  30
C        SUBROUTINE GDATA                                               GDAT  40
C                                                                       GDAT  50
C        PURPOSE                                                        GDAT  60
C           GENERATE INDEPENDENT VARIABLES UP TO THE M-TH POWER (THE    GDAT  70
C           HIGHEST DEGREE POLYNOMIAL SPECIFIED) AND COMPUTE MEANS,     GDAT  80
C           STANDARD DEVIATIONS, AND CORRELATION COEFFICIENTS.  THIS    GDAT  90
C           SUBROUTINE IS NORMALLY CALLED BEFORE SUBROUTINES ORDER,     GDAT 100
C           MINV AND MULTR IN THE PERFORMANCE OF A POLYNOMIAL           GDAT 110
C           REGRESSION.                                                 GDAT 120
C                                                                       GDAT 130
C        USAGE                                                          GDAT 140
C           CALL GDATA (N,M,X,XBAR,STD,D,SUMSQ)                         GDAT 150
C                                                                       GDAT 160
C        DESCRIPTION OF PARAMETERS                                      GDAT 170
C           N     - NUMBER OF OBSERVATIONS.                             GDAT 180
C           M     - THE HIGHEST DEGREE POLYNOMIAL TO BE FITTED.         GDAT 190
C           X     - INPUT MATRIX (N BY M+1) .  WHEN THE SUBROUTINE IS   GDAT 200
C                   CALLED, DATA FOR THE INDEPENDENT VARIABLE ARE       GDAT 210
C                   STORED IN THE FIRST COLUMN OF MATRIX X, AND DATA FORGDAT 220
C                   THE DEPENDENT VARIABLE ARE STORED IN THE LAST       GDAT 230
C                   COLUMN OF THE MATRIX.  UPON RETURNING TO THE        GDAT 240
C                   CALLING ROUTINE, GENERATED POWERS OF THE INDEPENDENTGDAT 250
C                   VARIABLE ARE STORED IN COLUMNS 2 THROUGH M.         GDAT 260
C           XBAR  - OUTPUT VECTOR OF LENGTH M+1 CONTAINING MEANS OF     GDAT 270
C                   INDEPENDENT AND DEPENDENT VARIABLES.                GDAT 280
C           STD   - OUTPUT VECTOR OF LENGTH M+1 CONTAINING STANDARD     GDAT 290
C                   DEVIATIONS OF INDEPENDENT AND DEPENDENT VARIABLES.  GDAT 300
C           D     - OUTPUT MATRIX (ONLY UPPER TRIANGULAR PORTION OF THE GDAT 310
C                   SYMMETRIC MATRIX OF M+1 BY M+1) CONTAINING CORRELA- GDAT 320
C                   TION COEFFICIENTS.  (STORAGE MODE OF 1)             GDAT 330
C           SUMSQ - OUTPUT VECTOR OF LENGTH M+1 CONTAINING SUMS OF      GDAT 340
C                   PRODUCTS OF DEVIATIONS FROM MEANS  OF INDEPENDENT   GDAT 350
C                   AND DEPENDENT VARIABLES.                            GDAT 360
C                                                                       GDAT 370
C        REMARKS                                                        GDAT 380
C           N MUST BE GREATER THAN M+1.                                 GDAT 390
C           IF M IS EQUAL TO 5 OR GREATER, SINGLE PRECISION MAY NOT BE  GDAT 400
C           SUFFICIENT TO GIVE SATISFACTORY COMPUTATIONAL RESULTS.      GDAT 410
C                                                                       GDAT 420
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  GDAT 430
C           NONE                                                        GDAT 440
C                                                                       GDAT 450
C        METHOD                                                         GDAT 460
C           REFER TO B. OSTLE, 'STATISTICS IN RESEARCH', THE IOWA STATE GDAT 470
C           COLLEGE PRESS, 1954, CHAPTER 6.                             GDAT 480
C                                                                       GDAT 490
C     ..................................................................GDAT 500
C                                                                       GDAT 510
      SUBROUTINE GDATA (N,M,X,XBAR,STD,D,SUMSQ)                         GDAT 520
      DIMENSION X(1),XBAR(1),STD(1),D(1),SUMSQ(1)                       GDAT 530
C                                                                       GDAT 540
C        ...............................................................GDAT 550
C                                                                       GDAT 560
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE  GDAT 570
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION      GDAT 580
C        STATEMENT WHICH FOLLOWS.                                       GDAT 590
C                                                                       GDAT 600
C     DOUBLE PRECISION X,XBAR,STD,D,SUMSQ,T1,T2                         GDAT 610
C                                                                       GDAT 620
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS    GDAT 630
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS      GDAT 640
C        ROUTINE.                                                       GDAT 650
C                                                                       GDAT 660
C        THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO      GDAT 670
C        CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  SQRT AND ABS IN   GDAT 680
C        STATEMENT 180 MUST BE CHANGED TO DSQRT AND DABS.               GDAT 690
C                                                                       GDAT 700
C        ...............................................................GDAT 710
C                                                                       GDAT 720
C     GENERATE INDEPENDENT VARIABLES                                    GDAT 730
C                                                                       GDAT 740
      IF(M-1) 105, 105, 90                                              GDAT 750
   90 L1=0                                                              GDAT 760
      DO 100 I=2,M                                                      GDAT 770
      L1=L1+N                                                           GDAT 780
      DO 100 J=1,N                                                      GDAT 790
      L=L1+J                                                            GDAT 800
      K=L-N                                                             GDAT 810
  100 X(L)=X(K)*X(J)                                                    GDAT 820
C                                                                       GDAT 830
C     CALCULATE MEANS                                                   GDAT 840
C                                                                       GDAT 850
  105 MM=M+1                                                            GDAT 860
      DF=N                                                              GDAT 870
      L=0                                                               GDAT 880
      DO 115 I=1,MM                                                     GDAT 890
      XBAR(I)=0.0                                                       GDAT 900
      DO 110 J=1,N                                                      GDAT 910
      L=L+1                                                             GDAT 920
  110 XBAR(I)=XBAR(I)+X(L)                                              GDAT 930
  115 XBAR(I)=XBAR(I)/DF                                                GDAT 940
C                                                                       GDAT 950
      DO 130 I=1,MM                                                     GDAT 960
  130 STD(I)=0.0                                                        GDAT 970
C                                                                       GDAT 980
C     CALCULATE SUMS OF CROSS-PRODUCTS OF DEVIATIONS                    GDAT 990
C                                                                       GDAT1000
      L=((MM+1)*MM)/2                                                   GDAT1010
      DO 150 I=1,L                                                      GDAT1020
  150 D(I)=0.0                                                          GDAT1030
      DO 170 K=1,N                                                      GDAT1040
      L=0                                                               GDAT1050
      DO 170 J=1,MM                                                     GDAT1060
      L2=N*(J-1)+K                                                      GDAT1070
      T2=X(L2)-XBAR(J)                                                  GDAT1080
      STD(J)=STD(J)+T2                                                  GDAT1090
      DO 170 I=1,J                                                      GDAT1100
      L1=N*(I-1)+K                                                      GDAT1110
      T1=X(L1)-XBAR(I)                                                  GDAT1120
      L=L+1                                                             GDAT1130
  170 D(L)=D(L)+T1*T2                                                   GDAT1140
      L=0                                                               GDAT1150
      DO 175 J=1,MM                                                     GDAT1160
      DO 175 I=1,J                                                      GDAT1170
      L=L+1                                                             GDAT1180
  175 D(L)=D(L)-STD(I)*STD(J)/DF                                        GDAT1190
      L=0                                                               GDAT1200
      DO 180 I=1,MM                                                     GDAT1210
      L=L+I                                                             GDAT1220
      SUMSQ(I)=D(L)                                                     GDAT1230
  180 STD(I)= SQRT( ABS(D(L)))                                          GDAT1240
C                                                                       GDAT1250
C     CALCULATE CORRELATION COEFFICIENTS                                GDAT1260
C                                                                       GDAT1270
      L=0                                                               GDAT1280
      DO 190 J=1,MM                                                     GDAT1290
      DO 190 I=1,J                                                      GDAT1300
      L=L+1                                                             GDAT1310
  190 D(L)=D(L)/(STD(I)*STD(J))                                         GDAT1320
C                                                                       GDAT1330
C     CALCULATE STANDARD DEVIATIONS                                     GDAT1340
C                                                                       GDAT1350
      DF=SQRT(DF-1.0)                                                   GDAT1360
      DO 200 I=1,MM                                                     GDAT1370
  200 STD(I)=STD(I)/DF                                                  GDAT1380
      RETURN                                                            GDAT1390
      END                                                               GDAT1400
