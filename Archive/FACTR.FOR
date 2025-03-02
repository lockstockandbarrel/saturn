C                                                                       FCTR  10
C     ..................................................................FCTR  20
C                                                                       FCTR  30
C        SUBROUTINE FACTR                                               FCTR  40
C                                                                       FCTR  50
C        PURPOSE                                                        FCTR  60
C           FACTORIZATION OF THE MATRIX A INTO A PRODUCT OF A LOWER     FCTR  70
C           TRIANGULAR MATRIX L AND AN UPPER TRIANGULAR MATRIX U.  L HASFCTR  80
C           UNIT DIAGONAL WHICH IS NOT STORED.                          FCTR  90
C                                                                       FCTR 100
C        USAGE                                                          FCTR 110
C           CALL FACTR(A,PER,N,IA,IER)                                  FCTR 120
C                                                                       FCTR 130
C        DESCRIPTION OF PARAMETERS                                      FCTR 140
C           A      MATRIX A                                             FCTR 150
C           PER    ONE DIMENSIONAL ARRAY WHERE PERMUTATIONS OF ROWS OF  FCTR 160
C                  THE MATRIX ARE STORED                                FCTR 170
C                  DIMENSION OF PER MUST BE GREATER THAN OR EQUAL TO N  FCTR 180
C           N      ORDER OF THE MATRIX A                                FCTR 190
C           IA     SIZE OF THE FIRST DIMENSION ASSIGNED TO THE ARRAY A  FCTR 200
C                  IN THE CALLING PROGRAM WHEN THE MATRIX IS IN DOUBLE  FCTR 210
C                  SUBSCRIPTED DATA STORAGE MODE.  IA=N WHEN THE MATRIX FCTR 220
C                  IS IN SSP VECTOR STORAGE MODE.                       FCTR 230
C           IER    ERROR INDICATOR WHICH IS ZERO IF THERE IS NO ERROR,  FCTR 240
C                  AND IS THREE IF THE PROCEDURE FAILS.                 FCTR 250
C                                                                       FCTR 260
C        REMARKS                                                        FCTR 270
C           THE ORIGINAL MATRIX, A,IS REPLACED BY THE TRIANGULAR FACTORSFCTR 280
C                                                                       FCTR 290
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  FCTR 300
C           NONE                                                        FCTR 310
C                                                                       FCTR 320
C        METHOD                                                         FCTR 330
C           SUCCESSIVE COMPUTATION OF THE COLUMNS OF L AND THE          FCTR 340
C           CORRESPONDING ROWS OF U.                                    FCTR 350
C                                                                       FCTR 360
C        REFERENCES                                                     FCTR 370
C           J. H. WILKINSON - THE ALGEBRAIC EIGENVALUE PROBLEM -        FCTR 380
C           CLARENDON PRESS, OXFORD, 1965. H. J. BOWDLER, R. S. MARTIN, FCTR 390
C           G. PETERS, AND J. H. WILKINSON - 'SOLUTION OF REAL AND      FCTR 400
C           COMPLEX SYSTEMS OF LINEAR EQUATIONS', NUMERISCHE MATHEMATIK,FCTR 410
C           VOL. 8, NO. 3, 1966, P. 217-234.                            FCTR 420
C                                                                       FCTR 430
C     ..................................................................FCTR 440
C                                                                       FCTR 450
      SUBROUTINE FACTR(A,PER,N,IA,IER)                                  FCTR 460
      DIMENSION A(1),PER(1)                                             FCTR 470
      DOUBLE PRECISION DP                                               FCTR 480
C                                                                       FCTR 490
C        COMPUTATION OF WEIGHTS FOR EQUILIBRATION                       FCTR 500
C                                                                       FCTR 510
      DO 20 I=1,N                                                       FCTR 520
      X=0.                                                              FCTR 530
      IJ=I                                                              FCTR 540
      DO 10 J=1,N                                                       FCTR 550
      IF (ABS(A(IJ))-X)10,10,5                                          FCTR 560
    5 X=ABS(A(IJ))                                                      FCTR 570
   10 IJ=IJ+IA                                                          FCTR 580
      IF (X) 110,110,20                                                 FCTR 590
   20 PER(I)=1./X                                                       FCTR 600
      I0=0                                                              FCTR 610
      DO 100 I=1,N                                                      FCTR 620
      IM1=I-1                                                           FCTR 630
      IP1=I+1                                                           FCTR 640
      IPIVOT=I                                                          FCTR 650
      X=0.                                                              FCTR 660
C                                                                       FCTR 670
C        COMPUTATION OF THE ITH COLUMN OF L                             FCTR 680
C                                                                       FCTR 690
      DO 50 K=I,N                                                       FCTR 700
      KI=I0+K                                                           FCTR 710
      DP=A(KI)                                                          FCTR 720
      IF (I-1) 110,40,25                                                FCTR 730
   25 KJ=K                                                              FCTR 740
      DO 30 J=1,IM1                                                     FCTR 750
      IJ=I0+J                                                           FCTR 760
      DP=DP-1.D0*A(KJ)*A(IJ)                                            FCTR 770
   30 KJ=KJ+IA                                                          FCTR 780
      A(KI)=DP                                                          FCTR 790
C                                                                       FCTR 800
C        SEARCH FOR EQUILIBRATED PIVOT                                  FCTR 810
C                                                                       FCTR 820
   40 IF (X-DABS(DP)*PER(K))45,50,50                                    FCTR 830
   45 IPIVOT=K                                                          FCTR 840
      X=DABS(DP)*PER(K)                                                 FCTR 850
   50 CONTINUE                                                          FCTR 860
      IF (X)110,110,55                                                  FCTR 870
C                                                                       FCTR 880
C        PERMUTATION OF ROWS IF REQUIRED                                FCTR 890
C                                                                       FCTR 900
   55 IF (IPIVOT-I) 110,70,57                                           FCTR 910
   57 KI=IPIVOT                                                         FCTR 920
      IJ=I                                                              FCTR 930
      DO 60 J=1,N                                                       FCTR 940
      X=A(IJ)                                                           FCTR 950
      A(IJ)=A(KI)                                                       FCTR 960
      A(KI)=X                                                           FCTR 970
      KI=KI+IA                                                          FCTR 980
   60 IJ=IJ+IA                                                          FCTR 990
      PER(IPIVOT)=PER(I)                                                FCTR1000
   70 PER(I)=IPIVOT                                                     FCTR1010
      IF (I-N) 72,100,100                                               FCTR1020
   72 IJ=I0+I                                                           FCTR1030
      X=A(IJ)                                                           FCTR1040
C                                                                       FCTR1050
C        COMPUTATION OF THE ITH ROW OF U                                FCTR1060
C                                                                       FCTR1070
      K0=I0+IA                                                          FCTR1080
      DO 90 K=IP1,N                                                     FCTR1090
      KI=I0+K                                                           FCTR1100
      A(KI)=A(KI)/X                                                     FCTR1110
      IF (I-1)110,90,75                                                 FCTR1120
   75 IJ=I                                                              FCTR1130
      KI=K0+I                                                           FCTR1140
      DP=A(KI)                                                          FCTR1150
      DO 80 J=1,IM1                                                     FCTR1160
      KJ=K0+J                                                           FCTR1170
      DP=DP-1.D0*A(IJ)*A(KJ)                                            FCTR1180
   80 IJ=IJ+IA                                                          FCTR1190
      A(KI)=DP                                                          FCTR1200
   90 K0=K0+IA                                                          FCTR1210
  100 I0=I0+IA                                                          FCTR1220
      IER=0                                                             FCTR1230
      RETURN                                                            FCTR1240
  110 IER=3                                                             FCTR1250
      RETURN                                                            FCTR1260
      END                                                               FCTR1270
