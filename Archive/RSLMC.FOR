C                                                                       RSLM  10
C     ..................................................................RSLM  20
C                                                                       RSLM  30
C        SUBROUTINE RSLMC                                               RSLM  40
C                                                                       RSLM  50
C        PURPOSE                                                        RSLM  60
C           SOLUTION OF A SYSTEM OF LINEAR EQUATIONS AX=B               RSLM  70
C                                                                       RSLM  80
C        USAGE                                                          RSLM  90
C           CALL RSLMC(A,AF,B,X,N,EPSI,IER,IA,V,PER)                    RSLM 100
C                                                                       RSLM 110
C        DESCRIPTION OF PARAMETERS                                      RSLM 120
C           A      INPUT MATRIX                                         RSLM 130
C           AF     ARRAY OF THE FACTORIZATION OF THE ORIGINAL MATRIX    RSLM 140
C           B      RIGHT HAND SIDE VECTOR                               RSLM 150
C           X      VECTOR CONTAINING THE SOLUTION ON RETURN             RSLM 160
C           N      ORDER OF THE SYSTEM                                  RSLM 170
C           EPSI   RELATIVE PRECISION INDICATOR(REQUIRED INPUT)         RSLM 180
C           IER    ERROR INDICATOR                                      RSLM 190
C                     =0 IF EACH COMPONENT OF X MEETS THE PRECISION EPSIRSLM 200
C                     =1 IF ONLY THE NORM OF X MEETS THIS PRECISION     RSLM 210
C                     =2 IF THE PRECISION IN THE NORM OF THE COMPUTED   RSLM 220
C                        SOLUTION IS LOWER THAN EPSI                    RSLM 230
C                     =3 IF THE SOLUTION OBTAINED HAS NO MEANING AT ALL RSLM 240
C                     =4 IF A DIAGONAL TERM OF THE UPPER TRIANGULAR     RSLM 250
C                        FACTOR IS ZERO                                 RSLM 260
C           IA     SIZE OF THE FIRST DIMENSION ASSIGNED TO THE ARRAY A  RSLM 270
C                  IN THE CALLING PROGRAM WHEN THE MATRIX IS IN DOUBLE  RSLM 280
C                  SUBSCRIPTED DATA STORAGE MODE.  IA=N WHEN THE MATRIX RSLM 290
C                  IS IN SSP VECTOR STORAGE MODE.                       RSLM 300
C           V      WORKING STORAGE VECTOR                               RSLM 310
C                  DIMENSION OF V MUST BE GREATER THAN OR EQUAL TO N    RSLM 320
C           PER    VECTOR WHERE PERMUTATIONS OF ROWS OF THE MATRIX ARE  RSLM 330
C                  STORED                                               RSLM 340
C                  DIMENSION OF PER MUST BE GREATER THAN OR EQUAL TO N  RSLM 350
C                                                                       RSLM 360
C        REMARKS                                                        RSLM 370
C           THE MATRIX OF THE SYSTEM MAY BE FACTORIZED BY THE SUBROUTINERSLM 380
C           FACTR IN THE ARRAY AF PRIOR TO ENTRY TO THIS SUBROUTINE.    RSLM 390
C           THE LOWER TRIANGULAR FACTOR MUST HAVE AN UNIT DIAGONAL.     RSLM 400
C           EPSI IS MODIFIED WHEN IER=2                                 RSLM 410
C                                                                       RSLM 420
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  RSLM 430
C           NONE                                                        RSLM 440
C                                                                       RSLM 450
C        METHOD                                                         RSLM 460
C           A TRIAL SOLUTION IS FIRST COMPUTED.  THEN CORRECTIONS ARE   RSLM 470
C           CALCULATED FROM RESIDUAL VECTORS.                           RSLM 480
C                                                                       RSLM 490
C        REFERENCES                                                     RSLM 500
C           J. H. WILKINSON - THE ALGEBRAIC EIGENVALUE PROBLEM -        RSLM 510
C           CLARENDON PRESS OXFORD, 1965.  H. J. BOWDLER, R. S. MARTIN, RSLM 520
C           G. PETERS, AND J. H. WILKINSON.  'SOLUTION OF REAL AND      RSLM 530
C           COMPLEX SYSTEMS OF LINEAR EQUATIONS', NUMERISCHE MATHEMATIK,RSLM 540
C           VOL. 8, NO. 3, 1966, 217-234.                               RSLM 550
C                                                                       RSLM 560
C     ..................................................................RSLM 570
C                                                                       RSLM 580
      SUBROUTINE RSLMC (A,AF,B,X,N,EPSI,IER,IA,V,PER)                   RSLM 590
      DIMENSION A(1),AF(1),B(1),X(1),V(1),PER(1)                        RSLM 600
      DOUBLE PRECISION DP                                               RSLM 610
C                                                                       RSLM 620
C        INITIALIZATION                                                 RSLM 630
C                                                                       RSLM 640
      D0=0.                                                             RSLM 650
      IER=0                                                             RSLM 660
      ITE=0                                                             RSLM 670
      DO 10 I=1,N                                                       RSLM 680
      V(I)=B(I)                                                         RSLM 690
   10 X(I)=0.                                                           RSLM 700
   20 ITE=ITE+1                                                         RSLM 710
C                                                                       RSLM 720
C        THE PERMUTATIONS OF ROWS OF A ARE APPLIED TO V                 RSLM 730
C                                                                       RSLM 740
      DO 30 I=1,N                                                       RSLM 750
      K=PER(I)                                                          RSLM 760
      IF (K-I)25,30,25                                                  RSLM 770
   25 D1=V(K)                                                           RSLM 780
      V(K)=V(I)                                                         RSLM 790
      V(I)=D1                                                           RSLM 800
   30 CONTINUE                                                          RSLM 810
C                                                                       RSLM 820
C        SOLUTION OF THE LOWER TRIANGULAR SYSTEM                        RSLM 830
C                                                                       RSLM 840
      DO 50 I=2,N                                                       RSLM 850
      IM1=I-1                                                           RSLM 860
      DP=V(I)                                                           RSLM 870
      IK=I                                                              RSLM 880
      DO 40 K=1,IM1                                                     RSLM 890
      DP=DP-1.D0*AF(IK)*V(K)                                            RSLM 900
   40 IK=IK+IA                                                          RSLM 910
   50 V(I)=DP                                                           RSLM 920
C                                                                       RSLM 930
C        SOLUTION OF THE UPPER TRIANGULAR SYSTEM                        RSLM 940
C                                                                       RSLM 950
      IF(AF(IK)) 58,54,58                                               RSLM 960
   54 IER=4                                                             RSLM 970
      GO TO 82                                                          RSLM 980
   58 V(N)=DP/AF(IK)                                                    RSLM 990
      DO 70 I=2,N                                                       RSLM1000
      IM1=N-I+1                                                         RSLM1010
      INF=IM1+1                                                         RSLM1020
      DP=V(IM1)                                                         RSLM1030
      IK=(IM1-1)*IA+IM1                                                 RSLM1040
      D1=AF(IK)                                                         RSLM1050
      DO 60 K=INF,N                                                     RSLM1060
      IK=IK+IA                                                          RSLM1070
   60 DP=DP-1.D0*AF(IK)*V(K)                                            RSLM1080
   70 V(IM1)=DP/D1                                                      RSLM1090
C                                                                       RSLM1100
C        TEST OF PRECISION                                              RSLM1110
C                                                                       RSLM1120
      D1=0.                                                             RSLM1130
      D2=0.                                                             RSLM1140
      KLE=0                                                             RSLM1150
      DO 80 I=1,N                                                       RSLM1160
      D1=D1+ABS(V(I))                                                   RSLM1170
      D2=D2+ABS(X(I))                                                   RSLM1180
      IF (ABS(V(I))-EPSI*ABS(X(I))) 80,80,75                            RSLM1190
   75 KLE=1                                                             RSLM1200
   80 CONTINUE                                                          RSLM1210
      IF (KLE)140,82,85                                                 RSLM1220
   82 RETURN                                                            RSLM1230
   85 IF (ITE-1)140,90,87                                               RSLM1240
C                                                                       RSLM1250
C        ITERATIONS ARE STOPPED WHEN THE NORM OF THE CORRECTION IS MORE RSLM1260
C        THAN HALF OF THE ONE OF THE FORMER                             RSLM1270
C                                                                       RSLM1280
   87 IF (D0-2.*D1)120,90,90                                            RSLM1290
   90 DO 95 I=1,N                                                       RSLM1300
   95 X(I)=X(I)+V(I)                                                    RSLM1310
      DO 110 I=1,N                                                      RSLM1320
      DP=B(I)                                                           RSLM1330
      IK=I                                                              RSLM1340
      DO 100 K=1,N                                                      RSLM1350
      DP=DP-1.D0*A(IK)*X(K)                                             RSLM1360
  100 IK=IK+IA                                                          RSLM1370
  110 V(I)=DP                                                           RSLM1380
      D0=D1                                                             RSLM1390
      GO TO 20                                                          RSLM1400
  120 IF(ITE-2)140,140,125                                              RSLM1410
  125 IF (D1-EPSI*D2)127,127,130                                        RSLM1420
  127 IER=1                                                             RSLM1430
      RETURN                                                            RSLM1440
  130 IER=2                                                             RSLM1450
      EPSI=D1/D2                                                        RSLM1460
      RETURN                                                            RSLM1470
  140 IER=3                                                             RSLM1480
      RETURN                                                            RSLM1490
      END                                                               RSLM1500
