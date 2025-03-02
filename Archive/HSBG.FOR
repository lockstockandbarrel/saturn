C                                                                       HSBG  10
C     ..................................................................HSBG  20
C                                                                       HSBG  30
C        SUBROUTINE HSBG                                                HSBG  40
C                                                                       HSBG  50
C        PURPOSE                                                        HSBG  60
C           TO REDUCE A REAL MATRIX INTO UPPER ALMOST TRIANGULAR FORM   HSBG  70
C                                                                       HSBG  80
C        USAGE                                                          HSBG  90
C           CALL HSBG(N,A,IA)                                           HSBG 100
C                                                                       HSBG 110
C        DESCRIPTION OF THE PARAMETERS                                  HSBG 120
C           N      ORDER OF THE MATRIX                                  HSBG 130
C           A      THE INPUT MATRIX, N BY N                             HSBG 140
C           IA     SIZE OF THE FIRST DIMENSION ASSIGNED TO THE ARRAY    HSBG 150
C                  A IN THE CALLING PROGRAM WHEN THE MATRIX IS IN       HSBG 160
C                  DOUBLE SUBSCRIPTED DATA STORAGE MODE.  IA=N WHEN     HSBG 170
C                  THE MATRIX IS IN SSP VECTOR STORAGE MODE.            HSBG 180
C                                                                       HSBG 190
C        REMARKS                                                        HSBG 200
C           THE HESSENBERG FORM REPLACES THE ORIGINAL MATRIX IN THE     HSBG 210
C           ARRAY A.                                                    HSBG 220
C                                                                       HSBG 230
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  HSBG 240
C           NONE                                                        HSBG 250
C                                                                       HSBG 260
C        METHOD                                                         HSBG 270
C           SIMILARITY TRANSFORMATIONS USING ELEMENTARY ELIMINATION     HSBG 280
C           MATRICES, WITH PARTIAL PIVOTING.                            HSBG 290
C                                                                       HSBG 300
C        REFERENCES                                                     HSBG 310
C           J.H. WILKINSON - THE ALGEBRAIC EIGENVALUE PROBLEM -         HSBG 320
C           CLARENDON PRESS, OXFORD, 1965.                              HSBG 330
C                                                                       HSBG 340
C     ..................................................................HSBG 350
C                                                                       HSBG 360
      SUBROUTINE HSBG(N,A,IA)                                           HSBG 370
      DIMENSION A(1)                                                    HSBG 380
      DOUBLE PRECISION S                                                HSBG 390
      L=N                                                               HSBG 400
      NIA=L*IA                                                          HSBG 410
      LIA=NIA-IA                                                        HSBG 420
C                                                                       HSBG 430
C        L IS THE ROW INDEX OF THE ELIMINATION                          HSBG 440
C                                                                       HSBG 450
   20 IF(L-3) 360,40,40                                                 HSBG 460
   40 LIA=LIA-IA                                                        HSBG 470
      L1=L-1                                                            HSBG 480
      L2=L1-1                                                           HSBG 490
C                                                                       HSBG 500
C        SEARCH FOR THE PIVOTAL ELEMENT IN THE LTH ROW                  HSBG 510
C                                                                       HSBG 520
      ISUB=LIA+L                                                        HSBG 530
      IPIV=ISUB-IA                                                      HSBG 540
      PIV=ABS(A(IPIV))                                                  HSBG 550
      IF(L-3) 90,90,50                                                  HSBG 560
   50 M=IPIV-IA                                                         HSBG 570
      DO 80 I=L,M,IA                                                    HSBG 580
      T=ABS(A(I))                                                       HSBG 590
      IF(T-PIV) 80,80,60                                                HSBG 600
   60 IPIV=I                                                            HSBG 610
      PIV=T                                                             HSBG 620
   80 CONTINUE                                                          HSBG 630
   90 IF(PIV) 100,320,100                                               HSBG 640
  100 IF(PIV-ABS(A(ISUB))) 180,180,120                                  HSBG 650
C                                                                       HSBG 660
C        INTERCHANGE THE COLUMNS                                        HSBG 670
C                                                                       HSBG 680
  120 M=IPIV-L                                                          HSBG 690
      DO 140 I=1,L                                                      HSBG 700
      J=M+I                                                             HSBG 710
      T=A(J)                                                            HSBG 720
      K=LIA+I                                                           HSBG 730
      A(J)=A(K)                                                         HSBG 740
  140 A(K)=T                                                            HSBG 750
C                                                                       HSBG 760
C        INTERCHANGE THE ROWS                                           HSBG 770
C                                                                       HSBG 780
      M=L2-M/IA                                                         HSBG 790
      DO 160 I=L1,NIA,IA                                                HSBG 800
      T=A(I)                                                            HSBG 810
      J=I-M                                                             HSBG 820
      A(I)=A(J)                                                         HSBG 830
  160 A(J)=T                                                            HSBG 840
C                                                                       HSBG 850
C        TERMS OF THE ELEMENTARY TRANSFORMATION                         HSBG 860
C                                                                       HSBG 870
  180 DO 200 I=L,LIA,IA                                                 HSBG 880
  200 A(I)=A(I)/A(ISUB)                                                 HSBG 890
C                                                                       HSBG 900
C        RIGHT TRANSFORMATION                                           HSBG 910
C                                                                       HSBG 920
      J=-IA                                                             HSBG 930
      DO 240 I=1,L2                                                     HSBG 940
      J=J+IA                                                            HSBG 950
      LJ=L+J                                                            HSBG 960
      DO 220 K=1,L1                                                     HSBG 970
      KJ=K+J                                                            HSBG 980
      KL=K+LIA                                                          HSBG 990
  220 A(KJ)=A(KJ)-A(LJ)*A(KL)                                           HSBG1000
  240 CONTINUE                                                          HSBG1010
C                                                                       HSBG1020
C        LEFT TRANSFORMATION                                            HSBG1030
C                                                                       HSBG1040
      K=-IA                                                             HSBG1050
      DO 300 I=1,N                                                      HSBG1060
      K=K+IA                                                            HSBG1070
      LK=K+L1                                                           HSBG1080
      S=A(LK)                                                           HSBG1090
      LJ=L-IA                                                           HSBG1100
      DO 280 J=1,L2                                                     HSBG1110
      JK=K+J                                                            HSBG1120
      LJ=LJ+IA                                                          HSBG1130
  280 S=S+A(LJ)*A(JK)*1.0D0                                             HSBG1140
  300 A(LK)=S                                                           HSBG1150
C                                                                       HSBG1160
C        SET THE LOWER PART OF THE MATRIX TO ZERO                       HSBG1170
C                                                                       HSBG1180
      DO 310 I=L,LIA,IA                                                 HSBG1190
  310 A(I)=0.0                                                          HSBG1200
  320 L=L1                                                              HSBG1210
      GO TO 20                                                          HSBG1220
  360 RETURN                                                            HSBG1230
      END                                                               HSBG1240
