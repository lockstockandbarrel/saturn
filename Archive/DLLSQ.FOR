C                                                                       DLLS  10
C     ..................................................................DLLS  20
C                                                                       DLLS  30
C        SUBROUTINE DLLSQ                                               DLLS  40
C                                                                       DLLS  50
C        PURPOSE                                                        DLLS  60
C           TO SOLVE LINEAR LEAST SQUARES PROBLEMS, I.E. TO MINIMIZE    DLLS  70
C           THE EUCLIDEAN NORM OF B-A*X, WHERE A IS A M BY N MATRIX     DLLS  80
C           WITH M NOT LESS THAN N. IN THE SPECIAL CASE M=N SYSTEMS OF  DLLS  90
C           LINEAR EQUATIONS MAY BE SOLVED.                             DLLS 100
C                                                                       DLLS 110
C        USAGE                                                          DLLS 120
C           CALL DLLSQ (A,B,M,N,L,X,IPIV,EPS,IER,AUX)                   DLLS 130
C                                                                       DLLS 140
C        DESCRIPTION OF PARAMETERS                                      DLLS 150
C           A      - DOUBLE PRECISION M BY N COEFFICIENT MATRIX         DLLS 160
C                    (DESTROYED).                                       DLLS 170
C           B      - DOUBLE PRECISION M BY L RIGHT HAND SIDE MATRIX     DLLS 180
C                    (DESTROYED).                                       DLLS 190
C           M      - ROW NUMBER OF MATRICES A AND B.                    DLLS 200
C           N      - COLUMN NUMBER OF MATRIX A, ROW NUMBER OF MATRIX X. DLLS 210
C           L      - COLUMN NUMBER OF MATRICES B AND X.                 DLLS 220
C           X      - DOUBLE PRECISION N BY L SOLUTION MATRIX.           DLLS 230
C           IPIV   - INTEGER OUTPUT VECTOR OF DIMENSION N WHICH         DLLS 240
C                    CONTAINS INFORMATIONS ON COLUMN INTERCHANGES       DLLS 250
C                    IN MATRIX A. (SEE REMARK NO.3).                    DLLS 260
C           EPS    - SINGLE PRECISION INPUT PARAMETER WHICH SPECIFIES   DLLS 270
C                    A RELATIVE TOLERANCE FOR DETERMINATION OF RANK OF  DLLS 280
C                    MATRIX A.                                          DLLS 290
C           IER    - A RESULTING ERROR PARAMETER.                       DLLS 300
C           AUX    - A DOUBLE PRECISION AUXILIARY STORAGE ARRAY OF      DLLS 310
C                    DIMENSION MAX(2*N,L). ON RETURN FIRST L LOCATIONS  DLLS 320
C                    OF AUX CONTAIN THE RESULTING LEAST SQUARES.        DLLS 330
C                                                                       DLLS 340
C        REMARKS                                                        DLLS 350
C           (1) NO ACTION BESIDES ERROR MESSAGE IER=-2 IN CASE          DLLS 360
C               M LESS THAN N.                                          DLLS 370
C           (2) NO ACTION BESIDES ERROR MESSAGE IER=-1 IN CASE          DLLS 380
C               OF A ZERO-MATRIX A.                                     DLLS 390
C           (3) IF RANK K OF MATRIX A IS FOUND TO BE LESS THAN N BUT    DLLS 400
C               GREATER THAN 0, THE PROCEDURE RETURNS WITH ERROR CODE   DLLS 410
C               IER=K INTO CALLING PROGRAM. THE LAST N-K ELEMENTS OF    DLLS 420
C               VECTOR IPIV DENOTE THE USELESS COLUMNS IN MATRIX A.     DLLS 430
C               THE REMAINING USEFUL COLUMNS FORM A BASE OF MATRIX A.   DLLS 440
C           (4) IF THE PROCEDURE WAS SUCCESSFUL, ERROR PARAMETER IER    DLLS 450
C               IS SET TO 0.                                            DLLS 460
C                                                                       DLLS 470
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DLLS 480
C           NONE                                                        DLLS 490
C                                                                       DLLS 500
C        METHOD                                                         DLLS 510
C           HOUSEHOLDER TRANSFORMATIONS ARE USED TO TRANSFORM MATRIX A  DLLS 520
C           TO UPPER TRIANGULAR FORM. AFTER HAVING APPLIED THE SAME     DLLS 530
C           TRANSFORMATION TO THE RIGHT HAND SIDE MATRIX B, AN          DLLS 540
C           APPROXIMATE SOLUTION OF THE PROBLEM IS COMPUTED BY          DLLS 550
C           BACK SUBSTITUTION. FOR REFERENCE, SEE                       DLLS 560
C           G. GOLUB, NUMERICAL METHODS FOR SOLVING LINEAR LEAST        DLLS 570
C           SQUARES PROBLEMS, NUMERISCHE MATHEMATIK, VOL.7,             DLLS 580
C           ISS.3 (1965), PP.206-216.                                   DLLS 590
C                                                                       DLLS 600
C     ..................................................................DLLS 610
C                                                                       DLLS 620
      SUBROUTINE DLLSQ(A,B,M,N,L,X,IPIV,EPS,IER,AUX)                    DLLS 630
C                                                                       DLLS 640
      DIMENSION A(1),B(1),X(1),IPIV(1),AUX(1)                           DLLS 650
      DOUBLE PRECISION A,B,X,AUX,PIV,H,SIG,BETA,TOL                     DLLS 660
C                                                                       DLLS 670
C     ERROR TEST                                                        DLLS 680
      IF(M-N)30,1,1                                                     DLLS 690
C                                                                       DLLS 700
C     GENERATION OF INITIAL VECTOR S(K) (K=1,2,...,N) IN STORAGE        DLLS 710
C     LOCATIONS AUX(K) (K=1,2,...,N)                                    DLLS 720
    1 PIV=0.D0                                                          DLLS 730
      IEND=0                                                            DLLS 740
      DO 4 K=1,N                                                        DLLS 750
      IPIV(K)=K                                                         DLLS 760
      H=0.D0                                                            DLLS 770
      IST=IEND+1                                                        DLLS 780
      IEND=IEND+M                                                       DLLS 790
      DO 2 I=IST,IEND                                                   DLLS 800
    2 H=H+A(I)*A(I)                                                     DLLS 810
      AUX(K)=H                                                          DLLS 820
      IF(H-PIV)4,4,3                                                    DLLS 830
    3 PIV=H                                                             DLLS 840
      KPIV=K                                                            DLLS 850
    4 CONTINUE                                                          DLLS 860
C                                                                       DLLS 870
C     ERROR TEST                                                        DLLS 880
      IF(PIV)31,31,5                                                    DLLS 890
C                                                                       DLLS 900
C     DEFINE TOLERANCE FOR CHECKING RANK OF A                           DLLS 910
    5 SIG=DSQRT(PIV)                                                    DLLS 920
      TOL=SIG*ABS(EPS)                                                  DLLS 930
C                                                                       DLLS 940
C                                                                       DLLS 950
C     DECOMPOSITION LOOP                                                DLLS 960
      LM=L*M                                                            DLLS 970
      IST=-M                                                            DLLS 980
      DO 21 K=1,N                                                       DLLS 990
      IST=IST+M+1                                                       DLLS1000
      IEND=IST+M-K                                                      DLLS1010
      I=KPIV-K                                                          DLLS1020
      IF(I)8,8,6                                                        DLLS1030
C                                                                       DLLS1040
C     INTERCHANGE K-TH COLUMN OF A WITH KPIV-TH IN CASE KPIV.GT.K       DLLS1050
    6 H=AUX(K)                                                          DLLS1060
      AUX(K)=AUX(KPIV)                                                  DLLS1070
      AUX(KPIV)=H                                                       DLLS1080
      ID=I*M                                                            DLLS1090
      DO 7 I=IST,IEND                                                   DLLS1100
      J=I+ID                                                            DLLS1110
      H=A(I)                                                            DLLS1120
      A(I)=A(J)                                                         DLLS1130
    7 A(J)=H                                                            DLLS1140
C                                                                       DLLS1150
C     COMPUTATION OF PARAMETER SIG                                      DLLS1160
    8 IF(K-1)11,11,9                                                    DLLS1170
    9 SIG=0.D0                                                          DLLS1180
      DO 10 I=IST,IEND                                                  DLLS1190
   10 SIG=SIG+A(I)*A(I)                                                 DLLS1200
      SIG=DSQRT(SIG)                                                    DLLS1210
C                                                                       DLLS1220
C     TEST ON SINGULARITY                                               DLLS1230
      IF(SIG-TOL)32,32,11                                               DLLS1240
C                                                                       DLLS1250
C     GENERATE CORRECT SIGN OF PARAMETER SIG                            DLLS1260
   11 H=A(IST)                                                          DLLS1270
      IF(H)12,13,13                                                     DLLS1280
   12 SIG=-SIG                                                          DLLS1290
C                                                                       DLLS1300
C     SAVE INTERCHANGE INFORMATION                                      DLLS1310
   13 IPIV(KPIV)=IPIV(K)                                                DLLS1320
      IPIV(K)=KPIV                                                      DLLS1330
C                                                                       DLLS1340
C     GENERATION OF VECTOR UK IN K-TH COLUMN OF MATRIX A AND OF         DLLS1350
C     PARAMETER BETA                                                    DLLS1360
      BETA=H+SIG                                                        DLLS1370
      A(IST)=BETA                                                       DLLS1380
      BETA=1.D0/(SIG*BETA)                                              DLLS1390
      J=N+K                                                             DLLS1400
      AUX(J)=-SIG                                                       DLLS1410
      IF(K-N)14,19,19                                                   DLLS1420
C                                                                       DLLS1430
C     TRANSFORMATION OF MATRIX A                                        DLLS1440
   14 PIV=0.D0                                                          DLLS1450
      ID=0                                                              DLLS1460
      JST=K+1                                                           DLLS1470
      KPIV=JST                                                          DLLS1480
      DO 18 J=JST,N                                                     DLLS1490
      ID=ID+M                                                           DLLS1500
      H=0.D0                                                            DLLS1510
      DO 15 I=IST,IEND                                                  DLLS1520
      II=I+ID                                                           DLLS1530
   15 H=H+A(I)*A(II)                                                    DLLS1540
      H=BETA*H                                                          DLLS1550
      DO 16 I=IST,IEND                                                  DLLS1560
      II=I+ID                                                           DLLS1570
   16 A(II)=A(II)-A(I)*H                                                DLLS1580
C                                                                       DLLS1590
C     UPDATING OF ELEMENT S(J) STORED IN LOCATION AUX(J)                DLLS1600
      II=IST+ID                                                         DLLS1610
      H=AUX(J)-A(II)*A(II)                                              DLLS1620
      AUX(J)=H                                                          DLLS1630
      IF(H-PIV)18,18,17                                                 DLLS1640
   17 PIV=H                                                             DLLS1650
      KPIV=J                                                            DLLS1660
   18 CONTINUE                                                          DLLS1670
C                                                                       DLLS1680
C     TRANSFORMATION OF RIGHT HAND SIDE MATRIX B                        DLLS1690
   19 DO 21 J=K,LM,M                                                    DLLS1700
      H=0.D0                                                            DLLS1710
      IEND=J+M-K                                                        DLLS1720
      II=IST                                                            DLLS1730
      DO 20 I=J,IEND                                                    DLLS1740
      H=H+A(II)*B(I)                                                    DLLS1750
   20 II=II+1                                                           DLLS1760
      H=BETA*H                                                          DLLS1770
      II=IST                                                            DLLS1780
      DO 21 I=J,IEND                                                    DLLS1790
      B(I)=B(I)-A(II)*H                                                 DLLS1800
   21 II=II+1                                                           DLLS1810
C     END OF DECOMPOSITION LOOP                                         DLLS1820
C                                                                       DLLS1830
C                                                                       DLLS1840
C     BACK SUBSTITUTION AND BACK INTERCHANGE                            DLLS1850
      IER=0                                                             DLLS1860
      I=N                                                               DLLS1870
      LN=L*N                                                            DLLS1880
      PIV=1.D0/AUX(2*N)                                                 DLLS1890
      DO 22 K=N,LN,N                                                    DLLS1900
      X(K)=PIV*B(I)                                                     DLLS1910
   22 I=I+M                                                             DLLS1920
      IF(N-1)26,26,23                                                   DLLS1930
   23 JST=(N-1)*M+N                                                     DLLS1940
      DO 25 J=2,N                                                       DLLS1950
      JST=JST-M-1                                                       DLLS1960
      K=N+N+1-J                                                         DLLS1970
      PIV=1.D0/AUX(K)                                                   DLLS1980
      KST=K-N                                                           DLLS1990
      ID=IPIV(KST)-KST                                                  DLLS2000
      IST=2-J                                                           DLLS2010
      DO 25 K=1,L                                                       DLLS2020
      H=B(KST)                                                          DLLS2030
      IST=IST+N                                                         DLLS2040
      IEND=IST+J-2                                                      DLLS2050
      II=JST                                                            DLLS2060
      DO 24 I=IST,IEND                                                  DLLS2070
      II=II+M                                                           DLLS2080
   24 H=H-A(II)*X(I)                                                    DLLS2090
      I=IST-1                                                           DLLS2100
      II=I+ID                                                           DLLS2110
      X(I)=X(II)                                                        DLLS2120
      X(II)=PIV*H                                                       DLLS2130
   25 KST=KST+M                                                         DLLS2140
C                                                                       DLLS2150
C                                                                       DLLS2160
C     COMPUTATION OF LEAST SQUARES                                      DLLS2170
   26 IST=N+1                                                           DLLS2180
      IEND=0                                                            DLLS2190
      DO 29 J=1,L                                                       DLLS2200
      IEND=IEND+M                                                       DLLS2210
      H=0.D0                                                            DLLS2220
      IF(M-N)29,29,27                                                   DLLS2230
   27 DO 28 I=IST,IEND                                                  DLLS2240
   28 H=H+B(I)*B(I)                                                     DLLS2250
      IST=IST+M                                                         DLLS2260
   29 AUX(J)=H                                                          DLLS2270
      RETURN                                                            DLLS2280
C                                                                       DLLS2290
C     ERROR RETURN IN CASE M LESS THAN N                                DLLS2300
   30 IER=-2                                                            DLLS2310
      RETURN                                                            DLLS2320
C                                                                       DLLS2330
C     ERROR RETURN IN CASE OF ZERO-MATRIX A                             DLLS2340
   31 IER=-1                                                            DLLS2350
      RETURN                                                            DLLS2360
C                                                                       DLLS2370
C     ERROR RETURN IN CASE OF RANK OF MATRIX A LESS THAN N              DLLS2380
   32 IER=K-1                                                           DLLS2390
      RETURN                                                            DLLS2400
      END                                                               DLLS2410
