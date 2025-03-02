C                                                                       LLSQ  10
C     ..................................................................LLSQ  20
C                                                                       LLSQ  30
C        SUBROUTINE LLSQ                                                LLSQ  40
C                                                                       LLSQ  50
C        PURPOSE                                                        LLSQ  60
C           TO SOLVE LINEAR LEAST SQUARES PROBLEMS, I.E. TO MINIMIZE    LLSQ  70
C           THE EUCLIDEAN NORM OF B-A*X, WHERE A IS A M BY N MATRIX     LLSQ  80
C           WITH M NOT LESS THAN N. IN THE SPECIAL CASE M=N SYSTEMS OF  LLSQ  90
C           LINEAR EQUATIONS MAY BE SOLVED.                             LLSQ 100
C                                                                       LLSQ 110
C        USAGE                                                          LLSQ 120
C           CALL LLSQ (A,B,M,N,L,X,IPIV,EPS,IER,AUX)                    LLSQ 130
C                                                                       LLSQ 140
C        DESCRIPTION OF PARAMETERS                                      LLSQ 150
C           A      - M BY N COEFFICIENT MATRIX (DESTROYED).             LLSQ 160
C           B      - M BY L RIGHT HAND SIDE MATRIX (DESTROYED).         LLSQ 170
C           M      - ROW NUMBER OF MATRICES A AND B.                    LLSQ 180
C           N      - COLUMN NUMBER OF MATRIX A, ROW NUMBER OF MATRIX X. LLSQ 190
C           L      - COLUMN NUMBER OF MATRICES B AND X.                 LLSQ 200
C           X      - N BY L SOLUTION MATRIX.                            LLSQ 210
C           IPIV   - INTEGER OUTPUT VECTOR OF DIMENSION N WHICH         LLSQ 220
C                    CONTAINS INFORMATIONS ON COLUMN INTERCHANGES       LLSQ 230
C                    IN MATRIX A. (SEE REMARK NO.3).                    LLSQ 240
C           EPS    - INPUT PARAMETER WHICH SPECIFIES A RELATIVE         LLSQ 250
C                    TOLERANCE FOR DETERMINATION OF RANK OF MATRIX A.   LLSQ 260
C           IER    - A RESULTING ERROR PARAMETER.                       LLSQ 270
C           AUX    - AUXILIARY STORAGE ARRAY OF DIMENSION MAX(2*N,L).   LLSQ 280
C                    ON RETURN FIRST L LOCATIONS OF AUX CONTAIN THE     LLSQ 290
C                    RESULTING LEAST SQUARES.                           LLSQ 300
C                                                                       LLSQ 310
C        REMARKS                                                        LLSQ 320
C           (1) NO ACTION BESIDES ERROR MESSAGE IER=-2 IN CASE          LLSQ 330
C               M LESS THAN N.                                          LLSQ 340
C           (2) NO ACTION BESIDES ERROR MESSAGE IER=-1 IN CASE          LLSQ 350
C               OF A ZERO-MATRIX A.                                     LLSQ 360
C           (3) IF RANK K OF MATRIX A IS FOUND TO BE LESS THAN N BUT    LLSQ 370
C               GREATER THAN 0, THE PROCEDURE RETURNS WITH ERROR CODE   LLSQ 380
C               IER=K INTO CALLING PROGRAM. THE LAST N-K ELEMENTS OF    LLSQ 390
C               VECTOR IPIV DENOTE THE USELESS COLUMNS IN MATRIX A.     LLSQ 400
C               THE REMAINING USEFUL COLUMNS FORM A BASE OF MATRIX A.   LLSQ 410
C           (4) IF THE PROCEDURE WAS SUCCESSFUL, ERROR PARAMETER IER    LLSQ 420
C               IS SET TO 0.                                            LLSQ 430
C                                                                       LLSQ 440
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  LLSQ 450
C           NONE                                                        LLSQ 460
C                                                                       LLSQ 470
C        METHOD                                                         LLSQ 480
C           HOUSEHOLDER TRANSFORMATIONS ARE USED TO TRANSFORM MATRIX A  LLSQ 490
C           TO UPPER TRIANGULAR FORM. AFTER HAVING APPLIED THE SAME     LLSQ 500
C           TRANSFORMATION TO THE RIGHT HAND SIDE MATRIX B, AN          LLSQ 510
C           APPROXIMATE SOLUTION OF THE PROBLEM IS COMPUTED BY          LLSQ 520
C           BACK SUBSTITUTION. FOR REFERENCE, SEE                       LLSQ 530
C           G. GOLUB, NUMERICAL METHODS FOR SOLVING LINEAR LEAST        LLSQ 540
C           SQUARES PROBLEMS, NUMERISCHE MATHEMATIK, VOL.7,             LLSQ 550
C           ISS.3 (1965), PP.206-216.                                   LLSQ 560
C                                                                       LLSQ 570
C     ..................................................................LLSQ 580
C                                                                       LLSQ 590
      SUBROUTINE LLSQ(A,B,M,N,L,X,IPIV,EPS,IER,AUX)                     LLSQ 600
C                                                                       LLSQ 610
      DIMENSION A(1),B(1),X(1),IPIV(1),AUX(1)                           LLSQ 620
C                                                                       LLSQ 630
C     ERROR TEST                                                        LLSQ 640
      IF(M-N)30,1,1                                                     LLSQ 650
C                                                                       LLSQ 660
C     GENERATION OF INITIAL VECTOR S(K) (K=1,2,...,N) IN STORAGE        LLSQ 670
C     LOCATIONS AUX(K) (K=1,2,...,N)                                    LLSQ 680
    1 PIV=0.                                                            LLSQ 690
      IEND=0                                                            LLSQ 700
      DO 4 K=1,N                                                        LLSQ 710
      IPIV(K)=K                                                         LLSQ 720
      H=0.                                                              LLSQ 730
      IST=IEND+1                                                        LLSQ 740
      IEND=IEND+M                                                       LLSQ 750
      DO 2 I=IST,IEND                                                   LLSQ 760
    2 H=H+A(I)*A(I)                                                     LLSQ 770
      AUX(K)=H                                                          LLSQ 780
      IF(H-PIV)4,4,3                                                    LLSQ 790
    3 PIV=H                                                             LLSQ 800
      KPIV=K                                                            LLSQ 810
    4 CONTINUE                                                          LLSQ 820
C                                                                       LLSQ 830
C     ERROR TEST                                                        LLSQ 840
      IF(PIV)31,31,5                                                    LLSQ 850
C                                                                       LLSQ 860
C     DEFINE TOLERANCE FOR CHECKING RANK OF A                           LLSQ 870
    5 SIG=SQRT(PIV)                                                     LLSQ 880
      TOL=SIG*ABS(EPS)                                                  LLSQ 890
C                                                                       LLSQ 900
C                                                                       LLSQ 910
C     DECOMPOSITION LOOP                                                LLSQ 920
      LM=L*M                                                            LLSQ 930
      IST=-M                                                            LLSQ 940
      DO 21 K=1,N                                                       LLSQ 950
      IST=IST+M+1                                                       LLSQ 960
      IEND=IST+M-K                                                      LLSQ 970
      I=KPIV-K                                                          LLSQ 980
      IF(I)8,8,6                                                        LLSQ 990
C                                                                       LLSQ1000
C     INTERCHANGE K-TH COLUMN OF A WITH KPIV-TH IN CASE KPIV.GT.K       LLSQ1010
    6 H=AUX(K)                                                          LLSQ1020
      AUX(K)=AUX(KPIV)                                                  LLSQ1030
      AUX(KPIV)=H                                                       LLSQ1040
      ID=I*M                                                            LLSQ1050
      DO 7 I=IST,IEND                                                   LLSQ1060
      J=I+ID                                                            LLSQ1070
      H=A(I)                                                            LLSQ1080
      A(I)=A(J)                                                         LLSQ1090
    7 A(J)=H                                                            LLSQ1100
C                                                                       LLSQ1110
C     COMPUTATION OF PARAMETER SIG                                      LLSQ1120
    8 IF(K-1)11,11,9                                                    LLSQ1130
    9 SIG=0.                                                            LLSQ1140
      DO 10 I=IST,IEND                                                  LLSQ1150
   10 SIG=SIG+A(I)*A(I)                                                 LLSQ1160
      SIG=SQRT(SIG)                                                     LLSQ1170
C                                                                       LLSQ1180
C     TEST ON SINGULARITY                                               LLSQ1190
      IF(SIG-TOL)32,32,11                                               LLSQ1200
C                                                                       LLSQ1210
C     GENERATE CORRECT SIGN OF PARAMETER SIG                            LLSQ1220
   11 H=A(IST)                                                          LLSQ1230
      IF(H)12,13,13                                                     LLSQ1240
   12 SIG=-SIG                                                          LLSQ1250
C                                                                       LLSQ1260
C     SAVE INTERCHANGE INFORMATION                                      LLSQ1270
   13 IPIV(KPIV)=IPIV(K)                                                LLSQ1280
      IPIV(K)=KPIV                                                      LLSQ1290
C                                                                       LLSQ1300
C     GENERATION OF VECTOR UK IN K-TH COLUMN OF MATRIX A AND OF         LLSQ1310
C     PARAMETER BETA                                                    LLSQ1320
      BETA=H+SIG                                                        LLSQ1330
      A(IST)=BETA                                                       LLSQ1340
      BETA=1./(SIG*BETA)                                                LLSQ1350
      J=N+K                                                             LLSQ1360
      AUX(J)=-SIG                                                       LLSQ1370
      IF(K-N)14,19,19                                                   LLSQ1380
C                                                                       LLSQ1390
C     TRANSFORMATION OF MATRIX A                                        LLSQ1400
   14 PIV=0.                                                            LLSQ1410
      ID=0                                                              LLSQ1420
      JST=K+1                                                           LLSQ1430
      KPIV=JST                                                          LLSQ1440
      DO 18 J=JST,N                                                     LLSQ1450
      ID=ID+M                                                           LLSQ1460
      H=0.                                                              LLSQ1470
      DO 15 I=IST,IEND                                                  LLSQ1480
      II=I+ID                                                           LLSQ1490
   15 H=H+A(I)*A(II)                                                    LLSQ1500
      H=BETA*H                                                          LLSQ1510
      DO 16 I=IST,IEND                                                  LLSQ1520
      II=I+ID                                                           LLSQ1530
   16 A(II)=A(II)-A(I)*H                                                LLSQ1540
C                                                                       LLSQ1550
C     UPDATING OF ELEMENT S(J) STORED IN LOCATION AUX(J)                LLSQ1560
      II=IST+ID                                                         LLSQ1570
      H=AUX(J)-A(II)*A(II)                                              LLSQ1580
      AUX(J)=H                                                          LLSQ1590
      IF(H-PIV)18,18,17                                                 LLSQ1600
   17 PIV=H                                                             LLSQ1610
      KPIV=J                                                            LLSQ1620
   18 CONTINUE                                                          LLSQ1630
C                                                                       LLSQ1640
C     TRANSFORMATION OF RIGHT HAND SIDE MATRIX B                        LLSQ1650
   19 DO 21 J=K,LM,M                                                    LLSQ1660
      H=0.                                                              LLSQ1670
      IEND=J+M-K                                                        LLSQ1680
      II=IST                                                            LLSQ1690
      DO 20 I=J,IEND                                                    LLSQ1700
      H=H+A(II)*B(I)                                                    LLSQ1710
   20 II=II+1                                                           LLSQ1720
      H=BETA*H                                                          LLSQ1730
      II=IST                                                            LLSQ1740
      DO 21 I=J,IEND                                                    LLSQ1750
      B(I)=B(I)-A(II)*H                                                 LLSQ1760
   21 II=II+1                                                           LLSQ1770
C     END OF DECOMPOSITION LOOP                                         LLSQ1780
C                                                                       LLSQ1790
C                                                                       LLSQ1800
C     BACK SUBSTITUTION AND BACK INTERCHANGE                            LLSQ1810
      IER=0                                                             LLSQ1820
      I=N                                                               LLSQ1830
      LN=L*N                                                            LLSQ1840
      PIV=1./AUX(2*N)                                                   LLSQ1850
      DO 22 K=N,LN,N                                                    LLSQ1860
      X(K)=PIV*B(I)                                                     LLSQ1870
   22 I=I+M                                                             LLSQ1880
      IF(N-1)26,26,23                                                   LLSQ1890
   23 JST=(N-1)*M+N                                                     LLSQ1900
      DO 25 J=2,N                                                       LLSQ1910
      JST=JST-M-1                                                       LLSQ1920
      K=N+N+1-J                                                         LLSQ1930
      PIV=1./AUX(K)                                                     LLSQ1940
      KST=K-N                                                           LLSQ1950
      ID=IPIV(KST)-KST                                                  LLSQ1960
      IST=2-J                                                           LLSQ1970
      DO 25 K=1,L                                                       LLSQ1980
      H=B(KST)                                                          LLSQ1990
      IST=IST+N                                                         LLSQ2000
      IEND=IST+J-2                                                      LLSQ2010
      II=JST                                                            LLSQ2020
      DO 24 I=IST,IEND                                                  LLSQ2030
      II=II+M                                                           LLSQ2040
   24 H=H-A(II)*X(I)                                                    LLSQ2050
      I=IST-1                                                           LLSQ2060
      II=I+ID                                                           LLSQ2070
      X(I)=X(II)                                                        LLSQ2080
      X(II)=PIV*H                                                       LLSQ2090
   25 KST=KST+M                                                         LLSQ2100
C                                                                       LLSQ2110
C                                                                       LLSQ2120
C     COMPUTATION OF LEAST SQUARES                                      LLSQ2130
   26 IST=N+1                                                           LLSQ2140
      IEND=0                                                            LLSQ2150
      DO 29 J=1,L                                                       LLSQ2160
      IEND=IEND+M                                                       LLSQ2170
      H=0.                                                              LLSQ2180
      IF(M-N)29,29,27                                                   LLSQ2190
   27 DO 28 I=IST,IEND                                                  LLSQ2200
   28 H=H+B(I)*B(I)                                                     LLSQ2210
      IST=IST+M                                                         LLSQ2220
   29 AUX(J)=H                                                          LLSQ2230
      RETURN                                                            LLSQ2240
C                                                                       LLSQ2250
C     ERROR RETURN IN CASE M LESS THAN N                                LLSQ2260
   30 IER=-2                                                            LLSQ2270
      RETURN                                                            LLSQ2280
C                                                                       LLSQ2290
C     ERROR RETURN IN CASE OF ZERO-MATRIX A                             LLSQ2300
   31 IER=-1                                                            LLSQ2310
      RETURN                                                            LLSQ2320
C                                                                       LLSQ2330
C     ERROR RETURN IN CASE OF RANK OF MATRIX A LESS THAN N              LLSQ2340
   32 IER=K-1                                                           LLSQ2350
      RETURN                                                            LLSQ2360
      END                                                               LLSQ2370
