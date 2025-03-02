C                                                                       DMAT  10
C     ..................................................................DMAT  20
C                                                                       DMAT  30
C        SUBROUTINE DMATX                                               DMAT  40
C                                                                       DMAT  50
C        PURPOSE                                                        DMAT  60
C           COMPUTE MEANS OF VARIABLES IN EACH GROUP AND A POOLED       DMAT  70
C           DISPERSION MATRIX FOR ALL THE GROUPS. NORMALLY THIS SUB-    DMAT  80
C           ROUTINE IS USED IN THE PERFORMANCE OF DISCRIMINANT ANALYSIS.DMAT  90
C                                                                       DMAT 100
C        USAGE                                                          DMAT 110
C           CALL DMATX (K,M,N,X,XBAR,D,CMEAN)                           DMAT 120
C                                                                       DMAT 130
C        DESCRIPTION OF PARAMETERS                                      DMAT 140
C           K     - NUMBER OF GROUPS                                    DMAT 150
C           M     - NUMBER OF VARIABLES (MUST BE THE SAME FOR ALL       DMAT 160
C                   GROUPS).                                            DMAT 170
C           N     - INPUT VECTOR OF LENGTH K CONTAINING SAMPLE SIZES OF DMAT 180
C                   GROUPS.                                             DMAT 190
C           X     - INPUT VECTOR CONTAINING DATA IN THE MANNER EQUIVA-  DMAT 200
C                   LENT TO A 3-DIMENSIONAL FORTRAN ARRAY, X(1,1,1),    DMAT 210
C                   X(2,1,1), X(3,1,1), ETC.  THE FIRST SUBSCRIPT IS    DMAT 220
C                   CASE NUMBER, THE SECOND SUBSCRIPT IS VARIABLE NUMBERDMAT 230
C                   AND THE THIRD SUBSCRIPT IS GROUP NUMBER.  THE       DMAT 240
C                   LENGTH OF VECTOR X IS EQUAL TO THE TOTAL NUMBER OF  DMAT 250
C                   DATA POINTS, T*M, WHERE T = N(1)+N(2)+...+N(K).     DMAT 260
C           XBAR  - OUTPUT MATRIX (M X K) CONTAINING MEANS OF VARIABLES DMAT 270
C                   IN K GROUPS.                                        DMAT 280
C           D     - OUTPUT MATRIX (M X M) CONTAINING POOLED DISPERSION. DMAT 290
C           CMEAN - WORKING VECTOR OF LENGTH M.                         DMAT 300
C                                                                       DMAT 310
C        REMARKS                                                        DMAT 320
C           THE NUMBER OF VARIABLES MUST BE GREATER THAN OR EQUAL TO    DMAT 330
C           THE NUMBER OF GROUPS.                                       DMAT 340
C                                                                       DMAT 350
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DMAT 360
C           NONE                                                        DMAT 370
C                                                                       DMAT 380
C        METHOD                                                         DMAT 390
C           REFER TO 'BMD COMPUTER PROGRAMS MANUAL', EDITED BY W. J.    DMAT 400
C           DIXON, UCLA, 1964, AND T. W. ANDERSON, 'INTRODUCTION TO     DMAT 410
C           MULTIVARIATE STATISTICAL ANALYSIS', JOHN WILEY AND SONS,    DMAT 420
C           1958, SECTION 6.6-6.8.                                      DMAT 430
C                                                                       DMAT 440
C     ..................................................................DMAT 450
C                                                                       DMAT 460
      SUBROUTINE DMATX (K,M,N,X,XBAR,D,CMEAN)                           DMAT 470
      DIMENSION N(1),X(1),XBAR(1),D(1),CMEAN(1)                         DMAT 480
C                                                                       DMAT 490
C        ...............................................................DMAT 500
C                                                                       DMAT 510
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE  DMAT 520
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION      DMAT 530
C        STATEMENT WHICH FOLLOWS.                                       DMAT 540
C                                                                       DMAT 550
C     DOUBLE PRECISION XBAR,D,CMEAN                                     DMAT 560
C                                                                       DMAT 570
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS    DMAT 580
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS      DMAT 590
C        ROUTINE.                                                       DMAT 600
C                                                                       DMAT 610
C        ...............................................................DMAT 620
C                                                                       DMAT 630
C     INITIALIZATION                                                    DMAT 640
C                                                                       DMAT 650
      MM=M*M                                                            DMAT 660
      DO 100 I=1,MM                                                     DMAT 670
  100 D(I)=0.0                                                          DMAT 680
C                                                                       DMAT 690
C     CALCULATE MEANS                                                   DMAT 700
C                                                                       DMAT 710
      N4=0                                                              DMAT 720
      L=0                                                               DMAT 730
      LM=0                                                              DMAT 740
      DO 160 NG=1,K                                                     DMAT 750
      N1=N(NG)                                                          DMAT 760
      FN=N1                                                             DMAT 770
      DO 130 J=1,M                                                      DMAT 780
      LM=LM+1                                                           DMAT 790
      XBAR(LM)=0.0                                                      DMAT 800
      DO 120 I=1,N1                                                     DMAT 810
      L=L+1                                                             DMAT 820
  120 XBAR(LM)=XBAR(LM)+X(L)                                            DMAT 830
  130 XBAR(LM)=XBAR(LM)/FN                                              DMAT 840
C                                                                       DMAT 850
C     CALCULATE SUMS OF CROSS-PRODUCTS OF DEVIATIONS                    DMAT 860
C                                                                       DMAT 870
      LMEAN=LM-M                                                        DMAT 880
      DO 150 I=1,N1                                                     DMAT 890
      LL=N4+I-N1                                                        DMAT 900
      DO 140 J=1,M                                                      DMAT 910
      LL=LL+N1                                                          DMAT 920
      N2=LMEAN+J                                                        DMAT 930
  140 CMEAN(J)=X(LL)-XBAR(N2)                                           DMAT 940
      LL=0                                                              DMAT 950
      DO 150 J=1,M                                                      DMAT 960
      DO 150 JJ=1,M                                                     DMAT 970
      LL=LL+1                                                           DMAT 980
  150 D(LL)=D(LL)+CMEAN(J)*CMEAN(JJ)                                    DMAT 990
  160 N4=N4+N1*M                                                        DMAT1000
C                                                                       DMAT1010
C     CALCULATE THE POOLED DISPERSION MATRIX                            DMAT1020
C                                                                       DMAT1030
      LL=-K                                                             DMAT1040
      DO 170 I=1,K                                                      DMAT1050
  170 LL=LL+N(I)                                                        DMAT1060
      FN=LL                                                             DMAT1070
      DO 180 I=1,MM                                                     DMAT1080
  180 D(I)=D(I)/FN                                                      DMAT1090
C                                                                       DMAT1100
      RETURN                                                            DMAT1110
      END                                                               DMAT1120
