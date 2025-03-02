C                                                                       MFGR  10
C     ..................................................................MFGR  20
C                                                                       MFGR  30
C        SUBROUTINE MFGR                                                MFGR  40
C                                                                       MFGR  50
C        PURPOSE                                                        MFGR  60
C           FOR A GIVEN M BY N MATRIX THE FOLLOWING CALCULATIONS        MFGR  70
C           ARE PERFORMED                                               MFGR  80
C           (1) DETERMINE RANK AND LINEARLY INDEPENDENT ROWS AND        MFGR  90
C               COLUMNS (BASIS).                                        MFGR 100
C           (2) FACTORIZE A SUBMATRIX OF MAXIMAL RANK.                  MFGR 110
C           (3) EXPRESS NON-BASIC ROWS IN TERMS OF BASIC ONES.          MFGR 120
C           (4) EXPRESS BASIC VARIABLES IN TERMS OF FREE ONES.          MFGR 130
C                                                                       MFGR 140
C        USAGE                                                          MFGR 150
C           CALL MFGR(A,M,N,EPS,IRANK,IROW,ICOL)                        MFGR 160
C                                                                       MFGR 170
C        DESCRIPTION OF PARAMETERS                                      MFGR 180
C           A      - GIVEN MATRIX WITH M ROWS AND N COLUMNS.            MFGR 190
C                    ON RETURN A CONTAINS THE FIVE SUBMATRICES          MFGR 200
C                    L, R, H, D, O.                                     MFGR 210
C           M      - NUMBER OF ROWS OF MATRIX A.                        MFGR 220
C           N      - NUMBER OF COLUMNS OF MATRIX A.                     MFGR 230
C           EPS    - TESTVALUE FOR ZERO AFFECTED BY ROUNDOFF NOISE.     MFGR 240
C           IRANK  - RESULTANT RANK OF GIVEN MATRIX.                    MFGR 250
C           IROW   - INTEGER VECTOR OF DIMENSION M CONTAINING THE       MFGR 260
C                    SUBSCRIPTS OF BASIC ROWS IN IROW(1),...,IROW(IRANK)MFGR 270
C           ICOL   - INTEGER VECTOR OF DIMENSION N CONTAINING THE       MFGR 280
C                    SUBSCRIPTS OF BASIC COLUMNS IN ICOL(1) UP TO       MFGR 290
C                    ICOL(IRANK).                                       MFGR 300
C                                                                       MFGR 310
C        REMARKS                                                        MFGR 320
C           THE LEFT HAND TRIANGULAR FACTOR IS NORMALIZED SUCH THAT     MFGR 330
C           THE DIAGONAL CONTAINS ALL ONES THUS ALLOWING TO STORE ONLY  MFGR 340
C           THE SUBDIAGONAL PART.                                       MFGR 350
C                                                                       MFGR 360
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  MFGR 370
C           NONE                                                        MFGR 380
C                                                                       MFGR 390
C        METHOD                                                         MFGR 400
C           GAUSSIAN ELIMINATION TECHNIQUE IS USED FOR CALCULATION      MFGR 410
C           OF THE TRIANGULAR FACTORS OF A GIVEN MATRIX.                MFGR 420
C           COMPLETE PIVOTING IS BUILT IN.                              MFGR 430
C           IN CASE OF A SINGULAR MATRIX ONLY THE TRIANGULAR FACTORS    MFGR 440
C           OF A SUBMATRIX OF MAXIMAL RANK ARE RETAINED.                MFGR 450
C           THE REMAINING PARTS OF THE RESULTANT MATRIX GIVE THE        MFGR 460
C           DEPENDENCIES OF ROWS AND THE SOLUTION OF THE HOMOGENEOUS    MFGR 470
C           MATRIX EQUATION A*X=0.                                      MFGR 480
C                                                                       MFGR 490
C     ..................................................................MFGR 500
C                                                                       MFGR 510
      SUBROUTINE MFGR(A,M,N,EPS,IRANK,IROW,ICOL)                        MFGR 520
C                                                                       MFGR 530
C        DIMENSIONED DUMMY VARIABLES                                    MFGR 540
      DIMENSION A(1),IROW(1),ICOL(1)                                    MFGR 550
C                                                                       MFGR 560
C       TEST OF SPECIFIED DIMENSIONS                                    MFGR 570
      IF(M)2,2,1                                                        MFGR 580
    1 IF(N)2,2,4                                                        MFGR 590
    2 IRANK=-1                                                          MFGR 600
    3 RETURN                                                            MFGR 610
C       RETURN IN CASE OF FORMAL ERRORS                                 MFGR 620
C                                                                       MFGR 630
C                                                                       MFGR 640
C        INITIALIZE COLUMN INDEX VECTOR                                 MFGR 650
C        SEARCH FIRST PIVOT ELEMENT                                     MFGR 660
    4 IRANK=0                                                           MFGR 670
      PIV=0.                                                            MFGR 680
      JJ=0                                                              MFGR 690
      DO 6 J=1,N                                                        MFGR 700
      ICOL(J)=J                                                         MFGR 710
      DO 6 I=1,M                                                        MFGR 720
      JJ=JJ+1                                                           MFGR 730
      HOLD=A(JJ)                                                        MFGR 740
      IF(ABS(PIV)-ABS(HOLD))5,6,6                                       MFGR 750
    5 PIV=HOLD                                                          MFGR 760
      IR=I                                                              MFGR 770
      IC=J                                                              MFGR 780
    6 CONTINUE                                                          MFGR 790
C                                                                       MFGR 800
C        INITIALIZE ROW INDEX VECTOR                                    MFGR 810
      DO 7 I=1,M                                                        MFGR 820
    7 IROW(I)=I                                                         MFGR 830
C                                                                       MFGR 840
C        SET UP INTERNAL TOLERANCE                                      MFGR 850
      TOL=ABS(EPS*PIV)                                                  MFGR 860
C                                                                       MFGR 870
C        INITIALIZE ELIMINATION LOOP                                    MFGR 880
      NM=N*M                                                            MFGR 890
      DO 19 NCOL=M,NM,M                                                 MFGR 900
C                                                                       MFGR 910
C        TEST FOR FEASIBILITY OF PIVOT ELEMENT                          MFGR 920
    8 IF(ABS(PIV)-TOL)20,20,9                                           MFGR 930
C                                                                       MFGR 940
C        UPDATE RANK                                                    MFGR 950
    9 IRANK=IRANK+1                                                     MFGR 960
C                                                                       MFGR 970
C        INTERCHANGE ROWS IF NECESSARY                                  MFGR 980
      JJ=IR-IRANK                                                       MFGR 990
      IF(JJ)12,12,10                                                    MFGR1000
   10 DO 11 J=IRANK,NM,M                                                MFGR1010
      I=J+JJ                                                            MFGR1020
      SAVE=A(J)                                                         MFGR1030
      A(J)=A(I)                                                         MFGR1040
   11 A(I)=SAVE                                                         MFGR1050
C                                                                       MFGR1060
C        UPDATE ROW INDEX VECTOR                                        MFGR1070
      JJ=IROW(IR)                                                       MFGR1080
      IROW(IR)=IROW(IRANK)                                              MFGR1090
      IROW(IRANK)=JJ                                                    MFGR1100
C                                                                       MFGR1110
C        INTERCHANGE COLUMNS IF NECESSARY                               MFGR1120
   12 JJ=(IC-IRANK)*M                                                   MFGR1130
      IF(JJ)15,15,13                                                    MFGR1140
   13 KK=NCOL                                                           MFGR1150
      DO 14 J=1,M                                                       MFGR1160
      I=KK+JJ                                                           MFGR1170
      SAVE=A(KK)                                                        MFGR1180
      A(KK)=A(I)                                                        MFGR1190
      KK=KK-1                                                           MFGR1200
   14 A(I)=SAVE                                                         MFGR1210
C                                                                       MFGR1220
C        UPDATE COLUMN INDEX VECTOR                                     MFGR1230
      JJ=ICOL(IC)                                                       MFGR1240
      ICOL(IC)=ICOL(IRANK)                                              MFGR1250
      ICOL(IRANK)=JJ                                                    MFGR1260
   15 KK=IRANK+1                                                        MFGR1270
      MM=IRANK-M                                                        MFGR1280
      LL=NCOL+MM                                                        MFGR1290
C                                                                       MFGR1300
C        TEST FOR LAST ROW                                              MFGR1310
      IF(MM)16,25,25                                                    MFGR1320
C                                                                       MFGR1330
C        TRANSFORM CURRENT SUBMATRIX AND SEARCH NEXT PIVOT              MFGR1340
   16 JJ=LL                                                             MFGR1350
      SAVE=PIV                                                          MFGR1360
      PIV=0.                                                            MFGR1370
      DO 19 J=KK,M                                                      MFGR1380
      JJ=JJ+1                                                           MFGR1390
      HOLD=A(JJ)/SAVE                                                   MFGR1400
      A(JJ)=HOLD                                                        MFGR1410
      L=J-IRANK                                                         MFGR1420
C                                                                       MFGR1430
C        TEST FOR LAST COLUMN                                           MFGR1440
      IF(IRANK-N)17,19,19                                               MFGR1450
   17 II=JJ                                                             MFGR1460
      DO 19 I=KK,N                                                      MFGR1470
      II=II+M                                                           MFGR1480
      MM=II-L                                                           MFGR1490
      A(II)=A(II)-HOLD*A(MM)                                            MFGR1500
      IF(ABS(A(II))-ABS(PIV))19,19,18                                   MFGR1510
   18 PIV=A(II)                                                         MFGR1520
      IR=J                                                              MFGR1530
      IC=I                                                              MFGR1540
   19 CONTINUE                                                          MFGR1550
C                                                                       MFGR1560
C        SET UP MATRIX EXPRESSING ROW DEPENDENCIES                      MFGR1570
   20 IF(IRANK-1)3,25,21                                                MFGR1580
   21 IR=LL                                                             MFGR1590
      DO 24 J=2,IRANK                                                   MFGR1600
      II=J-1                                                            MFGR1610
      IR=IR-M                                                           MFGR1620
      JJ=LL                                                             MFGR1630
      DO 23 I=KK,M                                                      MFGR1640
      HOLD=0.                                                           MFGR1650
      JJ=JJ+1                                                           MFGR1660
      MM=JJ                                                             MFGR1670
      IC=IR                                                             MFGR1680
      DO 22 L=1,II                                                      MFGR1690
      HOLD=HOLD+A(MM)*A(IC)                                             MFGR1700
      IC=IC-1                                                           MFGR1710
   22 MM=MM-M                                                           MFGR1720
   23 A(MM)=A(MM)-HOLD                                                  MFGR1730
   24 CONTINUE                                                          MFGR1740
C                                                                       MFGR1750
C        TEST FOR COLUMN REGULARITY                                     MFGR1760
   25 IF(N-IRANK)3,3,26                                                 MFGR1770
C                                                                       MFGR1780
C        SET UP MATRIX EXPRESSING BASIC VARIABLES IN TERMS OF FREE      MFGR1790
C       PARAMETERS (HOMOGENEOUS SOLUTION).                              MFGR1800
   26 IR=LL                                                             MFGR1810
      KK=LL+M                                                           MFGR1820
      DO 30 J=1,IRANK                                                   MFGR1830
      DO 29 I=KK,NM,M                                                   MFGR1840
      JJ=IR                                                             MFGR1850
      LL=I                                                              MFGR1860
      HOLD=0.                                                           MFGR1870
      II=J                                                              MFGR1880
   27 II=II-1                                                           MFGR1890
      IF(II)29,29,28                                                    MFGR1900
   28 HOLD=HOLD-A(JJ)*A(LL)                                             MFGR1910
      JJ=JJ-M                                                           MFGR1920
      LL=LL-1                                                           MFGR1930
      GOTO 27                                                           MFGR1940
   29 A(LL)=(HOLD-A(LL))/A(JJ)                                          MFGR1950
   30 IR=IR-1                                                           MFGR1960
      RETURN                                                            MFGR1970
      END                                                               MFGR1980
