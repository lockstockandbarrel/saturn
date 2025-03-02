C                                                                       DMGR  10
C     ..................................................................DMGR  20
C                                                                       DMGR  30
C        SUBROUTINE DMFGR                                               DMGR  40
C                                                                       DMGR  50
C        PURPOSE                                                        DMGR  60
C           FOR A GIVEN M BY N MATRIX THE FOLLOWING CALCULATIONS        DMGR  70
C           ARE PERFORMED                                               DMGR  80
C           (1) DETERMINE RANK AND LINEARLY INDEPENDENT ROWS AND        DMGR  90
C               COLUMNS (BASIS).                                        DMGR 100
C           (2) FACTORIZE A SUBMATRIX OF MAXIMAL RANK.                  DMGR 110
C           (3) EXPRESS NON-BASIC ROWS IN TERMS OF BASIC ONES.          DMGR 120
C           (4) EXPRESS BASIC VARIABLES IN TERMS OF FREE ONES.          DMGR 130
C                                                                       DMGR 140
C        USAGE                                                          DMGR 150
C           CALL DMFGR(A,M,N,EPS,IRANK,IROW,ICOL)                       DMGR 160
C                                                                       DMGR 170
C        DESCRIPTION OF PARAMETERS                                      DMGR 180
C           A      - DOUBLE PRECISION GIVEN MATRIX WITH M ROWS          DMGR 190
C                    AND N COLUMNS.                                     DMGR 200
C                    ON RETURN A CONTAINS THE TRIANGULAR FACTORS        DMGR 210
C                    OF A SUBMATRIX OF MAXIMAL RANK.                    DMGR 220
C           M      - NUMBER OF ROWS OF MATRIX A.                        DMGR 230
C           N      - NUMBER OF COLUMNS OF MATRIX A.                     DMGR 240
C           EPS    - SINGLE PRECISION TESTVALUE FOR ZERO AFFECTED BY    DMGR 250
C                    ROUNDOFF NOISE.                                    DMGR 260
C           IRANK  - RESULTANT RANK OF GIVEN MATRIX.                    DMGR 270
C           IROW   - INTEGER VECTOR OF DIMENSION M CONTAINING THE       DMGR 280
C                    SUBSCRIPTS OF BASIC ROWS IN IROW(1),...,IROW(IRANK)DMGR 290
C           ICOL   - INTEGER VECTOR OF DIMENSION N CONTAINING THE       DMGR 300
C                    SUBSCRIPTS OF BASIC COLUMNS IN ICOL(1) UP TO       DMGR 310
C                    ICOL(IRANK).                                       DMGR 320
C                                                                       DMGR 330
C        REMARKS                                                        DMGR 340
C           THE LEFT HAND TRIANGULAR FACTOR IS NORMALIZED SUCH THAT     DMGR 350
C           THE DIAGONAL CONTAINS ALL ONES THUS ALLOWING TO STORE ONLY  DMGR 360
C           THE SUBDIAGONAL PART.                                       DMGR 370
C                                                                       DMGR 380
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DMGR 390
C           NONE                                                        DMGR 400
C                                                                       DMGR 410
C        METHOD                                                         DMGR 420
C           GAUSSIAN ELIMINATION TECHNIQUE IS USED FOR CALCULATION      DMGR 430
C           OF THE TRIANGULAR FACTORS OF A GIVEN MATRIX.                DMGR 440
C           COMPLETE PIVOTING IS BUILT IN.                              DMGR 450
C           IN CASE OF A SINGULAR MATRIX ONLY THE TRIANGULAR FACTORS    DMGR 460
C           OF A SUBMATRIX OF MAXIMAL RANK ARE RETAINED.                DMGR 470
C           THE REMAINING PARTS OF THE RESULTANT MATRIX GIVE THE        DMGR 480
C           DEPENDENCIES OF ROWS AND THE SOLUTION OF THE HOMOGENEOUS    DMGR 490
C           MATRIX EQUATION A*X=0.                                      DMGR 500
C                                                                       DMGR 510
C     ..................................................................DMGR 520
C                                                                       DMGR 530
      SUBROUTINE DMFGR(A,M,N,EPS,IRANK,IROW,ICOL)                       DMGR 540
C                                                                       DMGR 550
C        DIMENSIONED DUMMY VARIABLES                                    DMGR 560
      DIMENSION A(1),IROW(1),ICOL(1)                                    DMGR 570
      DOUBLE PRECISION A,PIV,HOLD,SAVE                                  DMGR 580
C                                                                       DMGR 590
C        TEST OF SPECIFIED DIMENSIONS                                   DMGR 600
      IF(M)2,2,1                                                        DMGR 610
    1 IF(N)2,2,4                                                        DMGR 620
    2 IRANK=-1                                                          DMGR 630
    3 RETURN                                                            DMGR 640
C        RETURN IN CASE OF FORMAL ERRORS                                DMGR 650
C                                                                       DMGR 660
C                                                                       DMGR 670
C        INITIALIZE COLUMN INDEX VECTOR                                 DMGR 680
C        SEARCH FIRST PIVOT ELEMENT                                     DMGR 690
    4 IRANK=0                                                           DMGR 700
      PIV=0.D0                                                          DMGR 710
      JJ=0                                                              DMGR 720
      DO 6 J=1,N                                                        DMGR 730
      ICOL(J)=J                                                         DMGR 740
      DO 6 I=1,M                                                        DMGR 750
      JJ=JJ+1                                                           DMGR 760
      HOLD=A(JJ)                                                        DMGR 770
      IF(DABS(PIV)-DABS(HOLD))5,6,6                                     DMGR 780
    5 PIV=HOLD                                                          DMGR 790
      IR=I                                                              DMGR 800
      IC=J                                                              DMGR 810
    6 CONTINUE                                                          DMGR 820
C                                                                       DMGR 830
C        INITIALIZE ROW INDEX VECTOR                                    DMGR 840
      DO 7 I=1,M                                                        DMGR 850
    7 IROW(I)=I                                                         DMGR 860
C                                                                       DMGR 870
C        SET UP INTERNAL TOLERANCE                                      DMGR 880
      TOL=ABS(EPS*SNGL(PIV))                                            DMGR 890
C                                                                       DMGR 900
C        INITIALIZE ELIMINATION LOOP                                    DMGR 910
      NM=N*M                                                            DMGR 920
      DO 19 NCOL=M,NM,M                                                 DMGR 930
C                                                                       DMGR 940
C        TEST FOR FEASIBILITY OF PIVOT ELEMENT                          DMGR 950
    8 IF(ABS(SNGL(PIV))-TOL)20,20,9                                     DMGR 960
C                                                                       DMGR 970
C        UPDATE RANK                                                    DMGR 980
    9 IRANK=IRANK+1                                                     DMGR 990
C                                                                       DMGR1000
C        INTERCHANGE ROWS IF NECESSARY                                  DMGR1010
      JJ=IR-IRANK                                                       DMGR1020
      IF(JJ)12,12,10                                                    DMGR1030
   10 DO 11 J=IRANK,NM,M                                                DMGR1040
      I=J+JJ                                                            DMGR1050
      SAVE=A(J)                                                         DMGR1060
      A(J)=A(I)                                                         DMGR1070
   11 A(I)=SAVE                                                         DMGR1080
C                                                                       DMGR1090
C        UPDATE ROW INDEX VECTOR                                        DMGR1100
      JJ=IROW(IR)                                                       DMGR1110
      IROW(IR)=IROW(IRANK)                                              DMGR1120
      IROW(IRANK)=JJ                                                    DMGR1130
C                                                                       DMGR1140
C        INTERCHANGE COLUMNS IF NECESSARY                               DMGR1150
   12 JJ=(IC-IRANK)*M                                                   DMGR1160
      IF(JJ)15,15,13                                                    DMGR1170
   13 KK=NCOL                                                           DMGR1180
      DO 14 J=1,M                                                       DMGR1190
      I=KK+JJ                                                           DMGR1200
      SAVE=A(KK)                                                        DMGR1210
      A(KK)=A(I)                                                        DMGR1220
      KK=KK-1                                                           DMGR1230
   14 A(I)=SAVE                                                         DMGR1240
C                                                                       DMGR1250
C        UPDATE COLUMN INDEX VECTOR                                     DMGR1260
      JJ=ICOL(IC)                                                       DMGR1270
      ICOL(IC)=ICOL(IRANK)                                              DMGR1280
      ICOL(IRANK)=JJ                                                    DMGR1290
   15 KK=IRANK+1                                                        DMGR1300
      MM=IRANK-M                                                        DMGR1310
      LL=NCOL+MM                                                        DMGR1320
C                                                                       DMGR1330
C        TEST FOR LAST ROW                                              DMGR1340
      IF(MM)16,25,25                                                    DMGR1350
C                                                                       DMGR1360
C        TRANSFORM CURRENT SUBMATRIX AND SEARCH NEXT PIVOT              DMGR1370
   16 JJ=LL                                                             DMGR1380
      SAVE=PIV                                                          DMGR1390
      PIV=0.D0                                                          DMGR1400
      DO 19 J=KK,M                                                      DMGR1410
      JJ=JJ+1                                                           DMGR1420
      HOLD=A(JJ)/SAVE                                                   DMGR1430
      A(JJ)=HOLD                                                        DMGR1440
      L=J-IRANK                                                         DMGR1450
C                                                                       DMGR1460
C        TEST FOR LAST COLUMN                                           DMGR1470
      IF(IRANK-N)17,19,19                                               DMGR1480
   17 II=JJ                                                             DMGR1490
      DO 19 I=KK,N                                                      DMGR1500
      II=II+M                                                           DMGR1510
      MM=II-L                                                           DMGR1520
      A(II)=A(II)-HOLD*A(MM)                                            DMGR1530
      IF(DABS(A(II))-DABS(PIV))19,19,18                                 DMGR1540
   18 PIV=A(II)                                                         DMGR1550
      IR=J                                                              DMGR1560
      IC=I                                                              DMGR1570
   19 CONTINUE                                                          DMGR1580
C                                                                       DMGR1590
C        SET UP MATRIX EXPRESSING ROW DEPENDENCIES                      DMGR1600
   20 IF(IRANK-1)3,25,21                                                DMGR1610
   21 IR=LL                                                             DMGR1620
      DO 24 J=2,IRANK                                                   DMGR1630
      II=J-1                                                            DMGR1640
      IR=IR-M                                                           DMGR1650
      JJ=LL                                                             DMGR1660
      DO 23 I=KK,M                                                      DMGR1670
      HOLD=0.D0                                                         DMGR1680
      JJ=JJ+1                                                           DMGR1690
      MM=JJ                                                             DMGR1700
      IC=IR                                                             DMGR1710
      DO 22 L=1,II                                                      DMGR1720
      HOLD=HOLD+A(MM)*A(IC)                                             DMGR1730
      IC=IC-1                                                           DMGR1740
   22 MM=MM-M                                                           DMGR1750
   23 A(MM)=A(MM)-HOLD                                                  DMGR1760
   24 CONTINUE                                                          DMGR1770
C                                                                       DMGR1780
C        TEST FOR COLUMN REGULARITY                                     DMGR1790
   25 IF(N-IRANK)3,3,26                                                 DMGR1800
C                                                                       DMGR1810
C        SET UP MATRIX EXPRESSING BASIC VARIABLES IN TERMS OF FREE      DMGR1820
C        PARAMETERS (HOMOGENEOUS SOLUTION).                             DMGR1830
   26 IR=LL                                                             DMGR1840
      KK=LL+M                                                           DMGR1850
      DO 30 J=1,IRANK                                                   DMGR1860
      DO 29 I=KK,NM,M                                                   DMGR1870
      JJ=IR                                                             DMGR1880
      LL=I                                                              DMGR1890
      HOLD=0.D0                                                         DMGR1900
      II=J                                                              DMGR1910
   27 II=II-1                                                           DMGR1920
      IF(II)29,29,28                                                    DMGR1930
   28 HOLD=HOLD-A(JJ)*A(LL)                                             DMGR1940
      JJ=JJ-M                                                           DMGR1950
      LL=LL-1                                                           DMGR1960
      GOTO 27                                                           DMGR1970
   29 A(LL)=(HOLD-A(LL))/A(JJ)                                          DMGR1980
   30 IR=IR-1                                                           DMGR1990
      RETURN                                                            DMGR2000
      END                                                               DMGR2010
