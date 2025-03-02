C                                                                       DELB  10
C     ..................................................................DELB  20
C                                                                       DELB  30
C        SUBROUTINE DGELB                                               DELB  40
C                                                                       DELB  50
C        PURPOSE                                                        DELB  60
C           TO SOLVE A SYSTEM OF SIMULTANEOUS LINEAR EQUATIONS WITH A   DELB  70
C           COEFFICIENT MATRIX OF BAND STRUCTURE.                       DELB  80
C                                                                       DELB  90
C        USAGE                                                          DELB 100
C           CALL DGELB(R,A,M,N,MUD,MLD,EPS,IER)                         DELB 110
C                                                                       DELB 120
C        DESCRIPTION OF PARAMETERS                                      DELB 130
C           R      - DOUBLE PRECISION M BY N RIGHT HAND SIDE MATRIX     DELB 140
C                    (DESTROYED). ON RETURN R CONTAINS THE SOLUTION     DELB 150
C                    OF THE EQUATIONS.                                  DELB 160
C           A      - DOUBLE PRECISION M BY M COEFFICIENT MATRIX WITH    DELB 170
C                    BAND STRUCTURE (DESTROYED).                        DELB 180
C           M      - THE NUMBER OF EQUATIONS IN THE SYSTEM.             DELB 190
C           N      - THE NUMBER OF RIGHT HAND SIDE VECTORS.             DELB 200
C           MUD    - THE NUMBER OF UPPER CODIAGONALS (THAT MEANS        DELB 210
C                    CODIAGONALS ABOVE MAIN DIAGONAL).                  DELB 220
C           MLD    - THE NUMBER OF LOWER CODIAGONALS (THAT MEANS        DELB 230
C                    CODIAGONALS BELOW MAIN DIAGONAL).                  DELB 240
C           EPS    - SINGLE PRECISION INPUT CONSTANT WHICH IS USED AS   DELB 250
C                    RELATIVE TOLERANCE FOR TEST ON LOSS OF             DELB 260
C                    SIGNIFICANCE.                                      DELB 270
C           IER    - RESULTING ERROR PARAMETER CODED AS FOLLOWS         DELB 280
C                    IER=0  - NO ERROR,                                 DELB 290
C                    IER=-1 - NO RESULT BECAUSE OF WRONG INPUT PARAME-  DELB 300
C                             TERS M,MUD,MLD OR BECAUSE OF PIVOT ELEMENTDELB 310
C                             AT ANY ELIMINATION STEP EQUAL TO 0,       DELB 320
C                    IER=K  - WARNING DUE TO POSSIBLE LOSS OF SIGNIFI-  DELB 330
C                             CANCE INDICATED AT ELIMINATION STEP K+1,  DELB 340
C                             WHERE PIVOT ELEMENT WAS LESS THAN OR      DELB 350
C                             EQUAL TO THE INTERNAL TOLERANCE EPS TIMES DELB 360
C                             ABSOLUTELY GREATEST ELEMENT OF MATRIX A.  DELB 370
C                                                                       DELB 380
C        REMARKS                                                        DELB 390
C           BAND MATRIX A IS ASSUMED TO BE STORED ROWWISE IN THE FIRST  DELB 400
C           ME SUCCESSIVE STORAGE LOCATIONS OF TOTALLY NEEDED MA        DELB 410
C           STORAGE LOCATIONS, WHERE                                    DELB 420
C             MA=M*MC-ML*(ML+1)/2    AND    ME=MA-MU*(MU+1)/2    WITH   DELB 430
C             MC=MIN(M,1+MUD+MLD),  ML=MC-1-MLD,  MU=MC-1-MUD.          DELB 440
C           RIGHT HAND SIDE MATRIX R IS ASSUMED TO BE STORED COLUMNWISE DELB 450
C           IN N*M SUCCESSIVE STORAGE LOCATIONS. ON RETURN SOLUTION     DELB 460
C           MATRIX R IS STORED COLUMNWISE TOO.                          DELB 470
C           INPUT PARAMETERS M, MUD, MLD SHOULD SATISFY THE FOLLOWING   DELB 480
C           RESTRICTIONS     MUD NOT LESS THAN ZERO                     DELB 490
C                            MLD NOT LESS THAN ZERO                     DELB 500
C                            MUD+MLD NOT GREATER THAN 2*M-2.            DELB 510
C           NO ACTION BESIDES ERROR MESSAGE IER=-1 TAKES PLACE IF THESE DELB 520
C           RESTRICTIONS ARE NOT SATISFIED.                             DELB 530
C           THE PROCEDURE GIVES RESULTS IF THE RESTRICTIONS ON INPUT    DELB 540
C           PARAMETERS ARE SATISFIED AND IF PIVOT ELEMENTS AT ALL       DELB 550
C           ELIMINATION STEPS ARE DIFFERENT FROM 0. HOWEVER WARNING     DELB 560
C           IER=K - IF GIVEN - INDICATES POSSIBLE LOSS OF SIGNIFICANCE. DELB 570
C           IN CASE OF A WELL SCALED MATRIX A AND APPROPRIATE TOLERANCE DELB 580
C           EPS, IER=K MAY BE INTERPRETED THAT MATRIX A HAS THE RANK K. DELB 590
C           NO WARNING IS GIVEN IF MATRIX A HAS NO LOWER CODIAGONAL.    DELB 600
C                                                                       DELB 610
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DELB 620
C           NONE                                                        DELB 630
C                                                                       DELB 640
C        METHOD                                                         DELB 650
C           SOLUTION IS DONE BY MEANS OF GAUSS ELIMINATION WITH         DELB 660
C           COLUMN PIVOTING ONLY, IN ORDER TO PRESERVE BAND STRUCTURE   DELB 670
C           IN REMAINING COEFFICIENT MATRICES.                          DELB 680
C                                                                       DELB 690
C     ..................................................................DELB 700
C                                                                       DELB 710
      SUBROUTINE DGELB(R,A,M,N,MUD,MLD,EPS,IER)                         DELB 720
C                                                                       DELB 730
C                                                                       DELB 740
      DIMENSION R(1),A(1)                                               DELB 750
      DOUBLE PRECISION R,A,PIV,TB,TOL                                   DELB 760
C                                                                       DELB 770
C     TEST ON WRONG INPUT PARAMETERS                                    DELB 780
      IF(MLD)47,1,1                                                     DELB 790
    1 IF(MUD)47,2,2                                                     DELB 800
    2 MC=1+MLD+MUD                                                      DELB 810
      IF(MC+1-M-M)3,3,47                                                DELB 820
C                                                                       DELB 830
C     PREPARE INTEGER PARAMETERS                                        DELB 840
C        MC=NUMBER OF COLUMNS IN MATRIX A                               DELB 850
C        MU=NUMBER OF ZEROS TO BE INSERTED IN FIRST ROW OF MATRIX A     DELB 860
C        ML=NUMBER OF MISSING ELEMENTS IN LAST ROW OF MATRIX A          DELB 870
C        MR=INDEX OF LAST ROW IN MATRIX A WITH MC ELEMENTS              DELB 880
C        MZ=TOTAL NUMBER OF ZEROS TO BE INSERTED IN MATRIX A            DELB 890
C        MA=TOTAL NUMBER OF STORAGE LOCATIONS NECESSARY FOR MATRIX A    DELB 900
C        NM=NUMBER OF ELEMENTS IN MATRIX R                              DELB 910
    3 IF(MC-M)5,5,4                                                     DELB 920
    4 MC=M                                                              DELB 930
    5 MU=MC-MUD-1                                                       DELB 940
      ML=MC-MLD-1                                                       DELB 950
      MR=M-ML                                                           DELB 960
      MZ=(MU*(MU+1))/2                                                  DELB 970
      MA=M*MC-(ML*(ML+1))/2                                             DELB 980
      NM=N*M                                                            DELB 990
C                                                                       DELB1000
C     MOVE ELEMENTS BACKWARD AND SEARCH FOR ABSOLUTELY GREATEST ELEMENT DELB1010
C     (NOT NECESSARY IN CASE OF A MATRIX WITHOUT LOWER CODIAGONALS)     DELB1020
      IER=0                                                             DELB1030
      PIV=0.D0                                                          DELB1040
      IF(MLD)14,14,6                                                    DELB1050
    6 JJ=MA                                                             DELB1060
      J=MA-MZ                                                           DELB1070
      KST=J                                                             DELB1080
      DO 9 K=1,KST                                                      DELB1090
      TB=A(J)                                                           DELB1100
      A(JJ)=TB                                                          DELB1110
      TB=DABS(TB)                                                       DELB1120
      IF(TB-PIV)8,8,7                                                   DELB1130
    7 PIV=TB                                                            DELB1140
    8 J=J-1                                                             DELB1150
    9 JJ=JJ-1                                                           DELB1160
C                                                                       DELB1170
C     INSERT ZEROS IN FIRST MU ROWS (NOT NECESSARY IN CASE MZ=0)        DELB1180
      IF(MZ)14,14,10                                                    DELB1190
   10 JJ=1                                                              DELB1200
      J=1+MZ                                                            DELB1210
      IC=1+MUD                                                          DELB1220
      DO 13 I=1,MU                                                      DELB1230
      DO 12 K=1,MC                                                      DELB1240
      A(JJ)=0.D0                                                        DELB1250
      IF(K-IC)11,11,12                                                  DELB1260
   11 A(JJ)=A(J)                                                        DELB1270
      J=J+1                                                             DELB1280
   12 JJ=JJ+1                                                           DELB1290
   13 IC=IC+1                                                           DELB1300
C                                                                       DELB1310
C     GENERATE TEST VALUE FOR SINGULARITY                               DELB1320
   14 TOL=EPS*PIV                                                       DELB1330
C                                                                       DELB1340
C                                                                       DELB1350
C     START DECOMPOSITION LOOP                                          DELB1360
      KST=1                                                             DELB1370
      IDST=MC                                                           DELB1380
      IC=MC-1                                                           DELB1390
      DO 38 K=1,M                                                       DELB1400
      IF(K-MR-1)16,16,15                                                DELB1410
   15 IDST=IDST-1                                                       DELB1420
   16 ID=IDST                                                           DELB1430
      ILR=K+MLD                                                         DELB1440
      IF(ILR-M)18,18,17                                                 DELB1450
   17 ILR=M                                                             DELB1460
   18 II=KST                                                            DELB1470
C                                                                       DELB1480
C     PIVOT SEARCH IN FIRST COLUMN (ROW INDEXES FROM I=K UP TO I=ILR)   DELB1490
      PIV=0.D0                                                          DELB1500
      DO 22 I=K,ILR                                                     DELB1510
      TB=DABS(A(II))                                                    DELB1520
      IF(TB-PIV)20,20,19                                                DELB1530
   19 PIV=TB                                                            DELB1540
      J=I                                                               DELB1550
      JJ=II                                                             DELB1560
   20 IF(I-MR)22,22,21                                                  DELB1570
   21 ID=ID-1                                                           DELB1580
   22 II=II+ID                                                          DELB1590
C                                                                       DELB1600
C     TEST ON SINGULARITY                                               DELB1610
      IF(PIV)47,47,23                                                   DELB1620
   23 IF(IER)26,24,26                                                   DELB1630
   24 IF(PIV-TOL)25,25,26                                               DELB1640
   25 IER=K-1                                                           DELB1650
   26 PIV=1.D0/A(JJ)                                                    DELB1660
C                                                                       DELB1670
C     PIVOT ROW REDUCTION AND ROW INTERCHANGE IN RIGHT HAND SIDE R      DELB1680
      ID=J-K                                                            DELB1690
      DO 27 I=K,NM,M                                                    DELB1700
      II=I+ID                                                           DELB1710
      TB=PIV*R(II)                                                      DELB1720
      R(II)=R(I)                                                        DELB1730
   27 R(I)=TB                                                           DELB1740
C                                                                       DELB1750
C     PIVOT ROW REDUCTION AND ROW INTERCHANGE IN COEFFICIENT MATRIX A   DELB1760
      II=KST                                                            DELB1770
      J=JJ+IC                                                           DELB1780
      DO 28 I=JJ,J                                                      DELB1790
      TB=PIV*A(I)                                                       DELB1800
      A(I)=A(II)                                                        DELB1810
      A(II)=TB                                                          DELB1820
   28 II=II+1                                                           DELB1830
C                                                                       DELB1840
C     ELEMENT REDUCTION                                                 DELB1850
      IF(K-ILR)29,34,34                                                 DELB1860
   29 ID=KST                                                            DELB1870
      II=K+1                                                            DELB1880
      MU=KST+1                                                          DELB1890
      MZ=KST+IC                                                         DELB1900
      DO 33 I=II,ILR                                                    DELB1910
C                                                                       DELB1920
C     IN MATRIX A                                                       DELB1930
      ID=ID+MC                                                          DELB1940
      JJ=I-MR-1                                                         DELB1950
      IF(JJ)31,31,30                                                    DELB1960
   30 ID=ID-JJ                                                          DELB1970
   31 PIV=-A(ID)                                                        DELB1980
      J=ID+1                                                            DELB1990
      DO 32 JJ=MU,MZ                                                    DELB2000
      A(J-1)=A(J)+PIV*A(JJ)                                             DELB2010
   32 J=J+1                                                             DELB2020
      A(J-1)=0.D0                                                       DELB2030
C                                                                       DELB2040
C     IN MATRIX R                                                       DELB2050
      J=K                                                               DELB2060
      DO 33 JJ=I,NM,M                                                   DELB2070
      R(JJ)=R(JJ)+PIV*R(J)                                              DELB2080
   33 J=J+M                                                             DELB2090
   34 KST=KST+MC                                                        DELB2100
      IF(ILR-MR)36,35,35                                                DELB2110
   35 IC=IC-1                                                           DELB2120
   36 ID=K-MR                                                           DELB2130
      IF(ID)38,38,37                                                    DELB2140
   37 KST=KST-ID                                                        DELB2150
   38 CONTINUE                                                          DELB2160
C     END OF DECOMPOSITION LOOP                                         DELB2170
C                                                                       DELB2180
C                                                                       DELB2190
C     BACK SUBSTITUTION                                                 DELB2200
      IF(MC-1)46,46,39                                                  DELB2210
   39 IC=2                                                              DELB2220
      KST=MA+ML-MC+2                                                    DELB2230
      II=M                                                              DELB2240
      DO 45 I=2,M                                                       DELB2250
      KST=KST-MC                                                        DELB2260
      II=II-1                                                           DELB2270
      J=II-MR                                                           DELB2280
      IF(J)41,41,40                                                     DELB2290
   40 KST=KST+J                                                         DELB2300
   41 DO 43 J=II,NM,M                                                   DELB2310
      TB=R(J)                                                           DELB2320
      MZ=KST+IC-2                                                       DELB2330
      ID=J                                                              DELB2340
      DO 42 JJ=KST,MZ                                                   DELB2350
      ID=ID+1                                                           DELB2360
   42 TB=TB-A(JJ)*R(ID)                                                 DELB2370
   43 R(J)=TB                                                           DELB2380
      IF(IC-MC)44,45,45                                                 DELB2390
   44 IC=IC+1                                                           DELB2400
   45 CONTINUE                                                          DELB2410
   46 RETURN                                                            DELB2420
C                                                                       DELB2430
C                                                                       DELB2440
C     ERROR RETURN                                                      DELB2450
   47 IER=-1                                                            DELB2460
      RETURN                                                            DELB2470
      END                                                               DELB2480
