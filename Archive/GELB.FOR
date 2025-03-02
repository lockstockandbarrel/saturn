C                                                                       GELB  10
C     ..................................................................GELB  20
C                                                                       GELB  30
C        SUBROUTINE GELB                                                GELB  40
C                                                                       GELB  50
C        PURPOSE                                                        GELB  60
C           TO SOLVE A SYSTEM OF SIMULTANEOUS LINEAR EQUATIONS WITH A   GELB  70
C           COEFFICIENT MATRIX OF BAND STRUCTURE.                       GELB  80
C                                                                       GELB  90
C        USAGE                                                          GELB 100
C           CALL GELB(R,A,M,N,MUD,MLD,EPS,IER)                          GELB 110
C                                                                       GELB 120
C        DESCRIPTION OF PARAMETERS                                      GELB 130
C           R      - M BY N RIGHT HAND SIDE MATRIX (DESTROYED).         GELB 140
C                    ON RETURN R CONTAINS THE SOLUTION OF THE EQUATIONS.GELB 150
C           A      - M BY M COEFFICIENT MATRIX WITH BAND STRUCTURE      GELB 160
C                    (DESTROYED).                                       GELB 170
C           M      - THE NUMBER OF EQUATIONS IN THE SYSTEM.             GELB 180
C           N      - THE NUMBER OF RIGHT HAND SIDE VECTORS.             GELB 190
C           MUD    - THE NUMBER OF UPPER CODIAGONALS (THAT MEANS        GELB 200
C                    CODIAGONALS ABOVE MAIN DIAGONAL).                  GELB 210
C           MLD    - THE NUMBER OF LOWER CODIAGONALS (THAT MEANS        GELB 220
C                    CODIAGONALS BELOW MAIN DIAGONAL).                  GELB 230
C           EPS    - AN INPUT CONSTANT WHICH IS USED AS RELATIVE        GELB 240
C                    TOLERANCE FOR TEST ON LOSS OF SIGNIFICANCE.        GELB 250
C           IER    - RESULTING ERROR PARAMETER CODED AS FOLLOWS         GELB 260
C                    IER=0  - NO ERROR,                                 GELB 270
C                    IER=-1 - NO RESULT BECAUSE OF WRONG INPUT PARAME-  GELB 280
C                             TERS M,MUD,MLD OR BECAUSE OF PIVOT ELEMENTGELB 290
C                             AT ANY ELIMINATION STEP EQUAL TO 0,       GELB 300
C                    IER=K  - WARNING DUE TO POSSIBLE LOSS OF SIGNIFI-  GELB 310
C                             CANCE INDICATED AT ELIMINATION STEP K+1,  GELB 320
C                             WHERE PIVOT ELEMENT WAS LESS THAN OR      GELB 330
C                             EQUAL TO THE INTERNAL TOLERANCE EPS TIMES GELB 340
C                             ABSOLUTELY GREATEST ELEMENT OF MATRIX A.  GELB 350
C                                                                       GELB 360
C        REMARKS                                                        GELB 370
C           BAND MATRIX A IS ASSUMED TO BE STORED ROWWISE IN THE FIRST  GELB 380
C           ME SUCCESSIVE STORAGE LOCATIONS OF TOTALLY NEEDED MA        GELB 390
C           STORAGE LOCATIONS, WHERE                                    GELB 400
C             MA=M*MC-ML*(ML+1)/2    AND    ME=MA-MU*(MU+1)/2    WITH   GELB 410
C             MC=MIN(M,1+MUD+MLD),  ML=MC-1-MLD,  MU=MC-1-MUD.          GELB 420
C           RIGHT HAND SIDE MATRIX R IS ASSUMED TO BE STORED COLUMNWISE GELB 430
C           IN N*M SUCCESSIVE STORAGE LOCATIONS. ON RETURN SOLUTION     GELB 440
C           MATRIX R IS STORED COLUMNWISE TOO.                          GELB 450
C           INPUT PARAMETERS M, MUD, MLD SHOULD SATISFY THE FOLLOWING   GELB 460
C           RESTRICTIONS     MUD NOT LESS THAN ZERO                     GELB 470
C                            MLD NOT LESS THAN ZERO                     GELB 480
C                            MUD+MLD NOT GREATER THAN 2*M-2.            GELB 490
C           NO ACTION BESIDES ERROR MESSAGE IER=-1 TAKES PLACE IF THESE GELB 500
C           RESTRICTIONS ARE NOT SATISFIED.                             GELB 510
C           THE PROCEDURE GIVES RESULTS IF THE RESTRICTIONS ON INPUT    GELB 520
C           PARAMETERS ARE SATISFIED AND IF PIVOT ELEMENTS AT ALL       GELB 530
C           ELIMINATION STEPS ARE DIFFERENT FROM 0. HOWEVER WARNING     GELB 540
C           IER=K - IF GIVEN - INDICATES POSSIBLE LOSS OF SIGNIFICANCE. GELB 550
C           IN CASE OF A WELL SCALED MATRIX A AND APPROPRIATE TOLERANCE GELB 560
C           EPS, IER=K MAY BE INTERPRETED THAT MATRIX A HAS THE RANK K. GELB 570
C           NO WARNING IS GIVEN IF MATRIX A HAS NO LOWER CODIAGONAL.    GELB 580
C                                                                       GELB 590
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  GELB 600
C           NONE                                                        GELB 610
C                                                                       GELB 620
C        METHOD                                                         GELB 630
C           SOLUTION IS DONE BY MEANS OF GAUSS ELIMINATION WITH         GELB 640
C           COLUMN PIVOTING ONLY, IN ORDER TO PRESERVE BAND STRUCTURE   GELB 650
C           IN REMAINING COEFFICIENT MATRICES.                          GELB 660
C                                                                       GELB 670
C     ..................................................................GELB 680
C                                                                       GELB 690
      SUBROUTINE GELB(R,A,M,N,MUD,MLD,EPS,IER)                          GELB 700
C                                                                       GELB 710
C                                                                       GELB 720
      DIMENSION R(1),A(1)                                               GELB 730
C                                                                       GELB 740
C     TEST ON WRONG INPUT PARAMETERS                                    GELB 750
      IF(MLD)47,1,1                                                     GELB 760
    1 IF(MUD)47,2,2                                                     GELB 770
    2 MC=1+MLD+MUD                                                      GELB 780
      IF(MC+1-M-M)3,3,47                                                GELB 790
C                                                                       GELB 800
C     PREPARE INTEGER PARAMETERS                                        GELB 810
C        MC=NUMBER OF COLUMNS IN MATRIX A                               GELB 820
C        MU=NUMBER OF ZEROS TO BE INSERTED IN FIRST ROW OF MATRIX A     GELB 830
C        ML=NUMBER OF MISSING ELEMENTS IN LAST ROW OF MATRIX A          GELB 840
C        MR=INDEX OF LAST ROW IN MATRIX A WITH MC ELEMENTS              GELB 850
C        MZ=TOTAL NUMBER OF ZEROS TO BE INSERTED IN MATRIX A            GELB 860
C        MA=TOTAL NUMBER OF STORAGE LOCATIONS NECESSARY FOR MATRIX A    GELB 870
C        NM=NUMBER OF ELEMENTS IN MATRIX R                              GELB 880
    3 IF(MC-M)5,5,4                                                     GELB 890
    4 MC=M                                                              GELB 900
    5 MU=MC-MUD-1                                                       GELB 910
      ML=MC-MLD-1                                                       GELB 920
      MR=M-ML                                                           GELB 930
      MZ=(MU*(MU+1))/2                                                  GELB 940
      MA=M*MC-(ML*(ML+1))/2                                             GELB 950
      NM=N*M                                                            GELB 960
C                                                                       GELB 970
C     MOVE ELEMENTS BACKWARD AND SEARCH FOR ABSOLUTELY GREATEST ELEMENT GELB 980
C     (NOT NECESSARY IN CASE OF A MATRIX WITHOUT LOWER CODIAGONALS)     GELB 990
      IER=0                                                             GELB1000
      PIV=0.                                                            GELB1010
      IF(MLD)14,14,6                                                    GELB1020
    6 JJ=MA                                                             GELB1030
      J=MA-MZ                                                           GELB1040
      KST=J                                                             GELB1050
      DO 9 K=1,KST                                                      GELB1060
      TB=A(J)                                                           GELB1070
      A(JJ)=TB                                                          GELB1080
      TB=ABS(TB)                                                        GELB1090
      IF(TB-PIV)8,8,7                                                   GELB1100
    7 PIV=TB                                                            GELB1110
    8 J=J-1                                                             GELB1120
    9 JJ=JJ-1                                                           GELB1130
C                                                                       GELB1140
C     INSERT ZEROS IN FIRST MU ROWS (NOT NECESSARY IN CASE MZ=0)        GELB1150
      IF(MZ)14,14,10                                                    GELB1160
   10 JJ=1                                                              GELB1170
      J=1+MZ                                                            GELB1180
      IC=1+MUD                                                          GELB1190
      DO 13 I=1,MU                                                      GELB1200
      DO 12 K=1,MC                                                      GELB1210
      A(JJ)=0.                                                          GELB1220
      IF(K-IC)11,11,12                                                  GELB1230
   11 A(JJ)=A(J)                                                        GELB1240
      J=J+1                                                             GELB1250
   12 JJ=JJ+1                                                           GELB1260
   13 IC=IC+1                                                           GELB1270
C                                                                       GELB1280
C     GENERATE TEST VALUE FOR SINGULARITY                               GELB1290
   14 TOL=EPS*PIV                                                       GELB1300
C                                                                       GELB1310
C                                                                       GELB1320
C     START DECOMPOSITION LOOP                                          GELB1330
      KST=1                                                             GELB1340
      IDST=MC                                                           GELB1350
      IC=MC-1                                                           GELB1360
      DO 38 K=1,M                                                       GELB1370
      IF(K-MR-1)16,16,15                                                GELB1380
   15 IDST=IDST-1                                                       GELB1390
   16 ID=IDST                                                           GELB1400
      ILR=K+MLD                                                         GELB1410
      IF(ILR-M)18,18,17                                                 GELB1420
   17 ILR=M                                                             GELB1430
   18 II=KST                                                            GELB1440
C                                                                       GELB1450
C     PIVOT SEARCH IN FIRST COLUMN (ROW INDEXES FROM I=K UP TO I=ILR)   GELB1460
      PIV=0.                                                            GELB1470
      DO 22 I=K,ILR                                                     GELB1480
      TB=ABS(A(II))                                                     GELB1490
      IF(TB-PIV)20,20,19                                                GELB1500
   19 PIV=TB                                                            GELB1510
      J=I                                                               GELB1520
      JJ=II                                                             GELB1530
   20 IF(I-MR)22,22,21                                                  GELB1540
   21 ID=ID-1                                                           GELB1550
   22 II=II+ID                                                          GELB1560
C                                                                       GELB1570
C     TEST ON SINGULARITY                                               GELB1580
      IF(PIV)47,47,23                                                   GELB1590
   23 IF(IER)26,24,26                                                   GELB1600
   24 IF(PIV-TOL)25,25,26                                               GELB1610
   25 IER=K-1                                                           GELB1620
   26 PIV=1./A(JJ)                                                      GELB1630
C                                                                       GELB1640
C     PIVOT ROW REDUCTION AND ROW INTERCHANGE IN RIGHT HAND SIDE R      GELB1650
      ID=J-K                                                            GELB1660
      DO 27 I=K,NM,M                                                    GELB1670
      II=I+ID                                                           GELB1680
      TB=PIV*R(II)                                                      GELB1690
      R(II)=R(I)                                                        GELB1700
   27 R(I)=TB                                                           GELB1710
C                                                                       GELB1720
C     PIVOT ROW REDUCTION AND ROW INTERCHANGE IN COEFFICIENT MATRIX A   GELB1730
      II=KST                                                            GELB1740
      J=JJ+IC                                                           GELB1750
      DO 28 I=JJ,J                                                      GELB1760
      TB=PIV*A(I)                                                       GELB1770
      A(I)=A(II)                                                        GELB1780
      A(II)=TB                                                          GELB1790
   28 II=II+1                                                           GELB1800
C                                                                       GELB1810
C     ELEMENT REDUCTION                                                 GELB1820
      IF(K-ILR)29,34,34                                                 GELB1830
   29 ID=KST                                                            GELB1840
      II=K+1                                                            GELB1850
      MU=KST+1                                                          GELB1860
      MZ=KST+IC                                                         GELB1870
      DO 33 I=II,ILR                                                    GELB1880
C                                                                       GELB1890
C     IN MATRIX A                                                       GELB1900
      ID=ID+MC                                                          GELB1910
      JJ=I-MR-1                                                         GELB1920
      IF(JJ)31,31,30                                                    GELB1930
   30 ID=ID-JJ                                                          GELB1940
   31 PIV=-A(ID)                                                        GELB1950
      J=ID+1                                                            GELB1960
      DO 32 JJ=MU,MZ                                                    GELB1970
      A(J-1)=A(J)+PIV*A(JJ)                                             GELB1980
   32 J=J+1                                                             GELB1990
      A(J-1)=0.                                                         GELB2000
C                                                                       GELB2010
C     IN MATRIX R                                                       GELB2020
      J=K                                                               GELB2030
      DO 33 JJ=I,NM,M                                                   GELB2040
      R(JJ)=R(JJ)+PIV*R(J)                                              GELB2050
   33 J=J+M                                                             GELB2060
   34 KST=KST+MC                                                        GELB2070
      IF(ILR-MR)36,35,35                                                GELB2080
   35 IC=IC-1                                                           GELB2090
   36 ID=K-MR                                                           GELB2100
      IF(ID)38,38,37                                                    GELB2110
   37 KST=KST-ID                                                        GELB2120
   38 CONTINUE                                                          GELB2130
C     END OF DECOMPOSITION LOOP                                         GELB2140
C                                                                       GELB2150
C                                                                       GELB2160
C     BACK SUBSTITUTION                                                 GELB2170
      IF(MC-1)46,46,39                                                  GELB2180
   39 IC=2                                                              GELB2190
      KST=MA+ML-MC+2                                                    GELB2200
      II=M                                                              GELB2210
      DO 45 I=2,M                                                       GELB2220
      KST=KST-MC                                                        GELB2230
      II=II-1                                                           GELB2240
      J=II-MR                                                           GELB2250
      IF(J)41,41,40                                                     GELB2260
   40 KST=KST+J                                                         GELB2270
   41 DO 43 J=II,NM,M                                                   GELB2280
      TB=R(J)                                                           GELB2290
      MZ=KST+IC-2                                                       GELB2300
      ID=J                                                              GELB2310
      DO 42 JJ=KST,MZ                                                   GELB2320
      ID=ID+1                                                           GELB2330
   42 TB=TB-A(JJ)*R(ID)                                                 GELB2340
   43 R(J)=TB                                                           GELB2350
      IF(IC-MC)44,45,45                                                 GELB2360
   44 IC=IC+1                                                           GELB2370
   45 CONTINUE                                                          GELB2380
   46 RETURN                                                            GELB2390
C                                                                       GELB2400
C                                                                       GELB2410
C     ERROR RETURN                                                      GELB2420
   47 IER=-1                                                            GELB2430
      RETURN                                                            GELB2440
      END                                                               GELB2450
