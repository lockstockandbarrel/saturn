C                                                                       MPRC  10
C     ..................................................................MPRC  20
C                                                                       MPRC  30
C        SUBROUTINE MPRC                                                MPRC  40
C                                                                       MPRC  50
C        PURPOSE                                                        MPRC  60
C           TO PERMUTE THE ROWS OR COLUMNS OF A GIVEN MATRIX ACCORDING  MPRC  70
C           TO A GIVEN TRANSPOSITION VECTOR OR ITS INVERSE.  (SEE THE   MPRC  80
C           DISCUSSION ON PERMUTATIONS FOR DEFINITIONS AND NOTATION.)   MPRC  90
C                                                                       MPRC 100
C        USAGE                                                          MPRC 110
C           CALL MPRC(A,M,N,ITRA,INV,IROCO,IER)                         MPRC 120
C                                                                       MPRC 130
C        DESCRIPTION OF PARAMETERS                                      MPRC 140
C           A     - GIVEN M BY N MATRIX AND RESULTING PERMUTED MATRIX   MPRC 150
C           M     - NUMBER OF ROWS OF A                                 MPRC 160
C           N     - NUMBER OF COLUMNS OF A                              MPRC 170
C           ITRA  - GIVEN TRANSPOSITION VECTOR (DIMENSION M IF ROWS ARE MPRC 180
C                   PERMUTED, N IF COLUMNS ARE PERMUTED)                MPRC 190
C           INV   - INPUT PARAMETER                                     MPRC 200
C                   INV NON-ZERO  -  PERMUTE ACCORDING TO ITRA          MPRC 210
C                   INV    =   0  -  PERMUTE ACCORDING TO ITRA INVERSE  MPRC 220
C           IROCO - INPUT PARAMETER                                     MPRC 230
C                   IROCO NON-ZERO  -  PERMUTE THE COLUMNS OF A         MPRC 240
C                   IROCO    =   0  -  PERMUTE THE ROWS OF A            MPRC 250
C           IER   - RESULTING ERROR PARAMETER                           MPRC 260
C                   IER = -1  -  M AND N ARE NOT BOTH POSITIVE          MPRC 270
C                   IER =  0  -  NO ERROR                               MPRC 280
C                   IER =  1  -  ITRA IS NOT A TRANSPOSITION VECTOR ON  MPRC 290
C                                1,...,M IF ROWS ARE PERMUTED, 1,...,N  MPRC 300
C                                IF COLUMNS ARE PERMUTED                MPRC 310
C                                                                       MPRC 320
C        REMARKS                                                        MPRC 330
C           (1)  IF IER=-1 THERE IS NO COMPUTATION.                     MPRC 340
C           (2)  IF IER= 1, THEN COMPUTATION HAS BEEN UNSUCCESSFUL DUE  MPRC 350
C                TO ERROR, BUT THE MATRIX A WILL REFLECT THE ROW OR     MPRC 360
C                COLUMN INTERCHANGES PERFORMED BEFORE THE ERROR WAS     MPRC 370
C                DETECTED.                                              MPRC 380
C           (3)  THE MATRIX A IS ASSUMED TO BE STORED COLUMNWISE.       MPRC 390
C                                                                       MPRC 400
C        SUBROUTINES AND SUBPROGRAMS REQUIRED                           MPRC 410
C           NONE                                                        MPRC 420
C                                                                       MPRC 430
C        METHOD                                                         MPRC 440
C           THE ROWS OR COLUMNS ARE PERMUTED ELEMENTWISE, INTERCHANGING MPRC 450
C           ROW OR COLUMN 1 AND ITRA(1),...,ROW OR COLUMN K AND ITRA(K) MPRC 460
C           IN THAT ORDER IF INV=0, AND OTHERWISE INTERCHANGING ROW OR  MPRC 470
C           COLUMN K AND ITRA(K),...,ROW OR COLUMN 1 AND ITRA(1), WHERE MPRC 480
C           K IS M OR N DEPENDING ON WHETHER WE PERMUTE ROWS OR COLUMNS.MPRC 490
C                                                                       MPRC 500
C     ..................................................................MPRC 510
C                                                                       MPRC 520
      SUBROUTINE MPRC(A,M,N,ITRA,INV,IROCO,IER)                         MPRC 530
C                                                                       MPRC 540
C                                                                       MPRC 550
      DIMENSION A(1),ITRA(1)                                            MPRC 560
C                                                                       MPRC 570
C        TEST OF DIMENSIONS                                             MPRC 580
      IF(M)14,14,1                                                      MPRC 590
    1 IF(N)14,14,2                                                      MPRC 600
C                                                                       MPRC 610
C        DETERMINE WHICH ARE TO BE PERMUTED-THE ROWS OR THE COLUMNS     MPRC 620
    2 IF(IROCO)3,4,3                                                    MPRC 630
C                                                                       MPRC 640
C        INITIALIZE FOR COLUMN INTERCHANGES                             MPRC 650
    3 MM=M                                                              MPRC 660
      MMM=-1                                                            MPRC 670
      L=M                                                               MPRC 680
      LL=N                                                              MPRC 690
      GO TO 5                                                           MPRC 700
C                                                                       MPRC 710
C        INITIALIZE FOR ROW INTERCHANGES                                MPRC 720
    4 MM=1                                                              MPRC 730
      MMM=M                                                             MPRC 740
      L=N                                                               MPRC 750
      LL=M                                                              MPRC 760
C                                                                       MPRC 770
C        INITIALIZE LOOP OVER ALL ROWS OR COLUMNS                       MPRC 780
    5 IA=1                                                              MPRC 790
      ID=1                                                              MPRC 800
C                                                                       MPRC 810
C        TEST FOR INVERSE OPERATION                                     MPRC 820
      IF(INV)6,7,6                                                      MPRC 830
    6 IA=LL                                                             MPRC 840
      ID=-1                                                             MPRC 850
    7 DO 12 I=1,LL                                                      MPRC 860
      K=ITRA(IA)                                                        MPRC 870
      IF(K-IA)8,12,9                                                    MPRC 880
    8 IF(K)13,13,10                                                     MPRC 890
    9 IF(LL-K)13,10,10                                                  MPRC 900
C                                                                       MPRC 910
C        INITIALIZE ROW OR COLUMN INTERCHANGE                           MPRC 920
   10 IL=IA*MM                                                          MPRC 930
      K=K*MM                                                            MPRC 940
C                                                                       MPRC 950
C        PERFORM ROW OR COLUMN INTERCHANGE                              MPRC 960
      DO 11 J=1,L                                                       MPRC 970
      SAVE=A(IL)                                                        MPRC 980
      A(IL)=A(K)                                                        MPRC 990
      A(K)=SAVE                                                         MPRC1000
      K=K+MMM                                                           MPRC1010
   11 IL=IL+MMM                                                         MPRC1020
C                                                                       MPRC1030
C        ADDRESS NEXT INTERCHANGE STEP                                  MPRC1040
   12 IA=IA+ID                                                          MPRC1050
C                                                                       MPRC1060
C        NORMAL EXIT                                                    MPRC1070
      IER=0                                                             MPRC1080
      RETURN                                                            MPRC1090
C                                                                       MPRC1100
C        ERROR RETURN IN CASE ITRA IS NOT A TRANSPOSITION VECTOR        MPRC1110
   13 IER=1                                                             MPRC1120
      RETURN                                                            MPRC1130
C                                                                       MPRC1140
C        ERROR RETURN IN CASE OF ILLEGAL DIMENSIONS                     MPRC1150
   14 IER=-1                                                            MPRC1160
      RETURN                                                            MPRC1170
      END                                                               MPRC1180
