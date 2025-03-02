C                                                                       DMPR  10
C     ..................................................................DMPR  20
C                                                                       DMPR  30
C        SUBROUTINE DMPRC                                               DMPR  40
C                                                                       DMPR  50
C        PURPOSE                                                        DMPR  60
C           TO PERMUTE THE ROWS OR COLUMNS OF A GIVEN MATRIX ACCORDING  DMPR  70
C           TO A GIVEN TRANSPOSITION VECTOR OR ITS INVERSE.  (SEE THE   DMPR  80
C           DISCUSSION ON PERMUTATIONS FOR DEFINITIONS AND NOTATION.)   DMPR  90
C                                                                       DMPR 100
C        USAGE                                                          DMPR 110
C           CALL DMPRC(A,M,N,ITRA,INV,IROCO,IER)                        DMPR 120
C                                                                       DMPR 130
C        DESCRIPTION OF PARAMETERS                                      DMPR 140
C           A     - GIVEN DOUBLE PRECISION M BY N MATRIX AND RESULTING  DMPR 150
C                   PERMUTED MATRIX                                     DMPR 160
C           M     - NUMBER OF ROWS OF A                                 DMPR 170
C           N     - NUMBER OF COLUMNS OF A                              DMPR 180
C           ITRA  - GIVEN TRANSPOSITION VECTOR (DIMENSION M IF ROWS ARE DMPR 190
C                   PERMUTED, N IF COLUMNS ARE PERMUTED)                DMPR 200
C           INV   - INPUT PARAMETER                                     DMPR 210
C                   INV NON-ZERO  -  PERMUTE ACCORDING TO ITRA          DMPR 220
C                   INV    =   0  -  PERMUTE ACCORDING TO ITRA INVERSE  DMPR 230
C           IROCO - INPUT PARAMETER                                     DMPR 240
C                   IROCO NON-ZERO  -  PERMUTE THE COLUMNS OF A         DMPR 250
C                   IROCO    =   0  -  PERMUTE THE ROWS OF A            DMPR 260
C           IER   - RESULTING ERROR PARAMETER                           DMPR 270
C                   IER = -1  -  M AND N ARE NOT BOTH POSITIVE          DMPR 280
C                   IER =  0  -  NO ERROR                               DMPR 290
C                   IER =  1  -  ITRA IS NOT A TRANSPOSITION VECTOR ON  DMPR 300
C                                1,...,M IF ROWS ARE PERMUTED, 1,...,N  DMPR 310
C                                IF COLUMNS ARE PERMUTED                DMPR 320
C                                                                       DMPR 330
C        REMARKS                                                        DMPR 340
C           (1)  IF IER=-1 THERE IS NO COMPUTATION.                     DMPR 350
C           (2)  IF IER= 1, THEN COMPUTATION HAS BEEN UNSUCCESSFUL DUE  DMPR 360
C                TO ERROR, BUT THE MATRIX A WILL REFLECT THE ROW OR     DMPR 370
C                COLUMN INTERCHANGES PERFORMED BEFORE THE ERROR WAS     DMPR 380
C                DETECTED.                                              DMPR 390
C           (3)  THE MATRIX A IS ASSUMED TO BE STORED COLUMNWISE.       DMPR 400
C                                                                       DMPR 410
C        SUBROUTINES AND SUBPROGRAMS REQUIRED                           DMPR 420
C           NONE                                                        DMPR 430
C                                                                       DMPR 440
C        METHOD                                                         DMPR 450
C           THE ROWS OR COLUMNS ARE PERMUTED ELEMENTWISE, INTERCHANGING DMPR 460
C           ROW OR COLUMN 1 AND ITRA(1),...,ROW OR COLUMN K AND ITRA(K) DMPR 470
C           IN THAT ORDER IF INV=0, AND OTHERWISE INTERCHANGING ROW OR  DMPR 480
C           COLUMN K AND ITRA(K),...,ROW OR COLUMN 1 AND ITRA(1), WHERE DMPR 490
C           K IS M OR N DEPENDING ON WHETHER WE PERMUTE ROWS OR COLUMNS.DMPR 500
C                                                                       DMPR 510
C     ..................................................................DMPR 520
C                                                                       DMPR 530
      SUBROUTINE DMPRC(A,M,N,ITRA,INV,IROCO,IER)                        DMPR 540
C                                                                       DMPR 550
C                                                                       DMPR 560
      DIMENSION A(1),ITRA(1)                                            DMPR 570
      DOUBLE PRECISION A,SAVE                                           DMPR 580
C                                                                       DMPR 590
C        TEST OF DIMENSIONS                                             DMPR 600
      IF(M)14,14,1                                                      DMPR 610
    1 IF(N)14,14,2                                                      DMPR 620
C                                                                       DMPR 630
C        DETERMINE WHICH ARE TO BE PERMUTED-THE ROWS OR THE COLUMNS     DMPR 640
    2 IF(IROCO)3,4,3                                                    DMPR 650
C                                                                       DMPR 660
C        INITIALIZE FOR COLUMN INTERCHANGES                             DMPR 670
    3 MM=M                                                              DMPR 680
      MMM=-1                                                            DMPR 690
      L=M                                                               DMPR 700
      LL=N                                                              DMPR 710
      GO TO 5                                                           DMPR 720
C                                                                       DMPR 730
C        INITIALIZE FOR ROW INTERCHANGES                                DMPR 740
    4 MM=1                                                              DMPR 750
      MMM=M                                                             DMPR 760
      L=N                                                               DMPR 770
      LL=M                                                              DMPR 780
C                                                                       DMPR 790
C        INITIALIZE LOOP OVER ALL ROWS OR COLUMNS                       DMPR 800
    5 IA=1                                                              DMPR 810
      ID=1                                                              DMPR 820
C                                                                       DMPR 830
C        TEST FOR INVERSE OPERATION                                     DMPR 840
      IF(INV)6,7,6                                                      DMPR 850
    6 IA=LL                                                             DMPR 860
      ID=-1                                                             DMPR 870
    7 DO 12 I=1,LL                                                      DMPR 880
      K=ITRA(IA)                                                        DMPR 890
      IF(K-IA)8,12,9                                                    DMPR 900
    8 IF(K)13,13,10                                                     DMPR 910
    9 IF(LL-K)13,10,10                                                  DMPR 920
C                                                                       DMPR 930
C        INITIALIZE ROW OR COLUMN INTERCHANGE                           DMPR 940
   10 IL=IA*MM                                                          DMPR 950
      K=K*MM                                                            DMPR 960
C                                                                       DMPR 970
C        PERFORM ROW OR COLUMN INTERCHANGE                              DMPR 980
      DO 11 J=1,L                                                       DMPR 990
      SAVE=A(IL)                                                        DMPR1000
      A(IL)=A(K)                                                        DMPR1010
      A(K)=SAVE                                                         DMPR1020
      K=K+MMM                                                           DMPR1030
   11 IL=IL+MMM                                                         DMPR1040
C                                                                       DMPR1050
C        ADDRESS NEXT INTERCHANGE STEP                                  DMPR1060
   12 IA=IA+ID                                                          DMPR1070
C                                                                       DMPR1080
C        NORMAL EXIT                                                    DMPR1090
      IER=0                                                             DMPR1100
      RETURN                                                            DMPR1110
C                                                                       DMPR1120
C        ERROR RETURN IN CASE ITRA IS NOT A TRANSPOSITION VECTOR        DMPR1130
   13 IER=1                                                             DMPR1140
      RETURN                                                            DMPR1150
C                                                                       DMPR1160
C        ERROR RETURN IN CASE OF ILLEGAL DIMENSIONS                     DMPR1170
   14 IER=-1                                                            DMPR1180
      RETURN                                                            DMPR1190
      END                                                               DMPR1200
