C                                                                       QTES  10
C     ..................................................................QTES  20
C                                                                       QTES  30
C        SUBROUTINE QTEST                                               QTES  40
C                                                                       QTES  50
C        PURPOSE                                                        QTES  60
C           TEST WHETHER THREE OR MORE MATCHED GROUPS OF DICHOTOMOUS    QTES  70
C           DATA DIFFER SIGNIFICANTLY BY THE COCHRAN Q-TEST             QTES  80
C                                                                       QTES  90
C        USAGE                                                          QTES 100
C           CALL QTEST(A,N,M,Q,NDF)                                     QTES 110
C                                                                       QTES 120
C        DESCRIPTION OF PARAMETERS                                      QTES 130
C           A   - INPUT MATRIX, N BY M, OF DICHOTOMOUS DATA (0 AND 1)   QTES 140
C           N   - NUMBER OF SETS IN EACH GROUP                          QTES 150
C           M   - NUMBER OF GROUPS                                      QTES 160
C           Q   - COCHRAN Q STATISTIC (OUTPUT)                          QTES 170
C           NDF - NUMBER OF DEGREES OF FREEDOM (OUTPUT)                 QTES 180
C                                                                       QTES 190
C        REMARKS                                                        QTES 200
C           M MUST BE THREE OR GREATER                                  QTES 210
C                                                                       QTES 220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QTES 230
C           NONE                                                        QTES 240
C                                                                       QTES 250
C        METHOD                                                         QTES 260
C           DESCRIBED IN S. SIEGEL, 'NONPARAMETRIC STATISTICS FOR THE   QTES 270
C           BEHAVIORAL SCIENCES', MCGRAW-HILL, NEW YORK, 1956,          QTES 280
C           CHAPTER 7                                                   QTES 290
C                                                                       QTES 300
C     ..................................................................QTES 310
C                                                                       QTES 320
      SUBROUTINE QTEST(A,N,M,Q,NDF)                                     QTES 330
      DIMENSION A(1)                                                    QTES 340
C                                                                       QTES 350
C        COMPUTE SUM OF SQUARES OF ROW TOTALS, RSQ, AND GRAND TOTAL OF  QTES 360
C        ALL ELEMENTS, GD                                               QTES 370
C                                                                       QTES 380
      RSQ=0.0                                                           QTES 390
      GD=0.0                                                            QTES 400
      DO 20 I=1,N                                                       QTES 410
      TR=0.0                                                            QTES 420
      IJ=I-N                                                            QTES 430
      DO 10 J=1,M                                                       QTES 440
      IJ=IJ+N                                                           QTES 450
   10 TR=TR+A(IJ)                                                       QTES 460
      GD=GD+TR                                                          QTES 470
   20 RSQ=RSQ+TR*TR                                                     QTES 480
C                                                                       QTES 490
C        COMPUTE SUM OF SQUARES OF COLUMN TOTALS, CSQ                   QTES 500
C                                                                       QTES 510
      CSQ=0.0                                                           QTES 520
      IJ=0                                                              QTES 530
      DO 40 J=1,M                                                       QTES 540
      TC=0.0                                                            QTES 550
      DO 30 I=1,N                                                       QTES 560
      IJ=IJ+1                                                           QTES 570
   30 TC=TC+A(IJ)                                                       QTES 580
   40 CSQ=CSQ+TC*TC                                                     QTES 590
C                                                                       QTES 600
C        COMPUTE COCHRAN Q TEST VALUE                                   QTES 610
C                                                                       QTES 620
      FM=M                                                              QTES 630
      Q=(FM-1.0)*(FM*CSQ-GD*GD)/(FM*GD-RSQ)                             QTES 640
C                                                                       QTES 650
C        FIND DEGREES OF FREEDOM                                        QTES 660
C                                                                       QTES 670
      NDF=M-1                                                           QTES 680
      RETURN                                                            QTES 690
      END                                                               QTES 700
