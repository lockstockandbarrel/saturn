C                                                                       CHIS  10
C     ..................................................................CHIS  20
C                                                                       CHIS  30
C        SUBROUTINE CHISQ                                               CHIS  40
C                                                                       CHIS  50
C        PURPOSE                                                        CHIS  60
C           COMPUTE CHI-SQUARE FROM A CONTINGENCY TABLE                 CHIS  70
C                                                                       CHIS  80
C        USAGE                                                          CHIS  90
C           CALL CHISQ(A,N,M,CS,NDF,IERR,TR,TC)                         CHIS 100
C                                                                       CHIS 110
C        DESCRIPTION OF PARAMETERS                                      CHIS 120
C           A    - INPUT MATRIX, N BY M, CONTAINING CONTINGENCY TABLE   CHIS 130
C           N    - NUMBER OF ROWS IN A                                  CHIS 140
C           M    - NUMBER OF COLUMNS IN A                               CHIS 150
C           CS   - CHI-SQUARE (OUTPUT)                                  CHIS 160
C           NDF  - NUMBER OF DEGREES OF FREEDOM (OUTPUT)                CHIS 170
C           IERR - ERROR CODE (OUTPUT)                                  CHIS 180
C                    0 - NORMAL CASE                                    CHIS 190
C                    1 - EXPECTED VALUE IS LESS THAN 1.0 IN ONE OR      CHIS 200
C                        MORE CELLS                                     CHIS 210
C                    3 - NUMBER OF DEGREES OF FREEDOM IS ZERO           CHIS 220
C           TR   - WORK VECTOR OF LENGTH N                              CHIS 230
C           TC   - WORK VECTOR OF LENGTH M                              CHIS 240
C                                                                       CHIS 250
C        REMARKS                                                        CHIS 260
C           IF ONE OR MORE CELLS CONTAIN AN EXPECTED VALUE (I.E.,       CHIS 270
C           THEORETICAL VALUE) LESS THAN 1.0, CHI-SQUARE WILL BE        CHIS 280
C           COMPUTED, BUT ERROR CODE WILL BE SET TO 1.                  CHIS 290
C           SEE REFERENCE GIVEN BELOW.                                  CHIS 300
C           CHI-SQUARE IS SET TO ZERO IF EITHER N OR M IS ONE (ERROR    CHIS 310
C           CODE 3).                                                    CHIS 320
C                                                                       CHIS 330
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  CHIS 340
C           NONE                                                        CHIS 350
C                                                                       CHIS 360
C        METHOD                                                         CHIS 370
C           DESCRIBED IN S. SIEGEL, 'NONPARAMETRIC STATISTICS FOR THE   CHIS 380
C           BEHAVIORAL SCIENCES', MCGRAW-HILL, NEW YORK, 1956,          CHIS 390
C           CHAPTER 6 AND CHAPTER 8.                                    CHIS 400
C                                                                       CHIS 410
C     ..................................................................CHIS 420
C                                                                       CHIS 430
      SUBROUTINE CHISQ(A,N,M,CS,NDF,IERR,TR,TC)                         CHIS 440
      DIMENSION A(1),TR(1),TC(1)                                        CHIS 450
C                                                                       CHIS 460
      NM=N*M                                                            CHIS 470
      IERR=0                                                            CHIS 480
      CS=0.0                                                            CHIS 490
C                                                                       CHIS 500
C        FIND DEGREES OF FREEDOM                                        CHIS 510
C                                                                       CHIS 520
      NDF=(N-1)*(M-1)                                                   CHIS 530
      IF(NDF) 5,5,10                                                    CHIS 540
    5 IERR=3                                                            CHIS 550
      RETURN                                                            CHIS 560
C                                                                       CHIS 570
C        COMPUTE TOTALS OF ROWS                                         CHIS 580
C                                                                       CHIS 590
   10 DO 90 I=1,N                                                       CHIS 600
      TR(I)=0.0                                                         CHIS 610
      IJ=I-N                                                            CHIS 620
      DO 90 J=1,M                                                       CHIS 630
      IJ=IJ+N                                                           CHIS 640
   90 TR(I)=TR(I)+A(IJ)                                                 CHIS 650
C                                                                       CHIS 660
C        COMPUTE TOTALS OF COLUMNS                                      CHIS 670
C                                                                       CHIS 680
      IJ=0                                                              CHIS 690
      DO 100 J=1,M                                                      CHIS 700
      TC(J)=0.0                                                         CHIS 710
      DO 100 I=1,N                                                      CHIS 720
      IJ=IJ+1                                                           CHIS 730
  100 TC(J)=TC(J)+A(IJ)                                                 CHIS 740
C                                                                       CHIS 750
C        COMPUTE GRAND TOTAL                                            CHIS 760
C                                                                       CHIS 770
      GT=0.0                                                            CHIS 780
      DO 110 I=1,N                                                      CHIS 790
  110 GT=GT+TR(I)                                                       CHIS 800
C                                                                       CHIS 810
C        COMPUTE CHI SQUARE FOR 2 BY 2 TABLE (SPECIAL CASE)             CHIS 820
C                                                                       CHIS 830
      IF(NM-4) 130,120,130                                              CHIS 840
  120 CS=GT*(ABS(A(1)*A(4)-A(2)*A(3))-GT/2.0)**2  /(TC(1)*TC(2)*TR(1)   CHIS 850
     1*TR(2))                                                           CHIS 860
      RETURN                                                            CHIS 870
C                                                                       CHIS 880
C        COMPUTE CHI SQUARE FOR OTHER CONTINGENCY TABLES                CHIS 890
C                                                                       CHIS 900
  130 IJ=0                                                              CHIS 910
      DO 140 J=1,M                                                      CHIS 920
      DO 140 I=1,N                                                      CHIS 930
      IJ=IJ+1                                                           CHIS 940
      E=TR(I)*TC(J)/GT                                                  CHIS 950
      IF(E-1.0) 135, 140, 140                                           CHIS 960
  135 IERR=1                                                            CHIS 970
  140 CS=CS+(A(IJ)-E)*(A(IJ)-E)/E                                       CHIS 980
      RETURN                                                            CHIS 990
      END                                                               CHIS1000
