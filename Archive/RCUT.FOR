C                                                                       RCUT  10
C     ..................................................................RCUT  20
C                                                                       RCUT  30
C        SUBROUTINE RCUT                                                RCUT  40
C                                                                       RCUT  50
C        PURPOSE                                                        RCUT  60
C           PARTITION A MATRIX BETWEEN SPECIFIED ROWS TO FORM TWO       RCUT  70
C           RESULTANT MATRICES                                          RCUT  80
C                                                                       RCUT  90
C        USAGE                                                          RCUT 100
C           CALL RCUT (A,L,R,S,N,M,MS)                                  RCUT 110
C                                                                       RCUT 120
C        DESCRIPTION OF PARAMETERS                                      RCUT 130
C           A - NAME OF INPUT MATRIX                                    RCUT 140
C           L - ROW OF A ABOVE WHICH PARTITIONING TAKES PLACE           RCUT 150
C           R - NAME OF MATRIX TO BE FORMED FROM UPPER PORTION OF A     RCUT 160
C           S - NAME OF MATRIX TO BE FORMED FROM LOWER PORTION OF A     RCUT 170
C           N - NUMBER OF ROWS IN A                                     RCUT 180
C           M - NUMBER OF COLUMNS IN A                                  RCUT 190
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         RCUT 200
C                  0 - GENERAL                                          RCUT 210
C                  1 - SYMMETRIC                                        RCUT 220
C                  2 - DIAGONAL                                         RCUT 230
C                                                                       RCUT 240
C        REMARKS                                                        RCUT 250
C           MATRIX R CANNOT BE IN SAME LOCATION AS MATRIX A             RCUT 260
C           MATRIX S CANNOT BE IN SAME LOCATION AS MATRIX A             RCUT 270
C           MATRIX R CANNOT BE IN SAME LOCATION AS MATRIX S             RCUT 280
C           MATRIX R AND MATRIX S ARE ALWAYS GENERAL MATRICES           RCUT 290
C                                                                       RCUT 300
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  RCUT 310
C           LOC                                                         RCUT 320
C                                                                       RCUT 330
C        METHOD                                                         RCUT 340
C           ELEMENTS OF MATRIX A ABOVE ROW L ARE MOVED TO FORM MATRIX R RCUT 350
C           OF L-1 ROWS AND M COLUMNS. ELEMENTS OF MATRIX A IN ROW L    RCUT 360
C           AND BELOW ARE MOVED TO FORM MATRIX S OF N-L+1 ROWS AND M    RCUT 370
C           COLUMNS                                                     RCUT 380
C                                                                       RCUT 390
C     ..................................................................RCUT 400
C                                                                       RCUT 410
      SUBROUTINE RCUT(A,L,R,S,N,M,MS)                                   RCUT 420
      DIMENSION A(1),R(1),S(1)                                          RCUT 430
C                                                                       RCUT 440
      IR=0                                                              RCUT 450
      IS=0                                                              RCUT 460
      DO 70 J=1,M                                                       RCUT 470
      DO 70 I=1,N                                                       RCUT 480
C                                                                       RCUT 490
C        FIND LOCATION IN OUTPUT MATRIX AND SET TO ZERO                 RCUT 500
C                                                                       RCUT 510
      IF(I-L) 20,10,10                                                  RCUT 520
   10 IS=IS+1                                                           RCUT 530
      S(IS)=0.0                                                         RCUT 540
      GO TO 30                                                          RCUT 550
   20 IR=IR+1                                                           RCUT 560
      R(IR)=0.0                                                         RCUT 570
C                                                                       RCUT 580
C        LOCATE ELEMENT FOR ANY MATRIX STORAGE MODE                     RCUT 590
C                                                                       RCUT 600
   30 CALL LOC(I,J,IJ,N,M,MS)                                           RCUT 610
C                                                                       RCUT 620
C        TEST FOR ZERO ELEMENT IN DIAGONAL MATRIX                       RCUT 630
C                                                                       RCUT 640
      IF(IJ) 40,70,40                                                   RCUT 650
C                                                                       RCUT 660
C        DETERMINE WHETHER ABOVE OR BELOW L                             RCUT 670
C                                                                       RCUT 680
   40 IF(I-L) 60,50,50                                                  RCUT 690
   50 S(IS)=A(IJ)                                                       RCUT 700
      GO TO 70                                                          RCUT 710
   60 R(IR)=A(IJ)                                                       RCUT 720
   70 CONTINUE                                                          RCUT 730
      RETURN                                                            RCUT 740
      END                                                               RCUT 750
