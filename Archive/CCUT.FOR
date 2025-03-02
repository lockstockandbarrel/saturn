C                                                                       CCUT  10
C     ..................................................................CCUT  20
C                                                                       CCUT  30
C        SUBROUTINE CCUT                                                CCUT  40
C                                                                       CCUT  50
C        PURPOSE                                                        CCUT  60
C           PARTITION A MATRIX BETWEEN SPECIFIED COLUMNS TO FORM TWO    CCUT  70
C           RESULTANT MATRICES                                          CCUT  80
C                                                                       CCUT  90
C        USAGE                                                          CCUT 100
C           CALL CCUT (A,L,R,S,N,M,MS)                                  CCUT 110
C                                                                       CCUT 120
C        DESCRIPTION OF PARAMETERS                                      CCUT 130
C           A - NAME OF INPUT MATRIX                                    CCUT 140
C           L - COLUMN OF A TO THE LEFT OF WHICH PARTITIONING TAKES     CCUT 150
C               PLACE                                                   CCUT 160
C           R - NAME OF MATRIX TO BE FORMED FROM LEFT PORTION OF A      CCUT 170
C           S - NAME OF MATRIX TO BE FORMED FROM RIGHT PORTION OF A     CCUT 180
C           N - NUMBER OF ROWS IN A                                     CCUT 190
C           M - NUMBER OF COLUMNS IN A                                  CCUT 200
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         CCUT 210
C                  0 - GENERAL                                          CCUT 220
C                  1 - SYMMETRIC                                        CCUT 230
C                  2 - DIAGONAL                                         CCUT 240
C                                                                       CCUT 250
C        REMARKS                                                        CCUT 260
C           MATRIX R CANNOT BE IN SAME LOCATION AS MATRIX A             CCUT 270
C           MATRIX S CANNOT BE IN SAME LOCATION AS MATRIX A             CCUT 280
C           MATRIX R CANNOT BE IN SAME LOCATION AS MATRIX S             CCUT 290
C           MATRIX R AND MATRIX S ARE ALWAYS GENERAL MATRICES           CCUT 300
C                                                                       CCUT 310
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  CCUT 320
C           LOC                                                         CCUT 330
C                                                                       CCUT 340
C        METHOD                                                         CCUT 350
C           ELEMENTS OF MATRIX A TO THE LEFT OF COLUMN L ARE MOVED TO   CCUT 360
C           FORM MATRIX R OF N ROWS AND L-1 COLUMNS. ELEMENTS OF        CCUT 370
C           MATRIX A IN COLUMN L AND TO THE RIGHT OF L ARE MOVED TO FORMCCUT 380
C           MATRIX S OF N ROWS AND M-L+1 COLUMNS.                       CCUT 390
C                                                                       CCUT 400
C     ..................................................................CCUT 410
C                                                                       CCUT 420
      SUBROUTINE CCUT(A,L,R,S,N,M,MS)                                   CCUT 430
      DIMENSION A(1),R(1),S(1)                                          CCUT 440
C                                                                       CCUT 450
      IR=0                                                              CCUT 460
      IS=0                                                              CCUT 470
      DO 70 J=1,M                                                       CCUT 480
      DO 70 I=1,N                                                       CCUT 490
C                                                                       CCUT 500
C        FIND LOCATION IN OUTPUT MATRIX AND SET TO ZERO                 CCUT 510
C                                                                       CCUT 520
      IF(J-L) 20,10,10                                                  CCUT 530
   10 IS=IS+1                                                           CCUT 540
      S(IS)=0.0                                                         CCUT 550
      GO TO 30                                                          CCUT 560
   20 IR=IR+1                                                           CCUT 570
      R(IR)=0.0                                                         CCUT 580
C                                                                       CCUT 590
C        LOCATE ELEMENT FOR ANY MATRIX STORAGE MODE                     CCUT 600
C                                                                       CCUT 610
   30 CALL LOC(I,J,IJ,N,M,MS)                                           CCUT 620
C                                                                       CCUT 630
C        TEST FOR ZERO ELEMENT IN DIAGONAL MATRIX                       CCUT 640
C                                                                       CCUT 650
      IF(IJ) 40,70,40                                                   CCUT 660
C                                                                       CCUT 670
C        DETERMINE WHETHER RIGHT OR LEFT OF L                           CCUT 680
C                                                                       CCUT 690
   40 IF(J-L) 60,50,50                                                  CCUT 700
   50 S(IS)=A(IJ)                                                       CCUT 710
      GO TO 70                                                          CCUT 720
   60 R(IR)=A(IJ)                                                       CCUT 730
   70 CONTINUE                                                          CCUT 740
      RETURN                                                            CCUT 750
      END                                                               CCUT 760
