C                                                                       CTIE  10
C     ..................................................................CTIE  20
C                                                                       CTIE  30
C        SUBROUTINE CTIE                                                CTIE  40
C                                                                       CTIE  50
C        PURPOSE                                                        CTIE  60
C           ADJOIN TWO MATRICES WITH SAME ROW DIMENSION TO FORM ONE     CTIE  70
C           RESULTANT MATRIX (SEE METHOD)                               CTIE  80
C                                                                       CTIE  90
C        USAGE                                                          CTIE 100
C           CALL CTIE(A,B,R,N,M,MSA,MSB,L)                              CTIE 110
C                                                                       CTIE 120
C        DESCRIPTION OF PARAMETERS                                      CTIE 130
C           A - NAME OF FIRST INPUT MATRIX                              CTIE 140
C           B - NAME OF SECOND INPUT MATRIX                             CTIE 150
C           R - NAME OF OUTPUT MATRIX                                   CTIE 160
C           N - NUMBER OF ROWS IN A,B,R                                 CTIE 170
C           M - NUMBER OF COLUMNS IN A                                  CTIE 180
C           MSA - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         CTIE 190
C                  0 - GENERAL                                          CTIE 200
C                  1 - SYMMETRIC                                        CTIE 210
C                  2 - DIAGONAL                                         CTIE 220
C           MSB - SAME AS MSA EXCEPT FOR MATRIX B                       CTIE 230
C           L - NUMBER OF COLUMNS IN B                                  CTIE 240
C                                                                       CTIE 250
C        REMARKS                                                        CTIE 260
C           MATRIX R CANNOT BE IN THE SAME LOCATION AS MATRICES A OR B  CTIE 270
C           MATRIX R IS ALWAYS A GENERAL MATRIX                         CTIE 280
C           MATRIX A MUST HAVE THE SAME NUMBER OF ROWS AS MATRIX B      CTIE 290
C                                                                       CTIE 300
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  CTIE 310
C           LOC                                                         CTIE 320
C                                                                       CTIE 330
C        METHOD                                                         CTIE 340
C           MATRIX B IS ATTACHED TO THE RIGHT OF MATRIX A .             CTIE 350
C           THE RESULTANT MATRIX R CONTAINS N ROWS AND M+L COLUMNS      CTIE 360
C                                                                       CTIE 370
C     ..................................................................CTIE 380
C                                                                       CTIE 390
      SUBROUTINE CTIE(A,B,R,N,M,MSA,MSB,L)                              CTIE 400
      DIMENSION A(1),B(1),R(1)                                          CTIE 410
C                                                                       CTIE 420
      MM=M                                                              CTIE 430
      IR=0                                                              CTIE 440
      MSX=MSA                                                           CTIE 450
      DO 6 JJ=1,2                                                       CTIE 460
      DO 5 J=1,MM                                                       CTIE 470
      DO 5 I=1,N                                                        CTIE 480
      IR=IR+1                                                           CTIE 490
      R(IR)=0.0                                                         CTIE 500
C                                                                       CTIE 510
C        LOCATE ELEMENT FOR ANY MATRIX STORAGE MODE                     CTIE 520
C                                                                       CTIE 530
      CALL LOC(I,J,IJ,N,MM,MSX)                                         CTIE 540
C                                                                       CTIE 550
C        TEST FOR ZERO ELEMENT IN DIAGONAL MATRIX                       CTIE 560
C                                                                       CTIE 570
      IF(IJ) 2,5,2                                                      CTIE 580
C                                                                       CTIE 590
C        MOVE ELEMENT TO MATRIX R                                       CTIE 600
C                                                                       CTIE 610
    2 GO TO(3,4),JJ                                                     CTIE 620
    3 R(IR)=A(IJ)                                                       CTIE 630
      GO TO 5                                                           CTIE 640
    4 R(IR)=B(IJ)                                                       CTIE 650
    5 CONTINUE                                                          CTIE 660
C                                                                       CTIE 670
C        REPEAT ABOVE FOR MATRIX B                                      CTIE 680
C                                                                       CTIE 690
      MSX=MSB                                                           CTIE 700
      MM=L                                                              CTIE 710
    6 CONTINUE                                                          CTIE 720
      RETURN                                                            CTIE 730
      END                                                               CTIE 740
