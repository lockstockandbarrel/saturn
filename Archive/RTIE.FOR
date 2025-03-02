C                                                                       RTIE  10
C     ..................................................................RTIE  20
C                                                                       RTIE  30
C        SUBROUTINE RTIE                                                RTIE  40
C                                                                       RTIE  50
C        PURPOSE                                                        RTIE  60
C           ADJOIN TWO MATRICES WITH SAME COLUMN DIMENSION TO FORM ONE  RTIE  70
C           RESULTANT MATRIX (SEE METHOD)                               RTIE  80
C                                                                       RTIE  90
C        USAGE                                                          RTIE 100
C           CALL RTIE(A,B,R,N,M,MSA,MSB,L)                              RTIE 110
C                                                                       RTIE 120
C        DESCRIPTION OF PARAMETERS                                      RTIE 130
C           A - NAME OF FIRST INPUT MATRIX                              RTIE 140
C           B - NAME OF SECOND INPUT MATRIX                             RTIE 150
C           R - NAME OF OUTPUT MATRIX                                   RTIE 160
C           N - NUMBER OF ROWS IN A                                     RTIE 170
C           M - NUMBER OF COLUMNS IN A,B,R                              RTIE 180
C           MSA - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         RTIE 190
C                  0 - GENERAL                                          RTIE 200
C                  1 - SYMMETRIC                                        RTIE 210
C                  2 - DIAGONAL                                         RTIE 220
C           MSB - SAME AS MSA EXCEPT FOR MATRIX B                       RTIE 230
C           L - NUMBER OF ROWS IN B                                     RTIE 240
C                                                                       RTIE 250
C        REMARKS                                                        RTIE 260
C           MATRIX R CANNOT BE IN THE SAME LOCATION AS MATRICES A OR B  RTIE 270
C           MATRIX R IS ALWAYS A GENERAL MATRIX                         RTIE 280
C           MATRIX A MUST HAVE THE SAME NUMBER OF COLUMNS AS MATRIX B   RTIE 290
C                                                                       RTIE 300
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  RTIE 310
C           LOC                                                         RTIE 320
C                                                                       RTIE 330
C        METHOD                                                         RTIE 340
C           MATRIX B IS ATTACHED TO THE BOTTOM OF MATRIX A .            RTIE 350
C           THE RESULTANT MATRIX R CONTAINS N+L ROWS AND M COLUMNS.     RTIE 360
C                                                                       RTIE 370
C     ..................................................................RTIE 380
C                                                                       RTIE 390
      SUBROUTINE RTIE(A,B,R,N,M,MSA,MSB,L)                              RTIE 400
      DIMENSION A(1),B(1),R(1)                                          RTIE 410
C                                                                       RTIE 420
      NN=N                                                              RTIE 430
      IR=0                                                              RTIE 440
      NX=NN                                                             RTIE 450
      MSX=MSA                                                           RTIE 460
      DO 9 J=1,M                                                        RTIE 470
      DO 8 II=1,2                                                       RTIE 480
      DO 7 I=1,NN                                                       RTIE 490
      IR=IR+1                                                           RTIE 500
      R(IR)=0.0                                                         RTIE 510
C                                                                       RTIE 520
C        LOCATE ELEMENT FOR ANY MATRIX STORAGE MODE                     RTIE 530
C                                                                       RTIE 540
      CALL LOC(I,J,IJ,NN,M,MSX)                                         RTIE 550
C                                                                       RTIE 560
C        TEST FOR ZERO ELEMENT IN DIAGONAL MATRIX                       RTIE 570
C                                                                       RTIE 580
      IF(IJ) 2,7,2                                                      RTIE 590
C                                                                       RTIE 600
C        MOVE ELEMENT TO MATRIX R                                       RTIE 610
C                                                                       RTIE 620
    2 GO TO(3,4),II                                                     RTIE 630
    3 R(IR)=A(IJ)                                                       RTIE 640
      GO TO 7                                                           RTIE 650
    4 R(IR)=B(IJ)                                                       RTIE 660
    7 CONTINUE                                                          RTIE 670
C                                                                       RTIE 680
C        REPEAT ABOVE FOR MATRIX B                                      RTIE 690
C                                                                       RTIE 700
      MSX=MSB                                                           RTIE 710
    8 NN=L                                                              RTIE 720
C                                                                       RTIE 730
C        RESET FOR NEXT COLUMN                                          RTIE 740
C                                                                       RTIE 750
      MSX=MSA                                                           RTIE 760
    9 NN=NX                                                             RTIE 770
      RETURN                                                            RTIE 780
      END                                                               RTIE 790
