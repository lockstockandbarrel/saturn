C                                                                       MSTR  10
C     ..................................................................MSTR  20
C                                                                       MSTR  30
C        SUBROUTINE MSTR                                                MSTR  40
C                                                                       MSTR  50
C        PURPOSE                                                        MSTR  60
C           CHANGE STORAGE MODE OF A MATRIX                             MSTR  70
C                                                                       MSTR  80
C        USAGE                                                          MSTR  90
C           CALL MSTR(A,R,N,MSA,MSR)                                    MSTR 100
C                                                                       MSTR 110
C        DESCRIPTION OF PARAMETERS                                      MSTR 120
C           A - NAME OF INPUT MATRIX                                    MSTR 130
C           R - NAME OF OUTPUT MATRIX                                   MSTR 140
C           N - NUMBER OF ROWS AND COLUMNS IN A AND R                   MSTR 150
C           MSA - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         MSTR 160
C                  0 - GENERAL                                          MSTR 170
C                  1 - SYMMETRIC                                        MSTR 180
C                  2 - DIAGONAL                                         MSTR 190
C           MSR - SAME AS MSA EXCEPT FOR MATRIX R                       MSTR 200
C                                                                       MSTR 210
C        REMARKS                                                        MSTR 220
C           MATRIX R CANNOT BE IN THE SAME LOCATION AS MATRIX A         MSTR 230
C           MATRIX A MUST BE A SQUARE MATRIX                            MSTR 240
C                                                                       MSTR 250
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  MSTR 260
C           LOC                                                         MSTR 270
C                                                                       MSTR 280
C        METHOD                                                         MSTR 290
C           MATRIX A IS RESTRUCTURED TO FORM MATRIX R.                  MSTR 300
C            MSA MSR                                                    MSTR 310
C             0   0  MATRIX A IS MOVED TO MATRIX R                      MSTR 320
C             0   1  THE UPPER TRIANGLE ELEMENTS OF A GENERAL MATRIX    MSTR 330
C                    ARE USED TO FORM A SYMMETRIC MATRIX                MSTR 340
C             0   2  THE DIAGONAL ELEMENTS OF A GENERAL MATRIX ARE USED MSTR 350
C                    TO FORM A DIAGONAL MATRIX                          MSTR 360
C             1   0  A SYMMETRIC MATRIX IS EXPANDED TO FORM A GENERAL   MSTR 370
C                    MATRIX                                             MSTR 380
C             1   1  MATRIX A IS MOVED TO MATRIX R                      MSTR 390
C             1   2  THE DIAGONAL ELEMENTS OF A SYMMETRIC MATRIX ARE    MSTR 400
C                    USED TO FORM A DIAGONAL MATRIX                     MSTR 410
C             2   0  A DIAGONAL MATRIX IS EXPANDED BY INSERTING MISSING MSTR 420
C                    ZERO ELEMENTS TO FORM A GENERAL MATRIX             MSTR 430
C             2   1  A DIAGONAL MATRIX IS EXPANDED BY INSERTING MISSING MSTR 440
C                    ZERO ELEMENTS TO FORM A SYMMETRIC MATRIX           MSTR 450
C             2   2  MATRIX A IS MOVED TO MATRIX R                      MSTR 460
C                                                                       MSTR 470
C     ..................................................................MSTR 480
C                                                                       MSTR 490
      SUBROUTINE MSTR(A,R,N,MSA,MSR)                                    MSTR 500
      DIMENSION A(1),R(1)                                               MSTR 510
C                                                                       MSTR 520
      DO 20 I=1,N                                                       MSTR 530
      DO 20 J=1,N                                                       MSTR 540
C                                                                       MSTR 550
C        IF R IS GENERAL, FORM ELEMENT                                  MSTR 560
C                                                                       MSTR 570
      IF(MSR) 5,10,5                                                    MSTR 580
C                                                                       MSTR 590
C        IF IN LOWER TRIANGLE OF SYMMETRIC OR DIAGONAL R, BYPASS        MSTR 600
C                                                                       MSTR 610
    5 IF(I-J) 10,10,20                                                  MSTR 620
   10 CALL LOC(I,J,IR,N,N,MSR)                                          MSTR 630
C                                                                       MSTR 640
C        IF IN UPPER AND OFF DIAGONAL  OF DIAGONAL R, BYPASS            MSTR 650
C                                                                       MSTR 660
      IF(IR) 20,20,15                                                   MSTR 670
C                                                                       MSTR 680
C        OTHERWISE, FORM R(I,J)                                         MSTR 690
C                                                                       MSTR 700
   15 R(IR)=0.0                                                         MSTR 710
      CALL LOC(I,J,IA,N,N,MSA)                                          MSTR 720
C                                                                       MSTR 730
C        IF THERE IS NO A(I,J), LEAVE R(I,J) AT 0.0                     MSTR 740
C                                                                       MSTR 750
      IF(IA) 20,20,18                                                   MSTR 760
   18 R(IR)=A(IA)                                                       MSTR 770
   20 CONTINUE                                                          MSTR 780
      RETURN                                                            MSTR 790
      END                                                               MSTR 800
