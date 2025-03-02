C                                                                       MPRD  10
C     ..................................................................MPRD  20
C                                                                       MPRD  30
C        SUBROUTINE MPRD                                                MPRD  40
C                                                                       MPRD  50
C        PURPOSE                                                        MPRD  60
C           MULTIPLY TWO MATRICES TO FORM A RESULTANT MATRIX            MPRD  70
C                                                                       MPRD  80
C        USAGE                                                          MPRD  90
C           CALL MPRD(A,B,R,N,M,MSA,MSB,L)                              MPRD 100
C                                                                       MPRD 110
C        DESCRIPTION OF PARAMETERS                                      MPRD 120
C           A - NAME OF FIRST INPUT MATRIX                              MPRD 130
C           B - NAME OF SECOND INPUT MATRIX                             MPRD 140
C           R - NAME OF OUTPUT MATRIX                                   MPRD 150
C           N - NUMBER OF ROWS IN A AND R                               MPRD 160
C           M - NUMBER OF COLUMNS IN A AND ROWS IN B                    MPRD 170
C           MSA - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         MPRD 180
C                  0 - GENERAL                                          MPRD 190
C                  1 - SYMMETRIC                                        MPRD 200
C                  2 - DIAGONAL                                         MPRD 210
C           MSB - SAME AS MSA EXCEPT FOR MATRIX B                       MPRD 220
C           L - NUMBER OF COLUMNS IN B AND R                            MPRD 230
C                                                                       MPRD 240
C        REMARKS                                                        MPRD 250
C           MATRIX R CANNOT BE IN THE SAME LOCATION AS MATRICES A OR B  MPRD 260
C           NUMBER OF COLUMNS OF MATRIX A MUST BE EQUAL TO NUMBER OF ROWMPRD 270
C           OF MATRIX B                                                 MPRD 280
C                                                                       MPRD 290
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  MPRD 300
C           LOC                                                         MPRD 310
C                                                                       MPRD 320
C        METHOD                                                         MPRD 330
C           THE M BY L MATRIX B IS PREMULTIPLIED BY THE N BY M MATRIX A MPRD 340
C           AND THE RESULT IS STORED IN THE N BY L MATRIX R. THIS IS A  MPRD 350
C           ROW INTO COLUMN PRODUCT.                                    MPRD 360
C           THE FOLLOWING TABLE SHOWS THE STORAGE MODE OF THE OUTPUT    MPRD 370
C           MATRIX FOR ALL COMBINATIONS OF INPUT MATRICES               MPRD 380
C                         A                B                R           MPRD 390
C                      GENERAL          GENERAL          GENERAL        MPRD 400
C                      GENERAL          SYMMETRIC        GENERAL        MPRD 410
C                      GENERAL          DIAGONAL         GENERAL        MPRD 420
C                      SYMMETRIC        GENERAL          GENERAL        MPRD 430
C                      SYMMETRIC        SYMMETRIC        GENERAL        MPRD 440
C                      SYMMETRIC        DIAGONAL         GENERAL        MPRD 450
C                      DIAGONAL         GENERAL          GENERAL        MPRD 460
C                      DIAGONAL         SYMMETRIC        GENERAL        MPRD 470
C                      DIAGONAL         DIAGONAL         DIAGONAL       MPRD 480
C                                                                       MPRD 490
C     ..................................................................MPRD 500
C                                                                       MPRD 510
      SUBROUTINE MPRD(A,B,R,N,M,MSA,MSB,L)                              MPRD 520
      DIMENSION A(1),B(1),R(1)                                          MPRD 530
C                                                                       MPRD 540
C        SPECIAL CASE FOR DIAGONAL BY DIAGONAL                          MPRD 550
C                                                                       MPRD 560
      MS=MSA*10+MSB                                                     MPRD 570
      IF(MS-22) 30,10,30                                                MPRD 580
   10 DO 20 I=1,N                                                       MPRD 590
   20 R(I)=A(I)*B(I)                                                    MPRD 600
      RETURN                                                            MPRD 610
C                                                                       MPRD 620
C        ALL OTHER CASES                                                MPRD 630
C                                                                       MPRD 640
   30 IR=1                                                              MPRD 650
      DO 90 K=1,L                                                       MPRD 660
      DO 90 J=1,N                                                       MPRD 670
      R(IR)=0                                                           MPRD 680
      DO 80 I=1,M                                                       MPRD 690
      IF(MS) 40,60,40                                                   MPRD 700
   40 CALL LOC(J,I,IA,N,M,MSA)                                          MPRD 710
      CALL LOC(I,K,IB,M,L,MSB)                                          MPRD 720
      IF(IA) 50,80,50                                                   MPRD 730
   50 IF(IB) 70,80,70                                                   MPRD 740
   60 IA=N*(I-1)+J                                                      MPRD 750
      IB=M*(K-1)+I                                                      MPRD 760
   70 R(IR)=R(IR)+A(IA)*B(IB)                                           MPRD 770
   80 CONTINUE                                                          MPRD 780
   90 IR=IR+1                                                           MPRD 790
      RETURN                                                            MPRD 800
      END                                                               MPRD 810
