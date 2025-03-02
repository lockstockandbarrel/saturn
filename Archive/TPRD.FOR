C                                                                       TPRD  10
C     ..................................................................TPRD  20
C                                                                       TPRD  30
C        SUBROUTINE TPRD                                                TPRD  40
C                                                                       TPRD  50
C        PURPOSE                                                        TPRD  60
C           TRANSPOSE A MATRIX AND POSTMULTIPLY BY ANOTHER TO FORM      TPRD  70
C           A RESULTANT MATRIX                                          TPRD  80
C                                                                       TPRD  90
C        USAGE                                                          TPRD 100
C           CALL TPRD(A,B,R,N,M,MSA,MSB,L)                              TPRD 110
C                                                                       TPRD 120
C        DESCRIPTION OF PARAMETERS                                      TPRD 130
C           A - NAME OF FIRST INPUT MATRIX                              TPRD 140
C           B - NAME OF SECOND INPUT MATRIX                             TPRD 150
C           R - NAME OF OUTPUT MATRIX                                   TPRD 160
C           N - NUMBER OF ROWS IN A AND B                               TPRD 170
C           M - NUMBER OF COLUMNS IN A AND ROWS IN R                    TPRD 180
C           MSA - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         TPRD 190
C                  0 - GENERAL                                          TPRD 200
C                  1 - SYMMETRIC                                        TPRD 210
C                  2 - DIAGONAL                                         TPRD 220
C           MSB - SAME AS MSA EXCEPT FOR MATRIX B                       TPRD 230
C           L - NUMBER OF COLUMNS IN B AND R                            TPRD 240
C                                                                       TPRD 250
C        REMARKS                                                        TPRD 260
C           MATRIX R CANNOT BE IN THE SAME LOCATION AS MATRICES A OR B  TPRD 270
C                                                                       TPRD 280
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  TPRD 290
C           LOC                                                         TPRD 300
C                                                                       TPRD 310
C        METHOD                                                         TPRD 320
C           MATRIX TRANSPOSE OF A IS NOT ACTUALLY CALCULATED. INSTEAD,  TPRD 330
C           ELEMENTS IN MATRIX A ARE TAKEN COLUMNWISE RATHER THAN       TPRD 340
C           ROWWISE FOR MULTIPLICATION BY MATRIX B.                     TPRD 350
C           THE FOLLOWING TABLE SHOWS THE STORAGE MODE OF THE OUTPUT    TPRD 360
C           MATRIX FOR ALL COMBINATIONS OF INPUT MATRICES               TPRD 370
C                         A                B                R           TPRD 380
C                      GENERAL          GENERAL          GENERAL        TPRD 390
C                      GENERAL          SYMMETRIC        GENERAL        TPRD 400
C                      GENERAL          DIAGONAL         GENERAL        TPRD 410
C                      SYMMETRIC        GENERAL          GENERAL        TPRD 420
C                      SYMMETRIC        SYMMETRIC        GENERAL        TPRD 430
C                      SYMMETRIC        DIAGONAL         GENERAL        TPRD 440
C                      DIAGONAL         GENERAL          GENERAL        TPRD 450
C                      DIAGONAL         SYMMETRIC        GENERAL        TPRD 460
C                      DIAGONAL         DIAGONAL         DIAGONAL       TPRD 470
C                                                                       TPRD 480
C     ..................................................................TPRD 490
C                                                                       TPRD 500
      SUBROUTINE TPRD(A,B,R,N,M,MSA,MSB,L)                              TPRD 510
      DIMENSION A(1),B(1),R(1)                                          TPRD 520
C                                                                       TPRD 530
C        SPECIAL CASE FOR DIAGONAL BY DIAGONAL                          TPRD 540
C                                                                       TPRD 550
      MS=MSA*10+MSB                                                     TPRD 560
      IF(MS-22) 30,10,30                                                TPRD 570
   10 DO 20 I=1,N                                                       TPRD 580
   20 R(I)=A(I)*B(I)                                                    TPRD 590
      RETURN                                                            TPRD 600
C                                                                       TPRD 610
C        MULTIPLY TRANSPOSE OF A BY B                                   TPRD 620
C                                                                       TPRD 630
   30 IR=1                                                              TPRD 640
      DO 90 K=1,L                                                       TPRD 650
      DO 90 J=1,M                                                       TPRD 660
      R(IR)=0.0                                                         TPRD 670
      DO 80 I=1,N                                                       TPRD 680
      IF(MS) 40,60,40                                                   TPRD 690
   40 CALL LOC(I,J,IA,N,M,MSA)                                          TPRD 700
      CALL LOC(I,K,IB,N,L,MSB)                                          TPRD 710
      IF(IA) 50,80,50                                                   TPRD 720
   50 IF(IB) 70,80,70                                                   TPRD 730
   60 IA=N*(J-1)+I                                                      TPRD 740
      IB=N*(K-1)+I                                                      TPRD 750
   70 R(IR)=R(IR)+A(IA)*B(IB)                                           TPRD 760
   80 CONTINUE                                                          TPRD 770
   90 IR=IR+1                                                           TPRD 780
      RETURN                                                            TPRD 790
      END                                                               TPRD 800
