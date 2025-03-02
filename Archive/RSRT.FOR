C                                                                       RSRT  10
C     ..................................................................RSRT  20
C                                                                       RSRT  30
C        SUBROUTINE RSRT                                                RSRT  40
C                                                                       RSRT  50
C        PURPOSE                                                        RSRT  60
C           SORT ROWS OF A MATRIX                                       RSRT  70
C                                                                       RSRT  80
C        USAGE                                                          RSRT  90
C           CALL RSRT(A,B,R,N,M,MS)                                     RSRT 100
C                                                                       RSRT 110
C        DESCRIPTION OF PARAMETERS                                      RSRT 120
C           A - NAME OF INPUT MATRIX TO BE SORTED                       RSRT 130
C           B - NAME OF INPUT VECTOR WHICH CONTAINS SORTING KEY         RSRT 140
C           R - NAME OF SORTED OUTPUT MATRIX                            RSRT 150
C           N - NUMBER OF ROWS IN A AND R AND LENGTH OF B               RSRT 160
C           M - NUMBER OF COLUMNS IN A AND R                            RSRT 170
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         RSRT 180
C                  0 - GENERAL                                          RSRT 190
C                  1 - SYMMETRIC                                        RSRT 200
C                  2 - DIAGONAL                                         RSRT 210
C                                                                       RSRT 220
C        REMARKS                                                        RSRT 230
C           MATRIX R CANNOT BE IN THE SAME LOCATION AS MATRIX A         RSRT 240
C           MATRIX R IS ALWAYS A GENERAL MATRIX                         RSRT 250
C           N MUST BE GREATER THAN ONE.                                 RSRT 260
C                                                                       RSRT 270
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  RSRT 280
C           LOC                                                         RSRT 290
C                                                                       RSRT 300
C        METHOD                                                         RSRT 310
C           ROWS OF INPUT MATRIX A ARE SORTED TO FORM OUTPUT MATRIX R.  RSRT 320
C           THE SORTED ROW SEQUENCE IS DETERMINED BY THE VALUES OF      RSRT 330
C           ELEMENTS IN COLUMN VECTOR B. THE LOWEST VALUED ELEMENT IN   RSRT 340
C           B WILL CAUSE THE CORRESPONDING ROW OF A TO BE PLACED IN THE RSRT 350
C           FIRST ROW OF R. THE HIGHEST VALUED ELEMENT OF B WILL CAUSE  RSRT 360
C           THE CORRESPONDING ROW OF A TO BE PLACED IN THE LAST ROW OF  RSRT 370
C           R. IF DUPLICATE VALUES EXIST IN B, THE CORRESPONDING ROWS   RSRT 380
C           OF A ARE MOVED TO R IN THE SAME ORDER AS IN A.              RSRT 390
C                                                                       RSRT 400
C     ..................................................................RSRT 410
C                                                                       RSRT 420
      SUBROUTINE RSRT(A,B,R,N,M,MS)                                     RSRT 430
      DIMENSION A(1),B(1),R(1)                                          RSRT 440
C                                                                       RSRT 450
C        MOVE SORTING KEY VECTOR TO FIRST COLUMN OF OUTPUT MATRIX       RSRT 460
C        AND BUILD ORIGINAL SEQUENCE LIST IN SECOND COLUMN              RSRT 470
C                                                                       RSRT 480
      DO 10 I=1,N                                                       RSRT 490
      R(I)=B(I)                                                         RSRT 500
      I2=I+N                                                            RSRT 510
   10 R(I2)=I                                                           RSRT 520
C                                                                       RSRT 530
C        SORT ELEMENTS IN SORTING KEY VECTOR (ORIGINAL SEQUENCE LIST    RSRT 540
C        IS RESEQUENCED ACCORDINGLY)                                    RSRT 550
C                                                                       RSRT 560
      L=N+1                                                             RSRT 570
   20 ISORT=0                                                           RSRT 580
      L=L-1                                                             RSRT 590
      DO 40 I=2,L                                                       RSRT 600
      IF(R(I)-R(I-1)) 30,40,40                                          RSRT 610
   30 ISORT=1                                                           RSRT 620
      RSAVE=R(I)                                                        RSRT 630
      R(I)=R(I-1)                                                       RSRT 640
      R(I-1)=RSAVE                                                      RSRT 650
      I2=I+N                                                            RSRT 660
      SAVER=R(I2)                                                       RSRT 670
      R(I2)=R(I2-1)                                                     RSRT 680
      R(I2-1)=SAVER                                                     RSRT 690
   40 CONTINUE                                                          RSRT 700
      IF(ISORT) 20,50,20                                                RSRT 710
C                                                                       RSRT 720
C        MOVE ROWS FROM MATRIX A TO MATRIX R (NUMBER IN SECOND COLUMN   RSRT 730
C        OF R REPRESENTS ROW NUMBER OF MATRIX A TO BE MOVED)            RSRT 740
C                                                                       RSRT 750
   50 DO 80 I=1,N                                                       RSRT 760
C                                                                       RSRT 770
C        GET ROW NUMBER IN MATRIX A                                     RSRT 780
C                                                                       RSRT 790
      I2=I+N                                                            RSRT 800
      IN=R(I2)                                                          RSRT 810
C                                                                       RSRT 820
      IR=I-N                                                            RSRT 830
      DO 80 J=1,M                                                       RSRT 840
C                                                                       RSRT 850
C        LOCATE ELEMENT IN OUTPUT MATRIX                                RSRT 860
C                                                                       RSRT 870
      IR=IR+N                                                           RSRT 880
C                                                                       RSRT 890
C        LOCATE ELEMENT IN INPUT MATRIX                                 RSRT 900
C                                                                       RSRT 910
      CALL LOC(IN,J,IA,N,M,MS)                                          RSRT 920
C                                                                       RSRT 930
C        TEST FOR ZERO ELEMENT IN DIAGONAL MATRIX                       RSRT 940
C                                                                       RSRT 950
      IF(IA) 60,70,60                                                   RSRT 960
C                                                                       RSRT 970
C        MOVE ELEMENT TO OUTPUT MATRIX                                  RSRT 980
C                                                                       RSRT 990
   60 R(IR)=A(IA)                                                       RSRT1000
      GO TO 80                                                          RSRT1010
   70 R(IR)=0                                                           RSRT1020
   80 CONTINUE                                                          RSRT1030
      RETURN                                                            RSRT1040
      END                                                               RSRT1050
