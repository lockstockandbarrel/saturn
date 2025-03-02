C                                                                       CSRT  10
C     ..................................................................CSRT  20
C                                                                       CSRT  30
C        SUBROUTINE CSRT                                                CSRT  40
C                                                                       CSRT  50
C        PURPOSE                                                        CSRT  60
C           SORT COLUMNS OF A MATRIX                                    CSRT  70
C                                                                       CSRT  80
C        USAGE                                                          CSRT  90
C           CALL CSRT(A,B,R,N,M,MS)                                     CSRT 100
C                                                                       CSRT 110
C        DESCRIPTION OF PARAMETERS                                      CSRT 120
C           A - NAME OF INPUT MATRIX TO BE SORTED                       CSRT 130
C           B - NAME OF INPUT VECTOR WHICH CONTAINS SORTING KEY         CSRT 140
C           R - NAME OF SORTED OUTPUT MATRIX                            CSRT 150
C           N - NUMBER OF ROWS IN A AND R                               CSRT 160
C           M - NUMBER OF COLUMNS IN A AND R AND LENGTH OF B            CSRT 170
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         CSRT 180
C                  0 - GENERAL                                          CSRT 190
C                  1 - SYMMETRIC                                        CSRT 200
C                  2 - DIAGONAL                                         CSRT 210
C                                                                       CSRT 220
C        REMARKS                                                        CSRT 230
C           MATRIX R CANNOT BE IN THE SAME LOCATION AS MATRIX A         CSRT 240
C           MATRIX R IS ALWAYS A GENERAL MATRIX                         CSRT 250
C           M MUST BE GREATER THAN ONE.                                 CSRT 260
C                                                                       CSRT 270
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  CSRT 280
C           LOC                                                         CSRT 290
C           CCPY                                                        CSRT 300
C                                                                       CSRT 310
C        METHOD                                                         CSRT 320
C           COLUMNS OF INPUT MATRIX A ARE SORTED TO FORM OUTPUT MATRIX  CSRT 330
C           R. THE SORTED COLUMN SEQUENCE IS DETERMINED BY THE VALUES OFCSRT 340
C           ELEMENTS IN ROW VECTOR B. THE LOWEST VALUED ELEMENT IN      CSRT 350
C           B WILL CAUSE THE CORRESPONDING COLUMN OF A TO BE PLACED IN  CSRT 360
C           THE FIRST COLUMN OF R. THE HIGHEST VALUED ELEMENT OF B WILL CSRT 370
C           CAUSE THE CORRESPONDING ROW OF A TO BE PLACED IN THE LAST   CSRT 380
C           COLUMN OF R. IF DUPLICATE VALUES EXIST IN B, THE            CSRT 390
C           CORRESPONDING COLUMNS OF A ARE MOVED TO R IN THE SAME ORDER CSRT 400
C           AS IN A.                                                    CSRT 410
C                                                                       CSRT 420
C     ..................................................................CSRT 430
C                                                                       CSRT 440
      SUBROUTINE CSRT(A,B,R,N,M,MS)                                     CSRT 450
      DIMENSION A(1),B(1),R(1)                                          CSRT 460
C                                                                       CSRT 470
C        MOVE SORTING KEY VECTOR TO FIRST ROW OF OUTPUT MATRIX          CSRT 480
C        AND BUILD ORIGINAL SEQUENCE LIST IN SECOND ROW                 CSRT 490
C                                                                       CSRT 500
      IK=1                                                              CSRT 510
      DO 10 J=1,M                                                       CSRT 520
      R(IK)=B(J)                                                        CSRT 530
      R(IK+1)=J                                                         CSRT 540
   10 IK=IK+N                                                           CSRT 550
C                                                                       CSRT 560
C        SORT ELEMENTS IN SORTING KEY VECTOR (ORIGINAL SEQUENCE LIST    CSRT 570
C        IS RESEQUENCED ACCORDINGLY)                                    CSRT 580
C                                                                       CSRT 590
      L=M+1                                                             CSRT 600
   20 ISORT=0                                                           CSRT 610
      L=L-1                                                             CSRT 620
      IP=1                                                              CSRT 630
      IQ=N+1                                                            CSRT 640
      DO 50 J=2,L                                                       CSRT 650
      IF(R(IQ)-R(IP)) 30,40,40                                          CSRT 660
   30 ISORT=1                                                           CSRT 670
      RSAVE=R(IQ)                                                       CSRT 680
      R(IQ)=R(IP)                                                       CSRT 690
      R(IP)=RSAVE                                                       CSRT 700
      SAVER=R(IQ+1)                                                     CSRT 710
      R(IQ+1)=R(IP+1)                                                   CSRT 720
      R(IP+1)=SAVER                                                     CSRT 730
   40 IP=IP+N                                                           CSRT 740
      IQ=IQ+N                                                           CSRT 750
   50 CONTINUE                                                          CSRT 760
      IF(ISORT) 20,60,20                                                CSRT 770
C                                                                       CSRT 780
C        MOVE COLUMNS FROM MATRIX A TO MATRIX R (NUMBER IN SECOND ROW   CSRT 790
C        OF R REPRESENTS COLUMN NUMBER OF MATRIX A TO BE MOVED)         CSRT 800
C                                                                       CSRT 810
   60 IQ=-N                                                             CSRT 820
      DO 70 J=1,M                                                       CSRT 830
      IQ=IQ+N                                                           CSRT 840
C                                                                       CSRT 850
C        GET COLUMN NUMBER IN MATRIX A                                  CSRT 860
C                                                                       CSRT 870
      I2=IQ+2                                                           CSRT 880
      IN=R(I2)                                                          CSRT 890
C                                                                       CSRT 900
C        MOVE COLUMN                                                    CSRT 910
C                                                                       CSRT 920
      IR=IQ+1                                                           CSRT 930
      CALL CCPY(A,IN,R(IR),N,M,MS)                                      CSRT 940
   70 CONTINUE                                                          CSRT 950
      RETURN                                                            CSRT 960
      END                                                               CSRT 970
