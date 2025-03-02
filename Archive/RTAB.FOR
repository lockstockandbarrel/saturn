C                                                                       RTAB  10
C     ..................................................................RTAB  20
C                                                                       RTAB  30
C        SUBROUTINE RTAB                                                RTAB  40
C                                                                       RTAB  50
C        PURPOSE                                                        RTAB  60
C           TABULATE ROWS OF A MATRIX TO FORM A SUMMARY MATRIX          RTAB  70
C                                                                       RTAB  80
C        USAGE                                                          RTAB  90
C           CALL RTAB(A,B,R,S,N,M,MS,L)                                 RTAB 100
C                                                                       RTAB 110
C        DESCRIPTION OF PARAMETERS                                      RTAB 120
C           A - NAME OF INPUT MATRIX                                    RTAB 130
C           B - NAME OF INPUT VECTOR OF LENGTH N CONTAINING KEY         RTAB 140
C           R - NAME OF OUTPUT MATRIX CONTAINING SUMMARY OF ROW DATA.   RTAB 150
C               IT IS INITIALLY SET TO ZERO BY THIS SUBROUTINE.         RTAB 160
C           S - NAME OF OUTPUT VECTOR OF LENGTH L+1 CONTAINING COUNTS   RTAB 170
C           N - NUMBER OF ROWS IN A                                     RTAB 180
C           M - NUMBER OF COLUMNS IN A AND R                            RTAB 190
C           L - NUMBER OF ROWS IN R                                     RTAB 200
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         RTAB 210
C                  0 - GENERAL                                          RTAB 220
C                  1 - SYMMETRIC                                        RTAB 230
C                  2 - DIAGONAL                                         RTAB 240
C                                                                       RTAB 250
C        REMARKS                                                        RTAB 260
C           MATRIX R IS ALWAYS A GENERAL MATRIX                         RTAB 270
C                                                                       RTAB 280
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  RTAB 290
C           LOC                                                         RTAB 300
C           RADD                                                        RTAB 310
C                                                                       RTAB 320
C        METHOD                                                         RTAB 330
C           ROWS OF DATA IN MATRIX A ARE TABULATED BASED ON THE KEY     RTAB 340
C           CONTAINED IN VECTOR B. THE FLOATING POINT NUMBER IN B(I) IS RTAB 350
C           TRUNCATED TO FORM J. THE ITH ROW OF A IS ADDED TO THE JTH   RTAB 360
C           ROW OF R ELEMENT BY ELEMENT AND ONE IS ADDED TO S(J). IF J  RTAB 370
C           IS NOT BETWEEN ONE AND L, ONE IS ADDED TO S(L+1). THIS      RTAB 380
C           PROCEDURE IS REPEATED FOR EVERY ELEMENT IN VECTOR B.        RTAB 390
C           UPON COMPLETION, THE OUTPUT MATRIX R CONTAINS A SUMMARY OF  RTAB 400
C           ROW DATA AS SPECIFIED BY VECTOR B. EACH ELEMENT IN VECTOR S RTAB 410
C           CONTAINS A COUNT OF THE NUMBER OF ROWS OF A USED TO FORM THERTAB 420
C           CORRESPONDING ROW OF R. ELEMENT S(L+1) CONTAINS A COUNT OF  RTAB 430
C           THE NUMBER OF ROWS OF A NOT INCLUDED IN R AS A RESULT OF J  RTAB 440
C           BEING LESS THAN ONE OR GREATER THAN L.                      RTAB 450
C                                                                       RTAB 460
C     ..................................................................RTAB 470
C                                                                       RTAB 480
      SUBROUTINE RTAB(A,B,R,S,N,M,MS,L)                                 RTAB 490
      DIMENSION A(1),B(1),R(1),S(1)                                     RTAB 500
C                                                                       RTAB 510
C        CLEAR OUTPUT AREAS                                             RTAB 520
C                                                                       RTAB 530
      CALL LOC(M,L,IT,M,L,0)                                            RTAB 540
      DO 10 IR=1,IT                                                     RTAB 550
   10 R(IR)=0.0                                                         RTAB 560
      DO 20 IS=1,L                                                      RTAB 570
   20 S(IS)=0.0                                                         RTAB 580
      S(L+1)=0.0                                                        RTAB 590
C                                                                       RTAB 600
      DO 60 I=1,N                                                       RTAB 610
C                                                                       RTAB 620
C        TEST FOR THE KEY OUTSIDE THE RANGE                             RTAB 630
C                                                                       RTAB 640
      JR=B(I)                                                           RTAB 650
      IF (JR-1) 50,40,30                                                RTAB 660
   30 IF (JR-L) 40,40,50                                                RTAB 670
C                                                                       RTAB 680
C                                                                       RTAB 690
C        ADD ROW OF A TO ROW OF R AND 1 TO COUNT                        RTAB 700
C                                                                       RTAB 710
   40 CALL RADD (A,I,R,JR,N,M,MS,L)                                     RTAB 720
      S(JR)=S(JR)+1.0                                                   RTAB 730
      GO TO 60                                                          RTAB 740
C                                                                       RTAB 750
   50 S(L+1)=S(L+1)+1.0                                                 RTAB 760
   60 CONTINUE                                                          RTAB 770
      RETURN                                                            RTAB 780
      END                                                               RTAB 790
