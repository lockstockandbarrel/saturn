C                                                                       CTAB  10
C     ..................................................................CTAB  20
C                                                                       CTAB  30
C        SUBROUTINE CTAB                                                CTAB  40
C                                                                       CTAB  50
C        PURPOSE                                                        CTAB  60
C           TABULATE COLUMNS OF A MATRIX TO FORM A SUMMARY MATRIX       CTAB  70
C                                                                       CTAB  80
C        USAGE                                                          CTAB  90
C           CALL CTAB(A,B,R,S,N,M,MS,L)                                 CTAB 100
C                                                                       CTAB 110
C        DESCRIPTION OF PARAMETERS                                      CTAB 120
C           A - NAME OF INPUT MATRIX                                    CTAB 130
C           B - NAME OF INPUT VECTOR OF LENGTH M CONTAINING KEY         CTAB 140
C           R - NAME OF OUTPUT MATRIX CONTAINING SUMMARY OF COLUMN DATA.CTAB 150
C               IT IS INITIALLY SET TO ZERO BY THIS SUBROUTINE.         CTAB 160
C           S - NAME OF OUTPUT VECTOR OF LENGTH L+1 CONTAINING COUNTS   CTAB 170
C           N - NUMBER OF ROWS IN A AND R                               CTAB 180
C           M - NUMBER OF COLUMNS IN A                                  CTAB 190
C           L - NUMBER OF COLUMNS IN R                                  CTAB 200
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         CTAB 210
C                  0 - GENERAL                                          CTAB 220
C                  1 - SYMMETRIC                                        CTAB 230
C                  2 - DIAGONAL                                         CTAB 240
C                                                                       CTAB 250
C        REMARKS                                                        CTAB 260
C           MATRIX R IS ALWAYS A GENERAL MATRIX                         CTAB 270
C                                                                       CTAB 280
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  CTAB 290
C           LOC                                                         CTAB 300
C           CADD                                                        CTAB 310
C                                                                       CTAB 320
C        METHOD                                                         CTAB 330
C           COLUMNS OF DATA IN MATRIX A ARE TABULATED BASED ON THE KEY  CTAB 340
C           CONTAINED IN VECTOR B. THE FLOATING POINT NUMBER IN B(I) IS CTAB 350
C           TRUNCATED TO FORM J. THE ITH COLUMN OF A IS ADDED TO THE JTHCTAB 360
C           COLUMN OF MATRIX R AND ONE IS ADDED TO S(J). IF THE VALUE OFCTAB 370
C           J IS NOT BETWEEN 1 AND L, ONE IS ADDED TO S(L+1)            CTAB 380
C           UPON COMPLETION, THE OUTPUT MATRIX R CONTAINS A SUMMARY OF  CTAB 390
C           COLUMN DATA AS SPECIFIED BY VECTOR B. EACH ELEMENT IN VECTORCTAB 400
C           S CONTAINS A COUNT OF THE NUMBER OF COLUMNS OF A USED TO    CTAB 410
C           FORM R. ELEMENT S(L+1) CONTAINS THE NUMBER OF COLUMNS OF A  CTAB 420
C           NOT INCLUDED IN R AS A RESULT OF J BEING LESS THAN ONE OR   CTAB 430
C           GREATER THAN L.                                             CTAB 440
C                                                                       CTAB 450
C     ..................................................................CTAB 460
C                                                                       CTAB 470
      SUBROUTINE CTAB(A,B,R,S,N,M,MS,L)                                 CTAB 480
      DIMENSION A(1),B(1),R(1),S(1)                                     CTAB 490
C                                                                       CTAB 500
C        CLEAR OUTPUT AREAS                                             CTAB 510
C                                                                       CTAB 520
      CALL LOC(N,L,IT,N,L,0)                                            CTAB 530
      DO 10 IR=1,IT                                                     CTAB 540
   10 R(IR)=0.0                                                         CTAB 550
      DO 20 IS=1,L                                                      CTAB 560
   20 S(IS)=0.0                                                         CTAB 570
      S(L+1)=0.0                                                        CTAB 580
C                                                                       CTAB 590
      DO 60 I=1,M                                                       CTAB 600
C                                                                       CTAB 610
C        TEST FOR THE KEY OUTSIDE THE RANGE                             CTAB 620
C                                                                       CTAB 630
      JR=B(I)                                                           CTAB 640
      IF (JR-1) 50,40,30                                                CTAB 650
   30 IF (JR-L) 40,40,50                                                CTAB 660
C                                                                       CTAB 670
C                                                                       CTAB 680
C        ADD COLUMN OF A TO COLUMN OF R AND 1 TO COUNT                  CTAB 690
C                                                                       CTAB 700
   40 CALL CADD (A,I,R,JR,N,M,MS,L)                                     CTAB 710
      S(JR)=S(JR)+1.0                                                   CTAB 720
      GO TO 60                                                          CTAB 730
C                                                                       CTAB 740
   50 S(L+1)=S(L+1)+1.0                                                 CTAB 750
   60 CONTINUE                                                          CTAB 760
      RETURN                                                            CTAB 770
      END                                                               CTAB 780
