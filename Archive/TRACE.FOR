C                                                                       TRAC  10
C     ..................................................................TRAC  20
C                                                                       TRAC  30
C        SUBROUTINE TRACE                                               TRAC  40
C                                                                       TRAC  50
C        PURPOSE                                                        TRAC  60
C           COMPUTE CUMULATIVE PERCENTAGE OF EIGENVALUES GREATER THAN   TRAC  70
C           OR EQUAL TO A CONSTANT SPECIFIED BY THE USER.  THIS SUB-    TRAC  80
C           ROUTINE NORMALLY OCCURS IN A SEQUENCE OF CALLS TO SUB-      TRAC  90
C           ROUTINES CORRE, EIGEN, TRACE, LOAD, AND VARMX IN THE PER-   TRAC 100
C           FORMANCE OF A FACTOR ANALYSIS.                              TRAC 110
C                                                                       TRAC 120
C        USAGE                                                          TRAC 130
C           CALL TRACE (M,R,CON,K,D)                                    TRAC 140
C                                                                       TRAC 150
C        DESCRIPTION OF PARAMETERS                                      TRAC 160
C           M     - NUMBER OF VARIABLES. M MUST BE > OR = TO 1          TRAC 170
C           R     - INPUT MATRIX (SYMMETRIC AND STORED IN COMPRESSED    TRAC 180
C                   FORM WITH ONLY UPPER TRIANGLE BY COLUMN IN CORE)    TRAC 190
C                   CONTAINING EIGENVALUES IN DIAGONAL.  EIGENVALUES ARETRAC 200
C                   ARRANGED IN DESCENDING ORDER.  THE ORDER OF MATRIX RTRAC 210
C                   IS M BY M.  ONLY M*(M+1)/2 ELEMENTS ARE IN STORAGE. TRAC 220
C                   (STORAGE MODE OF 1)                                 TRAC 230
C           CON   - A CONSTANT USED TO DECIDE HOW MANY EIGENVALUES TO   TRAC 240
C                   RETAIN.  CUMULATIVE PERCENTAGE OF EIGENVALUES       TRAC 250
C                   WHICH ARE GREATER THAN OR EQUAL TO THIS VALUE IS    TRAC 260
C                   CALCULATED.                                         TRAC 270
C           K     - OUTPUT VARIABLE CONTAINING THE NUMBER OF EIGENVALUESTRAC 280
C                   GREATER THAN OR EQUAL TO CON.  (K IS THE NUMBER OF  TRAC 290
C                   FACTORS.)                                           TRAC 300
C           D     - OUTPUT VECTOR OF LENGTH M CONTAINING CUMULATIVE     TRAC 310
C                   PERCENTAGE OF EIGENVALUES WHICH ARE GREATER THAN    TRAC 320
C                   OR EQUAL TO CON.                                    TRAC 330
C        REMARKS                                                        TRAC 340
C           NONE                                                        TRAC 350
C                                                                       TRAC 360
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  TRAC 370
C           NONE                                                        TRAC 380
C                                                                       TRAC 390
C        METHOD                                                         TRAC 400
C           EACH EIGENVALUE GREATER THAN OR EQUAL TO CON IS DIVIDED BY MTRAC 410
C           AND THE RESULT IS ADDED TO THE PREVIOUS TOTAL TO OBTAIN     TRAC 420
C           THE CUMULATIVE PERCENTAGE FOR EACH EIGENVALUE.              TRAC 430
C                                                                       TRAC 440
C     ..................................................................TRAC 450
C                                                                       TRAC 460
      SUBROUTINE TRACE (M,R,CON,K,D)                                    TRAC 470
      DIMENSION R(1),D(1)                                               TRAC 480
C                                                                       TRAC 490
C        ...............................................................TRAC 500
C                                                                       TRAC 510
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE  TRAC 520
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION      TRAC 530
C        STATEMENT WHICH FOLLOWS.                                       TRAC 540
C                                                                       TRAC 550
C     DOUBLE PRECISION R,D                                              TRAC 560
C                                                                       TRAC 570
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS    TRAC 580
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS      TRAC 590
C        ROUTINE.                                                       TRAC 600
C                                                                       TRAC 610
C        ...............................................................TRAC 620
C                                                                       TRAC 630
      FM=M                                                              TRAC 640
      L=0                                                               TRAC 650
      DO 100 I=1,M                                                      TRAC 660
      L=L+I                                                             TRAC 670
  100 D(I)=R(L)                                                         TRAC 680
      K=0                                                               TRAC 690
C                                                                       TRAC 700
C     TEST WHETHER I-TH EIGENVALUE IS GREATER                           TRAC 710
C     THAN OR EQUAL TO THE CONSTANT                                     TRAC 720
C                                                                       TRAC 730
      DO 110 I=1,M                                                      TRAC 740
      IF(D(I)-CON) 120, 105, 105                                        TRAC 750
  105 K=K+1                                                             TRAC 760
  110 D(I)=D(I)/FM                                                      TRAC 770
C                                                                       TRAC 780
C     COMPUTE CUMULATIVE PERCENTAGE OF EIGENVALUES                      TRAC 790
C                                                                       TRAC 800
  120 DO 130 I=2,K                                                      TRAC 810
  130 D(I)=D(I)+D(I-1)                                                  TRAC 820
      RETURN                                                            TRAC 830
      END                                                               TRAC 840
