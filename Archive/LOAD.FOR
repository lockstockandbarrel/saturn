C                                                                       LOAD  10
C     ..................................................................LOAD  20
C                                                                       LOAD  30
C        SUBROUTINE LOAD                                                LOAD  40
C                                                                       LOAD  50
C        PURPOSE                                                        LOAD  60
C           COMPUTE A FACTOR MATRIX (LOADING) FROM EIGENVALUES AND      LOAD  70
C           ASSOCIATED EIGENVECTORS.  THIS SUBROUTINE NORMALLY OCCURS   LOAD  80
C           IN A SEQUENCE OF CALLS TO SUBROUTINES CORRE, EIGEN, TRACE,  LOAD  90
C           LOAD, AND VARMX IN THE PERFORMANCE OF A FACTOR ANALYSIS.    LOAD 100
C                                                                       LOAD 110
C        USAGE                                                          LOAD 120
C           CALL LOAD (M,K,R,V)                                         LOAD 130
C                                                                       LOAD 140
C        DESCRIPTION OF PARAMETERS                                      LOAD 150
C           M     - NUMBER OF VARIABLES.                                LOAD 160
C           K     - NUMBER OF FACTORS. K MUST BE GREATER THAN OR EQUAL  LOAD 170
C                   TO 1 AND LESS THAN OR EQUAL TO M.                   LOAD 171
C           R     - A MATRIX (SYMMETRIC AND STORED IN COMPRESSED FORM   LOAD 180
C                   WITH ONLY UPPER TRIANGLE BY COLUMN IN CORE) CON-    LOAD 190
C                   TAINING EIGENVALUES IN DIAGONAL.  EIGENVALUES ARE   LOAD 200
C                   ARRANGED IN DESCENDING ORDER, AND FIRST K           LOAD 210
C                   EIGENVALUES ARE USED BY THIS SUBROUTINE.  THE ORDER LOAD 220
C                   OF MATRIX R IS M BY M.  ONLY M*(M+1)/2 ELEMENTS ARE LOAD 230
C                   IN STORAGE.  (STORAGE MODE OF 1)                    LOAD 240
C           V     - WHEN THIS SUBROUTINE IS CALLED, MATRIX V (M X M)    LOAD 250
C                   CONTAINS EIGENVECTORS COLUMNWISE.  UPON RETURNING TOLOAD 260
C                   THE CALLING PROGRAM, MATRIX V CONTAINS A FACTOR     LOAD 270
C                   MATRIX (M X K).                                     LOAD 280
C                                                                       LOAD 290
C        REMARKS                                                        LOAD 300
C           NONE                                                        LOAD 310
C                                                                       LOAD 320
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  LOAD 330
C           NONE                                                        LOAD 340
C                                                                       LOAD 350
C        METHOD                                                         LOAD 360
C           NORMALIZED EIGENVECTORS ARE CONVERTED TO THE FACTOR PATTERN LOAD 370
C           BY MULTIPLYING THE ELEMENTS OF EACH VECTOR BY THE SQUARE    LOAD 380
C           ROOT OF THE CORRESPONDING EIGENVALUE.                       LOAD 390
C                                                                       LOAD 400
C     ..................................................................LOAD 410
C                                                                       LOAD 420
      SUBROUTINE LOAD (M,K,R,V)                                         LOAD 430
      DIMENSION R(1),V(1)                                               LOAD 440
C                                                                       LOAD 450
C        ...............................................................LOAD 460
C                                                                       LOAD 470
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE  LOAD 480
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION      LOAD 490
C        STATEMENT WHICH FOLLOWS.                                       LOAD 500
C                                                                       LOAD 510
C     DOUBLE PRECISION R,V,SQ                                           LOAD 520
C                                                                       LOAD 530
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS    LOAD 540
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS      LOAD 550
C        ROUTINE.                                                       LOAD 560
C                                                                       LOAD 570
C        THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO      LOAD 580
C        CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  SQRT IN STATEMENT LOAD 590
C        150 MUST BE CHANGED TO DSQRT.                                  LOAD 600
C                                                                       LOAD 610
C        ...............................................................LOAD 620
C                                                                       LOAD 630
      L=0                                                               LOAD 640
      JJ=0                                                              LOAD 650
      DO 160 J=1,K                                                      LOAD 660
      JJ=JJ+J                                                           LOAD 670
  150 SQ= SQRT(R(JJ))                                                   LOAD 680
      DO 160 I=1,M                                                      LOAD 690
      L=L+1                                                             LOAD 700
  160 V(L)=SQ*V(L)                                                      LOAD 710
      RETURN                                                            LOAD 720
      END                                                               LOAD 730
