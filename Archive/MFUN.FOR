C                                                                       MFUN  10
C     ..................................................................MFUN  20
C                                                                       MFUN  30
C        SUBROUTINE MFUN                                                MFUN  40
C                                                                       MFUN  50
C        PURPOSE                                                        MFUN  60
C           APPLY A FUNCTION TO EACH ELEMENT OF A MATRIX TO FORM A      MFUN  70
C           RESULTANT MATRIX                                            MFUN  80
C                                                                       MFUN  90
C        USAGE                                                          MFUN 100
C           CALL MFUN (A,F,R,N,M,MS)                                    MFUN 110
C           AN EXTERNAL STATEMENT MUST PRECEDE CALL STATEMENT IN ORDER  MFUN 120
C           TO IDENTIFY PARAMETER F AS THE NAME OF A FUNCTION           MFUN 130
C                                                                       MFUN 140
C        DESCRIPTION OF PARAMETERS                                      MFUN 150
C           A - NAME OF INPUT MATRIX                                    MFUN 160
C           F - NAME OF FORTRAN-FURNISHED OR USER FUNCTION SUBPROGRAM   MFUN 170
C           R - NAME OF OUTPUT MATRIX                                   MFUN 180
C           N - NUMBER OF ROWS IN MATRIX A AND R                        MFUN 190
C           M - NUMBER OF COLUMNS IN MATRIX A AND R                     MFUN 200
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A (AND R) MFUN 210
C                  0 - GENERAL                                          MFUN 220
C                  1 - SYMMETRIC                                        MFUN 230
C                  2 - DIAGONAL                                         MFUN 240
C                                                                       MFUN 250
C        REMARKS                                                        MFUN 260
C           PRECISION IS DEPENDENT UPON PRECISION OF FUNCTION USED      MFUN 270
C                                                                       MFUN 280
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  MFUN 290
C           LOC                                                         MFUN 300
C                                                                       MFUN 310
C        METHOD                                                         MFUN 320
C           FUNCTION F IS APPLIED TO EACH ELEMENT OF MATRIX A           MFUN 330
C           TO FORM MATRIX R                                            MFUN 340
C                                                                       MFUN 350
C     ..................................................................MFUN 360
C                                                                       MFUN 370
      SUBROUTINE MFUN(A,F,R,N,M,MS)                                     MFUN 380
      DIMENSION A(1),R(1)                                               MFUN 390
C                                                                       MFUN 400
C        COMPUTE VECTOR LENGTH, IT                                      MFUN 410
C                                                                       MFUN 420
      CALL LOC(N,M,IT,N,M,MS)                                           MFUN 430
C                                                                       MFUN 440
C        BUILD MATRIX R FOR ANY STORAGE MODE                            MFUN 450
C                                                                       MFUN 460
      DO 5 I=1,IT                                                       MFUN 470
    5 R(I)=F(A(I))                                                      MFUN 480
      RETURN                                                            MFUN 490
      END                                                               MFUN 500
