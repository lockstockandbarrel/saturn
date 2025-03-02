C                                                                       MATA  10
C     ..................................................................MATA  20
C                                                                       MATA  30
C        SUBROUTINE MATA                                                MATA  40
C                                                                       MATA  50
C        PURPOSE                                                        MATA  60
C           PREMULTIPLY A MATRIX BY ITS TRANSPOSE TO FORM A             MATA  70
C           SYMMETRIC MATRIX                                            MATA  80
C                                                                       MATA  90
C        USAGE                                                          MATA 100
C           CALL MATA(A,R,N,M,MS)                                       MATA 110
C                                                                       MATA 120
C        DESCRIPTION OF PARAMETERS                                      MATA 130
C           A  - NAME OF INPUT MATRIX                                   MATA 140
C           R  - NAME OF OUTPUT MATRIX                                  MATA 150
C           N  - NUMBER OF ROWS IN A                                    MATA 160
C           M  - NUMBER OF COLUMNS IN A. ALSO NUMBER OF ROWS AND        MATA 170
C                NUMBER OF COLUMNS OF R.                                MATA 180
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         MATA 190
C                  0 - GENERAL                                          MATA 200
C                  1 - SYMMETRIC                                        MATA 210
C                  2 - DIAGONAL                                         MATA 220
C                                                                       MATA 230
C        REMARKS                                                        MATA 240
C           MATRIX R CANNOT BE IN THE SAME LOCATION AS MATRIX A         MATA 250
C           MATRIX R IS ALWAYS A SYMMETRIC MATRIX WITH A STORAGE MODE=1 MATA 260
C                                                                       MATA 270
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  MATA 280
C           LOC                                                         MATA 290
C                                                                       MATA 300
C        METHOD                                                         MATA 310
C           CALCULATION OF (A TRANSPOSE A) RESULTS IN A SYMMETRIC MATRIXMATA 320
C           REGARDLESS OF THE STORAGE MODE OF THE INPUT MATRIX. THE     MATA 330
C           ELEMENTS OF MATRIX A ARE NOT CHANGED.                       MATA 340
C                                                                       MATA 350
C     ..................................................................MATA 360
C                                                                       MATA 370
      SUBROUTINE MATA(A,R,N,M,MS)                                       MATA 380
      DIMENSION A(1),R(1)                                               MATA 390
C                                                                       MATA 400
      DO 60 K=1,M                                                       MATA 410
      KX=(K*K-K)/2                                                      MATA 420
      DO 60 J=1,M                                                       MATA 430
      IF(J-K) 10,10,60                                                  MATA 440
   10 IR=J+KX                                                           MATA 450
      R(IR)=0                                                           MATA 460
      DO 60 I=1,N                                                       MATA 470
      IF(MS) 20,40,20                                                   MATA 480
   20 CALL LOC(I,J,IA,N,M,MS)                                           MATA 490
      CALL LOC(I,K,IB,N,M,MS)                                           MATA 500
      IF(IA) 30,60,30                                                   MATA 510
   30 IF(IB) 50,60,50                                                   MATA 520
   40 IA=N*(J-1)+I                                                      MATA 530
      IB=N*(K-1)+I                                                      MATA 540
   50 R(IR)=R(IR)+A(IA)*A(IB)                                           MATA 550
   60 CONTINUE                                                          MATA 560
      RETURN                                                            MATA 570
      END                                                               MATA 580
