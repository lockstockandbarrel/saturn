C                                                                       MTRA  10
C     ..................................................................MTRA  20
C                                                                       MTRA  30
C        SUBROUTINE MTRA                                                MTRA  40
C                                                                       MTRA  50
C        PURPOSE                                                        MTRA  60
C           TRANSPOSE A MATRIX                                          MTRA  70
C                                                                       MTRA  80
C        USAGE                                                          MTRA  90
C           CALL MTRA(A,R,N,M,MS)                                       MTRA 100
C                                                                       MTRA 110
C        DESCRIPTION OF PARAMETERS                                      MTRA 120
C           A - NAME OF MATRIX TO BE TRANSPOSED                         MTRA 130
C           R - NAME OF OUTPUT MATRIX                                   MTRA 140
C           N - NUMBER OF ROWS OF A AND COLUMNS OF R                    MTRA 150
C           M - NUMBER OF COLUMNS OF A AND ROWS OF R                    MTRA 160
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A (AND R) MTRA 170
C                  0 - GENERAL                                          MTRA 180
C                  1 - SYMMETRIC                                        MTRA 190
C                  2 - DIAGONAL                                         MTRA 200
C                                                                       MTRA 210
C        REMARKS                                                        MTRA 220
C           MATRIX R CANNOT BE IN THE SAME LOCATION AS MATRIX A         MTRA 230
C                                                                       MTRA 240
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  MTRA 250
C           MCPY                                                        MTRA 260
C                                                                       MTRA 270
C        METHOD                                                         MTRA 280
C           TRANSPOSE N BY M MATRIX A TO FORM M BY N MATRIX R BY MOVING MTRA 290
C           EACH ROW OF A INTO THE CORRESPONDING COLUMN OF R. IF MATRIX MTRA 300
C           A IS SYMMETRIC OR DIAGONAL, MATRIX R IS THE SAME AS A.      MTRA 310
C                                                                       MTRA 320
C     ..................................................................MTRA 330
C                                                                       MTRA 340
      SUBROUTINE MTRA(A,R,N,M,MS)                                       MTRA 350
      DIMENSION A(1),R(1)                                               MTRA 360
C                                                                       MTRA 370
C        IF MS IS 1 OR 2, COPY A                                        MTRA 380
C                                                                       MTRA 390
      IF(MS) 10,20,10                                                   MTRA 400
   10 CALL MCPY(A,R,N,N,MS)                                             MTRA 410
      RETURN                                                            MTRA 420
C                                                                       MTRA 430
C        TRANSPOSE GENERAL MATRIX                                       MTRA 440
C                                                                       MTRA 450
   20 IR=0                                                              MTRA 460
      DO 30 I=1,N                                                       MTRA 470
      IJ=I-N                                                            MTRA 480
      DO 30 J=1,M                                                       MTRA 490
      IJ=IJ+N                                                           MTRA 500
      IR=IR+1                                                           MTRA 510
   30 R(IR)=A(IJ)                                                       MTRA 520
      RETURN                                                            MTRA 530
      END                                                               MTRA 540
