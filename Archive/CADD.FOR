C                                                                       CADD  10
C     ..................................................................CADD  20
C                                                                       CADD  30
C        SUBROUTINE CADD                                                CADD  40
C                                                                       CADD  50
C        PURPOSE                                                        CADD  60
C           ADD COLUMN OF ONE MATRIX TO COLUMN OF ANOTHER MATRIX        CADD  70
C                                                                       CADD  80
C        USAGE                                                          CADD  90
C           CALL CADD(A,ICA,R,ICR,N,M,MS,L)                             CADD 100
C                                                                       CADD 110
C        DESCRIPTION OF PARAMETERS                                      CADD 120
C           A   - NAME OF INPUT MATRIX                                  CADD 130
C           ICA - COLUMN IN MATRIX A TO BE ADDED TO COLUMN ICR OF R     CADD 140
C           R   - NAME OF OUTPUT MATRIX                                 CADD 150
C           ICR - COLUMN IN MATRIX R WHERE SUMMATION IS DEVELOPED       CADD 160
C           N   - NUMBER OF ROWS IN A AND R                             CADD 170
C           M   - NUMBER OF COLUMNS IN A                                CADD 180
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         CADD 190
C                  0 - GENERAL                                          CADD 200
C                  1 - SYMMETRIC                                        CADD 210
C                  2 - DIAGONAL                                         CADD 220
C           L   - NUMBER OF COLUMNS IN R                                CADD 230
C                                                                       CADD 240
C        REMARKS                                                        CADD 250
C           MATRIX R MUST BE A GENERAL MATRIX                           CADD 260
C           MATRIX R CANNOT BE IN THE SAME LOCATION AS MATRIX A UNLESS  CADD 270
C           A IS GENERAL                                                CADD 280
C                                                                       CADD 290
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  CADD 300
C           LOC                                                         CADD 310
C                                                                       CADD 320
C        METHOD                                                         CADD 330
C           EACH ELEMENT OF COLUMN ICA OF MATRIX A IS ADDED TO          CADD 340
C           CORRESPONDING ELEMENT OF COLUMN ICR OF MATRIX R             CADD 350
C                                                                       CADD 360
C     ..................................................................CADD 370
C                                                                       CADD 380
      SUBROUTINE CADD(A,ICA,R,ICR,N,M,MS,L)                             CADD 390
      DIMENSION A(1),R(1)                                               CADD 400
C                                                                       CADD 410
      IR=N*(ICR-1)                                                      CADD 420
      DO 2 I=1,N                                                        CADD 430
      IR=IR+1                                                           CADD 440
C                                                                       CADD 450
C        LOCATE INPUT ELEMENT FOR ANY MATRIX STORAGE MODE               CADD 460
C                                                                       CADD 470
      CALL LOC(I,ICA,IA,N,M,MS)                                         CADD 480
C                                                                       CADD 490
C        TEST FOR ZERO ELEMENT IN DIAGONAL MATRIX                       CADD 500
C                                                                       CADD 510
      IF(IA) 1,2,1                                                      CADD 520
C                                                                       CADD 530
C        ADD ELEMENTS                                                   CADD 540
C                                                                       CADD 550
    1 R(IR)=R(IR)+A(IA)                                                 CADD 560
    2 CONTINUE                                                          CADD 570
      RETURN                                                            CADD 580
      END                                                               CADD 590
