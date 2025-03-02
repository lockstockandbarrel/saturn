C                                                                       RADD  10
C     ..................................................................RADD  20
C                                                                       RADD  30
C        SUBROUTINE RADD                                                RADD  40
C                                                                       RADD  50
C        PURPOSE                                                        RADD  60
C           ADD ROW OF ONE MATRIX TO ROW OF ANOTHER MATRIX              RADD  70
C                                                                       RADD  80
C        USAGE                                                          RADD  90
C           CALL RADD(A,IRA,R,IRR,N,M,MS,L)                             RADD 100
C                                                                       RADD 110
C        DESCRIPTION OF PARAMETERS                                      RADD 120
C           A   - NAME OF INPUT MATRIX                                  RADD 130
C           IRA - ROW IN MATRIX A TO BE ADDED TO ROW IRR OF MATRIX R    RADD 140
C           R   - NAME OF OUTPUT MATRIX                                 RADD 150
C           IRR - ROW IN MATRIX R WHERE SUMMATION IS DEVELOPED          RADD 160
C           N   - NUMBER OF ROWS IN A                                   RADD 170
C           M   - NUMBER OF COLUMNS IN A AND R                          RADD 180
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         RADD 190
C                  0 - GENERAL                                          RADD 200
C                  1 - SYMMETRIC                                        RADD 210
C                  2 - DIAGONAL                                         RADD 220
C           L   - NUMBER OF ROWS IN R                                   RADD 230
C                                                                       RADD 240
C        REMARKS                                                        RADD 250
C           MATRIX R MUST BE A GENERAL MATRIX                           RADD 260
C           MATRIX R CANNOT BE IN THE SAME LOCATION AS MATRIX A UNLESS  RADD 270
C           A IS GENERAL                                                RADD 280
C                                                                       RADD 290
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  RADD 300
C           LOC                                                         RADD 310
C                                                                       RADD 320
C        METHOD                                                         RADD 330
C           EACH ELEMENT OF ROW IRA OF MATRIX A IS ADDED TO             RADD 340
C           CORRESPONDING ELEMENT OF ROW IRR OF MATRIX R                RADD 350
C                                                                       RADD 360
C     ..................................................................RADD 370
C                                                                       RADD 380
      SUBROUTINE RADD(A,IRA,R,IRR,N,M,MS,L)                             RADD 390
      DIMENSION A(1),R(1)                                               RADD 400
C                                                                       RADD 410
      IR=IRR-L                                                          RADD 420
      DO 2 J=1,M                                                        RADD 430
      IR=IR+L                                                           RADD 440
C                                                                       RADD 450
C        LOCATE INPUT ELEMENT FOR ANY MATRIX STORAGE MODE               RADD 460
C                                                                       RADD 470
      CALL LOC(IRA,J,IA,N,M,MS)                                         RADD 480
C                                                                       RADD 490
C        TEST FOR ZERO ELEMENT IN DIAGONAL MATRIX                       RADD 500
C                                                                       RADD 510
      IF(IA) 1,2,1                                                      RADD 520
C                                                                       RADD 530
C        ADD ELEMENTS                                                   RADD 540
C                                                                       RADD 550
    1 R(IR)=R(IR)+A(IA)                                                 RADD 560
    2 CONTINUE                                                          RADD 570
      RETURN                                                            RADD 580
      END                                                               RADD 590
