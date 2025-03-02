C                                                                       XCPY  10
C     ..................................................................XCPY  20
C                                                                       XCPY  30
C        SUBROUTINE XCPY                                                XCPY  40
C                                                                       XCPY  50
C        PURPOSE                                                        XCPY  60
C           COPY A PORTION OF A MATRIX                                  XCPY  70
C                                                                       XCPY  80
C        USAGE                                                          XCPY  90
C           CALL XCPY(A,R,L,K,NR,MR,NA,MA,MS)                           XCPY 100
C                                                                       XCPY 110
C        DESCRIPTION OF PARAMETERS                                      XCPY 120
C           A  - NAME OF INPUT MATRIX                                   XCPY 130
C           R  - NAME OF OUTPUT MATRIX                                  XCPY 140
C           L  - ROW OF A WHERE FIRST ELEMENT OF R CAN BE FOUND         XCPY 150
C           K  - COLUMN OF A WHERE FIRST ELEMENT OF R CAN BE FOUND      XCPY 160
C           NR - NUMBER OF ROWS TO BE COPIED INTO R                     XCPY 170
C           MR - NUMBER OF COLUMNS TO BE COPIED INTO R                  XCPY 180
C           NA - NUMBER OF ROWS IN A                                    XCPY 190
C           MA - NUMBER OF COLUMNS IN A                                 XCPY 200
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         XCPY 210
C                  0 - GENERAL                                          XCPY 220
C                  1 - SYMMETRIC                                        XCPY 230
C                  2 - DIAGONAL                                         XCPY 240
C                                                                       XCPY 250
C        REMARKS                                                        XCPY 260
C           MATRIX R CANNOT BE IN THE SAME LOCATION AS MATRIX A         XCPY 270
C           MATRIX R IS ALWAYS A GENERAL MATRIX                         XCPY 280
C                                                                       XCPY 290
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  XCPY 300
C           LOC                                                         XCPY 310
C                                                                       XCPY 320
C        METHOD                                                         XCPY 330
C           MATRIX R IS FORMED BY COPYING A PORTION OF MATRIX A. THIS   XCPY 340
C           IS DONE BY EXTRACTING NR ROWS AND MR COLUMNS OF MATRIX A,   XCPY 350
C           STARTING WITH ELEMENT AT ROW L, COLUMN K                    XCPY 360
C                                                                       XCPY 370
C     ..................................................................XCPY 380
C                                                                       XCPY 390
      SUBROUTINE XCPY(A,R,L,K,NR,MR,NA,MA,MS)                           XCPY 400
      DIMENSION A(1),R(1)                                               XCPY 410
C                                                                       XCPY 420
C        INITIALIZE                                                     XCPY 430
C                                                                       XCPY 440
      IR=0                                                              XCPY 450
      L2=L+NR-1                                                         XCPY 460
      K2=K+MR-1                                                         XCPY 470
C                                                                       XCPY 480
      DO 5 J=K,K2                                                       XCPY 490
      DO 5 I=L,L2                                                       XCPY 500
      IR=IR+1                                                           XCPY 510
      R(IR)=0.0                                                         XCPY 520
C                                                                       XCPY 530
C        LOCATE ELEMENT FOR ANY MATRIX STORAGE MODE                     XCPY 540
C                                                                       XCPY 550
      CALL LOC(I,J,IA,NA,MA,MS)                                         XCPY 560
C                                                                       XCPY 570
C        TEST FOR ZERO ELEMENT IN DIAGONAL MATRIX                       XCPY 580
C                                                                       XCPY 590
      IF(IA) 4,5,4                                                      XCPY 600
    4 R(IR)=A(IA)                                                       XCPY 610
    5 CONTINUE                                                          XCPY 620
      RETURN                                                            XCPY 630
      END                                                               XCPY 640
