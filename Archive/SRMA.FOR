C                                                                       SRMA  10
C     ..................................................................SRMA  20
C                                                                       SRMA  30
C        SUBROUTINE SRMA                                                SRMA  40
C                                                                       SRMA  50
C        PURPOSE                                                        SRMA  60
C           MULTIPLY ROW OF MATRIX BY A SCALAR AND ADD TO ANOTHER ROW   SRMA  70
C           OF THE SAME MATRIX                                          SRMA  80
C                                                                       SRMA  90
C        USAGE                                                          SRMA 100
C           CALL SRMA(A,C,N,M,LA,LB)                                    SRMA 110
C                                                                       SRMA 120
C        DESCRIPTION OF PARAMETERS                                      SRMA 130
C           A  - NAME OF MATRIX                                         SRMA 140
C           C  - SCALAR                                                 SRMA 150
C           N  - NUMBER OF ROWS IN A                                    SRMA 160
C           M  - NUMBER OF COLUMNS IN A                                 SRMA 170
C           LA - ROW IN A TO BE MULTIPLIED BY SCALAR                    SRMA 180
C           LB - ROW IN A TO WHICH PRODUCT IS ADDED                     SRMA 190
C                IF 0 IS SPECIFIED, PRODUCT REPLACES ELEMENTS IN ROW LA SRMA 200
C                                                                       SRMA 210
C        REMARKS                                                        SRMA 220
C           MATRIX A MUST BE A GENERAL MATRIX                           SRMA 230
C                                                                       SRMA 240
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  SRMA 250
C           NONE                                                        SRMA 260
C                                                                       SRMA 270
C        METHOD                                                         SRMA 280
C           EACH ELEMENT OF ROW LA IS MULTIPLIED BY SCALAR C AND THE    SRMA 290
C           PRODUCT IS ADDED TO THE CORRESPONDING ELEMENT OF ROW LB.    SRMA 300
C           ROW LA REMAINS UNAFFECTED BY THE OPERATION.                 SRMA 310
C           IF PARAMETER LB CONTAINS ZERO, MULTIPLICATION BY THE SCALAR SRMA 320
C           IS PERFORMED AND THE PRODUCT REPLACES ELEMENTS IN ROW LA.   SRMA 330
C                                                                       SRMA 340
C     ..................................................................SRMA 350
C                                                                       SRMA 360
      SUBROUTINE SRMA(A,C,N,M,LA,LB)                                    SRMA 370
      DIMENSION A(1)                                                    SRMA 380
C                                                                       SRMA 390
      LAJ=LA-N                                                          SRMA 400
      LBJ=LB-N                                                          SRMA 410
      DO 3 J=1,M                                                        SRMA 420
C                                                                       SRMA 430
C        LOCATE ELEMENT IN BOTH ROWS                                    SRMA 440
C                                                                       SRMA 450
      LAJ=LAJ+N                                                         SRMA 460
      LBJ=LBJ+N                                                         SRMA 470
C                                                                       SRMA 480
C        CHECK LB FOR ZERO                                              SRMA 490
C                                                                       SRMA 500
      IF(LB) 1,2,1                                                      SRMA 510
C                                                                       SRMA 520
C        IF NOT, MULTIPLY BY CONSTANT AND ADD TO OTHER ROW              SRMA 530
C                                                                       SRMA 540
    1 A(LBJ)=A(LAJ)*C+A(LBJ)                                            SRMA 550
      GO TO 3                                                           SRMA 560
C                                                                       SRMA 570
C        OTHERWISE, MULTIPLY ROW BY CONSTANT                            SRMA 580
C                                                                       SRMA 590
    2 A(LAJ)=A(LAJ)*C                                                   SRMA 600
    3 CONTINUE                                                          SRMA 610
      RETURN                                                            SRMA 620
      END                                                               SRMA 630
