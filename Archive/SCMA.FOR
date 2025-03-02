C                                                                       SCMA  10
C     ..................................................................SCMA  20
C                                                                       SCMA  30
C        SUBROUTINE SCMA                                                SCMA  40
C                                                                       SCMA  50
C        PURPOSE                                                        SCMA  60
C           MULTIPLY COLUMN OF MATRIX BY A SCALAR AND ADD TO ANOTHER    SCMA  70
C           COLUMN OF THE SAME MATRIX                                   SCMA  80
C                                                                       SCMA  90
C        USAGE                                                          SCMA 100
C           CALL SCMA(A,C,N,LA,LB)                                      SCMA 110
C                                                                       SCMA 120
C        DESCRIPTION OF PARAMETERS                                      SCMA 130
C           A  - NAME OF MATRIX                                         SCMA 140
C           C  - SCALAR                                                 SCMA 150
C           N  - NUMBER OF ROWS IN A                                    SCMA 160
C           LA - COLUMN IN A TO BE MULTIPLIED BY SCALAR                 SCMA 170
C           LB - COLUMN IN A TO WHICH PRODUCT IS ADDED                  SCMA 180
C                IF 0 IS SPECIFIED, PRODUCT REPLACES ELEMENTS IN LA     SCMA 190
C                                                                       SCMA 200
C        REMARKS                                                        SCMA 210
C           MATRIX A MUST BE A GENERAL MATRIX                           SCMA 220
C                                                                       SCMA 230
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  SCMA 240
C           NONE                                                        SCMA 250
C                                                                       SCMA 260
C        METHOD                                                         SCMA 270
C           EACH ELEMENT OF COLUMN LA IS MULTIPLIED BY SCALAR C AND THE SCMA 280
C           PRODUCT IS ADDED TO THE CORRESPONDING ELEMENT OF COLUMN LB. SCMA 290
C           COLUMN LA REMAINS UNAFFECTED BY THE OPERATION.              SCMA 300
C           IF PARAMETER LB CONTAINS ZERO, MULTIPLICATION BY THE SCALAR SCMA 310
C           IS PERFORMED AND THE PRODUCT REPLACES ELEMENTS IN LA.       SCMA 320
C                                                                       SCMA 330
C     ..................................................................SCMA 340
C                                                                       SCMA 350
      SUBROUTINE SCMA(A,C,N,LA,LB)                                      SCMA 360
      DIMENSION A(1)                                                    SCMA 370
C                                                                       SCMA 380
C        LOCATE STARTING POINT OF BOTH COLUMNS                          SCMA 390
C                                                                       SCMA 400
      ILA=N*(LA-1)                                                      SCMA 410
      ILB=N*(LB-1)                                                      SCMA 420
C                                                                       SCMA 430
      DO 3 I=1,N                                                        SCMA 440
      ILA=ILA+1                                                         SCMA 450
      ILB=ILB+1                                                         SCMA 460
C                                                                       SCMA 470
C        CHECK LB FOR ZERO                                              SCMA 480
C                                                                       SCMA 490
      IF(LB) 1,2,1                                                      SCMA 500
C                                                                       SCMA 510
C        IF NOT MULTIPLY BY CONSTANT AND ADD TO SECOND COLUMN           SCMA 520
C                                                                       SCMA 530
    1 A(ILB)=A(ILA)*C+A(ILB)                                            SCMA 540
      GO TO 3                                                           SCMA 550
C                                                                       SCMA 560
C        OTHERWISE, MULTIPLY COLUMN BY CONSTANT                         SCMA 570
C                                                                       SCMA 580
    2 A(ILA)=A(ILA)*C                                                   SCMA 590
    3 CONTINUE                                                          SCMA 600
      RETURN                                                            SCMA 610
      END                                                               SCMA 620
