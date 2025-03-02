C                                                                       MSUB  10
C     ..................................................................MSUB  20
C                                                                       MSUB  30
C        SUBROUTINE MSUB                                                MSUB  40
C                                                                       MSUB  50
C        PURPOSE                                                        MSUB  60
C           SUBTRACT TWO MATRICES ELEMENT BY ELEMENT TO FORM RESULTANT  MSUB  70
C           MATRIX                                                      MSUB  80
C                                                                       MSUB  90
C        USAGE                                                          MSUB 100
C           CALL MSUB(A,B,R,N,M,MSA,MSB)                                MSUB 110
C                                                                       MSUB 120
C        DESCRIPTION OF PARAMETERS                                      MSUB 130
C           A - NAME OF INPUT MATRIX                                    MSUB 140
C           B - NAME OF INPUT MATRIX                                    MSUB 150
C           R - NAME OF OUTPUT MATRIX                                   MSUB 160
C           N - NUMBER OF ROWS IN A,B,R                                 MSUB 170
C           M - NUMBER OF COLUMNS IN A,B,R                              MSUB 180
C           MSA - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         MSUB 190
C                  0 - GENERAL                                          MSUB 200
C                  1 - SYMMETRIC                                        MSUB 210
C                  2 - DIAGONAL                                         MSUB 220
C           MSB - SAME AS MSA EXCEPT FOR MATRIX B                       MSUB 230
C                                                                       MSUB 240
C        REMARKS                                                        MSUB 250
C           NONE                                                        MSUB 260
C                                                                       MSUB 270
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  MSUB 280
C           LOC                                                         MSUB 290
C                                                                       MSUB 300
C        METHOD                                                         MSUB 310
C           STRUCTURE OF OUTPUT MATRIX IS FIRST DETERMINED. SUBTRACTION MSUB 320
C           OF MATRIX B ELEMENTS FROM CORRESPONDING MATRIX A ELEMENTS   MSUB 330
C           IS THEN PERFORMED.                                          MSUB 340
C           THE FOLLOWING TABLE SHOWS THE STORAGE MODE OF THE OUTPUT    MSUB 350
C           MATRIX FOR ALL COMBINATIONS OF INPUT MATRICES               MSUB 360
C                         A                B                 R          MSUB 370
C                      GENERAL          GENERAL          GENERAL        MSUB 380
C                      GENERAL          SYMMETRIC        GENERAL        MSUB 390
C                      GENERAL          DIAGONAL         GENERAL        MSUB 400
C                      SYMMETRIC        GENERAL          GENERAL        MSUB 410
C                      SYMMETRIC        SYMMETRIC        SYMMETRIC      MSUB 420
C                      SYMMETRIC        DIAGONAL         SYMMETRIC      MSUB 430
C                      DIAGONAL         GENERAL          GENERAL        MSUB 440
C                      DIAGONAL         SYMMETRIC        SYMMETRIC      MSUB 450
C                      DIAGONAL         DIAGONAL         DIAGONAL       MSUB 460
C                                                                       MSUB 470
C     ..................................................................MSUB 480
C                                                                       MSUB 490
      SUBROUTINE MSUB(A,B,R,N,M,MSA,MSB)                                MSUB 500
      DIMENSION A(1),B(1),R(1)                                          MSUB 510
C                                                                       MSUB 520
C        DETERMINE STORAGE MODE OF OUTPUT MATRIX                        MSUB 530
C                                                                       MSUB 540
      IF(MSA-MSB) 7,5,7                                                 MSUB 550
    5 CALL LOC(N,M,NM,N,M,MSA)                                          MSUB 560
      GO TO 100                                                         MSUB 570
    7 MTEST=MSA*MSB                                                     MSUB 580
      MSR=0                                                             MSUB 590
      IF(MTEST) 20,20,10                                                MSUB 600
   10 MSR=1                                                             MSUB 610
   20 IF(MTEST-2) 35,35,30                                              MSUB 620
   30 MSR=2                                                             MSUB 630
C                                                                       MSUB 640
C        LOCATE ELEMENTS AND PERFORM SUBTRACTION                        MSUB 650
C                                                                       MSUB 660
   35 DO 90 J=1,M                                                       MSUB 670
      DO 90 I=1,N                                                       MSUB 680
      CALL LOC(I,J,IJR,N,M,MSR)                                         MSUB 690
      IF(IJR) 40,90,40                                                  MSUB 700
   40 CALL LOC(I,J,IJA,N,M,MSA)                                         MSUB 710
      AEL=0.0                                                           MSUB 720
      IF(IJA) 50,60,50                                                  MSUB 730
   50 AEL=A(IJA)                                                        MSUB 740
   60 CALL LOC(I,J,IJB,N,M,MSB)                                         MSUB 750
      BEL=0.0                                                           MSUB 760
      IF(IJB) 70,80,70                                                  MSUB 770
   70 BEL=B(IJB)                                                        MSUB 780
   80 R(IJR)=AEL-BEL                                                    MSUB 790
   90 CONTINUE                                                          MSUB 800
      RETURN                                                            MSUB 810
C                                                                       MSUB 820
C        SUBTRACT MATRICES FOR OTHER CASES                              MSUB 830
C                                                                       MSUB 840
  100 DO 110 I=1,NM                                                     MSUB 850
  110 R(I)=A(I)-B(I)                                                    MSUB 860
      RETURN                                                            MSUB 870
      END                                                               MSUB 880
