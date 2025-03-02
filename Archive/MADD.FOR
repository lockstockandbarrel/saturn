C                                                                       MADD  10
C     ..................................................................MADD  20
C                                                                       MADD  30
C        SUBROUTINE MADD                                                MADD  40
C                                                                       MADD  50
C        PURPOSE                                                        MADD  60
C           ADD TWO MATRICES ELEMENT BY ELEMENT TO FORM RESULTANT       MADD  70
C           MATRIX                                                      MADD  80
C                                                                       MADD  90
C        USAGE                                                          MADD 100
C           CALL MADD(A,B,R,N,M,MSA,MSB)                                MADD 110
C                                                                       MADD 120
C        DESCRIPTION OF PARAMETERS                                      MADD 130
C           A - NAME OF INPUT MATRIX                                    MADD 140
C           B - NAME OF INPUT MATRIX                                    MADD 150
C           R - NAME OF OUTPUT MATRIX                                   MADD 160
C           N - NUMBER OF ROWS IN A,B,R                                 MADD 170
C           M - NUMBER OF COLUMNS IN A,B,R                              MADD 180
C           MSA - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         MADD 190
C                  0 - GENERAL                                          MADD 200
C                  1 - SYMMETRIC                                        MADD 210
C                  2 - DIAGONAL                                         MADD 220
C           MSB - SAME AS MSA EXCEPT FOR MATRIX B                       MADD 230
C                                                                       MADD 240
C        REMARKS                                                        MADD 250
C           NONE                                                        MADD 260
C                                                                       MADD 270
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  MADD 280
C           LOC                                                         MADD 290
C                                                                       MADD 300
C        METHOD                                                         MADD 310
C           STORAGE MODE OF OUTPUT MATRIX IS FIRST DETERMINED. ADDITION MADD 320
C           OF CORRESPONDING ELEMENTS IS THEN PERFORMED.                MADD 330
C           THE FOLLOWING TABLE SHOWS THE STORAGE MODE OF THE OUTPUT    MADD 340
C           MATRIX FOR ALL COMBINATIONS OF INPUT MATRICES               MADD 350
C                         A                B                 R          MADD 360
C                      GENERAL          GENERAL          GENERAL        MADD 370
C                      GENERAL          SYMMETRIC        GENERAL        MADD 380
C                      GENERAL          DIAGONAL         GENERAL        MADD 390
C                      SYMMETRIC        GENERAL          GENERAL        MADD 400
C                      SYMMETRIC        SYMMETRIC        SYMMETRIC      MADD 410
C                      SYMMETRIC        DIAGONAL         SYMMETRIC      MADD 420
C                      DIAGONAL         GENERAL          GENERAL        MADD 430
C                      DIAGONAL         SYMMETRIC        SYMMETRIC      MADD 440
C                      DIAGONAL         DIAGONAL         DIAGONAL       MADD 450
C                                                                       MADD 460
C     ..................................................................MADD 470
C                                                                       MADD 480
      SUBROUTINE MADD(A,B,R,N,M,MSA,MSB)                                MADD 490
      DIMENSION A(1),B(1),R(1)                                          MADD 500
C                                                                       MADD 510
C        DETERMINE STORAGE MODE OF OUTPUT MATRIX                        MADD 520
C                                                                       MADD 530
      IF(MSA-MSB) 7,5,7                                                 MADD 540
    5 CALL LOC(N,M,NM,N,M,MSA)                                          MADD 550
      GO TO 100                                                         MADD 560
    7 MTEST=MSA*MSB                                                     MADD 570
      MSR=0                                                             MADD 580
      IF(MTEST) 20,20,10                                                MADD 590
   10 MSR=1                                                             MADD 600
   20 IF(MTEST-2) 35,35,30                                              MADD 610
   30 MSR=2                                                             MADD 620
C                                                                       MADD 630
C        LOCATE ELEMENTS AND PERFORM ADDITION                           MADD 640
C                                                                       MADD 650
   35 DO 90 J=1,M                                                       MADD 660
      DO 90 I=1,N                                                       MADD 670
      CALL LOC(I,J,IJR,N,M,MSR)                                         MADD 680
      IF(IJR) 40,90,40                                                  MADD 690
   40 CALL LOC(I,J,IJA,N,M,MSA)                                         MADD 700
      AEL=0.0                                                           MADD 710
      IF(IJA) 50,60,50                                                  MADD 720
   50 AEL=A(IJA)                                                        MADD 730
   60 CALL LOC(I,J,IJB,N,M,MSB)                                         MADD 740
      BEL=0.0                                                           MADD 750
      IF(IJB) 70,80,70                                                  MADD 760
   70 BEL=B(IJB)                                                        MADD 770
   80 R(IJR)=AEL+BEL                                                    MADD 780
   90 CONTINUE                                                          MADD 790
      RETURN                                                            MADD 800
C                                                                       MADD 810
C        ADD MATRICES FOR OTHER CASES                                   MADD 820
C                                                                       MADD 830
  100 DO 110 I=1,NM                                                     MADD 840
  110 R(I)=A(I)+B(I)                                                    MADD 850
      RETURN                                                            MADD 860
      END                                                               MADD 870
