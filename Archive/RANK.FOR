C                                                                       RANK  10
C     ..................................................................RANK  20
C                                                                       RANK  30
C        SUBROUTINE RANK                                                RANK  40
C                                                                       RANK  50
C        PURPOSE                                                        RANK  60
C           RANK A VECTOR OF VALUES                                     RANK  70
C                                                                       RANK  80
C        USAGE                                                          RANK  90
C           CALL RANK(A,R,N)                                            RANK 100
C                                                                       RANK 110
C        DESCRIPTION OF PARAMETERS                                      RANK 120
C           A - INPUT VECTOR OF N VALUES                                RANK 130
C           R - OUTPUT VECTOR OF LENGTH N. SMALLEST VALUE IS RANKED 1,  RANK 140
C               LARGEST IS RANKED N. TIES ARE ASSIGNED AVERAGE OF TIED  RANK 150
C               RANKS                                                   RANK 160
C           N - NUMBER OF VALUES                                        RANK 170
C                                                                       RANK 180
C        REMARKS                                                        RANK 190
C           NONE                                                        RANK 200
C                                                                       RANK 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  RANK 220
C           NONE                                                        RANK 230
C                                                                       RANK 240
C        METHOD                                                         RANK 250
C           VECTOR IS SEARCHED FOR SUCCESSIVELY LARGER ELEMENTS. IF TIESRANK 260
C           OCCUR, THEY ARE LOCATED AND THEIR RANK VALUE COMPUTED.      RANK 270
C           FOR EXAMPLE, IF 2 VALUES ARE TIED FOR SIXTH RANK, THEY ARE  RANK 280
C           ASSIGNED A RANK OF 6.5 (=(6+7)/2)                           RANK 290
C                                                                       RANK 300
C     ..................................................................RANK 310
C                                                                       RANK 320
      SUBROUTINE RANK(A,R,N)                                            RANK 330
      DIMENSION A(1),R(1)                                               RANK 340
C                                                                       RANK 350
C        INITIALIZATION                                                 RANK 360
C                                                                       RANK 370
      DO 10 I=1,N                                                       RANK 380
   10 R(I)=0.0                                                          RANK 390
C                                                                       RANK 400
C        FIND RANK OF DATA                                              RANK 410
C                                                                       RANK 420
      DO 100 I=1,N                                                      RANK 430
C                                                                       RANK 440
C        TEST WHETHER DATA POINT IS ALREADY RANKED                      RANK 450
C                                                                       RANK 460
      IF(R(I)) 20, 20, 100                                              RANK 470
C                                                                       RANK 480
C        DATA POINT TO BE RANKED                                        RANK 490
C                                                                       RANK 500
   20 SMALL=0.0                                                         RANK 510
      EQUAL=0.0                                                         RANK 520
      X=A(I)                                                            RANK 530
      DO 50 J=1,N                                                       RANK 540
      IF(A(J)-X) 30, 40, 50                                             RANK 550
C        COUNT NUMBER OF DATA POINTS WHICH ARE SMALLER                  RANK 560
C                                                                       RANK 570
C                                                                       RANK 580
   30 SMALL=SMALL+1.0                                                   RANK 590
      GO TO 50                                                          RANK 600
C                                                                       RANK 610
C        COUNT NUMBER OF DATA POINTS WHICH ARE EQUAL                    RANK 620
C                                                                       RANK 630
   40 EQUAL=EQUAL+1.0                                                   RANK 640
      R(J)=-1.0                                                         RANK 650
   50 CONTINUE                                                          RANK 660
C                                                                       RANK 670
C        TEST FOR TIE                                                   RANK 680
C                                                                       RANK 690
      IF(EQUAL-1.0) 60, 60, 70                                          RANK 700
C                                                                       RANK 710
C        STORE RANK OF DATA POINT WHERE NO TIE                          RANK 720
C                                                                       RANK 730
   60 R(I)=SMALL+1.0                                                    RANK 740
      GO TO 100                                                         RANK 750
C                                                                       RANK 760
C        CALCULATE RANK OF TIED DATA POINTS                             RANK 770
C                                                                       RANK 780
   70 P=SMALL + (EQUAL + 1.0)*0.5                                       RANK 790
      DO 90 J=I,N                                                       RANK 800
      IF(R(J)+1.0) 90, 80, 90                                           RANK 810
   80 R(J)=P                                                            RANK 820
   90 CONTINUE                                                          RANK 830
  100 CONTINUE                                                          RANK 840
      RETURN                                                            RANK 850
      END                                                               RANK 860
