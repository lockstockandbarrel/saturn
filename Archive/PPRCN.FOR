C                                                                       PPRC  10
C     ..................................................................PPRC  20
C                                                                       PPRC  30
C        SUBROUTINE PPRCN                                               PPRC  40
C                                                                       PPRC  50
C        PURPOSE                                                        PPRC  60
C           TO COMPUTE, GIVEN TWO PERMUTATION VECTORS IP1 AND IP2, THE  PPRC  70
C           COMPOSITION IP2(IP1) AND THE CONJUGATE IP1(IP2(IP1 INVERSE))PPRC  80
C           OF IP2 BY IP1.  (SEE THE GENERAL DISCUSSION FOR DEFINITIONS PPRC  90
C           AND NOTATION.)                                              PPRC 100
C                                                                       PPRC 110
C        USAGE                                                          PPRC 120
C           CALL PPRCN(IP1,IP2,IP3,N,IPAR,IER)                          PPRC 130
C                                                                       PPRC 140
C        DESCRIPTION OF PARAMETERS                                      PPRC 150
C           IP1  - GIVEN PERMUTATION VECTOR (DIMENSION N)               PPRC 160
C           IP2  - GIVEN PERMUTATION VECTOR (DIMENSION N)               PPRC 170
C           IP3  - RESULTING PERMUTATION VECTOR (DIMENSION N)           PPRC 180
C           N    - DIMENSION OF VECTORS IP1, IP2 AND IP3                PPRC 190
C           IPAR - INPUT PARAMETER                                      PPRC 200
C                  IPAR NON-NEGATIVE - COMPUTE IP2(IP1)                 PPRC 210
C                  IPAR NEGATIVE     - COMPUTE IP1(IP2(IP1 INVERSE))    PPRC 220
C           IER  - RESULTING ERROR PARAMETER                            PPRC 230
C                  IER=-1  -  N IS NOT POSITIVE                         PPRC 240
C                  IER= 0  -  NO ERROR                                  PPRC 250
C                  IER= 1  -  IP1 AND IP2 ARE NOT BOTH PERMUTATION      PPRC 260
C                             VECTORS ON 1,...,N                        PPRC 270
C                                                                       PPRC 280
C        REMARKS                                                        PPRC 290
C           (1)  IF IER=-1 THERE HAS BEEN NO COMPUTATION.               PPRC 300
C           (2)  IF IER=1, THEN COMPUTATION HAS BEEN UNSUCCESSFUL DUE TOPPRC 310
C                ERROR AND THE PARTIAL RESULTS FOUND IN IP2 ARE USELESS.PPRC 320
C           (3)  IP3 CANNOT HAVE THE SAME STORAGE ALLOCATION AS IP1 OR  PPRC 330
C                IP2.                                                   PPRC 340
C                                                                       PPRC 350
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PPRC 360
C           PERM                                                        PPRC 370
C                                                                       PPRC 380
C        METHOD                                                         PPRC 390
C           SUBROUTINE PERM IS USED TO CHECK THAT IP1 AND IP2 ARE PERMU-PPRC 400
C           TATION VECTORS.  IF IP2(IP1) IS COMPUTED, IP3(I) IS SET TO  PPRC 410
C           IP2(IP1(I)) FOR I=1,...,N.  IF IP1(IP2(IP1 INVERSE)) IS     PPRC 420
C           COMPUTED, FIRST IP3 IS SET TO IP1 INVERSE BY SUBROUTINE PERMPPRC 430
C           AND THEN IP3(I) IS SET TO IP1(IP2(IP3(I))) FOR I=1,...,N.   PPRC 440
C                                                                       PPRC 450
C     ..................................................................PPRC 460
C                                                                       PPRC 470
      SUBROUTINE PPRCN(IP1,IP2,IP3,N,IPAR,IER)                          PPRC 480
C                                                                       PPRC 490
C                                                                       PPRC 500
      DIMENSION IP1(1),IP2(1),IP3(1)                                    PPRC 510
C                                                                       PPRC 520
C        CHECK THAT N IS POSITIVE AND THAT IP2 IS A PERMUTATION VECTOR  PPRC 530
      CALL PERM(IP2,IP3,N,-1,IER)                                       PPRC 540
C                                                                       PPRC 550
C        TEST IER TO SEE IF THERE IS AN ERROR                           PPRC 560
      IF(IER)7,1,7                                                      PPRC 570
C                                                                       PPRC 580
C        CHECK THAT IP1 IS A PERMUTATION VECTOR AND COMPUTE IP1 INVERSE PPRC 590
    1 CALL PERM(IP1,IP3,N,-1,IER)                                       PPRC 600
C                                                                       PPRC 610
C        TEST IER TO SEE IF THERE IS AN ERROR                           PPRC 620
      IF(IER)7,2,7                                                      PPRC 630
C                                                                       PPRC 640
C        TEST IPAR FOR THE DESIRED OPERATION                            PPRC 650
    2 IF(IPAR)3,5,5                                                     PPRC 660
C                                                                       PPRC 670
C        COMPUTE IP1(IP2(IP1 INVERSE))                                  PPRC 680
    3 DO 4 I=1,N                                                        PPRC 690
      K=IP3(I)                                                          PPRC 700
      J=IP2(K)                                                          PPRC 710
    4 IP3(I)=IP1(J)                                                     PPRC 720
      RETURN                                                            PPRC 730
C                                                                       PPRC 740
C        COMPUTE IP2(IP1)                                               PPRC 750
    5 DO 6 I=1,N                                                        PPRC 760
      K=IP1(I)                                                          PPRC 770
    6 IP3(I)=IP2(K)                                                     PPRC 780
    7 RETURN                                                            PPRC 790
      END                                                               PPRC 800
