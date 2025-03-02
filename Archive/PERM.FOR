C                                                                       PERM  10
C     ..................................................................PERM  20
C                                                                       PERM  30
C        SUBROUTINE PERM                                                PERM  40
C                                                                       PERM  50
C        PURPOSE                                                        PERM  60
C           TO COMPUTE THE PERMUTATION VECTOR THAT IS INVERSE TO A GIVENPERM  70
C           PERMUTATION VECTOR, THE PERMUTATION VECTOR THAT IS EQUIVA-  PERM  80
C           LENT TO A GIVEN TRANSPOSITION VECTOR AND A TRANSPOSITION    PERM  90
C           VECTOR THAT IS EQUIVALENT TO A GIVEN PERMUTATION VECTOR.    PERM 100
C           (SEE THE GENERAL DISCUSSION FOR DEFINITIONS AND NOTATION.)  PERM 110
C                                                                       PERM 120
C        USAGE                                                          PERM 130
C           CALL PERM(IP1,IP2,N,IPAR,IER)                               PERM 140
C                                                                       PERM 150
C        DESCRIPTION OF PARAMETERS                                      PERM 160
C           IP1  - GIVEN PERMUTATION OR TRANSPOSITION VECTOR            PERM 170
C                  (DIMENSION N)                                        PERM 180
C           IP2  - RESULTING PERMUTATION OR TRANSPOSITION VECTOR        PERM 190
C                  (DIMENSION N)                                        PERM 200
C           N    - DIMENSION OF VECTORS IP1 AND IP2                     PERM 210
C           IPAR - INPUT PARAMETER                                      PERM 220
C                  IPAR NEGATIVE - COMPUTE THE PERMUTATION VECTOR IP2   PERM 230
C                                  THAT IS THE INVERSE OF THE PERMUTA-  PERM 240
C                                  TION VECTOR IP1                      PERM 250
C                  IPAR  =  ZERO - COMPUTE THE PERMUTATION VECTOR IP2   PERM 260
C                                  THAT IS EQUIVALENT TO THE TRANSPOSI- PERM 270
C                                  TION VECTOR IP1                      PERM 280
C                  IPAR POSITIVE - COMPUTE A TRANSPOSITION VECTOR IP2   PERM 290
C                                  THAT IS EQUIVALENT TO THE PERMUTATIONPERM 300
C                                  VECTOR IP1                           PERM 310
C           IER  - RESULTING ERROR PARAMETER                            PERM 320
C                  IER=-1  -  N IS NOT POSITIVE                         PERM 330
C                  IER= 0  -  NO ERROR                                  PERM 340
C                  IER= 1  -  IP1 IS EITHER NOT A PERMUTATION VECTOR OR PERM 350
C                             NOT A TRANSPOSITION VECTOR ON 1,...,N,    PERM 360
C                             DEPENDING ON WHETHER IPAR IS NON-ZERO OR  PERM 370
C                             ZERO, RESPECTIVELY                        PERM 380
C                                                                       PERM 390
C        REMARKS                                                        PERM 400
C           (1)  IF IER=-1 THERE HAS BEEN NO COMPUTATION.               PERM 410
C           (2)  IF IER=1, THEN COMPUTATION HAS BEEN UNSUCCESSFUL DUE TOPERM 420
C                ERROR AND THE PARTIAL RESULTS FOUND IN IP2 ARE USELESS.PERM 430
C           (3)  IP2 CANNOT HAVE THE SAME STORAGE ALLOCATION AS IP1.    PERM 440
C                                                                       PERM 450
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PERM 460
C           NONE                                                        PERM 470
C                                                                       PERM 480
C        METHOD                                                         PERM 490
C           (1)  IPAR NEGATIVE - FOR EACH I, I=1,...,N, IP2(IP1(I)) IS  PERM 500
C                                SET TO I.                              PERM 510
C           (2)  IPAR  =  ZERO - INITIALLY IP2(I) IS SET TO I FOR       PERM 520
C                                I=1,...,N.  THEN, FOR I=1,...,N IN THATPERM 530
C                                ORDER, IP2(I) AND IP2(IP1(I)) ARE      PERM 540
C                                INTERCHANGED.                          PERM 550
C           (3)  IPAR POSITIVE - INITIALLY IP1 IS MOVED TO IP2.  THEN   PERM 560
C                                THE FOLLOWING TWO STEPS ARE REPEATED   PERM 570
C                                FOR I SUCCESSIVELY EQUAL TO 1,...,N.   PERM 580
C                                (A) FIND THE SMALLEST J GREATER THAN ORPERM 590
C                                    EQUAL TO I SUCH THAT IP2(J)=I.     PERM 600
C                                (B) SET IP2(J) TO IP2(I).              PERM 610
C                                                                       PERM 620
C     ..................................................................PERM 630
C                                                                       PERM 640
      SUBROUTINE PERM(IP1,IP2,N,IPAR,IER)                               PERM 650
C                                                                       PERM 660
C                                                                       PERM 670
      DIMENSION IP1(1),IP2(1)                                           PERM 680
C                                                                       PERM 690
C        TEST DIMENSION                                                 PERM 700
      IF(N)19,19,1                                                      PERM 710
C                                                                       PERM 720
C        TEST IPAR TO DETERMINE WHETHER IP1 IS TO BE INTERPRETED AS     PERM 730
C        A PERMUTATION VECTOR OR AS A TRANSPOSITION VECTOR              PERM 740
    1 IF(IPAR)2,13,2                                                    PERM 750
C                                                                       PERM 760
C        CHECK THAT IP1 IS A PERMUTATION VECTOR AND COMPUTE IP1 INVERSE PERM 770
    2 DO 3 I=1,N                                                        PERM 780
    3 IP2(I)=0                                                          PERM 790
      DO 6 I=1,N                                                        PERM 800
      K=IP1(I)                                                          PERM 810
      IF(K-N)4,5,20                                                     PERM 820
    4 IF(K)20,20,5                                                      PERM 830
    5 IF(IP2(K))20,6,20                                                 PERM 840
    6 IP2(K)=I                                                          PERM 850
C                                                                       PERM 860
C        TEST IPAR FOR THE DESIRED OPERATION                            PERM 870
      IF(IPAR)12,7,7                                                    PERM 880
C                                                                       PERM 890
C        COMPUTE TRANSPOSITION VECTOR IP2 FOR PERMUTATION VECTOR IP1    PERM 900
    7 DO 8 I=1,N                                                        PERM 910
    8 IP2(I)=IP1(I)                                                     PERM 920
      NN=N-1                                                            PERM 930
      IF(NN)12,12,9                                                     PERM 940
    9 DO 11 I=1,NN                                                      PERM 950
      DO 10 J=1,NN                                                      PERM 960
      IF(IP2(J)-I)10,11,10                                              PERM 970
   10 CONTINUE                                                          PERM 980
      J=N                                                               PERM 990
   11 IP2(J)=IP2(I)                                                     PERM1000
C                                                                       PERM1010
C        NORMAL RETURN - NO ERROR                                       PERM1020
   12 IER=0                                                             PERM1030
      RETURN                                                            PERM1040
C                                                                       PERM1050
C        COMPUTE PERMUTATION VECTOR IP2 FOR TRANSPOSITION VECTOR IP1    PERM1060
   13 DO 14 I=1,N                                                       PERM1070
   14 IP2(I)=I                                                          PERM1080
      DO 18 I=1,N                                                       PERM1090
      K=IP1(I)                                                          PERM1100
      IF(K-I)15,18,16                                                   PERM1110
   15 IF(K)20,20,17                                                     PERM1120
   16 IF(N-K)20,17,17                                                   PERM1130
   17 J=IP2(I)                                                          PERM1140
      IP2(I)=IP2(K)                                                     PERM1150
      IP2(K)=J                                                          PERM1160
   18 CONTINUE                                                          PERM1170
      GO TO 12                                                          PERM1180
C                                                                       PERM1190
C        ERROR RETURN - N IS NOT POSITIVE                               PERM1200
   19 IER=-1                                                            PERM1210
      RETURN                                                            PERM1220
C                                                                       PERM1230
C        ERROR RETURN - IP1 IS EITHER NOT A PERMUTATION VECTOR          PERM1240
C                       OR NOT A TRANSPOSITION VECTOR                   PERM1250
   20 IER=1                                                             PERM1260
      RETURN                                                            PERM1270
      END                                                               PERM1280
