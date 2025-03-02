C                                                                       POIN  10
C     ..................................................................POIN  20
C                                                                       POIN  30
C        SUBROUTINE POINT                                               POIN  40
C                                                                       POIN  50
C        PURPOSE                                                        POIN  60
C           TO COMPUTE THE POINT-BISERIAL CORRELATION COEFFICIENT       POIN  70
C           BETWEEN TWO VARIABLES, WHEN ONE OF THE VARIABLES IS A BINARYPOIN  80
C           VARIABLE AND ONE IS CONTINUOUS.  THIS IS A SPECIAL CASE OF  POIN  90
C           THE PEARSON PRODUCT-MOMENT CORRELATION COEFFICIENT.         POIN 100
C                                                                       POIN 110
C        USAGE                                                          POIN 120
C           CALL POINT (N,A,B,HI,ANS,IER)                               POIN 130
C                                                                       POIN 140
C        DESCRIPTION OF PARAMETERS                                      POIN 150
C           N   - NUMBER OF OBSERVATIONS                                POIN 160
C           A   - INPUT VECTOR OF LENGTH N CONTAINING THE CONTINUOUS    POIN 170
C                 VARIABLE                                              POIN 180
C           B   - INPUT VECTOR OF LENGTH N CONTAINING THE DICHOTOMOUS   POIN 190
C                 (BINARY) VARIABLE                                     POIN 200
C           HI  - INPUT NUMERICAL CODE TO INDICATE THE HIGHER CATEGORY. POIN 210
C                 ANY VALUE OF THE BINARY VARIABLE NOT LESS THAN HI WILLPOIN 220
C                 BE CLASSIFIED IN THE HIGHER OF THE TWO CATEGORIES.    POIN 230
C           ANS - OUTPUT VECTOR OF LENGTH 9 CONTAINING THE FOLLOWING    POIN 240
C                 RESULTS                                               POIN 250
C                    ANS(1)- MEAN OF VARIABLE A                         POIN 260
C                    ANS(2)- STANDARD DEVIATION OF VARIABLE A           POIN 270
C                    ANS(3)- NUMBER OF OBSERVATIONS IN THE HIGHER       POIN 280
C                            CATEGORY OF VARIABLE B                     POIN 290
C                    ANS(4)- NUMBER OF OBSERVATIONS IN THE LOWER        POIN 300
C                            CATEGORY OF VARIABLE B                     POIN 310
C                    ANS(5)- MEAN OF VARIABLE A FOR ONLY THOSE          POIN 320
C                            OBSERVATIONS IN THE HIGHER CATEGORY OF     POIN 330
C                            VARIABLE B                                 POIN 340
C                    ANS(6)- MEAN OF VARIABLE A FOR ONLY THOSE          POIN 350
C                            OBSERVATIONS IN THE LOWER CATEGORY OF      POIN 360
C                            VARIABLE B                                 POIN 370
C                    ANS(7)- POINT-BISERIAL CORRELATION COEFFICIENT     POIN 380
C                    ANS(8)- T-TEST FOR THE SIGNIFICANCE OF THE         POIN 390
C                            DIFFERENCE BETWEEN THE MEANS OF VARIABLE A POIN 400
C                            FOR THE HIGHER AND LOWER CATEGORIES        POIN 410
C                            RESPECTIVELY.                              POIN 420
C                    ANS(9)- DEGREES OF FREEDOM FOR THE T-TEST          POIN 430
C           IER- 1, IF ALL ELEMENTS OF B ARE NOT LESS THAN HI.          POIN 440
C                -1, IF ALL ELEMENTS OF B ARE LESS THAN HI.             POIN 450
C                0, OTHERWISE.  IF IER IS NON-ZERO, ANS(I), I=5,...,9,  POIN 460
C                IS SET TO 10**75.                                      POIN 470
C                                                                       POIN 480
C        REMARKS                                                        POIN 490
C           THE SYMBOLS USED TO IDENTFY THE VALUES OF THE TWO CATEGORIESPOIN 500
C           OF VARIABLE B MUST BE NUMERIC.  ALPHABETIC OR SPECIAL       POIN 510
C           CHARACTERS CANNOT BE USED.                                  POIN 520
C           THE T-TEST(ANS(8)) IS A TEST OF WHETHER THE POINT-BISERIAL  POIN 530
C           COEFFICIENT DIFFERS SIGNIFICANTLY FROM ZERO.                POIN 540
C                                                                       POIN 550
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  POIN 560
C           NONE                                                        POIN 570
C                                                                       POIN 580
C        METHOD                                                         POIN 590
C           REFER TO P. HORST, 'PSYCHOLOGICAL MEASUREMENT AND           POIN 600
C           PREDICTION', P. 91 (WADSWORTH, 1966).                       POIN 610
C                                                                       POIN 620
C     ..................................................................POIN 630
C                                                                       POIN 640
      SUBROUTINE POINT (N,A,B,HI,ANS,IER)                               POIN 650
C                                                                       POIN 660
      DIMENSION A(1),B(1),ANS(1)                                        POIN 670
C                                                                       POIN 680
C        COMPUTE MEAN AND STANDARD DEVIATION OF VARIABLE A              POIN 690
C                                                                       POIN 700
      IER=0                                                             POIN 710
      SUM=0.0                                                           POIN 720
      SUM2=0.0                                                          POIN 730
      DO 10 I=1,N                                                       POIN 740
      SUM=SUM+A(I)                                                      POIN 750
   10 SUM2=SUM2+A(I)*A(I)                                               POIN 760
      FN=N                                                              POIN 770
      ANS(1)=SUM/FN                                                     POIN 780
      ANS(2)=(SUM2-ANS(1)*SUM)/(FN-1.0)                                 POIN 790
      ANS(2)= SQRT(ANS(2))                                              POIN 800
C                                                                       POIN 810
C        FIND NUMBERS OF CASES IN THE HIGHER AND LOWER CATEGORIES       POIN 820
C                                                                       POIN 830
      P=0.0                                                             POIN 840
      SUM=0.0                                                           POIN 850
      SUM2=0.0                                                          POIN 860
      DO 30 I=1,N                                                       POIN 870
      IF(B(I)-HI) 20, 25, 25                                            POIN 880
   20 SUM2=SUM2+A(I)                                                    POIN 890
      GO TO 30                                                          POIN 900
   25 P=P+1.0                                                           POIN 910
      SUM=SUM+A(I)                                                      POIN 920
   30 CONTINUE                                                          POIN 930
C                                                                       POIN 940
      Q=FN-P                                                            POIN 950
      ANS(3)=P                                                          POIN 960
      ANS(4)=Q                                                          POIN 970
      IF (P) 35,35,40                                                   POIN 980
   35 IER=-1                                                            POIN 990
      GO TO 50                                                          POIN1000
   40 ANS(5)=SUM/P                                                      POIN1010
      IF (Q) 45,45,60                                                   POIN1020
   45 IER=1                                                             POIN1030
   50 DO 55 I=5,9                                                       POIN1040
   55 ANS(I)=1.E75                                                      POIN1050
      GO TO 65                                                          POIN1060
   60 ANS(6)=SUM2/Q                                                     POIN1070
C                                                                       POIN1080
C        COMPUTE THE POINT-BISERIAL CORRELATION                         POIN1090
C                                                                       POIN1100
      R=((ANS(5)-ANS(1))/ANS(2))* SQRT(P/Q)                             POIN1110
      ANS(7)=R                                                          POIN1120
C                                                                       POIN1130
C        COMPUTE T RATIO USED TO TEST THE HYPOTHESIS OF ZERO CORRELATIONPOIN1140
C                                                                       POIN1150
      T=R* SQRT((FN-2.0)/(1.0-R*R))                                     POIN1160
      ANS(8)=T                                                          POIN1170
C                                                                       POIN1180
C        COMPUTE DEGREES OF FREEDOM                                     POIN1190
C                                                                       POIN1200
      ANS(9)=FN-2                                                       POIN1210
C                                                                       POIN1220
   65 RETURN                                                            POIN1230
      END                                                               POIN1240
