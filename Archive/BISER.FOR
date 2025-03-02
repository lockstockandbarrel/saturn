C                                                                       BISE  10
C     ..................................................................BISE  20
C                                                                       BISE  30
C        SUBROUTINE BISER                                               BISE  40
C                                                                       BISE  50
C        PURPOSE                                                        BISE  60
C           TO COMPUTE THE BISERIAL CORRELATION COEFFICIENT BETWEEN TWO BISE  70
C           CONTINUOUS VARIABLES WHEN ONE OF THEM HAS BEEN ARTIFICIALLY BISE  80
C           DICHOTOMIZED.                                               BISE  90
C                                                                       BISE 100
C        USAGE                                                          BISE 110
C           CALL BISER (N,A,B,HI,ANS,IER)                               BISE 120
C                                                                       BISE 130
C        DESCRIPTION OF PARAMETERS                                      BISE 140
C           N   - NUMBER OF OBSERVATIONS                                BISE 150
C           A   - INPUT VECTOR OF LENGTH N CONTAINING THE CONTINUOUS    BISE 160
C                 VARIABLE                                              BISE 170
C           B   - INPUT VECTOR OF LENGTH N CONTAINING THE DICHOTOMIZED  BISE 180
C                 VARIABLE                                              BISE 190
C           HI  - INPUT - NUMERICAL CODE TO INDICATE THE HIGHER CATEGORYBISE 200
C                 OF THE DICHOTOMIZED VARIABLE.  ANY VALUE IN VECTOR B  BISE 210
C                 EQUAL TO OR GREATER THAN HI WILL BE CLASSIFIED INTO   BISE 220
C                 THE HIGHER CATEGORY.                                  BISE 230
C           ANS - OUTPUT VECTOR OF LENGTH 8 CONTAINING THE FOLLOWING    BISE 240
C                 ANS(1) - MEAN OF VARIABLE A                           BISE 250
C                 ANS(2) - STANDARD DEVIATION OF VARIABLE A             BISE 260
C                 ANS(3) - PROPORTION OF THE CASES IN THE HIGHER        BISE 270
C                          CATEGORY OF VARIABLE B                       BISE 280
C                 ANS(4) - PROPORTION OF THE CASES IN THE LOWER         BISE 290
C                          CATEGORY OF VARIABLE B                       BISE 300
C                 ANS(5) - MEAN OF VARIABLE A FOR THOSE CASES FALLING   BISE 310
C                          INTO THE HIGHER CATEGORY OF VARIABLE B       BISE 320
C                 ANS(6) - MEAN OF VARIABLE A FOR THOSE CASES FALLING   BISE 330
C                          INTO THE LOWER CATEGORY OF VARIABLE B        BISE 340
C                 ANS(7) - BISERIAL CORRELATION COEFFICIENT             BISE 350
C                 ANS(8) - STANDARD ERROR OF BISERIAL CORRELATION       BISE 360
C                          COEFFICIENT                                  BISE 370
C           IER -  1, IF NO CASES ARE IN THE LOWER CATEGORY OF VARIABLE BISE 380
C                 B.                                                    BISE 390
C                 -1, IF ALL CASES ARE IN THE LOWER CATEGORY OF         BISE 400
C                 VARIABLE B.                                           BISE 410
C                 0, OTHERWISE.                                         BISE 420
C                 IF IER IS NON-ZERO, ANS(I)=10**75,I=5,...,8.          BISE 430
C                                                                       BISE 440
C        REMARKS                                                        BISE 450
C           THE VALUES OF THE DICHOTOMIZED VARIABLE, B, MUST BE IN      BISE 460
C           NUMERIC FORM.  THEY CANNOR BE SPECIFIED BY MEANS OF         BISE 470
C           ALPHABETIC OR SPECIAL CHARACTERS.                           BISE 480
C                                                                       BISE 490
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  BISE 500
C           NDTRI                                                       BISE 510
C                                                                       BISE 520
C        METHOD                                                         BISE 530
C           REFER TO P. HORST, 'PSYCHOLOGICAL MEASUREMENT AND           BISE 540
C           PREDICTION', P.95-96 (WADSWORTH, 1966).                     BISE 550
C                                                                       BISE 560
C     ..................................................................BISE 570
C                                                                       BISE 580
      SUBROUTINE BISER (N,A,B,HI,ANS,IER)                               BISE 590
C                                                                       BISE 600
      DIMENSION A(1),B(1),ANS(1)                                        BISE 610
C                                                                       BISE 620
C        COMPUTE MEAN AND STANDARD DEVIATION OF VARIABLE A              BISE 630
C                                                                       BISE 640
      IER=0                                                             BISE 650
      SUM=0.0                                                           BISE 660
      SUM2=0.0                                                          BISE 670
      DO 10 I=1,N                                                       BISE 680
      SUM=SUM+A(I)                                                      BISE 690
   10 SUM2=SUM2+A(I)*A(I)                                               BISE 700
      FN=N                                                              BISE 710
      ANS(1)=SUM/FN                                                     BISE 720
      ANS(2)=(SUM2-ANS(1)*SUM)/(FN-1.0)                                 BISE 730
      ANS(2)= SQRT(ANS(2))                                              BISE 740
C                                                                       BISE 750
C        FIND PROPORTIONS OF CASES IN THE HIGHER AND LOWER CATEGORIES   BISE 760
C                                                                       BISE 770
      P=0.0                                                             BISE 780
      SUM=0.0                                                           BISE 790
      SUM2=0.0                                                          BISE 800
      DO 30 I=1,N                                                       BISE 810
      IF(B(I)-HI) 20, 25, 25                                            BISE 820
   20 SUM2=SUM2+A(I)                                                    BISE 830
      GO TO 30                                                          BISE 840
   25 P=P+1.0                                                           BISE 850
      SUM=SUM+A(I)                                                      BISE 860
   30 CONTINUE                                                          BISE 870
      ANS(4)=1.0                                                        BISE 880
      ANS(3)=0.0                                                        BISE 890
      Q=FN-P                                                            BISE 900
      IF (P) 35,35,40                                                   BISE 910
   35 IER=-1                                                            BISE 920
      GO TO 50                                                          BISE 930
   40 ANS(5)=SUM/P                                                      BISE 940
      IF (Q) 45,45,60                                                   BISE 950
   45 IER=1                                                             BISE 960
      ANS(4)=0.0                                                        BISE 970
      ANS(3)=1.0                                                        BISE 980
   50 DO 55 I=5,8                                                       BISE 990
   55 ANS(I)=1.E75                                                      BISE1000
      GO TO 65                                                          BISE1010
   60 ANS(6)=SUM2/Q                                                     BISE1020
      P=P/FN                                                            BISE1030
      Q=1.0-P                                                           BISE1040
C                                                                       BISE1050
C        FIND ORDINATE OF THE NORMAL DISTRIBUTION CURVE AT THE POINT OF BISE1060
C        DIVISION BETWEEN SEGMENTS CONTAINING P AND Q PROPORTIONS       BISE1070
C                                                                       BISE1080
      CALL NDTRI (Q,X,Y,ER)                                             BISE1090
C                                                                       BISE1100
C        COMPUTE THE BISERIAL COEFFICIENT OF CORRELATION                BISE1110
C                                                                       BISE1120
      R=((ANS(5)-ANS(1))/ANS(2))*(P/Y)                                  BISE1130
C                                                                       BISE1140
C        COMPUTE THE STANDARD ERROR OF R                                BISE1150
C                                                                       BISE1160
      ANS(8)=( SQRT(P*Q)/Y-R*R)/SQRT(FN)                                BISE1170
C                                                                       BISE1180
C        STORE RESULTS                                                  BISE1190
C                                                                       BISE1200
      ANS(3)=P                                                          BISE1210
      ANS(4)=Q                                                          BISE1220
      ANS(7)=R                                                          BISE1230
C                                                                       BISE1240
   65 RETURN                                                            BISE1250
      END                                                               BISE1260
C                                                                       BISE  10
C     ..................................................................BISE  20
C                                                                       BISE  30
C        SUBROUTINE BISER                                               BISE  40
C                                                                       BISE  50
C        PURPOSE                                                        BISE  60
C           TO COMPUTE THE BISERIAL CORRELATION COEFFICIENT BETWEEN TWO BISE  70
C           CONTINUOUS VARIABLES WHEN ONE OF THEM HAS BEEN ARTIFICIALLY BISE  80
C           DICHOTOMIZED.                                               BISE  90
C                                                                       BISE 100
C        USAGE                                                          BISE 110
C           CALL BISER (N,A,B,HI,ANS,IER)                               BISE 120
C                                                                       BISE 130
C        DESCRIPTION OF PARAMETERS                                      BISE 140
C           N   - NUMBER OF OBSERVATIONS                                BISE 150
C           A   - INPUT VECTOR OF LENGTH N CONTAINING THE CONTINUOUS    BISE 160
C                 VARIABLE                                              BISE 170
C           B   - INPUT VECTOR OF LENGTH N CONTAINING THE DICHOTOMIZED  BISE 180
C                 VARIABLE                                              BISE 190
C           HI  - INPUT - NUMERICAL CODE TO INDICATE THE HIGHER CATEGORYBISE 200
C                 OF THE DICHOTOMIZED VARIABLE.  ANY VALUE IN VECTOR B  BISE 210
C                 EQUAL TO OR GREATER THAN HI WILL BE CLASSIFIED INTO   BISE 220
C                 THE HIGHER CATEGORY.                                  BISE 230
C           ANS - OUTPUT VECTOR OF LENGTH 8 CONTAINING THE FOLLOWING    BISE 240
C                 ANS(1) - MEAN OF VARIABLE A                           BISE 250
C                 ANS(2) - STANDARD DEVIATION OF VARIABLE A             BISE 260
C                 ANS(3) - PROPORTION OF THE CASES IN THE HIGHER        BISE 270
C                          CATEGORY OF VARIABLE B                       BISE 280
C                 ANS(4) - PROPORTION OF THE CASES IN THE LOWER         BISE 290
C                          CATEGORY OF VARIABLE B                       BISE 300
C                 ANS(5) - MEAN OF VARIABLE A FOR THOSE CASES FALLING   BISE 310
C                          INTO THE HIGHER CATEGORY OF VARIABLE B       BISE 320
C                 ANS(6) - MEAN OF VARIABLE A FOR THOSE CASES FALLING   BISE 330
C                          INTO THE LOWER CATEGORY OF VARIABLE B        BISE 340
C                 ANS(7) - BISERIAL CORRELATION COEFFICIENT             BISE 350
C                 ANS(8) - STANDARD ERROR OF BISERIAL CORRELATION       BISE 360
C                          COEFFICIENT                                  BISE 370
C           IER -  1, IF NO CASES ARE IN THE LOWER CATEGORY OF VARIABLE BISE 380
C                 B.                                                    BISE 390
C                 -1, IF ALL CASES ARE IN THE LOWER CATEGORY OF         BISE 400
C                 VARIABLE B.                                           BISE 410
C                 0, OTHERWISE.                                         BISE 420
C                 IF IER IS NON-ZERO, ANS(I)=10**75,I=5,...,8.          BISE 430
C                                                                       BISE 440
C        REMARKS                                                        BISE 450
C           THE VALUES OF THE DICHOTOMIZED VARIABLE, B, MUST BE IN      BISE 460
C           NUMERIC FORM.  THEY CANNOR BE SPECIFIED BY MEANS OF         BISE 470
C           ALPHABETIC OR SPECIAL CHARACTERS.                           BISE 480
C                                                                       BISE 490
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  BISE 500
C           NDTRI                                                       BISE 510
C                                                                       BISE 520
C        METHOD                                                         BISE 530
C           REFER TO P. HORST, 'PSYCHOLOGICAL MEASUREMENT AND           BISE 540
C           PREDICTION', P.95-96 (WADSWORTH, 1966).                     BISE 550
C                                                                       BISE 560
C     ..................................................................BISE 570
C                                                                       BISE 580
      SUBROUTINE BISER (N,A,B,HI,ANS,IER)                               BISE 590
C                                                                       BISE 600
      DIMENSION A(1),B(1),ANS(1)                                        BISE 610
C                                                                       BISE 620
C        COMPUTE MEAN AND STANDARD DEVIATION OF VARIABLE A              BISE 630
C                                                                       BISE 640
      IER=0                                                             BISE 650
      SUM=0.0                                                           BISE 660
      SUM2=0.0                                                          BISE 670
      DO 10 I=1,N                                                       BISE 680
      SUM=SUM+A(I)                                                      BISE 690
   10 SUM2=SUM2+A(I)*A(I)                                               BISE 700
      FN=N                                                              BISE 710
      ANS(1)=SUM/FN                                                     BISE 720
      ANS(2)=(SUM2-ANS(1)*SUM)/(FN-1.0)                                 BISE 730
      ANS(2)= SQRT(ANS(2))                                              BISE 740
C                                                                       BISE 750
C        FIND PROPORTIONS OF CASES IN THE HIGHER AND LOWER CATEGORIES   BISE 760
C                                                                       BISE 770
      P=0.0                                                             BISE 780
      SUM=0.0                                                           BISE 790
      SUM2=0.0                                                          BISE 800
      DO 30 I=1,N                                                       BISE 810
      IF(B(I)-HI) 20, 25, 25                                            BISE 820
   20 SUM2=SUM2+A(I)                                                    BISE 830
      GO TO 30                                                          BISE 840
   25 P=P+1.0                                                           BISE 850
      SUM=SUM+A(I)                                                      BISE 860
   30 CONTINUE                                                          BISE 870
      ANS(4)=1.0                                                        BISE 880
      ANS(3)=0.0                                                        BISE 890
      Q=FN-P                                                            BISE 900
      IF (P) 35,35,40                                                   BISE 910
   35 IER=-1                                                            BISE 920
      GO TO 50                                                          BISE 930
   40 ANS(5)=SUM/P                                                      BISE 940
      IF (Q) 45,45,60                                                   BISE 950
   45 IER=1                                                             BISE 960
      ANS(4)=0.0                                                        BISE 970
      ANS(3)=1.0                                                        BISE 980
   50 DO 55 I=5,8                                                       BISE 990
   55 ANS(I)=1.E75                                                      BISE1000
      GO TO 65                                                          BISE1010
   60 ANS(6)=SUM2/Q                                                     BISE1020
      P=P/FN                                                            BISE1030
      Q=1.0-P                                                           BISE1040
C                                                                       BISE1050
C        FIND ORDINATE OF THE NORMAL DISTRIBUTION CURVE AT THE POINT OF BISE1060
C        DIVISION BETWEEN SEGMENTS CONTAINING P AND Q PROPORTIONS       BISE1070
C                                                                       BISE1080
      CALL NDTRI (Q,X,Y,ER)                                             BISE1090
C                                                                       BISE1100
C        COMPUTE THE BISERIAL COEFFICIENT OF CORRELATION                BISE1110
C                                                                       BISE1120
      R=((ANS(5)-ANS(1))/ANS(2))*(P/Y)                                  BISE1130
C                                                                       BISE1140
C        COMPUTE THE STANDARD ERROR OF R                                BISE1150
C                                                                       BISE1160
      ANS(8)=( SQRT(P*Q)/Y-R*R)/SQRT(FN)                                BISE1170
C                                                                       BISE1180
C        STORE RESULTS                                                  BISE1190
C                                                                       BISE1200
      ANS(3)=P                                                          BISE1210
      ANS(4)=Q                                                          BISE1220
      ANS(7)=R                                                          BISE1230
C                                                                       BISE1240
   65 RETURN                                                            BISE1250
      END                                                               BISE1260
