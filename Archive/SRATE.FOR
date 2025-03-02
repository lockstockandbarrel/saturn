C                                                                       SRAT  10
C     ..................................................................SRAT  20
C                                                                       SRAT  30
C        SUBROUTINE SRATE                                               SRAT  40
C                                                                       SRAT  50
C        PURPOSE                                                        SRAT  60
C           TO COMPUTE THE PROPORTION OF SUBJECTS SURVIVING, THE        SRAT  70
C           SURVIVAL RATES AND THE STANDARD ERRORS FOR SUCCESSIVELY     SRAT  80
C           REDUCED TIME PERIODS.  THE SURVIVAL RATE IS COMPUTED FOR    SRAT  90
C           EACH OF K PERIODS, WHERE K IS A CONSTANT TO BE SPECIFIED ANDSRAT 100
C           IS LESS THAN OR EQUAL TO N (WHERE N = TOTAL NUMBER OF       SRAT 110
C           PERIODS).                                                   SRAT 120
C                                                                       SRAT 130
C        USAGE                                                          SRAT 140
C           CALL SRATE (N,K,X,IE)                                       SRAT 150
C                                                                       SRAT 160
C        DESCRIPTION OF PARAMETERS                                      SRAT 170
C           N - THE TOTAL NUMBER OF PERIODS AFTER TREATMENT OR DIAGNOSISSRAT 180
C           K - THE SPECIFIED PERIOD UP TO WHICH SURVIVAL RATES ARE TO  SRAT 190
C               BE CALCULATED                                           SRAT 200
C           X - AN INPUT AND OUTPUT MATRIX (N X 9) CONTAINING THE       SRAT 210
C               FOLLOWING INFORMATION                                   SRAT 220
C               FOR INPUT--STORED IN THE N ROWS OF EACH COLUMN          SRAT 230
C                  COL 1 - NUMBER OF SUBJECTS ALIVE AT THE BEGINNING OF SRAT 240
C                          PERIOD                                       SRAT 250
C                  COL 2 - NUMBER OF SUBJECTS WHICH DIED DURING THE     SRAT 260
C                          PERIOD                                       SRAT 270
C                  COL 3 - NUMBER OF SUBJECTS LOST TO FOLLOW-UP DURING  SRAT 280
C                          THE PERIOD                                   SRAT 290
C                  COL 4 - NUMBER OF SUBJECTS WITHDRAWN ALIVE DURING THESRAT 300
C                          PERIOD                                       SRAT 310
C                    FOR OUTPUT--STORED IN THE FIRST K ROWS OF EACH     SRAT 320
C                    COLUMN                                             SRAT 330
C                  COL 5 - EFFECTIVE NUMBER EXPOSED TO THE RISK OF DYINGSRAT 340
C                  COL 6 - PROPORTION WHO DIED DURING THE PERIOD        SRAT 350
C                  COL 7 - PROPORTION WHO SURVIVED DURING THE PERIOD    SRAT 360
C                  COL 8 - SURVIAL RATE                                 SRAT 370
C                  COL 9 - STANDARD ERROR OF THE SURVIVAL RATE          SRAT 380
C           IE- 1, IF K IS NOT IN THE CLOSED INTERVAL (0,N).            SRAT 390
C               2, IF THE NUMBER OF SUBJECTS ALIVE AT THE BEGINNING     SRAT 400
C               OF PERIOD I IS LESS THAN THE SUM OF THOSE WHICH DIED,   SRAT 410
C               WERE LOST, OR WERE WITHDRAWN DURING PERIOD I (I=1,...N) SRAT 420
C               3, IF THE NUMBER OF SUBJECTS WHICH DIED, WERE LOST, OR  SRAT 430
C               WERE WITHDRAWN IN PERIOD I IS NOT EQUAL TO THE NUMBER   SRAT 440
C               ALIVE AT THE BEGINNING OF PERIOD I LESS THE NUMBER      SRAT 450
C               ALIVE AT THE BEGINNING OF PERIOD I + 1 (I=1,...N-1)     SRAT 460
C                                                                       SRAT 470
C        REMARKS                                                        SRAT 480
C           IF THE SUBJECTS IN A GIVEN GROUP ARE ALL DIAGNOSED OR       SRAT 490
C           TREATED AT THE SAME TIME, THE CONSTANT K MAY BE SET EQUAL TOSRAT 500
C           N .  IF THE SUBJECTS IN A GIVEN GOUP ENTER THE STUDY AT     SRAT 510
C           VARYING TIMES, K CAN BE NO GREATER THAN N-1.                SRAT 520
C                                                                       SRAT 530
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  SRAT 540
C           NONE                                                        SRAT 550
C                                                                       SRAT 560
C        METHOD                                                         SRAT 570
C           REFER TO S. J. CUTLER AND F. EDERER 'MAXIMUM UTILIZATION OF SRAT 580
C           THE LIFE TABLE METHOD IN ANALYZING SURVIVAL', JOURNAL OF    SRAT 590
C           CHRONIC DISEASES, DECEMBER, 1958.  PP 699-712.              SRAT 600
C                                                                       SRAT 610
C     ..................................................................SRAT 620
C                                                                       SRAT 630
      SUBROUTINE SRATE (N,K,X,IE)                                       SRAT 640
C                                                                       SRAT 650
      DIMENSION X(1)                                                    SRAT 660
C                                                                       SRAT 670
C        INITIALIZATION AND ERROR CHECKING                              SRAT 680
C                                                                       SRAT 690
      IE=0                                                              SRAT 700
      NP4=4*N+1                                                         SRAT 710
      NP9=NP4+NP4+N-2                                                   SRAT 720
      DO 1 I=NP4,NP9                                                    SRAT 730
    1 X(I)=0.0                                                          SRAT 740
      IF (K) 2,2,3                                                      SRAT 750
    2 IE=1                                                              SRAT 760
      GO TO 45                                                          SRAT 770
    3 IF(K-N) 4,4,2                                                     SRAT 780
    4 DO 9 I=1,N                                                        SRAT 790
      NP4=I+N                                                           SRAT 800
      NP9=NP4+N                                                         SRAT 810
      NP1=NP9+N                                                         SRAT 820
      IF(INT(X(I)-X(NP4)-X(NP9)-X(NP1)+.01)) 5,6,6                      SRAT 830
    5 IE=2                                                              SRAT 840
      GO TO 45                                                          SRAT 850
    6 IF(I-N) 7,9,9                                                     SRAT 860
    7 IF (INT(X(I+1)-X(I)+X(NP4)+X(NP9)+X(NP1)+.01)) 8,9,8              SRAT 870
    8 IE=3                                                              SRAT 880
        GO TO 45                                                        SRAT 890
    9 CONTINUE                                                          SRAT 900
   15 L1=0                                                              SRAT 910
      L2=L1+N                                                           SRAT 920
      L3=L2+N                                                           SRAT 930
      L4=L3+N                                                           SRAT 940
      L5=L4+N                                                           SRAT 950
      L6=L5+N                                                           SRAT 960
      L7=L6+N                                                           SRAT 970
      L8=L7+N                                                           SRAT 980
      L9=L8+N                                                           SRAT 990
      LD=L2                                                             SRAT1000
      LE=L5                                                             SRAT1010
      LQ=L6                                                             SRAT1020
      SUM=0.0                                                           SRAT1030
C                                                                       SRAT1040
      DO 40 I=1,K                                                       SRAT1050
C                                                                       SRAT1060
C        COMPUTE EFFECTIVE NUMBER EXPOSED TO RISK OF DYING              SRAT1070
C                                                                       SRAT1080
      L1=L1+1                                                           SRAT1090
      L3=L3+1                                                           SRAT1100
      L4=L4+1                                                           SRAT1110
      L5=L5+1                                                           SRAT1120
      X(L5)=X(L1)-(X(L3)+X(L4))/2.0                                     SRAT1130
C                                                                       SRAT1140
C        COMPUTE PROPORTION OF DYING                                    SRAT1150
C                                                                       SRAT1160
      L2=L2+1                                                           SRAT1170
      L6=L6+1                                                           SRAT1180
      X(L6)=X(L2)/X(L5)                                                 SRAT1190
C                                                                       SRAT1200
C        COMPUTE PROPORTION OF SURVIVING                                SRAT1210
C                                                                       SRAT1220
      L7=L7+1                                                           SRAT1230
      X(L7)=1.0-X(L6)                                                   SRAT1240
C                                                                       SRAT1250
C        COMPUTE SURVIVAL RATE                                          SRAT1260
C                                                                       SRAT1270
      L8=L8+1                                                           SRAT1280
      IF (I-1) 20, 20, 25                                               SRAT1290
   20 X(L8)=X(L7)                                                       SRAT1300
      GO TO 30                                                          SRAT1310
   25 X(L8)=X(L8-1)*X(L7)                                               SRAT1320
C                                                                       SRAT1330
C        COMPUTE STANDARD ERROR OF SURVIVAL RATE                        SRAT1340
C                                                                       SRAT1350
   30 L9=L9+1                                                           SRAT1360
      SUM=SUM+X(L6)/(X(L5)-X(L2))                                       SRAT1370
   40 X(L9)=X(L8)*SQRT(SUM)                                             SRAT1380
C                                                                       SRAT1390
   45 RETURN                                                            SRAT1400
      END                                                               SRAT1410
