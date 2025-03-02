C                                                                       BESK  10
C     ..................................................................BESK  20
C                                                                       BESK  30
C        SUBROUTINE BESK                                                BESK  40
C                                                                       BESK  50
C           COMPUTE THE K BESSEL FUNCTION FOR A GIVEN ARGUMENT AND ORDERBESK  60
C                                                                       BESK  70
C        USAGE                                                          BESK  80
C           CALL BESK(X,N,BK,IER)                                       BESK  90
C                                                                       BESK 100
C        DESCRIPTION OF PARAMETERS                                      BESK 110
C           X  -THE ARGUMENT OF THE K BESSEL FUNCTION DESIRED           BESK 120
C           N  -THE ORDER OF THE K BESSEL FUNCTION DESIRED              BESK 130
C           BK -THE RESULTANT K BESSEL FUNCTION                         BESK 140
C           IER-RESULTANT ERROR CODE WHERE                              BESK 150
C              IER=0  NO ERROR                                          BESK 160
C              IER=1  N IS NEGATIVE                                     BESK 170
C              IER=2  X IS ZERO OR NEGATIVE                             BESK 180
C              IER=3  X .GT. 170, MACHINE RANGE EXCEEDED                BESK 190
C              IER=4  BK .GT. 10**70                                    BESK 200
C                                                                       BESK 210
C        REMARKS                                                        BESK 220
C           N MUST BE GREATER THAN OR EQUAL TO ZERO                     BESK 230
C                                                                       BESK 240
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  BESK 250
C           NONE                                                        BESK 260
C                                                                       BESK 270
C        METHOD                                                         BESK 280
C           COMPUTES ZERO ORDER AND FIRST ORDER BESSEL FUNCTIONS USING  BESK 290
C           SERIES APPROXIMATIONS AND THEN COMPUTES N TH ORDER FUNCTION BESK 300
C           USING RECURRENCE RELATION.                                  BESK 310
C           RECURRENCE RELATION AND POLYNOMIAL APPROXIMATION TECHNIQUE  BESK 320
C           AS DESCRIBED BY A.J.M.HITCHCOCK,'POLYNOMIAL APPROXIMATIONS  BESK 330
C           TO BESSEL FUNCTIONS OF ORDER ZERO AND ONE AND TO RELATED    BESK 340
C           FUNCTIONS', M.T.A.C., V.11,1957,PP.86-88, AND G.N. WATSON,  BESK 350
C           'A TREATISE ON THE THEORY OF BESSEL FUNCTIONS', CAMBRIDGE   BESK 360
C           UNIVERSITY PRESS, 1958, P. 62                               BESK 370
C                                                                       BESK 380
C     ..................................................................BESK 390
C                                                                       BESK 400
      SUBROUTINE BESK(X,N,BK,IER)                                       BESK 410
      DIMENSION T(12)                                                   BESK 420
      BK=.0                                                             BESK 430
      IF(N)10,11,11                                                     BESK 440
   10 IER=1                                                             BESK 450
      RETURN                                                            BESK 460
   11 IF(X)12,12,20                                                     BESK 470
   12 IER=2                                                             BESK 480
      RETURN                                                            BESK 490
   20 IF(X-170.0)22,22,21                                               BESK 500
   21 IER=3                                                             BESK 510
      RETURN                                                            BESK 520
   22 IER=0                                                             BESK 530
      IF(X-1.)36,36,25                                                  BESK 540
   25 A=EXP(-X)                                                         BESK 550
      B=1./X                                                            BESK 560
      C=SQRT(B)                                                         BESK 570
      T(1)=B                                                            BESK 580
      DO 26 L=2,12                                                      BESK 590
   26 T(L)=T(L-1)*B                                                     BESK 600
      IF(N-1)27,29,27                                                   BESK 610
C                                                                       BESK 620
C     COMPUTE KO USING POLYNOMIAL APPROXIMATION                         BESK 630
C                                                                       BESK 640
   27 G0=A*(1.2533141-.1566642*T(1)+.08811128*T(2)-.09139095*T(3)       BESK 650
     2+.1344596*T(4)-.2299850*T(5)+.3792410*T(6)-.5247277*T(7)          BESK 660
     3+.5575368*T(8)-.4262633*T(9)+.2184518*T(10)-.06680977*T(11)       BESK 670
     4+.009189383*T(12))*C                                              BESK 680
      IF(N)20,28,29                                                     BESK 690
   28 BK=G0                                                             BESK 700
      RETURN                                                            BESK 710
C                                                                       BESK 720
C     COMPUTE K1 USING POLYNOMIAL APPROXIMATION                         BESK 730
C                                                                       BESK 740
   29 G1=A*(1.2533141+.4699927*T(1)-.1468583*T(2)+.1280427*T(3)         BESK 750
     2-.1736432*T(4)+.2847618*T(5)-.4594342*T(6)+.6283381*T(7)          BESK 760
     3-.6632295*T(8)+.5050239*T(9)-.2581304*T(10)+.07880001*T(11)       BESK 770
     4-.01082418*T(12))*C                                               BESK 780
      IF(N-1)20,30,31                                                   BESK 790
   30 BK=G1                                                             BESK 800
      RETURN                                                            BESK 810
C                                                                       BESK 820
C     FROM KO,K1 COMPUTE KN USING RECURRENCE RELATION                   BESK 830
C                                                                       BESK 840
   31 DO 35 J=2,N                                                       BESK 850
      GJ=2.*(FLOAT(J)-1.)*G1/X+G0                                       BESK 860
      IF(GJ-1.0E70)33,33,32                                             BESK 870
   32 IER=4                                                             BESK 880
      GO TO 34                                                          BESK 890
   33 G0=G1                                                             BESK 900
   35 G1=GJ                                                             BESK 910
   34 BK=GJ                                                             BESK 920
      RETURN                                                            BESK 930
   36 B=X/2.                                                            BESK 940
      A=.5772157+ALOG(B)                                                BESK 950
      C=B*B                                                             BESK 960
      IF(N-1)37,43,37                                                   BESK 970
C                                                                       BESK 980
C     COMPUTE KO USING SERIES EXPANSION                                 BESK 990
C                                                                       BESK1000
   37 G0=-A                                                             BESK1010
      X2J=1.                                                            BESK1020
      FACT=1.                                                           BESK1030
      HJ=.0                                                             BESK1040
      DO 40 J=1,6                                                       BESK1050
      RJ=1./FLOAT(J)                                                    BESK1060
      X2J=X2J*C                                                         BESK1070
      FACT=FACT*RJ*RJ                                                   BESK1080
      HJ=HJ+RJ                                                          BESK1090
   40 G0=G0+X2J*FACT*(HJ-A)                                             BESK1100
      IF(N)43,42,43                                                     BESK1110
   42 BK=G0                                                             BESK1120
      RETURN                                                            BESK1130
C                                                                       BESK1140
C     COMPUTE K1 USING SERIES EXPANSION                                 BESK1150
C                                                                       BESK1160
   43 X2J=B                                                             BESK1170
      FACT=1.                                                           BESK1180
      HJ=1.                                                             BESK1190
      G1=1./X+X2J*(.5+A-HJ)                                             BESK1200
      DO 50 J=2,8                                                       BESK1210
      X2J=X2J*C                                                         BESK1220
      RJ=1./FLOAT(J)                                                    BESK1230
      FACT=FACT*RJ*RJ                                                   BESK1240
      HJ=HJ+RJ                                                          BESK1250
   50 G1=G1+X2J*FACT*(.5+(A-HJ)*FLOAT(J))                               BESK1260
      IF(N-1)31,52,31                                                   BESK1270
   52 BK=G1                                                             BESK1280
      RETURN                                                            BESK1290
      END                                                               BESK1300
