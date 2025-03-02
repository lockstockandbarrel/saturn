C                                                                       DEL2  10
C     ..................................................................DEL2  20
C                                                                       DEL2  30
C        SUBROUTINE DELI2                                               DEL2  40
C                                                                       DEL2  50
C        PURPOSE                                                        DEL2  60
C           COMPUTES THE GENERALIZED ELLIPTIC INTEGRAL OF SECOND KIND   DEL2  70
C                                                                       DEL2  80
C        USAGE                                                          DEL2  90
C           CALL DELI2(R,X,CK,A,B)                                      DEL2 100
C                                                                       DEL2 110
C        DESCRIPTION OF PARAMETERS                                      DEL2 120
C           R     - RESULT VALUE IN DOUBLE PRECISION                    DEL2 130
C           X     - UPPER INTEGRATION BOUND (ARGUMENT OF ELLIPTIC       DEL2 140
C                   INTEGRAL OF SECOND KIND) IN DOUBLE PRECISION        DEL2 150
C           CK    - COMPLEMENTARY MODULUS IN DOUBLE PRECISION           DEL2 160
C           A     - DOUBLE PRECISION CONSTANT TERM IN NUMERATOR         DEL2 170
C           B     - DOUBLE PRECISION QUATRATIC TERM IN NUMERATOR        DEL2 180
C                                                                       DEL2 190
C        REMARKS                                                        DEL2 200
C           DOUBLE PRECISION MODULUS K = DSQRT(1.D0-CK*CK).             DEL2 210
C           SPECIAL CASES OF THE GENERALIZED ELLIPTIC INTEGRAL OF       DEL2 220
C           SECOND KIND ARE                                             DEL2 230
C           F(DATAN(X),K) OBTAINED WITH A=1.D0, B=1.D0                  DEL2 240
C           E(DATAN(X),K) OBTAINED WITH A=1.D0, B=CK*CK                 DEL2 250
C           B(DATAN(X),K) OBTAINED WITH A=1.D0, B=0.D0                  DEL2 260
C           D(DATAN(X),K) OBTAINED WITH A=0.D0, B=1.D0.                 DEL2 270
C                                                                       DEL2 280
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DEL2 290
C           NONE                                                        DEL2 300
C                                                                       DEL2 310
C        METHOD                                                         DEL2 320
C           DEFINITION                                                  DEL2 330
C           R=INTEGRAL((A+B*T*T)/(SQRT((1+T*T)*(1+(CK*T)**2))*(1+T*T)), DEL2 340
C                  SUMMED OVER T FROM 0 TO X).                          DEL2 350
C           EQUIVALENT IS THE DEFINITION                                DEL2 360
C           R=INTEGRAL((A+(B-A)*(SIN(T))**2)/SQRT(1-(K*SIN(T))**2),     DEL2 370
C                  SUMMED OVER T FROM 0 TO ATAN(X)).                    DEL2 380
C           EVALUATION                                                  DEL2 390
C           LANDENS TRANSFORMATION IS USED FOR CALCULATION.             DEL2 400
C           REFERENCE                                                   DEL2 410
C           R. BULIRSCH, NUMERICAL CALCULATION OF ELLIPTIC INTEGRALS ANDDEL2 420
C                  ELLIPTIC FUNCTIONS                                   DEL2 430
C                  HANDBOOK SERIES OF SPECIAL FUNCTIONS                 DEL2 440
C                  NUMERISCHE MATHEMATIK VOL. 7, 1965, PP. 78-90.       DEL2 450
C                                                                       DEL2 460
C     ..................................................................DEL2 470
C                                                                       DEL2 480
      SUBROUTINE DELI2(R,X,CK,A,B)                                      DEL2 490
C                                                                       DEL2 500
      DOUBLE PRECISION R,X,A,B,AN,AA,ANG,AANG,PIM,PIMA,ARI,AARI         DEL2 510
      DOUBLE PRECISION GEO,SGEO,C,D,P,CK                                DEL2 520
C                                                                       DEL2 530
C        TEST ARGUMENT                                                  DEL2 540
C                                                                       DEL2 550
      IF(X)2,1,2                                                        DEL2 560
    1 R=0.D0                                                            DEL2 570
      RETURN                                                            DEL2 580
C                                                                       DEL2 590
C        TEST MODULUS                                                   DEL2 600
C                                                                       DEL2 610
    2 C=0.D0                                                            DEL2 620
      D=0.5D0                                                           DEL2 630
      IF(CK)7,3,7                                                       DEL2 640
    3 R=DSQRT(1.D0+X*X)                                                 DEL2 650
      R=(A-B)*DABS(X)/R+B*DLOG(DABS(X)+R)                               DEL2 660
    4 R=R+C*(A-B)                                                       DEL2 670
C                                                                       DEL2 680
C        TEST SIGN OF ARGUMENT                                          DEL2 690
C                                                                       DEL2 700
      IF(X)5,6,6                                                        DEL2 710
    5 R=-R                                                              DEL2 720
    6 RETURN                                                            DEL2 730
C                                                                       DEL2 740
C        INITIALIZATION                                                 DEL2 750
C                                                                       DEL2 760
    7 AN=(B+A)*0.5D0                                                    DEL2 770
      AA=A                                                              DEL2 780
      R=B                                                               DEL2 790
      ANG=DABS(1.D0/X)                                                  DEL2 800
      PIM=0.D0                                                          DEL2 810
      ISI=0                                                             DEL2 820
      ARI=1.D0                                                          DEL2 830
      GEO=DABS(CK)                                                      DEL2 840
C                                                                       DEL2 850
C        LANDEN TRANSFORMATION                                          DEL2 860
C                                                                       DEL2 870
    8 R=AA*GEO+R                                                        DEL2 880
      SGEO=ARI*GEO                                                      DEL2 890
      AA=AN                                                             DEL2 900
      AARI=ARI                                                          DEL2 910
C                                                                       DEL2 920
C        ARITHMETIC MEAN                                                DEL2 930
C                                                                       DEL2 940
      ARI=GEO+ARI                                                       DEL2 950
C                                                                       DEL2 960
C        SUM OF SINE VALUES                                             DEL2 970
C                                                                       DEL2 980
      AN=(R/ARI+AA)*0.5D0                                               DEL2 990
      AANG=DABS(ANG)                                                    DEL21000
      ANG=-SGEO/ANG+ANG                                                 DEL21010
      PIMA=PIM                                                          DEL21020
      IF(ANG)10,9,11                                                    DEL21030
C                                                                       DEL21040
C        REPLACE 0 BY SMALL VALUE                                       DEL21050
C                                                                       DEL21060
    9 ANG=-1.D-17*AANG                                                  DEL21070
   10 PIM=PIM+3.1415926535897932                                        DEL21080
      ISI=ISI+1                                                         DEL21090
   11 AANG=ARI*ARI+ANG*ANG                                              DEL21100
      P=D/DSQRT(AANG)                                                   DEL21110
      IF(ISI-4)13,12,12                                                 DEL21120
   12 ISI=ISI-4                                                         DEL21130
   13 IF(ISI-2)15,14,14                                                 DEL21140
   14 P=-P                                                              DEL21150
   15 C=C+P                                                             DEL21160
      D=D*(AARI-GEO)*0.5D0/ARI                                          DEL21170
      IF(DABS(AARI-GEO)-1.D-9*AARI)17,17,16                             DEL21180
   16 SGEO=DSQRT(SGEO)                                                  DEL21190
C                                                                       DEL21200
C        GEOMETRIC MEAN                                                 DEL21210
C                                                                       DEL21220
      GEO=SGEO+SGEO                                                     DEL21230
      PIM=PIM+PIMA                                                      DEL21240
      ISI=ISI+ISI                                                       DEL21250
      GOTO 8                                                            DEL21260
C                                                                       DEL21270
C        ACCURACY WAS SUFFICIENT                                        DEL21280
C                                                                       DEL21290
   17 R=(DATAN(ARI/ANG)+PIM)*AN/ARI                                     DEL21300
      C=C+D*ANG/AANG                                                    DEL21310
      GOTO 4                                                            DEL21320
      END                                                               DEL21330
