C                                                                       DHAR  10
C     ..................................................................DHAR  20
C                                                                       DHAR  30
C        SUBROUTINE DHARM                                               DHAR  40
C                                                                       DHAR  50
C        PURPOSE                                                        DHAR  60
C           PERFORMS DISCRETE COMPLEX FOURIER TRANSFORMS ON A COMPLEX   DHAR  70
C           DOUBLE PRECISION,THREE DIMENSIONAL ARRAY                    DHAR  80
C                                                                       DHAR  90
C        USAGE                                                          DHAR 100
C           CALL DHARM(A,M,INV,S,IFSET,IFERR)                           DHAR 110
C                                                                       DHAR 120
C        DESCRIPTION OF PARAMETERS                                      DHAR 130
C           A     - A DOUBLE PRECISION VECTOR                           DHAR 140
C                   AS INPUT, A CONTAINS THE COMPLEX, 3-DIMENSIONAL     DHAR 150
C                   ARRAY TO BE TRANSFORMED.  THE REAL PART OF          DHAR 160
C                   A(I1,I2,I3) IS STORED IN VECTOR FASHION IN A CELL   DHAR 170
C                   WITH INDEX 2*(I3*N1*N2 + I2*N1 + I1) + 1 WHERE      DHAR 180
C                   NI = 2**M(I), I=1,2,3 AND I1 = 0,1,...,N1-1 ETC.    DHAR 190
C                   THE IMAGINARY PART IS IN THE CELL IMMEDIATELY       DHAR 200
C                   FOLLOWING.  NOTE THAT THE SUBSCRIPT I1 INCREASES    DHAR 210
C                   MOST RAPIDLY AND I3 INCREASES LEAST RAPIDLY.        DHAR 220
C                   AS OUTPUT, A CONTAINS THE COMPLEX FOURIER           DHAR 230
C                   TRANSFORM.  THE NUMBER OF CORE LOCATIONS OF         DHAR 240
C                   ARRAY A IS 2*(N1*N2*N3)                             DHAR 250
C           M     - A THREE CELL VECTOR WHICH DETERMINES THE SIZES      DHAR 260
C                   OF THE 3 DIMENSIONS OF THE ARRAY A.   THE SIZE,     DHAR 270
C                   NI, OF THE I DIMENSION OF A IS 2**M(I), I = 1,2,3   DHAR 280
C           INV   - A VECTOR WORK AREA FOR BIT AND INDEX MANIPULATION   DHAR 290
C                   OF DIMENSION ONE FOURTH OF THE QUANTITY             DHAR 300
C                   MAX(N1,N2,N3)                                       DHAR 310
C           S     - A DOUBLE PRECISION VECTOR WORK AREA FOR SINE TABLES DHAR 320
C                   WITH DIMENSION THE SAME AS INV                      DHAR 330
C           IFSET - AN OPTION PARAMETER WITH THE FOLLOWING SETTINGS     DHAR 340
C                      0    SET UP SINE AND INV TABLES ONLY             DHAR 350
C                      1    SET UP SINE AND INV TABLES ONLY AND         DHAR 360
C                           CALCULATE FOURIER TRANSFORM                 DHAR 370
C                     -1    SET UP SINE AND INV TABLES ONLY AND         DHAR 380
C                           CALCULATE INVERSE FOURIER TRANSFORM (FOR    DHAR 390
C                           THE MEANING OF INVERSE SEE THE EQUATIONS    DHAR 400
C                           UNDER METHOD BELOW)                         DHAR 410
C                      2    CALCULATE FOURIER TRANSFORM ONLY (ASSUME    DHAR 420
C                           SINE AND INV TABLES EXIST)                  DHAR 430
C                     -2    CALCULATE INVERSE FOURIER TRANSFORM ONLY    DHAR 440
C                           (ASSUME SINE AND INV TABLES EXIST)          DHAR 450
C           IFERR - ERROR INDICATOR.   WHEN IFSET IS 0,+1,-1,           DHAR 460
C                   IFERR = 1 MEANS THE MAXIMUM M(I) IS GREATER THAN    DHAR 470
C                   20, I=1,2,3   WHEN IFSET IS 2,-2 , IFERR = 1        DHAR 480
C                   MEANS THAT THE SINE AND INV TABLES ARE NOT LARGE    DHAR 490
C                   ENOUGH OR HAVE NOT BEEN COMPUTED .                  DHAR 500
C                   IF ON RETURN IFERR = 0 THEN NONE OF THE ABOVE       DHAR 510
C                   CONDITIONS ARE PRESENT                              DHAR 520
C                                                                       DHAR 530
C        REMARKS                                                        DHAR 540
C           THIS SUBROUTINE IS TO BE USED FOR COMPLEX, DOUBLE PRECISION,DHAR 550
C           3-DIMENSIONAL ARRAYS IN WHICH EACH DIMENSION IS A POWER OF  DHAR 560
C           2. THE MAXIMUM M(I) MUST NOT BE LESS THAN 3 OR GREATER THAN DHAR 570
C           20, I = 1,2,3.                                              DHAR 580
C                                                                       DHAR 590
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DHAR 600
C           NONE                                                        DHAR 610
C                                                                       DHAR 620
C        METHOD                                                         DHAR 630
C           FOR IFSET = +1, OR +2, THE FOURIER TRANSFORM OF COMPLEX     DHAR 640
C           ARRAY A IS OBTAINED.                                        DHAR 650
C                                                                       DHAR 660
C                  N1-1   N2-1   N3-1                L1   L2   L3       DHAR 670
C     X(J1,J2,J3)=SUM    SUM    SUM    A(K1,K2,K3)*W1  *W2  *W3         DHAR 680
C                  K1=0   K2=0   K3=0                                   DHAR 690
C                                                                       DHAR 700
C                  WHERE WI IS THE N(I) ROOT OF UNITY AND L1=K1*J1,     DHAR 710
C                        L2=K2*J2, L3=K3*J3                             DHAR 720
C                                                                       DHAR 730
C                                                                       DHAR 740
C           FOR IFSET = -1, OR -2, THE INVERSE FOURIER TRANSFORM A OF   DHAR 750
C           COMPLEX ARRAY X IS OBTAINED.                                DHAR 760
C                                                                       DHAR 770
C     A(K1,K2,K3)=                                                      DHAR 780
C               1      N1-1   N2-1   N3-1                -L1  -L2  -L3  DHAR 790
C           -------- *SUM    SUM    SUM    X(J1,J2,J3)*W1  *W2  *W3     DHAR 800
C           N1*N2*N3   J1=0   J2=0   J3=0                               DHAR 810
C                                                                       DHAR 820
C                                                                       DHAR 830
C           SEE J.W. COOLEY AND J.W. TUKEY, 'AN ALGORITHM FOR THE       DHAR 840
C           MACHINE CALCULATION OF COMPLEX FOURIER SERIES',             DHAR 850
C           MATHEMATICS OF COMPUTATIONS, VOL. 19 (APR. 1965), P. 297.   DHAR 860
C                                                                       DHAR 870
C     ..................................................................DHAR 880
C                                                                       DHAR 890
      SUBROUTINE DHARM(A,M,INV,S,IFSET,IFERR)                           DHAR 900
      DIMENSION A(1),INV(1),S(1),N(3),M(3),NP(3),W(2),W2(2),W3(2)       DHAR 910
      DOUBLE PRECISION A,R,W3,AWI,THETA,ROOT2,S,T,W,W2,FN,AWR           DHAR 920
      EQUIVALENCE (N1,N(1)),(N2,N(2)),(N3,N(3))                         DHAR 930
   10 IF( IABS(IFSET) - 1) 900,900,12                                   DHAR 940
   12 MTT=MAX0(M(1),M(2),M(3)) -2                                       DHAR 950
      ROOT2=DSQRT(2.0D0)                                                DHAR 960
      IF (MTT-MT ) 14,14,13                                             DHAR 970
   13 IFERR=1                                                           DHAR 980
      RETURN                                                            DHAR 990
C                   IFERR = 1 MEANS THE MAXIMUM M(I) IS GREATER THAN    DHARM470
C                  20 , I=1,2,3   WHEN IFSET IS 2,-2 , IFERR = 1        DHARM480
C                   MEANS THAT THE SINE AND INV TABLES ARE NOT LARGE    DHARM490
C                   ENOUGH OR HAVE NOT BEEN COMPUTED .                  DHARM500
C                   IF ON RETURN IFERR = 0 THEN NONE OF THE ABOVE       DHARM510
C                   CONDITIONS ARE PRESENT                              DHARM520
C           -------- *SUM    SUM    SUM    X(J1,J2,J3)*W1  *W2  *W3     DHARM800
   14 IFERR=0                                                           DHAR1000
      M1=M(1)                                                           DHAR1010
      M2=M(2)                                                           DHAR1020
      M3=M(3)                                                           DHAR1030
      N1=2**M1                                                          DHAR1040
      N2=2**M2                                                          DHAR1050
      N3=2**M3                                                          DHAR1060
   16 IF(IFSET) 18,18,20                                                DHAR1070
   18 NX= N1*N2*N3                                                      DHAR1080
      FN = NX                                                           DHAR1090
      DO 19 I = 1,NX                                                    DHAR1100
      A(2*I-1) = A(2*I-1)/FN                                            DHAR1110
   19 A(2*I) = -A(2*I)/FN                                               DHAR1120
   20 NP(1)=N1*2                                                        DHAR1130
      NP(2)= NP(1)*N2                                                   DHAR1140
      NP(3)=NP(2)*N3                                                    DHAR1150
      DO 250 ID=1,3                                                     DHAR1160
      IL = NP(3)-NP(ID)                                                 DHAR1170
      IL1 = IL+1                                                        DHAR1180
      MI = M(ID)                                                        DHAR1190
      IF (MI)250,250,30                                                 DHAR1200
   30 IDIF=NP(ID)                                                       DHAR1210
      KBIT=NP(ID)                                                       DHAR1220
      MEV = 2*(MI/2)                                                    DHAR1230
      IF (MI - MEV )60,60,40                                            DHAR1240
C                                                                       DHAR1250
C     M IS ODD. DO L=1 CASE                                             DHAR1260
   40 KBIT=KBIT/2                                                       DHAR1270
      KL=KBIT-2                                                         DHAR1280
      DO 50 I=1,IL1,IDIF                                                DHAR1290
      KLAST=KL+I                                                        DHAR1300
      DO 50 K=I,KLAST,2                                                 DHAR1310
      KD=K+KBIT                                                         DHAR1320
C                                                                       DHAR1330
C     DO ONE STEP WITH L=1,J=0                                          DHAR1340
C     A(K)=A(K)+A(KD)                                                   DHAR1350
C     A(KD)=A(K)-A(KD)                                                  DHAR1360
C                                                                       DHAR1370
      T=A(KD)                                                           DHAR1380
      A(KD)=A(K)-T                                                      DHAR1390
      A(K)=A(K)+T                                                       DHAR1400
      T=A(KD+1)                                                         DHAR1410
      A(KD+1)=A(K+1)-T                                                  DHAR1420
   50 A(K+1)=A(K+1)+T                                                   DHAR1430
      IF (MI - 1)250,250,52                                             DHAR1440
   52 LFIRST =3                                                         DHAR1450
C                                                                       DHAR1460
C     DEF - JLAST = 2**(L-2) -1                                         DHAR1470
      JLAST=1                                                           DHAR1480
      GO TO 70                                                          DHAR1490
C                                                                       DHAR1500
C     M IS EVEN                                                         DHAR1510
   60 LFIRST = 2                                                        DHAR1520
      JLAST=0                                                           DHAR1530
   70 DO 240 L=LFIRST,MI,2                                              DHAR1540
      JJDIF=KBIT                                                        DHAR1550
      KBIT=KBIT/4                                                       DHAR1560
      KL=KBIT-2                                                         DHAR1570
C                                                                       DHAR1580
C     DO FOR J=0                                                        DHAR1590
      DO 80 I=1,IL1,IDIF                                                DHAR1600
      KLAST=I+KL                                                        DHAR1610
      DO 80 K=I,KLAST,2                                                 DHAR1620
      K1=K+KBIT                                                         DHAR1630
      K2=K1+KBIT                                                        DHAR1640
      K3=K2+KBIT                                                        DHAR1650
C                                                                       DHAR1660
C     DO TWO STEPS WITH J=0                                             DHAR1670
C     A(K)=A(K)+A(K2)                                                   DHAR1680
C     A(K2)=A(K)-A(K2)                                                  DHAR1690
C     A(K1)=A(K1)+A(K3)                                                 DHAR1700
C     A(K3)=A(K1)-A(K3)                                                 DHAR1710
C                                                                       DHAR1720
C     A(K)=A(K)+A(K1)                                                   DHAR1730
C     A(K1)=A(K)-A(K1)                                                  DHAR1740
C     A(K2)=A(K2)+A(K3)*I                                               DHAR1750
C     A(K3)=A(K2)-A(K3)*I                                               DHAR1760
C                                                                       DHAR1770
      T=A(K2)                                                           DHAR1780
      A(K2)=A(K)-T                                                      DHAR1790
      A(K)=A(K)+T                                                       DHAR1800
      T=A(K2+1)                                                         DHAR1810
      A(K2+1)=A(K+1)-T                                                  DHAR1820
      A(K+1)=A(K+1)+T                                                   DHAR1830
C                                                                       DHAR1840
      T=A(K3)                                                           DHAR1850
      A(K3)=A(K1)-T                                                     DHAR1860
      A(K1)=A(K1)+T                                                     DHAR1870
      T=A(K3+1)                                                         DHAR1880
      A(K3+1)=A(K1+1)-T                                                 DHAR1890
      A(K1+1)=A(K1+1)+T                                                 DHAR1900
C                                                                       DHAR1910
      T=A(K1)                                                           DHAR1920
      A(K1)=A(K)-T                                                      DHAR1930
      A(K)=A(K)+T                                                       DHAR1940
      T=A(K1+1)                                                         DHAR1950
      A(K1+1)=A(K+1)-T                                                  DHAR1960
      A(K+1)=A(K+1)+T                                                   DHAR1970
C                                                                       DHAR1980
      R=-A(K3+1)                                                        DHAR1990
      T = A(K3)                                                         DHAR2000
      A(K3)=A(K2)-R                                                     DHAR2010
      A(K2)=A(K2)+R                                                     DHAR2020
      A(K3+1)=A(K2+1)-T                                                 DHAR2030
   80 A(K2+1)=A(K2+1)+T                                                 DHAR2040
      IF (JLAST) 235,235,82                                             DHAR2050
   82 JJ=JJDIF   +1                                                     DHAR2060
C                                                                       DHAR2070
C     DO FOR J=1                                                        DHAR2080
      ILAST= IL +JJ                                                     DHAR2090
      DO 85 I = JJ,ILAST,IDIF                                           DHAR2100
      KLAST = KL+I                                                      DHAR2110
      DO 85 K=I,KLAST,2                                                 DHAR2120
      K1 = K+KBIT                                                       DHAR2130
      K2 = K1+KBIT                                                      DHAR2140
      K3 = K2+KBIT                                                      DHAR2150
C                                                                       DHAR2160
C     LETTING W=(1+I)/ROOT2,W3=(-1+I)/ROOT2,W2=I,                       DHAR2170
C     A(K)=A(K)+A(K2)*I                                                 DHAR2180
C     A(K2)=A(K)-A(K2)*I                                                DHAR2190
C     A(K1)=A(K1)*W+A(K3)*W3                                            DHAR2200
C     A(K3)=A(K1)*W-A(K3)*W3                                            DHAR2210
C                                                                       DHAR2220
C     A(K)=A(K)+A(K1)                                                   DHAR2230
C     A(K1)=A(K)-A(K1)                                                  DHAR2240
C     A(K2)=A(K2)+A(K3)*I                                               DHAR2250
C     A(K3)=A(K2)-A(K3)*I                                               DHAR2260
C                                                                       DHAR2270
      R =-A(K2+1)                                                       DHAR2280
      T = A(K2)                                                         DHAR2290
      A(K2) = A(K)-R                                                    DHAR2300
      A(K) = A(K)+R                                                     DHAR2310
      A(K2+1)=A(K+1)-T                                                  DHAR2320
      A(K+1)=A(K+1)+T                                                   DHAR2330
C                                                                       DHAR2340
      AWR=A(K1)-A(K1+1)                                                 DHAR2350
      AWI = A(K1+1)+A(K1)                                               DHAR2360
      R=-A(K3)-A(K3+1)                                                  DHAR2370
      T=A(K3)-A(K3+1)                                                   DHAR2380
      A(K3)=(AWR-R)/ROOT2                                               DHAR2390
      A(K3+1)=(AWI-T)/ROOT2                                             DHAR2400
      A(K1)=(AWR+R)/ROOT2                                               DHAR2410
      A(K1+1)=(AWI+T)/ROOT2                                             DHAR2420
      T= A(K1)                                                          DHAR2430
      A(K1)=A(K)-T                                                      DHAR2440
      A(K)=A(K)+T                                                       DHAR2450
      T=A(K1+1)                                                         DHAR2460
      A(K1+1)=A(K+1)-T                                                  DHAR2470
      A(K+1)=A(K+1)+T                                                   DHAR2480
      R=-A(K3+1)                                                        DHAR2490
      T=A(K3)                                                           DHAR2500
      A(K3)=A(K2)-R                                                     DHAR2510
      A(K2)=A(K2)+R                                                     DHAR2520
      A(K3+1)=A(K2+1)-T                                                 DHAR2530
   85 A(K2+1)=A(K2+1)+T                                                 DHAR2540
      IF(JLAST-1) 235,235,90                                            DHAR2550
   90 JJ= JJ + JJDIF                                                    DHAR2560
C                                                                       DHAR2570
C     NOW DO THE REMAINING J'S                                          DHAR2580
      DO 230 J=2,JLAST                                                  DHAR2590
C                                                                       DHAR2600
C     FETCH W'S                                                         DHAR2610
C     DEF- W=W**INV(J), W2=W**2, W3=W**3                                DHAR2620
   96 I=INV(J+1)                                                        DHAR2630
   98 IC=NT-I                                                           DHAR2640
      W(1)=S(IC)                                                        DHAR2650
      W(2)=S(I)                                                         DHAR2660
      I2=2*I                                                            DHAR2670
      I2C=NT-I2                                                         DHAR2680
      IF(I2C)120,110,100                                                DHAR2690
C                                                                       DHAR2700
C     2*I IS IN FIRST QUADRANT                                          DHAR2710
  100 W2(1)=S(I2C)                                                      DHAR2720
      W2(2)=S(I2)                                                       DHAR2730
      GO TO 130                                                         DHAR2740
  110 W2(1)=0.                                                          DHAR2750
      W2(2)=1.                                                          DHAR2760
      GO TO 130                                                         DHAR2770
C                                                                       DHAR2780
C     2*I IS IN SECOND QUADRANT                                         DHAR2790
  120 I2CC = I2C+NT                                                     DHAR2800
      I2C=-I2C                                                          DHAR2810
      W2(1)=-S(I2C)                                                     DHAR2820
      W2(2)=S(I2CC)                                                     DHAR2830
  130 I3=I+I2                                                           DHAR2840
      I3C=NT-I3                                                         DHAR2850
      IF(I3C)160,150,140                                                DHAR2860
C                                                                       DHAR2870
C     I3 IN FIRST QUADRANT                                              DHAR2880
  140 W3(1)=S(I3C)                                                      DHAR2890
      W3(2)=S(I3)                                                       DHAR2900
      GO TO 200                                                         DHAR2910
  150 W3(1)=0.                                                          DHAR2920
      W3(2)=1.                                                          DHAR2930
      GO TO 200                                                         DHAR2940
C                                                                       DHAR2950
  160 I3CC=I3C+NT                                                       DHAR2960
      IF(I3CC)190,180,170                                               DHAR2970
C                                                                       DHAR2980
C     I3 IN SECOND QUADRANT                                             DHAR2990
  170 I3C=-I3C                                                          DHAR3000
      W3(1)=-S(I3C)                                                     DHAR3010
      W3(2)=S(I3CC)                                                     DHAR3020
      GO TO 200                                                         DHAR3030
  180 W3(1)=-1.                                                         DHAR3040
      W3(2)=0.                                                          DHAR3050
      GO TO 200                                                         DHAR3060
C                                                                       DHAR3070
C     3*I IN THIRD QUADRANT                                             DHAR3080
  190 I3CCC=NT+I3CC                                                     DHAR3090
      I3CC = -I3CC                                                      DHAR3100
      W3(1)=-S(I3CCC)                                                   DHAR3110
      W3(2)=-S(I3CC)                                                    DHAR3120
  200 ILAST=IL+JJ                                                       DHAR3130
      DO 220 I=JJ,ILAST,IDIF                                            DHAR3140
      KLAST=KL+I                                                        DHAR3150
      DO 220 K=I,KLAST,2                                                DHAR3160
      K1=K+KBIT                                                         DHAR3170
      K2=K1+KBIT                                                        DHAR3180
      K3=K2+KBIT                                                        DHAR3190
C                                                                       DHAR3200
C     DO TWO STEPS WITH J NOT 0                                         DHAR3210
C     A(K)=A(K)+A(K2)*W2                                                DHAR3220
C     A(K2)=A(K)-A(K2)*W2                                               DHAR3230
C     A(K1)=A(K1)*W+A(K3)*W3                                            DHAR3240
C     A(K3)=A(K1)*W-A(K3)*W3                                            DHAR3250
C                                                                       DHAR3260
C     A(K)=A(K)+A(K1)                                                   DHAR3270
C     A(K1)=A(K)-A(K1)                                                  DHAR3280
C     A(K2)=A(K2)+A(K3)*I                                               DHAR3290
C     A(K3)=A(K2)-A(K3)*I                                               DHAR3300
C                                                                       DHAR3310
      R=A(K2)*W2(1)-A(K2+1)*W2(2)                                       DHAR3320
      T=A(K2)*W2(2)+A(K2+1)*W2(1)                                       DHAR3330
      A(K2)=A(K)-R                                                      DHAR3340
      A(K)=A(K)+R                                                       DHAR3350
      A(K2+1)=A(K+1)-T                                                  DHAR3360
      A(K+1)=A(K+1)+T                                                   DHAR3370
C                                                                       DHAR3380
      R=A(K3)*W3(1)-A(K3+1)*W3(2)                                       DHAR3390
      T=A(K3)*W3(2)+A(K3+1)*W3(1)                                       DHAR3400
      AWR=A(K1)*W(1)-A(K1+1)*W(2)                                       DHAR3410
      AWI=A(K1)*W(2)+A(K1+1)*W(1)                                       DHAR3420
      A(K3)=AWR-R                                                       DHAR3430
      A(K3+1)=AWI-T                                                     DHAR3440
      A(K1)=AWR+R                                                       DHAR3450
      A(K1+1)=AWI+T                                                     DHAR3460
      T=A(K1)                                                           DHAR3470
      A(K1)=A(K)-T                                                      DHAR3480
      A(K)=A(K)+T                                                       DHAR3490
      T=A(K1+1)                                                         DHAR3500
      A(K1+1)=A(K+1)-T                                                  DHAR3510
      A(K+1)=A(K+1)+T                                                   DHAR3520
      R=-A(K3+1)                                                        DHAR3530
      T=A(K3)                                                           DHAR3540
      A(K3)=A(K2)-R                                                     DHAR3550
      A(K2)=A(K2)+R                                                     DHAR3560
      A(K3+1)=A(K2+1)-T                                                 DHAR3570
  220 A(K2+1)=A(K2+1)+T                                                 DHAR3580
C     END OF I AND K LOOPS                                              DHAR3590
C                                                                       DHAR3600
  230 JJ=JJDIF+JJ                                                       DHAR3610
C     END OF J-LOOP                                                     DHAR3620
C                                                                       DHAR3630
  235 JLAST=4*JLAST+3                                                   DHAR3640
  240 CONTINUE                                                          DHAR3650
C     END OF  L  LOOP                                                   DHAR3660
C                                                                       DHAR3670
  250 CONTINUE                                                          DHAR3680
C     END OF  ID  LOOP                                                  DHAR3690
C                                                                       DHAR3700
C     WE NOW HAVE THE COMPLEX FOURIER SUMS BUT THEIR ADDRESSES ARE      DHAR3710
C     BIT-REVERSED.  THE FOLLOWING ROUTINE PUTS THEM IN ORDER           DHAR3720
      NTSQ=NT*NT                                                        DHAR3730
      M3MT=M3-MT                                                        DHAR3740
  350 IF(M3MT) 370,360,360                                              DHAR3750
C                                                                       DHAR3760
C     M3 GR. OR EQ. MT                                                  DHAR3770
  360 IGO3=1                                                            DHAR3780
      N3VNT=N3/NT                                                       DHAR3790
      MINN3=NT                                                          DHAR3800
      GO TO 380                                                         DHAR3810
C                                                                       DHAR3820
C     M3 LESS THAN MT                                                   DHAR3830
  370 IGO3=2                                                            DHAR3840
      N3VNT=1                                                           DHAR3850
      NTVN3=NT/N3                                                       DHAR3860
      MINN3=N3                                                          DHAR3870
  380 JJD3 = NTSQ/N3                                                    DHAR3880
      M2MT=M2-MT                                                        DHAR3890
  450 IF (M2MT)470,460,460                                              DHAR3900
C                                                                       DHAR3910
C     M2 GR. OR EQ. MT                                                  DHAR3920
  460 IGO2=1                                                            DHAR3930
      N2VNT=N2/NT                                                       DHAR3940
      MINN2=NT                                                          DHAR3950
      GO TO 480                                                         DHAR3960
C                                                                       DHAR3970
C     M2 LESS THAN MT                                                   DHAR3980
  470 IGO2 = 2                                                          DHAR3990
      N2VNT=1                                                           DHAR4000
      NTVN2=NT/N2                                                       DHAR4010
      MINN2=N2                                                          DHAR4020
  480 JJD2=NTSQ/N2                                                      DHAR4030
      M1MT=M1-MT                                                        DHAR4040
  550 IF(M1MT)570,560,560                                               DHAR4050
C                                                                       DHAR4060
C     M1 GR. OR EQ. MT                                                  DHAR4070
  560 IGO1=1                                                            DHAR4080
      N1VNT=N1/NT                                                       DHAR4090
      MINN1=NT                                                          DHAR4100
      GO TO 580                                                         DHAR4110
C                                                                       DHAR4120
C     M1 LESS THAN MT                                                   DHAR4130
  570 IGO1=2                                                            DHAR4140
      N1VNT=1                                                           DHAR4150
      NTVN1=NT/N1                                                       DHAR4160
      MINN1=N1                                                          DHAR4170
  580 JJD1=NTSQ/N1                                                      DHAR4180
  600 JJ3=1                                                             DHAR4190
      J=1                                                               DHAR4200
      DO 880 JPP3=1,N3VNT                                               DHAR4210
      IPP3=INV(JJ3)                                                     DHAR4220
      DO 870 JP3=1,MINN3                                                DHAR4230
      GO TO (610,620),IGO3                                              DHAR4240
  610 IP3=INV(JP3)*N3VNT                                                DHAR4250
      GO TO 630                                                         DHAR4260
  620 IP3=INV(JP3)/NTVN3                                                DHAR4270
  630 I3=(IPP3+IP3)*N2                                                  DHAR4280
  700 JJ2=1                                                             DHAR4290
      DO 870 JPP2=1,N2VNT                                               DHAR4300
      IPP2=INV(JJ2)+I3                                                  DHAR4310
      DO 860 JP2=1,MINN2                                                DHAR4320
      GO TO (710,720),IGO2                                              DHAR4330
  710 IP2=INV(JP2)*N2VNT                                                DHAR4340
      GO TO 730                                                         DHAR4350
  720 IP2=INV(JP2)/NTVN2                                                DHAR4360
  730 I2=(IPP2+IP2)*N1                                                  DHAR4370
  800 JJ1=1                                                             DHAR4380
      DO 860 JPP1=1,N1VNT                                               DHAR4390
      IPP1=INV(JJ1)+I2                                                  DHAR4400
      DO 850 JP1=1,MINN1                                                DHAR4410
      GO TO (810,820),IGO1                                              DHAR4420
  810 IP1=INV(JP1)*N1VNT                                                DHAR4430
      GO TO 830                                                         DHAR4440
  820 IP1=INV(JP1)/NTVN1                                                DHAR4450
  830 I=2*(IPP1+IP1)+1                                                  DHAR4460
      IF (J-I) 840,850,850                                              DHAR4470
  840 T=A(I)                                                            DHAR4480
      A(I)=A(J)                                                         DHAR4490
      A(J)=T                                                            DHAR4500
      T=A(I+1)                                                          DHAR4510
      A(I+1)=A(J+1)                                                     DHAR4520
      A(J+1)=T                                                          DHAR4530
  850 J=J+2                                                             DHAR4550
  860 JJ1=JJ1+JJD1                                                      DHAR4560
C     END OF JPP1 AND JP2                                               DHAR4570
C                                                                       DHAR4580
  870 JJ2=JJ2+JJD2                                                      DHAR4590
C     END OF JPP2 AND JP3 LOOPS                                         DHAR4600
C                                                                       DHAR4610
  880 JJ3 = JJ3+JJD3                                                    DHAR4620
C     END OF JPP3 LOOP                                                  DHAR4630
C                                                                       DHAR4640
  890 IF(IFSET)891,895,895                                              DHAR4650
  891 DO 892 I = 1,NX                                                   DHAR4660
  892 A(2*I) = -A(2*I)                                                  DHAR4670
  895 RETURN                                                            DHAR4680
C                                                                       DHAR4690
C     THE FOLLOWING PROGRAM COMPUTES THE SIN AND INV TABLES.            DHAR4700
C                                                                       DHAR4710
  900 MT=MAX0(M(1),M(2),M(3)) -2                                        DHAR4720
      MT = MAX0(2,MT)                                                   DHAR4730
  904 IF (MT-18) 906,906,13                                             DHAR4740
  906 IFERR=0                                                           DHAR4770
      NT=2**MT                                                          DHAR4780
      NTV2=NT/2                                                         DHAR4790
C                                                                       DHAR4800
C     SET UP SIN TABLE                                                  DHAR4810
C     THETA=PIE/2**(L+1) FOR L=1                                        DHAR4820
  910 THETA=.7853981633974483                                           DHAR4830
C                                                                       DHAR4840
C     JSTEP=2**(MT-L+1) FOR L=1                                         DHAR4850
      JSTEP=NT                                                          DHAR4860
C                                                                       DHAR4870
C     JDIF=2**(MT-L) FOR L=1                                            DHAR4880
      JDIF=NTV2                                                         DHAR4890
      S(JDIF)=DSIN(THETA)                                               DHAR4900
      DO 950 L=2,MT                                                     DHAR4910
      THETA=THETA/2.0D0                                                 DHAR4920
      JSTEP2=JSTEP                                                      DHAR4930
      JSTEP=JDIF                                                        DHAR4940
      JDIF=JSTEP/2                                                      DHAR4950
      S(JDIF)=DSIN(THETA)                                               DHAR4960
      JC1=NT-JDIF                                                       DHAR4970
      S(JC1)=DCOS(THETA)                                                DHAR4980
      JLAST=NT-JSTEP2                                                   DHAR4990
      IF(JLAST - JSTEP) 950,920,920                                     DHAR5000
  920 DO 940 J=JSTEP,JLAST,JSTEP                                        DHAR5010
      JC=NT-J                                                           DHAR5020
      JD=J+JDIF                                                         DHAR5030
  940 S(JD)=S(J)*S(JC1)+S(JDIF)*S(JC)                                   DHAR5040
  950 CONTINUE                                                          DHAR5050
C                                                                       DHAR5060
C     SET UP INV(J) TABLE                                               DHAR5070
C                                                                       DHAR5080
  960 MTLEXP=NTV2                                                       DHAR5090
C                                                                       DHAR5100
C     MTLEXP=2**(MT-L). FOR L=1                                         DHAR5110
      LM1EXP=1                                                          DHAR5120
C                                                                       DHAR5130
C     LM1EXP=2**(L-1). FOR L=1                                          DHAR5140
      INV(1)=0                                                          DHAR5150
      DO 980 L=1,MT                                                     DHAR5160
      INV(LM1EXP+1) = MTLEXP                                            DHAR5170
      DO 970 J=2,LM1EXP                                                 DHAR5180
      JJ=J+LM1EXP                                                       DHAR5190
  970 INV(JJ)=INV(J)+MTLEXP                                             DHAR5200
      MTLEXP=MTLEXP/2                                                   DHAR5210
  980 LM1EXP=LM1EXP*2                                                   DHAR5220
  982 IF(IFSET)12,895,12                                                DHAR5230
      END                                                               DHAR5240
