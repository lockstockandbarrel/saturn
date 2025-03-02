C                                                                       HARM  10
C     ..................................................................HARM  20
C                                                                       HARM  30
C        SUBROUTINE HARM                                                HARM  40
C                                                                       HARM  50
C        PURPOSE                                                        HARM  60
C           PERFORMS DISCRETE COMPLEX FOURIER TRANSFORMS ON A COMPLEX   HARM  70
C           THREE DIMENSIONAL ARRAY                                     HARM  80
C                                                                       HARM  90
C        USAGE                                                          HARM 100
C           CALL HARM (A,M,INV,S,IFSET,IFERR)                           HARM 110
C                                                                       HARM 120
C        DESCRIPTION OF PARAMETERS                                      HARM 130
C           A     - AS INPUT, A CONTAINS THE COMPLEX, 3-DIMENSIONAL     HARM 140
C                   ARRAY TO BE TRANSFORMED.  THE REAL PART OF          HARM 150
C                   A(I1,I2,I3) IS STORED IN VECTOR FASHION IN A CELL   HARM 160
C                   WITH INDEX 2*(I3*N1*N2 + I2*N1 + I1) + 1 WHERE      HARM 170
C                   NI = 2**M(I), I=1,2,3 AND I1 = 0,1,...,N1-1 ETC.    HARM 180
C                   THE IMAGINARY PART IS IN THE CELL IMMEDIATELY       HARM 190
C                   FOLLOWING.  NOTE THAT THE SUBSCRIPT I1 INCREASES    HARM 200
C                   MOST RAPIDLY AND I3 INCREASES LEAST RAPIDLY.        HARM 210
C                   AS OUTPUT, A CONTAINS THE COMPLEX FOURIER           HARM 220
C                   TRANSFORM.  THE NUMBER OF CORE LOCATIONS OF         HARM 230
C                   ARRAY A IS 2*(N1*N2*N3)                             HARM 240
C           M     - A THREE CELL VECTOR WHICH DETERMINES THE SIZES      HARM 250
C                   OF THE 3 DIMENSIONS OF THE ARRAY A.   THE SIZE,     HARM 260
C                   NI, OF THE I DIMENSION OF A IS 2**M(I), I = 1,2,3   HARM 270
C           INV   - A VECTOR WORK AREA FOR BIT AND INDEX MANIPULATION   HARM 280
C                   OF DIMENSION ONE FOURTH OF THE QUANTITY             HARM 290
C                   MAX(N1,N2,N3)                                       HARM 300
C           S     - A VECTOR WORK AREA FOR SINE TABLES WITH DIMENSION   HARM 310
C                   THE SAME AS INV                                     HARM 320
C           IFSET - AN OPTION PARAMETER WITH THE FOLLOWING SETTINGS     HARM 330
C                      0    SET UP SINE AND INV TABLES ONLY             HARM 340
C                      1    SET UP SINE AND INV TABLES ONLY AND         HARM 350
C                           CALCULATE FOURIER TRANSFORM                 HARM 360
C                     -1    SET UP SINE AND INV TABLES ONLY AND         HARM 370
C                           CALCULATE INVERSE FOURIER TRANSFORM (FOR    HARM 380
C                           THE MEANING OF INVERSE SEE THE EQUATIONS    HARM 390
C                           UNDER METHOD BELOW)                         HARM 400
C                      2    CALCULATE FOURIER TRANSFORM ONLY (ASSUME    HARM 410
C                           SINE AND INV TABLES EXIST)                  HARM 420
C                     -2    CALCULATE INVERSE FOURIER TRANSFORM ONLY    HARM 430
C                           (ASSUME SINE AND INV TABLES EXIST)          HARM 440
C           IFERR - ERROR INDICATOR.   WHEN IFSET IS 0,+1,-1,           HARM 450
C                   IFERR = 1 MEANS THE MAXIMUM M(I) IS GREATER THAN    HARM 460
C                  20 , I=1,2,3   WHEN IFSET IS 2,-2 , IFERR = 1        HARM 470
C                   MEANS THAT THE SINE AND INV TABLES ARE NOT LARGE    HARM 480
C                   ENOUGH OR HAVE NOT BEEN COMPUTED .                  HARM 490
C                   IF ON RETURN IFERR = 0 THEN NONE OF THE ABOVE       HARM 500
C                   CONDITIONS ARE PRESENT                              HARM 510
C                                                                       HARM 520
C        REMARKS                                                        HARM 530
C           THIS SUBROUTINE IS TO BE USED FOR COMPLEX, 3-DIMENSIONAL    HARM 540
C           ARRAYS IN WHICH EACH DIMENSION IS A POWER OF 2.  THE        HARM 550
C           MAXIMUM M(I) MUST NOT BE LESS THAN 3 OR GREATER THAN 20,    HARM 560
C           I = 1,2,3                                                   HARM 570
C                                                                       HARM 580
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  HARM 590
C           NONE                                                        HARM 600
C                                                                       HARM 610
C        METHOD                                                         HARM 620
C           FOR IFSET = +1, OR +2, THE FOURIER TRANSFORM OF COMPLEX     HARM 630
C           ARRAY A IS OBTAINED.                                        HARM 640
C                                                                       HARM 650
C                  N1-1   N2-1   N3-1                L1   L2   L3       HARM 660
C     X(J1,J2,J3)=SUM    SUM    SUM    A(K1,K2,K3)*W1  *W2  *W3         HARM 670
C                  K1=0   K2=0   K3=0                                   HARM 680
C                                                                       HARM 690
C                  WHERE WI IS THE N(I) ROOT OF UNITY AND L1=K1*J1,     HARM 700
C                        L2=K2*J2, L3=K3*J3                             HARM 710
C                                                                       HARM 720
C                                                                       HARM 730
C           FOR IFSET = -1, OR -2, THE INVERSE FOURIER TRANSFORM A OF   HARM 740
C           COMPLEX ARRAY X IS OBTAINED.                                HARM 750
C                                                                       HARM 760
C     A(K1,K2,K3)=                                                      HARM 770
C               1      N1-1   N2-1   N3-1                -L1  -L2  -L3  HARM 780
C           -------- *SUM    SUM    SUM    X(J1,J2,J3)*W1  *W2  *W3     HARM 790
C           N1*N2*N3   J1=0   J2=0   J3=0                               HARM 800
C                                                                       HARM 810
C                                                                       HARM 820
C           SEE J.W. COOLEY AND J.W. TUKEY, 'AN ALGORITHM FOR THE       HARM 830
C           MACHINE CALCULATION OF COMPLEX FOURIER SERIES',             HARM 840
C           MATHEMATICS OF COMPUTATIONS, VOL. 19 (APR. 1965), P. 297.   HARM 850
C                                                                       HARM 860
C     ..................................................................HARM 870
C                                                                       HARM 880
      SUBROUTINE HARM(A,M,INV,S,IFSET, IFERR)                           HARM 890
      DIMENSION A(1),INV(1),S(1),N(3),M(3),NP(3),W(2),W2(2),W3(2)       HARM 900
      EQUIVALENCE (N1,N(1)),(N2,N(2)),(N3,N(3))                         HARM 910
   10 IF( IABS(IFSET) - 1) 900,900,12                                   HARM 920
   12 MTT=MAX0(M(1),M(2),M(3)) -2                                       HARM 930
      ROOT2 = SQRT(2.)                                                  HARM 940
      IF (MTT-MT ) 14,14,13                                             HARM 950
   13 IFERR=1                                                           HARM 960
      RETURN                                                            HARM 970
   14 IFERR=0                                                           HARM 980
      M1=M(1)                                                           HARM 990
      M2=M(2)                                                           HARM1000
      M3=M(3)                                                           HARM1010
      N1=2**M1                                                          HARM1020
      N2=2**M2                                                          HARM1030
      N3=2**M3                                                          HARM1040
   16 IF(IFSET) 18,18,20                                                HARM1050
   18 NX= N1*N2*N3                                                      HARM1060
      FN = NX                                                           HARM1070
      DO 19 I = 1,NX                                                    HARM1080
      A(2*I-1) = A(2*I-1)/FN                                            HARM1090
   19 A(2*I) = -A(2*I)/FN                                               HARM1100
   20 NP(1)=N1*2                                                        HARM1110
      NP(2)= NP(1)*N2                                                   HARM1120
      NP(3)=NP(2)*N3                                                    HARM1130
      DO 250 ID=1,3                                                     HARM1140
      IL = NP(3)-NP(ID)                                                 HARM1150
      IL1 = IL+1                                                        HARM1160
      MI = M(ID)                                                        HARM1170
      IF (MI)250,250,30                                                 HARM1180
   30 IDIF=NP(ID)                                                       HARM1190
      KBIT=NP(ID)                                                       HARM1200
      MEV = 2*(MI/2)                                                    HARM1210
      IF (MI - MEV )60,60,40                                            HARM1220
C                                                                       HARM1230
C     M IS ODD. DO L=1 CASE                                             HARM1240
   40 KBIT=KBIT/2                                                       HARM1250
      KL=KBIT-2                                                         HARM1260
      DO 50 I=1,IL1,IDIF                                                HARM1270
      KLAST=KL+I                                                        HARM1280
      DO 50 K=I,KLAST,2                                                 HARM1290
      KD=K+KBIT                                                         HARM1300
C                                                                       HARM1310
C     DO ONE STEP WITH L=1,J=0                                          HARM1320
C     A(K)=A(K)+A(KD)                                                   HARM1330
C     A(KD)=A(K)-A(KD)                                                  HARM1340
C                                                                       HARM1350
      T=A(KD)                                                           HARM1360
      A(KD)=A(K)-T                                                      HARM1370
      A(K)=A(K)+T                                                       HARM1380
      T=A(KD+1)                                                         HARM1390
      A(KD+1)=A(K+1)-T                                                  HARM1400
   50 A(K+1)=A(K+1)+T                                                   HARM1410
      IF (MI - 1)250,250,52                                             HARM1420
   52 LFIRST =3                                                         HARM1430
C                                                                       HARM1440
C     DEF - JLAST = 2**(L-2) -1                                         HARM1450
      JLAST=1                                                           HARM1460
      GO TO 70                                                          HARM1470
C                                                                       HARM1480
C     M IS EVEN                                                         HARM1490
   60 LFIRST = 2                                                        HARM1500
      JLAST=0                                                           HARM1510
   70 DO 240 L=LFIRST,MI,2                                              HARM1520
      JJDIF=KBIT                                                        HARM1530
      KBIT=KBIT/4                                                       HARM1540
      KL=KBIT-2                                                         HARM1550
C                                                                       HARM1560
C     DO FOR J=0                                                        HARM1570
      DO 80 I=1,IL1,IDIF                                                HARM1580
      KLAST=I+KL                                                        HARM1590
      DO 80 K=I,KLAST,2                                                 HARM1600
      K1=K+KBIT                                                         HARM1610
      K2=K1+KBIT                                                        HARM1620
      K3=K2+KBIT                                                        HARM1630
C                                                                       HARM1640
C     DO TWO STEPS WITH J=0                                             HARM1650
C     A(K)=A(K)+A(K2)                                                   HARM1660
C     A(K2)=A(K)-A(K2)                                                  HARM1670
C     A(K1)=A(K1)+A(K3)                                                 HARM1680
C     A(K3)=A(K1)-A(K3)                                                 HARM1690
C                                                                       HARM1700
C     A(K)=A(K)+A(K1)                                                   HARM1710
C     A(K1)=A(K)-A(K1)                                                  HARM1720
C     A(K2)=A(K2)+A(K3)*I                                               HARM1730
C     A(K3)=A(K2)-A(K3)*I                                               HARM1740
C                                                                       HARM1750
      T=A(K2)                                                           HARM1760
      A(K2)=A(K)-T                                                      HARM1770
      A(K)=A(K)+T                                                       HARM1780
      T=A(K2+1)                                                         HARM1790
      A(K2+1)=A(K+1)-T                                                  HARM1800
      A(K+1)=A(K+1)+T                                                   HARM1810
C                                                                       HARM1820
      T=A(K3)                                                           HARM1830
      A(K3)=A(K1)-T                                                     HARM1840
      A(K1)=A(K1)+T                                                     HARM1850
      T=A(K3+1)                                                         HARM1860
      A(K3+1)=A(K1+1)-T                                                 HARM1870
      A(K1+1)=A(K1+1)+T                                                 HARM1880
C                                                                       HARM1890
      T=A(K1)                                                           HARM1900
      A(K1)=A(K)-T                                                      HARM1910
      A(K)=A(K)+T                                                       HARM1920
      T=A(K1+1)                                                         HARM1930
      A(K1+1)=A(K+1)-T                                                  HARM1940
      A(K+1)=A(K+1)+T                                                   HARM1950
C                                                                       HARM1960
      R=-A(K3+1)                                                        HARM1970
      T = A(K3)                                                         HARM1980
      A(K3)=A(K2)-R                                                     HARM1990
      A(K2)=A(K2)+R                                                     HARM2000
      A(K3+1)=A(K2+1)-T                                                 HARM2010
   80 A(K2+1)=A(K2+1)+T                                                 HARM2020
      IF (JLAST) 235,235,82                                             HARM2030
   82 JJ=JJDIF   +1                                                     HARM2040
C                                                                       HARM2050
C     DO FOR J=1                                                        HARM2060
      ILAST= IL +JJ                                                     HARM2070
      DO 85 I = JJ,ILAST,IDIF                                           HARM2080
      KLAST = KL+I                                                      HARM2090
      DO 85 K=I,KLAST,2                                                 HARM2100
      K1 = K+KBIT                                                       HARM2110
      K2 = K1+KBIT                                                      HARM2120
      K3 = K2+KBIT                                                      HARM2130
C                                                                       HARM2140
C     LETTING W=(1+I)/ROOT2,W3=(-1+I)/ROOT2,W2=I,                       HARM2150
C     A(K)=A(K)+A(K2)*I                                                 HARM2160
C     A(K2)=A(K)-A(K2)*I                                                HARM2170
C     A(K1)=A(K1)*W+A(K3)*W3                                            HARM2180
C     A(K3)=A(K1)*W-A(K3)*W3                                            HARM2190
C                                                                       HARM2200
C     A(K)=A(K)+A(K1)                                                   HARM2210
C     A(K1)=A(K)-A(K1)                                                  HARM2220
C     A(K2)=A(K2)+A(K3)*I                                               HARM2230
C     A(K3)=A(K2)-A(K3)*I                                               HARM2240
C                                                                       HARM2250
      R =-A(K2+1)                                                       HARM2260
      T = A(K2)                                                         HARM2270
      A(K2) = A(K)-R                                                    HARM2280
      A(K) = A(K)+R                                                     HARM2290
      A(K2+1)=A(K+1)-T                                                  HARM2300
      A(K+1)=A(K+1)+T                                                   HARM2310
C                                                                       HARM2320
      AWR=A(K1)-A(K1+1)                                                 HARM2330
      AWI = A(K1+1)+A(K1)                                               HARM2340
      R=-A(K3)-A(K3+1)                                                  HARM2350
      T=A(K3)-A(K3+1)                                                   HARM2360
      A(K3)=(AWR-R)/ROOT2                                               HARM2370
      A(K3+1)=(AWI-T)/ROOT2                                             HARM2380
      A(K1)=(AWR+R)/ROOT2                                               HARM2390
      A(K1+1)=(AWI+T)/ROOT2                                             HARM2400
      T= A(K1)                                                          HARM2410
      A(K1)=A(K)-T                                                      HARM2420
      A(K)=A(K)+T                                                       HARM2430
      T=A(K1+1)                                                         HARM2440
      A(K1+1)=A(K+1)-T                                                  HARM2450
      A(K+1)=A(K+1)+T                                                   HARM2460
      R=-A(K3+1)                                                        HARM2470
      T=A(K3)                                                           HARM2480
      A(K3)=A(K2)-R                                                     HARM2490
      A(K2)=A(K2)+R                                                     HARM2500
      A(K3+1)=A(K2+1)-T                                                 HARM2510
   85 A(K2+1)=A(K2+1)+T                                                 HARM2520
      IF(JLAST-1) 235,235,90                                            HARM2530
   90 JJ= JJ + JJDIF                                                    HARM2540
C                                                                       HARM2550
C     NOW DO THE REMAINING J'S                                          HARM2560
      DO 230 J=2,JLAST                                                  HARM2570
C                                                                       HARM2580
C     FETCH W'S                                                         HARM2590
C     DEF- W=W**INV(J), W2=W**2, W3=W**3                                HARM2600
   96 I=INV(J+1)                                                        HARM2610
   98 IC=NT-I                                                           HARM2620
      W(1)=S(IC)                                                        HARM2630
      W(2)=S(I)                                                         HARM2640
      I2=2*I                                                            HARM2650
      I2C=NT-I2                                                         HARM2660
      IF(I2C)120,110,100                                                HARM2670
C                                                                       HARM2680
C     2*I IS IN FIRST QUADRANT                                          HARM2690
  100 W2(1)=S(I2C)                                                      HARM2700
      W2(2)=S(I2)                                                       HARM2710
      GO TO 130                                                         HARM2720
  110 W2(1)=0.                                                          HARM2730
      W2(2)=1.                                                          HARM2740
      GO TO 130                                                         HARM2750
C                                                                       HARM2760
C     2*I IS IN SECOND QUADRANT                                         HARM2770
  120 I2CC = I2C+NT                                                     HARM2780
      I2C=-I2C                                                          HARM2790
      W2(1)=-S(I2C)                                                     HARM2800
      W2(2)=S(I2CC)                                                     HARM2810
  130 I3=I+I2                                                           HARM2820
      I3C=NT-I3                                                         HARM2830
      IF(I3C)160,150,140                                                HARM2840
C                                                                       HARM2850
C     I3 IN FIRST QUADRANT                                              HARM2860
  140 W3(1)=S(I3C)                                                      HARM2870
      W3(2)=S(I3)                                                       HARM2880
      GO TO 200                                                         HARM2890
  150 W3(1)=0.                                                          HARM2900
      W3(2)=1.                                                          HARM2910
      GO TO 200                                                         HARM2920
C                                                                       HARM2930
  160 I3CC=I3C+NT                                                       HARM2940
      IF(I3CC)190,180,170                                               HARM2950
C                                                                       HARM2960
C     I3 IN SECOND QUADRANT                                             HARM2970
  170 I3C=-I3C                                                          HARM2980
      W3(1)=-S(I3C)                                                     HARM2990
      W3(2)=S(I3CC)                                                     HARM3000
      GO TO 200                                                         HARM3010
  180 W3(1)=-1.                                                         HARM3020
      W3(2)=0.                                                          HARM3030
      GO TO 200                                                         HARM3040
C                                                                       HARM3050
C     3*I IN THIRD QUADRANT                                             HARM3060
  190 I3CCC=NT+I3CC                                                     HARM3070
      I3CC = -I3CC                                                      HARM3080
      W3(1)=-S(I3CCC)                                                   HARM3090
      W3(2)=-S(I3CC)                                                    HARM3100
  200 ILAST=IL+JJ                                                       HARM3110
      DO 220 I=JJ,ILAST,IDIF                                            HARM3120
      KLAST=KL+I                                                        HARM3130
      DO 220 K=I,KLAST,2                                                HARM3140
      K1=K+KBIT                                                         HARM3150
      K2=K1+KBIT                                                        HARM3160
      K3=K2+KBIT                                                        HARM3170
C                                                                       HARM3180
C     DO TWO STEPS WITH J NOT 0                                         HARM3190
C     A(K)=A(K)+A(K2)*W2                                                HARM3200
C     A(K2)=A(K)-A(K2)*W2                                               HARM3210
C     A(K1)=A(K1)*W+A(K3)*W3                                            HARM3220
C     A(K3)=A(K1)*W-A(K3)*W3                                            HARM3230
C                                                                       HARM3240
C     A(K)=A(K)+A(K1)                                                   HARM3250
C     A(K1)=A(K)-A(K1)                                                  HARM3260
C     A(K2)=A(K2)+A(K3)*I                                               HARM3270
C     A(K3)=A(K2)-A(K3)*I                                               HARM3280
C                                                                       HARM3290
      R=A(K2)*W2(1)-A(K2+1)*W2(2)                                       HARM3300
      T=A(K2)*W2(2)+A(K2+1)*W2(1)                                       HARM3310
      A(K2)=A(K)-R                                                      HARM3320
      A(K)=A(K)+R                                                       HARM3330
      A(K2+1)=A(K+1)-T                                                  HARM3340
      A(K+1)=A(K+1)+T                                                   HARM3350
C                                                                       HARM3360
      R=A(K3)*W3(1)-A(K3+1)*W3(2)                                       HARM3370
      T=A(K3)*W3(2)+A(K3+1)*W3(1)                                       HARM3380
      AWR=A(K1)*W(1)-A(K1+1)*W(2)                                       HARM3390
      AWI=A(K1)*W(2)+A(K1+1)*W(1)                                       HARM3400
      A(K3)=AWR-R                                                       HARM3410
      A(K3+1)=AWI-T                                                     HARM3420
      A(K1)=AWR+R                                                       HARM3430
      A(K1+1)=AWI+T                                                     HARM3440
      T=A(K1)                                                           HARM3450
      A(K1)=A(K)-T                                                      HARM3460
      A(K)=A(K)+T                                                       HARM3470
      T=A(K1+1)                                                         HARM3480
      A(K1+1)=A(K+1)-T                                                  HARM3490
      A(K+1)=A(K+1)+T                                                   HARM3500
      R=-A(K3+1)                                                        HARM3510
      T=A(K3)                                                           HARM3520
      A(K3)=A(K2)-R                                                     HARM3530
      A(K2)=A(K2)+R                                                     HARM3540
      A(K3+1)=A(K2+1)-T                                                 HARM3550
  220 A(K2+1)=A(K2+1)+T                                                 HARM3560
C     END OF I AND K LOOPS                                              HARM3570
C                                                                       HARM3580
  230 JJ=JJDIF+JJ                                                       HARM3590
C     END OF J-LOOP                                                     HARM3600
C                                                                       HARM3610
  235 JLAST=4*JLAST+3                                                   HARM3620
  240 CONTINUE                                                          HARM3630
C     END OF  L  LOOP                                                   HARM3640
C                                                                       HARM3650
  250 CONTINUE                                                          HARM3660
C     END OF  ID  LOOP                                                  HARM3670
C                                                                       HARM3680
C     WE NOW HAVE THE COMPLEX FOURIER SUMS BUT THEIR ADDRESSES ARE      HARM3690
C     BIT-REVERSED.  THE FOLLOWING ROUTINE PUTS THEM IN ORDER           HARM3700
      NTSQ=NT*NT                                                        HARM3710
      M3MT=M3-MT                                                        HARM3720
  350 IF(M3MT) 370,360,360                                              HARM3730
C                                                                       HARM3740
C     M3 GR. OR EQ. MT                                                  HARM3750
  360 IGO3=1                                                            HARM3760
      N3VNT=N3/NT                                                       HARM3770
      MINN3=NT                                                          HARM3780
      GO TO 380                                                         HARM3790
C                                                                       HARM3800
C     M3 LESS THAN MT                                                   HARM3810
  370 IGO3=2                                                            HARM3820
      N3VNT=1                                                           HARM3830
      NTVN3=NT/N3                                                       HARM3840
      MINN3=N3                                                          HARM3850
  380 JJD3 = NTSQ/N3                                                    HARM3860
      M2MT=M2-MT                                                        HARM3870
  450 IF (M2MT)470,460,460                                              HARM3880
C                                                                       HARM3890
C     M2 GR. OR EQ. MT                                                  HARM3900
  460 IGO2=1                                                            HARM3910
      N2VNT=N2/NT                                                       HARM3920
      MINN2=NT                                                          HARM3930
      GO TO 480                                                         HARM3940
C                                                                       HARM3950
C     M2 LESS THAN MT                                                   HARM3960
  470 IGO2 = 2                                                          HARM3970
      N2VNT=1                                                           HARM3980
      NTVN2=NT/N2                                                       HARM3990
      MINN2=N2                                                          HARM4000
  480 JJD2=NTSQ/N2                                                      HARM4010
      M1MT=M1-MT                                                        HARM4020
  550 IF(M1MT)570,560,560                                               HARM4030
C                                                                       HARM4040
C     M1 GR. OR EQ. MT                                                  HARM4050
  560 IGO1=1                                                            HARM4060
      N1VNT=N1/NT                                                       HARM4070
      MINN1=NT                                                          HARM4080
      GO TO 580                                                         HARM4090
C                                                                       HARM4100
C     M1 LESS THAN MT                                                   HARM4110
  570 IGO1=2                                                            HARM4120
      N1VNT=1                                                           HARM4130
      NTVN1=NT/N1                                                       HARM4140
      MINN1=N1                                                          HARM4150
  580 JJD1=NTSQ/N1                                                      HARM4160
  600 JJ3=1                                                             HARM4170
      J=1                                                               HARM4180
      DO 880 JPP3=1,N3VNT                                               HARM4190
      IPP3=INV(JJ3)                                                     HARM4200
      DO 870 JP3=1,MINN3                                                HARM4210
      GO TO (610,620),IGO3                                              HARM4220
  610 IP3=INV(JP3)*N3VNT                                                HARM4230
      GO TO 630                                                         HARM4240
  620 IP3=INV(JP3)/NTVN3                                                HARM4250
  630 I3=(IPP3+IP3)*N2                                                  HARM4260
  700 JJ2=1                                                             HARM4270
      DO 870 JPP2=1,N2VNT                                               HARM4280
      IPP2=INV(JJ2)+I3                                                  HARM4290
      DO 860 JP2=1,MINN2                                                HARM4300
      GO TO (710,720),IGO2                                              HARM4310
  710 IP2=INV(JP2)*N2VNT                                                HARM4320
      GO TO 730                                                         HARM4330
  720 IP2=INV(JP2)/NTVN2                                                HARM4340
  730 I2=(IPP2+IP2)*N1                                                  HARM4350
  800 JJ1=1                                                             HARM4360
      DO 860 JPP1=1,N1VNT                                               HARM4370
      IPP1=INV(JJ1)+I2                                                  HARM4380
      DO 850 JP1=1,MINN1                                                HARM4390
      GO TO (810,820),IGO1                                              HARM4400
  810 IP1=INV(JP1)*N1VNT                                                HARM4410
      GO TO 830                                                         HARM4420
  820 IP1=INV(JP1)/NTVN1                                                HARM4430
  830 I=2*(IPP1+IP1)+1                                                  HARM4440
      IF (J-I) 840,850,850                                              HARM4450
  840 T=A(I)                                                            HARM4460
      A(I)=A(J)                                                         HARM4470
      A(J)=T                                                            HARM4480
      T=A(I+1)                                                          HARM4490
      A(I+1)=A(J+1)                                                     HARM4500
      A(J+1)=T                                                          HARM4510
  850 J=J+2                                                             HARM4530
  860 JJ1=JJ1+JJD1                                                      HARM4540
C     END OF JPP1 AND JP2                                               HARM4550
C                                                                       HARM4560
  870 JJ2=JJ2+JJD2                                                      HARM4570
C     END OF JPP2 AND JP3 LOOPS                                         HARM4580
C                                                                       HARM4590
  880 JJ3 = JJ3+JJD3                                                    HARM4600
C     END OF JPP3 LOOP                                                  HARM4610
C                                                                       HARM4620
  890 IF(IFSET)891,895,895                                              HARM4630
  891 DO 892 I = 1,NX                                                   HARM4640
  892 A(2*I) = -A(2*I)                                                  HARM4650
  895 RETURN                                                            HARM4660
C                                                                       HARM4670
C     THE FOLLOWING PROGRAM COMPUTES THE SIN AND INV TABLES.            HARM4680
C                                                                       HARM4690
  900 MT=MAX0(M(1),M(2),M(3)) -2                                        HARM4700
      MT = MAX0(2,MT)                                                   HARM4710
  904 IF (MT-18) 906,906,13                                             HARM4720
  906 IFERR=0                                                           HARM4750
      NT=2**MT                                                          HARM4760
      NTV2=NT/2                                                         HARM4770
C                                                                       HARM4780
C     SET UP SIN TABLE                                                  HARM4790
C     THETA=PIE/2**(L+1) FOR L=1                                        HARM4800
  910 THETA=.7853981634                                                 HARM4810
C                                                                       HARM4820
C     JSTEP=2**(MT-L+1) FOR L=1                                         HARM4830
      JSTEP=NT                                                          HARM4840
C                                                                       HARM4850
C     JDIF=2**(MT-L) FOR L=1                                            HARM4860
      JDIF=NTV2                                                         HARM4880
      S(JDIF)=SIN(THETA)                                                HARM4890
      DO 950 L=2,MT                                                     HARM4900
      THETA=THETA/2.                                                    HARM4910
      JSTEP2=JSTEP                                                      HARM4920
      JSTEP=JDIF                                                        HARM4930
      JDIF=JSTEP/2                                                      HARM4940
      S(JDIF)=SIN(THETA)                                                HARM4950
      JC1=NT-JDIF                                                       HARM4960
      S(JC1)=COS(THETA)                                                 HARM4970
      JLAST=NT-JSTEP2                                                   HARM4980
      IF(JLAST - JSTEP) 950,920,920                                     HARM4990
  920 DO 940 J=JSTEP,JLAST,JSTEP                                        HARM5000
      JC=NT-J                                                           HARM5010
      JD=J+JDIF                                                         HARM5020
  940 S(JD)=S(J)*S(JC1)+S(JDIF)*S(JC)                                   HARM5030
  950 CONTINUE                                                          HARM5040
C                                                                       HARM5050
C     SET UP INV(J) TABLE                                               HARM5060
C                                                                       HARM5070
  960 MTLEXP=NTV2                                                       HARM5080
C                                                                       HARM5090
C     MTLEXP=2**(MT-L). FOR L=1                                         HARM5100
      LM1EXP=1                                                          HARM5110
C                                                                       HARM5120
C     LM1EXP=2**(L-1). FOR L=1                                          HARM5130
      INV(1)=0                                                          HARM5140
      DO 980 L=1,MT                                                     HARM5150
      INV(LM1EXP+1) = MTLEXP                                            HARM5160
      DO 970 J=2,LM1EXP                                                 HARM5170
      JJ=J+LM1EXP                                                       HARM5180
  970 INV(JJ)=INV(J)+MTLEXP                                             HARM5190
      MTLEXP=MTLEXP/2                                                   HARM5200
  980 LM1EXP=LM1EXP*2                                                   HARM5210
  982 IF(IFSET)12,895,12                                                HARM5220
      END                                                               HARM5230
