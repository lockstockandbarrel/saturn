C                                                                       JELF  10
C     ..................................................................JELF  20
C                                                                       JELF  30
C        SUBROUTINE JELF                                                JELF  40
C                                                                       JELF  50
C        PURPOSE                                                        JELF  60
C           COMPUTES THE THREE JACOBIAN ELLIPTIC FUNCTIONS SN, CN, DN.  JELF  70
C                                                                       JELF  80
C        USAGE                                                          JELF  90
C           CALL JELF(SN,CN,DN,X,SCK)                                   JELF 100
C                                                                       JELF 110
C        DESCRIPTION OF PARAMETERS                                      JELF 120
C           SN    - RESULT VALUE SN(X)                                  JELF 130
C           CN    - RESULT VALUE CN(X)                                  JELF 140
C           DN    - RESULT VALUE DN(X)                                  JELF 150
C           X     - ARGUMENT OF JACOBIAN ELLIPTIC FUNCTIONS             JELF 160
C           SCK   - SQUARE OF COMPLEMENTARY MODULUS                     JELF 170
C                                                                       JELF 180
C        REMARKS                                                        JELF 190
C           NONE                                                        JELF 200
C                                                                       JELF 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  JELF 220
C           NONE                                                        JELF 230
C                                                                       JELF 240
C        METHOD                                                         JELF 250
C           DEFINITION                                                  JELF 260
C           X=INTEGRAL(1/SQRT((1-T*T)*(1-(K*T)**2)), SUMMED OVER        JELF 270
C           T FROM 0 TO SN), WHERE K=SQRT(1-SCK).                       JELF 280
C           SN*SN + CN*CN = 1                                           JELF 290
C           (K*SN)**2 + DN**2 = 1.                                      JELF 300
C           EVALUATION                                                  JELF 310
C           CALCULATION IS DONE USING THE PROCESS OF THE ARITHMETIC     JELF 320
C           GEOMETRIC MEAN TOGETHER WITH GAUSS DESCENDING TRANSFORMATIONJELF 330
C           BEFORE INVERSION OF THE INTEGRAL TAKES PLACE.               JELF 340
C           REFERENCE                                                   JELF 350
C           R. BULIRSCH, NUMERICAL CALCULATION OF ELLIPTIC INTEGRALS ANDJELF 360
C                  ELLIPTIC FUNCTIOMS.                                  JELF 370
C                  HANDBOOK SERIES OF SPECIAL FUNCTIONS                 JELF 380
C                  NUMERISCHE MATHEMATIK VOL. 7, 1965, PP. 78-90.       JELF 390
C                                                                       JELF 400
C     ..................................................................JELF 410
C                                                                       JELF 420
      SUBROUTINE JELF(SN,CN,DN,X,SCK)                                   JELF 430
C                                                                       JELF 440
C                                                                       JELF 450
      DIMENSION ARI(12),GEO(12)                                         JELF 460
C     TEST MODULUS                                                      JELF 470
      CM=SCK                                                            JELF 480
      Y=X                                                               JELF 490
      IF(SCK)3,1,4                                                      JELF 500
    1 D=EXP(X)                                                          JELF 510
      A=1./D                                                            JELF 520
      B=A+D                                                             JELF 530
      CN=2./B                                                           JELF 540
      DN=CN                                                             JELF 550
      SN=TANH(X)                                                        JELF 560
C        DEGENERATE CASE SCK=0 GIVES RESULTS                            JELF 570
C           CN X = DN X = 1/COSH X                                      JELF 580
C           SN X = TANH X                                               JELF 590
    2 RETURN                                                            JELF 600
C        JACOBIS MODULUS TRANSFORMATION                                 JELF 610
    3 D=1.-SCK                                                          JELF 620
      CM=-SCK/D                                                         JELF 630
      D=SQRT(D)                                                         JELF 640
      Y=D*X                                                             JELF 650
    4 A=1.                                                              JELF 660
      DN=1.                                                             JELF 670
      DO 6 I=1,12                                                       JELF 680
      L=I                                                               JELF 690
      ARI(I)=A                                                          JELF 700
      CM=SQRT(CM)                                                       JELF 710
      GEO(I)=CM                                                         JELF 720
      C=(A+CM)*.5                                                       JELF 730
      IF(ABS(A-CM)-1.E-4*A)7,7,5                                        JELF 740
    5 CM=A*CM                                                           JELF 750
    6 A=C                                                               JELF 760
C                                                                       JELF 770
C     START BACKWARD RECURSION                                          JELF 780
    7 Y=C*Y                                                             JELF 790
      SN=SIN(Y)                                                         JELF 800
      CN=COS(Y)                                                         JELF 810
      IF(SN)8,13,8                                                      JELF 820
    8 A=CN/SN                                                           JELF 830
      C=A*C                                                             JELF 840
      DO 9 I=1,L                                                        JELF 850
      K=L-I+1                                                           JELF 860
      B=ARI(K)                                                          JELF 870
      A=C*A                                                             JELF 880
      C=DN*C                                                            JELF 890
      DN=(GEO(K)+A)/(B+A)                                               JELF 900
    9 A=C/B                                                             JELF 910
      A=1./SQRT(C*C+1.)                                                 JELF 920
      IF(SN)10,11,11                                                    JELF 930
   10 SN=-A                                                             JELF 940
      GOTO 12                                                           JELF 950
   11 SN=A                                                              JELF 960
   12 CN=C*SN                                                           JELF 970
   13 IF(SCK)14,2,2                                                     JELF 980
   14 A=DN                                                              JELF 990
      DN=CN                                                             JELF1000
      CN=A                                                              JELF1010
      SN=SN/D                                                           JELF1020
      RETURN                                                            JELF1030
      END                                                               JELF1040
