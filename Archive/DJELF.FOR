C                                                                       DJEL  10
C     ..................................................................DJEL  20
C                                                                       DJEL  30
C        SUBROUTINE DJELF                                               DJEL  40
C                                                                       DJEL  50
C        PURPOSE                                                        DJEL  60
C           COMPUTES THE THREE JACOBIAN ELLIPTIC FUNCTIONS SN, CN, DN.  DJEL  70
C                                                                       DJEL  80
C        USAGE                                                          DJEL  90
C           CALL DJELF(SN,CN,DN,X,SCK)                                  DJEL 100
C                                                                       DJEL 110
C        DESCRIPTION OF PARAMETERS                                      DJEL 120
C           SN    - RESULT VALUE SN(X) IN DOUBLE PRECISION              DJEL 130
C           CN    - RESULT VALUE CN(X) IN DOUBLE PRECISION              DJEL 140
C           DN    - RESULT VALUE DN(X) IN DOUBLE PRECISION              DJEL 150
C           X     - DOUBLE PRECISION ARGUMENT OF JACOBIAN ELLIPTIC      DJEL 160
C                   FUNCTIONS                                           DJEL 170
C           SCK   - SQUARE OF COMPLEMENTARY MODULUS IN DOUBLE PRECISION DJEL 180
C                                                                       DJEL 190
C        REMARKS                                                        DJEL 200
C           NONE                                                        DJEL 210
C                                                                       DJEL 220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DJEL 230
C           NONE                                                        DJEL 240
C                                                                       DJEL 250
C        METHOD                                                         DJEL 260
C           DEFINITION                                                  DJEL 270
C           X=INTEGRAL(1/SQRT((1-T*T)*(1-(K*T)**2)), SUMMED OVER        DJEL 280
C           T FROM 0 TO SN), WHERE K=SQRT(1-SCK).                       DJEL 290
C           SN*SN + CN*CN = 1                                           DJEL 300
C           (K*SN)**2 + DN**2 = 1.                                      DJEL 310
C           EVALUATION                                                  DJEL 320
C           CALCULATION IS DONE USING THE PROCESS OF THE ARITHMETIC     DJEL 330
C           GEOMETRIC MEAN TOGETHER WITH GAUSS DESCENDING TRANSFORMATIONDJEL 340
C           BEFORE INVERSION OF THE INTEGRAL TAKES PLACE.               DJEL 350
C           REFERENCE                                                   DJEL 360
C           R. BULIRSCH, NUMERICAL CALCULATION OF ELLIPTIC INTEGRALS ANDDJEL 370
C                  ELLIPTIC FUNCTIOMS.                                  DJEL 380
C                  HANDBOOK SERIES OF SPECIAL FUNCTIONS                 DJEL 390
C                  NUMERISCHE MATHEMATIK VOL. 7, 1965, PP. 78-90.       DJEL 400
C                                                                       DJEL 410
C     ..................................................................DJEL 420
C                                                                       DJEL 430
      SUBROUTINE DJELF(SN,CN,DN,X,SCK)                                  DJEL 440
C                                                                       DJEL 450
      DIMENSION ARI(12),GEO(12)                                         DJEL 460
      DOUBLE PRECISION SN,CN,DN,X,SCK,ARI,GEO,CM,Y,A,B,C,D              DJEL 470
C                                                                       DJEL 480
C        TEST MODULUS                                                   DJEL 490
C                                                                       DJEL 500
      CM=SCK                                                            DJEL 510
      Y=X                                                               DJEL 520
      IF(SCK)3,1,4                                                      DJEL 530
    1 D=DEXP(X)                                                         DJEL 540
      A=1.D0/D                                                          DJEL 550
      B=A+D                                                             DJEL 560
      CN=2.D0/B                                                         DJEL 570
      DN=CN                                                             DJEL 580
      A=(D-A)/2.D0                                                      DJEL 590
      SN=A*CN                                                           DJEL 600
C        DEGENERATE CASE SCK=0 GIVES RESULTS                            DJEL 610
C           CN X = DN X = 1/COSH X                                      DJEL 620
C           SN X = TANH X                                               DJEL 630
    2 RETURN                                                            DJEL 640
C                                                                       DJEL 650
C        JACOBIS MODULUS TRANSFORMATION                                 DJEL 660
C                                                                       DJEL 670
    3 D=1.D0-SCK                                                        DJEL 680
      CM=-SCK/D                                                         DJEL 690
      D=DSQRT(D)                                                        DJEL 700
      Y=D*X                                                             DJEL 710
    4 A=1.D0                                                            DJEL 720
      DN=1.D0                                                           DJEL 730
      DO 6 I=1,12                                                       DJEL 740
      L=I                                                               DJEL 750
      ARI(I)=A                                                          DJEL 760
      CM=DSQRT(CM)                                                      DJEL 770
      GEO(I)=CM                                                         DJEL 780
      C=(A+CM)*.5D0                                                     DJEL 790
      IF(DABS(A-CM)-1.D-9*A)7,7,5                                       DJEL 800
    5 CM=A*CM                                                           DJEL 810
    6 A=C                                                               DJEL 820
C                                                                       DJEL 830
C        START BACKWARD RECURSION                                       DJEL 840
C                                                                       DJEL 850
    7 Y=C*Y                                                             DJEL 860
      SN=DSIN(Y)                                                        DJEL 870
      CN=DCOS(Y)                                                        DJEL 880
      IF(SN)8,13,8                                                      DJEL 890
    8 A=CN/SN                                                           DJEL 900
      C=A*C                                                             DJEL 910
      DO 9 I=1,L                                                        DJEL 920
      K=L-I+1                                                           DJEL 930
      B=ARI(K)                                                          DJEL 940
      A=C*A                                                             DJEL 950
      C=DN*C                                                            DJEL 960
      DN=(GEO(K)+A)/(B+A)                                               DJEL 970
    9 A=C/B                                                             DJEL 980
      A=1.D0/DSQRT(C*C+1.D0)                                            DJEL 990
      IF(SN)10,11,11                                                    DJEL1000
   10 SN=-A                                                             DJEL1010
      GOTO 12                                                           DJEL1020
   11 SN=A                                                              DJEL1030
   12 CN=C*SN                                                           DJEL1040
   13 IF(SCK)14,2,2                                                     DJEL1050
   14 A=DN                                                              DJEL1060
      DN=CN                                                             DJEL1070
      CN=A                                                              DJEL1080
      SN=SN/D                                                           DJEL1090
      RETURN                                                            DJEL1100
      END                                                               DJEL1110
