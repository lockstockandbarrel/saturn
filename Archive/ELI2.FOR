C                                                                       ELI2  10
C     ..................................................................ELI2  20
C                                                                       ELI2  30
C        SUBROUTINE ELI2                                                ELI2  40
C                                                                       ELI2  50
C        PURPOSE                                                        ELI2  60
C           COMPUTES THE GENERALIZED ELLIPTIC INTEGRAL OF SECOND KIND   ELI2  70
C                                                                       ELI2  80
C        USAGE                                                          ELI2  90
C           CALL ELI2(R,X,CK,A,B)                                       ELI2 100
C                                                                       ELI2 110
C        DESCRIPTION OF PARAMETERS                                      ELI2 120
C           R     - RESULT VALUE                                        ELI2 130
C           X     - UPPER INTEGRATION BOUND (ARGUMENT OF ELLIPTIC       ELI2 140
C                   INTEGRAL OF SECOND KIND)                            ELI2 150
C           CK    - COMPLEMENTARY MODULUS                               ELI2 160
C           A     - CONSTANT TERM IN NUMERATOR                          ELI2 170
C           B     - QUADRATIC TERM IN NUMERATOR                         ELI2 180
C                                                                       ELI2 190
C        REMARKS                                                        ELI2 200
C           MODULUS K = SQRT(1.-CK*CK).                                 ELI2 210
C           SPECIAL CASES OF THE GENERALIZED ELLIPTIC INTEGRAL OF       ELI2 220
C           SECOND KIND ARE                                             ELI2 230
C           F(ATAN(X),K) OBTAINED WITH A=1., B=1.                       ELI2 240
C           E(ATAN(X),K) OBTAINED WITH A=1., B=CK*CK.                   ELI2 250
C           B(ATAN(X),K) OBTAINED WITH A=1., B=0.                       ELI2 260
C           D(ATAN(X),K) OBTAINED WITH A=0., B=1.                       ELI2 270
C                                                                       ELI2 280
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  ELI2 290
C           NONE                                                        ELI2 300
C                                                                       ELI2 310
C        METHOD                                                         ELI2 320
C           DEFINITION                                                  ELI2 330
C           R=INTEGRAL((A+B*T*T)/(SQRT((1+T*T)*(1+(CK*T)**2))*(1+T*T)), ELI2 340
C                  SUMMED OVER T FROM 0 TO X).                          ELI2 350
C           EQUIVALENT IS THE DEFINITION                                ELI2 360
C           R=INTEGRAL((A+(B-A)*(SIN(T))**2)/SQRT(1-(K*SIN(T))**2),     ELI2 370
C                  SUMMED OVER T FROM 0 TO ATAN(X)).                    ELI2 380
C           EVALUATION                                                  ELI2 390
C           LANDENS TRANSFORMATION IS USED FOR CALCULATION.             ELI2 400
C           REFERENCE                                                   ELI2 410
C           R. BULIRSCH, NUMERICAL CALCULATION OF ELLIPTIC INTEGRALS ANDELI2 420
C                  ELLIPTIC FUNCTIONS                                   ELI2 430
C                  HANDBOOK SERIES OF SPECIAL FUNCTIONS                 ELI2 440
C                  NUMERISCHE MATHEMATIK VOL. 7, 1965, PP. 78-90.       ELI2 450
C                                                                       ELI2 460
C     ..................................................................ELI2 470
C                                                                       ELI2 480
      SUBROUTINE ELI2(R,X,CK,A,B)                                       ELI2 490
C        TEST ARGUMENT                                                  ELI2 500
      IF(X)2,1,2                                                        ELI2 510
    1 R=0.                                                              ELI2 520
      RETURN                                                            ELI2 530
C        TEST MODULUS                                                   ELI2 540
    2 C=0.                                                              ELI2 550
      D=0.5                                                             ELI2 560
      IF(CK)7,3,7                                                       ELI2 570
    3 R=SQRT(1.+X*X)                                                    ELI2 580
      R=(A-B)*ABS(X)/R+B*ALOG(ABS(X)+R)                                 ELI2 590
C        TEST SIGN OF ARGUMENT                                          ELI2 600
    4 R=R+C*(A-B)                                                       ELI2 610
      IF(X)5,6,6                                                        ELI2 620
    5 R=-R                                                              ELI2 630
    6 RETURN                                                            ELI2 640
C        INITIALIZATION                                                 ELI2 650
    7 AN=(B+A)*0.5                                                      ELI2 660
      AA=A                                                              ELI2 670
      R=B                                                               ELI2 680
      ANG=ABS(1./X)                                                     ELI2 690
      PIM=0.                                                            ELI2 700
      ISI=0                                                             ELI2 710
      ARI=1.                                                            ELI2 720
      GEO=ABS(CK)                                                       ELI2 730
C        LANDEN TRANSFORMATION                                          ELI2 740
    8 R=AA*GEO+R                                                        ELI2 750
      SGEO=ARI*GEO                                                      ELI2 760
      AA=AN                                                             ELI2 770
      AARI=ARI                                                          ELI2 780
C        ARITHMETIC MEAN                                                ELI2 790
      ARI=GEO+ARI                                                       ELI2 800
C        SUM OF SINE VALUES                                             ELI2 810
      AN=(R/ARI+AA)*0.5                                                 ELI2 820
      AANG=ABS(ANG)                                                     ELI2 830
      ANG=-SGEO/ANG+ANG                                                 ELI2 840
      PIMA=PIM                                                          ELI2 850
      IF(ANG)10,9,11                                                    ELI2 860
    9 ANG=-1.E-8*AANG                                                   ELI2 870
   10 PIM=PIM+3.1415927                                                 ELI2 880
      ISI=ISI+1                                                         ELI2 890
   11 AANG=ARI*ARI+ANG*ANG                                              ELI2 900
      P=D/SQRT(AANG)                                                    ELI2 910
      IF(ISI-4)13,12,12                                                 ELI2 920
   12 ISI=ISI-4                                                         ELI2 930
   13 IF(ISI-2)15,14,14                                                 ELI2 940
   14 P=-P                                                              ELI2 950
   15 C=C+P                                                             ELI2 960
      D=D*(AARI-GEO)*0.5/ARI                                            ELI2 970
      IF(ABS(AARI-GEO)-1.E-4*AARI)17,17,16                              ELI2 980
   16 SGEO=SQRT(SGEO)                                                   ELI2 990
C        GEOMETRIC MEAN                                                 ELI21000
      GEO=SGEO+SGEO                                                     ELI21010
      PIM=PIM+PIMA                                                      ELI21020
      ISI=ISI+ISI                                                       ELI21030
      GOTO 8                                                            ELI21040
C        ACCURACY WAS SUFFICIENT                                        ELI21050
   17 R=(ATAN(ARI/ANG)+PIM)*AN/ARI                                      ELI21060
      C=C+D*ANG/AANG                                                    ELI21070
      GOTO 4                                                            ELI21080
      END                                                               ELI21090
