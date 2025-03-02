C                                                                       RHAR  10
C     ..................................................................RHAR  20
C                                                                       RHAR  30
C        SUBROUTINE RHARM                                               RHAR  40
C                                                                       RHAR  50
C        PURPOSE                                                        RHAR  60
C           FINDS THE FOURIER COEFFICIENTS OF ONE DIMENSIONAL REAL DATA RHAR  70
C                                                                       RHAR  80
C        USAGE                                                          RHAR  90
C           CALL RHARM (A,M,INV,S,IFERR)                                RHAR 100
C                                                                       RHAR 110
C        DESCRIPTION OF PARAMETERS                                      RHAR 120
C           A     - AS INPUT, CONTAINS ONE DIMENSIONAL REAL DATA. A IS  RHAR 130
C                   2*N+4 CORE LOCATIONS, WHERE N = 2**M. 2*N REAL      RHAR 140
C                   NUMBERS ARE PUT INTO THE FIRST 2*N CORE LOCATIONS   RHAR 150
C                   OF A                                                RHAR 160
C                   AS OUTPUT, A CONTAINS THE FOURIER COEFFICIENTS      RHAR 170
C                   A0/2,B0=0,A1,B1,A2,B2,...,AN/2,BN=0 RESPECTIVELY IN RHAR 180
C                   THE FIRST 2N+2 CORE LOCATIONS OF A                  RHAR 190
C           M     - AN INTEGER WHICH DETERMINES THE SIZE OF THE VECTOR  RHAR 200
C                   A. THE SIZE OF A IS 2*(2**M) + 4                    RHAR 210
C           INV   - A VECTOR WORK AREA FOR BIT AND INDEX MANIPULATION OFRHAR 220
C                   DIMENSION ONE EIGHTH THE NUMBER OF REAL INPUT, VIZ.,RHAR 230
C                   (1/8)*2*(2**M)                                      RHAR 240
C           S     - A VECTOR WORK AREA FOR SINE TABLES WITH DIMENSION   RHAR 250
C                   THE SAME AS INV                                     RHAR 260
C           IFERR - A RETURNED VALUE OF 1 MEANS THAT M IS LESS THAN 3 ORRHAR 270
C                   GREATER THAN 20. OTHERWISE IFERR IS SET = 0         RHAR 280
C                                                                       RHAR 290
C        REMARKS                                                        RHAR 300
C           THIS SUBROUTINE GIVES THE FOURIER COEFFICIENTS OF 2*(2**M)  RHAR 310
C           REAL POINTS. SEE SUBROUTINE HARM FOR THREE DIMENSIONAL,     RHAR 320
C           COMPLEX FOURIER TRANSFORMS                                  RHAR 330
C                                                                       RHAR 340
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  RHAR 350
C           HARM                                                        RHAR 360
C                                                                       RHAR 370
C        METHOD                                                         RHAR 380
C           THE FOURIER COEFFICIENTS A0,B0=0,A1,B1,...,AN,BN=0 ARE      RHAR 390
C           OBTAINED FOR INPUT XJ, J=0,1,2,...,2N-1 FOR THE FOLLOWING   RHAR 400
C           EQUATION (PI = 3.14159...)                                  RHAR 410
C                                                                       RHAR 420
C                 N-1                                               J   RHAR 430
C     XJ=(1/2)A0+SUM (AK*COS(PI*J*K/N)+BK*SIN(PI*J*K/N))+(1/2)AN(-1)    RHAR 440
C                 K=1                                                   RHAR 450
C                                                                       RHAR 460
C        SEE REFERENCE UNDER SUBROUTINE HARM                            RHAR 470
C                                                                       RHAR 480
C     ..................................................................RHAR 490
C                                                                       RHAR 500
      SUBROUTINE RHARM(A,M,INV,S,IFERR)                                 RHAR 510
      DIMENSION A(1),L(3),INV(1),S(1)                                   RHAR 520
      IFSET=1                                                           RHAR 530
      L(1)=M                                                            RHAR 540
      L(2)=0                                                            RHAR 550
      L(3)=0                                                            RHAR 560
      NTOT=2**M                                                         RHAR 570
      NTOT2 = 2*NTOT                                                    RHAR 580
      FN = NTOT                                                         RHAR 590
      DO   3 I = 2,NTOT2,2                                              RHAR 600
   3  A(I) = -A(I)                                                      RHAR 610
      DO   6 I = 1,NTOT2                                                RHAR 620
   6  A(I) = A(I)/FN                                                    RHAR 630
      CALL HARM(A,L,INV,S ,IFSET,IFERR)                                 RHAR 640
C                                                                       RHAR 650
C     MOVE LAST HALF OF A(J)S DOWN ONE SLOT AND ADD A(N) AT BOTTOM TO   RHAR 660
C     GIVE ARRAY FOR A1PRIME AND A2PRIME CALCULATION                    RHAR 670
C                                                                       RHAR 680
   21 DO  52 I=1,NTOT,2                                                 RHAR 690
      J0=NTOT2+2-I                                                      RHAR 700
      A(J0)=A(J0-2)                                                     RHAR 710
   52 A(J0+1)=A(J0-1)                                                   RHAR 720
      A(NTOT2+3)=A(1)                                                   RHAR 730
      A(NTOT2+4)=A(2)                                                   RHAR 740
C                                                                       RHAR 750
C     CALCULATE A1PRIMES AND STORE IN FIRST N SLOTS                     RHAR 760
C     CALCULATE A2PRIMES AND STORE IN SECOND N SLOTS IN REVERSE ORDER   RHAR 770
      K0=NTOT+1                                                         RHAR 780
      DO 104 I=1,K0,2                                                   RHAR 790
      K1=NTOT2-I+4                                                      RHAR 800
      AP1RE=.5*(A(I)+A(K1))                                             RHAR 810
      AP2RE=-.5*(A(I+1)+A(K1+1))                                        RHAR 820
      AP1IM=.5*(-A(I+1)+A(K1+1))                                        RHAR 830
      AP2IM=-.5*(A(I)-A(K1))                                            RHAR 840
      A(I)=AP1RE                                                        RHAR 850
      A(I+1)=AP1IM                                                      RHAR 860
      A(K1)=AP2RE                                                       RHAR 870
  104 A(K1+1)=AP2IM                                                     RHAR 880
      NTO = NTOT/2                                                      RHAR 890
  110 NT=NTO+1                                                          RHAR 900
      DEL=3.1415927/FLOAT(NTOT)                                         RHAR 910
      SS=SIN(DEL)                                                       RHAR 920
      SC=COS(DEL)                                                       RHAR 930
      SI=0.0                                                            RHAR 940
      CO=1.0                                                            RHAR 950
C                                                                       RHAR 960
C     COMPUTE C(J)S FOR J=0 THRU J=N                                    RHAR 970
  114 DO 116 I=1,NT                                                     RHAR 980
      K6=NTOT2-2*I+5                                                    RHAR 990
      AP2RE=A(K6)*CO+A(K6+1)*SI                                         RHAR1000
      AP2IM=-A(K6)*SI+A(K6+1)*CO                                        RHAR1010
      CIRE=.5*(A(2*I-1)+AP2RE)                                          RHAR1020
      CIIM=.5*(A(2*I)+AP2IM)                                            RHAR1030
      CNIRE=.5*(A(2*I-1)-AP2RE)                                         RHAR1040
      CNIIM=.5*(A(2*I)-AP2IM)                                           RHAR1050
      A(2*I-1)=CIRE                                                     RHAR1060
      A(2*I)=CIIM                                                       RHAR1070
      A(K6)=CNIRE                                                       RHAR1080
      A(K6+1)=-CNIIM                                                    RHAR1090
      SIS=SI                                                            RHAR1100
      SI=SI*SC+CO*SS                                                    RHAR1110
  116 CO=CO*SC-SIS*SS                                                   RHAR1120
C                                                                       RHAR1130
C     SHIFT C(J)S FOR J=N/2+1 TO J=N UP ONE SLOT                        RHAR1140
      DO 117 I=1,NTOT,2                                                 RHAR1150
      K8=NTOT+4+I                                                       RHAR1160
      A(K8-2)=A(K8)                                                     RHAR1170
  117 A(K8-1)=A(K8+1)                                                   RHAR1180
      DO 500 I=3,NTOT2,2                                                RHAR1190
      A(I) = 2. * A(I)                                                  RHAR1200
  500 A(I + 1) = -2. * A(I + 1)                                         RHAR1210
      RETURN                                                            RHAR1220
      END                                                               RHAR1230
