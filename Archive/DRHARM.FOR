C                                                                       DRAR  10
C     ..................................................................DRAR  20
C                                                                       DRAR  30
C        SUBROUTINE DRHARM                                              DRAR  40
C                                                                       DRAR  50
C        PURPOSE                                                        DRAR  60
C           FINDS THE FOURIER COEFFICIENTS OF ONE DIMENSIONAL DOUBLE    DRAR  70
C           PRECISION REAL DATA                                         DRAR  80
C                                                                       DRAR  90
C        USAGE                                                          DRAR 100
C           CALL DRHARM(A,M,INV,S,IFERR)                                DRAR 110
C                                                                       DRAR 120
C        DESCRIPTION OF PARAMETERS                                      DRAR 130
C           A     - A DOUBLE PRECISION VECTOR                           DRAR 140
C                   AS INPUT, CONTAINS ONE DIMENSIONAL REAL DATA. A IS  DRAR 150
C                   2*N+4 CORE LOCATIONS, WHERE N = 2**M. 2*N REAL      DRAR 160
C                   NUMBERS ARE PUT INTO THE FIRST 2*N CORE LOCATIONS   DRAR 170
C                   OF A                                                DRAR 180
C                   AS OUTPUT, A CONTAINS THE FOURIER COEFFICIENTS      DRAR 190
C                   A0/2,B0=0,A1,B1,A2,B2,...,AN/2,BN=0 RESPECTIVELY IN DRAR 200
C                   THE FIRST 2N+2 CORE LOCATIONS OF A                  DRAR 210
C           M     - AN INTEGER WHICH DETERMINES THE SIZE OF THE VECTOR  DRAR 220
C                   A. THE SIZE OF A IS 2*(2**M) + 4                    DRAR 230
C           INV   - A VECTOR WORK AREA FOR BIT AND INDEX MANIPULATION OFDRAR 240
C                   DIMENSION ONE EIGHTH THE NUMBER OF REAL INPUT, VIZ.,DRAR 250
C                   (1/8)*2*(2**M)                                      DRAR 260
C           S     - A DOUBLE PRECISION VECTOR WORK AREA FOR SINE TABLES DRAR 270
C                   WITH DIMENSION THE SAME AS INV                      DRAR 280
C           IFERR - A RETURNED VALUE OF 1 MEANS THAT M IS LESS THAN 3 ORDRAR 290
C                   GREATER THAN 20. OTHERWISE IFERR IS SET = 0         DRAR 300
C                                                                       DRAR 310
C        REMARKS                                                        DRAR 320
C           THIS SUBROUTINE GIVES THE FOURIER COEFFICIENTS OF 2*(2**M)  DRAR 330
C           REAL POINTS. SEE SUBROUTINE DHARM FOR THREE DIMENSIONAL,    DRAR 340
C           DOUBLE PRECISION, COMPLEX FOURIER TRANSFORMS.               DRAR 350
C                                                                       DRAR 360
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DRAR 370
C           DHARM                                                       DRAR 380
C                                                                       DRAR 390
C        METHOD                                                         DRAR 400
C           THE FOURIER COEFFICIENTS A0,B0=0,A1,B1,...,AN,BN=0 ARE      DRAR 410
C           OBTAINED FOR INPUT XJ, J=0,1,2,...,2N-1 FOR THE FOLLOWING   DRAR 420
C           EQUATION (PI = 3.14159...)                                  DRAR 430
C                                                                       DRAR 440
C                 N-1                                               J   DRAR 450
C     XJ=(1/2)A0+SUM (AK*COS(PI*J*K/N)+BK*SIN(PI*J*K/N))+(1/2)AN(-1)    DRAR 460
C                 K=1                                                   DRAR 470
C                                                                       DRAR 480
C           SEE REFERENCE UNDER SUBROUTINE DHARM                        DRAR 490
C                                                                       DRAR 500
C     ..................................................................DRAR 510
C                                                                       DRAR 520
      SUBROUTINE DRHARM(A,M,INV,S,IFERR)                                DRAR 530
      DIMENSION A(1),L(3),INV(1),S(1)                                   DRAR 540
      DOUBLE PRECISION A,SI,AP1IM,FN,CO,CIRE,AP2IM,S,SS,DEL,CIIM,AP1RE, DRAR 550
     1 CNIRE,SC,SIS,AP2RE,CNIIM                                         DRAR 560
      IFSET=1                                                           DRAR 570
      L(1)=M                                                            DRAR 580
      L(2)=0                                                            DRAR 590
      L(3)=0                                                            DRAR 600
      NTOT=2**M                                                         DRAR 610
      NTOT2 = 2*NTOT                                                    DRAR 620
      FN = NTOT                                                         DRAR 630
      DO   3 I = 2,NTOT2,2                                              DRAR 640
   3  A(I) = -A(I)                                                      DRAR 650
      DO   6 I = 1,NTOT2                                                DRAR 660
   6  A(I) = A(I)/FN                                                    DRAR 670
      CALL DHARM(A,L,INV,S,IFSET,IFERR)                                 DRAR 680
C                                                                       DRAR 690
C     MOVE LAST HALF OF A(J)S DOWN ONE SLOT AND ADD A(N) AT BOTTOM TO   DRAR 700
C     GIVE ARRAY FOR A1PRIME AND A2PRIME CALCULATION                    DRAR 710
C                                                                       DRAR 720
   21 DO  52 I=1,NTOT,2                                                 DRAR 730
      J0=NTOT2+2-I                                                      DRAR 740
      A(J0)=A(J0-2)                                                     DRAR 750
   52 A(J0+1)=A(J0-1)                                                   DRAR 760
      A(NTOT2+3)=A(1)                                                   DRAR 770
      A(NTOT2+4)=A(2)                                                   DRAR 780
C                                                                       DRAR 790
C     CALCULATE A1PRIMES AND STORE IN FIRST N SLOTS                     DRAR 800
C     CALCULATE A2PRIMES AND STORE IN SECOND N SLOTS IN REVERSE ORDER   DRAR 810
      K0=NTOT+1                                                         DRAR 820
      DO 104 I=1,K0,2                                                   DRAR 830
      K1=NTOT2-I+4                                                      DRAR 840
      AP1RE=.5*(A(I)+A(K1))                                             DRAR 850
      AP2RE=-.5*(A(I+1)+A(K1+1))                                        DRAR 860
      AP1IM=.5*(-A(I+1)+A(K1+1))                                        DRAR 870
      AP2IM=-.5*(A(I)-A(K1))                                            DRAR 880
      A(I)=AP1RE                                                        DRAR 890
      A(I+1)=AP1IM                                                      DRAR 900
      A(K1)=AP2RE                                                       DRAR 910
  104 A(K1+1)=AP2IM                                                     DRAR 920
      NTO = NTOT/2                                                      DRAR 930
  110 NT=NTO+1                                                          DRAR 940
      DEL=3.141592653589793/DFLOAT(NTOT)                                DRAR 950
      SS=DSIN(DEL)                                                      DRAR 960
      SC=DCOS(DEL)                                                      DRAR 970
      SI=0.0                                                            DRAR 980
      CO=1.0                                                            DRAR 990
C                                                                       DRAR1000
C     COMPUTE C(J)S FOR J=0 THRU J=N                                    DRAR1010
  114 DO 116 I=1,NT                                                     DRAR1020
      K6=NTOT2-2*I+5                                                    DRAR1030
      AP2RE=A(K6)*CO+A(K6+1)*SI                                         DRAR1040
      AP2IM=-A(K6)*SI+A(K6+1)*CO                                        DRAR1050
      CIRE=.5*(A(2*I-1)+AP2RE)                                          DRAR1060
      CIIM=.5*(A(2*I)+AP2IM)                                            DRAR1070
      CNIRE=.5*(A(2*I-1)-AP2RE)                                         DRAR1080
      CNIIM=.5*(A(2*I)-AP2IM)                                           DRAR1090
      A(2*I-1)=CIRE                                                     DRAR1100
      A(2*I)=CIIM                                                       DRAR1110
      A(K6)=CNIRE                                                       DRAR1120
      A(K6+1)=-CNIIM                                                    DRAR1130
      SIS=SI                                                            DRAR1140
      SI=SI*SC+CO*SS                                                    DRAR1150
  116 CO=CO*SC-SIS*SS                                                   DRAR1160
C                                                                       DRAR1170
C     SHIFT C(J)S FOR J=N/2+1 TO J=N UP ONE SLOT                        DRAR1180
      DO 117 I=1,NTOT,2                                                 DRAR1190
      K8=NTOT+4+I                                                       DRAR1200
      A(K8-2)=A(K8)                                                     DRAR1210
  117 A(K8-1)=A(K8+1)                                                   DRAR1220
      DO 500 I=3,NTOT2,2                                                DRAR1230
      A(I) = 2. * A(I)                                                  DRAR1240
  500 A(I + 1) = -2. * A(I + 1)                                         DRAR1250
      RETURN                                                            DRAR1260
      END                                                               DRAR1270
