C                                                                       DAPC  10
C     ..................................................................DAPC  20
C                                                                       DAPC  30
C        SUBROUTINE DAPCH                                               DAPC  40
C                                                                       DAPC  50
C        PURPOSE                                                        DAPC  60
C           SET UP NORMAL EQUATIONS OF LEAST SQUARES FIT IN TERMS OF    DAPC  70
C           CHEBYSHEV POLYNOMIALS FOR A GIVEN DISCRETE FUNCTION         DAPC  80
C                                                                       DAPC  90
C        USAGE                                                          DAPC 100
C           CALL DAPCH(DATI,N,IP,XD,X0,WORK,IER)                        DAPC 110
C                                                                       DAPC 120
C        DESCRIPTION OF PARAMETERS                                      DAPC 130
C           DATI  - VECTOR OF DIMENSION 3*N (OR DIMENSION 2*N+1)        DAPC 140
C                   CONTAINING THE GIVEN ARGUMENTS, FOLLOWED BY THE     DAPC 150
C                   FUNCTION VALUES AND N (RESPECTIVELY 1) WEIGHT       DAPC 160
C                   VALUES. THE CONTENT OF VECTOR DATI REMAINS          DAPC 170
C                   UNCHANGED.                                          DAPC 180
C           N     - NUMBER OF GIVEN POINTS                              DAPC 190
C           IP    - DIMENSION OF LEAST SQUARES FIT, I.E. NUMBER OF      DAPC 200
C                   CHEBYSHEV POLYNOMIALS USED AS FUNDAMENTAL FUNCTIONS DAPC 210
C                   IP SHOULD NOT EXCEED N                              DAPC 220
C           XD    - RESULTANT MULTIPLICATIVE CONSTANT FOR LINEAR        DAPC 230
C                   TRANSFORMATION OF ARGUMENT RANGE                    DAPC 240
C           X0    - RESULTANT ADDITIVE CONSTANT FOR LINEAR              DAPC 250
C                   TRANSFORMATION OF ARGUMENT RANGE                    DAPC 260
C           WORK  - WORKING STORAGE OF DIMENSION (IP+1)*(IP+2)/2        DAPC 270
C                   ON RETURN WORK CONTAINS THE SYMMETRIC COEFFICIENT   DAPC 280
C                   MATRIX OF THE NORMAL EQUATIONS IN COMPRESSED FORM   DAPC 290
C                   FOLLOWED IMMEDIATELY BY RIGHT HAND SIDE             DAPC 300
C                   AND SQUARE SUM OF FUNCTION VALUES                   DAPC 310
C           IER   - RESULTING ERROR PARAMETER                           DAPC 320
C                   IER =-1 MEANS FORMAL ERRORS IN DIMENSION            DAPC 330
C                   IER = 0 MEANS NO ERRORS                             DAPC 340
C                   IER = 1 MEANS COINCIDING ARGUMENTS                  DAPC 350
C        NOTE: DATI,WORK,XD,X0,XA,XE,XM,DF,T,SUM ARE DOUBLE PRECISION   DAPC 360
C        REMARKS                                                        DAPC 370
C           NO WEIGHTS ARE USED IF THE VALUE OF DATI(2*N+1) IS          DAPC 380
C           NOT POSITIVE.                                               DAPC 390
C           EXECUTION OF SUBROUTINE APCH IS A PREPARATORY STEP FOR      DAPC 400
C           CALCULATION OF LEAST SQUARES FITS IN CHEBYSHEV POLYNOMIALS  DAPC 410
C           IT SHOULD BE FOLLOWED BY EXECUTION OF SUBROUTINE APFS       DAPC 420
C                                                                       DAPC 430
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DAPC 440
C           NONE                                                        DAPC 450
C                                                                       DAPC 460
C        METHOD                                                         DAPC 470
C           THE LEAST SQUARE FIT IS DETERMINED USING CHEBYSHEV          DAPC 480
C           POLYNOMIALS AS FUNDAMENTAL FUNCTION SYSTEM.                 DAPC 490
C           THE METHOD IS DISCUSSED IN THE ARTICLE                      DAPC 500
C           A.T.BERZTISS, LEAST SQUARES FITTING TO IRREGULARLY SPACED   DAPC 510
C           DATA, SIAM REVIEW, VOL.6, ISS.3, 1964, PP. 203-227.         DAPC 520
C                                                                       DAPC 530
C     ..................................................................DAPC 540
C                                                                       DAPC 550
      SUBROUTINE DAPCH(DATI,N,IP,XD,X0,WORK,IER)                        DAPC 560
C                                                                       DAPC 570
        DOUBLE PRECISION DATI,WORK,XD,X0,XA,XE,XM,DF,T,SUM              DAPC 580
C       DIMENSIONED DUMMY VARIABLES                                     DAPC 590
      DIMENSION DATI(1),WORK(1)                                         DAPC 600
C                                                                       DAPC 610
C        CHECK FOR FORMAL ERRORS IN SPECIFIED DIMENSIONS                DAPC 620
      IF(N-1)19,20,1                                                    DAPC 630
    1 IF(IP)19,19,2                                                     DAPC 640
C                                                                       DAPC 650
C        SEARCH SMALLEST AND LARGEST ARGUMENT                           DAPC 660
    2 IF(IP-N)3,3,19                                                    DAPC 670
    3 XA=DATI(1)                                                        DAPC 680
      X0=XA                                                             DAPC 690
      XE=0.D0                                                           DAPC 700
      DO 7 I=1,N                                                        DAPC 710
      XM=DATI(I)                                                        DAPC 720
      IF(XA-XM)5,5,4                                                    DAPC 730
    4 XA=XM                                                             DAPC 740
    5 IF(X0-XM)6,7,7                                                    DAPC 750
    6 X0=XM                                                             DAPC 760
    7 CONTINUE                                                          DAPC 770
C                                                                       DAPC 780
C        INITIALIZE CALCULATION OF NORMAL EQUATIONS                     DAPC 790
      XD=X0-XA                                                          DAPC 800
      M=(IP*(IP+1))/2                                                   DAPC 810
      IEND=M+IP+1                                                       DAPC 820
      MT2=IP+IP                                                         DAPC 830
      MT2M=MT2-1                                                        DAPC 840
C                                                                       DAPC 850
C        SET WORKING STORAGE AND RIGHT HAND SIDE TO ZERO                DAPC 860
      DO 8 I=1,IP                                                       DAPC 870
      J=MT2-I                                                           DAPC 880
      WORK(J)=0.D0                                                      DAPC 890
      WORK(I)=0.D0                                                      DAPC 900
      K=M+I                                                             DAPC 910
    8 WORK(K)=0.D0                                                      DAPC 920
C                                                                       DAPC 930
C        CHECK FOR DEGENERATE ARGUMENT RANGE                            DAPC 940
      IF(XD)20,20,9                                                     DAPC 950
C                                                                       DAPC 960
C        CALCULATE CONSTANTS FOR REDUCTION OF ARGUMENTS                 DAPC 970
    9 X0=-(X0+XA)/XD                                                    DAPC 980
      XD=2.D0/XD                                                        DAPC 990
      SUM=0.D0                                                          DAPC1000
C                                                                       DAPC1010
C        START GREAT LOOP OVER ALL GIVEN POINTS                         DAPC1020
      DO 15 I=1,N                                                       DAPC1030
      T=DATI(I)*XD+X0                                                   DAPC1040
      J=I+N                                                             DAPC1050
      DF=DATI(J)                                                        DAPC1060
C                                                                       DAPC1070
C        CALCULATE AND STORE VALUES OF CHEBYSHEV POLYNOMIALS            DAPC1080
C        FOR ARGUMENT T                                                 DAPC1090
      XA=1.D0                                                           DAPC1100
      XM=T                                                              DAPC1110
      IF(DATI(2*N+1))11,11,10                                           DAPC1120
   10 J=J+N                                                             DAPC1130
      XA=DATI(J)                                                        DAPC1140
      XM=T*XA                                                           DAPC1150
   11 T=T+T                                                             DAPC1160
      SUM=SUM+DF*DF*XA                                                  DAPC1170
      DF=DF+DF                                                          DAPC1180
      J=1                                                               DAPC1190
   12 K=M+J                                                             DAPC1200
      WORK(K)=WORK(K)+DF*XA                                             DAPC1210
   13 WORK(J)=WORK(J)+XA                                                DAPC1220
      IF(J-MT2M)14,15,15                                                DAPC1230
   14 J=J+1                                                             DAPC1240
      XE=T*XM-XA                                                        DAPC1250
      XA=XM                                                             DAPC1260
      XM=XE                                                             DAPC1270
      IF(J-IP)12,12,13                                                  DAPC1280
   15 CONTINUE                                                          DAPC1290
      WORK(IEND)=SUM+SUM                                                DAPC1300
C                                                                       DAPC1310
C        CALCULATE MATRIX OF NORMAL EQUATIONS                           DAPC1320
      LL=M                                                              DAPC1330
      KK=MT2M                                                           DAPC1340
      JJ=1                                                              DAPC1350
      K=KK                                                              DAPC1360
      DO 18 J=1,M                                                       DAPC1370
      WORK(LL)=WORK(K)+WORK(JJ)                                         DAPC1380
      LL=LL-1                                                           DAPC1390
      IF(K-JJ)16,16,17                                                  DAPC1400
   16 KK=KK-2                                                           DAPC1410
      K=KK                                                              DAPC1420
      JJ=1                                                              DAPC1430
      GOTO 18                                                           DAPC1440
   17 JJ=JJ+1                                                           DAPC1450
      K=K-1                                                             DAPC1460
   18 CONTINUE                                                          DAPC1470
      IER=0                                                             DAPC1480
      RETURN                                                            DAPC1490
C                                                                       DAPC1500
C        ERROR RETURN IN CASE OF FORMAL ERRORS                          DAPC1510
   19 IER=-1                                                            DAPC1520
      RETURN                                                            DAPC1530
C                                                                       DAPC1540
C        ERROR RETURN IN CASE OF COINCIDING ARGUMENTS                   DAPC1550
   20 IER=1                                                             DAPC1560
      RETURN                                                            DAPC1570
      END                                                               DAPC1580
