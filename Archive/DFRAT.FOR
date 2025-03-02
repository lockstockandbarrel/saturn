C                                                                       DFRA  10
C     ..................................................................DFRA  20
C                                                                       DFRA  30
C        SUBROUTINE DFRAT                                               DFRA  40
C                                                                       DFRA  50
C        PURPOSE                                                        DFRA  60
C           DFRAT IS USED FOR HANDLING OF DATA AND FUNDAMENTAL FUNCTIONSDFRA  70
C           WITH RATIONAL APPROXIMATION. IT IS A SUBSTANTIAL PART OF    DFRA  80
C           RATIONAL APPROXIMATION AND HAS NO MEANING INDEPENDENTLY     DFRA  90
C                                                                       DFRA 100
C        USAGE                                                          DFRA 110
C           CALL DFRAT(I,N,M,P,DATI,WGT,IER)                            DFRA 120
C                                                                       DFRA 130
C        DESCRIPTION OF PARAMETERS                                      DFRA 140
C           I     - SUBSCRIPT OF CURRENT DATA POINT                     DFRA 150
C           N     - NUMBER OF ALL DATA POINTS                           DFRA 160
C           M     - NUMBER OF FUNDAMENTAL FUNCTIONS USED                DFRA 170
C           P     - ARRAY OF DIMENSION M+1 AT LEAST, WHICH CONTAINS     DFRA 180
C                   ON RETURN THE VALUES OF THE M FUNDAMENTAL           DFRA 190
C                   FUNCTIONS, FOLLOWED BY CURRENT FUNCTION VALUE       DFRA 200
C                   P MUST BE OF DOUBLE PRECISION                       DFRA 210
C           DATI  - ARRAY CONTAINING GIVEN N ARGUMENTS, FOLLOWED        DFRA 220
C                   BY N FUNCTION VALUES AND FINALLY BY 1 RESPECTIVELY  DFRA 230
C                   N WEIGHT VALUES                                     DFRA 240
C                   DATI MUST BE OF DOUBLE PRECISION                    DFRA 250
C           WGT   - RESULTANT WEIGHT FACTOR USED FOR I-TH TERM          DFRA 260
C                   WGT MUST BE OF DOUBLE PRECISION                     DFRA 270
C           IER   - RESULTANT ERROR PARAMETER, COMBINED WITH INPUT      DFRA 280
C                   VALUES FOR CONTROL                                  DFRA 290
C                   IER(2) MEANS DIMENSION OF NUMERATOR                 DFRA 300
C                   IER(3) MEANS DIMENSION OF DENOMINATOR               DFRA 310
C                   IER(1) IS USED AS RESULTANT ERROR PARAMETER,        DFRA 320
C                   IER(1) = 0 IN CASE OF NO ERRORS                     DFRA 330
C                   IER(1) = 1 OTHERWISE (ZERO DENOMINATOR)             DFRA 340
C                                                                       DFRA 350
C        REMARKS                                                        DFRA 360
C           VECTOR IER IS USED FOR COMMUNICATION BETWEEN DARAT AND DFRATDFRA 370
C                                                                       DFRA 380
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DFRA 390
C           DCNP                                                        DFRA 400
C                                                                       DFRA 410
C        METHOD                                                         DFRA 420
C           CF. MATHEMATICAL DESCRIPTION OF SUBROUTINE ARAT             DFRA 430
C                                                                       DFRA 440
C     ..................................................................DFRA 450
C                                                                       DFRA 460
      SUBROUTINE DFRAT(I,N,M,P,DATI,WGT,IER)                            DFRA 470
C                                                                       DFRA 480
C                                                                       DFRA 490
C        DIMENSIONED DUMMY VARIABLES                                    DFRA 500
      DIMENSION P(1),DATI(1),IER(1)                                     DFRA 510
      DOUBLE PRECISION P,DATI,WGT,T,F,FNUM,FDEN                         DFRA 520
C                                                                       DFRA 530
C        INITIALIZATION                                                 DFRA 540
      IP=IER(2)                                                         DFRA 550
      IQ=IER(3)                                                         DFRA 560
      IQM1=IQ-1                                                         DFRA 570
      IPQ=IP+IQ                                                         DFRA 580
C                                                                       DFRA 590
C        LOOK UP ARGUMENT, FUNCTION VALUE AND WEIGHT                    DFRA 600
C        LOOK UP NUMERATOR AND DENOMINATOR                              DFRA 610
      T=DATI(I)                                                         DFRA 620
      J=I+N                                                             DFRA 630
      F=DATI(J)                                                         DFRA 640
      FNUM=P(J)                                                         DFRA 650
      J=J+N                                                             DFRA 660
      WGT=1.D0                                                          DFRA 670
      IF(DATI(2*N+1))2,2,1                                              DFRA 680
    1 WGT=DATI(J)                                                       DFRA 690
    2 FDEN=P(J)                                                         DFRA 700
C                                                                       DFRA 710
C        CALCULATE FUNCTION VALUE USED                                  DFRA 720
      F=F*FDEN-FNUM                                                     DFRA 730
C                                                                       DFRA 740
C        CHECK FOR ZERO DENOMINATOR                                     DFRA 750
      IF(FDEN)4,3,4                                                     DFRA 760
C                                                                       DFRA 770
C        ERROR RETURN IN CASE OF ZERO DENOMINATOR                       DFRA 780
    3 IER(1)=1                                                          DFRA 790
      RETURN                                                            DFRA 800
C                                                                       DFRA 810
C        CALCULATE WEIGHT FACTORS USED                                  DFRA 820
    4 WGT=WGT/(FDEN*FDEN)                                               DFRA 830
      FNUM=-FNUM/FDEN                                                   DFRA 840
C                                                                       DFRA 850
C        CALCULATE FUNDAMENTAL FUNCTIONS                                DFRA 860
      J=IQM1                                                            DFRA 870
      IF(IP-IQ)6,6,5                                                    DFRA 880
    5 J=IP-1                                                            DFRA 890
    6 CALL DCNP(P(IQ),T,J)                                              DFRA 900
C                                                                       DFRA 910
C        STORE VALUES OF DENOMINATOR FUNDAMENTAL FUNCTIONS              DFRA 920
    7 IF(IQM1)10,10,8                                                   DFRA 930
    8 DO 9 II=1,IQM1                                                    DFRA 940
      J=II+IQ                                                           DFRA 950
    9 P(II)=P(J)*FNUM                                                   DFRA 960
C                                                                       DFRA 970
C        STORE FUNCTION VALUE                                           DFRA 980
   10 P(IPQ)=F                                                          DFRA 990
C                                                                       DFRA1000
C        NORMAL RETURN                                                  DFRA1010
      IER(1)=0                                                          DFRA1020
      RETURN                                                            DFRA1030
      END                                                               DFRA1040
