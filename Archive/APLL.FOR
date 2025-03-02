C                                                                       APLL  10
C     ..................................................................APLL  20
C                                                                       APLL  30
C        SUBROUTINE APLL                                                APLL  40
C                                                                       APLL  50
C        PURPOSE                                                        APLL  60
C           SET UP NORMAL EQUATIONS FOR A LINEAR LEAST SQUARES FIT      APLL  70
C           TO A GIVEN DISCRETE FUNCTION                                APLL  80
C                                                                       APLL  90
C        USAGE                                                          APLL 100
C           CALL APLL(FFCT,N,IP,P,WORK,DATI,IER)                        APLL 110
C           SUBROUTINE FFCT REQUIRES AN EXTERNAL STATEMENT              APLL 120
C                                                                       APLL 130
C        DESCRIPTION OF PARAMETERS                                      APLL 140
C           FFCT  - USER CODED SUBROUTINE WHICH MUST BE DECLARED        APLL 150
C                   EXTERNAL IN THE MAIN PROGRAM. IT IS CALLED          APLL 160
C                   CALL FFCT(I,N,IP,P,DATI,WGT,IER) AND RETURNS        APLL 170
C                   THE VALUES OF THE FUNDAMENTAL FUNCTIONS FOR         APLL 180
C                   THE I-TH ARGUMENT IN P(1) UP TO P(IP)               APLL 190
C                   FOLLOWED BY THE I-TH FUNCTION VALUE IN P(IP+1)      APLL 200
C                   N IS THE NUMBER OF ALL POINTS                       APLL 210
C                   DATI IS A DUMMY PARAMETER WHICH IS USED AS ARRAY    APLL 220
C                   NAME. THE GIVEN DATA SET MAY BE ALLOCATED IN DATI   APLL 230
C                   WGT IS THE WEIGHT FACTOR FOR THE I-TH POINT         APLL 240
C                   IER IS USED AS RESULTANT ERROR PARAMETER IN FFCT    APLL 250
C           N     - NUMBER OF GIVEN POINTS                              APLL 260
C           IP    - NUMBER OF FUNDAMENTAL FUNCTIONS USED FOR LEAST      APLL 270
C                   SQUARES FIT                                         APLL 280
C                   IP SHOULD NOT EXCEED N                              APLL 290
C           P     - WORKING STORAGE OF DIMENSION IP+1, WHICH            APLL 300
C                   IS USED AS INTERFACE BETWEEN APLL AND THE USER      APLL 310
C                   CODED SUBROUTINE FFCT                               APLL 320
C           WORK  - WORKING STORAGE OF DIMENSION (IP+1)*(IP+2)/2.       APLL 330
C                   ON RETURN WORK CONTAINS THE SYMMETRIC COEFFICIENT   APLL 340
C                   MATRIX OF THE NORMAL EQUATIONS IN COMPRESSED FORM,  APLL 350
C                   I.E. UPPER TRINGULAR PART ONLY STORED COLUMNWISE.   APLL 360
C                   THE FOLLOWING IP POSITIONS CONTAIN THE RIGHT        APLL 370
C                   HAND SIDE AND WORK((IP+1)*(IP+2)/2) CONTAINS        APLL 380
C                   THE WEIGHTED SQUARE SUM OF THE FUNCTION VALUES      APLL 390
C           DATI  - DUMMY ENTRY TO COMMUNICATE AN ARRAY NAME BETWEEN    APLL 400
C                   MAIN LINE AND SUBROUTINE FFCT.                      APLL 410
C           IER   - RESULTING ERROR PARAMETER                           APLL 420
C                   IER =-1 MEANS FORMAL ERRORS IN SPECIFIED DIMENSIONS APLL 430
C                   IER = 0 MEANS NO ERRORS                             APLL 440
C                   IER = 1 MEANS ERROR IN EXTERNAL SUBROUTINE FFCT     APLL 450
C                                                                       APLL 460
C        REMARKS                                                        APLL 470
C           TO ALLOW FOR EASY COMMUNICATION OF INTEGER VALUES           APLL 480
C           BETWEEN MAINLINE AND EXTERNAL SUBROUTINE FFCT, THE ERROR    APLL 490
C           PARAMETER IER IS TREATED AS A VECTOR OF DIMENSION 1 WITHIN  APLL 500
C           SUBROUTINE APLL. ADDITIONAL COMPONENTS OF IER MAY BE        APLL 510
C           INTRODUCED BY THE USER FOR COMMUNICATION BACK AND FORTH.    APLL 520
C           IN THIS CASE, HOWEVER, THE USER MUST SPECIFY IER AS A       APLL 530
C           VECTOR IN HIS MAINLINE.                                     APLL 540
C           EXECUTION OF SUBROUTINE APLL IS A PREPARATORY STEP FOR      APLL 550
C           CALCULATION OF THE LINEAR LEAST SQUARES FIT.                APLL 560
C           NORMALLY IT IS FOLLOWED BY EXECUTION OF SUBROUTINE APFS     APLL 570
C                                                                       APLL 580
C       SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                   APLL 590
C           THE EXTERNAL SUBROUTINE FFCT MUST BE FURNISHED BY THE USER  APLL 600
C                                                                       APLL 610
C        METHOD                                                         APLL 620
C           HANDLING OF THE GIVEN DATA SET (ARGUMENTS,FUNCTION VALUES   APLL 630
C           AND WEIGHTS) IS COMPLETELY LEFT TO THE USER                 APLL 640
C           ESSENTIALLY HE HAS THREE CHOICES                            APLL 650
C           (1) THE I-TH VALUES OF ARGUMENT, FUNCTION VALUE AND WEIGHT  APLL 660
C               ARE CALCULATED WITHIN SUBROUTINE FFCT FOR GIVEN I.      APLL 670
C           (2) THE I-TH VALUES OF ARGUMENT, FUNCTION VALUE AND WEIGHT  APLL 680
C               ARE DETERMINED BY TABLE LOOK UP. THE STORAGE LOCATIONS  APLL 690
C               REQUIRED ARE ALLOCATED WITHIN THE DUMMY ARRAY DATI      APLL 700
C               (POSSIBLY IN P TOO, IN EXCESS OF THE SPECIFIED IP + 1   APLL 710
C               LOCATIONS).                                             APLL 720
C               ANOTHER POSSIBILITY WOULD BE TO USE COMMON AS INTERFACE APLL 730
C               BETWEEN MAIN LINE AND SUBROUTINE FFCT AND TO ALLOCATE   APLL 740
C               STORAGE FOR THE DATA SET IN COMMON.                     APLL 750
C           (3) THE I-TH VALUES OF ARGUMENT, FUNCTION VALUE AND WEIGHT  APLL 760
C               ARE READ IN FROM AN EXTERNAL DEVICE. THIS MAY BE EASILY APLL 770
C               ACCOMPLISHED SINCE I IS USED STRICTLY INCREASING FROM   APLL 780
C               ONE UP TO N WITHIN APLL                                 APLL 790
C                                                                       APLL 800
C     ..................................................................APLL 810
C                                                                       APLL 820
      SUBROUTINE APLL(FFCT,N,IP,P,WORK,DATI,IER)                        APLL 830
C                                                                       APLL 840
C                                                                       APLL 850
C        DIMENSIONED DUMMY VARIABLES                                    APLL 860
      DIMENSION P(1),WORK(1),DATI(1),IER(1)                             APLL 870
C                                                                       APLL 880
C        CHECK FOR FORMAL ERRORS IN SPECIFIED DIMENSIONS                APLL 890
      IF(N)10,10,1                                                      APLL 900
    1 IF(IP)10,10,2                                                     APLL 910
    2 IF(N-IP)10,3,3                                                    APLL 920
C                                                                       APLL 930
C        SET WORKING STORAGE AND RIGHT HAND SIDE TO ZERO                APLL 940
    3 IPP1=IP+1                                                         APLL 950
      M=IPP1*(IP+2)/2                                                   APLL 960
      IER(1)=0                                                          APLL 970
      DO 4 I=1,M                                                        APLL 980
    4 WORK(I)=0.                                                        APLL 990
C                                                                       APLL1000
C        START GREAT LOOP OVER ALL GIVEN POINTS                         APLL1010
      DO 8 I=1,N                                                        APLL1020
      CALL FFCT(I,N,IP,P,DATI,WGT,IER)                                  APLL1030
      IF(IER(1))9,5,9                                                   APLL1040
    5 J=0                                                               APLL1050
      DO 7 K=1,IPP1                                                     APLL1060
      AUX=P(K)*WGT                                                      APLL1070
      DO 6 L=1,K                                                        APLL1080
      J=J+1                                                             APLL1090
    6 WORK(J)=WORK(J)+P(L)*AUX                                          APLL1100
    7 CONTINUE                                                          APLL1110
    8 CONTINUE                                                          APLL1120
C                                                                       APLL1130
C        NORMAL RETURN                                                  APLL1140
    9 RETURN                                                            APLL1150
C                                                                       APLL1160
C        ERROR RETURN IN CASE OF FORMAL ERRORS                          APLL1170
   10 IER(1)=-1                                                         APLL1180
      RETURN                                                            APLL1190
      END                                                               APLL1200
