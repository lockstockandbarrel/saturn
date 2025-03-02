C                                                                       DAPL  10
C     ..................................................................DAPL  20
C                                                                       DAPL  30
C        SUBROUTINE DAPLL                                               DAPL  40
C                                                                       DAPL  50
C        PURPOSE                                                        DAPL  60
C           SET UP NORMAL EQUATIONS FOR A LINEAR LEAST SQUARES FIT      DAPL  70
C           TO A GIVEN DISCRETE FUNCTION                                DAPL  80
C                                                                       DAPL  90
C        USAGE                                                          DAPL 100
C           CALL DAPLL(FFCT,N,IP,P,WORK,DATI,IER)                       DAPL 110
C           SUBROUTINE FFCT REQUIRES AN EXTERNAL STATEMENT              DAPL 120
C                                                                       DAPL 130
C        DESCRIPTION OF PARAMETERS                                      DAPL 140
C           FFCT  - USER CODED SUBROUTINE WHICH MUST BE DECLARED        DAPL 150
C                   EXTERNAL IN THE MAIN PROGRAM. IT IS CALLED          DAPL 160
C                   CALL FFCT(I,N,IP,P,DATI,WGT,IER) AND RETURNS        DAPL 170
C                   THE VALUES OF THE FUNDAMENTAL FUNCTIONS FOR         DAPL 180
C                   THE I-TH ARGUMENT IN P(1) UP TO P(IP)               DAPL 190
C                   FOLLOWED BY THE I-TH FUNCTION VALUE IN P(IP+1)      DAPL 200
C                   N IS THE NUMBER OF ALL POINTS                       DAPL 210
C                   DATI IS A DUMMY PARAMETER WHICH IS USED AS ARRAY    DAPL 220
C                   NAME. THE GIVEN DATA SET MAY BE ALLOCATED IN DATI   DAPL 230
C                   WGT IS THE WEIGHT FACTOR FOR THE I-TH POINT         DAPL 240
C                   IER IS USED AS RESULTANT ERROR PARAMETER IN FFCT    DAPL 250
C           N     - NUMBER OF GIVEN POINTS                              DAPL 260
C           IP    - NUMBER OF FUNDAMENTAL FUNCTIONS USED FOR LEAST      DAPL 270
C                   SQUARES FIT                                         DAPL 280
C                   IP SHOULD NOT EXCEED N                              DAPL 290
C           P     - WORKING STORAGE OF DIMENSION IP+1, WHICH            DAPL 300
C                   IS USED AS INTERFACE BETWEEN APLL AND THE USER      DAPL 310
C                   CODED SUBROUTINE FFCT                               DAPL 320
C           WORK  - WORKING STORAGE OF DIMENSION (IP+1)*(IP+2)/2.       DAPL 330
C                   ON RETURN WORK CONTAINS THE SYMMETRIC COEFFICIENT   DAPL 340
C                   MATRIX OF THE NORMAL EQUATIONS IN COMPRESSED FORM,  DAPL 350
C                   I.E. UPPER TRINGULAR PART ONLY STORED COLUMNWISE.   DAPL 360
C                   THE FOLLOWING IP POSITIONS CONTAIN THE RIGHT        DAPL 370
C                   HAND SIDE AND WORK((IP+1)*(IP+2)/2) CONTAINS        DAPL 380
C                   THE WEIGHTED SQUARE SUM OF THE FUNCTION VALUES      DAPL 390
C           DATI  - DUMMY ENTRY TO COMMUNICATE AN ARRAY NAME BETWEEN    DAPL 400
C                   MAIN LINE AND SUBROUTINE FFCT.                      DAPL 410
C           IER   - RESULTING ERROR PARAMETER                           DAPL 420
C                   IER =-1 MEANS FORMAL ERRORS IN SPECIFIED DIMENSIONS DAPL 430
C                   IER = 0 MEANS NO ERRORS                             DAPL 440
C                   IER = 1 MEANS ERROR IN EXTERNAL SUBROUTINE FFCT     DAPL 450
C        NOTE: P,WORK,DATI,WGT AND AUX ARE DOUBLE PRECISION             DAPL 460
C        REMARKS                                                        DAPL 470
C           TO ALLOW FOR EASY COMMUNICATION OF INTEGER VALUES           DAPL 480
C           BETWEEN MAINLINE AND EXTERNAL SUBROUTINE FFCT, THE ERROR    DAPL 490
C           PARAMETER IER IS TREATED AS A VECTOR OF DIMENSION 1 WITHIN  DAPL 500
C           SUBROUTINE APLL. ADDITIONAL COMPONENTS OF IER MAY BE        DAPL 510
C           INTRODUCED BY THE USER FOR COMMUNICATION BACK AND FORTH.    DAPL 520
C           IN THIS CASE, HOWEVER, THE USER MUST SPECIFY IER AS A       DAPL 530
C           VECTOR IN HIS MAINLINE.                                     DAPL 540
C           EXECUTION OF SUBROUTINE APLL IS A PREPARATORY STEP FOR      DAPL 550
C           CALCULATION OF THE LINEAR LEAST SQUARES FIT.                DAPL 560
C           NORMALLY IT IS FOLLOWED BY EXECUTION OF SUBROUTINE APFS     DAPL 570
C                                                                       DAPL 580
C       SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                   DAPL 590
C           THE EXTERNAL SUBROUTINE FFCT MUST BE FURNISHED BY THE USER  DAPL 600
C                                                                       DAPL 610
C        METHOD                                                         DAPL 620
C           HANDLING OF THE GIVEN DATA SET (ARGUMENTS,FUNCTION VALUES   DAPL 630
C           AND WEIGHTS) IS COMPLETELY LEFT TO THE USER                 DAPL 640
C           ESSENTIALLY HE HAS THREE CHOICES                            DAPL 650
C           (1) THE I-TH VALUES OF ARGUMENT, FUNCTION VALUE AND WEIGHT  DAPL 660
C               ARE CALCULATED WITHIN SUBROUTINE FFCT FOR GIVEN I.      DAPL 670
C           (2) THE I-TH VALUES OF ARGUMENT, FUNCTION VALUE AND WEIGHT  DAPL 680
C               ARE DETERMINED BY TABLE LOOK UP. THE STORAGE LOCATIONS  DAPL 690
C               REQUIRED ARE ALLOCATED WITHIN THE DUMMY ARRAY DATI      DAPL 700
C               (POSSIBLY IN P TOO, IN EXCESS OF THE SPECIFIED IP + 1   DAPL 710
C               LOCATIONS).                                             DAPL 720
C               ANOTHER POSSIBILITY WOULD BE TO USE COMMON AS INTERFACE DAPL 730
C               BETWEEN MAIN LINE AND SUBROUTINE FFCT AND TO ALLOCATE   DAPL 740
C               STORAGE FOR THE DATA SET IN COMMON.                     DAPL 750
C           (3) THE I-TH VALUES OF ARGUMENT, FUNCTION VALUE AND WEIGHT  DAPL 760
C               ARE READ IN FROM AN EXTERNAL DEVICE. THIS MAY BE EASILY DAPL 770
C               ACCOMPLISHED SINCE I IS USED STRICTLY INCREASING FROM   DAPL 780
C               ONE UP TO N WITHIN APLL                                 DAPL 790
C                                                                       DAPL 800
C     ..................................................................DAPL 810
C                                                                       DAPL 820
      SUBROUTINE DAPLL(FFCT,N,IP,P,WORK,DATI,IER)                       DAPL 830
C                                                                       DAPL 840
      DOUBLE PRECISION P,WORK,DATI,WGT,AUX                              DAPL 850
C        DIMENSIONED DUMMY VARIABLES                                    DAPL 860
      DIMENSION P(1),WORK(1),DATI(1),IER(1)                             DAPL 870
C                                                                       DAPL 880
C        CHECK FOR FORMAL ERRORS IN SPECIFIED DIMENSIONS                DAPL 890
      IF(N)10,10,1                                                      DAPL 900
    1 IF(IP)10,10,2                                                     DAPL 910
    2 IF(N-IP)10,3,3                                                    DAPL 920
C                                                                       DAPL 930
C        SET WORKING STORAGE AND RIGHT HAND SIDE TO ZERO                DAPL 940
    3 IPP1=IP+1                                                         DAPL 950
      M=IPP1*(IP+2)/2                                                   DAPL 960
      IER(1)=0                                                          DAPL 970
      DO 4 I=1,M                                                        DAPL 980
    4 WORK(I)=0.0D0                                                     DAPL 990
C                                                                       DAPL1000
C        START GREAT LOOP OVER ALL GIVEN POINTS                         DAPL1010
      DO 8 I=1,N                                                        DAPL1020
      CALL FFCT(I,N,IP,P,DATI,WGT,IER)                                  DAPL1030
      IF(IER(1))9,5,9                                                   DAPL1040
    5 J=0                                                               DAPL1050
      DO 7 K=1,IPP1                                                     DAPL1060
      AUX=P(K)*WGT                                                      DAPL1070
      DO 6 L=1,K                                                        DAPL1080
      J=J+1                                                             DAPL1090
    6 WORK(J)=WORK(J)+P(L)*AUX                                          DAPL1100
    7 CONTINUE                                                          DAPL1110
    8 CONTINUE                                                          DAPL1120
C                                                                       DAPL1130
C        NORMAL RETURN                                                  DAPL1140
    9 RETURN                                                            DAPL1150
C                                                                       DAPL1160
C        ERROR RETURN IN CASE OF FORMAL ERRORS                          DAPL1170
   10 IER(1)=-1                                                         DAPL1180
      RETURN                                                            DAPL1190
      END                                                               DAPL1200
