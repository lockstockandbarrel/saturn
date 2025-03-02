C                                                                       DQAT  10
C     ..................................................................DQAT  20
C                                                                       DQAT  30
C        SUBROUTINE DQATR                                               DQAT  40
C                                                                       DQAT  50
C                                                                       DQAT  60
C        PURPOSE                                                        DQAT  70
C           TO COMPUTE AN APPROXIMATION FOR INTEGRAL(FCT(X), SUMMED     DQAT  80
C           OVER X FROM XL TO XU).                                      DQAT  90
C                                                                       DQAT 100
C        USAGE                                                          DQAT 110
C           CALL DQATR (XL,XU,EPS,NDIM,FCT,Y,IER,AUX)                   DQAT 120
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT.               DQAT 130
C                                                                       DQAT 140
C        DESCRIPTION OF PARAMETERS                                      DQAT 150
C           XL     - DOUBLE PRECISION LOWER BOUND OF THE INTERVAL.      DQAT 160
C           XU     - DOUBLE PRECISION UPPER BOUND OF THE INTERVAL.      DQAT 170
C           EPS    - SINGLE PRECISION UPPER BOUND OF THE ABSOLUTE ERROR.DQAT 180
C           NDIM   - THE DIMENSION OF THE AUXILIARY STORAGE ARRAY AUX.  DQAT 190
C                    NDIM-1 IS THE MAXIMAL NUMBER OF BISECTIONS OF      DQAT 200
C                    THE INTERVAL (XL,XU).                              DQAT 210
C           FCT    - THE NAME OF THE EXTERNAL DOUBLE PRECISION FUNCTION DQAT 220
C                    SUBPROGRAM USED.                                   DQAT 230
C           Y      - RESULTING DOUBLE PRECISION APPROXIMATION FOR THE   DQAT 240
C                    INTEGRAL VALUE.                                    DQAT 250
C           IER    - A RESULTING ERROR PARAMETER.                       DQAT 260
C           AUX    - AUXILIARY DOUBLE PRECISION STORAGE ARRAY WITH      DQAT 270
C                    DIMENSION NDIM.                                    DQAT 280
C                                                                       DQAT 290
C        REMARKS                                                        DQAT 300
C           ERROR PARAMETER IER IS CODED IN THE FOLLOWING FORM          DQAT 310
C           IER=0  - IT WAS POSSIBLE TO REACH THE REQUIRED ACCURACY.    DQAT 320
C                    NO ERROR.                                          DQAT 330
C           IER=1  - IT IS IMPOSSIBLE TO REACH THE REQUIRED ACCURACY    DQAT 340
C                    BECAUSE OF ROUNDING ERRORS.                        DQAT 350
C           IER=2  - IT WAS IMPOSSIBLE TO CHECK ACCURACY BECAUSE NDIM   DQAT 360
C                    IS LESS THAN 5, OR THE REQUIRED ACCURACY COULD NOT DQAT 370
C                    BE REACHED WITHIN NDIM-1 STEPS. NDIM SHOULD BE     DQAT 380
C                    INCREASED.                                         DQAT 390
C                                                                       DQAT 400
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DQAT 410
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DQAT 420
C           MUST BE CODED BY THE USER. ITS DOUBLE PRECISION ARGUMENT X  DQAT 430
C           SHOULD NOT BE DESTROYED.                                    DQAT 440
C                                                                       DQAT 450
C        METHOD                                                         DQAT 460
C           EVALUATION OF Y IS DONE BY MEANS OF TRAPEZOIDAL RULE IN     DQAT 470
C           CONNECTION WITH ROMBERGS PRINCIPLE. ON RETURN Y CONTAINS    DQAT 480
C           THE BEST POSSIBLE APPROXIMATION OF THE INTEGRAL VALUE AND   DQAT 490
C           VECTOR AUX THE UPWARD DIAGONAL OF ROMBERG SCHEME.           DQAT 500
C           COMPONENTS AUX(I) (I=1,2,...,IEND, WITH IEND LESS THAN OR   DQAT 510
C           EQUAL TO NDIM) BECOME APPROXIMATIONS TO INTEGRAL VALUE WITH DQAT 520
C           DECREASING ACCURACY BY MULTIPLICATION WITH (XU-XL).         DQAT 530
C           FOR REFERENCE, SEE                                          DQAT 540
C           (1) FILIPPI, DAS VERFAHREN VON ROMBERG-STIEFEL-BAUER ALS    DQAT 550
C               SPEZIALFALL DES ALLGEMEINEN PRINZIPS VON RICHARDSON,    DQAT 560
C               MATHEMATIK-TECHNIK-WIRTSCHAFT, VOL.11, ISS.2 (1964),    DQAT 570
C               PP.49-54.                                               DQAT 580
C           (2) BAUER, ALGORITHM 60, CACM, VOL.4, ISS.6 (1961), PP.255. DQAT 590
C                                                                       DQAT 600
C     ..................................................................DQAT 610
C                                                                       DQAT 620
      SUBROUTINE DQATR(XL,XU,EPS,NDIM,FCT,Y,IER,AUX)                    DQAT 630
C                                                                       DQAT 640
C                                                                       DQAT 650
      DIMENSION AUX(1)                                                  DQAT 660
      DOUBLE PRECISION AUX,XL,XU,X,Y,H,HH,HD,P,Q,SM,FCT                 DQAT 670
C                                                                       DQAT 680
C     PREPARATIONS OF ROMBERG-LOOP                                      DQAT 690
      AUX(1)=.5D0*(FCT(XL)+FCT(XU))                                     DQAT 700
      H=XU-XL                                                           DQAT 710
      IF(NDIM-1)8,8,1                                                   DQAT 720
    1 IF(H)2,10,2                                                       DQAT 730
C                                                                       DQAT 740
C     NDIM IS GREATER THAN 1 AND H IS NOT EQUAL TO 0.                   DQAT 750
    2 HH=H                                                              DQAT 760
      E=EPS/DABS(H)                                                     DQAT 770
      DELT2=0.                                                          DQAT 780
      P=1.D0                                                            DQAT 790
      JJ=1                                                              DQAT 800
      DO 7 I=2,NDIM                                                     DQAT 810
      Y=AUX(1)                                                          DQAT 820
      DELT1=DELT2                                                       DQAT 830
      HD=HH                                                             DQAT 840
      HH=.5D0*HH                                                        DQAT 850
      P=.5D0*P                                                          DQAT 860
      X=XL+HH                                                           DQAT 870
      SM=0.D0                                                           DQAT 880
      DO 3 J=1,JJ                                                       DQAT 890
      SM=SM+FCT(X)                                                      DQAT 900
    3 X=X+HD                                                            DQAT 910
      AUX(I)=.5D0*AUX(I-1)+P*SM                                         DQAT 920
C     A NEW APPROXIMATION OF INTEGRAL VALUE IS COMPUTED BY MEANS OF     DQAT 930
C     TRAPEZOIDAL RULE.                                                 DQAT 940
C                                                                       DQAT 950
C     START OF ROMBERGS EXTRAPOLATION METHOD.                           DQAT 960
      Q=1.D0                                                            DQAT 970
      JI=I-1                                                            DQAT 980
      DO 4 J=1,JI                                                       DQAT 990
      II=I-J                                                            DQAT1000
      Q=Q+Q                                                             DQAT1010
      Q=Q+Q                                                             DQAT1020
    4 AUX(II)=AUX(II+1)+(AUX(II+1)-AUX(II))/(Q-1.D0)                    DQAT1030
C     END OF ROMBERG-STEP                                               DQAT1040
C                                                                       DQAT1050
      DELT2=DABS(Y-AUX(1))                                              DQAT1060
      IF(I-5)7,5,5                                                      DQAT1070
    5 IF(DELT2-E)10,10,6                                                DQAT1080
    6 IF(DELT2-DELT1)7,11,11                                            DQAT1090
    7 JJ=JJ+JJ                                                          DQAT1100
    8 IER=2                                                             DQAT1110
    9 Y=H*AUX(1)                                                        DQAT1120
      RETURN                                                            DQAT1130
   10 IER=0                                                             DQAT1140
      GO TO 9                                                           DQAT1150
   11 IER=1                                                             DQAT1160
      Y=H*Y                                                             DQAT1170
      RETURN                                                            DQAT1180
      END                                                               DQAT1190
