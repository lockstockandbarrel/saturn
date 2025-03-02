C                                                                       QATR  10
C     ..................................................................QATR  20
C                                                                       QATR  30
C        SUBROUTINE QATR                                                QATR  40
C                                                                       QATR  50
C        PURPOSE                                                        QATR  60
C           TO COMPUTE AN APPROXIMATION FOR INTEGRAL(FCT(X), SUMMED     QATR  70
C           OVER X FROM XL TO XU).                                      QATR  80
C                                                                       QATR  90
C        USAGE                                                          QATR 100
C           CALL QATR (XL,XU,EPS,NDIM,FCT,Y,IER,AUX)                    QATR 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT.               QATR 120
C                                                                       QATR 130
C        DESCRIPTION OF PARAMETERS                                      QATR 140
C           XL     - THE LOWER BOUND OF THE INTERVAL.                   QATR 150
C           XU     - THE UPPER BOUND OF THE INTERVAL.                   QATR 160
C           EPS    - THE UPPER BOUND OF THE ABSOLUTE ERROR.             QATR 170
C           NDIM   - THE DIMENSION OF THE AUXILIARY STORAGE ARRAY AUX.  QATR 180
C                    NDIM-1 IS THE MAXIMAL NUMBER OF BISECTIONS OF      QATR 190
C                    THE INTERVAL (XL,XU).                              QATR 200
C           FCT    - THE NAME OF THE EXTERNAL FUNCTION SUBPROGRAM USED. QATR 210
C           Y      - THE RESULTING APPROXIMATION FOR THE INTEGRAL VALUE.QATR 220
C           IER    - A RESULTING ERROR PARAMETER.                       QATR 230
C           AUX    - AN AUXILIARY STORAGE ARRAY WITH DIMENSION NDIM.    QATR 240
C                                                                       QATR 250
C        REMARKS                                                        QATR 260
C           ERROR PARAMETER IER IS CODED IN THE FOLLOWING FORM          QATR 270
C           IER=0  - IT WAS POSSIBLE TO REACH THE REQUIRED ACCURACY.    QATR 280
C                    NO ERROR.                                          QATR 290
C           IER=1  - IT IS IMPOSSIBLE TO REACH THE REQUIRED ACCURACY    QATR 300
C                    BECAUSE OF ROUNDING ERRORS.                        QATR 310
C           IER=2  - IT WAS IMPOSSIBLE TO CHECK ACCURACY BECAUSE NDIM   QATR 320
C                    IS LESS THAN 5, OR THE REQUIRED ACCURACY COULD NOT QATR 330
C                    BE REACHED WITHIN NDIM-1 STEPS. NDIM SHOULD BE     QATR 340
C                    INCREASED.                                         QATR 350
C                                                                       QATR 360
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QATR 370
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE CODED BY    QATR 380
C           THE USER. ITS ARGUMENT X SHOULD NOT BE DESTROYED.           QATR 390
C                                                                       QATR 400
C        METHOD                                                         QATR 410
C           EVALUATION OF Y IS DONE BY MEANS OF TRAPEZOIDAL RULE IN     QATR 420
C           CONNECTION WITH ROMBERGS PRINCIPLE. ON RETURN Y CONTAINS    QATR 430
C           THE BEST POSSIBLE APPROXIMATION OF THE INTEGRAL VALUE AND   QATR 440
C           VECTOR AUX THE UPWARD DIAGONAL OF ROMBERG SCHEME.           QATR 450
C           COMPONENTS AUX(I) (I=1,2,...,IEND, WITH IEND LESS THAN OR   QATR 460
C           EQUAL TO NDIM) BECOME APPROXIMATIONS TO INTEGRAL VALUE WITH QATR 470
C           DECREASING ACCURACY BY MULTIPLICATION WITH (XU-XL).         QATR 480
C           FOR REFERENCE, SEE                                          QATR 490
C           (1) FILIPPI, DAS VERFAHREN VON ROMBERG-STIEFEL-BAUER ALS    QATR 500
C               SPEZIALFALL DES ALLGEMEINEN PRINZIPS VON RICHARDSON,    QATR 510
C               MATHEMATIK-TECHNIK-WIRTSCHAFT, VOL.11, ISS.2 (1964),    QATR 520
C               PP.49-54.                                               QATR 530
C           (2) BAUER, ALGORITHM 60, CACM, VOL.4, ISS.6 (1961), PP.255. QATR 540
C                                                                       QATR 550
C     ..................................................................QATR 560
C                                                                       QATR 570
      SUBROUTINE QATR(XL,XU,EPS,NDIM,FCT,Y,IER,AUX)                     QATR 580
C                                                                       QATR 590
C                                                                       QATR 600
      DIMENSION AUX(1)                                                  QATR 610
C                                                                       QATR 620
C     PREPARATIONS OF ROMBERG-LOOP                                      QATR 630
      AUX(1)=.5*(FCT(XL)+FCT(XU))                                       QATR 640
      H=XU-XL                                                           QATR 650
      IF(NDIM-1)8,8,1                                                   QATR 660
    1 IF(H)2,10,2                                                       QATR 670
C                                                                       QATR 680
C     NDIM IS GREATER THAN 1 AND H IS NOT EQUAL TO 0.                   QATR 690
    2 HH=H                                                              QATR 700
      E=EPS/ABS(H)                                                      QATR 710
      DELT2=0.                                                          QATR 720
      P=1.                                                              QATR 730
      JJ=1                                                              QATR 740
      DO 7 I=2,NDIM                                                     QATR 750
      Y=AUX(1)                                                          QATR 760
      DELT1=DELT2                                                       QATR 770
      HD=HH                                                             QATR 780
      HH=.5*HH                                                          QATR 790
      P=.5*P                                                            QATR 800
      X=XL+HH                                                           QATR 810
      SM=0.                                                             QATR 820
      DO 3 J=1,JJ                                                       QATR 830
      SM=SM+FCT(X)                                                      QATR 840
    3 X=X+HD                                                            QATR 850
      AUX(I)=.5*AUX(I-1)+P*SM                                           QATR 860
C     A NEW APPROXIMATION OF INTEGRAL VALUE IS COMPUTED BY MEANS OF     QATR 870
C     TRAPEZOIDAL RULE.                                                 QATR 880
C                                                                       QATR 890
C     START OF ROMBERGS EXTRAPOLATION METHOD.                           QATR 900
      Q=1.                                                              QATR 910
      JI=I-1                                                            QATR 920
      DO 4 J=1,JI                                                       QATR 930
      II=I-J                                                            QATR 940
      Q=Q+Q                                                             QATR 950
      Q=Q+Q                                                             QATR 960
    4 AUX(II)=AUX(II+1)+(AUX(II+1)-AUX(II))/(Q-1.)                      QATR 970
C     END OF ROMBERG-STEP                                               QATR 980
C                                                                       QATR 990
      DELT2=ABS(Y-AUX(1))                                               QATR1000
      IF(I-5)7,5,5                                                      QATR1010
    5 IF(DELT2-E)10,10,6                                                QATR1020
    6 IF(DELT2-DELT1)7,11,11                                            QATR1030
    7 JJ=JJ+JJ                                                          QATR1040
    8 IER=2                                                             QATR1050
    9 Y=H*AUX(1)                                                        QATR1060
      RETURN                                                            QATR1070
   10 IER=0                                                             QATR1080
      GO TO 9                                                           QATR1090
   11 IER=1                                                             QATR1100
      Y=H*Y                                                             QATR1110
      RETURN                                                            QATR1120
      END                                                               QATR1130
