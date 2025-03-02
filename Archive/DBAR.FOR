C                                                                       DBAR  10
C     ..................................................................DBAR  20
C                                                                       DBAR  30
C     SUBROUTINE DBAR                                                   DBAR  40
C                                                                       DBAR  50
C     PURPOSE                                                           DBAR  60
C        TO COMPUTE, AT A GIVEN POINT X, AN APPROXIMATION Z TO THE      DBAR  70
C        DERIVATIVE OF AN ANALYTICALLY GIVEN FUNCTION FCT THAT IS 11-   DBAR  80
C        TIMES DIFFERENTIABLE IN A DOMAIN CONTAINING A CLOSED INTERVAL -DBAR  90
C        THE SET OF T BETWEEN X AND X+H (H POSITIVE OR NEGATIVE) - USINGDBAR 100
C        FUNCTION VALUES ONLY ON THAT INTERVAL.                         DBAR 110
C                                                                       DBAR 120
C      USAGE                                                            DBAR 130
C        CALL DBAR(X,H,IH,FCT,Z)                                        DBAR 140
C        PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                   DBAR 150
C                                                                       DBAR 160
C     DESCRIPTION OF PARAMETERS                                         DBAR 170
C        X   - THE POINT AT WHICH THE DERIVATIVE IS TO BE COMPUTED      DBAR 180
C        H   - THE NUMBER THAT DEFINES THE CLOSED INTERVAL WHOSE END-   DBAR 190
C              POINTS ARE X AND X+H (SEE PURPOSE)                       DBAR 200
C        IH  - INPUT PARAMETER (SEE REMARKS AND METHOD)                 DBAR 210
C              IH NON-ZERO - THE SUBROUTINE GENERATES THE INTERNAL      DBAR 220
C                            VALUE HH                                   DBAR 230
C              IH    =   0 - THE INTERNAL VALUE HH IS SET TO H          DBAR 240
C        FCT - THE NAME OF THE EXTERNAL FUNCTION SUBPROGRAM THAT WILL   DBAR 250
C              GENERATE THE NECESSARY FUNCTION VALUES                   DBAR 260
C        Z   - RESULTING DERIVATIVE VALUE                               DBAR 270
C                                                                       DBAR 280
C     REMARKS                                                           DBAR 290
C        (1)  IF H = 0, THEN THERE IS NO COMPUTATION.                   DBAR 300
C        (2)  THE (MAGNITUDE OF THE) INTERNAL VALUE HH, WHICH IS DETER- DBAR 310
C             MINED ACCORDING TO IH, IS THE MAXIMUM STEP-SIZE USED IN   DBAR 320
C             THE COMPUTATION OF THE ONE-SIDED DIVIDED DIFFERENCES (SEE DBAR 330
C             METHOD.)  IF IH IS NON-ZERO, THEN THE SUBROUTINE GENERATESDBAR 340
C             HH ACCORDING TO CRITERIA THAT BALANCE ROUND-OFF AND TRUN- DBAR 350
C             CATION ERROR.  HH ALWAYS HAS THE SAME SIGN AS H AND IT IS DBAR 360
C             ALWAYS LESS THAN OR EQUAL TO THE MAGNITUDE OF H IN AB-    DBAR 370
C             SOLUTE VALUE, SO THAT ALL COMPUTATION OCCURS IN THE CLOSEDDBAR 380
C             INTERVAL DETERMINED BY H.                                 DBAR 390
C                                                                       DBAR 400
C     SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                     DBAR 410
C        THE EXTERNAL FUNCTION SUBPROGRAM FCT(T) MUST BE FURNISHED BY   DBAR 420
C        THE USER.                                                      DBAR 430
C                                                                       DBAR 440
C     METHOD                                                            DBAR 450
C        THE COMPUTATION OF Z IS BASED ON RICHARDSON'S AND ROMBERG'S    DBAR 460
C        EXTRAPOLATION METHOD AS APPLIED TO THE SEQUENCE OF ONE-SIDED   DBAR 470
C        DIVIDED DIFFERENCES ASSOCIATED WITH THE POINT PAIRS            DBAR 480
C        (X,X+(K*HH)/10)K=1,...,10.  (SEE FILLIPI, S. AND ENGELS, H.,   DBAR 490
C        ALTES UND NEUES ZUR NUMERISCHEN DIFFERENTIATION, ELECTRONISCHE DBAR 500
C        DATENVERARBEITUNG, ISS. 2 (1966), PP. 57-65.)                  DBAR 510
C                                                                       DBAR 520
C     ..................................................................DBAR 530
C                                                                       DBAR 540
      SUBROUTINE DBAR(X,H,IH,FCT,Z)                                     DBAR 550
C                                                                       DBAR 560
C                                                                       DBAR 570
      DIMENSION AUX(10)                                                 DBAR 580
C                                                                       DBAR 590
C        NO ACTION IN CASE OF ZERO INTERVAL LENGTH                      DBAR 600
      IF(H)1,17,1                                                       DBAR 610
C                                                                       DBAR 620
C        GENERATE STEPSIZE HH FOR DIVIDED DIFFERENCES                   DBAR 630
    1 C=ABS(H)                                                          DBAR 640
      B=H                                                               DBAR 650
      D=X                                                               DBAR 660
      D=FCT(D)                                                          DBAR 670
      IF(IH)2,9,2                                                       DBAR 680
    2 HH=.5                                                             DBAR 690
      IF(C-HH)3,4,4                                                     DBAR 700
    3 HH=B                                                              DBAR 710
    4 HH=SIGN(HH,B)                                                     DBAR 720
      Z=ABS((FCT(X+HH)-D)/HH)                                           DBAR 730
      A=ABS(D)                                                          DBAR 740
      HH=1.                                                             DBAR 750
      IF(A-1.)6,6,5                                                     DBAR 760
    5 HH=HH*A                                                           DBAR 770
    6 IF(Z-1.)8,8,7                                                     DBAR 780
    7 HH=HH/Z                                                           DBAR 790
    8 IF(HH-C)10,10,9                                                   DBAR 800
    9 HH=B                                                              DBAR 810
   10 HH=SIGN(HH,B)                                                     DBAR 820
C                                                                       DBAR 830
C        INITIALIZE DIFFERENTIATION LOOP                                DBAR 840
      Z=(FCT(X+HH)-D)/HH                                                DBAR 850
      J=10                                                              DBAR 860
      JJ=J-1                                                            DBAR 870
      AUX(J)=Z                                                          DBAR 880
      DH=HH/FLOAT(J)                                                    DBAR 890
      DZ=1.E75                                                          DBAR 900
C                                                                       DBAR 910
C        START DIFFERENTIATION LOOP                                     DBAR 920
   11 J=J-1                                                             DBAR 930
      C=J                                                               DBAR 940
      HH=C*DH                                                           DBAR 950
      AUX(J)=(FCT(X+HH)-D)/HH                                           DBAR 960
C                                                                       DBAR 970
C        INITIALIZE EXTRAPOLATION LOOP                                  DBAR 980
      D2=1.E75                                                          DBAR 990
      B=0.                                                              DBAR1000
      A=1./C                                                            DBAR1010
C                                                                       DBAR1020
C        START EXTRAPOLATION LOOP                                       DBAR1030
      DO 12 I=J,JJ                                                      DBAR1040
      D1=D2                                                             DBAR1050
      B=B+A                                                             DBAR1060
      HH=(AUX(I)-AUX(I+1))/B                                            DBAR1070
      AUX(I+1)=AUX(I)+HH                                                DBAR1080
C                                                                       DBAR1090
C        TEST ON OSCILLATING INCREMENTS                                 DBAR1100
      D2=ABS(HH)                                                        DBAR1110
      IF(D2-D1)12,13,13                                                 DBAR1120
   12 CONTINUE                                                          DBAR1130
C        END OF EXTRAPOLATION LOOP                                      DBAR1140
C                                                                       DBAR1150
C        UPDATE RESULT VALUE Z                                          DBAR1160
      I=JJ+1                                                            DBAR1170
      GO TO 14                                                          DBAR1180
   13 D2=D1                                                             DBAR1190
      JJ=I                                                              DBAR1200
   14 IF(D2-DZ)15,16,16                                                 DBAR1210
   15 DZ=D2                                                             DBAR1220
      Z=AUX(I)                                                          DBAR1230
   16 IF(J-1)17,17,11                                                   DBAR1240
C        END OF DIFFERENTIATION LOOP                                    DBAR1250
C                                                                       DBAR1260
   17 RETURN                                                            DBAR1270
      END                                                               DBAR1280
