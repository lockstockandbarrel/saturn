C                                                                       DDBA  10
C     ..................................................................DDBA  20
C                                                                       DDBA  30
C     SUBROUTINE DDBAR                                                  DDBA  40
C                                                                       DDBA  50
C     PURPOSE                                                           DDBA  60
C        TO COMPUTE, AT A GIVEN POINT X, AN APPROXIMATION Z TO THE      DDBA  70
C        DERIVATIVE OF AN ANALYTICALLY GIVEN FUNCTION FCT THAT IS 11-   DDBA  80
C        TIMES DIFFERENTIABLE IN A DOMAIN CONTAINING A CLOSED INTERVAL -DDBA  90
C        THE SET OF T BETWEEN X AND X+H (H POSITIVE OR NEGATIVE) - USINGDDBA 100
C        FUNCTION VALUES ONLY ON THAT INTERVAL.                         DDBA 110
C                                                                       DDBA 120
C      USAGE                                                            DDBA 130
C        CALL DDBAR(X,H,IH,FCT,Z,)                                      DDBA 140
C        PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                   DDBA 150
C                                                                       DDBA 160
C     DESCRIPTION OF PARAMETERS                                         DDBA 170
C        X   - THE POINT AT WHICH THE DERIVATIVE IS TO BE COMPUTED      DDBA 180
C              X IS IN DOUBLE PRECISION                                 DDBA 190
C        H   - THE NUMBER THAT DEFINES THE CLOSED INTERVAL WHOSE END-   DDBA 200
C              POINTS ARE X AND X+H (SEE PURPOSE)                       DDBA 210
C              H IS IN SINGLE PRECISION                                 DDBA 220
C        IH  - INPUT PARAMETER (SEE REMARKS AND METHOD)                 DDBA 230
C              IH NON-ZERO - THE SUBROUTINE GENERATES THE INTERNAL      DDBA 240
C                            VALUE HH                                   DDBA 250
C              IH    =   0 - THE INTERNAL VALUE HH IS SET TO H          DDBA 260
C        FCT - THE NAME OF THE EXTERNAL DOUBLE PRECISION FUNCTION       DDBA 270
C              SUBPROGRAM THAT WILL GENERATE THE NECESSARY FUNCTION     DDBA 280
C              VALUES.                                                  DDBA 290
C        Z   - RESULTING DERIVATIVE VALUE - DOUBLE PRECISION            DDBA 300
C                                                                       DDBA 310
C     REMARKS                                                           DDBA 320
C        (1)  IF H = 0, THEN THERE IS NO COMPUTATION.                   DDBA 330
C        (2)  THE (MAGNITUDE OF THE) INTERNAL VALUE HH, WHICH IS DETER- DDBA 340
C             MINED ACCORDING TO IH, IS THE MAXIMUM STEP-SIZE USED IN   DDBA 350
C             THE COMPUTATION OF THE ONE-SIDED DIVIDED DIFFERENCES (SEE DDBA 360
C             METHOD.)  IF IH IS NON-ZERO, THEN THE SUBROUTINE GENERATESDDBA 370
C             HH ACCORDING TO CRITERIA THAT BALANCE ROUND-OFF AND TRUN- DDBA 380
C             CATION ERROR.  HH ALWAYS HAS THE SAME SIGN AS H AND IT IS DDBA 390
C             ALWAYS LESS THAN OR EQUAL TO THE MAGNITUDE OF H IN AB-    DDBA 400
C             SOLUTE VALUE, SO THAT ALL COMPUTATION OCCURS IN THE CLOSEDDDBA 410
C             INTERVAL DETERMINED BY H.                                 DDBA 420
C                                                                       DDBA 430
C     SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                     DDBA 440
C        THE EXTERNAL FUNCTION SUBPROGRAM FCT(T) MUST BE FURNISHED BY   DDBA 450
C        THE USER. FCT(T) IS IN DOUBLE PRECISION                        DDBA 460
C                                                                       DDBA 470
C     METHOD                                                            DDBA 480
C        THE COMPUTATION OF Z IS BASED ON RICHARDSON'S AND ROMBERG'S    DDBA 490
C        EXTRAPOLATION METHOD AS APPLIED TO THE SEQUENCE OF ONE-SIDED   DDBA 500
C        DIVIDED DIFFERENCES ASSOCIATED WITH THE POINT PAIRS            DDBA 510
C        (X,X+(K*HH)/10)K=1,...,10.  (SEE FILLIPI, S. AND ENGELS, H.,   DDBA 520
C        ALTES UND NEUES ZUR NUMERISCHEN DIFFERENTIATION, ELECTRONISCHE DDBA 530
C        DATENVERARBEITUNG, ISS. 2 (1966), PP. 57-65.)                  DDBA 540
C                                                                       DDBA 550
C     ..................................................................DDBA 560
C                                                                       DDBA 570
      SUBROUTINE DDBAR(X,H,IH,FCT,Z)                                    DDBA 580
C                                                                       DDBA 590
C                                                                       DDBA 600
      DIMENSION AUX(10)                                                 DDBA 610
      DOUBLE PRECISION X,FCT,Z,AUX,A,B,C,D,DH,HH                        DDBA 620
C                                                                       DDBA 630
C        NO ACTION IN CASE OF ZERO INTERVAL LENGTH                      DDBA 640
      IF(H)1,17,1                                                       DDBA 650
C                                                                       DDBA 660
C        GENERATE STEPSIZE HH FOR DIVIDED DIFFERENCES                   DDBA 670
    1 C=ABS(H)                                                          DDBA 680
      B=H                                                               DDBA 690
      D=X                                                               DDBA 700
      D=FCT(D)                                                          DDBA 710
      IF(IH)2,9,2                                                       DDBA 720
    2 HH=.5D-2                                                          DDBA 730
      IF(C-HH)3,4,4                                                     DDBA 740
    3 HH=B                                                              DDBA 750
    4 HH=DSIGN(HH,B)                                                    DDBA 760
      Z=DABS((FCT(X+HH)-D)/HH)                                          DDBA 770
      A=DABS(D)                                                         DDBA 780
      HH=1.D-2                                                          DDBA 790
      IF(A-1.D0)6,6,5                                                   DDBA 800
    5 HH=HH*A                                                           DDBA 810
    6 IF(Z-1.D0)8,8,7                                                   DDBA 820
    7 HH=HH/Z                                                           DDBA 830
    8 IF(HH-C)10,10,9                                                   DDBA 840
    9 HH=B                                                              DDBA 850
   10 HH=DSIGN(HH,B)                                                    DDBA 860
C                                                                       DDBA 870
C        INITIALIZE DIFFERENTIATION LOOP                                DDBA 880
      Z=(FCT(X+HH)-D)/HH                                                DDBA 890
      J=10                                                              DDBA 900
      JJ=J-1                                                            DDBA 910
      AUX(J)=Z                                                          DDBA 920
      DH=HH/DFLOAT(J)                                                   DDBA 930
      DZ=1.E75                                                          DDBA 940
C                                                                       DDBA 950
C        START DIFFERENTIATION LOOP                                     DDBA 960
   11 J=J-1                                                             DDBA 970
      C=J                                                               DDBA 980
      HH=C*DH                                                           DDBA 990
      AUX(J)=(FCT(X+HH)-D)/HH                                           DDBA1000
C                                                                       DDBA1010
C        INITIALIZE EXTRAPOLATION LOOP                                  DDBA1020
      D2=1.E75                                                          DDBA1030
      B=0.D0                                                            DDBA1040
      A=1.D0/C                                                          DDBA1050
C                                                                       DDBA1060
C        START EXTRAPOLATION LOOP                                       DDBA1070
      DO 12 I=J,JJ                                                      DDBA1080
      D1=D2                                                             DDBA1090
      B=B+A                                                             DDBA1100
      HH=(AUX(I)-AUX(I+1))/B                                            DDBA1110
      AUX(I+1)=AUX(I)+HH                                                DDBA1120
C                                                                       DDBA1130
C        TEST ON OSCILLATING INCREMENTS                                 DDBA1140
      D2=DABS(HH)                                                       DDBA1150
      IF(D2-D1)12,13,13                                                 DDBA1160
   12 CONTINUE                                                          DDBA1170
C        END OF EXTRAPOLATION LOOP                                      DDBA1180
C                                                                       DDBA1190
C        UPDATE RESULT VALUE Z                                          DDBA1200
      I=JJ+1                                                            DDBA1210
      GO TO 14                                                          DDBA1220
   13 D2=D1                                                             DDBA1230
      JJ=I                                                              DDBA1240
   14 IF(D2-DZ)15,16,16                                                 DDBA1250
   15 DZ=D2                                                             DDBA1260
      Z=AUX(I)                                                          DDBA1270
   16 IF(J-1)17,17,11                                                   DDBA1280
C        END OF DIFFERENTIATION LOOP                                    DDBA1290
C                                                                       DDBA1300
   17 RETURN                                                            DDBA1310
      END                                                               DDBA1320
