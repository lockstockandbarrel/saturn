C                                                                       DDCA  10
C     ..................................................................DDCA  20
C                                                                       DDCA  30
C     SUBROUTINE DDCAR                                                  DDCA  40
C                                                                       DDCA  50
C     PURPOSE                                                           DDCA  60
C        TO COMPUTE, AT A GIVEN POINT X, AN APPROXIMATION Z TO THE      DDCA  70
C        DERIVATIVE OF AN ANALYTICALLY GIVEN FUNCTION FCT THAT IS 11-   DDCA  80
C        TIMES DIFFERENTIABLE IN A DOMAIN CONTAINING A CLOSED, 2-SIDED  DDCA  90
C        SYMMETRIC INTERVAL OF RADIUS ABSOLUTE H ABOUT X, USING FUNCTIONDDCA 100
C        VALUES ONLY ON THAT CLOSED INTERVAL.                           DDCA 110
C                                                                       DDCA 120
C     USAGE                                                             DDCA 130
C        CALL DDCAR(X,H,IH,FCT,Z)                                       DDCA 140
C        PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                   DDCA 150
C                                                                       DDCA 160
C     DESCRIPTION OF PARAMETERS                                         DDCA 170
C        X   - THE POINT AT WHICH THE DERIVATIVE IS TO BE COMPUTED      DDCA 180
C              X IS IN DOUBLE PRECISION.                                DDCA 190
C        H   - THE NUMBER WHOSE ABSOLUTE VALUE DEFINES THE CLOSED,      DDCA 200
C              SYMMETRIC 2-SIDED INTERVAL ABOUT X (SEE PURPOSE)         DDCA 210
C              H IS IN SINGLE PRECISION                                 DDCA 220
C        IH  - INPUT PARAMETER (SEE REMARKS AND METHOD)                 DDCA 230
C              IH NON-ZERO - THE SUBROUTINE GENERATES THE INTERNAL      DDCA 240
C                            VALUE HH                                   DDCA 250
C              IH    =   0 - THE INTERNAL VALUE HH IS SET TO ABSOLUTE H DDCA 260
C        FCT - THE NAME OF THE EXTERNAL DOUBLE PRECISION FUNCTION       DDCA 270
C              SUBPROGRAM THAT WILL GENERATE THE NECESSARY FUNCTION     DDCA 280
C              VALUES.                                                  DDCA 290
C        Z   - RESULTING DERIVATIVE VALUE - DOUBLE PRECISION            DDCA 300
C                                                                       DDCA 310
C     REMARKS                                                           DDCA 320
C        (1)  IF H = 0, THEN THERE IS NO COMPUTATION.                   DDCA 330
C        (2)  THE INTERNAL VALUE HH, WHICH IS DETERMINED ACCORDING TO   DDCA 340
C             IH, IS THE MAXIMUM STEP-SIZE USED IN THE COMPUTATION OF   DDCA 350
C             THE CENTRAL DIVIDED DIFFERENCES (SEE METHOD.)  IF IH IS   DDCA 360
C             NON-ZERO, THEN THE SUBROUTINE GENERATES HH ACCORDING TO   DDCA 370
C             CRITERIA THAT BALANCE ROUND-OFF AND TRUNCATION ERROR.  HH DDCA 380
C             IS ALWAYS LESS THAN OR EQUAL TO ABSOLUTE H IN ABSOLUTE    DDCA 390
C             VALUE, SO THAT ALL COMPUTATION OCCURS WITHIN A RADIUS     DDCA 400
C             ABSOLUTE H OF X.                                          DDCA 410
C                                                                       DDCA 420
C     SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                     DDCA 430
C        THE EXTERNAL FUNCTION SUBPROGRAM FCT(T) MUST BE FURNISHED BY   DDCA 440
C        THE USER. FCT(T) IS IN DOUBLE PRECISION                        DDCA 450
C                                                                       DDCA 460
C     METHOD                                                            DDCA 470
C        THE COMPUTATION OF Z IS BASED ON RICHARDSON'S AND ROMBERG'S    DDCA 480
C        EXTRAPOLATION METHOD AS APPLIED TO THE SEQUENCE OF CENTRAL     DDCA 490
C        DIVIDED DIFFERENCES ASSOCIATED WITH THE POINT PAIRS            DDCA 500
C        (X-(K*HH)/5,X+(K*HH)/5) K=1,...,5.  (SEE FILLIPI, S. AND       DDCA 510
C        ENGELS, H., ALTES UND NEUES ZUR NUMERISCHEN DIFFERENTIATION,   DDCA 520
C        ELECTRONISCHE DATENVERARBEITUNG, ISS. 2 (1966), PP. 57-65.)    DDCA 530
C                                                                       DDCA 540
C     ..................................................................DDCA 550
C                                                                       DDCA 560
      SUBROUTINE DDCAR(X,H,IH,FCT,Z)                                    DDCA 570
C                                                                       DDCA 580
C                                                                       DDCA 590
      DIMENSION AUX(5)                                                  DDCA 600
      DOUBLE PRECISION X,FCT,Z,AUX,A,B,C,DH,HH                          DDCA 610
C                                                                       DDCA 620
C        NO ACTION IN CASE OF ZERO INTERVAL LENGTH                      DDCA 630
      IF(H)1,17,1                                                       DDCA 640
C                                                                       DDCA 650
C        GENERATE STEPSIZE HH FOR DIVIDED DIFFERENCES                   DDCA 660
    1 C=ABS(H)                                                          DDCA 670
      IF(IH)2,9,2                                                       DDCA 680
    2 HH=.5D-2                                                          DDCA 690
      IF(C-HH)3,4,4                                                     DDCA 700
    3 HH=C                                                              DDCA 710
    4 A=FCT(X+HH)                                                       DDCA 720
      B=FCT(X-HH)                                                       DDCA 730
      Z=DABS((A-B)/(HH+HH))                                             DDCA 740
      A=.5D0*DABS(A+B)                                                  DDCA 750
      HH=.5D-2                                                          DDCA 760
      IF(A-1.D0)6,6,5                                                   DDCA 770
    5 HH=HH*A                                                           DDCA 780
    6 IF(Z-1.D0)8,8,7                                                   DDCA 790
    7 HH=HH/Z                                                           DDCA 800
    8 IF(HH-C)10,10,9                                                   DDCA 810
    9 HH=C                                                              DDCA 820
C                                                                       DDCA 830
C        INITIALIZE DIFFERENTIATION LOOP                                DDCA 840
   10 Z=(FCT(X+HH)-FCT(X-HH))/(HH+HH)                                   DDCA 850
      J=5                                                               DDCA 860
      JJ=J-1                                                            DDCA 870
      AUX(J)=Z                                                          DDCA 880
      DH=HH/DFLOAT(J)                                                   DDCA 890
      DZ=1.E75                                                          DDCA 900
C                                                                       DDCA 910
C        START DIFFERENTIATION LOOP                                     DDCA 920
   11 J=J-1                                                             DDCA 930
      C=J                                                               DDCA 940
      HH=C*DH                                                           DDCA 950
      AUX(J)=(FCT(X+HH)-FCT(X-HH))/(HH+HH)                              DDCA 960
C                                                                       DDCA 970
C        INITIALIZE EXTRAPOLATION LOOP                                  DDCA 980
      D2=1.E75                                                          DDCA 990
      B=0.D0                                                            DDCA1000
      A=1.D0/C                                                          DDCA1010
C                                                                       DDCA1020
C        START EXTRAPOLATION LOOP                                       DDCA1030
      DO 12 I=J,JJ                                                      DDCA1040
      D1=D2                                                             DDCA1050
      B=B+A                                                             DDCA1060
      HH=(AUX(I)-AUX(I+1))/(B*(2.D0+B))                                 DDCA1070
      AUX(I+1)=AUX(I)+HH                                                DDCA1080
C                                                                       DDCA1090
C        TEST ON OSCILLATING INCREMENTS                                 DDCA1100
      D2=DABS(HH)                                                       DDCA1110
      IF(D2-D1)12,13,13                                                 DDCA1120
   12 CONTINUE                                                          DDCA1130
C        END OF EXTRAPOLATION LOOP                                      DDCA1140
C                                                                       DDCA1150
C        UPDATE RESULT VALUE Z                                          DDCA1160
      I=JJ+1                                                            DDCA1170
      GO TO 14                                                          DDCA1180
   13 D2=D1                                                             DDCA1190
      JJ=I                                                              DDCA1200
   14 IF(D2-DZ)15,16,16                                                 DDCA1210
   15 DZ=D2                                                             DDCA1220
      Z=AUX(I)                                                          DDCA1230
   16 IF(J-1)17,17,11                                                   DDCA1240
C        END OF DIFFERENTIATION LOOP                                    DDCA1250
C                                                                       DDCA1260
   17 RETURN                                                            DDCA1270
      END                                                               DDCA1280
