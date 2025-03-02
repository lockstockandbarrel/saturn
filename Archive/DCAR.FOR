C                                                                       DCAR  10
C     ..................................................................DCAR  20
C                                                                       DCAR  30
C        SUBROUTINE DCAR                                                DCAR  40
C                                                                       DCAR  50
C     PURPOSE                                                           DCAR  60
C        TO COMPUTE, AT A GIVEN POINT X, AN APPROXIMATION Z TO THE      DCAR  70
C        DERIVATIVE OF AN ANALYTICALLY GIVEN FUNCTION FCT THAT IS 11-   DCAR  80
C        TIMES DIFFERENTIABLE IN A DOMAIN CONTAINING A CLOSED, 2-SIDED  DCAR  90
C        SYMMETRIC INTERVAL OF RADIUS ABSOLUTE H ABOUT X, USING FUNCTIONDCAR 100
C        VALUES ONLY ON THAT CLOSED INTERVAL.                           DCAR 110
C                                                                       DCAR 120
C     USAGE                                                             DCAR 130
C           CALL DCAR (X,H,IH,FCT,Z)                                    DCAR 140
C        PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                   DCAR 150
C                                                                       DCAR 160
C     DESCRIPTION OF PARAMETERS                                         DCAR 170
C        X   - THE POINT AT WHICH THE DERIVATIVE IS TO BE COMPUTED      DCAR 180
C        H   - THE NUMBER WHOSE ABSOLUTE VALUE DEFINES THE CLOSED,      DCAR 190
C              SYMMETRIC 2-SIDED INTERVAL ABOUT X (SEE PURPOSE)         DCAR 200
C        IH  - INPUT PARAMETER (SEE REMARKS AND METHOD)                 DCAR 210
C              IH NON-ZERO - THE SUBROUTINE GENERATES THE INTERNAL      DCAR 220
C                            VALUE HH                                   DCAR 230
C              IH    =   0 - THE INTERNAL VALUE HH IS SET TO ABSOLUTE H DCAR 240
C        FCT - THE NAME OF THE EXTERNAL FUNCTION SUBPROGRAM THAT WILL   DCAR 250
C              GENERATE THE NECESSARY FUNCTION VALUES                   DCAR 260
C        Z   - RESULTING DERIVATIVE VALUE                               DCAR 270
C                                                                       DCAR 280
C     REMARKS                                                           DCAR 290
C        (1)  IF H = 0, THEN THERE IS NO COMPUTATION.                   DCAR 300
C        (2)  THE INTERNAL VALUE HH, WHICH IS DETERMINED ACCORDING TO   DCAR 310
C             IH, IS THE MAXIMUM STEP-SIZE USED IN THE COMPUTATION OF   DCAR 320
C             THE CENTRAL DIVIDED DIFFERENCES (SEE METHOD.)  IF IH IS   DCAR 330
C             NON-ZERO, THEN THE SUBROUTINE GENERATES HH ACCORDING TO   DCAR 340
C             CRITERIA THAT BALANCE ROUND-OFF AND TRUNCATION ERROR.  HH DCAR 350
C             IS ALWAYS LESS THAN OR EQUAL TO ABSOLUTE H IN ABSOLUTE    DCAR 360
C             VALUE, SO THAT ALL COMPUTATION OCCURS WITHIN A RADIUS     DCAR 370
C             ABSOLUTE H OF X.                                          DCAR 380
C                                                                       DCAR 390
C     SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                     DCAR 400
C        THE EXTERNAL FUNCTION SUBPROGRAM FCT(T) MUST BE FURNISHED BY   DCAR 410
C        THE USER.                                                      DCAR 420
C                                                                       DCAR 430
C     METHOD                                                            DCAR 440
C        THE COMPUTATION OF Z IS BASED ON RICHARDSON'S AND ROMBERG'S    DCAR 450
C        EXTRAPOLATION METHOD AS APPLIED TO THE SEQUENCE OF CENTRAL     DCAR 460
C        DIVIDED DIFFERENCES ASSOCIATED WITH THE POINT PAIRS            DCAR 470
C        (X-(K*HH)/5,X+(K*HH)/5) K=1,...,5.  (SEE FILLIPI, S. AND       DCAR 480
C        ENGELS, H., ALTES UND NEUES ZUR NUMERISCHEN DIFFERENTIATION,   DCAR 490
C        ELECTRONISCHE DATENVERARBEITUNG, ISS. 2 (1966), PP. 57-65.)    DCAR 500
C                                                                       DCAR 510
C     ..................................................................DCAR 520
C                                                                       DCAR 530
      SUBROUTINE DCAR(X,H,IH,FCT,Z)                                     DCAR 540
C                                                                       DCAR 550
C                                                                       DCAR 560
      DIMENSION AUX(5)                                                  DCAR 570
C                                                                       DCAR 580
C        NO ACTION IN CASE OF ZERO INTERVAL LENGTH                      DCAR 590
      IF(H)1,17,1                                                       DCAR 600
C                                                                       DCAR 610
C        GENERATE STEPSIZE HH FOR DIVIDED DIFFERENCES                   DCAR 620
    1 C=ABS(H)                                                          DCAR 630
      IF(IH)2,9,2                                                       DCAR 640
    2 HH=.5                                                             DCAR 650
      IF(C-HH)3,4,4                                                     DCAR 660
    3 HH=C                                                              DCAR 670
    4 A=FCT(X+HH)                                                       DCAR 680
      B=FCT(X-HH)                                                       DCAR 690
      Z=ABS((A-B)/(HH+HH))                                              DCAR 700
      A=.5*ABS(A+B)                                                     DCAR 710
      HH=.5                                                             DCAR 720
      IF(A-1.)6,6,5                                                     DCAR 730
    5 HH=HH*A                                                           DCAR 740
    6 IF(Z-1.)8,8,7                                                     DCAR 750
    7 HH=HH/Z                                                           DCAR 760
    8 IF(HH-C)10,10,9                                                   DCAR 770
    9 HH=C                                                              DCAR 780
C                                                                       DCAR 790
C        INITIALIZE DIFFERENTIATION LOOP                                DCAR 800
   10 Z=(FCT(X+HH)-FCT(X-HH))/(HH+HH)                                   DCAR 810
      J=5                                                               DCAR 820
      JJ=J-1                                                            DCAR 830
      AUX(J)=Z                                                          DCAR 840
      DH=HH/FLOAT(J)                                                    DCAR 850
      DZ=1.E75                                                          DCAR 860
C                                                                       DCAR 870
C        START DIFFERENTIATION LOOP                                     DCAR 880
   11 J=J-1                                                             DCAR 890
      C=J                                                               DCAR 900
      HH=C*DH                                                           DCAR 910
      AUX(J)=(FCT(X+HH)-FCT(X-HH))/(HH+HH)                              DCAR 920
C                                                                       DCAR 930
C        INITIALIZE EXTRAPOLATION LOOP                                  DCAR 940
      D2=1.E75                                                          DCAR 950
      B=0.                                                              DCAR 960
      A=1./C                                                            DCAR 970
C                                                                       DCAR 980
C        START EXTRAPOLATION LOOP                                       DCAR 990
      DO 12 I=J,JJ                                                      DCAR1000
      D1=D2                                                             DCAR1010
      B=B+A                                                             DCAR1020
      HH=(AUX(I)-AUX(I+1))/(B*(2.+B))                                   DCAR1030
      AUX(I+1)=AUX(I)+HH                                                DCAR1040
C                                                                       DCAR1050
C        TEST ON OSCILLATING INCREMENTS                                 DCAR1060
      D2=ABS(HH)                                                        DCAR1070
      IF(D2-D1)12,13,13                                                 DCAR1080
   12 CONTINUE                                                          DCAR1090
C        END OF EXTRAPOLATION LOOP                                      DCAR1100
C                                                                       DCAR1110
C        UPDATE RESULT VALUE Z                                          DCAR1120
      I=JJ+1                                                            DCAR1130
      GO TO 14                                                          DCAR1140
   13 D2=D1                                                             DCAR1150
      JJ=I                                                              DCAR1160
   14 IF(D2-DZ)15,16,16                                                 DCAR1170
   15 DZ=D2                                                             DCAR1180
      Z=AUX(I)                                                          DCAR1190
   16 IF(J-1)17,17,11                                                   DCAR1200
C        END OF DIFFERENTIATION LOOP                                    DCAR1210
C                                                                       DCAR1220
   17 RETURN                                                            DCAR1230
      END                                                               DCAR1240
