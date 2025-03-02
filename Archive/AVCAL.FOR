C                                                                       AVCA  10
C     ..................................................................AVCA  20
C                                                                       AVCA  30
C        SUBROUTINE AVCAL                                               AVCA  40
C                                                                       AVCA  50
C        PURPOSE                                                        AVCA  60
C           PERFORM THE CALCULUS OF A FACTORIAL EXPERIMENT USING        AVCA  70
C           OPERATOR SIGMA AND OPERATOR DELTA.  THIS SUBROUTINE IS      AVCA  80
C           PRECEDED BY SUBROUTINE ADVAT AND FOLLOWED BY SUBROUTINE     AVCA  90
C           MEANQ IN THE PERFORMANCE OF ANALYSIS OF VARIANCE FOR A      AVCA 100
C           COMPLETE FACTORIAL DESIGN.                                  AVCA 110
C                                                                       AVCA 120
C        USAGE                                                          AVCA 130
C           CALL AVCAL (K,LEVEL,X,L,ISTEP,LASTS)                        AVCA 140
C                                                                       AVCA 150
C        DESCRIPTION OF PARAMETERS                                      AVCA 160
C           K     - NUMBER OF VARIABLES (FACTORS). K MUST BE .GT. ONE.  AVCA 170
C           LEVEL - INPUT VECTOR OF LENGTH K CONTAINING LEVELS (CATE-   AVCA 180
C                   GORIES) WITHIN EACH VARIABLE.                       AVCA 190
C           X     - INPUT VECTOR CONTAINING DATA.  DATA HAVE BEEN PLACEDAVCA 200
C                   IN VECTOR X BY SUBROUTINE AVDAT.  THE LENGTH OF X   AVCA 210
C                   IS (LEVEL(1)+1)*(LEVEL(2)+1)*...*(LEVEL(K)+1).      AVCA 220
C           L     - THE POSITION IN VECTOR X WHERE THE LAST INPUT DATA  AVCA 230
C                   IS LOCATED.  L HAS BEEN CALCULATED BY SUBROUTINE    AVCA 240
C                   AVDAT.                                              AVCA 250
C           ISTEP - INPUT VECTOR OF LENGTH K CONTAINING STORAGE CONTROL AVCA 260
C                   STEPS WHICH HAVE BEEN CALCULATED BY SUBROUTINE      AVCA 270
C                   AVDAT.                                              AVCA 280
C           LASTS - WORKING VECTOR OF LENGTH K.                         AVCA 290
C                                                                       AVCA 300
C        REMARKS                                                        AVCA 310
C           THIS SUBROUTINE MUST FOLLOW SUBROUTINE AVDAT.               AVCA 320
C                                                                       AVCA 330
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  AVCA 340
C           NONE                                                        AVCA 350
C                                                                       AVCA 360
C        METHOD                                                         AVCA 370
C           THE METHOD IS BASED ON THE TECHNIQUE DISCUSSED BY H. O.     AVCA 380
C           HARTLEY IN 'MATHEMATICAL METHODS FOR DIGITAL COMPUTERS',    AVCA 390
C           EDITED BY A. RALSTON AND H. WILF, JOHN WILEY AND SONS,      AVCA 400
C           1962, CHAPTER 20.                                           AVCA 410
C                                                                       AVCA 420
C     ..................................................................AVCA 430
C                                                                       AVCA 440
      SUBROUTINE AVCAL (K,LEVEL,X,L,ISTEP,LASTS)                        AVCA 450
      DIMENSION LEVEL(1),X(1),ISTEP(1),LASTS(1)                         AVCA 460
C                                                                       AVCA 470
C        ...............................................................AVCA 480
C                                                                       AVCA 490
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE  AVCA 500
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION      AVCA 510
C        STATEMENT WHICH FOLLOWS.                                       AVCA 520
C                                                                       AVCA 530
C     DOUBLE PRECISION X,SUM                                            AVCA 540
C                                                                       AVCA 550
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS    AVCA 560
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS      AVCA 570
C        ROUTINE.                                                       AVCA 580
C                                                                       AVCA 590
C        ...............................................................AVCA 600
C                                                                       AVCA 610
C     CALCULATE THE LAST DATA POSITION OF EACH FACTOR                   AVCA 620
C                                                                       AVCA 630
      LASTS(1)=L+1                                                      AVCA 640
      DO 145 I=2,K                                                      AVCA 650
  145 LASTS(I)=LASTS(I-1)+ISTEP(I)                                      AVCA 660
C                                                                       AVCA 670
C     PERFORM CALCULUS OF OPERATION                                     AVCA 680
C                                                                       AVCA 690
  150 DO 175 I=1,K                                                      AVCA 700
      L=1                                                               AVCA 710
      LL=1                                                              AVCA 720
      SUM=0.0                                                           AVCA 730
      NN=LEVEL(I)                                                       AVCA 740
      FN=NN                                                             AVCA 750
      INCRE=ISTEP(I)                                                    AVCA 760
      LAST=LASTS(I)                                                     AVCA 770
C                                                                       AVCA 780
C     SIGMA OPERATION                                                   AVCA 790
C                                                                       AVCA 800
  155 DO 160 J=1,NN                                                     AVCA 810
      SUM=SUM+X(L)                                                      AVCA 820
  160 L=L+INCRE                                                         AVCA 830
      X(L)=SUM                                                          AVCA 840
C                                                                       AVCA 850
C     DELTA OPERATION                                                   AVCA 860
C                                                                       AVCA 870
      DO 165 J=1,NN                                                     AVCA 880
      X(LL)=FN*X(LL)-SUM                                                AVCA 890
  165 LL=LL+INCRE                                                       AVCA 900
      SUM=0.0                                                           AVCA 910
      IF(L-LAST) 167, 175, 175                                          AVCA 920
  167 IF(L-LAST+INCRE) 168, 168, 170                                    AVCA 930
  168 L=L+INCRE                                                         AVCA 940
      LL=LL+INCRE                                                       AVCA 950
      GO TO 155                                                         AVCA 960
  170 L=L+INCRE+1-LAST                                                  AVCA 970
      LL=LL+INCRE+1-LAST                                                AVCA 980
      GO TO 155                                                         AVCA 990
  175 CONTINUE                                                          AVCA1000
      RETURN                                                            AVCA1010
      END                                                               AVCA1020
