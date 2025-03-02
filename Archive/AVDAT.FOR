C                                                                       AVDA  10
C     ..................................................................AVDA  20
C                                                                       AVDA  30
C        SUBROUTINE AVDAT                                               AVDA  40
C                                                                       AVDA  50
C        PURPOSE                                                        AVDA  60
C           PLACE DATA FOR ANALYSIS OF VARIANCE IN PROPERLY DISTRIBUTED AVDA  70
C           POSITIONS OF STORAGE.  THIS SUBROUTINE IS NORMALLY FOLLOWED AVDA  80
C           BY CALLS TO AVCAL AND MEANQ SUBROUTINES IN THE PERFORMANCE  AVDA  90
C           OF ANALYSIS OF VARIANCE FOR A COMPLETE FACTORIAL DESIGN.    AVDA 100
C                                                                       AVDA 110
C        USAGE                                                          AVDA 120
C           CALL AVDAT (K,LEVEL,N,X,L,ISTEP,KOUNT)                      AVDA 130
C                                                                       AVDA 140
C        DESCRIPTION OF PARAMETERS                                      AVDA 150
C           K     - NUMBER OF VARIABLES (FACTORS). K MUST BE .GT. ONE.  AVDA 160
C           LEVEL - INPUT VECTOR OF LENGTH K CONTAINING LEVELS (CATE-   AVDA 170
C                   GORIES) WITHIN EACH VARIABLE.                       AVDA 180
C           N     - TOTAL NUMBER OF DATA POINTS READ IN.                AVDA 190
C           X     - WHEN THE SUBROUTINE IS CALLED, THIS VECTOR CONTAINS AVDA 200
C                   DATA IN LOCATIONS X(1) THROUGH X(N).  UPON RETURNINGAVDA 210
C                   TO THE CALLING ROUTINE, THE VECTOR CONTAINS THE DATAAVDA 220
C                   IN PROPERLY REDISTRIBUTED LOCATIONS OF VECTOR X.    AVDA 230
C                   THE LENGTH OF VECTOR X IS CALCULATED BY (1) ADDING  AVDA 240
C                   ONE TO EACH LEVEL OF VARIABLE AND (2) OBTAINING THE AVDA 250
C                   CUMULATIVE PRODUCT OF ALL LEVELS.  (THE LENGTH OF   AVDA 260
C                   X = (LEVEL(1)+1)*(LEVEL(2)+1)*...*(LEVEL(K)+1).)    AVDA 270
C           L     - OUTPUT VARIABLE CONTAINING THE POSITION IN VECTOR X AVDA 280
C                   WHERE THE LAST INPUT DATA IS STORED.                AVDA 290
C           ISTEP - OUTPUT VECTOR OF LENGTH K CONTAINING CONTROL STEPS  AVDA 300
C                   WHICH ARE USED TO LOCATE DATA IN PROPER POSITIONS   AVDA 310
C                   OF VECTOR X.                                        AVDA 320
C           KOUNT - WORKING VECTOR OF LENGTH K.                         AVDA 330
C                                                                       AVDA 340
C        REMARKS                                                        AVDA 350
C           INPUT DATA MUST BE ARRANGED IN THE FOLLOWING MANNER.        AVDA 360
C           CONSIDER THE 3-VARIABLE ANALYSIS OF VARIANCE DESIGN, WHERE  AVDA 370
C           ONE VARIABLE HAS 3 LEVELS AND THE OTHER TWO VARIABLES HAVE  AVDA 380
C           2 LEVELS.  THE DATA MAY BE REPRESENTED IN THE FORM X(I,J,K),AVDA 390
C           I=1,2,3  J=1,2  K=1,2.  IN ARRANGING DATA, THE INNER        AVDA 400
C           SUBSCRIPT, NAMELY I, CHANGES FIRST.  WHEN I=3, THE NEXT     AVDA 410
C           INNER SUBSCRIPT, J, CHANGES AND SO ON UNTIL I=3, J=2, AND   AVDA 420
C           K=2.                                                        AVDA 430
C                                                                       AVDA 440
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  AVDA 450
C           NONE                                                        AVDA 460
C                                                                       AVDA 470
C        METHOD                                                         AVDA 480
C           THE METHOD IS BASED ON THE TECHNIQUE DISCUSSED BY H. O.     AVDA 490
C           HARTLEY IN 'MATHEMATICAL METHODS FOR DIGITAL COMPUTERS',    AVDA 500
C           EDITED BY A. RALSTON AND H. WILF, JOHN WILEY AND SONS,      AVDA 510
C           1962, CHAPTER 20.                                           AVDA 520
C                                                                       AVDA 530
C     ..................................................................AVDA 540
C                                                                       AVDA 550
      SUBROUTINE AVDAT (K,LEVEL,N,X,L,ISTEP,KOUNT)                      AVDA 560
      DIMENSION LEVEL(1),X(1),ISTEP(1),KOUNT(1)                         AVDA 570
C                                                                       AVDA 580
C        ...............................................................AVDA 590
C                                                                       AVDA 600
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE  AVDA 610
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION      AVDA 620
C        STATEMENT WHICH FOLLOWS.                                       AVDA 630
C                                                                       AVDA 640
C     DOUBLE PRECISION X                                                AVDA 650
C                                                                       AVDA 660
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS    AVDA 670
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS      AVDA 680
C        ROUTINE.                                                       AVDA 690
C                                                                       AVDA 700
C        ...............................................................AVDA 710
C                                                                       AVDA 720
C     CALCULATE TOTAL DATA AREA REQUIRED                                AVDA 730
C                                                                       AVDA 740
      M=LEVEL(1)+1                                                      AVDA 750
      DO 105 I=2,K                                                      AVDA 760
  105 M=M*(LEVEL(I)+1)                                                  AVDA 770
C                                                                       AVDA 780
C     MOVE DATA TO THE UPPER PART OF THE ARRAY X                        AVDA 790
C     FOR THE PURPOSE OF REARRANGEMENT                                  AVDA 800
C                                                                       AVDA 810
      N1=M+1                                                            AVDA 820
      N2=N+1                                                            AVDA 830
      DO 107 I=1,N                                                      AVDA 840
      N1=N1-1                                                           AVDA 850
      N2=N2-1                                                           AVDA 860
  107 X(N1)=X(N2)                                                       AVDA 870
C                                                                       AVDA 880
C     CALCULATE MULTIPLIERS TO BE USED IN FINDING STORAGE LOCATIONS FOR AVDA 890
C     INPUT DATA                                                        AVDA 900
C                                                                       AVDA 910
      ISTEP(1)=1                                                        AVDA 920
      DO 110 I=2,K                                                      AVDA 930
  110 ISTEP(I)=ISTEP(I-1)*(LEVEL(I-1)+1)                                AVDA 940
      DO 115 I=1,K                                                      AVDA 950
  115 KOUNT(I)=1                                                        AVDA 960
C                                                                       AVDA 970
C     PLACE DATA IN PROPER LOCATIONS                                    AVDA 980
C                                                                       AVDA 990
      N1=N1-1                                                           AVDA1000
      DO 135 I=1,N                                                      AVDA1010
      L=KOUNT(1)                                                        AVDA1020
      DO 120 J=2,K                                                      AVDA1030
  120 L=L+ISTEP(J)*(KOUNT(J)-1)                                         AVDA1040
      N1=N1+1                                                           AVDA1050
      X(L)=X(N1)                                                        AVDA1060
      DO 130 J=1,K                                                      AVDA1070
      IF(KOUNT(J)-LEVEL(J)) 124, 125, 124                               AVDA1080
  124 KOUNT(J)=KOUNT(J)+1                                               AVDA1090
      GO TO 135                                                         AVDA1100
  125 KOUNT(J)=1                                                        AVDA1110
  130 CONTINUE                                                          AVDA1120
  135 CONTINUE                                                          AVDA1130
      RETURN                                                            AVDA1140
      END                                                               AVDA1150
