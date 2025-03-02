C                                                                       DDGT  10
C     ..................................................................DDGT  20
C                                                                       DDGT  30
C        SUBROUTINE DDGT3                                               DDGT  40
C                                                                       DDGT  50
C        PURPOSE                                                        DDGT  60
C           TO COMPUTE A VECTOR OF DERIVATIVE VALUES GIVEN VECTORS OF   DDGT  70
C           ARGUMENT VALUES AND CORRESPONDING FUNCTION VALUES.          DDGT  80
C                                                                       DDGT  90
C        USAGE                                                          DDGT 100
C           CALL DDGT3(X,Y,Z,NDIM,IER)                                  DDGT 110
C                                                                       DDGT 120
C        DESCRIPTION OF PARAMETERS                                      DDGT 130
C           X     -  GIVEN VECTOR OF DOUBLE PRECISION ARGUMENT VALUES   DDGT 140
C                    (DIMENSION NDIM)                                   DDGT 150
C           Y     -  GIVEN VECTOR OF DOUBLE PRECISION FUNCTION VALUES   DDGT 160
C                    CORRESPONDING TO X (DIMENSION NDIM)                DDGT 170
C           Z     -  RESULTING VECTOR OF DOUBLE PRECISION DERIVATIVE    DDGT 180
C                    VALUES (DIMENSION NDIM)                            DDGT 190
C           NDIM  -  DIMENSION OF VECTORS X,Y AND Z                     DDGT 200
C           IER   -  RESULTING ERROR PARAMETER                          DDGT 210
C                    IER  = -1  - NDIM IS LESS THAN 3                   DDGT 220
C                    IER  =  0  - NO ERROR                              DDGT 230
C                    IER POSITIVE  - X(IER) = X(IER-1) OR X(IER) =      DDGT 240
C                                    X(IER-2)                           DDGT 250
C                                                                       DDGT 260
C        REMARKS                                                        DDGT 270
C           (1)   IF IER = -1,2,3, THEN THERE IS NO COMPUTATION.        DDGT 280
C           (2)   IF IER =  4,...,N, THEN THE DERIVATIVE VALUES Z(1)    DDGT 290
C                 ,..., Z(IER-1) HAVE BEEN COMPUTED.                    DDGT 300
C           (3)   Z CAN HAVE THE SAME STORAGE ALLOCATION AS X OR Y.  IF DDGT 310
C                 X OR Y IS DISTINCT FROM Z, THEN IT IS NOT DESTROYED.  DDGT 320
C                                                                       DDGT 330
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DDGT 340
C           NONE                                                        DDGT 350
C                                                                       DDGT 360
C        METHOD                                                         DDGT 370
C           EXCEPT AT THE ENDPOINTS X(1) AND X(NDIM), Z(I) IS THE       DDGT 380
C           DERIVATIVE AT X(I) OF THE LAGRANGIAN INTERPOLATION          DDGT 390
C           POLYNOMIAL OF DEGREE 2 RELEVANT TO THE 3 SUCCESSIVE POINTS  DDGT 400
C           (X(I+K),Y(I+K)) K = -1,0,1. (SEE HILDEBRAND, F.B.,          DDGT 410
C           INTRODUCTION TO NUMERICAL ANALYSIS, MC GRAW-HILL, NEW YORK/ DDGT 420
C           TORONTO/LONDON, 1956, PP. 64-68.)                           DDGT 430
C                                                                       DDGT 440
C     ..................................................................DDGT 450
C                                                                       DDGT 460
      SUBROUTINE DDGT3(X,Y,Z,NDIM,IER)                                  DDGT 470
C                                                                       DDGT 480
C                                                                       DDGT 490
      DIMENSION X(1),Y(1),Z(1)                                          DDGT 500
      DOUBLE PRECISION X,Y,Z,DY1,DY2,DY3,A,B                            DDGT 510
C                                                                       DDGT 520
C        TEST OF DIMENSION AND ERROR EXIT IN CASE NDIM IS LESS THAN 3   DDGT 530
      IER=-1                                                            DDGT 540
      IF(NDIM-3)8,1,1                                                   DDGT 550
C                                                                       DDGT 560
C        PREPARE DIFFERENTIATION LOOP                                   DDGT 570
    1 A=X(1)                                                            DDGT 580
      B=Y(1)                                                            DDGT 590
      I=2                                                               DDGT 600
      DY2=X(2)-A                                                        DDGT 610
      IF(DY2)2,9,2                                                      DDGT 620
    2 DY2=(Y(2)-B)/DY2                                                  DDGT 630
C                                                                       DDGT 640
C        START DIFFERENTIATION LOOP                                     DDGT 650
      DO 6 I=3,NDIM                                                     DDGT 660
      A=X(I)-A                                                          DDGT 670
      IF(A)3,9,3                                                        DDGT 680
    3 A=(Y(I)-B)/A                                                      DDGT 690
      B=X(I)-X(I-1)                                                     DDGT 700
      IF(B)4,9,4                                                        DDGT 710
    4 DY1=DY2                                                           DDGT 720
      DY2=(Y(I)-Y(I-1))/B                                               DDGT 730
      DY3=A                                                             DDGT 740
      A=X(I-1)                                                          DDGT 750
      B=Y(I-1)                                                          DDGT 760
      IF(I-3)5,5,6                                                      DDGT 770
    5 Z(1)=DY1+DY3-DY2                                                  DDGT 780
    6 Z(I-1)=DY1+DY2-DY3                                                DDGT 790
C        END OF DIFFERENTIATION LOOP                                    DDGT 800
C                                                                       DDGT 810
C        NORMAL EXIT                                                    DDGT 820
      IER=0                                                             DDGT 830
      I=NDIM                                                            DDGT 840
    7 Z(I)=DY2+DY3-DY1                                                  DDGT 850
    8 RETURN                                                            DDGT 860
C                                                                       DDGT 870
C        ERROR EXIT IN CASE OF IDENTICAL ARGUMENTS                      DDGT 880
    9 IER=I                                                             DDGT 890
      I=I-1                                                             DDGT 900
      IF(I-2)8,8,7                                                      DDGT 910
      END                                                               DDGT 920
