C                                                                       DGT3  10
C     ..................................................................DGT3  20
C                                                                       DGT3  30
C        SUBROUTINE DGT3                                                DGT3  40
C                                                                       DGT3  50
C        PURPOSE                                                        DGT3  60
C           TO COMPUTE A VECTOR OF DERIVATIVE VALUES GIVEN VECTORS OF   DGT3  70
C           ARGUMENT VALUES AND CORRESPONDING FUNCTION VALUES.          DGT3  80
C                                                                       DGT3  90
C        USAGE                                                          DGT3 100
C           CALL DGT3(X,Y,Z,NDIM,IER)                                   DGT3 110
C                                                                       DGT3 120
C        DESCRIPTION OF PARAMETERS                                      DGT3 130
C           X     -  GIVEN VECTOR OF ARGUMENT VALUES (DIMENSION NDIM)   DGT3 140
C           Y     -  GIVEN VECTOR OF FUNCTION VALUES CORRESPONDING TO X DGT3 150
C                    (DIMENSION NDIM)                                   DGT3 160
C           Z     -  RESULTING VECTOR OF DERIVATIVE VALUES (DIMENSION   DGT3 170
C                    NDIM)                                              DGT3 180
C           NDIM  -  DIMENSION OF VECTORS X,Y AND Z                     DGT3 190
C           IER   -  RESULTING ERROR PARAMETER                          DGT3 200
C                    IER  = -1  - NDIM IS LESS THAN 3                   DGT3 210
C                    IER  =  0  - NO ERROR                              DGT3 220
C                    IER POSITIVE  - X(IER) = X(IER-1) OR X(IER) =      DGT3 230
C                                    X(IER-2)                           DGT3 240
C                                                                       DGT3 250
C        REMARKS                                                        DGT3 260
C           (1)   IF IER = -1,2,3, THEN THERE IS NO COMPUTATION.        DGT3 270
C           (2)   IF IER =  4,...,N, THEN THE DERIVATIVE VALUES Z(1)    DGT3 280
C                 ,..., Z(IER-1) HAVE BEEN COMPUTED.                    DGT3 290
C           (3)   Z CAN HAVE THE SAME STORAGE ALLOCATION AS X OR Y.  IF DGT3 300
C                 X OR Y IS DISTINCT FROM Z, THEN IT IS NOT DESTROYED.  DGT3 310
C                                                                       DGT3 320
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DGT3 330
C           NONE                                                        DGT3 340
C                                                                       DGT3 350
C        METHOD                                                         DGT3 360
C           EXCEPT AT THE ENDPOINTS X(1) AND X(NDIM), Z(I) IS THE       DGT3 370
C           DERIVATIVE AT X(I) OF THE LAGRANGIAN INTERPOLATION          DGT3 380
C           POLYNOMIAL OF DEGREE 2 RELEVANT TO THE 3 SUCCESSIVE POINTS  DGT3 390
C           (X(I+K),Y(I+K)) K = -1,0,1. (SEE HILDEBRAND, F.B.,          DGT3 400
C           INTRODUCTION TO NUMERICAL ANALYSIS, MC GRAW-HILL, NEW YORK/ DGT3 410
C           TORONTO/LONDON, 1956, PP. 64-68.)                           DGT3 420
C                                                                       DGT3 430
C     ..................................................................DGT3 440
C                                                                       DGT3 450
      SUBROUTINE DGT3(X,Y,Z,NDIM,IER)                                   DGT3 460
C                                                                       DGT3 470
C                                                                       DGT3 480
      DIMENSION X(1),Y(1),Z(1)                                          DGT3 490
C                                                                       DGT3 500
C        TEST OF DIMENSION AND ERROR EXIT IN CASE NDIM IS LESS THAN 3   DGT3 510
      IER=-1                                                            DGT3 520
      IF(NDIM-3)8,1,1                                                   DGT3 530
C                                                                       DGT3 540
C        PREPARE DIFFERENTIATION LOOP                                   DGT3 550
    1 A=X(1)                                                            DGT3 560
      B=Y(1)                                                            DGT3 570
      I=2                                                               DGT3 580
      DY2=X(2)-A                                                        DGT3 590
      IF(DY2)2,9,2                                                      DGT3 600
    2 DY2=(Y(2)-B)/DY2                                                  DGT3 610
C                                                                       DGT3 620
C        START DIFFERENTIATION LOOP                                     DGT3 630
      DO 6 I=3,NDIM                                                     DGT3 640
      A=X(I)-A                                                          DGT3 650
      IF(A)3,9,3                                                        DGT3 660
    3 A=(Y(I)-B)/A                                                      DGT3 670
      B=X(I)-X(I-1)                                                     DGT3 680
      IF(B)4,9,4                                                        DGT3 690
    4 DY1=DY2                                                           DGT3 700
      DY2=(Y(I)-Y(I-1))/B                                               DGT3 710
      DY3=A                                                             DGT3 720
      A=X(I-1)                                                          DGT3 730
      B=Y(I-1)                                                          DGT3 740
      IF(I-3)5,5,6                                                      DGT3 750
    5 Z(1)=DY1+DY3-DY2                                                  DGT3 760
    6 Z(I-1)=DY1+DY2-DY3                                                DGT3 770
C        END DIFFERENTIATION LOOP                                       DGT3 780
C                                                                       DGT3 790
C        NORMAL EXIT                                                    DGT3 800
      IER=0                                                             DGT3 810
      I=NDIM                                                            DGT3 820
    7 Z(I)=DY2+DY3-DY1                                                  DGT3 830
    8 RETURN                                                            DGT3 840
C                                                                       DGT3 850
C        ERROR EXIT IN CASE OF IDENTICAL ARGUMENTS                      DGT3 860
    9 IER=I                                                             DGT3 870
      I=I-1                                                             DGT3 880
      IF(I-2)8,8,7                                                      DGT3 890
      END                                                               DGT3 900
