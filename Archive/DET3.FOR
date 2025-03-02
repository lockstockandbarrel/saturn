C                                                                       DET3  10
C     ..................................................................DET3  20
C                                                                       DET3  30
C        SUBROUTINE DET3                                                DET3  40
C                                                                       DET3  50
C        PURPOSE                                                        DET3  60
C           TO COMPUTE A VECTOR OF DERIVATIVE VALUES GIVEN A VECTOR OF  DET3  70
C           FUNCTION VALUES WHOSE ENTRIES CORRESPOND TO EQUIDISTANTLY   DET3  80
C           SPACED ARGUMENT VALUES.                                     DET3  90
C                                                                       DET3 100
C        USAGE                                                          DET3 110
C           CALL DET3(H,Y,Z,NDIM,IER)                                   DET3 120
C                                                                       DET3 130
C        DESCRIPTION OF PARAMETERS                                      DET3 140
C           H     -  CONSTANT DIFFERENCE BETWEEN SUCCESSIVE ARGUMENT    DET3 150
C                    VALUES (H IS POSITIVE IF THE ARGUMENT VALUES       DET3 160
C                    INCREASE AND NEGATIVE OTHERWISE)                   DET3 170
C           Y     -  GIVEN VECTOR OF FUNCTION VALUES (DIMENSION NDIM)   DET3 180
C           Z     -  RESULTING VECTOR OF DERIVATIVE VALUES (DIMENSION   DET3 190
C                    NDIM)                                              DET3 200
C           NDIM  -  DIMENSION OF VECTORS Y AND Z                       DET3 210
C           IER   -  RESULTING ERROR PARAMETER                          DET3 220
C                    IER = -1  - NDIM IS LESS THAN 3                    DET3 230
C                    IER =  0  - NO ERROR                               DET3 240
C                    IER =  1  - H = 0                                  DET3 250
C                                                                       DET3 260
C        REMARKS                                                        DET3 270
C           (1)   IF IER = -1,1, THEN THERE IS NO COMPUTATION.          DET3 280
C           (2)   Z CAN HAVE THE SAME STORAGE ALLOCATION AS Y. IF Y IS  DET3 290
C                 DISTINCT FROM Z, THEN IT IS NOT DESTROYED.            DET3 300
C                                                                       DET3 310
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DET3 320
C           NONE                                                        DET3 330
C                                                                       DET3 340
C        METHOD                                                         DET3 350
C           IF X IS THE (SUPPRESSED) VECTOR OF ARGUMENT VALUES, THEN    DET3 360
C           EXCEPT AT THE ENDPOINTS X(1) AND X(NDIM), Z(I) IS THE       DET3 370
C           DERIVATIVE AT X(I) OF THE LAGRANGIAN INTERPOLATION          DET3 380
C           POLYNOMIAL OF DEGREE 2 RELEVANT TO THE 3 SUCCESSIVE POINTS  DET3 390
C           (X(I+K),Y(I+K)) K = -1,0,1.  (SEE HILDEBRAND, F.B.,         DET3 400
C           INTRODUCTION TO NUMERICAL ANALYSIS, MC-GRAW-HILL, NEW YORK/ DET3 410
C           TORONTO/LONDON, 1956, PP.82-84.)                            DET3 420
C                                                                       DET3 430
C     ..................................................................DET3 440
C                                                                       DET3 450
      SUBROUTINE DET3(H,Y,Z,NDIM,IER)                                   DET3 460
C                                                                       DET3 470
C                                                                       DET3 480
      DIMENSION Y(1),Z(1)                                               DET3 490
C                                                                       DET3 500
C        TEST OF DIMENSION                                              DET3 510
      IF(NDIM-3)4,1,1                                                   DET3 520
C                                                                       DET3 530
C        TEST OF STEPSIZE                                               DET3 540
    1 IF(H)2,5,2                                                        DET3 550
C                                                                       DET3 560
C        PREPARE DIFFERENTIATION LOOP                                   DET3 570
    2 HH=.5/H                                                           DET3 580
      YY=Y(NDIM-2)                                                      DET3 590
      B=Y(2)+Y(2)                                                       DET3 600
      B=HH*(B+B-Y(3)-Y(1)-Y(1)-Y(1))                                    DET3 610
C                                                                       DET3 620
C        START DIFFERENTIATION LOOP                                     DET3 630
      DO 3 I=3,NDIM                                                     DET3 640
      A=B                                                               DET3 650
      B=HH*(Y(I)-Y(I-2))                                                DET3 660
    3 Z(I-2)=A                                                          DET3 670
C        END OF DIFFERENTIATION LOOP                                    DET3 680
C                                                                       DET3 690
C        NORMAL EXIT                                                    DET3 700
      IER=0                                                             DET3 710
      A=Y(NDIM-1)+Y(NDIM-1)                                             DET3 720
      Z(NDIM)=HH*(Y(NDIM)+Y(NDIM)+Y(NDIM)-A-A+YY)                       DET3 730
      Z(NDIM-1)=B                                                       DET3 740
      RETURN                                                            DET3 750
C                                                                       DET3 760
C        ERROR EXIT IN CASE NDIM IS LESS THAN 3                         DET3 770
    4 IER=-1                                                            DET3 780
      RETURN                                                            DET3 790
C                                                                       DET3 800
C        ERROR EXIT IN CASE OF ZERO STEPSIZE                            DET3 810
    5 IER=1                                                             DET3 820
      RETURN                                                            DET3 830
      END                                                               DET3 840
