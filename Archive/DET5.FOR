C                                                                       DET5  10
C     ..................................................................DET5  20
C                                                                       DET5  30
C        SUBROUTINE DET5                                                DET5  40
C                                                                       DET5  50
C        PURPOSE                                                        DET5  60
C           TO COMPUTE A VECTOR OF DERIVATIVE VALUES GIVEN A VECTOR OF  DET5  70
C           FUNCTION VALUES WHOSE ENTRIES CORRESPOND TO EQUIDISTANTLY   DET5  80
C           SPACED ARGUMENT VALUES.                                     DET5  90
C                                                                       DET5 100
C        USAGE                                                          DET5 110
C           CALL DET5(H,Y,Z,NDIM,IER)                                   DET5 120
C                                                                       DET5 130
C        DESCRIPTION OF PARAMETERS                                      DET5 140
C           H     -  CONSTANT DIFFERENCE BETWEEN SUCCESSIVE ARGUMENT    DET5 150
C                    VALUES (H IS POSITIVE IF THE ARGUMENT VALUES       DET5 160
C                    INCREASE AND NEGATIVE OTHERWISE)                   DET5 170
C           Y     -  GIVEN VECTOR OF FUNCTION VALUES (DIMENSION NDIM)   DET5 180
C           Z     -  RESULTING VECTOR OF DERIVATIVE VALUES (DIMENSION   DET5 190
C                    NDIM)                                              DET5 200
C           NDIM  -  DIMENSION OF VECTORS Y AND Z                       DET5 210
C           IER   -  RESULTING ERROR PARAMETER                          DET5 220
C                    IER = -1  - NDIM IS LESS THAN 5                    DET5 230
C                    IER =  0  - NO ERROR                               DET5 240
C                    IER =  1  - H = 0                                  DET5 250
C                                                                       DET5 260
C        REMARKS                                                        DET5 270
C           (1)   IF IER = -1,1, THEN THERE IS NO COMPUTATION.          DET5 280
C           (2)   Z CAN HAVE THE SAME STORAGE ALLOCATION AS Y.  IF Y IS DET5 290
C                 DISTINCT FROM Z, THEN IT IS NOT DESTROYED.            DET5 300
C                                                                       DET5 310
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DET5 320
C           NONE                                                        DET5 330
C                                                                       DET5 340
C        METHOD                                                         DET5 350
C           IF X IS THE (SUPPRESSED) VECTOR OF ARGUMENT VALUES, THEN    DET5 360
C           EXCEPT AT THE POINTS X(1),X(2),X(NDIM-1) AND X(NDIM), Z(I)  DET5 370
C           IS THE DERIVATIVE AT X(I) OF THE LAGRANGIAN INTERPOLATION   DET5 380
C           POLYNOMIAL OF DEGREE 4 RELEVANT TO THE 5 SUCCESSIVE POINTS  DET5 390
C           (X(I+K),Y(I+K)) K = -2,-1,...,2.  (SEE HILDEBRAND, F.B.,    DET5 400
C           INTRODUCTION TO NUMERICAL ANALYSIS, MC GRAW-HILL, NEW YORK/ DET5 410
C           TORONTO/LONDON, 1956, PP. 82-84.)                           DET5 420
C                                                                       DET5 430
C     ..................................................................DET5 440
C                                                                       DET5 450
      SUBROUTINE DET5(H,Y,Z,NDIM,IER)                                   DET5 460
C                                                                       DET5 470
C                                                                       DET5 480
      DIMENSION Y(1),Z(1)                                               DET5 490
C                                                                       DET5 500
C        TEST OF DIMENSION                                              DET5 510
      IF(NDIM-5)4,1,1                                                   DET5 520
C                                                                       DET5 530
C        TEST OF STEPSIZE                                               DET5 540
    1 IF(H)2,5,2                                                        DET5 550
C                                                                       DET5 560
C        PREPARE DIFFERENTIATION LOOP                                   DET5 570
    2 HH=.08333333/H                                                    DET5 580
      YY=Y(NDIM-4)                                                      DET5 590
      B=HH*(-25.*Y(1)+48.*Y(2)-36.*Y(3)+16.*Y(4)-3.*Y(5))               DET5 600
      C=HH*(-3.*Y(1)-10.*Y(2)+18.*Y(3)-6.*Y(4)+Y(5))                    DET5 610
C                                                                       DET5 620
C        START DIFFERENTIATION LOOP                                     DET5 630
      DO 3 I=5,NDIM                                                     DET5 640
      A=B                                                               DET5 650
      B=C                                                               DET5 660
      C=HH*(Y(I-4)-Y(I)+8.*(Y(I-1)-Y(I-3)))                             DET5 670
    3 Z(I-4)=A                                                          DET5 680
C        END OF DIFFERENTIATION LOOP                                    DET5 690
C                                                                       DET5 700
C        NORMAL EXIT                                                    DET5 710
      IER=0                                                             DET5 720
      A=HH*(-YY+6.*Y(NDIM-3)-18.*Y(NDIM-2)+10.*Y(NDIM-1)+3.*Y(NDIM))    DET5 730
     0Z(NDIM)=HH*(3.*YY-16.*Y(NDIM-3)+36.*Y(NDIM-2)-48.*Y(NDIM-1)       DET5 740
     1            +25.*Y(NDIM))                                         DET5 750
      Z(NDIM-1)=A                                                       DET5 760
      Z(NDIM-2)=C                                                       DET5 770
      Z(NDIM-3)=B                                                       DET5 780
      RETURN                                                            DET5 790
C                                                                       DET5 800
C        ERROR EXIT IN CASE NDIM IS LESS THAN 5                         DET5 810
    4 IER=-1                                                            DET5 820
      RETURN                                                            DET5 830
C                                                                       DET5 840
C        ERROR EXIT IN CASE OF ZERO STEPSIZE                            DET5 850
    5 IER=1                                                             DET5 860
      RETURN                                                            DET5 870
      END                                                               DET5 880
