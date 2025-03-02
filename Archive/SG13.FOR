C                                                                       SG13  10
C     ..................................................................SG13  20
C                                                                       SG13  30
C        SUBROUTINE SG13                                                SG13  40
C                                                                       SG13  50
C        PURPOSE                                                        SG13  60
C           TO COMPUTE A VECTOR OF SMOOTHED FUNCTION VALUES GIVEN       SG13  70
C           VECTORS OF ARGUMENT VALUES AND CORRESPONDING FUNCTION       SG13  80
C           VALUES.                                                     SG13  90
C                                                                       SG13 100
C        USAGE                                                          SG13 110
C           CALL SG13(X,Y,Z,NDIM,IER)                                   SG13 120
C                                                                       SG13 130
C        DESCRIPTION OF PARAMETERS                                      SG13 140
C           X     -  GIVEN VECTOR OF ARGUMENT VALUES (DIMENSION NDIM)   SG13 150
C           Y     -  GIVEN VECTOR OF FUNCTION VALUES CORRESPONDING TO X SG13 160
C                    (DIMENSION NDIM)                                   SG13 170
C           Z     -  RESULTING VECTOR OF SMOOTHED FUNCTION VALUES       SG13 180
C                    (DIMENSION NDIM)                                   SG13 190
C           NDIM  -  DIMENSION OF VECTORS X,Y,AND Z                     SG13 200
C           IER   -  RESULTING ERROR PARAMETER                          SG13 210
C                    IER = -1  - NDIM IS LESS THAN 3                    SG13 220
C                    IER =  0  - NO ERROR                               SG13 230
C                                                                       SG13 240
C        REMARKS                                                        SG13 250
C           (1)  IF IER=-1 THERE HAS BEEN NO COMPUTATION.               SG13 260
C           (2)   Z CAN HAVE THE SAME STORAGE ALLOCATION AS X OR Y.  IF SG13 270
C                 X OR Y IS DISTINCT FROM Z, THEN IT IS NOT DESTROYED.  SG13 280
C                                                                       SG13 290
C        SUBROUTINES AND SUBPROGRAMS REQUIRED                           SG13 300
C           NONE                                                        SG13 310
C                                                                       SG13 320
C        METHOD                                                         SG13 330
C           EXCEPT AT THE ENDPOINTS X(1) AND X(NDIM), EACH SMOOTHED     SG13 340
C           VALUE Z(I) IS OBTAINED BY EVALUATING AT X(I) THE LEAST-     SG13 350
C           SQUARES POLYNOMIAL OF DEGREE 1 RELEVANT TO THE 3 SUCCESSIVE SG13 360
C           POINTS (X(I+K),Y(I+K)) K = -1,0,1.(SEE HILDEBRAND, F.B.,    SG13 370
C           INTRODUCTION  TO NUMERICAL ANALYSIS, MC GRAW-HILL, NEW YORK/SG13 380
C           TORONTO/LONDON, 1956, PP.258-311.)                          SG13 390
C                                                                       SG13 400
C     ..................................................................SG13 410
C                                                                       SG13 420
      SUBROUTINE SG13(X,Y,Z,NDIM,IER)                                   SG13 430
C                                                                       SG13 440
C                                                                       SG13 450
      DIMENSION X(1),Y(1),Z(1)                                          SG13 460
C                                                                       SG13 470
C        TEST OF DIMENSION                                              SG13 480
      IF(NDIM-3)7,1,1                                                   SG13 490
C                                                                       SG13 500
C        START LOOP                                                     SG13 510
    1 DO 6 I=3,NDIM                                                     SG13 520
      XM=.3333333*(X(I-2)+X(I-1)+X(I))                                  SG13 530
      YM=.3333333*(Y(I-2)+Y(I-1)+Y(I))                                  SG13 540
      T1=X(I-2)-XM                                                      SG13 550
      T2=X(I-1)-XM                                                      SG13 560
      T3=X(I)-XM                                                        SG13 570
      XM=T1*T1+T2*T2+T3*T3                                              SG13 580
      IF(XM)3,3,2                                                       SG13 590
    2 XM=(T1*(Y(I-2)-YM)+T2*(Y(I-1)-YM)+T3*(Y(I)-YM))/XM                SG13 600
C                                                                       SG13 610
C        CHECK FIRST POINT                                              SG13 620
    3 IF(I-3)4,4,5                                                      SG13 630
    4 H=XM*T1+YM                                                        SG13 640
    5 Z(I-2)=H                                                          SG13 650
    6 H=XM*T2+YM                                                        SG13 660
C        END OF LOOP                                                    SG13 670
C                                                                       SG13 680
C        UPDATE LAST TWO COMPONENTS                                     SG13 690
      Z(NDIM-1)=H                                                       SG13 700
      Z(NDIM)=XM*T3+YM                                                  SG13 710
      IER=0                                                             SG13 720
      RETURN                                                            SG13 730
C                                                                       SG13 740
C        ERROR EXIT IN CASE NDIM IS LESS THAN 3                         SG13 750
    7 IER=-1                                                            SG13 760
      RETURN                                                            SG13 770
      END                                                               SG13 780
