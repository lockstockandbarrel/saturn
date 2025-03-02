C                                                                       DSG1  10
C     ..................................................................DSG1  20
C                                                                       DSG1  30
C        SUBROUTINE DSG13                                               DSG1  40
C                                                                       DSG1  50
C        PURPOSE                                                        DSG1  60
C           TO COMPUTE A VECTOR OF SMOOTHED FUNCTION VALUES GIVEN       DSG1  70
C           VECTORS OF ARGUMENT VALUES AND CORRESPONDING FUNCTION       DSG1  80
C           VALUES.                                                     DSG1  90
C                                                                       DSG1 100
C        USAGE                                                          DSG1 110
C           CALL DSG13(X,Y,Z,NDIM,IER)                                  DSG1 120
C                                                                       DSG1 130
C        DESCRIPTION OF PARAMETERS                                      DSG1 140
C           X     -  GIVEN VECTOR OF DOUBLE PRECISION ARGUMENT VALUES   DSG1 150
C                    (DIMENSION NDIM)                                   DSG1 160
C           Y     -  GIVEN VECTOR OF DOUBLE PRECISION FUNCTION VALUES   DSG1 170
C                    CORRESPONDING TO X (DIMENSION NDIM)                DSG1 180
C           Z     -  RESULTING VECTOR OF DOUBLE PRECISION SMOOTHED      DSG1 190
C                    FUNCTION VALUES (DIMENSION NDIM)                   DSG1 200
C           NDIM  -  DIMENSION OF VECTORS X,Y,AND Z                     DSG1 210
C           IER   -  RESULTING ERROR PARAMETER                          DSG1 220
C                    IER = -1  - NDIM IS LESS THAN 3                    DSG1 230
C                    IER =  0  - NO ERROR                               DSG1 240
C                                                                       DSG1 250
C        REMARKS                                                        DSG1 260
C           (1)  IF IER=-1 THERE HAS BEEN NO COMPUTATION.               DSG1 270
C           (2)   Z CAN HAVE THE SAME STORAGE ALLOCATION AS X OR Y.  IF DSG1 280
C                 X OR Y IS DISTINCT FROM Z, THEN IT IS NOT DESTROYED.  DSG1 290
C                                                                       DSG1 300
C        SUBROUTINES AND SUBPROGRAMS REQUIRED                           DSG1 310
C           NONE                                                        DSG1 320
C                                                                       DSG1 330
C        METHOD                                                         DSG1 340
C           EXCEPT AT THE ENDPOINTS X(1) AND X(NDIM), EACH SMOOTHED     DSG1 350
C           VALUE Z(I) IS OBTAINED BY EVALUATING AT X(I) THE LEAST-     DSG1 360
C           SQUARES POLYNOMIAL OF DEGREE 1 RELEVANT TO THE 3 SUCCESSIVE DSG1 370
C           POINTS (X(I+K),Y(I+K)) K = -1,0,1.(SEE HILDEBRAND, F.B.,    DSG1 380
C           INTRODUCTION  TO NUMERICAL ANALYSIS, MC GRAW-HILL, NEW YORK/DSG1 390
C           TORONTO/LONDON, 1956, PP.258-311.)                          DSG1 400
C                                                                       DSG1 410
C     ..................................................................DSG1 420
C                                                                       DSG1 430
      SUBROUTINE DSG13(X,Y,Z,NDIM,IER)                                  DSG1 440
C                                                                       DSG1 450
C                                                                       DSG1 460
      DIMENSION X(1),Y(1),Z(1)                                          DSG1 470
      DOUBLE PRECISION X,Y,Z,XM,YM,T1,T2,T3,H                           DSG1 480
C                                                                       DSG1 490
C        TEST OF DIMENSION                                              DSG1 500
      IF(NDIM-3)7,1,1                                                   DSG1 510
C                                                                       DSG1 520
C        START LOOP                                                     DSG1 530
    1 DO 6 I=3,NDIM                                                     DSG1 540
      XM=.33333333333333333D0*(X(I-2)+X(I-1)+X(I))                      DSG1 550
      YM=.33333333333333333D0*(Y(I-2)+Y(I-1)+Y(I))                      DSG1 560
      T1=X(I-2)-XM                                                      DSG1 570
      T2=X(I-1)-XM                                                      DSG1 580
      T3=X(I)-XM                                                        DSG1 590
      XM=T1*T1+T2*T2+T3*T3                                              DSG1 600
      IF(XM)3,3,2                                                       DSG1 610
    2 XM=(T1*(Y(I-2)-YM)+T2*(Y(I-1)-YM)+T3*(Y(I)-YM))/XM                DSG1 620
C                                                                       DSG1 630
C        CHECK FIRST POINT                                              DSG1 640
    3 IF(I-3)4,4,5                                                      DSG1 650
    4 H=XM*T1+YM                                                        DSG1 660
    5 Z(I-2)=H                                                          DSG1 670
    6 H=XM*T2+YM                                                        DSG1 680
C        END OF LOOP                                                    DSG1 690
C                                                                       DSG1 700
C        UPDATE LAST TWO COMPONENTS                                     DSG1 710
      Z(NDIM-1)=H                                                       DSG1 720
      Z(NDIM)=XM*T3+YM                                                  DSG1 730
      IER=0                                                             DSG1 740
      RETURN                                                            DSG1 750
C                                                                       DSG1 760
C        ERROR EXIT IN CASE NDIM IS LESS THAN 3                         DSG1 770
    7 IER=-1                                                            DSG1 780
      RETURN                                                            DSG1 790
      END                                                               DSG1 800
