C                                                                       DSE3  10
C     ..................................................................DSE3  20
C                                                                       DSE3  30
C        SUBROUTINE DSE35                                               DSE3  40
C                                                                       DSE3  50
C        PURPOSE                                                        DSE3  60
C           TO COMPUTE A VECTOR OF SMOOTHED FUNCTION VALUES GIVEN A     DSE3  70
C           VECTOR OF FUNCTION VALUES WHOSE ENTRIES CORRESPOND TO       DSE3  80
C           EQUIDISTANTLY SPACED ARGUMENT VALUES.                       DSE3  90
C                                                                       DSE3 100
C        USAGE                                                          DSE3 110
C           CALL DSE35(Y,Z,NDIM,IER)                                    DSE3 120
C                                                                       DSE3 130
C        DESCRIPTION OF PARAMETERS                                      DSE3 140
C           Y     -  GIVEN VECTOR OF DOUBLE PRECISION FUNCTION VALUES   DSE3 150
C                    (DIMENSION NDIM)                                   DSE3 160
C           Z     -  RESULTING VECTOR OF DOUBLE PRECISION SMOOTHED      DSE3 170
C                    FUNCTION VALUES (DIMENSION NDIM)                   DSE3 180
C           NDIM  -  DIMENSION OF VECTORS Y AND Z                       DSE3 190
C           IER   -  RESULTING ERROR PARAMETER                          DSE3 200
C                    IER = -1  - NDIM IS LESS THAN 5                    DSE3 210
C                    IER =  0  - NO ERROR                               DSE3 220
C                                                                       DSE3 230
C        REMARKS                                                        DSE3 240
C           (1)  IF IER=-1 THERE HAS BEEN NO COMPUTATION.               DSE3 250
C           (2)   Z CAN HAVE THE SAME STORAGE ALLOCATION AS Y.  IF Y IS DSE3 260
C                 DISTINCT FROM Z, THEN IT IS NOT DESTROYED.            DSE3 270
C                                                                       DSE3 280
C        SUBROUTINE AND FUNCTION SUBPROGRAMS REQUIRED                   DSE3 290
C           NONE                                                        DSE3 300
C                                                                       DSE3 310
C        METHOD                                                         DSE3 320
C           IF X IS THE (SUPPRESSED) VECTOR OF ARGUMENT VALUES, THEN    DSE3 330
C           EXCEPT AT THE POINTS X(1),X(2),X(NDIM-1) AND X(NDIM), EACH  DSE3 340
C           SMOOTHED VALUE Z(I) IS OBTAINED BY EVALUATING AT X(I) THE   DSE3 350
C           LEAST-SQUARES POLYNOMIAL OF DEGREE 3 RELEVANT TO THE 5      DSE3 360
C           SUCCESSIVE POINTS (X(I+K),Y(I+K)) K = -2,-1,...,2.  (SEE    DSE3 370
C           HILDEBRAND, F.B., INTRODUCTION TO NUMERICAL ANALYSIS,       DSE3 380
C           MC GRAW-HILL, NEW YORK/TORONTO/LONDON, 1956, PP. 295-302.)  DSE3 390
C                                                                       DSE3 400
C     ..................................................................DSE3 410
C                                                                       DSE3 420
      SUBROUTINE DSE35(Y,Z,NDIM,IER)                                    DSE3 430
C                                                                       DSE3 440
C                                                                       DSE3 450
      DIMENSION Y(1),Z(1)                                               DSE3 460
      DOUBLE PRECISION Y,Z,A,B,C,D                                      DSE3 470
C                                                                       DSE3 480
C        TEST OF DIMENSION                                              DSE3 490
      IF(NDIM-5)4,1,1                                                   DSE3 500
C                                                                       DSE3 510
C        PREPARE LOOP                                                   DSE3 520
    1 B=Y(1)                                                            DSE3 530
      C=Y(2)                                                            DSE3 540
C                                                                       DSE3 550
C        START LOOP                                                     DSE3 560
      DO 3 I=5,NDIM                                                     DSE3 570
      A=B                                                               DSE3 580
      B=C                                                               DSE3 590
      C=Y(I-2)                                                          DSE3 600
C                                                                       DSE3 610
C        GENERATE FOURTH CENTRAL DIFFERENCE                             DSE3 620
      D=C-B-Y(I-1)                                                      DSE3 630
      D=D+D+C                                                           DSE3 640
      D=D+D+A+Y(I)                                                      DSE3 650
C                                                                       DSE3 660
C        CHECK FIRST TWO COMPONENTS                                     DSE3 670
      IF(I-5)2,2,3                                                      DSE3 680
    2 Z(1)=A-.014285714285714286D0*D                                    DSE3 690
      Z(2)=B+.057142857142857143D0*D                                    DSE3 700
    3 Z(I-2)=C-.08571428571428571D0*D                                   DSE3 710
C        END OF LOOP                                                    DSE3 720
C                                                                       DSE3 730
C        UPDATE LAST TWO COMPONENTS                                     DSE3 740
      Z(NDIM-1)=Y(NDIM-1)+.057142857142857143D0*D                       DSE3 750
      Z(NDIM)=Y(NDIM)-.014285714285714286D0*D                           DSE3 760
      IER=0                                                             DSE3 770
      RETURN                                                            DSE3 780
C                                                                       DSE3 790
C        ERROR EXIT IN CASE NDIM IS LESS THAN 5                         DSE3 800
    4 IER=-1                                                            DSE3 810
      RETURN                                                            DSE3 820
      END                                                               DSE3 830
