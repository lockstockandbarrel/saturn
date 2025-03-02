C                                                                       SE35  10
C     ..................................................................SE35  20
C                                                                       SE35  30
C        SUBROUTINE SE35                                                SE35  40
C                                                                       SE35  50
C        PURPOSE                                                        SE35  60
C           TO COMPUTE A VECTOR OF SMOOTHED FUNCTION VALUES GIVEN A     SE35  70
C           VECTOR OF FUNCTION VALUES WHOSE ENTRIES CORRESPOND TO       SE35  80
C           EQUIDISTANTLY SPACED ARGUMENT VALUES.                       SE35  90
C                                                                       SE35 100
C        USAGE                                                          SE35 110
C           CALL SE35(Y,Z,NDIM,IER)                                     SE35 120
C                                                                       SE35 130
C        DESCRIPTION OF PARAMETERS                                      SE35 140
C           Y     -  GIVEN VECTOR OF FUNCTION VALUES (DIMENSION NDIM)   SE35 150
C           Z     -  RESULTING VECTOR OF SMOOTHED FUNCTION VALUES       SE35 160
C                    (DIMENSION NDIM)                                   SE35 170
C           NDIM  -  DIMENSION OF VECTORS Y AND Z                       SE35 180
C           IER   -  RESULTING ERROR PARAMETER                          SE35 190
C                    IER = -1  - NDIM IS LESS THAN 5                    SE35 200
C                    IER =  0  - NO ERROR                               SE35 210
C                                                                       SE35 220
C        REMARKS                                                        SE35 230
C           (1)  IF IER=-1 THERE HAS BEEN NO COMPUTATION.               SE35 240
C           (2)   Z CAN HAVE THE SAME STORAGE ALLOCATION AS Y.  IF Y IS SE35 250
C                 DISTINCT FROM Z, THEN IT IS NOT DESTROYED.            SE35 260
C                                                                       SE35 270
C        SUBROUTINE AND FUNCTION SUBPROGRAMS REQUIRED                   SE35 280
C           NONE                                                        SE35 290
C                                                                       SE35 300
C        METHOD                                                         SE35 310
C           IF X IS THE (SUPPRESSED) VECTOR OF ARGUMENT VALUES, THEN    SE35 320
C           EXCEPT AT THE POINTS X(1),X(2),X(NDIM-1) AND X(NDIM), EACH  SE35 330
C           SMOOTHED VALUE Z(I) IS OBTAINED BY EVALUATING AT X(I) THE   SE35 340
C           LEAST-SQUARES POLYNOMIAL OF DEGREE 3 RELEVANT TO THE 5      SE35 350
C           SUCCESSIVE POINTS (X(I+K),Y(I+K)) K = -2,-1,...,2.  (SEE    SE35 360
C           HILDEBRAND, F.B., INTRODUCTION TO NUMERICAL ANALYSIS,       SE35 370
C           MC GRAW-HILL, NEW YORK/TORONTO/LONDON, 1956, PP. 295-302.)  SE35 380
C                                                                       SE35 390
C     ..................................................................SE35 400
C                                                                       SE35 410
      SUBROUTINE SE35(Y,Z,NDIM,IER)                                     SE35 420
C                                                                       SE35 430
C                                                                       SE35 440
      DIMENSION Y(1),Z(1)                                               SE35 450
C                                                                       SE35 460
C        TEST OF DIMENSION                                              SE35 470
      IF(NDIM-5)4,1,1                                                   SE35 480
C                                                                       SE35 490
C        PREPARE LOOP                                                   SE35 500
    1 B=Y(1)                                                            SE35 510
      C=Y(2)                                                            SE35 520
C                                                                       SE35 530
C        START LOOP                                                     SE35 540
      DO 3 I=5,NDIM                                                     SE35 550
      A=B                                                               SE35 560
      B=C                                                               SE35 570
      C=Y(I-2)                                                          SE35 580
C                                                                       SE35 590
C        GENERATE FOURTH CENTRAL DIFFERENCE                             SE35 600
      D=C-B-Y(I-1)                                                      SE35 610
      D=D+D+C                                                           SE35 620
      D=D+D+A+Y(I)                                                      SE35 630
C                                                                       SE35 640
C        CHECK FIRST TWO COMPONENTS                                     SE35 650
      IF(I-5)2,2,3                                                      SE35 660
    2 Z(1)=A-.01428571*D                                                SE35 670
      Z(2)=B+.05714286*D                                                SE35 680
    3 Z(I-2)=C-.08571429*D                                              SE35 690
C        END OF LOOP                                                    SE35 700
C                                                                       SE35 710
C        UPDATE LAST TWO COMPONENTS                                     SE35 720
      Z(NDIM-1)=Y(NDIM-1)+.05714286*D                                   SE35 730
      Z(NDIM)=Y(NDIM)-.01428571*D                                       SE35 740
      IER=0                                                             SE35 750
      RETURN                                                            SE35 760
C                                                                       SE35 770
C        ERROR EXIT IN CASE NDIM IS LESS THAN 5                         SE35 780
    4 IER=-1                                                            SE35 790
      RETURN                                                            SE35 800
      END                                                               SE35 810
