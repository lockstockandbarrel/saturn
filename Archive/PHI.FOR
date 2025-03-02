C                                                                       PHI   10
C     ..................................................................PHI   20
C                                                                       PHI   30
C        SUBROUTINE PHI                                                 PHI   40
C                                                                       PHI   50
C        PURPOSE                                                        PHI   60
C           TO COMPUTE THE PHI COEFFICIENT BETWEEN TWO VARIABLES  WHICH PHI   70
C           ARE DICHOTOMOUS.                                            PHI   80
C                                                                       PHI   90
C        USAGE                                                          PHI  100
C           CALL PHI (N,U,V,HU,HV,P,CH,XP,IE)                           PHI  110
C                                                                       PHI  120
C        DESCRIPTION OF PARAMETERS                                      PHI  130
C           N  - NUMBER OF OBSERVATIONS                                 PHI  140
C           U  - INPUT VECTOR OF LENGTH N CONTAINING THE FIRST DICHOTO- PHI  150
C                MOUS VARIABLE                                          PHI  160
C           V  - INPUT VECTOR OF LENGTH N CONTAINING THE SECOND DICHOTO-PHI  170
C                MOUS VARIABLE                                          PHI  180
C           HU - INPUT NUMERICAL CODE WHICH INDICATES THE HIGHER        PHI  190
C                CATEGORY OF THE FIRST VARIABLE.  ANY OBSERVATION IN    PHI  200
C                VECTOR U WHICH HAS A VALUE EQUAL TO OR GREATER THAN HU PHI  210
C                WILL BE CLASSIFIED IN THE HIGHER CATEGORY.             PHI  220
C           HV - INPUT NUMERICAL CODE FOR VECTOR V, SIMILAR TO HU       PHI  230
C           P  - PHI COEFFICIENT COMPUTED                               PHI  240
C           CH - CHI-SQUARE COMPUTED AS A FUNCTION OF PHI COEFFICIENT   PHI  250
C                (DEGREES OF FREEDOM FOR CHI-SQUARE = 1)                PHI  260
C           XP - COMPUTED VALUE OF THE MAXIMAL PHI COEFFICIENT THAT     PHI  270
C                CAN BE ATTAINED IN THE PROBLEM                         PHI  280
C           IE - IF IE IS NON-ZERO, SOME CELL IN THE 2 BY 2 TABLE IS    PHI  290
C                NULL.  IF SO, P, CH, AND XP ARE SET TO 10**75.         PHI  300
C                                                                       PHI  310
C        REMARKS                                                        PHI  320
C           VARIABLES U AND V MUST BE SPECIFIED NUMERIC.                PHI  330
C           THE PHI COEFFICIENT IS A SPECIAL CASE OF THE                PHI  340
C           PEARSON PRODUCT-MOMENT CORRELATION WHEN BOTH VARIABLES ARE  PHI  350
C           BINARY.                                                     PHI  360
C                                                                       PHI  370
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PHI  380
C           NONE                                                        PHI  390
C                                                                       PHI  400
C        METHOD                                                         PHI  410
C           REFER TO P. HORST, 'PYSCHOLOGICAL MEASUREMENT AND           PHI  420
C           PREDICTION', P. 94 (WADSWORTH, 1966).                       PHI  430
C                                                                       PHI  440
C     ..................................................................PHI  450
C                                                                       PHI  460
      SUBROUTINE PHI (N,U,V,HU,HV,P,CH,XP,IE)                           PHI  470
C                                                                       PHI  480
      DIMENSION U(1),V(1)                                               PHI  490
C                                                                       PHI  500
C        CONSTRUCT A 2X2 CONTINGENCY TABLE                              PHI  510
C                                                                       PHI  520
      IE=0                                                              PHI  530
      A=0.0                                                             PHI  540
      B=0.0                                                             PHI  550
      C=0.0                                                             PHI  560
      D=0.0                                                             PHI  570
C                                                                       PHI  580
      DO 40 I=1,N                                                       PHI  590
      IF(U(I)-HU) 10,25,25                                              PHI  600
   10 IF(V(I)-HV) 15,20,20                                              PHI  610
   15 D=D+1.0                                                           PHI  620
      GO TO 40                                                          PHI  630
   20 B=B+1.0                                                           PHI  640
      GO TO 40                                                          PHI  650
   25 IF(V(I)-HV) 30,35,35                                              PHI  660
   30 C=C+1.0                                                           PHI  670
      GO TO 40                                                          PHI  680
   35 A=A+1.0                                                           PHI  690
   40 CONTINUE                                                          PHI  700
      IF(A) 100,100,41                                                  PHI  710
   41 IF(B) 100,100,42                                                  PHI  720
   42 IF(C) 100,100,43                                                  PHI  730
   43 IF(D) 100,100,44                                                  PHI  740
C                                                                       PHI  750
C        COMPUTE THE PHI COEFFICIENT                                    PHI  760
C                                                                       PHI  770
   44 P=(A*D-B*C)/ SQRT((A+B)*(C+D)*(A+C)*(B+D))                        PHI  780
C                                                                       PHI  790
C        COMPUTE CHI-SQURE                                              PHI  800
C                                                                       PHI  810
      T=N                                                               PHI  820
      CH=T*P*P                                                          PHI  830
C                                                                       PHI  840
C        COMPUTE THE MAXIMAL PHI COEFFICIENT                            PHI  850
C                                                                       PHI  860
      P1=(A+C)/T                                                        PHI  870
      P2=(B+D)/T                                                        PHI  880
      P3=(A+B)/T                                                        PHI  890
      P4=(C+D)/T                                                        PHI  900
      IF(P1-P2) 75, 45, 45                                              PHI  910
   45 IF(P3-P4) 65, 50, 50                                              PHI  920
   50 IF(P1-P3) 60, 55, 55                                              PHI  930
   55 XP=SQRT((P3/P4)*(P2/P1))                                          PHI  940
      GO TO 95                                                          PHI  950
   60 XP=SQRT((P1/P2)*(P4/P3))                                          PHI  960
      GO TO 95                                                          PHI  970
   65 IF(P1-P4) 70, 55, 55                                              PHI  980
   70 XP=SQRT((P2/P1)*(P3/P4))                                          PHI  990
      GO TO 95                                                          PHI 1000
   75 IF(P3-P4) 90, 80, 80                                              PHI 1010
   80 IF(P2-P3) 60, 85, 85                                              PHI 1020
   85 XP=SQRT((P4/P3)*(P1/P2))                                          PHI 1030
      GO TO 95                                                          PHI 1040
   90 IF(P2-P4) 70, 85, 85                                              PHI 1050
C                                                                       PHI 1060
   95 RETURN                                                            PHI 1070
  100 IE=1                                                              PHI 1080
      P=1.E75                                                           PHI 1090
      CH=1.E75                                                          PHI 1100
      XP=1.E75                                                          PHI 1110
      GO TO 95                                                          PHI 1120
      END                                                               PHI 1130
