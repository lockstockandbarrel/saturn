C                                                                       RK2   10
C     ..................................................................RK2   20
C                                                                       RK2   30
C        SUBROUTINE RK2                                                 RK2   40
C                                                                       RK2   50
C        PURPOSE                                                        RK2   60
C           INTEGRATES A FIRST ORDER DIFFERENTIAL EQUATION              RK2   70
C           DY/DX=FUN(X,Y) AND PRODUCES A TABLE OF INTEGRATED VALUES    RK2   80
C                                                                       RK2   90
C        USAGE                                                          RK2  100
C           CALL RK2(FUN,H,XI,YI,K,N,VEC)                               RK2  110
C                                                                       RK2  120
C        DESCRIPTION OF PARAMETERS                                      RK2  130
C           FUN-USER-SUPPLIED FUNCTION SUBPROGRAM WITH ARGUMENTS X,Y    RK2  140
C               WHICH GIVES DY/DX                                       RK2  150
C           H  -STEP SIZE                                               RK2  160
C           XI -INITIAL VALUE OF X                                      RK2  170
C           YI -INITIAL VALUE OF Y WHERE YI=Y(XI)                       RK2  180
C           K  -THE INTERVAL AT WHICH COMPUTED VALUES ARE TO BE STORED  RK2  190
C           N  -THE NUMBER OF VALUES TO BE STORED                       RK2  200
C           VEC-THE RESULTANT VECTOR OF LENGTH N IN WHICH COMPUTED      RK2  210
C               VALUES OF Y ARE TO BE STORED                            RK2  220
C                                                                       RK2  230
C        REMARKS                                                        RK2  240
C           NONE                                                        RK2  250
C                                                                       RK2  260
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  RK2  270
C           FUN - USER-SUPPLIED FUNCTION SUBPROGRAM FOR DY/DX           RK2  280
C           CALLING PROGRAM MUST HAVE FORTRAN EXTERNAL STATEMENT        RK2  290
C           CONTAINING NAMES OF FUNCTION SUBPROGRAMS LISTED IN CALL TO  RK2  300
C           RK2                                                         RK2  310
C                                                                       RK2  320
C        METHOD                                                         RK2  330
C           FOURTH ORDER RUNGE-KUTTA INTEGRATION ON A RECURSIVE BASIS ASRK2  340
C           SHOWN IN F.B. HILDEBRAND, 'INTRODUCTION TO NUMERICAL        RK2  350
C           ANALYSIS', MCGRAW-HILL, NEW YORK, 1956                      RK2  360
C                                                                       RK2  370
C     ..................................................................RK2  380
C                                                                       RK2  390
      SUBROUTINE RK2(FUN,H,XI,YI,K,N,VEC)                               RK2  400
C                                                                       RK2  410
C        ...............................................................RK2  420
C                                                                       RK2  430
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE  RK2  440
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION      RK2  450
C        STATEMENT WHICH FOLLOWS.                                       RK2  460
C                                                                       RK2  470
C     DOUBLE PRECISION H,XI,YI,VEC,H2,Y,X,T1,T2,T3,T4,FUN               RK2  480
C                                                                       RK2  490
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS    RK2  500
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS      RK2  510
C        ROUTINE.                                                       RK2  520
C                                                                       RK2  530
C        USER FUNCTION SUBPROGRAM, FUN, MUST BE IN DOUBLE PRECISION.    RK2  540
C                                                                       RK2  550
C        ...............................................................RK2  560
C                                                                       RK2  570
      DIMENSION VEC(1)                                                  RK2  580
      H2=H/2.                                                           RK2  590
      Y=YI                                                              RK2  600
      X=XI                                                              RK2  610
      DO 2 I=1,N                                                        RK2  620
      DO 1 J=1,K                                                        RK2  630
      T1=H*FUN(X,Y)                                                     RK2  640
      T2=H*FUN(X+H2,Y+T1/2.)                                            RK2  650
      T3=H*FUN(X+H2,Y+T2/2.)                                            RK2  660
      T4=H*FUN(X+H,Y+T3)                                                RK2  670
      Y= Y+(T1+2.*T2+2.*T3+T4)/6.                                       RK2  680
    1 X=X+H                                                             RK2  690
    2 VEC(I)=Y                                                          RK2  700
      RETURN                                                            RK2  710
      END                                                               RK2  720
