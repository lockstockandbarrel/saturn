C                                                                       TIE   10
C     ..................................................................TIE   20
C                                                                       TIE   30
C        SUBROUTINE TIE                                                 TIE   40
C                                                                       TIE   50
C        PURPOSE                                                        TIE   60
C           CALCULATE CORRECTION FACTOR DUE TO TIES                     TIE   70
C                                                                       TIE   80
C        USAGE                                                          TIE   90
C           CALL TIE(R,N,KT,T)                                          TIE  100
C                                                                       TIE  110
C        DESCRIPTION OF PARAMETERS                                      TIE  120
C           R  - INPUT VECTOR OF RANKS OF LENGTH N CONTAINING VALUES    TIE  130
C                1 TO N                                                 TIE  140
C           N  - NUMBER OF RANKED VALUES                                TIE  150
C           KT - INPUT CODE FOR CALCULATION OF CORRECTION FACTOR        TIE  160
C                      1   SOLVE EQUATION 1                             TIE  170
C                      2   SOLVE EQUATION 2                             TIE  180
C           T  - CORRECTION FACTOR (OUTPUT)                             TIE  190
C                    EQUATION 1   T=SUM(CT**3-CT)/12                    TIE  200
C                    EQUATION 2   T=SUM(CT*(CT-1)/2)                    TIE  210
C                  WHERE CT IS THE NUMBER OF OBSERVATIONS TIED FOR A    TIE  220
C                        GIVEN RANK                                     TIE  230
C                                                                       TIE  240
C        REMARKS                                                        TIE  250
C           NONE                                                        TIE  260
C                                                                       TIE  270
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  TIE  280
C           NONE                                                        TIE  290
C                                                                       TIE  300
C        METHOD                                                         TIE  310
C           VECTOR IS SEARCHED FOR SUCCESSIVELY LARGER RANKS. TIES ARE  TIE  320
C           COUNTED AND CORRECTION FACTOR 1 OR 2 SUMMED.                TIE  330
C                                                                       TIE  340
C     ..................................................................TIE  350
C                                                                       TIE  360
      SUBROUTINE TIE(R,N,KT,T)                                          TIE  370
      DIMENSION R(1)                                                    TIE  380
C                                                                       TIE  390
C        INITIALIZATION                                                 TIE  400
C                                                                       TIE  410
      T=0.0                                                             TIE  420
      Y=0.0                                                             TIE  430
    5 X=1.0E38                                                          TIE  440
      IND=0                                                             TIE  450
C                                                                       TIE  460
C        FIND NEXT LARGEST RANK                                         TIE  470
C                                                                       TIE  480
      DO 30 I=1,N                                                       TIE  490
      IF(R(I)-Y) 30,30,10                                               TIE  500
   10 IF(R(I)-X) 20,30,30                                               TIE  510
   20 X=R(I)                                                            TIE  520
      IND=IND+1                                                         TIE  530
   30 CONTINUE                                                          TIE  540
C                                                                       TIE  550
C        IF ALL RANKS HAVE BEEN TESTED, RETURN                          TIE  560
C                                                                       TIE  570
      IF(IND) 90,90,40                                                  TIE  580
   40 Y=X                                                               TIE  590
      CT=0.0                                                            TIE  600
C                                                                       TIE  610
C        COUNT TIES                                                     TIE  620
C                                                                       TIE  630
      DO 60 I=1,N                                                       TIE  640
      IF(R(I)-X) 60,50,60                                               TIE  650
   50 CT=CT+1.0                                                         TIE  660
   60 CONTINUE                                                          TIE  670
C                                                                       TIE  680
C        CALCULATE CORRECTION FACTOR                                    TIE  690
C                                                                       TIE  700
      IF(CT) 70,5,70                                                    TIE  710
   70 IF(KT-1) 75,80,75                                                 TIE  720
   75 T=T+CT*(CT-1.)/2.0                                                TIE  730
      GO TO 5                                                           TIE  740
   80 T=T+(CT*CT*CT-CT)/12.0                                            TIE  750
      GO TO 5                                                           TIE  760
   90 RETURN                                                            TIE  770
      END                                                               TIE  780
