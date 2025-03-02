C                                                                       SMO   10
C     ..................................................................SMO   20
C                                                                       SMO   30
C        SUBROUTINE SMO                                                 SMO   40
C                                                                       SMO   50
C        PURPOSE                                                        SMO   60
C           TO SMOOTH OR FILTER SERIES A BY WEIGHTS W.                  SMO   70
C                                                                       SMO   80
C        USAGE                                                          SMO   90
C           CALL SMO (A,N,W,M,L,R)                                      SMO  100
C                                                                       SMO  110
C        DESCRIPTION OF PARAMETERS                                      SMO  120
C           A - INPUT VECTOR OF LENGTH N CONTAINING TIME SERIES DATA.   SMO  130
C           N - LENGTH OF SERIES A.                                     SMO  140
C           W - INPUT VECTOR OF LENGTH M CONTAINING WEIGHTS.            SMO  150
C           M - NUMBER OF ITEMS IN WEIGHT VECTOR.  M MUST BE AN ODD     SMO  160
C               INTEGER.  (IF M IS AN EVEN INTEGER, ANY FRACTION        SMO  170
C               RESULTING FROM THE CALCULATION OF (L*(M-1))/2 IN (1)    SMO  180
C               AND (2) BELOW WILL BE TRUNCATED.)                       SMO  190
C           L - SELECTION INTEGER.  FOR EXAMPLE, L=12 MEANS THAT WEIGHTSSMO  200
C               ARE APPLIED TO EVERY 12-TH ITEM OF A.  L=1 APPLIES      SMO  210
C               WEIGHTS TO SUCCESSIVE ITEMS OF A.  FOR MONTHLY DATA,    SMO  220
C               L=12 GIVES YEAR-TO-YEAR AVERAGES AND L=1 GIVES MONTH-TO-SMO  230
C               MONTH AVERAGES.                                         SMO  240
C           R - OUTPUT VECTOR OF LENGTH N.  FROM IL TO IH ELEMENTS OF   SMO  250
C               THE VECTOR R ARE FILLED WITH THE SMOOTHED SERIES AND    SMO  260
C               OTHER ELEMENTS WITH ZERO, WHERE                         SMO  270
C                    IL=(L*(M-1))/2+1  ................ (1)             SMO  280
C                    IH=N-(L*(M-1))/2  ................ (2)             SMO  290
C                                                                       SMO  300
C        REMARKS                                                        SMO  310
C           N MUST BE GREATER THAN OR EQUAL TO THE PRODUCT OF L*M.      SMO  320
C                                                                       SMO  330
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  SMO  340
C           NONE                                                        SMO  350
C                                                                       SMO  360
C        METHOD                                                         SMO  370
C           REFER TO THE ARTICLE 'FORTRAN SUBROUTINES FOR TIME SERIES   SMO  380
C           ANALYSIS', BY J. R. HEALY AND B. P. BOGERT, COMMUNICATIONS  SMO  390
C           OF ACM, V.6, NO.1, JANUARY, 1963.                           SMO  400
C                                                                       SMO  410
C     ..................................................................SMO  420
C                                                                       SMO  430
      SUBROUTINE SMO (A,N,W,M,L,R)                                      SMO  440
      DIMENSION A(1),W(1),R(1)                                          SMO  450
C                                                                       SMO  460
C     INITIALIZATION                                                    SMO  470
C                                                                       SMO  480
      DO 110 I=1,N                                                      SMO  490
  110 R(I)=0.0                                                          SMO  500
      IL=(L*(M-1))/2+1                                                  SMO  510
      IH=N-(L*(M-1))/2                                                  SMO  520
C                                                                       SMO  530
C     SMOOTH SERIES A BY WEIGHTS W                                      SMO  540
C                                                                       SMO  550
      DO 120 I=IL,IH                                                    SMO  560
      K=I-IL+1                                                          SMO  570
      DO 120 J=1,M                                                      SMO  580
      IP=(J*L)-L+K                                                      SMO  590
  120 R(I)=R(I)+A(IP)*W(J)                                              SMO  600
      RETURN                                                            SMO  610
      END                                                               SMO  620
