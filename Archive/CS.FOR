C                                                                       CS    10
C     ..................................................................CS    20
C                                                                       CS    30
C        SUBROUTINE CS                                                  CS    40
C                                                                       CS    50
C        PURPOSE                                                        CS    60
C           COMPUTES THE FRESNEL INTEGRALS.                             CS    70
C                                                                       CS    80
C        USAGE                                                          CS    90
C           CALL CS (C,S,X)                                             CS   100
C                                                                       CS   110
C        DESCRIPTION OF PARAMETERS                                      CS   120
C           C     - THE RESULTANT VALUE C(X).                           CS   130
C           S     - THE RESULTANT VALUE S(X).                           CS   140
C           X     - THE ARGUMENT OF FRESNEL INTEGRALS                   CS   150
C                   IF X IS NEGATIVE, THE ABSOLUTE VALUE IS USED.       CS   160
C                                                                       CS   170
C        REMARKS                                                        CS   180
C           THE ARGUMENT VALUE X REMAINS UNCHANGED.                     CS   190
C                                                                       CS   200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  CS   210
C           NONE                                                        CS   220
C                                                                       CS   230
C        METHOD                                                         CS   240
C           DEFINITION                                                  CS   250
C           C(X)=INTEGRAL(COS(T)/SQRT(2*LI*T) SUMMED OVER T FROM 0 TO X)CS   260
C           S(X)=INTEGRAL(SIN(T)/SQRT(I*LI*T) SUMMED OVER T FROM 0 TO X)CS   270
C           EVALUATION                                                  CS   280
C           USING DIFFERENT APPROXIMATIONS FOR X LESS THAN 4 AND X      CS   290
C           GREATER THAN 4.                                             CS   300
C           REFERENCE                                                   CS   310
C           'COMPUTATION OF FRESNEL INTEGRALS' BY BOERSMA,              CS   320
C           MATHEMATICAL TABLES AND OTHER AIDS TO COMPUTATION, VOL. 14, CS   330
C           1960, NO. 72, P. 380.                                       CS   340
C                                                                       CS   350
C     ..................................................................CS   360
C                                                                       CS   370
      SUBROUTINE CS(C,S,X)                                              CS   380
      Z=ABS(X)                                                          CS   390
      IF(Z-4.)1,1,2                                                     CS   400
    1 C=SQRT(Z)                                                         CS   410
      S=Z*C                                                             CS   420
      Z=(4.-Z)*(4.+Z)                                                   CS   430
      C=C*((((((5.100785E-11*Z+5.244297E-9)*Z+5.451182E-7)*Z            CS   440
     1+3.273308E-5)*Z+1.020418E-3)*Z+1.102544E-2)*Z+1.840965E-1)        CS   450
      S=S*(((((6.677681E-10*Z+5.883158E-8)*Z+5.051141E-6)*Z             CS   460
     1+2.441816E-4)*Z+6.121320E-3)*Z+8.026490E-2)                       CS   470
      RETURN                                                            CS   480
    2 D=COS(Z)                                                          CS   490
      S=SIN(Z)                                                          CS   500
      Z=4./Z                                                            CS   510
      A=(((((((8.768258E-4*Z-4.169289E-3)*Z+7.970943E-3)*Z-6.792801E-3) CS   520
     1*Z-3.095341E-4)*Z+5.972151E-3)*Z-1.606428E-5)*Z-2.493322E-2)*Z    CS   530
     2-4.444091E-9                                                      CS   540
      B=((((((-6.633926E-4*Z+3.401409E-3)*Z-7.271690E-3)*Z+7.428246E-3) CS   550
     1*Z-4.027145E-4)*Z-9.314910E-3)*Z-1.207998E-6)*Z+1.994711E-1       CS   560
      Z=SQRT(Z)                                                         CS   570
      C=0.5+Z*(D*A+S*B)                                                 CS   580
      S=0.5+Z*(S*A-D*B)                                                 CS   590
      RETURN                                                            CS   600
      END                                                               CS   610
