C                                                                       I0    10
C     ..................................................................I0    20
C                                                                       I0    30
C        SUBROUTINE I0                                                  I0    40
C                                                                       I0    50
C        PURPOSE                                                        I0    60
C            COMPUTE THE MODIFIED BESSEL FUNCTION I OF ORDER ZERO       I0    70
C                                                                       I0    80
C        USAGE                                                          I0    90
C            CALL I0(X,RI0)                                             I0   100
C                                                                       I0   110
C        DESCRIPTION OF PARAMETERS                                      I0   120
C            X    -GIVEN ARGUMENT OF THE BESSEL FUNCTION I OF ORDER 0   I0   130
C            RI0  -RESULTANT VALUE OF THE BESSEL FUNCTION I OF ORDER 0  I0   140
C                                                                       I0   150
C        REMARKS                                                        I0   160
C            LARGE VALUES OF THE ARGUMENT MAY CAUSE OVERFLOW IN THE     I0   170
C            BUILTIN EXP-FUNCTION                                       I0   180
C                                                                       I0   190
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  I0   200
C           NONE                                                        I0   210
C                                                                       I0   220
C        METHOD                                                         I0   230
C           POLYNOMIAL APPROXIMATIONS GIVEN BY E.E. ALLEN ARE USED FOR  I0   240
C           CALCULATION.                                                I0   250
C           FOR REFERENCE SEE                                           I0   260
C           M. ABRAMOWITZ AND I.A. STEGUN,'HANDBOOK OF MATHEMATICAL     I0   270
C           FUNCTIONS', U.S. DEPARTMENT OF COMMERCE, NATIONAL BUREAU OF I0   280
C           STANDARDS APPLIED MATHEMATICS SERIES, 1966, P.378.          I0   290
C                                                                       I0   300
C     ..................................................................I0   310
C                                                                       I0   320
      SUBROUTINE I0(X,RI0)                                              I0   330
      RI0=ABS(X)                                                        I0   340
      IF(RI0-3.75)1,1,2                                                 I0   350
    1 Z=X*X*7.111111E-2                                                 I0   360
      RI0=((((( 4.5813E-3*Z+3.60768E-2)*Z+2.659732E-1)*Z+1.206749E0)*Z  I0   370
     1+3.089942E0)*Z+3.515623E0)*Z+1.                                   I0   380
      RETURN                                                            I0   390
    2 Z=3.75/RI0                                                        I0   400
      RI0= EXP(RI0)/SQRT(RI0)*((((((((3.92377E-3*Z-1.647633E-2)*Z       I0   410
     1+2.635537E-2)*Z-2.057706E-2)*Z+9.16281E-3)*Z-1.57565E-3)*Z        I0   420
     2+2.25319E-3)*Z+1.328592E-2)*Z+3.989423E-1)                        I0   430
      RETURN                                                            I0   440
      END                                                               I0   450
