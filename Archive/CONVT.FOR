C                                                                       CONV  10
C     ..................................................................CONV  20
C                                                                       CONV  30
C        SUBROUTINE CONVT                                               CONV  40
C                                                                       CONV  50
C        PURPOSE                                                        CONV  60
C           CONVERT NUMBERS FROM SINGLE PRECISION TO DOUBLE PRECISION   CONV  70
C           OR FROM DOUBLE PRECISION TO SINGLE PRECISION.               CONV  80
C                                                                       CONV  90
C        USAGE                                                          CONV 100
C           CALL CONVT (N,M,MODE,S,D,MS)                                CONV 110
C                                                                       CONV 120
C        DESCRIPTION OF PARAMETERS                                      CONV 130
C           N    - NUMBER OF ROWS IN MATRICES S AND D.                  CONV 140
C           M    - NUMBER OF COLUMNS IN MATRICES S AND D.               CONV 150
C           MODE - CODE INDICATING TYPE OF CONVERSION                   CONV 160
C                    1 - FROM SINGLE PRECISION TO DOUBLE PRECISION      CONV 170
C                    2 - FROM DOUBLE PRECISION TO SINGLE PRECISION      CONV 180
C           S    - IF MODE=1, THIS MATRIX CONTAINS SINGLE PRECISION     CONV 190
C                  NUMBERS AS INPUT.  IF MODE=2, IT CONTAINS SINGLE     CONV 200
C                  PRECISION NUMBERS AS OUTPUT.  THE SIZE OF MATRIX S   CONV 210
C                  IS N BY M.                                           CONV 220
C           D    - IF MODE=1, THIS MATRIX CONTAINS DOUBLE PRECISION     CONV 230
C                  NUMBERS AS OUTPUT.  IF MODE=2, IT CONTAINS DOUBLE    CONV 240
C                  PRECISION NUMBERS AS INPUT.  THE SIZE OF MATRIX D IS CONV 250
C                  N BY M.                                              CONV 260
C           MS   - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX          CONV 270
C                    0 - GENERAL                                        CONV 280
C                    1 - SYMMETRIC                                      CONV 290
C                    2 - DIAGONAL                                       CONV 300
C                                                                       CONV 310
C        REMARKS                                                        CONV 320
C           MATRIX D CANNOT BE IN THE SAME LOCATION AS MATRIX S.        CONV 330
C           MATRIX D MUST BE DEFINED BY A DOUBLE PRECISION STATEMENT IN CONV 340
C           THE CALLING PROGRAM.                                        CONV 350
C                                                                       CONV 360
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  CONV 370
C           NONE                                                        CONV 380
C                                                                       CONV 390
C        METHOD                                                         CONV 400
C           ACCORDING TO THE TYPE OF CONVERSION INDICATED IN MODE, THIS CONV 410
C           SUBROUTINE COPIES NUMBERS FROM MATRIX S TO MATRIX D OR FROM CONV 420
C           MATRIX D TO MATRIX S.                                       CONV 430
C                                                                       CONV 440
C     ..................................................................CONV 450
C                                                                       CONV 460
      SUBROUTINE CONVT (N,M,MODE,S,D,MS)                                CONV 470
      DIMENSION S(1),D(1)                                               CONV 480
      DOUBLE PRECISION D                                                CONV 490
C                                                                       CONV 500
C        FIND STORAGE MODE OF MATRIX AND NUMBER OF DATA POINTS          CONV 510
C                                                                       CONV 520
      IF(MS-1) 2, 4, 6                                                  CONV 530
    2 NM=N*M                                                            CONV 540
      GO TO 8                                                           CONV 550
    4 NM=((N+1)*N)/2                                                    CONV 560
      GO TO 8                                                           CONV 570
    6 NM=N                                                              CONV 580
C                                                                       CONV 590
C        TEST TYPE OF CONVERSION                                        CONV 600
C                                                                       CONV 610
    8 IF(MODE-1) 10, 10, 20                                             CONV 620
C                                                                       CONV 630
C        SINGLE PRECISION TO DOUBLE PRECISION                           CONV 640
C                                                                       CONV 650
   10 DO 15 L=1,NM                                                      CONV 660
   15 D(L)=S(L)                                                         CONV 670
      GO TO 30                                                          CONV 680
C                                                                       CONV 690
C        DOUBLE PRECISION TO SINGLE PRECISION                           CONV 700
C                                                                       CONV 710
   20 DO 25 L=1,NM                                                      CONV 720
   25 S(L)=D(L)                                                         CONV 730
C                                                                       CONV 740
   30 RETURN                                                            CONV 750
      END                                                               CONV 760
