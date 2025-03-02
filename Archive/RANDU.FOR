C                                                                       RAND  10
C     ..................................................................RAND  20
C                                                                       RAND  30
C        SUBROUTINE RANDU                                               RAND  40
C                                                                       RAND  50
C        PURPOSE                                                        RAND  60
C           COMPUTES UNIFORMLY DISTRIBUTED RANDOM REAL NUMBERS BETWEEN  RAND  70
C           0 AND 1.0 AND RANDOM INTEGERS BETWEEN ZERO AND              RAND  80
C           2**31. EACH ENTRY USES AS INPUT AN INTEGER RANDOM NUMBER    RAND  90
C           AND PRODUCES A NEW INTEGER AND REAL RANDOM NUMBER.          RAND 100
C                                                                       RAND 110
C        USAGE                                                          RAND 120
C           CALL RANDU(IX,IY,YFL)                                       RAND 130
C                                                                       RAND 140
C        DESCRIPTION OF PARAMETERS                                      RAND 150
C           IX - FOR THE FIRST ENTRY THIS MUST CONTAIN ANY ODD INTEGER  RAND 160
C                NUMBER WITH NINE OR LESS DIGITS. AFTER THE FIRST ENTRY,RAND 170
C                IX SHOULD BE THE PREVIOUS VALUE OF IY COMPUTED BY THIS RAND 180
C                SUBROUTINE.                                            RAND 190
C           IY - A RESULTANT INTEGER RANDOM NUMBER REQUIRED FOR THE NEXTRAND 200
C                ENTRY TO THIS SUBROUTINE. THE RANGE OF THIS NUMBER IS  RAND 210
C                BETWEEN ZERO AND 2**31                                 RAND 220
C           YFL- THE RESULTANT UNIFORMLY DISTRIBUTED, FLOATING POINT,   RAND 230
C                RANDOM NUMBER IN THE RANGE 0 TO 1.0                    RAND 240
C                                                                       RAND 250
C        REMARKS                                                        RAND 260
C           THIS SUBROUTINE IS SPECIFIC TO SYSTEM/360 AND WILL PRODUCE  RAND 270
C           2**29 TERMS BEFORE REPEATING.  THE REFERENCE BELOW DISCUSSESRAND 280
C           SEEDS (65539 HERE), RUN PROBLEMS, AND PROBLEMS CONCERNING   RAND 290
C           RANDOM DIGITS USING THIS GENERATION SCHEME.  MACLAREN AND   RAND 300
C           MARSAGLIA, JACM 12, P. 83-89, DISCUSS CONGRUENTIAL          RAND 310
C           GENERATION METHODS AND TESTS.  THE USE OF TWO GENERATORS OF RAND 320
C           THE RANDU TYPE, ONE FILLING A TABLE AND ONE PICKING FROM THERAND 330
C           TABLE, IS OF BENEFIT IN SOME CASES.  65549 HAS BEEN         RAND 340
C           SUGGESTED AS A SEED WHICH HAS BETTER STATISTICAL PROPERTIES RAND 350
C           FOR HIGH ORDER BITS OF THE GENERATED DEVIATE.               RAND 360
C           SEEDS SHOULD BE CHOSEN IN ACCORDANCE WITH THE DISCUSSION    RAND 370
C           GIVEN IN THE REFERENCE BELOW.  ALSO, IT SHOULD BE NOTED THATRAND 380
C           IF FLOATING POINT RANDOM NUMBERS ARE DESIRED,AS ARE         RAND 390
C           AVAILABLE FROM RANDU, THE RANDOM CHARACTERISTICS OF THE     RAND 400
C           FLOATING POINT DEVIATES ARE MODIFIED AND IN FACT THESE      RAND 410
C           DEVIATES HAVE HIGH PROBABILITY OF HAVING A TRAILING LOW     RAND 420
C           ORDER ZERO BIT IN THEIR FRACTIONAL PART.                    RAND 430
C                                                                       RAND 440
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  RAND 450
C           NONE                                                        RAND 460
C                                                                       RAND 470
C        METHOD                                                         RAND 480
C           POWER RESIDUE METHOD DISCUSSED IN IBM MANUAL C20-8011,      RAND 490
C           RANDOM NUMBER GENERATION AND TESTING                        RAND 500
C                                                                       RAND 510
C     ..................................................................RAND 520
C                                                                       RAND 530
      SUBROUTINE RANDU(IX,IY,YFL)                                       RAND 540
      IY=IX*65539                                                       RAND 550
      IF(IY)5,6,6                                                       RAND 560
    5 IY=IY+2147483647+1                                                RAND 570
    6 YFL=IY                                                            RAND 580
      YFL=YFL*.4656613E-9                                               RAND 590
      RETURN                                                            RAND 600
      END                                                               RAND 610
