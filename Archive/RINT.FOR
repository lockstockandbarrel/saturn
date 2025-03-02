C                                                                       RINT  10
C     ..................................................................RINT  20
C                                                                       RINT  30
C        SUBROUTINE RINT                                                RINT  40
C                                                                       RINT  50
C        PURPOSE                                                        RINT  60
C           INTERCHANGE TWO ROWS OF A MATRIX                            RINT  70
C                                                                       RINT  80
C        USAGE                                                          RINT  90
C           CALL RINT(A,N,M,LA,LB)                                      RINT 100
C                                                                       RINT 110
C        DESCRIPTION OF PARAMETERS                                      RINT 120
C           A  - NAME OF MATRIX                                         RINT 130
C           N  - NUMBER OF ROWS IN A                                    RINT 140
C           M  - NUMBER OF COLUMNS IN A                                 RINT 150
C           LA - ROW TO BE INTERCHANGED WITH ROW LB                     RINT 160
C           LB - ROW TO BE INTERCHANGED WITH ROW LA                     RINT 170
C                                                                       RINT 180
C        REMARKS                                                        RINT 190
C           MATRIX A MUST BE A GENERAL MATRIX                           RINT 200
C                                                                       RINT 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  RINT 220
C           NONE                                                        RINT 230
C                                                                       RINT 240
C        METHOD                                                         RINT 250
C           EACH ELEMENT OF ROW LA IS INTERCHANGED WITH CORRESPONDING   RINT 260
C           ELEMENT OF ROW LB                                           RINT 270
C                                                                       RINT 280
C     ..................................................................RINT 290
C                                                                       RINT 300
      SUBROUTINE RINT(A,N,M,LA,LB)                                      RINT 310
      DIMENSION A(1)                                                    RINT 320
C                                                                       RINT 330
      LAJ=LA-N                                                          RINT 340
      LBJ=LB-N                                                          RINT 350
      DO 3 J=1,M                                                        RINT 360
C                                                                       RINT 370
C        LOCATE ELEMENTS IN BOTH ROWS                                   RINT 380
C                                                                       RINT 390
      LAJ=LAJ+N                                                         RINT 400
      LBJ=LBJ+N                                                         RINT 410
C                                                                       RINT 420
C        INTERCHANGE ELEMENTS                                           RINT 430
C                                                                       RINT 440
      SAVE=A(LAJ)                                                       RINT 450
      A(LAJ)=A(LBJ)                                                     RINT 460
    3 A(LBJ)=SAVE                                                       RINT 470
      RETURN                                                            RINT 480
      END                                                               RINT 490
