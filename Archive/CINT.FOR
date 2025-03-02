C                                                                       CINT  10
C     ..................................................................CINT  20
C                                                                       CINT  30
C        SUBROUTINE CINT                                                CINT  40
C                                                                       CINT  50
C        PURPOSE                                                        CINT  60
C           INTERCHANGE TWO COLUMNS OF A MATRIX                         CINT  70
C                                                                       CINT  80
C        USAGE                                                          CINT  90
C           CALL CINT(A,N,LA,LB)                                        CINT 100
C                                                                       CINT 110
C        DESCRIPTION OF PARAMETERS                                      CINT 120
C           A  - NAME OF MATRIX                                         CINT 130
C           N  - NUMBER OF ROWS IN A                                    CINT 140
C           LA - COLUMN TO BE INTERCHANGED WITH COLUMN LB               CINT 150
C           LB - COLUMN TO BE INTERCHANGED WITH COLUMN LA               CINT 160
C                                                                       CINT 170
C        REMARKS                                                        CINT 180
C           MATRIX A MUST BE A GENERAL MATRIX                           CINT 190
C                                                                       CINT 200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  CINT 210
C           NONE                                                        CINT 220
C                                                                       CINT 230
C        METHOD                                                         CINT 240
C           EACH ELEMENT OF COLUMN LA IS INTERCHANGED WITH CORRESPONDINGCINT 250
C           ELEMENT OF COLUMN LB                                        CINT 260
C                                                                       CINT 270
C     ..................................................................CINT 280
C                                                                       CINT 290
      SUBROUTINE CINT(A,N,LA,LB)                                        CINT 300
      DIMENSION A(1)                                                    CINT 310
C                                                                       CINT 320
C        LOCATE STARTING POINT OF BOTH COLUMNS                          CINT 330
C                                                                       CINT 340
      ILA=N*(LA-1)                                                      CINT 350
      ILB=N*(LB-1)                                                      CINT 360
C                                                                       CINT 370
      DO 3 I=1,N                                                        CINT 380
      ILA=ILA+1                                                         CINT 390
      ILB=ILB+1                                                         CINT 400
C                                                                       CINT 410
C        INTERCHANGE ELEMENTS                                           CINT 420
C                                                                       CINT 430
      SAVE=A(ILA)                                                       CINT 440
      A(ILA)=A(ILB)                                                     CINT 450
    3 A(ILB)=SAVE                                                       CINT 460
      RETURN                                                            CINT 470
      END                                                               CINT 480
