C                                                                       ARRA  10
C     ..................................................................ARRA  20
C                                                                       ARRA  30
C        SUBROUTINE ARRAY                                               ARRA  40
C                                                                       ARRA  50
C        PURPOSE                                                        ARRA  60
C           CONVERT DATA ARRAY FROM SINGLE TO DOUBLE DIMENSION OR VICE  ARRA  70
C           VERSA.  THIS SUBROUTINE IS USED TO LINK THE USER PROGRAM    ARRA  80
C           WHICH HAS DOUBLE DIMENSION ARRAYS AND THE SSP SUBROUTINES   ARRA  90
C           WHICH OPERATE ON ARRAYS OF DATA IN A VECTOR FASHION.        ARRA 100
C                                                                       ARRA 110
C        USAGE                                                          ARRA 120
C           CALL ARRAY (MODE,I,J,N,M,S,D)                               ARRA 130
C                                                                       ARRA 140
C        DESCRIPTION OF PARAMETERS                                      ARRA 150
C           MODE - CODE INDICATING TYPE OF CONVERSION                   ARRA 160
C                    1 - FROM SINGLE TO DOUBLE DIMENSION                ARRA 170
C                    2 - FROM DOUBLE TO SINGLE DIMENSION                ARRA 180
C           I    - NUMBER OF ROWS IN ACTUAL DATA MATRIX                 ARRA 190
C           J    - NUMBER OF COLUMNS IN ACTUAL DATA MATRIX              ARRA 200
C           N    - NUMBER OF ROWS SPECIFIED FOR THE MATRIX D IN         ARRA 210
C                  DIMENSION STATEMENT                                  ARRA 220
C           M    - NUMBER OF COLUMNS SPECIFIED FOR THE MATRIX D IN      ARRA 230
C                  DIMENSION STATEMENT                                  ARRA 240
C           S    - IF MODE=1, THIS VECTOR IS INPUT WHICH CONTAINS THE   ARRA 250
C                  ELEMENTS OF A DATA MATRIX OF SIZE I BY J. COLUMN I+1 ARRA 260
C                  OF DATA MATRIX FOLLOWS COLUMN I, ETC. IF MODE=2,     ARRA 270
C                  THIS VECTOR IS OUTPUT REPRESENTING A DATA MATRIX OF  ARRA 280
C                  SIZE I BY J CONTAINING ITS COLUMNS CONSECUTIVELY.    ARRA 290
C                  THE LENGTH OF S IS IJ, WHERE IJ=I*J.                 ARRA 300
C           D    - IF MODE=1, THIS MATRIX OF SIZE N BY M IS OUTPUT,     ARRA 310
C                  CONTAINING A DATA MATRIX OF SIZE I BY J IN THE FIRST ARRA 320
C                  I ROWS AND J COLUMNS. IF MODE=2, THIS N BY M MATRIX  ARRA 330
C                  IS INPUT CONTAINING A DATA MATRIX OF SIZE I BY J IN  ARRA 340
C                  THE FIRST I ROWS AND J COLUMNS.                      ARRA 350
C                                                                       ARRA 360
C        REMARKS                                                        ARRA 370
C           VECTOR S CAN BE IN THE SAME LOCATION AS MATRIX D.  VECTOR S ARRA 380
C           IS REFERRED AS A MATRIX IN OTHER SSP ROUTINES, SINCE IT     ARRA 390
C           CONTAINS A DATA MATRIX.                                     ARRA 400
C           THIS SUBROUTINE CONVERTS ONLY GENERAL DATA MATRICES (STORAGEARRA 410
C           MODE OF 0).                                                 ARRA 420
C                                                                       ARRA 430
C        SUBROUTINES AND FUNCTION SUBROUTINES REQUIRED                  ARRA 440
C           NONE                                                        ARRA 450
C                                                                       ARRA 460
C        METHOD                                                         ARRA 470
C           REFER TO THE DISCUSSION ON VARIABLE DATA SIZE IN THE SECTIONARRA 480
C           DESCRIBING OVERALL RULES FOR USAGE IN THIS MANUAL.          ARRA 490
C                                                                       ARRA 500
C     ..................................................................ARRA 510
C                                                                       ARRA 520
      SUBROUTINE ARRAY (MODE,I,J,N,M,S,D)                               ARRA 530
      DIMENSION S(1),D(1)                                               ARRA 540
C                                                                       ARRA 550
      NI=N-I                                                            ARRA 560
C                                                                       ARRA 570
C        TEST TYPE OF CONVERSION                                        ARRA 580
C                                                                       ARRA 590
      IF(MODE-1) 100, 100, 120                                          ARRA 600
C                                                                       ARRA 610
C        CONVERT FROM SINGLE TO DOUBLE DIMENSION                        ARRA 620
C                                                                       ARRA 630
  100 IJ=I*J+1                                                          ARRA 640
      NM=N*J+1                                                          ARRA 650
      DO 110 K=1,J                                                      ARRA 660
      NM=NM-NI                                                          ARRA 670
      DO 110 L=1,I                                                      ARRA 680
      IJ=IJ-1                                                           ARRA 690
      NM=NM-1                                                           ARRA 700
  110 D(NM)=S(IJ)                                                       ARRA 710
      GO TO 140                                                         ARRA 720
C                                                                       ARRA 730
C        CONVERT FROM DOUBLE TO SINGLE DIMENSION                        ARRA 740
C                                                                       ARRA 750
  120 IJ=0                                                              ARRA 760
      NM=0                                                              ARRA 770
      DO 130 K=1,J                                                      ARRA 780
      DO 125 L=1,I                                                      ARRA 790
      IJ=IJ+1                                                           ARRA 800
      NM=NM+1                                                           ARRA 810
  125 S(IJ)=D(NM)                                                       ARRA 820
  130 NM=NM+NI                                                          ARRA 830
C                                                                       ARRA 840
  140 RETURN                                                            ARRA 850
      END                                                               ARRA 860
