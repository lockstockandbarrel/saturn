C                                                                       SUBS  10
C     ..................................................................SUBS  20
C                                                                       SUBS  30
C        SUBROUTINE SUBST                                               SUBS  40
C                                                                       SUBS  50
C        PURPOSE                                                        SUBS  60
C           DERIVE A SUBSET VECTOR INDICATING WHICH OBSERVATIONS IN A   SUBS  70
C           SET HAVE SATISFIED CERTAIN CONDITIONS ON THE VARIABLES.     SUBS  80
C                                                                       SUBS  90
C        USAGE                                                          SUBS 100
C           CALL SUBST (A,C,R,B,S,NO,NV,NC)                             SUBS 110
C           PARAMETER B MUST BE DEFINED BY AN EXTERNAL STATEMENT IN THE SUBS 120
C           CALLING PROGRAM                                             SUBS 130
C                                                                       SUBS 140
C        DESCRIPTION OF PARAMETERS                                      SUBS 150
C           A  - OBSERVATION MATRIX, NO BY NV                           SUBS 160
C           C  - INPUT MATRIX, 3 BY NC, OF CONDITIONS TO BE CONSIDERED. SUBS 170
C                THE FIRST ELEMENT OF EACH COLUMN OF C REPRESENTS THE   SUBS 180
C                NUMBER OF THE VARIABLE (COLUMN OF THE MATRIX A) TO BE  SUBS 190
C                TESTED, THE SECOND ELEMENT OF EACH COLUMN IS A         SUBS 200
C                RELATIONAL CODE AS FOLLOWS                             SUBS 210
C                     1. FOR LT (LESS THAN)                             SUBS 220
C                     2. FOR LE (LESS THAN OR EQUAL TO)                 SUBS 230
C                     3. FOR EQ (EQUAL TO)                              SUBS 240
C                     4. FOR NE (NOT EQUAL TO)                          SUBS 250
C                     5. FOR GE (GREATER THAN OR EQUAL TO)              SUBS 260
C                     6. FOR GT (GREATER THAN)                          SUBS 270
C                THE THIRD ELEMENT OF EACH COLUMN IS A QUANTITY TO BE   SUBS 280
C                USED FOR COMPARISON WITH THE OBSERVATION VALUES. FOR   SUBS 290
C                EXAMPLE, THE FOLLOWING COLUMN IN C                     SUBS 300
C                          2.                                           SUBS 310
C                          5.                                           SUBS 320
C                         92.5                                          SUBS 330
C                CAUSES THE SECOND VARIABLE TO BE TESTED FOR GREATER    SUBS 340
C                THAN OR EQUAL TO 92.5                                  SUBS 350
C           R  - WORKING VECTOR USED TO STORE INTERMEDIATE RESULTS OF   SUBS 360
C                ABOVE TESTS ON A SINGLE OBSERVATION. IF CONDITION IS   SUBS 370
C                SATISFIED, R(I) IS SET TO 1. IF IT IS NOT, R(I) IS SET SUBS 380
C                TO 0. VECTOR LENGTH IS NC.                             SUBS 390
C           B  - NAME OF SUBROUTINE TO BE SUPPLIED BY THE USER. IT      SUBS 400
C                CONSISTS OF A BOOLEAN EXPRESSION LINKING THE           SUBS 410
C                INTERMEDIATE VALUES STORED IN VECTOR R. THE BOOLEAN    SUBS 420
C                OPERATORS ARE '*' FOR'AND', '+' FOR 'OR'. EXAMPLE      SUBS 430
C                     SUBROUTINE BOOL(R,T)                              SUBS 440
C                     DIMENSION R(3)                                    SUBS 450
C                     T=R(1)*(R(2)+R(3))                                SUBS 460
C                     RETURN                                            SUBS 470
C                     END                                               SUBS 480
C                THE ABOVE EXPRESSION IS TESTED FOR                     SUBS 490
C                     R(1).AND.(R(2).OR.R(3))                           SUBS 500
C           S  - OUTPUT VECTOR INDICATING, FOR EACH OBSERVATION,        SUBS 510
C                WHETHER OR NOT PROPOSITION B IS SATISFIED. IF IT IS,   SUBS 520
C                S(I) IS NON-ZERO. IF IT IS NOT, S(I) IS ZERO. VECTOR   SUBS 530
C                LENGTH IS NO.                                          SUBS 540
C           NO - NUMBER OF OBSERVATIONS. NO MUST BE > OR = TO 1.        SUBS 550
C           NV - NUMBER OF VARIABLES. NV MUST BE > OR = TO 1.           SUBS 560
C           NC - NUMBER OF BASIC CONDITIONS TO BE SATISFIED. NC MUST BE SUBS 570
C                GREATER THAN OR EQUAL TO 1.                            SUBS 571
C                                                                       SUBS 580
C        REMARKS                                                        SUBS 590
C           NONE                                                        SUBS 600
C                                                                       SUBS 610
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  SUBS 620
C           B  THE NAME OF ACTUAL SUBROUTINE SUPPLIED BY THE USER MAY   SUBS 630
C              BE DIFFERENT (E.G., BOOL), BUT SUBROUTINE SUBST ALWAYS   SUBS 640
C              CALLS IT AS B.  IN ORDER FOR SUBROUTINE SUBST TO DO THIS,SUBS 650
C              THE NAME OF THE USER-SUPPLIED SUBROUTINE MUST BE         SUBS 660
C              DEFINED BY AN EXTERNAL STATEMENT IN THE CALLING PROGRAM. SUBS 670
C              THE NAME MUST ALSO BE LISTED IN THE ''CALL SUBST''       SUBS 680
C              STATEMENT.  (SEE USAGE ABOVE)                            SUBS 690
C                                                                       SUBS 700
C        METHOD                                                         SUBS 710
C           THE FOLLOWING IS DONE FOR EACH OBSERVATION.                 SUBS 720
C           CONDITION MATRIX IS ANALYZED TO DETERMINE WHICH VARIABLES   SUBS 730
C           ARE TO BE EXAMINED. INTERMEDIATE VECTOR R IS FORMED. THE    SUBS 740
C           BOOLEAN EXPRESSION (IN SUBROUTINE B) IS THEN EVALUATED TO   SUBS 750
C           DERIVE THE ELEMENT IN SUBSET VECTOR S CORRESPONDING TO THE  SUBS 760
C           OBSERVATION.                                                SUBS 770
C                                                                       SUBS 780
C     ..................................................................SUBS 790
C                                                                       SUBS 800
      SUBROUTINE SUBST(A,C,R,B,S,NO,NV,NC)                              SUBS 810
      DIMENSION A(1),C(1),R(1),S(1)                                     SUBS 820
C                                                                       SUBS 830
      DO 9 I=1,NO                                                       SUBS 840
      IQ=I-NO                                                           SUBS 850
      K=-2                                                              SUBS 860
      DO 8 J=1,NC                                                       SUBS 870
C                                                                       SUBS 880
C        CLEAR R VECTOR                                                 SUBS 890
C                                                                       SUBS 900
      R(J)=0.0                                                          SUBS 910
C                                                                       SUBS 920
C         LOCATE ELEMENT IN OBSERVATION MATRIX AND RELATIONAL CODE      SUBS 930
C                                                                       SUBS 940
      K=K+3                                                             SUBS 950
      IZ=C(K)                                                           SUBS 960
      IA=IQ+IZ*NO                                                       SUBS 970
      IGO=C(K+1)                                                        SUBS 980
C                                                                       SUBS 990
C         FORM R VECTOR                                                 SUBS1000
C                                                                       SUBS1010
      Q=A(IA)-C(K+2)                                                    SUBS1020
      GO TO(1,2,3,4,5,6),IGO                                            SUBS1030
    1 IF(Q) 7,8,8                                                       SUBS1040
    2 IF(Q) 7,7,8                                                       SUBS1050
    3 IF(Q) 8,7,8                                                       SUBS1060
    4 IF(Q) 7,8,7                                                       SUBS1070
    5 IF(Q) 8,7,7                                                       SUBS1080
    6 IF(Q) 8,8,7                                                       SUBS1090
    7 R(J)=1.0                                                          SUBS1100
    8 CONTINUE                                                          SUBS1110
C                                                                       SUBS1120
C        CALCULATE S VECTOR                                             SUBS1130
C                                                                       SUBS1140
    9 CALL B(R,S(I))                                                    SUBS1150
      RETURN                                                            SUBS1160
      END                                                               SUBS1170
