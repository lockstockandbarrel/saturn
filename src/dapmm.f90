!*==dapmm.f90 processed by SPAG 8.01RF 00:34  2 Mar 2025
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
!     ..................................................................
!
!        SUBROUTINE DAPMM
!
!        PURPOSE
!           APPROXIMATE A FUNCTION TABULATED IN N POINTS BY ANY LINEAR
!           COMBINATION OF M GIVEN CONTINUOUS FUNCTIONS IN THE SENSE
!           OF CHEBYSHEV.
!
!        USAGE
!           CALL DAPMM(FCT,N,M,TOP,IHE,PIV,T,ITER,IER)
!           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT IN THE
!           CALLING PROGRAM.
!
!        DESCRIPTION OF PARAMETERS
!           FCT    - NAME OF SUBROUTINE TO BE SUPPLIED BY THE USER.
!                    IT COMPUTES VALUES OF M GIVEN FUNCTIONS FOR
!                    ARGUMENT VALUE X.
!                    USAGE
!                       CALL FCT(Y,X,K)
!                    DESCRIPTION OF PARAMETERS
!                       Y   - RESULT VECTOR OF DIMENSION M CONTAINING
!                             THE VALUES OF GIVEN CONTINUOUS FUNCTIONS
!                             FOR GIVEN ARGUMENT X
!                       X   - ARGUMENT VALUE
!                       K   - AN INTEGER VALUE WHICH IS EQUAL TO M-1
!                    REMARKS
!                       IF APPROXIMATION BY NORMAL CHEBYSHEV, SHIFTED
!                       CHEBYSHEV, LEGENDRE, LAGUERRE, HERMITE POLYNO-
!                       MIALS IS DESIRED SUBROUTINES CNP, CSP, LEP,
!                       LAP, HEP, RESPECTIVELY FROM SSP COULD BE USED.
!           N      - NUMBER OF DATA POINTS DEFINING THE FUNCTION WHICH
!                    IS TO BE APPROXIMATED
!           M      - NUMBER OF GIVEN CONTINUOUS FUNCTIONS FROM WHICH
!                    THE APPROXIMATING FUNCTION IS CONSTRUCTED.
!           TOP    - VECTOR OF DIMENSION 3*N.
!                    ON ENTRY IT MUST CONTAIN FROM TOP(1) UP TO TOP(N)
!                    THE GIVEN N FUNCTION VALUES AND FROM TOP(N+1) UP
!                    TO TOP(2*N) THE CORRESPONDING NODES
!                    ON RETURN TOP CONTAINS FROM TOP(1) UP TO TOP(N)
!                    THE ERRORS AT THOSE N NODES.
!                    OTHER VALUES OF TOP ARE SCRATCH.
!           IHE    - INTEGER VECTOR OF DIMENSION 3*M+4*N+6
!           PIV    - VECTOR OF DIMENSION 3*M+6.
!                    ON RETURN PIV CONTAINS AT PIV(1) UP TO PIV(M) THE
!                    RESULTING COEFFICIENTS OF LINEAR APPROXIMATION.
!           T      - AUXILIARY VECTOR OF DIMENSION (M+2)*(M+2)
!           ITER   - RESULTANT INTEGER WHICH SPECIFIES THE NUMBER OF
!                    ITERATIONS NEEDED
!           IER    - RESULTANT ERROR PARAMETER CODED IN THE FOLLOWING
!                    FORM
!                     IER=0  - NO ERROR
!                     IER=1  - THE NUMBER OF ITERATIONS HAS REACHED
!                              THE INTERNAL MAXIMUM N+M
!                     IER=-1 - NO RESULT BECAUSE OF WRONG INPUT PARA-
!                              METER M OR N OR SINCE AT SOME ITERATION
!                              NO SUITABLE PIVOT COULD BE FOUND
!        NOTE: DSUM,TOP,PIV,T,SAVE,HELP,RECI,TOL ARE DOUBLE PRECISION
!        REMARKS
!           NO ACTION BESIDES ERROR MESSAGE IN CASE M LESS THAN 1 OR
!           N LESS THAN 2.
!
!        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
!           THE EXTERNAL SUBROUTINE FCT MUST BE FURNISHED BY THE USER.
!
!        METHOD
!           THE PROBLEM OF APPROXIMATION A TABULATED FUNCTION BY ANY
!           LINEAR COMBINATION OF GIVEN FUNCTIONS IN THE SENSE OF
!           CHEBYSHEV (I.E. TO MINIMIZE THE MAXIMUM ERROR) IS TRANS-
!           FORMED INTO A LINEAR PROGRAMMING PROBLEM. APMM USES A
!           REVISED SIMPLEX METHOD TO SOLVE A CORRESPONDING DUAL
!           PROBLEM. FOR REFERENCE, SEE
!           I.BARRODALE/A.YOUNG, ALGORITHMS FOR BEST L-SUB-ONE AND
!           L-SUB-INFINITY, LINEAR APPROXIMATIONS ON A DISCRETE SET,
!           NUMERISCHE MATHEMATIK, VOL.8, ISS.3 (1966), PP.295-306.
!
!     ..................................................................
!
SUBROUTINE dapmm(fct,N,M,Top,Ihe,Piv,T,Iter,Ier)
   IMPLICIT NONE
!*** Start of declarations inserted by SPAG
   INTEGER i , ido , Ier , Ihe , ilab , ind , ipiv , ise , Iter , j , k , l , ll , M , N , nan
!*** End of declarations inserted by SPAG
!
   DOUBLE PRECISION dsum , Top , Piv , T , save , help , repi , tol
   DIMENSION Top(1) , Ihe(1) , Piv(1) , T(1)
   INTEGER :: spag_nextblock_1
   spag_nextblock_1 = 1
   SPAG_DispatchLoop_1: DO
      SELECT CASE (spag_nextblock_1)
      CASE (1)
!
!        TEST ON WRONG INPUT PARAMETERS N AND M
         Ier = -1
         IF ( N<=1 ) RETURN
         IF ( M<=0 ) RETURN
!
!        INITIALIZE CHARACTERISTIC VECTORS FOR THE TABLEAU
         Ier = 0
!
!        PREPARE TOP-ROW TOP
         DO i = 1 , N
            k = i + N
            j = k + N
            Top(j) = Top(k)
            Top(k) = -Top(i)
         ENDDO
!
!        PREPARE INVERSE TRANSFORMATION MATRIX T
         l = M + 2
         ll = l*l
         DO i = 1 , ll
            T(i) = 0.
         ENDDO
         k = 1
         j = l + 1
         DO i = 1 , l
            T(k) = 1.
            k = k + j
         ENDDO
!
!        PREPARE INDEX-VECTOR IHE
         DO i = 1 , l
            k = i + l
            j = k + l
            Ihe(i) = 0
            Ihe(k) = i
            Ihe(j) = 1 - i
         ENDDO
         nan = N + N
         k = l + l + l
         j = k + nan
         DO i = 1 , nan
            k = k + 1
            Ihe(k) = i
            j = j + 1
            Ihe(j) = i
         ENDDO
!
!        SET COUNTER ITER FOR ITERATION-STEPS
         Iter = -1
         spag_nextblock_1 = 2
      CASE (2)
         Iter = Iter + 1
!
!        TEST FOR MAXIMUM ITERATION-STEPS
         IF ( N+M<=Iter ) THEN
            Ier = 1
            spag_nextblock_1 = 11
            CYCLE SPAG_DispatchLoop_1
         ELSE
!
!        DETERMINE THE COLUMN WITH THE MOST POSITIVE ELEMENT IN TOP
            ise = 0
            ipiv = 0
            k = l + l + l
            save = 0.0D0
!
!        START TOP-LOOP
            DO i = 1 , nan
               ido = k + i
               help = Top(i)
               IF ( help>save ) THEN
                  save = help
                  ipiv = i
               ENDIF
               IF ( Ihe(ido)==0 ) ise = i
            ENDDO
!        END OF TOP-LOOP
!
!        IS OPTIMAL TABLEAU REACHED
            IF ( ipiv<=0 ) THEN
               spag_nextblock_1 = 11
               CYCLE SPAG_DispatchLoop_1
            ENDIF
!
!        DETERMINE THE PIVOT-ELEMENT FOR THE COLUMN CHOSEN UPOVE
            ilab = 1
            ind = 0
            j = ise
            IF ( j<=0 ) j = ipiv
            spag_nextblock_1 = 6
            CYCLE SPAG_DispatchLoop_1
         ENDIF
      CASE (3)
!
!        TRANSFER K-TH COLUMN FROM T TO PIV
         k = (k-1)*l
         DO i = 1 , l
            j = l + i
            k = k + 1
            Piv(j) = T(k)
         ENDDO
         spag_nextblock_1 = 4
      CASE (4)
!
!        IS ANOTHER COLUMN NEEDED FOR SEARCH FOR PIVOT-ELEMENT
         IF ( ise<=0 ) THEN
!
!        SEARCH PIVOT-ELEMENT PIV(IND)
            save = 1.D38
            ido = 0
            k = l + 1
            ll = l + l
            ind = 0
!
!        START PIVOT-LOOP
            DO i = k , ll
               j = i + l
               help = Piv(i)
               IF ( help>0 ) THEN
                  help = -help
                  IF ( ise/=0 ) THEN
                     help = -Piv(j)/help
                  ELSEIF ( Ihe(j)==0 ) THEN
                     ido = i
                     CYCLE
                  ENDIF
                  IF ( help<save ) THEN
                     save = help
                     ind = i
                  ENDIF
               ENDIF
            ENDDO
!        END OF PIVOT-LOOP
!
!        TEST FOR SUITABLE PIVOT-ELEMENT
            IF ( ind<=0 ) THEN
               IF ( ido<=0 ) THEN
!
!        SET ERROR PARAMETER IER=-1 SINCE NO SUITABLE PIVOT IS FOUND
                  Ier = -1
                  spag_nextblock_1 = 11
                  CYCLE SPAG_DispatchLoop_1
               ELSE
                  ind = ido
               ENDIF
            ENDIF
!        PIVOT-ELEMENT IS STORED IN PIV(IND)
!
!        COMPUTE THE RECIPROCAL OF THE PIVOT-ELEMENT REPI
            repi = 1.0D0/Piv(ind)
            ind = ind - l
!
!        UPDATE THE TOP-ROW TOP OF THE TABLEAU
            ilab = 0
            save = -Top(ipiv)*repi
            Top(ipiv) = save
!
!        INITIALIZE J AS COUNTER FOR TOP-LOOP
            j = nan
         ELSE
            ise = -ise
!
!        TRANSFER COLUMNS IN PIV
            j = l + 1
            ido = l + l
            DO i = j , ido
               k = i + l
               Piv(k) = Piv(i)
            ENDDO
            j = ipiv
            spag_nextblock_1 = 6
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         spag_nextblock_1 = 5
      CASE (5)
         IF ( j==ipiv ) THEN
            spag_nextblock_1 = 9
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         spag_nextblock_1 = 6
      CASE (6)
         k = 0
!
!        SEARCH COLUMN IN TRANSFORMATION-MATRIX T
         DO i = 1 , l
            IF ( Ihe(i)==j ) THEN
               k = i
               IF ( ilab>0 ) THEN
                  spag_nextblock_1 = 3
                  CYCLE SPAG_DispatchLoop_1
               ENDIF
               spag_nextblock_1 = 7
               CYCLE SPAG_DispatchLoop_1
            ENDIF
         ENDDO
!
!        GENERATE COLUMN USING SUBROUTINE FCT AND TRANSFORMATION-MATRIX
         i = l + l + l + nan + j
         i = Ihe(i) - N
         IF ( i<=0 ) THEN
            i = i + N
            k = 1
         ENDIF
         i = i + nan
!
!        CALL SUBROUTINE FCT
         CALL fct(Piv,Top(i),M-1)
!
!        PREPARE THE CALLED VECTOR PIV
         dsum = 0.D0
         ido = M
         DO i = 1 , M
            help = Piv(ido)
            IF ( k<=0 ) help = -help
            dsum = dsum + dble(help)
            Piv(ido+1) = help
            ido = ido - 1
         ENDDO
         Piv(l) = -dsum
         Piv(1) = 1.0D0
!
!        TRANSFORM VECTOR PIV WITH ROWS OF MATRIX T
         ido = ind
         IF ( ilab>0 ) THEN
            k = 1
            ido = k
         ENDIF
         DO
            dsum = 0.D0
            help = 0.0D0
!
!        START MULTIPLICATION-LOOP
            DO i = 1 , l
               dsum = dsum + dble(Piv(i)*T(ido))
               TOL=ABS(SNGL(DSUM))
               IF ( tol>help ) help = tol
               ido = ido + l
            ENDDO
!        END OF MULTIPLICATION-LOOP
!
            tol = 1.D-5*help
            IF ( abs(sngl(dsum))<=tol ) dsum = 0.D0
            IF ( ilab<=0 ) THEN
               spag_nextblock_1 = 8
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            i = k + l
            Piv(i) = dsum
!
!        TEST FOR LAST COLUMN-TERM
            k = k + 1
            IF ( k>l ) THEN
               spag_nextblock_1 = 4
               CYCLE SPAG_DispatchLoop_1
            ENDIF
            ido = k
         ENDDO
         spag_nextblock_1 = 7
      CASE (7)
         i = (k-1)*l + ind
         dsum = T(i)
         spag_nextblock_1 = 8
      CASE (8)
!
!        COMPUTE NEW TOP-ELEMENT
         dsum = dsum*dble(save)
         tol = 1.D-5*abs(sngl(dsum))
         Top(j) = Top(j) + sngl(dsum)
         IF ( dabs(Top(j))<=tol ) Top(j) = 0.D0
         spag_nextblock_1 = 9
      CASE (9)
!
!        TEST FOR LAST TOP-TERM
         j = j - 1
         IF ( j>0 ) THEN
            spag_nextblock_1 = 5
            CYCLE SPAG_DispatchLoop_1
         ENDIF
!        END OF TOP-LOOP
!
!        TRANSFORM PIVOT-COLUMN
         i = ind + l
         Piv(i) = -1.
         DO i = 1 , l
            j = i + l
            Piv(i) = -Piv(j)*repi
         ENDDO
!
!        UPDATE TRANSFORMATION-MATRIX T
         j = 0
         DO i = 1 , l
            ido = j + ind
            save = T(ido)
            T(ido) = 0.D0
            DO k = 1 , l
               ise = k + j
               T(ise) = T(ise) + save*Piv(k)
            ENDDO
            j = j + l
         ENDDO
!
!        UPDATE INDEX-VECTOR IHE
!        INITIALIZE CHARACTERISTICS
         j = 0
         k = 0
         ise = 0
         ido = 0
!
!        START QUESTION-LOOP
         DO i = 1 , l
            ll = i + l
            ilab = Ihe(ll)
            IF ( Ihe(i)==ipiv ) THEN
               ise = i
               j = ilab
            ENDIF
            IF ( ilab==ind ) THEN
               ido = i
               k = Ihe(i)
            ENDIF
         ENDDO
!        END OF QUESTION-LOOP
!
!        START MODIFICATION
         IF ( k<=0 ) THEN
            Ihe(ido) = ipiv
            IF ( ise<=0 ) THEN
               spag_nextblock_1 = 10
               CYCLE SPAG_DispatchLoop_1
            ENDIF
         ELSEIF ( ind/=j ) THEN
            ll = l + l + l + nan
            k = k + ll
            i = ipiv + ll
            ilab = Ihe(k)
            Ihe(k) = Ihe(i)
            Ihe(i) = ilab
            IF ( ise<=0 ) THEN
               spag_nextblock_1 = 10
               CYCLE SPAG_DispatchLoop_1
            ENDIF
         ELSE
            Ihe(ise) = 0
            spag_nextblock_1 = 10
            CYCLE SPAG_DispatchLoop_1
         ENDIF
         ido = ido + l
         i = ise + l
         Ihe(ido) = j
         Ihe(i) = ind
         Ihe(ise) = 0
         spag_nextblock_1 = 10
      CASE (10)
         ll = l + l
         j = ll + ind
         i = ll + l + ipiv
         ilab = Ihe(i)
         Ihe(i) = Ihe(j)
!        END OF MODIFICATION
!
         Ihe(j) = ilab
         spag_nextblock_1 = 2
         CYCLE SPAG_DispatchLoop_1
      CASE (11)
!
!        EVALUATE FINAL TABLEAU
!        COMPUTE SAVE AS MAXIMUM ERROR OF APPROXIMATION AND
!        HELP AS ADDITIVE CONSTANCE FOR RESULTING COEFFICIENTS
         save = 0.D0
         help = 0.D0
         k = l + l + l
         DO i = 1 , nan
            ido = k + i
            j = Ihe(ido)
            IF ( j<0 ) THEN
            ELSEIF ( j==0 ) THEN
               save = -Top(i)
            ELSE
               CYCLE
            ENDIF
            IF ( M+j+1==0 ) help = Top(i)
         ENDDO
!
!        PREPARE T,TOP,PIV
         T(1) = save
         ido = nan + 1
         j = nan + N
         DO i = ido , j
            Top(i) = save
         ENDDO
         DO i = 1 , M
            Piv(i) = help
         ENDDO
!
!        COMPUTE COEFFICIENTS OF RESULTING POLYNOMIAL IN PIV(1) UP TO PI
!        AND CALCULATE ERRORS AT GIVEN NODES IN TOP(1) UP TO TOP(N)
         DO i = 1 , nan
            ido = k + i
            j = Ihe(ido)
            IF ( j<0 ) THEN
               j = -j
               Piv(j) = help - Top(i)
            ELSEIF ( j/=0 ) THEN
               IF ( j<=N ) THEN
                  j = j + nan
                  Top(j) = save + Top(i)
               ENDIF
            ENDIF
         ENDDO
         DO i = 1 , N
            ido = nan + i
            Top(i) = Top(ido)
         ENDDO
         EXIT SPAG_DispatchLoop_1
      END SELECT
   ENDDO SPAG_DispatchLoop_1
END SUBROUTINE dapmm
