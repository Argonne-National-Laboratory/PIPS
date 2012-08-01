C COPYRIGHT (c) 2002 ENSEEIHT-IRIT, Toulouse, France and
C  Council for the Central Laboratory of the Research Councils.
C  Version 1.0.0 July 2004
C  Version 1.0.1 March 2008  Comments reflowed with length < 73
C  AUTHOR Daniel Ruiz (Daniel.Ruiz@enseeiht.fr)
C *** Copyright (c) 2004  Council for the Central Laboratory of the
C     Research Councils and Ecole Nationale Superieure
C     d'Electrotechnique, d'Electronique, d'Informatique,
C     d'Hydraulique et des Telecommunications de Toulouse.          ***
C *** Although every effort has been made to ensure robustness and  ***
C *** reliability of the subroutines in this MC77 suite, we         ***
C *** disclaim any liability arising through the use or misuse of   ***
C *** any of the subroutines.                                       ***

      SUBROUTINE MC77I(ICNTL, CNTL)
C
C  Purpose
C  =======
C
C  The components of the array ICNTL control the action of MC77A/AD.
C  Default values for these are set in this subroutine.
C
C  Parameters
C  ==========
C
      INTEGER LICNTL, LCNTL
      PARAMETER ( LICNTL=10, LCNTL=10 )
      INTEGER ICNTL(LICNTL)
      REAL CNTL(LCNTL)
C
C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I
C
C    ICNTL(1) has default value 6.
C     It is the output stream for error messages. If it
C     is negative, these messages will be suppressed.
C
C    ICNTL(2) has default value 6.
C     It is the output stream for warning messages.
C     If it is negative, these messages are suppressed.
C
C    ICNTL(3) has default value -1.
C     It is the output stream for monitoring printing.
C     If it is negative, these messages are suppressed.
C
C    ICNTL(4) has default value 0.
C     If left at the default value, the incoming data is checked for
C     out-of-range indices and duplicates, in which case the driver
C     routine will exit with an error.  Setting ICNTL(4) to any
C     other value will avoid the checks but is likely to cause problems
C     later if out-of-range indices or duplicates are present.
C     The user should only set ICNTL(4) nonzero if the data is
C     known to be in range without duplicates.
C
C    ICNTL(5) has default value 0.
C     If left at the default value, it indicates that A contains some
C     negative entries, and it is necessary that their absolute values
C     be computed internally.  Otherwise, the values in the input
C     matrix A will be considered as non-negative.
C
C    ICNTL(6) has default value 0.
C     If nonzero, the input matrix A is symmetric and the user must
C     only supply the lower triangular part of A in the appropriate
C     format.  Entries in the upper triangular part of a symmetric
C     matrix will be considered as out-of-range, and are also checked
C     when ICNTL(4) is 0.
C
C    ICNTL(7) has a default value of 10.
C     It specifies the maximum number of scaling iterations that
C     may be performed.
C     Note that iteration "0", corresponding to the initial
C     normalization of the data, is always performed.
C     Restriction:  ICNTL(7) > 0
C                   (otherwise, the driver stops with an error)
C     ( ... In future release : Restriction:  ICNTL(7) >= 0 ... )
C
C    ICNTL(8) to ICNTL(15) are not currently used by MC77A/AD but are
C     set to zero in this routine.
C
C
C    CNTL(1) has a default value of 0.
C     It specifies the tolerance value when to stopping the iterations,
C     that is it is the desired value such that all row and column norms
C     in the scaled matrix lie between (1 +/- CNTL(1)).
C     If CNTL(1) is less than or equal to 0, tolerance is not checked,
C     and the algorithm will stop when the maximum number of iterations
C     given in ICNTL(7) is reached.
C
C    CNTL(2) has a default value of 1.
C     It is used in conjunction with parameter JOB set to -1,
C     to specify a REAL value for the power of norm under consideration.
C     Restriction:  CNTL(2) >= 1
C                   (otherwise, the driver stops with an error)
C
C    CNTL(3) to CNTL(10) are not currently used by MC77A/AD but are
C     set to zero in this routine.

C Initialization of the ICNTL array.
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = -1
      ICNTL(4) = 0
      ICNTL(5) = 0
      ICNTL(6) = 0
      ICNTL(7) = 10
C     Currently unused control variables:
      DO 10 I = 8,LICNTL
        ICNTL(I) = 0
   10 CONTINUE

C Initialization of the CNTL array.
      CNTL(1) = ZERO
      CNTL(2) = ONE
C     Currently unused control variables:
      DO 20 I = 3,LCNTL
        CNTL(I) = ZERO
   20 CONTINUE

      RETURN
      END

C**********************************************************************
C***           DRIVER FOR THE HARWELL-BOEING SPARSE FORMAT          ***
C**********************************************************************
      SUBROUTINE MC77A(JOB,M,N,NNZ,JCST,IRN,A,IW,LIW,DW,LDW,
     &                  ICNTL,CNTL,INFO,RINFO)
C
C  Purpose
C  =======
C
C This subroutine computes scaling factors D(i), i=1..M, and E(j),
C j=1..N, for an MxN sparse matrix A = {a_ij} so that the scaled matrix
C D^{-1}*A*E^{-1} has both rows and columns norms all equal or close to
C 1. The rectangular case (M/=N) is, for the moment, only allowed for
C equilibration in the infinity norm, a particular case for which the
C convergence is ensured in all cases and trivial to monitor. See [1].
C
C  Parameters
C  ==========
C
      INTEGER LICNTL, LCNTL, LINFO, LRINFO
      PARAMETER ( LICNTL=10, LCNTL=10, LINFO=10, LRINFO=10 )
      INTEGER ICNTL(LICNTL),INFO(LINFO)
      REAL CNTL(LCNTL),RINFO(LRINFO)
C
      INTEGER JOB,M,N,NNZ,LIW,LDW
      INTEGER JCST(N+1),IRN(NNZ),IW(LIW)
      REAL A(NNZ),DW(LDW)
C
C JOB is an INTEGER variable which must be set by the user to
C control the action. It is not altered by the subroutine.
C Possible values for JOB are:
C   0 Equilibrate the infinity norm of rows and columns in matrix A.
C   1 Equilibrate the one norm of rows and columns in matrix A.
C   p Equilibrate the p-th norm (p>=2) of rows and columns in matrix A.
C  -1 Equilibrate the p-th norm of rows and columns in matrix A, with
C     a strictly positive REAL parameter p given in CNTL(2).
C   Restriction: JOB >= -1.
C
C M is an INTEGER variable which must be set by the user to the
C   number of rows in matrix A. It is not altered by the subroutine.
C   Restriction: M >= 1.
C
C N is an INTEGER variable which must be set by the user to the
C   number of columns in matrix A. It is not altered by the subroutine.
C   Restriction: N >= 1.
C             If ICNTL(6) /= 0 (symmetric matrix), N must be equal to M.
C
C NNZ is an INTEGER variable which must be set by the user to the number
C   of entries in the matrix. It is not altered by the subroutine.
C   Restriction: NNZ >= 0.
C
C JCST is an INTEGER array of length N+1.
C   JCST(J), J=1..N, must be set by the user to the position in array
C   IRN of the first row index of an entry in column J.
C   JCST(N+1) must be set to NNZ+1.
C   The array JCST is not altered by the subroutine.
C
C IRN is an INTEGER array of length NNZ.
C   IRN(K), K=1..NNZ, must be set by the user to hold the row indices of
C   the entries in the matrix.
C   The array IRN is not altered by the subroutine.
C   Restrictions:
C     The entries in A belonging to column J must be stored contiguously
C     in the positions JCST(J)..JCST(J+1)-1. The ordering of the row
C     indices within each column is unimportant.
C     Out-of-range indices and duplicates are not allowed.
C
C A is a REAL (REAL in the D-version) array of length NNZ.
C   The user must set A(K), K=1..NNZ, to the numerical value of the
C   entry that corresponds to IRN(K).
C   It is not altered by the subroutine.
C
C LIW is an INTEGER variable that must be set by the user to
C   the length of array IW. It is not altered by the subroutine.
C   Restriction:
C     If ICNTL(6) == 0:  LIW >= M+N
C     If ICNTL(6) /= 0:  LIW >= M
C
C IW is an INTEGER array of length LIW that is used for workspace.
C
C LDW is an INTEGER variable that must be set by the user to the
C   length of array DW. It is not altered by the subroutine.
C   Restriction:
C     If ICNTL(5) = 0 and ICNTL(6) == 0:  LDW >= NNZ + 2*(M+N)
C     If ICNTL(5) = 1 and ICNTL(6) == 0:  LDW >= 2*(M+N)
C     If ICNTL(5) = 0 and ICNTL(6) /= 0:  LDW >= NNZ + 2*M
C     If ICNTL(5) = 1 and ICNTL(6) /= 0:  LDW >= 2*M
C
C DW is a REAL (REAL in the D-version) array of length LDW
C   that need not be set by the user.
C   On return, DW(i) contains d_i, i=1..M, the diagonal entries in
C   the row-scaling matrix D, and when the input matrix A is
C   unsymmetric, DW(M+j) contains e_j, j=1..N, the diagonal entries in
C   the column-scaling matrix E.
C
C ICNTL is an INTEGER array of length 10. Its components control the
C   output of MC77A/AD and must be set by the user before calling
C   MC77A/AD. They are not altered by the subroutine.
C   See MC77I/ID for details.
C
C CNTL is a REAL (REAL in the D-version) array of length 10.
C   Its components control the output of MC77A/AD and must be set by the
C   user before calling MC77A/AD. They are not altered by the
C   subroutine.
C   See MC77I/ID for details.
C
C INFO is an INTEGER array of length 10 which need not be set by the
C   user. INFO(1) is set non-negative to indicate success. A negative
C   value is returned if an error occurred, a positive value if a
C   warning occurred. INFO(2) holds further information on the error.
C   On exit from the subroutine, INFO(1) will take one of the
C   following values:
C    0 : successful entry (for structurally nonsingular matrix).
C   +1 : The maximum number of iterations in ICNTL(7) has been reached.
C   -1 : M < 1.  Value of M held in INFO(2).
C   -2 : N < 1.  Value of N held in INFO(2).
C   -3 : ICNTL(6) /= 0 (input matrix is symmetric) and N /= M.
C        Value of (N-M) held in INFO(2).
C   -4 : NNZ < 1. Value of NNZ held in INFO(2).
C   -5 : the defined length LIW violates the restriction on LIW.
C        Value of LIW required given by INFO(2).
C   -6 : the defined length LDW violates the restriction on LDW.
C        Value of LDW required given by INFO(2).
C   -7 : entries are found whose row indices are out of range. INFO(2)
C        contains the index of a column in which such an entry is found.
C   -8 : repeated entries are found. INFO(2) contains the index
C        of a column in which such an entry is found.
C   -9 : An entry corresponding to an element in the upper triangular
C        part of the matrix is found in the input data (ICNTL(6)} \= 0).
C  -10 : ICNTL(7) is out of range and INFO(2) contains its value.
C  -11 : CNTL(2) is out of range.
C  -12 : JOB < -1.  Value of JOB held in INFO(2).
C  -13 : JOB /= 0 and N /= M. This is a restriction of the algorithm in
C        its present stage, e.g. for rectangular matrices, only the
C        scaling in infinity norm is allowed. Value of JOB held in
C        INFO(2).
C   INFO(3) returns the number of iterations performed.
C   INFO(4) to INFO(10) are not currently used and are set to zero by
C        the routine.
C
C RINFO is a REAL (REAL in the D-version) array of length 10
C which need not be set by the user.
C   When CNTL(1) is strictly positive, RINFO(1) will contain on exit the
C   maximum distance of all row norms to 1, and RINFO(2) will contain on
C   exit the maximum distance of all column norms to 1.
C   If ICNTL(6) /= 0 (indicating that the input matrix is symmetric),
C   RINFO(2) is not used since the row-scaling equals the column scaling
C   in this case.
C   RINFO(3) to RINFO(10) are not currently used and are set to zero.
C
C References:
C  [1]  D. R. Ruiz, (2001),
C       "A scaling algorithm to equilibrate both rows and columns norms
C       in matrices",
C       Technical Report RAL-TR-2001-034, RAL, Oxfordshire, England,
C             and Report RT/APO/01/4, ENSEEIHT-IRIT, Toulouse, France.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      REAL THRESH, DP
      INTEGER I, J, K, MAXIT, CHECK, SETUP
C External routines and functions
      EXTERNAL MC77N,MC77O,MC77P,MC77Q
C Intrinsic functions
      INTRINSIC ABS,MAX

C Reset informational output values
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = 0
      RINFO(1) = ZERO
      RINFO(2) = ZERO
C Check value of JOB
      IF (JOB.LT.-1) THEN
        INFO(1) = -12
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
        GO TO 99
      ENDIF
C Check bad values in ICNTL
Cnew  IF (ICNTL(7).LT.0) THEN
      IF (ICNTL(7).LE.0) THEN
        INFO(1) = -10
        INFO(2) = 7
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9002) INFO(1),7,ICNTL(7)
        GO TO 99
      ENDIF
C Check bad values in CNTL
      IF ((JOB.EQ.-1) .AND. (CNTL(2).LT.ONE)) THEN
        INFO(1) = -11
        INFO(2) = 2
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9008) INFO(1),2,CNTL(2)
        GO TO 99
      ENDIF
C Check value of M
      IF (M.LT.1) THEN
        INFO(1) = -1
        INFO(2) = M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'M',M
        GO TO 99
      ENDIF
C Check value of N
      IF (N.LT.1) THEN
        INFO(1) = -2
        INFO(2) = N
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
        GO TO 99
      ENDIF
C Symmetric case: M must be equal to N
      IF ( (ICNTL(6).NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -3
        INFO(2) = N-M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9003) INFO(1),(N-M)
        GO TO 99
      ENDIF
C For rectangular matrices, only scaling in infinity norm is allowed
      IF ( (JOB.NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -13
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9009) INFO(1),JOB,M,N
        GO TO 99
      ENDIF
C Check value of NNZ
      IF (NNZ.LT.1) THEN
        INFO(1) = -4
        INFO(2) = NNZ
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'NNZ',NNZ
        GO TO 99
      ENDIF
C Check LIW
      K = M
      IF (ICNTL(6).EQ.0)  K = M+N
      IF (LIW.LT.K) THEN
        INFO(1) = -5
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
        GO TO 99
      ENDIF
C Check LDW
      K = 2*M
      IF (ICNTL(6).EQ.0)  K = 2*(M+N)
      IF ( (ICNTL(5).EQ.0) .OR. (JOB.GE.2) .OR. (JOB.EQ.-1) )
     &   K = K + NNZ
      IF (LDW.LT.K) THEN
        INFO(1) = -6
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
        GO TO 99
      ENDIF
C Check row indices. Use IW(1:M) as workspace
C     Harwell-Boeing sparse format :
      IF (ICNTL(4).EQ.0) THEN
        DO 3 I = 1,M
          IW(I) = 0
    3   CONTINUE
        DO 5 J = 1,N
          DO 4 K = JCST(J),JCST(J+1)-1
            I = IRN(K)
C Check for row indices that are out of range
            IF (I.LT.1 .OR. I.GT.M) THEN
              INFO(1) = -7
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),J,I
              GO TO 99
            ENDIF
            IF (ICNTL(6).NE.0 .AND. I.LT.J) THEN
              INFO(1) = -9
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),J,I
              GO TO 99
            ENDIF
C Check for repeated row indices within a column
            IF (IW(I).EQ.J) THEN
              INFO(1) = -8
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9007) INFO(1),J,I
              GO TO 99
            ELSE
              IW(I) = J
            ENDIF
    4     CONTINUE
    5   CONTINUE
      ENDIF

C Print diagnostics on input
C ICNTL(9) could be set by the user to monitor the level of diagnostics
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9020) JOB,M,N,NNZ,ICNTL(7),CNTL(1)
C       Harwell-Boeing sparse format :
        IF (ICNTL(9).LT.1) THEN
          WRITE(ICNTL(3),9021) (JCST(J),J=1,N+1)
          WRITE(ICNTL(3),9023) (IRN(J),J=1,NNZ)
          WRITE(ICNTL(3),9024) (A(J),J=1,NNZ)
        ENDIF
      ENDIF

C Set components of INFO to zero
      DO 10 I=1,10
        INFO(I)  = 0
        RINFO(I) = ZERO
   10 CONTINUE

C Set the convergence threshold and the maximum number of iterations
      THRESH = MAX( ZERO, CNTL(1) )
      MAXIT  = ICNTL(7)
C Check for the particular inconsistency between MAXIT and THRESH.
Cnew  IF ( (MAXIT.LE.0) .AND. (CNTL(1).LE.ZERO) )  THEN
Cnew     INFO(1) = -14
Cnew     IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9010) INFO(1),MAXIT,CNTL(1)
Cnew     GO TO 99
Cnew  ENDIF

C   ICNTL(8) could be set by the user to indicate the frequency at which
C   convergence should be checked when CNTL(1) is strictly positive.
C   The default value used here 1 (e.g. check convergence every
C   iteration).
C     CHECK  = ICNTL(8)
      CHECK  = 1
      IF (CNTL(1).LE.ZERO)  CHECK = 0

C Prepare temporary data in working arrays as apropriate
      K = 2*M
      IF (ICNTL(6).EQ.0)  K = 2*(M+N)
      SETUP = 0
      IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
        SETUP = 1
C       Copy ABS(A(J))**JOB into DW(K+J), J=1..NNZ
        IF (JOB.EQ.-1) THEN
          DP = CNTL(2)
        ELSE
          DP = REAL(JOB)
        ENDIF
        DO 20 J = 1,NNZ
          DW(K+J) = ABS(A(J))**DP
   20   CONTINUE
C       Reset the Threshold to take into account the Hadamard p-th power
        THRESH = ONE - (ABS(ONE-THRESH))**DP
      ELSE
        IF (ICNTL(5).EQ.0) THEN
          SETUP = 1
C         Copy ABS(A(J)) into DW(K+J), J=1..NNZ
          DO 30 J = 1,NNZ
            DW(K+J) = ABS(A(J))
   30     CONTINUE
        ENDIF
      ENDIF

C Begin the computations (input matrix in Harwell-Boeing SPARSE format)
C     We consider first the un-symmetric case :
      IF (ICNTL(6).EQ.0)  THEN
C
C       Equilibrate matrix A in the infinity norm
        IF (JOB.EQ.0) THEN
C         IW(1:M+N), DW(M+N:2*(M+N)) are workspaces
          IF (SETUP.EQ.1) THEN
            CALL MC77N(M,N,NNZ,JCST,IRN,DW(2*(M+N)+1),DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77N(M,N,NNZ,JCST,IRN,A,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
C
C       Equilibrate matrix A in the p-th norm, 1 <= p < +infinity
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77O(M,N,NNZ,JCST,IRN,DW(2*(M+N)+1),DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77O(M,N,NNZ,JCST,IRN,A,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
C         Reset the scaling factors to their Hadamard p-th root, if
C         required and re-compute the actual value of the final error
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 40 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+N+I)**DP) )
   40       CONTINUE
            RINFO(2) = ZERO
            DO 50 J = 1,N
              DW(M+J) = DW(M+J)**DP
              IF (IW(M+J).NE.0)
     &           RINFO(2) = MAX( RINFO(2), ABS(ONE-DW(2*M+N+J)**DP) )
   50       CONTINUE
          ENDIF
        ENDIF
C
C     We treat then the symmetric (packed) case :
      ELSE
C
C       Equilibrate matrix A in the infinity norm
        IF (JOB.EQ.0) THEN
C         IW(1:M), DW(M:2*M) are workspaces
          IF (SETUP.EQ.1) THEN
            CALL MC77P(M,NNZ,JCST,IRN,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77P(M,NNZ,JCST,IRN,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
C
C       Equilibrate matrix A in the p-th norm, 1 <= p < +infinity
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77Q(M,NNZ,JCST,IRN,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77Q(M,NNZ,JCST,IRN,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
C         Reset the scaling factors to their Hadamard p-th root, if
C         required and re-compute the actual value of the final error
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 60 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+I)**DP) )
   60       CONTINUE
          ENDIF
        ENDIF
        RINFO(2) = RINFO(1)
C
      ENDIF
C     End of the un-symmetric/symmetric IF statement
C     The Harwell-Boeing sparse case has been treated


C Print diagnostics on output
C ICNTL(9) could be set by the user to monitor the level of diagnostics
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9030) (INFO(J),J=1,3),(RINFO(J),J=1,2)
        IF (ICNTL(9).LT.2) THEN
          WRITE(ICNTL(3),9031) (DW(I),I=1,M)
          IF (ICNTL(6).EQ.0)
     &      WRITE(ICNTL(3),9032) (DW(M+J),J=1,N)
        ENDIF
      ENDIF

C Return from subroutine.
   99 CONTINUE
      RETURN

 9001 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3,
     &        ' because ',(A),' = ',I10)
 9002 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Bad input control flag.',
     &        ' Value of ICNTL(',I1,') = ',I8)
 9003 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Input matrix is symmetric and N /= M.',
     &        ' Value of (N-M) = ',I8)
 9004 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Column ',I8,
     &        ' contains an entry with invalid row index ',I8)
 9007 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Column ',I8,
     &        ' contains two or more entries with row index ',I8)
 9008 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Bad input REAL control parameter.',
     &        ' Value of CNTL(',I1,') = ',1PD14.4)
 9009 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
     &        '        Only scaling in infinity norm is allowed',
     &        ' for rectangular matrices'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8)
C9010 FORMAT (' ****** Error in MC77A/AD. INFO(1) = ',I3/
Cnew &        '        Algorithm will loop indefinitely !!!',
Cnew &        '     Check for inconsistency'/
Cnew &        ' Maximum number of iterations = ',I8/
Cnew &        ' Threshold for convergence    = ',1PD14.4)
 9020 FORMAT (' ****** Input parameters for MC77A/AD:'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8/' NNZ = ',I8/
     &        ' Max n.b. Iters.  = ',I8/
     &        ' Cvgce. Threshold = ',1PD14.4)
 9021 FORMAT (' JCST(1:N+1) = ',8I8/(15X,8I8))
 9023 FORMAT (' IRN(1:NNZ)  = ',8I8/(15X,8I8))
 9024 FORMAT (' A(1:NNZ)    = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for MC77A/AD:'/
     &        ' INFO(1:3)   = ',3(2X,I8)/
     &        ' RINFO(1:2)  = ',2(1PD14.4))
 9031 FORMAT (' DW(1:M)     = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9032 FORMAT (' DW(M+1:M+N) = ',4(1PD14.4)/(15X,4(1PD14.4)))
c9031 FORMAT (' DW(1:M)     = ',5(F11.3)/(15X,5(F11.3)))
c9032 FORMAT (' DW(M+1:M+N) = ',5(F11.3)/(15X,5(F11.3)))
      END

C**********************************************************************
      SUBROUTINE MC77N(M,N,NNZ,JCST,IRN,A,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
C
C *********************************************************************
C ***           Infinity norm  ---  The un-symmetric case           ***
C ***                  Harwell-Boeing SPARSE format                 ***
C *********************************************************************
      INTEGER M,N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCST(N+1),IRN(NNZ),IW(M),JW(N)
      REAL THRESH,ERR(2)
      REAL A(NNZ),D(M),E(N),DW(M),EW(N)

C Variables are described in MC77A/AD
C The un-symmetric matrix is stored in Harwell-Boeing sparse format.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J, K
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO

C Initialisations ...
      DO 5 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
        D(I)  = ONE
    5 CONTINUE
      DO 10 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
        E(J)  = ONE
   10 CONTINUE

C Now, we compute the initial infinity-norm of each row and column.
      DO 20 J=1,N
        DO 30 K=JCST(J),JCST(J+1)-1
          I = IRN(K)
          IF (EW(J).LT.A(K)) THEN
            EW(J) = A(K)
            JW(J) = I
          ENDIF
          IF (DW(I).LT.A(K)) THEN
            DW(I) = A(K)
            IW(I) = J
          ENDIF
   30   CONTINUE
   20 CONTINUE

      DO 40 I=1,M
        IF (IW(I).GT.0) D(I) = SQRT(DW(I))
   40 CONTINUE
      DO 45 J=1,N
        IF (JW(J).GT.0) E(J) = SQRT(EW(J))
   45 CONTINUE

      DO 50 J=1,N
        I = JW(J)
        IF (I.GT.0) THEN
          IF (IW(I).EQ.J) THEN
            IW(I) = -IW(I)
            JW(J) = -JW(J)
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 I=1,M
        IF (IW(I).GT.0)  GOTO 99
   60 CONTINUE
      DO 65 J=1,N
        IF (JW(J).GT.0)  GOTO 99
   65 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE

        DO 120 J=1,N
          IF (JW(J).GT.0) THEN
            DO 130 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              S = A(K) / (D(I)*E(J))
              IF (EW(J).LT.S) THEN
                EW(J) = S
                JW(J) = I
              ENDIF
              IF (IW(I).GT.0) THEN
                IF (DW(I).LT.S) THEN
                  DW(I) = S
                  IW(I) = J
                ENDIF
              ENDIF
  130       CONTINUE
          ELSE
            DO 140 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              IF (IW(I).GT.0) THEN
                S = A(K) / (D(I)*E(J))
                IF (DW(I).LT.S) THEN
                  DW(I) = S
                  IW(I) = J
                ENDIF
              ENDIF
  140       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 I=1,M
          IF (IW(I).GT.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).GT.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE

        DO 160 J=1,N
          I = JW(J)
          IF (I.GT.0) THEN
            IF (IW(I).EQ.J) THEN
              IW(I) = -IW(I)
              JW(J) = -JW(J)
            ENDIF
          ENDIF
  160   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in arrays
C       DW and EW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).GT.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).GT.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) )  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77O(M,N,NNZ,JCST,IRN,A,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
C
C *********************************************************************
C ***              One norm  ---  The un-symmetric case             ***
C ***                  Harwell-Boeing SPARSE format                 ***
C *********************************************************************
      INTEGER M,N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCST(N+1),IRN(NNZ),IW(M),JW(N)
      REAL THRESH,ERR(2)
      REAL A(NNZ),D(M),E(N),DW(M),EW(N)

C Variables are described in MC77A/AD
C The un-symmetric matrix is stored in Harwell-Boeing sparse format.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J, K
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO

C Initialisations ...
      DO 10 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
        D(I)  = ONE
   10 CONTINUE
      DO 15 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
        E(J)  = ONE
   15 CONTINUE

C Now, we compute the initial one-norm of each row and column.
      DO 20 J=1,N
        DO 30 K=JCST(J),JCST(J+1)-1
          IF (A(K).GT.ZERO) THEN
            I = IRN(K)
            EW(J) = EW(J) + A(K)
            IF (JW(J).EQ.0) THEN
               JW(J) = K
            ELSE
               JW(J) = -1
            ENDIF
            DW(I) = DW(I) + A(K)
            IF (IW(I).EQ.0) THEN
               IW(I) = K
            ELSE
               IW(I) = -1
            ENDIF
          ENDIF
   30   CONTINUE
   20 CONTINUE

      DO 40 K=1,M
        IF (IW(K).NE.0) D(K) = SQRT(DW(K))
   40 CONTINUE
      DO 45 K=1,N
        IF (JW(K).NE.0) E(K) = SQRT(EW(K))
   45 CONTINUE

      DO 50 J=1,N
        K = JW(J)
        IF (K.GT.0) THEN
          I = IRN(K)
          IF (IW(I).EQ.K) THEN
            IW(I) = 0
            JW(J) = 0
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 K=1,M
        IF ( (IW(K).NE.0) .OR. (JW(K).NE.0) )  GOTO 99
   60 CONTINUE
C Since we enforce M=N in the case of the one-norm, the tests can be
C done in one shot. For future relases, (M/=N) we should incorporate
C instead :
Cnew  DO 60 I=1,M
Cnew    IF (IW(I).NE.0)  GOTO 99
Cne60 CONTINUE
Cnew  DO 65 J=1,N
Cnew    IF (JW(J).NE.0)  GOTO 99
Cne65 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE

        DO 120 J=1,N
          IF (JW(J).NE.0) THEN
            DO 130 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              S = A(K) / (D(I)*E(J))
              EW(J) = EW(J) + S
              DW(I) = DW(I) + S
  130       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 I=1,M
          IF (IW(I).NE.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).NE.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in arrays
C       DW and EW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).NE.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).NE.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) ) GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77P(N,NNZ,JCST,IRN,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
C
C *********************************************************************
C ***         Infinity norm  ---  The symmetric-packed case         ***
C ***                  Harwell-Boeing SPARSE format                 ***
C *********************************************************************
      INTEGER N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCST(N+1),IRN(NNZ),IJW(N)
      REAL THRESH,ERR
      REAL A(NNZ),DE(N),DEW(N)

C Variables are described in MC77A/AD
C The lower triangular part of the symmetric matrix is stored
C in Harwell-Boeing symmetric-packed sparse format.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J, K
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR = ZERO

C Initialisations ...
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE

C Now, we compute the initial infinity-norm of each row and column.
      DO 20 J=1,N
        DO 30 K=JCST(J),JCST(J+1)-1
          I = IRN(K)
          IF (DEW(J).LT.A(K)) THEN
            DEW(J) = A(K)
            IJW(J) = I
          ENDIF
          IF (DEW(I).LT.A(K)) THEN
            DEW(I) = A(K)
            IJW(I) = J
          ENDIF
   30   CONTINUE
   20 CONTINUE

      DO 40 K=1,N
        IF (IJW(K).GT.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE

      DO 50 J=1,N
        I = IJW(J)
        IF (I.GT.0) THEN
          IF (IJW(I).EQ.J) THEN
            IJW(I) = -IJW(I)
            IF (I.NE.J) IJW(J) = -IJW(J)
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 K=1,N
        IF (IJW(K).GT.0)  GOTO 99
   60 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE

        DO 120 J=1,N
          IF (IJW(J).GT.0) THEN
            DO 130 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              S = A(K) / (DE(I)*DE(J))
              IF (DEW(J).LT.S) THEN
                DEW(J) = S
                IJW(J) = I
              ENDIF
              IF (IJW(I).GT.0) THEN
                IF (DEW(I).LT.S) THEN
                  DEW(I) = S
                  IJW(I) = J
                ENDIF
              ENDIF
  130       CONTINUE
          ELSE
            DO 140 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              IF (IJW(I).GT.0) THEN
                S = A(K) / (DE(I)*DE(J))
                IF (DEW(I).LT.S) THEN
                  DEW(I) = S
                  IJW(I) = J
                ENDIF
              ENDIF
  140       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 K=1,N
          IF (IJW(K).GT.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE

        DO 160 J=1,N
          I = IJW(J)
          IF (I.GT.0) THEN
            IF (IJW(I).EQ.J) THEN
              IJW(I) = -IJW(I)
              IF (I.NE.J) IJW(J) = -IJW(J)
            ENDIF
          ENDIF
  160   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in array
C       DEW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).GT.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77Q(N,NNZ,JCST,IRN,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
C
C *********************************************************************
C ***            One norm  ---  The symmetric-packed case           ***
C ***                  Harwell-Boeing SPARSE format                 ***
C *********************************************************************
      INTEGER N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCST(N+1),IRN(NNZ),IJW(N)
      REAL THRESH,ERR
      REAL A(NNZ),DE(N),DEW(N)

C Variables are described in MC77A/AD
C The lower triangular part of the symmetric matrix is stored
C in Harwell-Boeing symmetric-packed sparse format.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J, K
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR = ZERO

C Initialisations ...
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE

C Now, we compute the initial one-norm of each row and column.
      DO 20 J=1,N
        DO 30 K=JCST(J),JCST(J+1)-1
          IF (A(K).GT.ZERO) THEN
            DEW(J) = DEW(J) + A(K)
            IF (IJW(J).EQ.0) THEN
               IJW(J) = K
            ELSE
               IJW(J) = -1
            ENDIF
            I = IRN(K)
            IF (I.NE.J) THEN
              DEW(I) = DEW(I) + A(K)
              IJW(I) = IJW(I) + 1
              IF (IJW(I).EQ.0) THEN
                 IJW(I) = K
              ELSE
                 IJW(I) = -1
              ENDIF
            ENDIF
          ENDIF
   30   CONTINUE
   20 CONTINUE

      DO 40 K=1,N
        IF (IJW(K).NE.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE

      DO 50 J=1,N
        K = IJW(J)
        IF (K.GT.0) THEN
          I = IRN(K)
          IF (IJW(I).EQ.K) THEN
            IJW(I) = 0
            IJW(J) = 0
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 K=1,N
        IF (IJW(K).NE.0)  GOTO 99
   60 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE

        DO 120 J=1,N
          IF (IJW(J).NE.0) THEN
            DO 130 K=JCST(J),JCST(J+1)-1
              I = IRN(K)
              S = A(K) / (DE(I)*DE(J))
              DEW(J) = DEW(J) + S
              IF (I.NE.J)  DEW(I) = DEW(I) + S
  130       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 K=1,N
          IF (IJW(K).NE.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in array
C       DEW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).NE.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
C***              DRIVER FOR THE GENERAL SPARSE FORMAT              ***
C**********************************************************************
      SUBROUTINE MC77B(JOB,M,N,NNZ,IRN,JCN,A,IW,LIW,DW,LDW,
     &                  ICNTL,CNTL,INFO,RINFO)
C
C  Purpose
C  =======
C
C This subroutine computes scaling factors D(i), i=1..M, and E(j),
C j=1..N, for an MxN sparse matrix A = {a_ij} so that the scaled matrix
C D^{-1}*A*E^{-1} has both rows and columns norms all equal or close to
C 1. The rectangular case (M/=N) is, for the moment, only allowed for
C equilibration in the infinity norm, a particular case for which the
C convergence is ensured in all cases and trivial to monitor. See [1].
C
C  Parameters
C  ==========
C
      INTEGER LICNTL, LCNTL, LINFO, LRINFO
      PARAMETER ( LICNTL=10, LCNTL=10, LINFO=10, LRINFO=10 )
      INTEGER ICNTL(LICNTL),INFO(LINFO)
      REAL CNTL(LCNTL),RINFO(LRINFO)
C
      INTEGER JOB,M,N,NNZ,LIW,LDW
      INTEGER JCN(NNZ),IRN(NNZ),IW(LIW)
      REAL A(NNZ),DW(LDW)
C
C JOB is an INTEGER variable which must be set by the user to
C control the action. It is not altered by the subroutine.
C Possible values for JOB are:
C   0 Equilibrate the infinity norm of rows and columns in matrix A.
C   1 Equilibrate the one norm of rows and columns in matrix A.
C   p Equilibrate the p-th norm (p>=2) of rows and columns in matrix A.
C  -1 Equilibrate the p-th norm of rows and columns in matrix A, with
C     a strictly positive REAL parameter p given in CNTL(2).
C   Restriction: JOB >= -1.
C
C M is an INTEGER variable which must be set by the user to the
C   number of rows in matrix A. It is not altered by the subroutine.
C   Restriction: M >= 1.
C
C N is an INTEGER variable which must be set by the user to the
C   number of columns in matrix A. It is not altered by the subroutine.
C   Restriction: N >= 1.
C            If ICNTL(6) /= 0 (symmetric matrix), N must be equal to M.
C
C NNZ is an INTEGER variable which must be set by the user to the number
C   of entries in the matrix. It is not altered by the subroutine.
C   Restriction: NNZ >= 1.
C
C IRN is an INTEGER array of length NNZ.
C   IRN(K), K=1..NNZ, must be set by the user to hold the row indices
C   of the entries in the matrix.
C   The array IRN is not altered by the subroutine.
C   Restriction:
C     Out-of-range indices and duplicates are not allowed.
C
C JCN is an INTEGER array of NNZ.
C   JCN(J), J=1..NNZ, must be set by the user to hold the column indices
C   of the entries in the matrix.
C   The array JCN is not altered by the subroutine.
C   Restriction:
C     Out-of-range indices and duplicates are not allowed.
C
C A is a REAL (REAL in the D-version) array of length NNZ.
C   The user must set A(K), K=1..NNZ, to the numerical value of the
C   entry that corresponds to IRN(K).
C   It is not altered by the subroutine.
C
C LIW is an INTEGER variable that must be set by the user to
C   the length of array IW. It is not altered by the subroutine.
C   Restriction:
C     If ICNTL(6) == 0:  LIW >= M+N
C     If ICNTL(6) /= 0:  LIW >= M
C
C IW is an INTEGER array of length LIW that is used for workspace.
C
C LDW is an INTEGER variable that must be set by the user to the
C   length of array DW. It is not altered by the subroutine.
C   Restriction:
C     If ICNTL(5) = 0 and ICNTL(6) == 0:  LDW >= NNZ + 2*(M+N)
C     If ICNTL(5) = 1 and ICNTL(6) == 0:  LDW >= 2*(M+N)
C     If ICNTL(5) = 0 and ICNTL(6) /= 0:  LDW >= NNZ + 2*M
C     If ICNTL(5) = 1 and ICNTL(6) /= 0:  LDW >= 2*M
C
C DW is a REAL (REAL in the D-version) array of length LDW
C   that need not be set by the user.
C   On return, DW(i) contains d_i, i=1..M, the diagonal entries in
C   the row-scaling matrix D, and when the input matrix A is
C   unsymmetric, DW(M+j) contains e_j, j=1..N, the diagonal entries in
C   the column-scaling matrix E.
C
C ICNTL is an INTEGER array of length 10. Its components control the
C   output of MC77A/AD and must be set by the user before calling
C   MC77A/AD. They are not altered by the subroutine.
C   See MC77I/ID for details.
C
C CNTL is a REAL (REAL in the D-version) array of length 10.
C   Its components control the output of MC77A/AD and must be set by the
C   user before calling MC77A/AD. They are not altered by the
C   subroutine.
C   See MC77I/ID for details.
C
C INFO is an INTEGER array of length 10 which need not be set by the
C   user. INFO(1) is set non-negative to indicate success. A negative
C   value is returned if an error occurred, a positive value if a
C   warning occurred. INFO(2) holds further information on the error.
C   On exit from the subroutine, INFO(1) will take one of the
C   following values:
C    0 : successful entry (for structurally nonsingular matrix).
C   +1 : The maximum number of iterations in ICNTL(7) has been reached.
C   -1 : M < 1.  Value of M held in INFO(2).
C   -2 : N < 1.  Value of N held in INFO(2).
C   -3 : ICNTL(6) /= 0 (input matrix is symmetric) and N /= M.
C        Value of (N-M) held in INFO(2).
C   -4 : NNZ < 1. Value of NNZ held in INFO(2).
C   -5 : the defined length LIW violates the restriction on LIW.
C        Value of LIW required given by INFO(2).
C   -6 : the defined length LDW violates the restriction on LDW.
C        Value of LDW required given by INFO(2).
C   -7 : entries are found whose row indices are out of range.
C        INFO(2) contains an index pointing to a value in
C        IRN and JCN in which such an entry is found.
C   -8 : repeated entries are found.
C        INFO(2) contains an index pointing to a value in
C        IRN and JCN in which such an entry is found.
C   -9 : An entry corresponding to an element in the upper triangular
C        part of the matrix is found in the input data (ICNTL(6)} \= 0).
C  -10 : ICNTL(7) is out of range and INFO(2) contains its value.
C  -11 : CNTL(2) is out of range.
C  -12 : JOB < -1.  Value of JOB held in INFO(2).
C  -13 : JOB /= 0 and N /= M. This is a restriction of the algorithm in
C        its present stage, e.g. for rectangular matrices, only the
C        scaling in infinity norm is allowed. Value of JOB held in
C        INFO(2).
C   INFO(3) returns the number of iterations performed.
C   INFO(4) to INFO(10) are not currently used and are set to zero by
C        the routine.
C
C RINFO is a REAL (REAL in the D-version) array of length 10
C which need not be set by the user.
C   When CNTL(1) is strictly positive, RINFO(1) will contain on exit the
C   maximum distance of all row norms to 1, and RINFO(2) will contain on
C   exit the maximum distance of all column norms to 1.
C   If ICNTL(6) /= 0 (indicating that the input matrix is symmetric),
C   RINFO(2) is not used since the row-scaling equals the column scaling
C   in this case.
C   RINFO(3) to RINFO(10) are not currently used and are set to zero.
C
C References:
C  [1]  D. R. Ruiz, (2001),
C       "A scaling algorithm to equilibrate both rows and columns norms
C       in matrices",
C       Technical Report RAL-TR-2001-034, RAL, Oxfordshire, England,
C             and Report RT/APO/01/4, ENSEEIHT-IRIT, Toulouse, France.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      REAL THRESH, DP
      INTEGER I, J, K, MAXIT, CHECK, SETUP
C External routines and functions
      EXTERNAL MC77R,MC77S,MC77T,MC77U
C Intrinsic functions
      INTRINSIC ABS,MAX

C Reset informational output values
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = 0
      RINFO(1) = ZERO
      RINFO(2) = ZERO
C Check value of JOB
      IF (JOB.LT.-1) THEN
        INFO(1) = -12
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
        GO TO 99
      ENDIF
C Check bad values in ICNTL
Cnew  IF (ICNTL(7).LT.0) THEN
      IF (ICNTL(7).LE.0) THEN
        INFO(1) = -10
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9002) INFO(1),7,ICNTL(7)
        GO TO 99
      ENDIF
C Check bad values in CNTL
      IF ((JOB.EQ.-1) .AND. (CNTL(2).LT.ONE)) THEN
        INFO(1) = -11
        INFO(2) = 2
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9008) INFO(1),2,CNTL(2)
        GO TO 99
      ENDIF
C Check value of M
      IF (M.LT.1) THEN
        INFO(1) = -1
        INFO(2) = M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'M',M
        GO TO 99
      ENDIF
C Check value of N
      IF (N.LT.1) THEN
        INFO(1) = -2
        INFO(2) = N
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
        GO TO 99
      ENDIF
C Symmetric case: M must be equal to N
      IF ( (ICNTL(6).NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -3
        INFO(2) = N-M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9003) INFO(1),(N-M)
        GO TO 99
      ENDIF
C For rectangular matrices, only scaling in infinity norm is allowed
      IF ( (JOB.NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -13
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9009) INFO(1),JOB,M,N
        GO TO 99
      ENDIF
C Check value of NNZ
      IF (NNZ.LT.1) THEN
        INFO(1) = -4
        INFO(2) = NNZ
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'NNZ',NNZ
        GO TO 99
      ENDIF
C Check LIW
      K = M
      IF (ICNTL(6).EQ.0)  K = M+N
      IF (LIW.LT.K) THEN
        INFO(1) = -5
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
        GO TO 99
      ENDIF
C Check LDW
      K = 2*M
      IF (ICNTL(6).EQ.0)  K = 2*(M+N)
      IF ( (ICNTL(5).EQ.0) .OR. (JOB.GE.2) .OR. (JOB.EQ.-1) )
     &   K = K + NNZ
      IF (LDW.LT.K) THEN
        INFO(1) = -6
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
        GO TO 99
      ENDIF
C Check row indices. Use IW(1:M) as workspace
C     General sparse format :
      IF (ICNTL(4).EQ.0) THEN
        DO 6 K = 1,NNZ
          I = IRN(K)
          J = JCN(K)
C Check for row indices that are out of range
          IF (I.LT.1 .OR. I.GT.M .OR. J.LT.1 .OR. J.GT.N) THEN
            INFO(1) = -7
            INFO(2) = K
            IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),K,I,J
            GO TO 99
          ENDIF
          IF (ICNTL(6).NE.0 .AND. I.LT.J) THEN
            INFO(1) = -9
            INFO(2) = K
            IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),K,I,J
            GO TO 99
          ENDIF
    6   CONTINUE
C Check for repeated row indices within a column
        DO 7 I = 1,M
          IW(I) = 0
    7   CONTINUE
        DO 9 J = 1,N
          DO 8 K = 1,NNZ
          IF (JCN(K).EQ.J) THEN
            I = IRN(K)
            IF (IW(I).EQ.J) THEN
              INFO(1) = -8
              INFO(2) = K
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9007) INFO(1),K,I,J
              GO TO 99
            ELSE
              IW(I) = J
            ENDIF
          ENDIF
    8     CONTINUE
    9   CONTINUE
      ENDIF

C Print diagnostics on input
C ICNTL(9) could be set by the user to monitor the level of diagnostics
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9020) JOB,M,N,NNZ,ICNTL(7),CNTL(1)
C       General sparse format :
        IF (ICNTL(9).LT.1) THEN
          WRITE(ICNTL(3),9022) (JCN(J),J=1,NNZ)
          WRITE(ICNTL(3),9023) (IRN(J),J=1,NNZ)
          WRITE(ICNTL(3),9024) (A(J),J=1,NNZ)
        ENDIF
      ENDIF

C Set components of INFO to zero
      DO 10 I=1,10
        INFO(I)  = 0
        RINFO(I) = ZERO
   10 CONTINUE

C Set the convergence threshold and the maximum number of iterations
      THRESH = MAX( ZERO, CNTL(1) )
      MAXIT  = ICNTL(7)
C Check for the particular inconsistency between MAXIT and THRESH.
Cnew  IF ( (MAXIT.LE.0) .AND. (CNTL(1).LE.ZERO) )  THEN
Cnew     INFO(1) = -14
Cnew     IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9010) INFO(1),MAXIT,CNTL(1)
Cnew     GO TO 99
Cnew  ENDIF

C   ICNTL(8) could be set by the user to indicate the frequency at which
C   convergence should be checked when CNTL(1) is strictly positive. The
C   default value used here 1 (e.g. check convergence every iteration).
C     CHECK  = ICNTL(8)
      CHECK  = 1
      IF (CNTL(1).LE.ZERO)  CHECK = 0

C Prepare temporary data in working arrays as apropriate
      K = 2*M
      IF (ICNTL(6).EQ.0)  K = 2*(M+N)
      SETUP = 0
      IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
        SETUP = 1
C       Copy ABS(A(J))**JOB into DW(K+J), J=1..NNZ
        IF (JOB.EQ.-1) THEN
          DP = CNTL(2)
        ELSE
          DP = REAL(JOB)
        ENDIF
        DO 20 J = 1,NNZ
          DW(K+J) = ABS(A(J))**DP
   20   CONTINUE
C       Reset the Threshold to take into account the Hadamard p-th power
        THRESH = ONE - (ABS(ONE-THRESH))**DP
      ELSE
        IF (ICNTL(5).EQ.0) THEN
          SETUP = 1
C         Copy ABS(A(J)) into DW(K+J), J=1..NNZ
          DO 30 J = 1,NNZ
            DW(K+J) = ABS(A(J))
   30     CONTINUE
        ENDIF
      ENDIF

C Begin the computations (input matrix in general SPARSE format) ...
C     We consider first the un-symmetric case :
      IF (ICNTL(6).EQ.0)  THEN
C
C       Equilibrate matrix A in the infinity norm
        IF (JOB.EQ.0) THEN
C         IW(1:M+N), DW(M+N:2*(M+N)) are workspaces
          IF (SETUP.EQ.1) THEN
            CALL MC77R(M,N,NNZ,JCN,IRN,DW(2*(M+N)+1),DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77R(M,N,NNZ,JCN,IRN,A,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
C
C       Equilibrate matrix A in the p-th norm, 1 <= p < +infinity
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77S(M,N,NNZ,JCN,IRN,DW(2*(M+N)+1),DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77S(M,N,NNZ,JCN,IRN,A,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
C         Reset the scaling factors to their Hadamard p-th root, if
C         required and re-compute the actual value of the final error
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 40 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+N+I)**DP) )
   40       CONTINUE
            RINFO(2) = ZERO
            DO 50 J = 1,N
              DW(M+J) = DW(M+J)**DP
              IF (IW(M+J).NE.0)
     &           RINFO(2) = MAX( RINFO(2), ABS(ONE-DW(2*M+N+J)**DP) )
   50       CONTINUE
          ENDIF
        ENDIF
C
C     We treat then the symmetric (packed) case :
      ELSE
C
C       Equilibrate matrix A in the infinity norm
        IF (JOB.EQ.0) THEN
C         IW(1:M), DW(M:2*M) are workspaces
          IF (SETUP.EQ.1) THEN
            CALL MC77T(M,NNZ,JCN,IRN,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77T(M,NNZ,JCN,IRN,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
C
C       Equilibrate matrix A in the p-th norm, 1 <= p < +infinity
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77U(M,NNZ,JCN,IRN,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77U(M,NNZ,JCN,IRN,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
C         Reset the scaling factors to their Hadamard p-th root, if
C         required and re-compute the actual value of the final error
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 60 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+I)**DP) )
   60       CONTINUE
          ENDIF
        ENDIF
        RINFO(2) = RINFO(1)
C
      ENDIF
C     End of the un-symmetric/symmetric IF statement
C     The general sparse case has been treated


C Print diagnostics on output
C ICNTL(9) could be set by the user to monitor the level of diagnostics
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9030) (INFO(J),J=1,3),(RINFO(J),J=1,2)
        IF (ICNTL(9).LT.2) THEN
          WRITE(ICNTL(3),9031) (DW(I),I=1,M)
          IF (ICNTL(6).EQ.0)
     &      WRITE(ICNTL(3),9032) (DW(M+J),J=1,N)
        ENDIF
      ENDIF

C Return from subroutine.
   99 CONTINUE
      RETURN

 9001 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3,
     &        ' because ',(A),' = ',I10)
 9002 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Bad input control flag.',
     &        ' Value of ICNTL(',I1,') = ',I8)
 9003 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Input matrix is symmetric and N /= M.',
     &        ' Value of (N-M) = ',I8)
 9004 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Entry ',I8,
     &        ' has invalid row index ',I8, ' or column index ',I8)
 9007 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Duplicate entry ',I8, '   with row index ',I8/
     &        '                                 and column index ',I8)
 9008 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Bad input REAL control parameter.',
     &        ' Value of CNTL(',I1,') = ',1PD14.4)
 9009 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
     &        '        Only scaling in infinity norm is allowed',
     &        ' for rectangular matrices'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8)
C9010 FORMAT (' ****** Error in MC77B/BD. INFO(1) = ',I3/
Cnew &        '        Algorithm will loop indefinitely !!!',
Cnew &        '     Check for inconsistency'/
Cnew &        ' Maximum number of iterations = ',I8/
Cnew &        ' Threshold for convergence    = ',1PD14.4)
 9020 FORMAT (' ****** Input parameters for MC77B/BD:'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8/' NNZ = ',I8/
     &        ' Max n.b. Iters.  = ',I8/
     &        ' Cvgce. Threshold = ',1PD14.4)
 9022 FORMAT (' JCN(1:NNZ)  = ',8I8/(15X,8I8))
 9023 FORMAT (' IRN(1:NNZ)  = ',8I8/(15X,8I8))
 9024 FORMAT (' A(1:NNZ)    = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for MC77B/BD:'/
     &        ' INFO(1:3)   = ',3(2X,I8)/
     &        ' RINFO(1:2)  = ',2(1PD14.4))
 9031 FORMAT (' DW(1:M)     = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9032 FORMAT (' DW(M+1:M+N) = ',4(1PD14.4)/(15X,4(1PD14.4)))
c9031 FORMAT (' DW(1:M)     = ',5(F11.3)/(15X,5(F11.3)))
c9032 FORMAT (' DW(M+1:M+N) = ',5(F11.3)/(15X,5(F11.3)))
      END

C**********************************************************************
      SUBROUTINE MC77R(M,N,NNZ,JCN,IRN,A,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
C
C *********************************************************************
C ***           Infinity norm  ---  The un-symmetric case           ***
C ***                     General SPARSE format                     ***
C *********************************************************************
      INTEGER M,N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCN(NNZ),IRN(NNZ),IW(M),JW(N)
      REAL THRESH,ERR(2)
      REAL A(NNZ),D(M),E(N),DW(M),EW(N)

C Variables are described in MC77A/AD
C The un-symmetric matrix is stored in general sparse format.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J, K
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO

C Initialisations ...
      DO 5 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
        D(I)  = ONE
    5 CONTINUE
      DO 10 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
        E(J)  = ONE
   10 CONTINUE

C Now, we compute the initial infinity-norm of each row and column.
      DO 20 K=1,NNZ
        I = IRN(K)
        J = JCN(K)
        IF (EW(J).LT.A(K)) THEN
          EW(J) = A(K)
          JW(J) = I
        ENDIF
        IF (DW(I).LT.A(K)) THEN
          DW(I) = A(K)
          IW(I) = J
        ENDIF
   20 CONTINUE

      DO 40 I=1,M
        IF (IW(I).GT.0) D(I) = SQRT(DW(I))
   40 CONTINUE
      DO 45 J=1,N
        IF (JW(J).GT.0) E(J) = SQRT(EW(J))
   45 CONTINUE

      DO 50 J=1,N
        I = JW(J)
        IF (I.GT.0) THEN
          IF (IW(I).EQ.J) THEN
            IW(I) = -IW(I)
            JW(J) = -JW(J)
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 I=1,M
        IF (IW(I).GT.0)  GOTO 99
   60 CONTINUE
      DO 65 J=1,N
        IF (JW(J).GT.0)  GOTO 99
   65 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE

        DO 120 K=1,NNZ
          I = IRN(K)
          J = JCN(K)
          IF (JW(J).GT.0) THEN
            S = A(K) / (D(I)*E(J))
            IF (EW(J).LT.S) THEN
              EW(J) = S
              JW(J) = I
            ENDIF
            IF (IW(I).GT.0) THEN
              IF (DW(I).LT.S) THEN
                DW(I) = S
                IW(I) = J
              ENDIF
            ENDIF
          ELSE
            IF (IW(I).GT.0) THEN
              S = A(K) / (D(I)*E(J))
              IF (DW(I).LT.S) THEN
                DW(I) = S
                IW(I) = J
              ENDIF
            ENDIF
          ENDIF
  120   CONTINUE

        DO 150 I=1,M
          IF (IW(I).GT.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).GT.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE

        DO 160 J=1,N
          I = JW(J)
          IF (I.GT.0) THEN
            IF (IW(I).EQ.J) THEN
              IW(I) = -IW(I)
              JW(J) = -JW(J)
            ENDIF
          ENDIF
  160   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in arrays
C       DW and EW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).GT.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).GT.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) )  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77S(M,N,NNZ,JCN,IRN,A,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
C
C *********************************************************************
C ***              One norm  ---  The un-symmetric case             ***
C ***                     General SPARSE format                     ***
C *********************************************************************
      INTEGER M,N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCN(NNZ),IRN(NNZ),IW(M),JW(N)
      REAL THRESH,ERR(2)
      REAL A(NNZ),D(M),E(N),DW(M),EW(N)

C Variables are described in MC77A/AD
C The un-symmetric matrix is stored in general sparse format.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J, K
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO

C Initialisations ...
      DO 10 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
        D(I)  = ONE
   10 CONTINUE
      DO 15 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
        E(J)  = ONE
   15 CONTINUE

C Now, we compute the initial one-norm of each row and column.
      DO 20 K=1,NNZ
        IF (A(K).GT.ZERO) THEN
          I = IRN(K)
          J = JCN(K)
          EW(J) = EW(J) + A(K)
          IF (JW(J).EQ.0) THEN
             JW(J) = K
          ELSE
             JW(J) = -1
          ENDIF
          DW(I) = DW(I) + A(K)
          IF (IW(I).EQ.0) THEN
             IW(I) = K
          ELSE
             IW(I) = -1
          ENDIF
        ENDIF
   20 CONTINUE

      DO 40 I=1,M
        IF (IW(I).NE.0) D(I) = SQRT(DW(I))
   40 CONTINUE
      DO 45 J=1,N
        IF (JW(J).NE.0) E(J) = SQRT(EW(J))
   45 CONTINUE

      DO 50 J=1,N
        K = JW(J)
        IF (K.GT.0) THEN
          I = IRN(K)
          IF (IW(I).EQ.K) THEN
            IW(I) = 0
            JW(J) = 0
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 K=1,M
        IF ( (IW(K).NE.0) .OR. (JW(K).NE.0) )  GOTO 99
   60 CONTINUE
C Since we enforce M=N in the case of the one-norm, the tests can be
C done in one shot. For future relases, (M/=N) we should incorporate
C instead :
Cnew  DO 60 I=1,M
Cnew    IF (IW(I).NE.0)  GOTO 99
Cne60 CONTINUE
Cnew  DO 65 J=1,N
Cnew    IF (JW(J).NE.0)  GOTO 99
Cne65 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE

        DO 120 K=1,NNZ
          J = JCN(K)
          I = IRN(K)
          S = A(K) / (D(I)*E(J))
          EW(J) = EW(J) + S
          DW(I) = DW(I) + S
  120   CONTINUE

        DO 150 I=1,M
          IF (IW(I).NE.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).NE.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in arrays
C       DW and EW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).NE.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).NE.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) ) GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77T(N,NNZ,JCN,IRN,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
C
C *********************************************************************
C ***         Infinity norm  ---  The symmetric-packed case         ***
C ***                     General SPARSE format                     ***
C *********************************************************************
      INTEGER N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCN(NNZ),IRN(NNZ),IJW(N)
      REAL THRESH,ERR
      REAL A(NNZ),DE(N),DEW(N)

C Variables are described in MC77A/AD
C The lower triangular part of the symmetric matrix is stored
C in general symmetric-packed sparse format.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J, K
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR = ZERO

C Initialisations ...
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE

C Now, we compute the initial infinity-norm of each row and column.
      DO 20 K=1,NNZ
        I = IRN(K)
        J = JCN(K)
        IF (DEW(J).LT.A(K)) THEN
          DEW(J) = A(K)
          IJW(J) = I
        ENDIF
        IF (DEW(I).LT.A(K)) THEN
          DEW(I) = A(K)
          IJW(I) = J
        ENDIF
   20 CONTINUE

      DO 40 K=1,N
        IF (IJW(K).GT.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE

      DO 50 J=1,N
        I = IJW(J)
        IF (I.GT.0) THEN
          IF (IJW(I).EQ.J) THEN
            IJW(I) = -IJW(I)
            IF (I.NE.J) IJW(J) = -IJW(J)
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 K=1,N
        IF (IJW(K).GT.0)  GOTO 99
   60 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE

        DO 120 K=1,NNZ
          I = IRN(K)
          J = JCN(K)
          IF (IJW(J).GT.0) THEN
            S = A(K) / (DE(I)*DE(J))
            IF (DEW(J).LT.S) THEN
              DEW(J) = S
              IJW(J) = I
            ENDIF
            IF (IJW(I).GT.0) THEN
              IF (DEW(I).LT.S) THEN
                DEW(I) = S
                IJW(I) = J
              ENDIF
            ENDIF
          ELSE
            IF (IJW(I).GT.0) THEN
              S = A(K) / (DE(I)*DE(J))
              IF (DEW(I).LT.S) THEN
                DEW(I) = S
                IJW(I) = J
              ENDIF
            ENDIF
          ENDIF
  120   CONTINUE

        DO 150 K=1,N
          IF (IJW(K).GT.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE

        DO 160 J=1,N
          I = IJW(J)
          IF (I.GT.0) THEN
            IF (IJW(I).EQ.J) THEN
              IJW(I) = -IJW(I)
              IF (I.NE.J) IJW(J) = -IJW(J)
            ENDIF
          ENDIF
  160   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in array
C       DEW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).GT.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77U(N,NNZ,JCN,IRN,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
C
C *********************************************************************
C ***            One norm  ---  The symmetric-packed case           ***
C ***                     General SPARSE format                     ***
C *********************************************************************
      INTEGER N,NNZ,MAXIT,NITER,CHECK,INFO
      INTEGER JCN(NNZ),IRN(NNZ),IJW(N)
      REAL THRESH,ERR
      REAL A(NNZ),DE(N),DEW(N)

C Variables are described in MC77A/AD
C The lower triangular part of the symmetric matrix is stored
C in general symmetric-packed sparse format.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J, K
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR = ZERO

C Initialisations ...
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE

C Now, we compute the initial one-norm of each row and column.
      DO 20 K=1,NNZ
        IF (A(K).GT.ZERO) THEN
          J = JCN(K)
          DEW(J) = DEW(J) + A(K)
          IF (IJW(J).EQ.0) THEN
             IJW(J) = K
          ELSE
             IJW(J) = -1
          ENDIF
          I = IRN(K)
          IF (I.NE.J) THEN
            DEW(I) = DEW(I) + A(K)
            IF (IJW(I).EQ.0) THEN
               IJW(I) = K
            ELSE
               IJW(I) = -1
            ENDIF
          ENDIF
        ENDIF
   20 CONTINUE

      DO 40 K=1,N
        IF (IJW(K).NE.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE

      DO 50 J=1,N
        K = IJW(J)
        IF (K.GT.0) THEN
          I = IRN(K)
          IF (IJW(I).EQ.K) THEN
            IJW(I) = 0
            IJW(J) = 0
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 K=1,N
        IF (IJW(K).NE.0)  GOTO 99
   60 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE

        DO 120 K=1,NNZ
          J = JCN(K)
          I = IRN(K)
          S = A(K) / (DE(I)*DE(J))
          DEW(J) = DEW(J) + S
          IF (I.NE.J)  DEW(I) = DEW(I) + S
  120   CONTINUE

        DO 150 K=1,N
          IF (IJW(K).NE.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in array
C       DEW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).NE.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
C***             DRIVER FOR THE CASE OF DENSE MATRICES              ***
C**********************************************************************
      SUBROUTINE MC77C(JOB,M,N,A,LDA,IW,LIW,DW,LDW,
     &                  ICNTL,CNTL,INFO,RINFO)
C
C  Purpose
C  =======
C
C This subroutine computes scaling factors D(i), i=1..M, and E(j),
C j=1..N, for an MxN sparse matrix A = {a_ij} so that the scaled matrix
C D^{-1}*A*E^{-1} has both rows and columns norms all equal or close to
C 1. The rectangular case (M/=N) is, for the moment, only allowed for
C equilibration in the infinity norm, a particular case for which the
C convergence is ensured in all cases and trivial to monitor. See [1].
C
C  Parameters
C  ==========
C
      INTEGER LICNTL, LCNTL, LINFO, LRINFO
      PARAMETER ( LICNTL=10, LCNTL=10, LINFO=10, LRINFO=10 )
      INTEGER ICNTL(LICNTL),INFO(LINFO)
      REAL CNTL(LCNTL),RINFO(LRINFO)
C
      INTEGER JOB,M,N,LDA,LIW,LDW
      INTEGER IW(LIW)
      REAL A(LDA,*),DW(LDW)
C
C JOB is an INTEGER variable which must be set by the user to
C control the action. It is not altered by the subroutine.
C Possible values for JOB are:
C   0 Equilibrate the infinity norm of rows and columns in matrix A.
C   1 Equilibrate the one norm of rows and columns in matrix A.
C   p Equilibrate the p-th norm (p>=2) of rows and columns in matrix A.
C  -1 Equilibrate the p-th norm of rows and columns in matrix A, with
C     a strictly positive REAL parameter p given in CNTL(2).
C   Restriction: JOB >= -1.
C
C M is an INTEGER variable which must be set by the user to the
C   number of rows in matrix A. It is not altered by the subroutine.
C   Restriction: M >= 1.
C
C N is an INTEGER variable which must be set by the user to the
C   number of columns in matrix A. It is not altered by the subroutine.
C   Restriction: N >= 1.
C             If ICNTL(6) /= 0 (symmetric matrix), N must be equal to M.
C
C A is a REAL (REAL in the D-version) array containing
C   the numerical values A(i,j), i=1..M, j=1..N, of the input matrix A.
C   It is not altered by the subroutine.
C
C LIW is an INTEGER variable that must be set by the user to
C   the length of array IW. It is not altered by the subroutine.
C   Restriction:
C     If ICNTL(6) == 0:  LIW >= M+N
C     If ICNTL(6) /= 0:  LIW >= M
C
C IW is an INTEGER array of length LIW that is used for workspace.
C
C LDW is an INTEGER variable that must be set by the user to the
C   length of array DW. It is not altered by the subroutine.
C   Restriction:
C     If ICNTL(5) = 0 and ICNTL(6) == 0:  LDW >= (LDA*N) + 2*(M+N)
C     If ICNTL(5) = 1 and ICNTL(6) == 0:  LDW >= 2*(M+N)
C     If ICNTL(5) = 0 and ICNTL(6) /= 0:  LDW >= M*(M+1)/2 + 2*M
C     If ICNTL(5) = 1 and ICNTL(6) /= 0:  LDW >= 2*M
C
C DW is a REAL (REAL in the D-version) array of length LDW
C   that need not be set by the user.
C   On return, DW(i) contains d_i, i=1..M, the diagonal entries in
C   the row-scaling matrix D, and when the input matrix A is
C   unsymmetric, DW(M+j) contains e_j, j=1..N, the diagonal entries in
C   the column-scaling matrix E.
C
C ICNTL is an INTEGER array of length 10. Its components control the
C   output of MC77A/AD and must be set by the user before calling
C   MC77A/AD. They are not altered by the subroutine.
C   See MC77I/ID for details.
C
C CNTL is a REAL (REAL in the D-version) array of length 10.
C   Its components control the output of MC77A/AD and must be set by the
C   user before calling MC77A/AD. They are not altered by the
C   subroutine.
C   See MC77I/ID for details.
C
C INFO is an INTEGER array of length 10 which need not be set by the
C   user. INFO(1) is set non-negative to indicate success. A negative
C   value is returned if an error occurred, a positive value if a
C   warning occurred. INFO(2) holds further information on the error.
C   On exit from the subroutine, INFO(1) will take one of the
C   following values:
C    0 : successful entry (for structurally nonsingular matrix).
C   +1 : The maximum number of iterations in ICNTL(7) has been reached.
C   -1 : M < 1.  Value of M held in INFO(2).
C   -2 : N < 1.  Value of N held in INFO(2).
C   -3 : ICNTL(6) /= 0 (input matrix is symmetric) and N /= M.
C        Value of (N-M) held in INFO(2).
C   -5 : the defined length LIW violates the restriction on LIW.
C        Value of LIW required given by INFO(2).
C   -6 : the defined length LDW violates the restriction on LDW.
C        Value of LDW required given by INFO(2).
C  -10 : ICNTL(7) is out of range and INFO(2) contains its value.
C  -11 : CNTL(2) is out of range.
C  -12 : JOB < -1.  Value of JOB held in INFO(2).
C  -13 : JOB /= 0 and N /= M. This is a restriction of the algorithm in
C        its present stage, e.g. for rectangular matrices, only the
C        scaling in infinity norm is allowed. Value of JOB held in
C        INFO(2).
C   INFO(3) returns the number of iterations performed.
C   INFO(4) to INFO(10) are not currently used and are set to zero by
C        the routine.
C
C RINFO is a REAL (REAL in the D-version) array of length 10
C which need not be set by the user.
C   When CNTL(1) is strictly positive, RINFO(1) will contain on exit the
C   maximum distance of all row norms to 1, and RINFO(2) will contain on
C   exit the maximum distance of all column norms to 1.
C   If ICNTL(6) /= 0 (indicating that the input matrix is symmetric),
C   RINFO(2) is not used since the row-scaling equals the column scaling
C   in this case.
C   RINFO(3) to RINFO(10) are not currently used and are set to zero.
C
C References:
C  [1]  D. R. Ruiz, (2001),
C       "A scaling algorithm to equilibrate both rows and columns norms
C       in matrices",
C       Technical Report RAL-TR-2001-034, RAL, Oxfordshire, England,
C             and Report RT/APO/01/4, ENSEEIHT-IRIT, Toulouse, France.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      REAL THRESH, DP
      INTEGER I, J, K, NNZ, MAXIT, CHECK, SETUP
C External routines and functions
      EXTERNAL MC77J,MC77K,MC77L,MC77M
C Intrinsic functions
      INTRINSIC ABS,MAX

C Reset informational output values
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = 0
      RINFO(1) = ZERO
      RINFO(2) = ZERO
C Check value of JOB
      IF (JOB.LT.-1) THEN
        INFO(1) = -12
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
        GO TO 99
      ENDIF
C Check bad values in ICNTL
Cnew  IF (ICNTL(7).LT.0) THEN
      IF (ICNTL(7).LE.0) THEN
        INFO(1) = -10
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9002) INFO(1),7,ICNTL(7)
        GO TO 99
      ENDIF
C Check bad values in CNTL
      IF ((JOB.EQ.-1) .AND. (CNTL(2).LT.ONE)) THEN
        INFO(1) = -11
        INFO(2) = 2
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9008) INFO(1),2,CNTL(2)
        GO TO 99
      ENDIF
C Check value of M
      IF (M.LT.1) THEN
        INFO(1) = -1
        INFO(2) = M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'M',M
        GO TO 99
      ENDIF
C Check value of N
      IF (N.LT.1) THEN
        INFO(1) = -2
        INFO(2) = N
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
        GO TO 99
      ENDIF
C Symmetric case: M must be equal to N
      IF ( (ICNTL(6).NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -3
        INFO(2) = N-M
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9003) INFO(1),(N-M)
        GO TO 99
      ENDIF
C For rectangular matrices, only scaling in infinity norm is allowed
      IF ( (JOB.NE.0) .AND. (N.NE.M) ) THEN
        INFO(1) = -13
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9009) INFO(1),JOB,M,N
        GO TO 99
      ENDIF
C Check LDA
      IF (ICNTL(6).EQ.0) THEN
        IF (M.GT.LDA) THEN
          INFO(1) = -4
          INFO(2) = LDA-M
          IF (ICNTL(1).GE.0)
     &      WRITE(ICNTL(1),9001) INFO(1),'LDA < M; (LDA-M)',INFO(2)
          GO TO 99
        ENDIF
      ELSE
        IF (((M*(M+1))/2).GT.LDA) THEN
          INFO(1) = -4
          INFO(2) = LDA-((M*(M+1))/2)
          IF (ICNTL(1).GE.0)
     &      WRITE(ICNTL(1),9001) INFO(1),
     &        'LDA < M*(M+1)/2; (LDA-(M*(M+1)/2))',INFO(2)
          GO TO 99
        ENDIF
      ENDIF
C Check LIW
      K = M
      IF (ICNTL(6).EQ.0)  K = M+N
      IF (LIW.LT.K) THEN
        INFO(1) = -5
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
        GO TO 99
      ENDIF
C Check LDW
      K = 2*M
      NNZ = (M*(M+1))/2
      IF (ICNTL(6).EQ.0)  THEN
        K = 2*(M+N)
        NNZ = LDA*N
      ENDIF
      IF ( (ICNTL(5).EQ.0) .OR. (JOB.GE.2) .OR. (JOB.EQ.-1) )
     &   K = K + NNZ
      IF (LDW.LT.K) THEN
        INFO(1) = -6
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
        GO TO 99
      ENDIF

C Print diagnostics on input
C ICNTL(9) could be set by the user to monitor the level of diagnostics
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9020) JOB,M,N,ICNTL(7),CNTL(1)
        IF (ICNTL(9).LT.1) THEN
          K = 1
          DO 5 I=1,M
C           General case :
            IF (ICNTL(6).EQ.0)
     &        WRITE(ICNTL(3),9021) I,(A(I,J),J=1,N)
C           Dense symmetric packed format :
            IF (ICNTL(6).NE.0) THEN
              WRITE(ICNTL(3),9022) I,(A(J,1),J=K,K+M-I)
              K = K+M-I+1
            ENDIF
    5     CONTINUE
        ENDIF
      ENDIF

C Set components of INFO to zero
      DO 10 I=1,10
        INFO(I)  = 0
        RINFO(I) = ZERO
   10 CONTINUE

C Set the convergence threshold and the maximum number of iterations
      THRESH = MAX( ZERO, CNTL(1) )
      MAXIT  = ICNTL(7)
C Check for the particular inconsistency between MAXIT and THRESH.
Cnew  IF ( (MAXIT.LE.0) .AND. (CNTL(1).LE.ZERO) )  THEN
Cnew     INFO(1) = -14
Cnew     IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9010) INFO(1),MAXIT,CNTL(1)
Cnew     GO TO 99
Cnew  ENDIF

C   ICNTL(8) could be set by the user to indicate the frequency at which
C   convergence should be checked when CNTL(1) is strictly positive.
C   The default value used here 1 (e.g. check convergence every
C   iteration).
C     CHECK  = ICNTL(8)
      CHECK  = 1
      IF (CNTL(1).LE.ZERO)  CHECK = 0

C Prepare temporary data in working arrays as apropriate
      K = 2*M
      IF (ICNTL(6).EQ.0)  K = 2*(M+N)
      SETUP = 0
      IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
        SETUP = 1
C       Copy ABS(A(J))**JOB into DW(K+J), J=1..NNZ
        IF (JOB.EQ.-1) THEN
          DP = CNTL(2)
        ELSE
          DP = REAL(JOB)
        ENDIF
        IF (ICNTL(6).EQ.0) THEN
          DO 25 J = 1,N
            DO 20 I = 1,M
              DW(K+(J-1)*LDA+I) = ABS(A(I,J))**DP
   20       CONTINUE
   25     CONTINUE
        ELSE
          NNZ = (M*(M+1))/2
          DO 30 J = 1,NNZ
            DW(K+J) = ABS(A(J,1))**DP
   30     CONTINUE
        ENDIF
C       Reset the Threshold to take into account the Hadamard p-th power
        THRESH = ONE - (ABS(ONE-THRESH))**DP
      ELSE
        IF (ICNTL(5).EQ.0) THEN
          SETUP = 1
C         Copy ABS(A(J)) into DW(K+J), J=1..NNZ
          IF (ICNTL(6).EQ.0) THEN
            DO 40 J = 1,N
              DO 35 I = 1,M
                DW(K+(J-1)*LDA+I) = ABS(A(I,J))
   35         CONTINUE
   40       CONTINUE
          ELSE
            NNZ = (M*(M+1))/2
            DO 45 J = 1,NNZ
              DW(K+J) = ABS(A(J,1))
   45       CONTINUE
          ENDIF
        ENDIF
      ENDIF

C Begin the computations (input matrix in DENSE format) ...
C     We consider first the un-symmetric case :
      IF (ICNTL(6).EQ.0)  THEN
C
C       Equilibrate matrix A in the infinity norm
        IF (JOB.EQ.0) THEN
C         IW(1:M+N), DW(M+N:2*(M+N)) are workspaces
          IF (SETUP.EQ.1) THEN
            CALL MC77J(M,N,DW(2*(M+N)+1),LDA,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77J(M,N,A,LDA,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
C
C       Equilibrate matrix A in the p-th norm, 1 <= p < +infinity
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77K(M,N,DW(2*(M+N)+1),LDA,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ELSE
            CALL MC77K(M,N,A,LDA,DW(1),DW(M+1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),IW(M+1),DW(M+N+1),DW(2*M+N+1),INFO(1))
          ENDIF
C         Reset the scaling factors to their Hadamard p-th root, if
C         required and re-compute the actual value of the final error
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 50 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+N+I)**DP) )
   50       CONTINUE
            RINFO(2) = ZERO
            DO 55 J = 1,N
              DW(M+J) = DW(M+J)**DP
              IF (IW(M+J).NE.0)
     &           RINFO(2) = MAX( RINFO(2), ABS(ONE-DW(2*M+N+J)**DP) )
   55       CONTINUE
          ENDIF
        ENDIF
C
C     We treat then the symmetric (packed) case :
      ELSE
C
C       Equilibrate matrix A in the infinity norm
        IF (JOB.EQ.0) THEN
C         IW(1:M), DW(M:2*M) are workspaces
          IF (SETUP.EQ.1) THEN
            CALL MC77L(M,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77L(M,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
C
C       Equilibrate matrix A in the p-th norm, 1 <= p < +infinity
        ELSE
          IF (SETUP.EQ.1) THEN
            CALL MC77M(M,DW(2*M+1),DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ELSE
            CALL MC77M(M,A,DW(1),
     &                  CHECK,THRESH,RINFO(1),MAXIT,INFO(3),
     &                  IW(1),DW(M+1),INFO(1))
          ENDIF
C         Reset the scaling factors to their Hadamard p-th root, if
C         required and re-compute the actual value of the final error
          IF ( (JOB.GE.2) .OR. (JOB.EQ.-1) ) THEN
            IF (JOB.EQ.-1) THEN
              DP = ONE / CNTL(2)
            ELSE
              DP = ONE / REAL(JOB)
            ENDIF
            RINFO(1) = ZERO
            DO 60 I = 1,M
              DW(I) = DW(I)**DP
              IF (IW(I).NE.0)
     &           RINFO(1) = MAX( RINFO(1), ABS(ONE-DW(M+I)**DP) )
   60       CONTINUE
          ENDIF
        ENDIF
        RINFO(2) = RINFO(1)
C
      ENDIF
C     End of the un-symmetric/symmetric IF statement


C Print diagnostics on output
C ICNTL(9) could be set by the user to monitor the level of diagnostics
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9030) (INFO(J),J=1,3),(RINFO(J),J=1,2)
        IF (ICNTL(9).LT.2) THEN
          WRITE(ICNTL(3),9031) (DW(I),I=1,M)
          IF (ICNTL(6).EQ.0)
     &      WRITE(ICNTL(3),9032) (DW(M+J),J=1,N)
        ENDIF
      ENDIF

C Return from subroutine.
   99 CONTINUE
      RETURN

 9001 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3,
     &        ' because ',(A),' = ',I10)
 9002 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        Bad input control flag.',
     &        ' Value of ICNTL(',I1,') = ',I8)
 9003 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        Input matrix is symmetric and N /= M.',
     &        ' Value of (N-M) = ',I8)
 9004 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        LDW too small, must be at least ',I8)
 9008 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        Bad input REAL control parameter.',
     &        ' Value of CNTL(',I1,') = ',1PD14.4)
 9009 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
     &        '        Only scaling in infinity norm is allowed',
     &        ' for rectangular matrices'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8)
C9010 FORMAT (' ****** Error in MC77C/CD. INFO(1) = ',I3/
Cnew &        '        Algorithm will loop indefinitely !!!',
Cnew &        '     Check for inconsistency'/
Cnew &        ' Maximum number of iterations = ',I8/
Cnew &        ' Threshold for convergence    = ',1PD14.4)
 9020 FORMAT (' ****** Input parameters for MC77C/CD:'/
     &        ' JOB = ',I8/' M   = ',I8/' N   = ',I8/
     &        ' Max n.b. Iters.  = ',I8/
     &        ' Cvgce. Threshold = ',1PD14.4)
 9021 FORMAT (' ROW     (I) = ',I8/
     &        ' A(I,1:N)    = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9022 FORMAT (' ROW     (I) = ',I8/
     &        ' A(I,I:N)    = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for MC77C/CD:'/
     &        ' INFO(1:3)   = ',3(2X,I8)/
     &        ' RINFO(1:2)  = ',2(1PD14.4))
 9031 FORMAT (' DW(1:M)     = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9032 FORMAT (' DW(M+1:M+N) = ',4(1PD14.4)/(15X,4(1PD14.4)))
c9031 FORMAT (' DW(1:M)     = ',5(F11.3)/(15X,5(F11.3)))
c9032 FORMAT (' DW(M+1:M+N) = ',5(F11.3)/(15X,5(F11.3)))
      END

C**********************************************************************
      SUBROUTINE MC77J(M,N,A,LDA,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
C
C *********************************************************************
C ***           Infinity norm  ---  The un-symmetric case           ***
C ***                         DENSE format                          ***
C *********************************************************************
      INTEGER M,N,LDA,MAXIT,NITER,CHECK,INFO
      INTEGER IW(M),JW(N)
      REAL THRESH,ERR(2)
      REAL A(LDA,N),D(M),E(N),DW(M),EW(N)

C Variables are described in MC77B/BD
C The input matrix A is dense and un-symmetric.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO

C Initialisations ...
      DO 5 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
        D(I)  = ONE
    5 CONTINUE
      DO 10 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
        E(J)  = ONE
   10 CONTINUE

C Now, we compute the initial infinity-norm of each row and column.
      DO 20 J=1,N
        DO 30 I=1,M
          IF (EW(J).LT.A(I,J)) THEN
            EW(J) = A(I,J)
            JW(J) = I
          ENDIF
          IF (DW(I).LT.A(I,J)) THEN
            DW(I) = A(I,J)
            IW(I) = J
          ENDIF
   30   CONTINUE
   20 CONTINUE

      DO 40 I=1,M
        IF (IW(I).GT.0) D(I) = SQRT(DW(I))
   40 CONTINUE
      DO 45 J=1,N
        IF (JW(J).GT.0) E(J) = SQRT(EW(J))
   45 CONTINUE

      DO 50 J=1,N
        I = JW(J)
        IF (I.GT.0) THEN
          IF (IW(I).EQ.J) THEN
            IW(I) = -IW(I)
            JW(J) = -JW(J)
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 I=1,M
        IF (IW(I).GT.0)  GOTO 99
   60 CONTINUE
      DO 65 J=1,N
        IF (JW(J).GT.0)  GOTO 99
   65 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE

        DO 120 J=1,N
          IF (JW(J).GT.0) THEN
            DO 130 I=1,M
              S = A(I,J) / (D(I)*E(J))
              IF (EW(J).LT.S) THEN
                EW(J) = S
                JW(J) = I
              ENDIF
              IF (IW(I).GT.0) THEN
                IF (DW(I).LT.S) THEN
                  DW(I) = S
                  IW(I) = J
                ENDIF
              ENDIF
  130       CONTINUE
          ELSE
            DO 140 I=1,M
              IF (IW(I).GT.0) THEN
                S = A(I,J) / (D(I)*E(J))
                IF (DW(I).LT.S) THEN
                  DW(I) = S
                  IW(I) = J
                ENDIF
              ENDIF
  140       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 I=1,M
          IF (IW(I).GT.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).GT.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE

        DO 160 J=1,N
          I = JW(J)
          IF (I.GT.0) THEN
            IF (IW(I).EQ.J) THEN
              IW(I) = -IW(I)
              JW(J) = -JW(J)
            ENDIF
          ENDIF
  160   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in arrays
C       DW and EW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).GT.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).GT.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) )  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77K(M,N,A,LDA,D,E,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IW,JW,DW,EW,INFO)
C
C *********************************************************************
C ***              One norm  ---  The un-symmetric case             ***
C ***                         DENSE format                          ***
C *********************************************************************
      INTEGER M,N,LDA,MAXIT,NITER,CHECK,INFO
      INTEGER IW(M),JW(N)
      REAL THRESH,ERR(2)
      REAL A(LDA,N),D(M),E(N),DW(M),EW(N)

C Variables are described in MC77B/BD
C The input matrix A is dense and un-symmetric.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J, K
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR(1) = ZERO
      ERR(2) = ZERO

C Initialisations ...
      DO 10 I = 1, M
        IW(I) = 0
        DW(I) = ZERO
        D(I)  = ONE
   10 CONTINUE
      DO 15 J = 1, N
        JW(J) = 0
        EW(J) = ZERO
        E(J)  = ONE
   15 CONTINUE

C Now, we compute the initial one-norm of each row and column.
      DO 20 J=1,N
        DO 30 I=1,M
          IF (A(I,J).GT.ZERO) THEN
            DW(I) = DW(I) + A(I,J)
            IW(I) = IW(I) + 1
            EW(J) = EW(J) + A(I,J)
            JW(J) = JW(J) + 1
          ENDIF
   30   CONTINUE
   20 CONTINUE

      DO 40 I=1,M
        IF (IW(I).GT.0) D(I) = SQRT(DW(I))
   40 CONTINUE
      DO 45 J=1,N
        IF (JW(J).GT.0) E(J) = SQRT(EW(J))
   45 CONTINUE

      DO 60 K=1,M
        IF ( (IW(K).GT.0) .OR. (JW(K).GT.0) )  GOTO 99
   60 CONTINUE
C Since we enforce M=N in the case of the one-norm, the tests can be
C done in one shot. For future relases, (M/=N) we should incorporate
C instead :
Cnew  DO 60 I=1,M
Cnew    IF (IW(I).GT.0)  GOTO 99
Cne60 CONTINUE
Cnew  DO 65 J=1,N
Cnew    IF (JW(J).GT.0)  GOTO 99
Cne65 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 I=1,M
          DW(I) = ZERO
  110   CONTINUE
        DO 115 J=1,N
          EW(J) = ZERO
  115   CONTINUE

        DO 120 J=1,N
          IF (JW(J).GT.0) THEN
            DO 130 I=1,M
              S = A(I,J) / (D(I)*E(J))
              EW(J) = EW(J) + S
              DW(I) = DW(I) + S
  130       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 I=1,M
          IF (IW(I).GT.0) D(I) = D(I)*SQRT(DW(I))
  150   CONTINUE
        DO 155 J=1,N
          IF (JW(J).GT.0) E(J) = E(J)*SQRT(EW(J))
  155   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in arrays
C       DW and EW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR(1) = ZERO
        DO 170 I=1,M
          IF (IW(I).GT.0) ERR(1) = MAX( ERR(1), ABS(ONE-DW(I)) )
  170   CONTINUE
        ERR(2) = ZERO
        DO 175 J=1,N
          IF (JW(J).GT.0) ERR(2) = MAX( ERR(2), ABS(ONE-EW(J)) )
  175   CONTINUE
        IF ( (ERR(1).LT.THRESH) .AND. (ERR(2).LT.THRESH) ) GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77L(N,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
C
C *********************************************************************
C ***         Infinity norm  ---  The symmetric-packed case         ***
C ***                         DENSE format                          ***
C *********************************************************************
      INTEGER N,MAXIT,NITER,CHECK,INFO
      INTEGER IJW(N)
      REAL THRESH,ERR
      REAL A(*),DE(N),DEW(N)

C Variables are described in MC77B/BD
C The lower triangular part of the dense symmetric matrix A
C is stored contiguously column by column.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J, K
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR = ZERO

C Initialisations ...
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE

C Now, we compute the initial infinity-norm of each row and column.
      K = 1
      DO 20 J=1,N
        DO 30 I=J,N
          IF (DEW(J).LT.A(K)) THEN
            DEW(J) = A(K)
            IJW(J) = I
          ENDIF
          IF (DEW(I).LT.A(K)) THEN
            DEW(I) = A(K)
            IJW(I) = J
          ENDIF
          K = K+1
   30   CONTINUE
   20 CONTINUE

      DO 40 K=1,N
        IF (IJW(K).GT.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE

      DO 50 J=1,N
        I = IJW(J)
        IF (I.GT.0) THEN
          IF (IJW(I).EQ.J) THEN
            IJW(I) = -IJW(I)
            IF (I.NE.J)  IJW(J) = -IJW(J)
          ENDIF
        ENDIF
   50 CONTINUE

      DO 60 K=1,N
        IF (IJW(K).GT.0)  GOTO 99
   60 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE

        K = 1
        DO 120 J=1,N
          IF (IJW(J).GT.0) THEN
            DO 130 I=J,N
              S = A(K) / (DE(I)*DE(J))
              IF (DEW(J).LT.S) THEN
                DEW(J) = S
                IJW(J) = I
              ENDIF
              IF (IJW(I).GT.0) THEN
                IF (DEW(I).LT.S) THEN
                  DEW(I) = S
                  IJW(I) = J
                ENDIF
              ENDIF
              K = K+1
  130       CONTINUE
          ELSE
            DO 140 I=J,N
              IF (IJW(I).GT.0) THEN
                S = A(K) / (DE(I)*DE(J))
                IF (DEW(I).LT.S) THEN
                  DEW(I) = S
                  IJW(I) = J
                ENDIF
              ENDIF
              K = K+1
  140       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 K=1,N
          IF (IJW(K).GT.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE

        DO 160 J=1,N
          I = IJW(J)
          IF (I.GT.0) THEN
            IF (IJW(I).EQ.J) THEN
              IJW(I) = -IJW(I)
              IF (I.NE.J)  IJW(J) = -IJW(J)
            ENDIF
          ENDIF
  160   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in array
C       DEW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).GT.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END

C**********************************************************************
      SUBROUTINE MC77M(N,A,DE,
     &                  CHECK,THRESH,ERR,MAXIT,NITER,
     &                  IJW,DEW,INFO)
C
C *********************************************************************
C ***            One norm  ---  The symmetric-packed case           ***
C ***                         DENSE format                          ***
C *********************************************************************
      INTEGER N,MAXIT,NITER,CHECK,INFO
      INTEGER IJW(N)
      REAL THRESH,ERR
      REAL A(*),DE(N),DEW(N)

C Variables are described in MC77B/BD
C The lower triangular part of the dense symmetric matrix A
C is stored contiguously column by column.

C Local variables and parameters
      REAL ONE, ZERO
      PARAMETER ( ONE=1.0E+00, ZERO=0.0E+00 )
      INTEGER I, J, K
      REAL S
C Intrinsic functions
      INTRINSIC SQRT, MAX, ABS

      INFO = 0
      NITER = 0
      ERR = ZERO

C Initialisations ...
      DO 10 K = 1, N
        IJW(K) = 0
        DEW(K) = ZERO
        DE(K)  = ONE
   10 CONTINUE

C Now, we compute the initial one-norm of each row and column.
      K = 1
      DO 20 J=1,N
        DO 30 I=J,N
          IF (A(K).GT.ZERO) THEN
            DEW(J) = DEW(J) + A(K)
            IJW(J) = IJW(J) + 1
            IF (I.NE.J) THEN
              DEW(I) = DEW(I) + A(K)
              IJW(I) = IJW(I) + 1
            ENDIF
          ENDIF
          K = K+1
   30   CONTINUE
   20 CONTINUE

      DO 40 K=1,N
        IF (IJW(K).GT.0) DE(K) = SQRT(DEW(K))
   40 CONTINUE

      DO 60 K=1,N
        IF (IJW(K).GT.0)  GOTO 99
   60 CONTINUE
      GOTO 200


C Then, iterate on the normalisation of rows and columns.
   99 NITER = NITER + 1
      IF ( (NITER.GT.MAXIT) .AND. (MAXIT.GE.0) )  THEN
         IF (CHECK.GT.0)  INFO = 1
         NITER = NITER - 1
         GOTO 100
      ENDIF

        DO 110 K=1,N
          DEW(K) = ZERO
  110   CONTINUE

        K = 1
        DO 120 J=1,N
          IF (IJW(J).GT.0) THEN
            DO 130 I=J,N
              S = A(K) / (DE(I)*DE(J))
              DEW(J) = DEW(J) + S
              IF (I.NE.J)  DEW(I) = DEW(I) + S
              K = K+1
  130       CONTINUE
          ENDIF
  120   CONTINUE

        DO 150 K=1,N
          IF (IJW(K).GT.0) DE(K) = DE(K)*SQRT(DEW(K))
  150   CONTINUE

C Stopping criterion :
        IF (CHECK.LE.0) GOTO 99
C       IF (MOD(NITER,CHECK).NE.0) GOTO 99

C N.B.  the test is performed on the basis of the values in array
C       DEW which, in fact, correspond to the row and column
C       norms of the scaled matrix at the previous iteration ...
  100   IF (INFO.NE.0)  GOTO 200
        ERR = ZERO
        DO 170 K=1,N
          IF (IJW(K).GT.0) ERR = MAX( ERR, ABS(ONE-DEW(K)) )
  170   CONTINUE
        IF (ERR.LT.THRESH)  GOTO 200

      IF (CHECK.GT.0)  GOTO 99

  200 RETURN
      END
