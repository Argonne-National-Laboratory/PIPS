C *******************************************************************
C COPYRIGHT (c) 1999 CCLRC Council for the Central Laboratory
*                    of the Research Councils
C All rights reserved.
C
C None of the comments in this Copyright notice between the lines
C of asterisks shall be removed or altered in any way.
C
C This Package is intended for compilation without modification,
C so most of the embedded comments have been removed.
C
C ALL USE IS SUBJECT TO LICENCE. For full details of a UK ACADEMIC
C Licence, see http://hsl.rl.ac.uk/acuk/cou.html
C
C Please note that for a UK ACADEMIC Licence:
C
C 1. The Packages may only be used for academic research or teaching
C    purposes by the Licensee, and must not be copied by the Licensee for
C    use by any other persons. Use of the Packages in any commercial
C    application shall be subject to prior written agreement between
C    Hyprotech UK Limited and the Licensee on suitable terms and
C    conditions, which will include financial conditions.
C 2. All information on the Package is provided to the Licensee on the
C    understanding that the details thereof are confidential.
C 3. All publications issued by the Licensee that include results obtained
C    with the help of one or more of the Packages shall acknowledge the
C    use of the Packages. The Licensee will notify the Numerical Analysis
C    Group at Rutherford Appleton Laboratory of any such publication.
C 4. The Packages may be modified by or on behalf of the Licensee
C    for such use in research applications but at no time shall such
C    Packages or modifications thereof become the property of the
C    Licensee. The Licensee shall make available free of charge to the
C    copyright holder for any purpose all information relating to
C    any modification.
C 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
C    direct or consequential loss or damage whatsoever arising out of
C    the use of Packages by the Licensee.
C *******************************************************************
C
C Original date 13 September 1999
C 01/11/00  Entries in IW initialized to zero in MA57OD to avoid copy
C           of unassigned variables by MA57ED.
C           AINPUT and IINPUT reset in call to MA57ED.
C 06/02/01  Default values for ICNTL(12) and ICNTL(13) changed.
C           Control for direct addressing in solve changed to be
C           on number of rows and number columns in block pivot.
C           Several comments changed as consequence.
C           INFO(31) added to record number of block pivots.
C           Subtroutines MA57XD and MA57YD added for efficiency when only
C           one rhs (equivalent to MA57QD and MA57RD resp).
C 04/07/01  Use of MC41 changed to use of MC71.
C 26/10/01  Printing controls corrected to ensure ICNTL(5) is used and
C           unit number always checked for being positive before printing.
C           Text and comments changed to reflect that D inverse is held in
C           factors and text for solution changed from Right-hand side
C           to solution.
C           Option of choosing two 1 x 1 pivots when 2 x 2 fails
C           removed.
C           MC47B/BD given remaining length in KEEP to avoid compresses
C 20/12/01  INFO(1) initilaized to zero in MA57ED
C 06/12/02  The test for convergence of iterative refinement changed to
C           avoid any problem with comparisons of numbers held in
C           registers.
C 25/03/03  MC50 (AMD with dense row protection) and MA27 (minimum degree)
C           added.  Invoked by ICNTL(6) equal to 2 and 3,respectively.
C           Routines MA57H/HD, MA57V/VD, and MA57Z/ZD have been added to
C           duplicate routines MA27H/HD, MA27G/GD, and MA27U/UD from
C           MA57 and MC50B/BD is another internal routine of MA57.
C           ICNTL(14) has been added to control density of rows regarded
C           as dense by the MC50 and MA27 orderings.
C 24/05/04  Statment functions in MA57U/UD replaced by in-line code.

C 12th July 2004 Version 1.0.0. Version numbering added.

C 20/07/04  Several changes incorporated for HSL 2004 code.
C           Removed unused INT,ABS from MA57U/UD
C           INFO(32), INFO(33), and INFO(34) added
C           INFO(32): number of zeros in the triangle of the factors
C           INFO(33): number of zeros in the rectangle of the factors
C           INFO(34): number of zero columns in rectangle of the factors
C           Static pivoting available (controlled by CNTL(4), CNTL(5))
C           Scaling using symmetrized MC64 (ICNTL(15))
C           Links to METIS_NODEND ordering


C 31st July 2004 Version 2.0.0 established at HSL 2004 release.
C 1st Sept  2004 Version 2.1.0. Default changed to static pivoting off.
C 10th Sept 2004 Version 2.2.0. Defaults for ICNTL(6), ICNTL(9) and
C           CNTL(5) changed. Scaling factors (optionally) printed.
C  4th Nov  2004 Version 2.2.1. Change to assembly of reals in MA57OD
C           leading to more efficient code at suggestion of Stephane
C           Pralet.

      SUBROUTINE MA57ID(CNTL, ICNTL)
C****************************************************************
      DOUBLE PRECISION    CNTL(5)
      INTEGER             ICNTL(20)
      INTEGER I
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C===============================================
C===============================================
      CNTL(1)   = 0.01D0
      CNTL(2)   = 1.0D-20
      CNTL(3)   = 0.5D0
      CNTL(4) = ZERO
      CNTL(5) = ZERO
      ICNTL(1)  = 6
      ICNTL(2)  = 6
      ICNTL(3)  = 6
      ICNTL(4)  = -1
      ICNTL(5)  = 2
      ICNTL(6)  = 2
      ICNTL(7)  = 1
      ICNTL(8)  = 0
      ICNTL(9)  = 10
      ICNTL(10) = 0
      ICNTL(11) = 16
      ICNTL(12) = 16
      ICNTL(13) = 10
      ICNTL(14) = 100
      ICNTL(15) = 1
      DO 110 I=16,20
        ICNTL(I) = 0
  110 CONTINUE
      RETURN
      END
CCC
      SUBROUTINE MA57AD(N,NE,IRN,JCN,LKEEP,KEEP,IWORK,ICNTL,INFO,RINFO)
      INTEGER N,NE,IRN(NE),JCN(NE),IWORK(5*N),LKEEP,KEEP(LKEEP),
     *        ICNTL(20),INFO(40)
      DOUBLE PRECISION RINFO(20)
C**** Still to be updated
      INTRINSIC MIN
      EXTERNAL MA57GD,MC47BD,MC50BD,MA57VD,MA57HD,MA57JD,MA57KD,
     *         MA57LD,MA57MD,MA57ND
      INTEGER I,IL,IN,IPE,IRNPRM,COUNT,FILS,FRERE,HOLD,IFCT,INVP,IPS,
     +        IW,IWFR,K,LDIAG,LP,LW,LROW,MAP,EXPNE,
     +        MP,NCMPA,NEMIN,NODE,NST,NSTEPS,NV,PERM,
     +        IW1,IW2,IW3,IW4,IW5,NSTK,ND,NELIM,NZE,ALENB
      INTEGER METOPT(8),METFTN
      DOUBLE PRECISION ZERO,THRESH
      PARAMETER (ZERO=0.0D0)
      LP = ICNTL(1)
      MP = ICNTL(3)
      LDIAG = ICNTL(5)
      DO 10 I = 1,40
        INFO(I) = 0
   10 CONTINUE
      DO 11 I = 1,20
        RINFO(I) = ZERO
   11 CONTINUE
      IF (N.LT.1)  GO TO 20
      IF (NE.LT.0) GO TO 30
      IF (LKEEP.LT.5*N+NE+MAX(N,NE)+42) GO TO 40
      IF (ICNTL(6).EQ.1) THEN
        DO 12 I = 1,N
          IWORK(I) = 0
   12   CONTINUE
        DO 14 I=1,N
          K = KEEP(I)
          IF (K.LE.0 .OR. K.GT.N) GO TO 80
          IF (IWORK(K).NE.0) GO TO 80
          IWORK(K) = I
   14   CONTINUE
      ENDIF
      IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
        WRITE(MP,99980) N,NE,(ICNTL(I),I=1,7),ICNTL(12),ICNTL(15)
99980 FORMAT (//'Entering analysis phase (MA57AD) with ...'/
     1      'N         Order of matrix                     =',I12/
     2      'NE        Number of entries                   =',I12/
     6      'ICNTL(1)  Stream for errors                   =',I12/
     7      ' --- (2)  Stream for warnings                 =',I12/
     8      ' --- (3)  Stream for monitoring               =',I12/
     9      ' --- (4)  Stream for statistics               =',I12/
     1      ' --- (5)  Level of diagnostic printing        =',I12/
     2      ' --- (6)  Flag for input pivot order          =',I12/
     2      ' --- (7)  Numerical pivoting control (st est) =',I12/
     2      ' --- (12) Node amalgamation parameter         =',I12/
     2      ' --- (15) Scaling control (storage estimate)  =',I12)
        K = MIN(10,NE)
        IF (LDIAG.GE.4) K = NE
        WRITE (MP,'(/A/(3(I6,A,2I8,A)))') ' Matrix entries:',
     +        (I,': (',IRN(I),JCN(I),')',I=1,K)
        IF (K.LT.NE) WRITE (MP,'(A)') '     . . .'
        IF (ICNTL(6).EQ.1) THEN
          K = MIN(10,N)
          IF (LDIAG.GE.4) K = N
          WRITE (MP,'(A,10I6:/(7X,10I6))') ' KEEP =', (KEEP(I),I=1,K)
          IF (K.LT.N) WRITE (MP,'(7X,A)') '     . . .'
        END IF
      END IF
      IW1 = 1
      IW2 = IW1 + N
      IW3 = IW2 + N
      IW4 = IW3 + N
      IW5 = IW4 + N
      FILS  = IW1
      FRERE = IW2
      ND    = IW3
      NELIM = IW4
      NV    = IW5
      PERM = 1
      NSTEPS = PERM + N
      EXPNE  = NSTEPS + 1
      HOLD   = EXPNE + 1
      LROW = HOLD + 40
      NODE = LROW + N
      NSTK = NODE + N
      MAP  = NSTK + N
      IRNPRM = MAP + MAX(N,NE)
      INVP  = NODE
      IW    = NODE
      IPE   = LROW
      IFCT  = MAP
      IPS   = MAP
      COUNT = NSTK
      KEEP(HOLD) = 0
      IF (ICNTL(6).NE.1) THEN
      IF (ICNTL(6) .NE. 3) THEN
        CALL MA57GD(N,NE,IRN,JCN,KEEP(IFCT),KEEP(IPE),KEEP(COUNT),
     +              KEEP(IW),IWFR,ICNTL,INFO)
        IF (ICNTL(6).EQ.4) THEN
          KEEP(IPE+N) = IWFR
          METFTN    = 1
          METOPT(1) = 0
CCCCCC
          CALL METIS_NODEND(N,KEEP(IPE),KEEP(IFCT),METFTN,METOPT,
     *                      KEEP(NSTK),KEEP(PERM))
          IF (KEEP(PERM).EQ.-1) GO TO 90
          GO TO 111
        ENDIF
      LW = LKEEP-IFCT+1
        IF (ICNTL(6).EQ.2) THEN
          CALL MC50BD(ICNTL(14),N,LW,KEEP(IPE),IWFR,KEEP(COUNT),
     +                KEEP(IFCT),IWORK(NV),
     +                KEEP(INVP),KEEP(PERM),NCMPA,IWORK(IW1),
     +                IWORK(IW2),IWORK(IW3),IWORK(IW4))
        ELSE
          CALL MC47BD(N,LW,KEEP(IPE),IWFR,KEEP(COUNT),
     +                KEEP(IFCT),IWORK(NV),
     +                KEEP(INVP),KEEP(PERM),NCMPA,IWORK(IW1),
     +                IWORK(IW2),IWORK(IW3),IWORK(IW4))
        ENDIF
        INFO(13) = NCMPA
        ELSE
      LW = LKEEP-IFCT+1
         CALL MA57VD(N,NE,IRN,JCN,KEEP(IFCT),LW,KEEP(IPE),IWORK(IW1),
     *               IWORK(IW2),IWFR,ICNTL,INFO)
         THRESH = FLOAT(ICNTL(14))/100.0
         CALL MA57HD(N,KEEP(IPE),KEEP(IFCT),LW,IWFR,IWORK(NV),
     *               IWORK(IW1),IWORK(IW2),IWORK(IW3),IWORK(IW4),
     +               2139062143,INFO(13),THRESH)
      DO 110 I = 1,N
        IF (IWORK(NV+I-1).NE.0) GO TO 110
        IN = I
  105   IL = IN
        IN = - KEEP(IPE+IL-1)
        IF (IWORK(NV+IN-1).EQ.0) GO TO 105
        KEEP(IPE+I-1) = -IN
  110 CONTINUE
        ENDIF
      ENDIF
  111 IF (ICNTL(6).EQ.1 .OR. ICNTL(6).EQ.4) THEN
        CALL MA57JD(N,NE,IRN,JCN,KEEP(PERM),KEEP(IFCT),KEEP(IPE),
     +              KEEP(COUNT),IWORK(IW1),IWFR,ICNTL,INFO)
        LW = 2*NE
        CALL MA57KD(N,KEEP(IPE),KEEP(IFCT),LW,IWFR,KEEP(PERM),
     +              KEEP(INVP),IWORK(NV),IWORK(IW1),NCMPA)
        INFO(13) = NCMPA
      END IF
      NEMIN = ICNTL(12)
      CALL MA57LD(N,KEEP(IPE),IWORK(NV),KEEP(IPS),IWORK(NELIM),
     +            KEEP(NSTK),KEEP(NODE),KEEP(PERM),
     +            KEEP(NSTEPS),IWORK(FILS),IWORK(FRERE),IWORK(ND),
     +            NEMIN,KEEP(IRNPRM))
      NST = KEEP(NSTEPS)
      CALL MA57MD(N,NE,IRN,JCN,KEEP(MAP),KEEP(IRNPRM),
     +            KEEP(LROW),KEEP(PERM),
     +            IWORK(IW2),IWORK(IW5))
      KEEP(EXPNE) = IWORK(IW5)
      CALL MA57ND(N,KEEP(LROW),KEEP(NSTK),IWORK(NELIM),
     +            IWORK(ND),NST,IWORK(IW1),IWORK(IW2),
     +            INFO,RINFO)
      ALENB    = 1
      IF (ICNTL(7).EQ.4) ALENB = ALENB + N + 5
      IF (ICNTL(15).EQ.1) ALENB = ALENB + N
      INFO(9)  = MAX(INFO(9)+ALENB,ALENB+KEEP(EXPNE)+1)
      INFO(11) = MAX(INFO(11)+ALENB,ALENB+KEEP(EXPNE)+1)
      INFO(10) = MAX(INFO(10),KEEP(EXPNE)+N+5)
      INFO(12) = MAX(INFO(12),KEEP(EXPNE)+N+5)
      IF (ICNTL(15).EQ.1) THEN
        INFO(9) = MAX(INFO(9),ALENB+3*KEEP(EXPNE)+3*N)
        INFO(11) = MAX(INFO(11),ALENB+3*KEEP(EXPNE)+3*N)
        INFO(10) = MAX(INFO(10),3*KEEP(EXPNE)+5*N+1)
        INFO(12) = MAX(INFO(12),3*KEEP(EXPNE)+5*N+1)
      ENDIF
      IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
        NZE = KEEP(EXPNE)
        WRITE (MP,99999) INFO(1),NZE,
     *                  (INFO(I),I=3,13),(RINFO(I),I=1,2)
99999 FORMAT (/'Leaving analysis phase (MA57AD) with ...'/
     1    'INFO(1)  Error indicator                      =',I12/
     2    'Number of entries in matrix with diagonal     =',I12/
     2    'INFO(3)  Number of out-of-range indices       =',I12/
     2    'INFO(4)  Number of off-diagonal duplicates    =',I12/
     2    'INFO(5)  Forecast real storage for factors    =',I12/
     3    '----(6)  Forecast integer storage for factors =',I12/
     3    '----(7)  Forecast maximum front size          =',I12/
     4    '----(8)  Number of nodes in assembly tree     =',I12/
     5    '----(9)  Size of FACT without compress        =',I12/
     6    '----(10) Size of IFACT without compress       =',I12/
     5    '----(11) Size of FACT with compress           =',I12/
     5    '----(12) Size of IFACT with compress          =',I12/
     5    '----(13) Number of compresses                 =',I12/
     9    'RINFO(1) Forecast additions for assembly      =',1P,D12.5/
     9    'RINFO(2) Forecast ops for elimination         =',1P,D12.5)
        K = MIN(10,N)
        IF (LDIAG.GE.4) K = N
        WRITE (MP,'(/A/(5I12))')  'Permutation array:',
     +                  (KEEP(I),I=1,K)
        IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
        WRITE (MP,'(/A/(5I12))')
     +        'Number of entries in rows of permuted matrix:',
     +        (KEEP(LROW+I-1),I=1,K)
        IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
        K = MIN(10,NZE)
        IF (LDIAG.GE.4) K = NZE
        WRITE (MP,'(/A/(5I12))')
     *        'Column indices of permuted matrix:',
     *                           (KEEP(IRNPRM+I-1),I=1,K)
        IF (K.LT.NZE) WRITE (MP,'(16X,A)') '     . . .'
        K = MIN(10,N)
        IF (LDIAG.GE.4) K = N
        WRITE (MP,'(/A/(5I12))')
     +    'Tree nodes at which variables eliminated:',
     +    (KEEP(NODE+I-1),I=1,K)
        IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
        K = MIN(10,NE)
        IF (LDIAG.GE.4) K = NE
        WRITE (MP,'(/A/(5I12))') 'Map array:',
     *                               (KEEP(I),I=MAP,MAP+K-1)
        IF (K.LT.NE) WRITE (MP,'(16X,A)') ' . . .'
      END IF
      RETURN
   20 INFO(1) = -1
      INFO(2) = N
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
     +    '**** Error return from MA57AD ****  INFO(1) =',INFO(1),
     +    'N has value ',INFO(2)
      RETURN
   30 INFO(1) = -2
      INFO(2) = NE
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
     +    '**** Error return from MA57AD ****  INFO(1) =',INFO(1),
     +    'NE has value',INFO(2)
       RETURN
   40 INFO(1) = -15
      INFO(2) = LKEEP
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10/A,I10)')
     +    '**** Error return from MA57AD ****  INFO(1) =',INFO(1),
     +    'LKEEP has value    ',INFO(2),
     +    'Should be at least ',5*N+NE+MAX(N,NE)+42
       RETURN
   80 INFO(1) = -9
      INFO(2) = I
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A/A,I10,A)')
     +    '**** Error return from MA57AD ****  INFO(1) =',INFO(1),
     +    'Invalid permutation supplied in KEEP',
     +    'Component',INFO(2),' is faulty'
      RETURN
   90 INFO(1) = -18
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A)')
     +    '**** Error return from MA57AD ****  INFO(1) =',INFO(1),
     +    'MeTiS ordering requested but MeTiS not linked'
      END
CCC
C--------------------------------------------------------------------
C-             Copyright Rutherford Appleton Laboratory
C--------------------------------------------------------------------
      SUBROUTINE MA57BD(N, NE, A, FACT, LFACT, IFACT, LIFACT,
     * LKEEP, KEEP, PPOS, ICNTL, CNTL, INFO, RINFO)
      INTEGER N,NE,LFACT,LIFACT,LKEEP
      DOUBLE PRECISION A(NE),FACT(LFACT)
      DOUBLE PRECISION RINFO(20)
      DOUBLE PRECISION CNTL(5)
      INTEGER ICNTL(20), IFACT(LIFACT)
      INTEGER   INFO(40), KEEP(LKEEP), PPOS(N)
      INTEGER EXPNE,HOLD,I,IRNPRM,K,LDIAG,LLFACT,LP,LROW,MAP,MM1,MM2,MP
      INTEGER J,JJ,KK,ISCALE,NUM,NE64,IDUP,IMAT,IPT,JLOOP,JNEW,NN,ISING
      INTEGER NSTEPS,NODE,NSTK,PERM,INEW,ALENB,BIGA
      DOUBLE PRECISION ONE,ZERO,RINF,FD05AD,FCT,SMAX,SMIN,REPS
      PARAMETER (ONE = 1.0D0, ZERO=0.0D0)
C?? To identify bug
      INTRINSIC MIN
      EXTERNAL MA57OD,MA57UD,FD05AD,MC34AD,MC64WD
      RINF = FD05AD(5)
      REPS = FD05AD(1)
      LP     = ICNTL(1)
      MP     = ICNTL(3)
      LDIAG  = ICNTL(5)
C??
      IF (N.LE.0)  GO TO 25
      IF (NE.LT.0) GO TO 30
      IF (LKEEP.LT.5*N+NE+MAX(N,NE)+42) GO TO 40
      IF (ICNTL(7).LT.1 .OR. ICNTL(7).GT.4) GO TO 35
      NSTEPS = KEEP(N+1)
      EXPNE  = KEEP(N+2)
      PERM = 1
      HOLD = PERM + N + 2
      LROW = HOLD + 40
      NODE = LROW + N
      NSTK = NODE + N
      MAP  = NSTK + N
      IRNPRM = MAP + MAX(NE,N)
      BIGA = LFACT
      LLFACT = LFACT - 1
      IF (ICNTL(15).EQ.1) THEN
        ISCALE = LLFACT - N + 1
        LLFACT = ISCALE - 1
      ENDIF
      IF (ICNTL(7).EQ.4) THEN
        LLFACT = LLFACT - N - 5
        MM1 = LLFACT+6
        MM2 = LLFACT+1
      ELSE
        MM1 = 1
        MM2 = 1
      ENDIF
      IF (LIFACT.LT.EXPNE+N+5)  GO TO 95
      IF (LLFACT.LT.EXPNE+1)   GO TO 85
      ALENB = 1
      IF (ICNTL(7).EQ.4)  ALENB = ALENB + N + 5
      IF (ICNTL(15).EQ.1) ALENB = ALENB + N
      IF (ICNTL(15).EQ.1)  THEN
        IF (LFACT .LT. ALENB + 3*EXPNE  + 3*N) GO TO 85
        IF (LIFACT .LT. 3*EXPNE + 5*N + 1) GO TO 95
      ENDIF
C*****************************
      IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
        WRITE (MP,99999)
99999 FORMAT (//'Entering factorization phase (MA57BD) with ...')
        IF (KEEP(HOLD).GT.0) WRITE (MP,99998)
99998 FORMAT ('Re-entry call after call to MA57ED')
        WRITE (MP,99997) N,NE,EXPNE,(ICNTL(I),I=1,5),ICNTL(7),ICNTL(8),
     +         ICNTL(11),ICNTL(15),LFACT, LIFACT, NSTEPS,
     +         CNTL(1), CNTL(2), CNTL(4), CNTL(5)
99997 FORMAT ('N       Order of input matrix               =',I12/
     2        'NE      Entries in input matrix             =',I12/
     2        '        Entries in input matrix (inc diags) =',I12/
     6        'ICNTL(1)  Stream for errors                 =',I12/
     7        ' --- (2)  Stream for warnings               =',I12/
     8        ' --- (3)  Stream for monitoring             =',I12/
     9        ' --- (4)  Stream for statistics             =',I12/
     1        ' --- (5)  Level of diagnostic printing      =',I12/
     1        ' --- (7)  Numerical pivoting control        =',I12/
     1        ' --- (8)  Restart or discard factors        =',I12/
     1        ' --- (11) Block size for Level 3 BLAS       =',I12/
     1        ' --- (15) Scaling control (1 on)            =',I12/
     4        'LFACT   Size of real working space          =',I12/
     5        'LIFACT  Size of integer working space       =',I12/
     7        '        Number nodes in assembly tree       =',I12/
     9        'CNTL(1) Value of threshold parameter        =',D12.5/
     9        'CNTL(2) Threshold for zero pivot            =',D12.5/
     9        'CNTL(4) Control for value of static pivots  =',D12.5/
     9        'CNTL(5) Control for number delayed pivots   =',D12.5)
        K = MIN(10,NE)
        IF (LDIAG.GE.4) K = NE
        IF (NE.GT.0) THEN
          WRITE (MP,'(/A/(3(I6,A,1P,D16.8,A)))') 'Matrix entries:',
     +     (I,': (',A(I),')',I=1,K)
          IF (K.LT.NE) WRITE (MP,'(A)') '     . . .'
        END IF
        K = MIN(10,N)
        IF (LDIAG.GE.4) K = N
        WRITE (MP,'(/A/(5I12))')  'Permutation array:',
     +                    (KEEP(I),I=1,K)
        IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
        WRITE (MP,'(/A/(5I12))')
     +          'Number of entries in rows of permuted matrix:',
     +          (KEEP(LROW+I-1),I=1,K)
        IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
        WRITE (MP,'(/A/(5I12))')
     +    'Tree nodes at which variables eliminated:',
     +    (KEEP(NODE+I-1),I=1,K)
        IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
        K = MIN(10,NSTEPS)
        IF (LDIAG.GE.4) K = NSTEPS
        IF (K.GT.0) WRITE (MP,'(/A/(5I12))')
     +     'Number of assemblies at each tree node:',
     +     (KEEP(NSTK+I-1),I=1,K)
        IF (K.LT.NSTEPS) WRITE (MP,'(16X,A)') ' . . .'
        K = MIN(10,NE)
        IF (LDIAG.GE.4) K = NE
        WRITE (MP,'(/A/(5I12))') 'Map array:',
     *                               (KEEP(I),I=MAP,MAP+K-1)
        IF (K.LT.NE) WRITE (MP,'(16X,A)') ' . . .'
        K = MIN(10,EXPNE)
        IF (LDIAG.GE.4) K = EXPNE
        WRITE (MP,'(/A/(5I12))')
     *          'Column indices of permuted matrix:',
     *                             (KEEP(IRNPRM+I-1),I=1,K)
        IF (K.LT.EXPNE) WRITE (MP,'(16X,A)') '     . . .'
      ENDIF
      IF (KEEP(HOLD) .GT. 0) GO TO 22
C***************************************************
C***************************************************
C?? For the moment to handle missing diagonals
      DO 19 K = 1,EXPNE
        FACT(LLFACT-EXPNE+K) = ZERO
   19 CONTINUE
      FACT(BIGA) = ZERO
      DO 20 K = 1,NE
        FACT(BIGA) = MAX(FACT(BIGA),ABS(A(K)))
        FACT(KEEP(MAP+K-1)+LLFACT-EXPNE) = A(K)
   20 CONTINUE
      RINFO(18) = FACT(BIGA)
      DO 21 K = 1,EXPNE
        IFACT(LIFACT-EXPNE+K) = KEEP(IRNPRM+K-1)
   21 CONTINUE
      DO 23 I = 1,N
        PPOS(KEEP(PERM+I-1)) = I
   23 CONTINUE
      IF (ICNTL(15).EQ.1) THEN
        IPT = 1
        IDUP = IPT+N+1
        IMAT = IDUP+N
        ISING = IMAT + MAX(NE,EXPNE)
        DO 4444 I = 1,N
          IFACT(IDUP+I-1) = 0
 4444   CONTINUE
C9999   CONTINUE
        IFACT(IPT) = 1
        KK = 1
        K = 1
        DO 3333 J = 1,N
          DO 2222 JJ = 1,KEEP(LROW+J-1)
            I = KEEP(PERM+IFACT(LIFACT-EXPNE+K)-1)
            IF (IFACT(IDUP+I-1).GE.IFACT(IPT+J-1)) THEN
              FACT(IFACT(IDUP+I-1)) =
     &          FACT(IFACT(IDUP+I-1)) + FACT(LLFACT-EXPNE+K)
            ELSE
              IF (FACT(LLFACT-EXPNE+K).NE.ZERO) THEN
                IFACT(IDUP+I-1) = KK
                FACT(KK) = FACT(LLFACT-EXPNE+K)
                IFACT(IMAT-1+KK) = I
                KK = KK+1
              ENDIF
            ENDIF
            K = K + 1
 2222     CONTINUE
          IFACT(IPT+J) = KK
 3333   CONTINUE
        CALL MC34AD(N,IFACT(IMAT),IFACT(IPT),.TRUE.,FACT,KEEP(PERM))
        NE64 = IFACT(IPT+N)-1
        DO 75 J = 1,N
          FCT = ZERO
          DO 60 K = IFACT(IPT+J-1),IFACT(IPT+J)-1
            FACT(K) = ABS(FACT(K))
            IF (FACT(K).GT.FCT) FCT = FACT(K)
   60     CONTINUE
          FACT(NE64+2*N+J) = FCT
          IF (FCT.NE.ZERO) THEN
            FCT = LOG(FCT)
          ELSE
            FCT = RINF/N
          ENDIF
          DO 70 K = IFACT(IPT+J-1),IFACT(IPT+J)-1
            IF (FACT(K).NE.ZERO) THEN
              FACT(K) = FCT - LOG(FACT(K))
            ELSE
              FACT(K) = RINF/N
            ENDIF
   70     CONTINUE
   75   CONTINUE
        CALL MC64WD(N,NE64,IFACT(IPT),IFACT(IMAT),FACT,KEEP(PERM),NUM,
     &      IFACT(IDUP),IFACT(IMAT+NE64),IFACT(IMAT+NE64+N),
     &      IFACT(IMAT+NE64+2*N),IFACT(IMAT+NE64+3*N),
     &      FACT(NE64+1),FACT(NE64+N+1))
        IF (NUM.EQ.N) THEN
          DO 80 J = 1,N
            IF (FACT(NE64+2*N+J).NE.ZERO) THEN
              FACT(NE64+N+J) = FACT(NE64+N+J) - LOG(FACT(NE64+2*N+J))
            ELSE
              FACT(NE64+N+J) = ZERO
            ENDIF
   80     CONTINUE
          DO 5555 I=1,N
            FACT(ISCALE+PPOS(I)-1) =
     &        SQRT(EXP(FACT(NE64+I)+FACT(NE64+N+I)))
 5555     CONTINUE
        ELSE
        K = 0
        DO 3501 I = 1,N
          IF (KEEP(PERM+I-1).LT.0) THEN
            PPOS(I) = -PPOS(I)
            IFACT(ISING+I-1) = 0
          ELSE
            K = K + 1
            IFACT(ISING+I-1) = K
          ENDIF
 3501   CONTINUE
        DO 3502 I = 1,N
          KEEP(PERM+ABS(PPOS(I))-1) = I
 3502   CONTINUE
        DO 3503 I = 1,N
          IFACT(IDUP+I-1) = 0
 3503   CONTINUE
        IFACT(IPT) = 1
        KK = 1
        K = 1
        JNEW = 0
        NN = N
        DO 3505 J = 1,N
          IF (PPOS(J).LT.0) THEN
            NN = NN - 1
            K = K + KEEP(LROW+J-1)
            GO TO 3505
          ENDIF
          JNEW = JNEW + 1
          DO 3504 JJ = 1,KEEP(LROW+J-1)
            I = KEEP(PERM+IFACT(LIFACT-EXPNE+K)-1)
            IF (PPOS(I).GT.0) THEN
              IF (IFACT(IDUP+I-1).GE.IFACT(IPT+J-1)) THEN
                FACT(IFACT(IDUP+I-1)) =
     &            FACT(IFACT(IDUP+I-1)) + FACT(LLFACT-EXPNE+K)
              ELSE
                IF (FACT(LLFACT-EXPNE+K).NE.ZERO) THEN
                  IFACT(IDUP+I-1) = KK
                  FACT(KK) = FACT(LLFACT-EXPNE+K)
                  IFACT(IMAT-1+KK) = IFACT(ISING+I-1)
                  KK = KK+1
                ENDIF
              ENDIF
            ENDIF
            K = K + 1
 3504     CONTINUE
          IFACT(IPT+JNEW) = KK
 3505   CONTINUE
      NE64 = IFACT(IPT+NN)-1
        CALL MC34AD(NN,IFACT(IMAT),IFACT(IPT),.TRUE.,FACT,KEEP(PERM))
        NE64 = IFACT(IPT+NN)-1
        DO 3508 J = 1,NN
          FCT = ZERO
          DO 3506 K = IFACT(IPT+J-1),IFACT(IPT+J)-1
            FACT(K) = ABS(FACT(K))
            IF (FACT(K).GT.FCT) FCT = FACT(K)
 3506     CONTINUE
          FACT(NE64+2*N+J) = FCT
          IF (FCT.NE.ZERO) THEN
            FCT = LOG(FCT)
          ELSE
            FCT = RINF/NN
          ENDIF
          DO 3507 K = IFACT(IPT+J-1),IFACT(IPT+J)-1
            IF (FACT(K).NE.ZERO) THEN
              FACT(K) = FCT - LOG(FACT(K))
            ELSE
              FACT(K) = RINF/NN
            ENDIF
 3507     CONTINUE
 3508   CONTINUE
        CALL MC64WD(NN,NE64,IFACT(IPT),IFACT(IMAT),FACT,KEEP(PERM),NUM,
     &      IFACT(IDUP),IFACT(IMAT+NE64),IFACT(IMAT+NE64+N),
     &      IFACT(IMAT+NE64+2*N),IFACT(IMAT+NE64+3*N),
     &      FACT(NE64+1),FACT(NE64+N+1))
        DO 3509 J = 1,NN
            IF (FACT(NE64+2*N+J).NE.ZERO) THEN
              FACT(NE64+N+J) = FACT(NE64+N+J) - LOG(FACT(NE64+2*N+J))
            ELSE
              FACT(NE64+N+J) = ZERO
            ENDIF
 3509     CONTINUE
          K=0
          DO 3510 I=1,N
            IF (PPOS(I).LT.0) THEN
              K = K + 1
              FACT(ISCALE-PPOS(I)-1) = ZERO
            ELSE
              FACT(ISCALE+PPOS(I)-1) =
     &          SQRT(EXP(FACT(NE64+I-K)+FACT(NE64+N+I-K)))
            ENDIF
 3510     CONTINUE
          DO 3516 I = 1,N
            KEEP(PERM+ABS(PPOS(I))-1) = I
 3516     CONTINUE
          K = 1
          DO 3514 JJ = 1,N
            J = PPOS(JJ)
            IF (J.GT.0) THEN
              DO 3511 JLOOP = 1,KEEP(LROW+JJ-1)
                I = IFACT(LIFACT-EXPNE+K)
                INEW = KEEP(PERM+I-1)
                IF (PPOS(INEW).LT.0)
     &            FACT(ISCALE+I-1) = MAX(FACT(ISCALE+I-1),
     &                 ABS(FACT(LLFACT-EXPNE+K))*FACT(ISCALE+J-1))
                K = K + 1
 3511         CONTINUE
            ELSE
              DO 3512 JLOOP = 1,KEEP(LROW+JJ-1)
                I = IFACT(LIFACT-EXPNE+K)
                INEW = KEEP(PERM+I-1)
                IF (I .NE. -J)  THEN
                FACT(ISCALE-J-1) =
     &              MAX(FACT(ISCALE-J-1),
     &              ABS(FACT(LLFACT-EXPNE+K))*FACT(ISCALE+I-1))
                ENDIF
                K = K + 1
 3512         CONTINUE
            ENDIF
 3514     CONTINUE
          DO 3513 I = 1,N
            INEW = KEEP(PERM+I-1)
            IF (PPOS(INEW) .LT. 0) THEN
              PPOS(INEW) = - PPOS(INEW)
              IF (FACT(ISCALE+I-1) .EQ. ZERO) THEN
                FACT(ISCALE+I-1) = ONE
              ELSE
                FACT(ISCALE+I-1) = ONE/FACT(ISCALE+I-1)
              ENDIF
            ENDIF
 3513     CONTINUE
        ENDIF
C8888     CONTINUE
          SMAX = FACT(ISCALE)
          SMIN = FACT(ISCALE)
          DO 5566 I = 1,N
            SMAX = MAX(SMAX,FACT(ISCALE+I-1))
            SMIN = MIN(SMIN,FACT(ISCALE+I-1))
 5566     CONTINUE
          RINFO(16) = SMIN
          RINFO(17) = SMAX
          K = 1
          FACT(BIGA) = ZERO
          DO 6666 JJ = 1,N
            J = PPOS(JJ)
            DO 7777 JLOOP = 1,KEEP(LROW+JJ-1)
              I = IFACT(LIFACT-EXPNE+K)
              FACT(LLFACT-EXPNE+K) =
     &          FACT(ISCALE+I-1)*FACT(LLFACT-EXPNE+K)*FACT(ISCALE+J-1)
              FACT(BIGA) = MAX(FACT(BIGA), ABS(FACT(LLFACT-EXPNE+K)))
              K = K + 1
 7777       CONTINUE
 6666     CONTINUE
CPRINT
C6661     CONTINUE
C6663     CONTINUE
      ENDIF
C**********************************
C**********************************
   22 CALL MA57OD(N, EXPNE, FACT, LLFACT, IFACT, LIFACT, KEEP(LROW),
     *            PPOS,
     *            NSTEPS, KEEP(NSTK), KEEP(NODE), FACT(MM1),
     *            FACT(MM2),
     *            KEEP(PERM),
     *            CNTL, ICNTL,
     *            INFO, RINFO, KEEP(HOLD), FACT(BIGA))
      IF (INFO(1).EQ.10 .OR. INFO(1).EQ.11) THEN
        IF (LDIAG.GT.2 .AND. MP.GE.0)  THEN
          IF (INFO(1).EQ.10) WRITE (MP,99982) INFO(1)
99982 FORMAT (/'Leaving factorization phase (MA57BD) with ...'/
     1  'Factorization suspended because of lack of real space'/
     1  'INFO (1) = ',I3)
          IF (INFO(1).EQ.11) WRITE (MP,99983) INFO(1)
99983 FORMAT (/'Leaving factorization phase (MA57BD) with ...'/
     1  'Factorization suspended because of lack of integer space'/
     1  'INFO (1) = ',I3)
        ENDIF
        RETURN
      ENDIF
      DO 24 I = 1,N
        KEEP(PERM+PPOS(I)-1) = I
   24 CONTINUE
        INFO(17) = ALENB + INFO(17)
        INFO(19) = ALENB + INFO(19)
      IF (ICNTL(15).EQ.1) THEN
        INFO(17) = MAX(INFO(17),ALENB + 3*EXPNE+3*N)
        INFO(19) = MAX(INFO(19),ALENB + 3*EXPNE+3*N)
        INFO(18) = MAX(INFO(18),3*EXPNE+5*N+1)
        INFO(20) = MAX(INFO(18),3*EXPNE+5*N+1)
      ENDIF
      IF (INFO(1).LT.0) RETURN
      GO TO 100
C************************
C************************
   25 INFO(1) = -1
      INFO(2) =  N
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
     +    '**** Error return from MA57BD ****  INFO(1) =',INFO(1),
     +    'N has value ',INFO(2)
      RETURN
   30 INFO(1) = -2
      INFO(2) = NE
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
     +    '**** Error return from MA57BD ****  INFO(1) =',INFO(1),
     +    'NE has value',INFO(2)
      RETURN
   40 INFO(1) = -15
      INFO(2) = LKEEP
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10/A,I10)')
     +    '**** Error return from MA57BD ****  INFO(1) =',INFO(1),
     +    'LKEEP has value    ',INFO(2),
     +    'Should be at least ',5*N+NE+MAX(N,NE)+42
      RETURN
   35 INFO(1) = -10
      INFO(2) = ICNTL(7)
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
     +    '**** Error return from MA57BD ****  INFO(1) =',INFO(1),
     +    'ICNTL(7) has value',ICNTL(7)
      RETURN
   85 INFO(1) = -3
      INFO(2) = LFACT
      INFO(17) = ALENB + EXPNE + 1
      IF (ICNTL(15).EQ.1) INFO(17) = ALENB + 3*EXPNE + 3*N
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
     +    '**** Error return from MA57BD ****  INFO(1) =',INFO(1),
     +    'Insufficient real space in FACT, LFACT = ',INFO(2)
      RETURN
   95 INFO(1) = -4
      INFO(2) = LIFACT
      INFO(18) = EXPNE+N+5
      IF (ICNTL(15).EQ.1) INFO(18) = 3*EXPNE + 5*N + 1
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
     +    '**** Error return from MA57BD ****  INFO(1) =',INFO(1),
     +    'Insufficient integer space in IFACT, LIFACT = ',INFO(2)
      RETURN
C****************
C****************
 100  IF (LDIAG.LE.2 .OR. MP.LT.0) RETURN
      WRITE (MP,99980) INFO(1), INFO(2),
     *    (INFO(I),I=14,25),INFO(28),INFO(29)
      WRITE (MP,99984) (INFO(I),I=31,35),RINFO(3), RINFO(4),
     *                 RINFO(5), RINFO(18)
99980 FORMAT (/'Leaving factorization phase (MA57BD) with ...'/
     1  'INFO (1)                                      =',I12/
     2  ' --- (2)                                      =',I12/
     3  ' --- (14) Number of entries in factors        =',I12/
     4  ' --- (15) Real storage for factors            =',I12/
     5  ' --- (16) Integer storage for factors         =',I12/
     6  ' --- (17) Min LFACT with compresses           =',I12/
     7  ' --- (18) Min LIFACT with compresses          =',I12/
     8  ' --- (19) Min LFACT without compresses        =',I12/
     9  ' --- (20) Min LIFACT without compresses       =',I12/
     *  ' --- (21) Order of largest frontal matrix     =',I12/
     1  ' --- (22) Number of 2x2 pivots                =',I12/
     2  ' --- (23) Number of delayed pivots            =',I12/
     3  ' --- (24) Number of negative eigenvalues      =',I12/
     4  ' --- (25) Rank of factorization               =',I12/
     5  ' --- (28) Number compresses on real data      =',I12/
     6  ' --- (29) Number compresses on integer data   =',I12)
      IF (ICNTL(15).EQ.1) WRITE (MP,99985) RINFO(16),RINFO(17)
99985 FORMAT (
     1  'RINFO(16) Minimum value of scaling factor     =  ',1PD10.3/
     2  '-----(17) Maximum value of scaling factor     =  ',1PD10.3)
99984 FORMAT (
     7  ' --- (31) Number of block pivots in factors   =',I12/
     7  ' --- (32) Number of zeros factors triangle    =',I12/
     7  ' --- (33) Number of zeros factors rectangle   =',I12/
     7  ' --- (34) Number of zero cols factors rect    =',I12/
     7  ' --- (35) Number of static pivots             =',I12/
     1  'RINFO(3)  Operations during node assembly     =  ',1PD10.3/
     2  '-----(4)  Operations during node elimination  =  ',1PD10.3/
     3  '-----(5)  Extra operations because of BLAS    =  ',1PD10.3/
     3  '-----(18) Largest modulus of entry in matrix  =  ',1PD10.3)
      IF (INFO(27).GT.0) WRITE (MP,99981) INFO(27),RINFO(14),RINFO(15)
99981 FORMAT (/'Matrix modification performed'/
     1  'INFO (27) Step at which matrix first modified =',I12/
     2  'RINFO(14) Maximum value added to diagonal     =  ',1PD10.3/
     2  'RINFO(15) Smallest pivot in modified matrix   =  ',1PD10.3)
      CALL MA57UD(FACT,LFACT,IFACT,LIFACT,ICNTL)
      IF (ICNTL(15).NE.1) RETURN
      K = MIN(10,N)
      IF (LDIAG.GE.4) K = N
      WRITE (MP,'(/A/(5D12.5))')  'Scaling factors:',
     +                    (FACT(ISCALE+I-1),I=1,K)
      IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
      END
      SUBROUTINE MA57CD(JOB,N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,W,
     *                  LW,IW1,ICNTL,INFO)
      INTEGER JOB,N,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
      DOUBLE PRECISION W(LW),RHS(LRHS,NRHS)
      INTEGER IW1(N),ICNTL(20),INFO(40)
      INTRINSIC MIN
      EXTERNAL MA57QD,MA57RD,MA57SD,MA57TD,MA57UD,MA57XD,MA57YD
      DOUBLE PRECISION SCALE,ONE
      PARAMETER (ONE = 1.0D0)
      INTEGER I,J,K,LDIAG,LLW,LP,MP,ISCALE
      LP = ICNTL(1)
      MP = ICNTL(3)
      LDIAG = ICNTL(5)
      IF (N.LE.0) THEN
        INFO(1) = -1
        INFO(2) = N
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
     +    '**** Error return from MA57CD ****  INFO(1) =',INFO(1),
     +    'N has value',N
        GOTO 500
      ENDIF
      IF (NRHS.LT.1) THEN
        INFO(1) = -16
        INFO(2) = NRHS
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I10,A)')
     +    '**** Error return from MA57CD ****  INFO(1) =',INFO(1),
     +    'value of NRHS =',NRHS,' is less than 1'
        GOTO 500
      ENDIF
      IF (LRHS.LT.N) THEN
        INFO(1) = -11
        INFO(2) = LRHS
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I10,A,I10)')
     +    '**** Error return from MA57CD ****  INFO(1) =',INFO(1),
     +    'value of LRHS =',LRHS,' is less than N=',N
        GOTO 500
      ENDIF
      IF (LW.LT.N*NRHS) THEN
        INFO(1) = -17
        INFO(2) = N*NRHS
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I10,A,I10)')
     +    '**** Error return from MA57CD ****  INFO(1) =',INFO(1),
     +    'value of LW =',LW,' is less than', N*NRHS
        GOTO 500
      ENDIF
      IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
        WRITE (MP,99999) JOB,N,(ICNTL(I),I=1,5),LFACT,LIFACT,NRHS,
     +         LRHS,LW,ICNTL(13)
99999 FORMAT(/'Entering solution phase (MA57CD) with ...'/
     +    'JOB       Control on coefficient matrix       =',I12/
     +    'N         Order of matrix                     =',I12/
     6    'ICNTL(1)  Stream for errors                   =',I12/
     7    ' --- (2)  Stream for warnings                 =',I12/
     8    ' --- (3)  Stream for monitoring               =',I12/
     9    ' --- (4)  Stream for statistics               =',I12/
     1    ' --- (5)  Level of diagnostic printing        =',I12/
     +    'LFACT     Length of array FACT                =',I12/
     +    'LIFACT    Length of array IFACT               =',I12/
     +    'NRHS      Number of right-hand sides          =',I12/
     +    'LRHS      Leading dimension of RHS array      =',I12/
     +    'LW        Leading dimension of work array     =',I12/
     +    'ICNTL(13) Threshold for Level 2 and 3 BLAS    =',I12)
        CALL MA57UD(FACT,LFACT,IFACT,LIFACT,ICNTL)
        IF (ICNTL(15).EQ.1) THEN
          ISCALE = LFACT-N
          K = MIN(10,N)
          IF (LDIAG.GE.4) K = N
          WRITE (MP,'(/A/(5D12.5))')  'Scaling factors:',
     +                        (FACT(ISCALE+I-1),I=1,K)
          IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
        ENDIF
        K = MIN(10,N)
        IF (LDIAG.GE.4) K = N
        DO 10 J = 1,NRHS
          WRITE(MP,'(/A,I10)') 'Right-hand side',J
          WRITE (MP,'((1P,5D13.3))') (RHS(I,J),I=1,K)
          IF (K.LT.N) WRITE (MP,'(A)') '     . . .'
   10   CONTINUE
      END IF
      LLW = LW/NRHS
      IF (ICNTL(15).EQ.1) THEN
        ISCALE = LFACT-N
        DO 5555 I = 1, N
          SCALE = FACT(ISCALE+I-1)
          IF (JOB.GE.4) SCALE = ONE/FACT(ISCALE+I-1)
          DO 4444 J = 1, NRHS
            RHS(I,J) = SCALE*RHS(I,J)
 4444     CONTINUE
 5555   CONTINUE
      ENDIF
      IF (JOB.LE.2) THEN
        IF (NRHS.EQ.1) THEN
          CALL MA57XD(N,FACT,LFACT,IFACT,LIFACT,RHS,LRHS,
     *                W,LLW,IW1,ICNTL)
        ELSE
          CALL MA57QD(N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *                W,LLW,IW1,ICNTL)
        ENDIF
        IF (JOB.EQ.2) GO TO 15
        IF (NRHS.EQ.1) THEN
          CALL MA57YD(N,FACT,LFACT,IFACT,LIFACT,RHS,LRHS,
     *                W,LLW,IW1,ICNTL)
        ELSE
          CALL MA57RD(N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *                W,LLW,IW1,ICNTL)
        ENDIF
      ENDIF
      IF (JOB.EQ.3)
     *  CALL MA57SD(FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *              W,LLW,ICNTL)
      IF (JOB.GE.4)
     *  CALL MA57TD(N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *              W,LLW,IW1,ICNTL)
   15 IF (ICNTL(15).EQ.1) THEN
        ISCALE = LFACT-N
        DO 6666 I = 1, N
          SCALE = FACT(ISCALE+I-1)
          IF (JOB.EQ.2) SCALE = ONE/FACT(ISCALE+I-1)
          DO 7777 J = 1, NRHS
            RHS(I,J) = SCALE*RHS(I,J)
 7777     CONTINUE
 6666   CONTINUE
      ENDIF
      IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
        WRITE (MP,'(//A)')
     *       'Leaving solution phase (MA57CD) with ...'
        DO 20 J = 1,NRHS
          WRITE(MP,'(/A,I10)') 'Solution       ',J
          WRITE (MP,'(1P,5D13.3)') (RHS(I,J),I=1,K)
          IF (K.LT.N) WRITE (MP,'(A)') '     . . .'
   20   CONTINUE
      ENDIF
  500 RETURN
      END
      SUBROUTINE MA57QD(N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *                  W,LW,IW1,ICNTL)
      INTEGER N,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
      DOUBLE PRECISION W(LW,NRHS),RHS(LRHS,NRHS)
      INTEGER IW1(N),ICNTL(20)
      INTRINSIC ABS
      EXTERNAL DGEMM,DTPSV
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      INTEGER APOS,I,IBLK,II,IPIV,IRHS,IWPOS,J,J1,J2,K,
     +        NCOLS,NROWS
      DOUBLE PRECISION W1
      APOS = 1
      IWPOS = 4
      DO 270 IBLK = 1,IFACT(3)
        IW1(IBLK) = IWPOS
        NCOLS = IFACT(IWPOS)
        NROWS = IFACT(IWPOS+1)
        IWPOS = IWPOS + 2
        IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN
          DO 10 I = 1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            DO 11 J = 1,NRHS
              W(I,J) = RHS(II,J)
   11       CONTINUE
   10     CONTINUE
          DO 12 J = 1,NRHS
            CALL DTPSV('L','N','U',NROWS,FACT(APOS),W(1,J),1)
   12     CONTINUE
          APOS = APOS + (NROWS* (NROWS+1))/2
          IF (NCOLS.GT.NROWS) CALL DGEMM('N','N',NCOLS-NROWS,NRHS,NROWS,
     +                                  ONE,FACT(APOS),NCOLS-NROWS,
     +                                  W,LW,ONE,W(NROWS+1,1),LW)
          APOS = APOS + NROWS* (NCOLS-NROWS)
          DO 35 I = 1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            DO 36 J = 1,NRHS
              RHS(II,J) = W(I,J)
   36       CONTINUE
   35     CONTINUE
        ELSE
        J1 = IWPOS
        J2 = IWPOS + NROWS - 1
        DO 130 IPIV = 1,NROWS
          APOS = APOS + 1
          DO 101 II = 1,NRHS
            W1 = RHS(ABS(IFACT(J1)),II)
            K = APOS
            DO 100 J = J1+1,J2
              IRHS = ABS(IFACT(J))
              RHS(IRHS,II) = RHS(IRHS,II) - FACT(K)*W1
              K = K + 1
  100       CONTINUE
  101     CONTINUE
          APOS = K
          J1 = J1 + 1
  130   CONTINUE
        J2 = IWPOS + NCOLS - 1
        DO 136 IPIV = 1,NROWS
          DO 135 II = 1,NRHS
            K = APOS
            W1 = RHS(ABS(IFACT(IWPOS+IPIV-1)),II)
            DO 133 J = J1,J2
              IRHS = ABS(IFACT(J))
              RHS(IRHS,II) = RHS(IRHS,II) + W1*FACT(K)
              K = K + 1
  133       CONTINUE
  135     CONTINUE
          APOS = K
  136   CONTINUE
      END IF
      IWPOS = IWPOS + NCOLS
  270 CONTINUE
      END
      SUBROUTINE MA57RD(N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *                  W,LW,IW1,ICNTL)
      INTEGER N,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
      DOUBLE PRECISION W(LW,NRHS),RHS(LRHS,NRHS)
      INTEGER IW1(N),ICNTL(20)
      INTRINSIC ABS
      EXTERNAL DGEMM,DTPSV
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      INTEGER APOS,APOS2,I,IBLK,II,IPIV,IRHS,IRHS1,
     +        IRHS2,IWPOS,J,JPIV,J1,J2,K,KK,LROW,NCOLS,NROWS
      DOUBLE PRECISION W1
      APOS = IFACT(1)
      APOS2 = IFACT(2)
      DO 380 IBLK = IFACT(3),1,-1
        IWPOS = IW1(IBLK)
        NCOLS = ABS(IFACT(IWPOS))
        NROWS = ABS(IFACT(IWPOS+1))
        APOS = APOS - NROWS* (NCOLS-NROWS)
        IWPOS = IWPOS + 2
        IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN
          DO 5 I = NROWS + 1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            DO 3 J = 1,NRHS
              W(I,J) = RHS(II,J)
    3       CONTINUE
    5     CONTINUE
          DO 10 IPIV = NROWS,1,-1
            IRHS = ABS(IFACT(IWPOS+IPIV-1))
            APOS = APOS - (NROWS+1-IPIV)
            DO 9 J = 1,NRHS
              W(IPIV,J) = RHS(IRHS,J)*FACT(APOS)
    9       CONTINUE
   10     CONTINUE
          JPIV = -1
          DO 20 IPIV = NROWS,1,-1
            IRHS = IFACT(IWPOS+IPIV-1)
            IF (IRHS.LT.0) THEN
              IRHS1 = -IFACT(IWPOS+IPIV-1+JPIV)
              DO 19 J = 1,NRHS
                W(IPIV,J) = RHS(IRHS1,J)*FACT(APOS2) + W(IPIV,J)
   19         CONTINUE
              IF (JPIV.EQ.1) APOS2 = APOS2 - 1
              JPIV = -JPIV
            END IF
   20     CONTINUE
          K = NCOLS - NROWS
          IF (K.GT.0) CALL DGEMM('T','N',NROWS,NRHS,K,ONE,
     +                           FACT(APOS+(NROWS*(NROWS+1))/2),K,
     +                           W(NROWS+1,1),LW,ONE,W,LW)
          DO 22 J = 1,NRHS
            CALL DTPSV('L','T','U',NROWS,FACT(APOS),W(1,J),1)
   22     CONTINUE
          DO 60 I = 1,NROWS
            II = ABS(IFACT(IWPOS+I-1))
            DO 59 J = 1,NRHS
              RHS(II,J) = W(I,J)
   59       CONTINUE
   60     CONTINUE
        ELSE
          J1 = IWPOS
          J2 = IWPOS + NCOLS - 1
          JPIV = -1
          DO 210 IPIV = NROWS,1,-1
            IRHS = IFACT(IWPOS+IPIV-1)
            LROW = NROWS + 1 - IPIV
            IF (IRHS.GT.0) THEN
              APOS = APOS - LROW
              DO 65 J = 1,NRHS
                RHS(IRHS,J) = RHS(IRHS,J)*FACT(APOS)
   65         CONTINUE
            ELSE
              IF (JPIV.EQ.-1) THEN
                IRHS1 = -IFACT(IWPOS+IPIV-2)
                IRHS2 = -IRHS
                APOS = APOS - LROW - LROW - 1
                DO 68 J = 1,NRHS
                  W1 = RHS(IRHS1,J)*FACT(APOS) +
     +                 RHS(IRHS2,J)*FACT(APOS2)
                  RHS(IRHS2,J) = RHS(IRHS1,J)*FACT(APOS2) +
     +                           RHS(IRHS2,J)*FACT(APOS+LROW+1)
                  RHS(IRHS1,J) = W1
   68           CONTINUE
                APOS2 = APOS2 - 1
              END IF
              JPIV = -JPIV
            END IF
  210     CONTINUE
          APOS = APOS + (NROWS* (NROWS+1))/2
          KK = APOS
          J1 = IWPOS + NROWS
          DO 220 IPIV = 1,NROWS
            IRHS = ABS(IFACT(IWPOS+IPIV-1))
            DO 218 II = 1,NRHS
              W1 = RHS(IRHS,II)
              K = KK
              DO 215 J = J1,J2
                W1 = W1 + FACT(K)*RHS(ABS(IFACT(J)),II)
                K = K + 1
  215         CONTINUE
              RHS(IRHS,II) = W1
  218       CONTINUE
            KK = K
  220     CONTINUE
          J2 = IWPOS + NROWS - 1
          DO 260 IPIV = 1,NROWS
            IRHS = ABS(IFACT(J1-1))
            APOS = APOS - IPIV
            DO 240 II = 1,NRHS
              W1 = RHS(IRHS,II)
              K = APOS + 1
              DO 230 J = J1,J2
                W1 = W1 - FACT(K)*RHS(ABS(IFACT(J)),II)
                K = K + 1
  230         CONTINUE
              RHS(IRHS,II) = W1
  240       CONTINUE
            J1 = J1 - 1
  260     CONTINUE
        END IF
  380 CONTINUE
      END
      SUBROUTINE MA57UD(FACT,LFACT,IFACT,LIFACT,ICNTL)
      INTEGER LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),ICNTL(20)
      INTRINSIC MIN,SIGN
      CHARACTER*72 LINE
      INTEGER APOS,APOS2,IBLK,ILINE,IROW,IWPOS,J,JPIV,J1,J2,K,
     +        LDIAG,LEN,MP,NBLK,NCOLS,NROWS
      CHARACTER*1 PM(-2:2)
      DATA PM/'*','-','.','+','.'/
      DOUBLE PRECISION ZERO,TINY,FD05AD
      PARAMETER (ZERO=0.0D0)
      EXTERNAL FD05AD
      MP = ICNTL(3)
      LDIAG = ICNTL(5)
      TINY = FD05AD(4)
      APOS2 = IFACT(1)
      NBLK = IFACT(3)
      IF (LDIAG.EQ.3) NBLK = MIN(1,NBLK)
      LEN = 12
      IF (LDIAG.EQ.5) LEN = 1
      IF (LEN.EQ.12) THEN
        IF (NBLK.EQ.IFACT(3)) THEN
          WRITE (MP,'(/A)')
     +      'For each block, the following information is provided:'
        ELSE
          WRITE (MP,'(/A,A)') 'For the first block only,',
     +      ' the following information is provided:'
        END IF
      END IF
      IF (LEN.EQ.12) WRITE (MP,'(A)')
     +    '   1. Block number, number of rows, number of columns',
     +    '   2. List of indices for the pivot, each negated if part of'
     +    ,'      a 2x2 pivot',
     +    '   3. The factorized block pivot',
     +    '      It has the form',
     +    '            -1  T',
     +    '        L  D   L ',
     +    '                         -1    T',
     +    '      and is printed as D and L  packed together.',
     +    '   4. List of indices for the non-pivot columns',
     +    '   5. The non-pivot part as rectangular block by rows'
      IWPOS = 4
      APOS = 1
      DO 300 IBLK = 1,NBLK
        NCOLS = IFACT(IWPOS)
        NROWS = IFACT(IWPOS+1)
        IWPOS = IWPOS + 2
        WRITE (MP,'(/4(A,I6))') 'Block pivot',IBLK,' with',NROWS,
     +        ' rows and', NCOLS,' columns'
        IF (LEN.EQ.12) WRITE (MP,'(6I12)')
     +                       (IFACT(K),K=IWPOS,IWPOS+NROWS-1)
        IF (LEN.EQ.1) WRITE (MP,'(72A1)') (PM(SIGN(1,IFACT(K))),
     +      K=IWPOS,IWPOS+NROWS-1)
        JPIV = 0
        DO 30 IROW = 1,NROWS
          IF (JPIV.EQ.1) THEN
            JPIV = 0
          ELSE
            IF (IFACT(IWPOS+IROW-1).LT.0) JPIV = 1
          END IF
          ILINE = 1
          DO 10 J = 1,IROW - 1
            WRITE (LINE(ILINE:ILINE+LEN-1),'(A)') ' '
            ILINE = ILINE + LEN
            IF (ILINE.GT.72) THEN
              WRITE (MP,'(A)') LINE
              ILINE = 1
            END IF
   10     CONTINUE
          DO 20 J = IROW,NROWS
            IF (LEN.EQ.12) WRITE (LINE(ILINE:ILINE+11),
     +          '(1P,D12.4)') FACT(APOS)
            IF (LEN.EQ.1) THEN
               IF (FACT(APOS).EQ.ZERO) THEN
                  WRITE (LINE(ILINE:ILINE),'(A)') '.'
               ELSE
                  WRITE (LINE(ILINE:ILINE),'(A)') '*'
               END IF
            END IF
            APOS = APOS + 1
            IF (J.EQ.IROW+1) THEN
              IF (JPIV.EQ.1) THEN
                IF (LEN.EQ.12) WRITE (LINE(ILINE:ILINE+11),
     +              '(1P,D12.4)') FACT(APOS2)
                IF (LEN.EQ.1) THEN
                    IF (FACT(APOS2).EQ.ZERO) THEN
                       WRITE (LINE(ILINE:ILINE),'(A)') '.'
                    ELSE
                       WRITE (LINE(ILINE:ILINE),'(A)') '*'
                    END IF
                END IF
                APOS2 = APOS2 + 1
              END IF
            END IF
            ILINE = ILINE + LEN
            IF (ILINE.GT.72) THEN
              WRITE (MP,'(A)') LINE
              ILINE = 1
            END IF
   20     CONTINUE
          IF (ILINE.GT.1) THEN
            LINE(ILINE:) = ' '
            WRITE (MP,'(A)') LINE
          END IF
   30   CONTINUE
        IWPOS = IWPOS + NROWS
        IF (LEN.EQ.12) WRITE (MP,'(6I12)') (IFACT(K),K=IWPOS,
     +      IWPOS+NCOLS-NROWS-1)
        IF (LEN.EQ.1) WRITE (MP,'(72A1)') (PM(SIGN(1,IFACT(K))),
     +      K=IWPOS,IWPOS+NCOLS-NROWS-1)
        IWPOS = IWPOS + NCOLS - NROWS
        DO 280 IROW = 1,NROWS
          J1 = NROWS
          J2 = NCOLS
          ILINE = 1
          DO 110 J = J1 + 1,J2
            IF (LEN.EQ.12) WRITE (LINE(ILINE:ILINE+11),
     +          '(1P,D12.4)') FACT(APOS)
            IF (LEN.EQ.1) THEN
               IF (FACT(APOS).EQ.ZERO) THEN
                  WRITE (LINE(ILINE:ILINE),'(A)') '.'
               ELSE
                  WRITE (LINE(ILINE:ILINE),'(A)') '*'
               END IF
            END IF
            APOS = APOS + 1
            ILINE = ILINE + LEN
            IF (ILINE.GT.72) THEN
              WRITE (MP,'(A)') LINE
              ILINE = 1
            END IF
  110     CONTINUE
          IF (ILINE.GT.1) THEN
            LINE(ILINE:) = ' '
            WRITE (MP,'(A)') LINE
          END IF
  280   CONTINUE
  300 CONTINUE
      END
      SUBROUTINE MA57SD(FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *                  W,LW,ICNTL)
      INTEGER LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
      DOUBLE PRECISION W(LW,NRHS),RHS(LRHS,NRHS)
      INTEGER ICNTL(20)
      INTRINSIC ABS
      EXTERNAL DGEMM,DTPSV
      INTEGER APOS,APOS2,I,IBLK,II,IPIV,IRHS,IRHS1,
     +        IRHS2,IWPOS,J,JPIV,NCOLS,NROWS
      DOUBLE PRECISION W1
      APOS = 1
      APOS2 = IFACT(1)
      IWPOS = 4
      DO 380 IBLK = 1,IFACT(3)
        NCOLS = IFACT(IWPOS)
        NROWS = IFACT(IWPOS+1)
        IWPOS = IWPOS + 2
        IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN
          DO 10 IPIV = 1,NROWS
            IRHS = ABS(IFACT(IWPOS+IPIV-1))
            DO 9 J = 1,NRHS
              W(IPIV,J) = RHS(IRHS,J)*FACT(APOS)
    9       CONTINUE
            APOS = APOS + (NROWS+1-IPIV)
   10     CONTINUE
          JPIV = 1
          DO 20 IPIV = 1,NROWS
            IRHS = IFACT(IWPOS+IPIV-1)
            IF (IRHS.LT.0) THEN
              IRHS1 = -IFACT(IWPOS+IPIV-1+JPIV)
              DO 19 J = 1,NRHS
                W(IPIV,J) = RHS(IRHS1,J)*FACT(APOS2) + W(IPIV,J)
   19         CONTINUE
              IF (JPIV.EQ.-1) APOS2 = APOS2 + 1
              JPIV = -JPIV
            END IF
   20     CONTINUE
          DO 60 I = 1,NROWS
            II = ABS(IFACT(IWPOS+I-1))
            DO 59 J = 1,NRHS
              RHS(II,J) = W(I,J)
   59       CONTINUE
   60     CONTINUE
        ELSE
          JPIV = 1
          DO 210 IPIV = 1,NROWS
            IRHS = IFACT(IWPOS+IPIV-1)
            IF (IRHS.GT.0) THEN
              DO 65 J = 1,NRHS
                RHS(IRHS,J) = RHS(IRHS,J)*FACT(APOS)
   65         CONTINUE
              APOS = APOS + NROWS - IPIV + 1
            ELSE
              IF (JPIV.EQ.1) THEN
                IRHS1 = -IRHS
                IRHS2 = -IFACT(IWPOS+IPIV)
                DO 68 J = 1,NRHS
                  W1 = RHS(IRHS1,J)*FACT(APOS) +
     +                 RHS(IRHS2,J)*FACT(APOS2)
                  RHS(IRHS2,J) = RHS(IRHS1,J)*FACT(APOS2) +
     +                           RHS(IRHS2,J)*FACT(APOS+NROWS-IPIV+1)
                  RHS(IRHS1,J) = W1
   68           CONTINUE
                APOS2 = APOS2 + 1
              END IF
              JPIV = -JPIV
              APOS = APOS + NROWS - IPIV + 1
            END IF
  210     CONTINUE
        END IF
        IWPOS = IWPOS + NCOLS
        APOS = APOS + NROWS*(NCOLS-NROWS)
  380 CONTINUE
      END
      SUBROUTINE MA57TD(N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *                  W,LW,IW1,ICNTL)
      INTEGER N,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
      DOUBLE PRECISION W(LW,NRHS),RHS(LRHS,NRHS)
      INTEGER IW1(N),ICNTL(20)
      INTRINSIC ABS
      EXTERNAL DGEMM,DTPSV
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      INTEGER APOS,I,IBLK,II,IPIV,IRHS,
     +        IWPOS,J,J1,J2,K,KK,NCOLS,NROWS
      DOUBLE PRECISION W1
      APOS = IFACT(1)
      IWPOS = 4
      DO 10 I = 1,IFACT(3)-1
        IW1(I) = IWPOS
        IWPOS = IWPOS + ABS(IFACT(IWPOS))+2
   10 CONTINUE
      IW1(IFACT(3)) = IWPOS
      DO 380 IBLK = IFACT(3),1,-1
        IWPOS = IW1(IBLK)
        NCOLS = ABS(IFACT(IWPOS))
        NROWS = ABS(IFACT(IWPOS+1))
        APOS = APOS - NROWS* (NCOLS-NROWS)
        IWPOS = IWPOS + 2
        IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN
          DO 5 I = 1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            DO 3 J = 1,NRHS
              W(I,J) = RHS(II,J)
    3       CONTINUE
    5     CONTINUE
          K = NCOLS - NROWS
          IF (K.GT.0) CALL DGEMM('T','N',NROWS,NRHS,K,ONE,
     +                           FACT(APOS),K,
     +                           W(NROWS+1,1),LW,ONE,W,LW)
          APOS = APOS-(NROWS*(NROWS+1))/2
          DO 22 J = 1,NRHS
            CALL DTPSV('L','T','U',NROWS,FACT(APOS),W(1,J),1)
   22     CONTINUE
          DO 60 I = 1,NROWS
            II = ABS(IFACT(IWPOS+I-1))
            DO 59 J = 1,NRHS
              RHS(II,J) = W(I,J)
   59       CONTINUE
   60     CONTINUE
        ELSE
          J1 = IWPOS
          J2 = IWPOS + NCOLS - 1
          KK = APOS
          J1 = IWPOS + NROWS
          DO 220 IPIV = 1,NROWS
            IRHS = ABS(IFACT(IWPOS+IPIV-1))
            DO 218 II = 1,NRHS
              W1 = RHS(IRHS,II)
              K = KK
              DO 215 J = J1,J2
                W1 = W1 + FACT(K)*RHS(ABS(IFACT(J)),II)
                K = K + 1
  215         CONTINUE
              RHS(IRHS,II) = W1
  218       CONTINUE
            KK = K
  220     CONTINUE
          J2 = IWPOS + NROWS - 1
          DO 260 IPIV = 1,NROWS
            IRHS = ABS(IFACT(J1-1))
            APOS = APOS - IPIV
            DO 240 II = 1,NRHS
              W1 = RHS(IRHS,II)
              K = APOS + 1
              DO 230 J = J1,J2
                W1 = W1 - FACT(K)*RHS(ABS(IFACT(J)),II)
                K = K + 1
  230         CONTINUE
              RHS(IRHS,II) = W1
  240       CONTINUE
            J1 = J1 - 1
  260     CONTINUE
        END IF
  380 CONTINUE
      END
      SUBROUTINE MA57DD(JOB,N,NE,A,IRN,JCN,FACT,LFACT,IFACT,LIFACT,
     *                  RHS,X,RESID,W,IW,ICNTL,CNTL,INFO,RINFO)
      INTEGER JOB,N,NE
      DOUBLE PRECISION A(NE)
      INTEGER IRN(NE),JCN(NE),LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT)
      DOUBLE PRECISION RHS(N),X(N),RESID(N),W(N,*)
      INTEGER IW(N),ICNTL(20)
      DOUBLE PRECISION CNTL(5)
      INTEGER INFO(40)
      DOUBLE PRECISION RINFO(20)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.D0,ONE=1.0D0)
      DOUBLE PRECISION COND(2),CTAU,DXMAX,ERROR,OLDOMG(2),OLDOM2,
     *                 OMEGA(2),OM2,TAU
      INTEGER I,ICNTLC(20),ITER,J,K,KASE,KK,LDIAG,LP,MP,KEEP71(5)
      LOGICAL LCOND(2)
      INTRINSIC MIN
      EXTERNAL MA57CD,MA57UD,FD05AD,MC71AD
      DOUBLE PRECISION EPS,FD05AD
      LP = ICNTL(1)
      MP = ICNTL(3)
      LDIAG = ICNTL(5)
      INFO(1) = 0
      IF (N.LE.0) THEN
        INFO(1) = -1
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I12)')
     +    '**** Error return from MA57DD ****  INFO(1) =',INFO(1),
     +    'N has value',N
        INFO(2) = N
        GOTO 500
      ENDIF
      IF (NE.LT.0) THEN
        INFO(1) = -2
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I12)')
     +    '**** Error return from MA57DD ****  INFO(1) =',INFO(1),
     +    'NE has value',NE
        INFO(2) = NE
        GOTO 500
      ENDIF
      IF (ICNTL(9).LT.1) THEN
        INFO(1) = -13
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I12)')
     +    '**** Error return from MA57DD ****  INFO(1) =',INFO(1),
     +    'ICNTL(9) has value',ICNTL(9)
        INFO(2) = ICNTL(9)
        GOTO 500
      ENDIF
      IF (JOB.LT.0 .OR. JOB.GT.4 .OR. (ICNTL(9).GT.1 .AND.
     *    (JOB.NE.0 .AND. JOB.NE.2)))  THEN
        INFO(1) = -12
        INFO(2) = JOB
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I12)')
     +    '**** Error return from MA57DD ****  INFO(1) =',INFO(1),
     +    'JOB has value',JOB
        IF (ICNTL(9).GT.1 .AND. LDIAG.GT.0 .AND. LP.GE.0)
     +    WRITE (LP,'(A,I3)') 'and ICNTL(9) =',ICNTL(9)
        GOTO 500
      ENDIF
      IF (NE.EQ.0) THEN
        IF (JOB.NE.3) THEN
          DO 8 I = 1,N
            RESID(I) = ZERO
  8       CONTINUE
        ENDIF
        DO 9 I = 1,N
          X(I) = ZERO
  9     CONTINUE
        INFO(30)=0
        DO 10 I = 6,13
          RINFO(I) = ZERO
 10     CONTINUE
        GO TO 500
      ENDIF
      IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
        WRITE (MP,99999) JOB,N,NE,(ICNTL(I),I=1,5),LFACT,LIFACT,
     +   ICNTL(9),ICNTL(10),ICNTL(13),CNTL(3)
99999 FORMAT(/'Entering iterative refinement solution phase ',
     +  '(MA57DD) with ...'/
     +  'JOB       Control for coefficient matrix      =',I12/
     +  'N         Order of matrix                     =',I12/
     +  'NE        Number of entries in matrix         =',I12/
     6  'ICNTL(1)  Stream for errors                   =',I12/
     7  ' --- (2)  Stream for warnings                 =',I12/
     8  ' --- (3)  Stream for monitoring               =',I12/
     9  ' --- (4)  Stream for statistics               =',I12/
     1  ' --- (5)  Level of diagnostic printing        =',I12/
     +  'LFACT     Length of array FACT                =',I12/
     +  'LIFACT    Length of array IFACT               =',I12/
     +  'ICNTL(9)  Number steps iterative refinement   =',I12/
     +  'ICNTL(10) Control for error analysis          =',I12/
     +  'ICNTL(13) Threshold for Level 2 and 3 BLAS    =',I12/
     +  'CNTL(3)   Convergence test for IR             =',1P,D12.4)
        CALL MA57UD(FACT,LFACT,IFACT,LIFACT,ICNTL)
        K = MIN(10,N)
        IF (LDIAG.GE.4) K = N
        WRITE(MP,'(/A)') 'Right-hand side'
        WRITE (MP,'((4X, 1P,5D13.3))') (RHS(I),I=1,K)
        IF (K.LT.N) WRITE (MP,'(A)') '     . . .'
      END IF
      DO 15 I=1,5
        ICNTLC(I) = ICNTL(I)
   15 CONTINUE
      ICNTLC(13) = ICNTL(13)
      ICNTLC(15) = ICNTL(15)
      ICNTLC(3) = -1
      IF (JOB.LE.2) THEN
        IF (JOB .LE. 1) THEN
          DO 14 I = 1,N
            X(I) = RHS(I)
            RESID(I) = RHS(I)
   14     CONTINUE
          CALL MA57CD(1,N,FACT,LFACT,IFACT,LIFACT,1,X,N,W,N,IW,
     +                ICNTLC,INFO)
        ELSE
          DO 13 I = 1,N
            RESID(I) = RHS(I)
   13     CONTINUE
        ENDIF
        IF (ICNTL(9).EQ.1) THEN
          DO 16 KK = 1,NE
            I = IRN(KK)
            J = JCN(KK)
            IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 16
            RESID(J) = RESID(J) - A(KK)*X(I)
            IF (I.NE.J) RESID(I) = RESID(I) - A(KK)*X(J)
   16     CONTINUE
          IF (JOB.EQ.0) GO to 340
        ELSE
          DO 18 I = 1,N
            W(I,1) = ZERO
            W(I,3) = ZERO
   18     CONTINUE
          DO 17 KK = 1,NE
            I = IRN(KK)
            J = JCN(KK)
            IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 17
            RESID(J) = RESID(J) - A(KK)*X(I)
            W(J,1) = W(J,1) + ABS(A(KK)*X(I))
            W(J,3) = W(J,3) + ABS(A(KK))
            IF (I.NE.J) THEN
              RESID(I) = RESID(I) - A(KK)*X(J)
              W(I,1) = W(I,1) + ABS(A(KK)*X(J))
              W(I,3) = W(I,3) + ABS(A(KK))
            ENDIF
   17     CONTINUE
        DXMAX = ZERO
        DO 221 I = 1,N
          DXMAX = MAX(DXMAX,ABS(X(I)))
  221   CONTINUE
      EPS = FD05AD(1)
        CTAU = 1000.*EPS
          OMEGA(1) = ZERO
          OMEGA(2) = ZERO
          DO 231 I = 1,N
            TAU = (W(I,3)*DXMAX+ABS(RHS(I)))*N*CTAU
            IF ((W(I,1)+ABS(RHS(I))).GT.TAU) THEN
              OMEGA(1) = MAX(OMEGA(1),ABS(RESID(I))/
     +                   (W(I,1)+ABS(RHS(I))))
              IW(I) = 1
            ELSE
              IF (TAU.GT.ZERO) THEN
                OMEGA(2) = MAX(OMEGA(2),ABS(RESID(I))/
     +                     (W(I,1)+W(I,3)*DXMAX))
              END IF
              IW(I) = 2
            END IF
  231     CONTINUE
          OM2 = OMEGA(1) + OMEGA(2)
          ITER = 0
          IF (OM2.LE.EPS) THEN
            GO TO 270
          ENDIF
          DO 251 I = 1,N
            W(I,2) = X(I)
  251     CONTINUE
          OLDOMG(1) = OMEGA(1)
          OLDOMG(2) = OMEGA(2)
          OLDOM2 = OM2
        ENDIF
      ENDIF
      DO 260 ITER = 1,ICNTL(9)
        CALL MA57CD(1,N,FACT,LFACT,IFACT,LIFACT,1,RESID,N,W,N,IW,
     +              ICNTLC,INFO)
        DO 141 I = 1,N
          X(I) = X(I) + RESID(I)
  141   CONTINUE
        IF (JOB.LT.4 .AND. ICNTL(9).EQ.1) GO TO 340
        IF (ICNTL(9).EQ.1) THEN
          DO 151 I = 1,N
            RESID(I) = RHS(I)
  151     CONTINUE
          DO 181 KK = 1,NE
            I = IRN(KK)
            J = JCN(KK)
            IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 181
            RESID(J) = RESID(J) - A(KK)*X(I)
            IF (I.NE.J) RESID(I) = RESID(I) - A(KK)*X(J)
  181     CONTINUE
          GO TO 340
        ELSE
          DO 153 I = 1,N
            RESID(I) = RHS(I)
            W(I,1) = ZERO
  153     CONTINUE
          DO 183 KK = 1,NE
            I = IRN(KK)
            J = JCN(KK)
            IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 183
            RESID(J) = RESID(J) - A(KK)*X(I)
            W(J,1) = W(J,1) + ABS(A(KK)*X(I))
            IF (I.NE.J) THEN
              RESID(I) = RESID(I) - A(KK)*X(J)
              W(I,1) = W(I,1) + ABS(A(KK)*X(J))
            ENDIF
  183     CONTINUE
        ENDIF
        DXMAX = ZERO
        DO 220 I = 1,N
          DXMAX = MAX(DXMAX,ABS(X(I)))
  220   CONTINUE
        OMEGA(1) = ZERO
        OMEGA(2) = ZERO
        DO 230 I = 1,N
          TAU = (W(I,3)*DXMAX+ABS(RHS(I)))*N*CTAU
          IF ((W(I,1)+ABS(RHS(I))).GT.TAU) THEN
            OMEGA(1) = MAX(OMEGA(1),ABS(RESID(I))/
     +                 (W(I,1)+ABS(RHS(I))))
            IW(I) = 1
          ELSE
            IF (TAU.GT.ZERO) THEN
              OMEGA(2) = MAX(OMEGA(2),ABS(RESID(I))/
     +                   (W(I,1)+W(I,3)*DXMAX))
            END IF
            IW(I) = 2
          END IF
  230   CONTINUE
        OM2 = OMEGA(1) + OMEGA(2)
        IF ((OM2+ONE).LE.ONE) THEN
          GO TO 270
        ENDIF
        IF (OM2.GT.OLDOM2*CNTL(3)) THEN
          IF (OM2.GT.OLDOM2) THEN
            OMEGA(1) = OLDOMG(1)
            OMEGA(2) = OLDOMG(2)
            DO 240 I = 1,N
              X(I) = W(I,2)
  240       CONTINUE
          END IF
          GO TO 270
        ELSE
          DO 250 I = 1,N
            W(I,2) = X(I)
  250     CONTINUE
          OLDOMG(1) = OMEGA(1)
          OLDOMG(2) = OMEGA(2)
          OLDOM2 = OM2
        END IF
  260 CONTINUE
      INFO(1) = -8
      IF (LP.GE.0 .AND. LDIAG.GE.1) WRITE (LP,9170) INFO(1),ICNTL(9)
 9170 FORMAT ('Error return from MA57D/DD because of ','nonconvergence',
     +       ' of iterative refinement'/'Error INFO(1) = ',I2,'  with',
     +       ' ICNTL','(9) = ',I10)
  270 RINFO(6)  = OMEGA(1)
      RINFO(7)  = OMEGA(2)
      RINFO(8) = ZERO
      DO 271 I=1,N
        RINFO(8) = MAX(RINFO(8),W(I,3))
  271 CONTINUE
      RINFO(9) = DXMAX
      RINFO(10) = ZERO
      DO 272 I=1,N
        RINFO(10) = MAX(RINFO(10),ABS(RESID(I)))
  272 CONTINUE
      IF (RINFO(8)*RINFO(9).NE.ZERO)
     *RINFO(10) = RINFO(10)/(RINFO(8)*RINFO(9))
      INFO(30) = ITER
      IF (INFO(1).LT.0) GO TO 340
      IF (ICNTL(10).LE.0) GO TO 340
      LCOND(1) = .FALSE.
      LCOND(2) = .FALSE.
      ERROR    = ZERO
      DO 280 I = 1,N
        IF (IW(I).EQ.1) THEN
          W(I,1) = W(I,1) + ABS(RHS(I))
          W(I,2) = ZERO
          LCOND(1) = .TRUE.
        ELSE
          W(I,2) = W(I,1) + W(I,3)*DXMAX
          W(I,1) = ZERO
          LCOND(2) = .TRUE.
        END IF
  280 CONTINUE
      DO 330 K = 1,2
        IF (LCOND(K)) THEN
          KASE = 0
          DO 310 KK = 1,40
            CALL MC71AD(N,KASE,W(1,3),COND(K),W(1,4),IW,KEEP71)
            IF (KASE.EQ.0) GO TO 320
            IF (KASE.EQ.1) THEN
              CALL MA57CD(1,N,FACT,LFACT,IFACT,LIFACT,1,W(1,3),
     *                    N,W(1,4),N,IW,ICNTLC,INFO)
              DO 290 I = 1,N
                W(I,3) = W(I,K)*W(I,3)
  290         CONTINUE
            END IF
            IF (KASE.EQ.2) THEN
              DO 300 I = 1,N
                W(I,3) = W(I,K)*W(I,3)
  300         CONTINUE
              CALL MA57CD(1,N,FACT,LFACT,IFACT,LIFACT,1,W(1,3),N,
     *                    W(1,4),N,IW,ICNTLC,INFO)
            END IF
  310     CONTINUE
          INFO(1) = -14
          IF (LP.GE.0 .AND. LDIAG.GE.1) WRITE (LP,9160)
 9160 FORMAT ('Error return from MA57D/DD because of ','error in MC71',
     +       'A/AD'/'Error not calculated')
  320     IF (DXMAX.GT.ZERO) COND(K) = COND(K)/DXMAX
          ERROR = ERROR + OMEGA(K)*COND(K)
        ELSE
          COND(K) = ZERO
        ENDIF
  330 CONTINUE
      RINFO(11)  = COND(1)
      RINFO(12)  = COND(2)
      RINFO(13)  = ERROR
 340  IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
        WRITE(MP,99980) INFO(1)
99980 FORMAT (/'Leaving iterative refinement solution phase ',
     +  '(MA57DD) with ...'/
     1      'INFO (1)                                      =',I12/)
        IF (INFO(1).LT.0) GO TO 500
        IF (ICNTL(9).GT.1) THEN
          WRITE(MP,99981) INFO(30),(RINFO(I),I=6,10)
99981     FORMAT(
     1     'INFO(30)  Number steps iterative ref   =',I10/
     1     'RINFO(6)  Backward errors  (OMEGA(1))  =',1PD10.3/
     2     '-----(7)  Backward errors  (OMEGA(2))  =',1PD10.3/
     3     '-----(8)  Infinity norm of matrix      =',1PD10.3/
     4     '-----(9)  Infinity norm of solution    =',1PD10.3/
     5     '-----(10) Norm of scaled residuals     =',1PD10.3)
          IF (ICNTL(10).GT.0) WRITE(MP,99979) (RINFO(I),I=11,13)
99979       FORMAT (
     1       'RINFO(11) Condition number (COND(1))   =',1PD10.3/
     1       'RINFO(12) Condition number (COND(2))   =',1PD10.3/
     1       'RINFO(13) Error in solution            =',1PD10.3)
          WRITE(MP,'(/A,I10)') 'Residual'
          K=MIN(N,10)
          IF (LDIAG.GE.4) K = N
          WRITE (MP,'(1P,5D13.3)') (RESID(I),I=1,K)
          IF (K.LT.N) WRITE (MP,'(A)') '     . . .'
        ELSE
          IF (JOB.GE.1 .AND. JOB.LE.3) THEN
            WRITE(MP,'(/A,I10)') 'Correction to solution'
          ELSE
            WRITE(MP,'(/A,I10)') 'Residual'
          ENDIF
          K=MIN(N,10)
          IF (LDIAG.GE.4) K = N
          WRITE (MP,'(1P,5D13.3)') (RESID(I),I=1,K)
          IF (K.LT.N) WRITE (MP,'(A)') '     . . .'
        END IF
        K=MIN(N,10)
        IF (LDIAG.GE.4) K = N
        WRITE(MP,'(/A,I10)') 'Solution'
        WRITE (MP,'(1P,5D13.3)') (X(I),I=1,K)
        IF (K.LT.N) WRITE (MP,'(A)') '     . . .'
      END IF
 500  RETURN
      END
      SUBROUTINE MA57ED(N,IC,KEEP,FACT,LFACT,NEWFAC,LNEW,
     *                  IFACT,LIFACT,NEWIFC,LINEW,INFO)
      INTEGER N,IC,KEEP(*),LFACT,LNEW,LIFACT,LINEW,INFO(40)
      DOUBLE PRECISION FACT(LFACT),NEWFAC(LNEW)
      INTEGER IFACT(LIFACT),NEWIFC(LINEW)
      INTEGER APOSBB,ASTK,HOLD,I,ISTK,IWPOS,MOVE,NFRONT
      HOLD = N + 3
      INFO(1) = 0
      INFO(2) = 0
      IF (IC.GE.1) THEN
        IF (LINEW.LE.LIFACT) THEN
          INFO(1) = -7
          INFO(2) = LINEW
          RETURN
        ENDIF
        IWPOS = KEEP(HOLD+7)
        ISTK  = KEEP(HOLD+14)
        NFRONT = KEEP(HOLD+23)
        DO 10 I = 1,IWPOS+NFRONT-1
          NEWIFC(I) = IFACT(I)
   10   CONTINUE
        MOVE = LINEW - LIFACT
        DO 20 I = ISTK+1,LIFACT
          NEWIFC(I+MOVE) = IFACT(I)
   20   CONTINUE
          KEEP(HOLD+13) = KEEP(HOLD+13) + MOVE
          KEEP(HOLD+14) = ISTK + MOVE
          KEEP(HOLD+18) = KEEP(HOLD+18) + MOVE
      ENDIF
      IF (IC.NE.1) THEN
        IF (LNEW.LE.LFACT) THEN
          INFO(1) = -7
          INFO(2) = LNEW
          RETURN
        ENDIF
        APOSBB = KEEP(HOLD+9)
        ASTK   = KEEP(HOLD+15)
        DO 60 I = 1, APOSBB-1
          NEWFAC(I) = FACT(I)
   60   CONTINUE
        MOVE = LNEW - LFACT
        DO 70 I = ASTK+1,LFACT
          NEWFAC(I+MOVE) = FACT(I)
   70   CONTINUE
        KEEP(HOLD+12) = KEEP(HOLD+12) + MOVE
        KEEP(HOLD+15) = ASTK + MOVE
        KEEP(HOLD+19) = KEEP(HOLD+19) + MOVE
      ENDIF
      RETURN
      END
      SUBROUTINE MA57GD(N,NE,IRN,JCN,IW,IPE,COUNT,FLAG,IWFR,ICNTL,INFO)
      INTEGER N,NE,IRN(NE),JCN(NE),IW(NE*2+N),IPE(N),COUNT(N),
     +        FLAG(N),IWFR,ICNTL(20),INFO(40)
      INTRINSIC MAX,MIN
      INTEGER I,J,K,L,LDIAG,WP
      WP = ICNTL(2)
      LDIAG = ICNTL(5)
      IF (WP.LT.0) LDIAG = 0
      INFO(1) = 0
      INFO(3) = 0
      DO 10 I = 1,N
        FLAG(I) = 0
        COUNT(I) = 0
   10 CONTINUE
      DO 20 K = 1,NE
        I = IRN(K)
        J = JCN(K)
        IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) THEN
          INFO(3) = INFO(3) + 1
          INFO(1) = 1
          IF (INFO(3).EQ.1 .AND. LDIAG.GT.1) WRITE (WP,'(2A,I2)')
     +        '*** Warning message from subroutine MA57AD ***',
     +        ' INFO(1) =',INFO(1)
          IF (INFO(3).LE.10 .AND. LDIAG.GT.1) WRITE (WP,'(3(I10,A))')
     +         K,'th entry (in row',I,' and column',J,') ignored'
        ELSE IF (I.NE.J) THEN
          COUNT(I) = COUNT(I) + 1
          COUNT(J) = COUNT(J) + 1
        END IF
   20 CONTINUE
      IPE(1) = COUNT(1)+1
      DO 30 I = 2,N
        IPE(I) = IPE(I-1) + COUNT(I)
   30 CONTINUE
      DO 40 K = 1,NE
        I = IRN(K)
        J = JCN(K)
        IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N .OR. I.EQ.J) GO TO 40
        IPE(I) = IPE(I) - 1
        IW(IPE(I)) = J
        IPE(J) = IPE(J) - 1
        IW(IPE(J)) = I
   40 CONTINUE
      INFO(4) = 0
      IWFR = 1
      DO 60 I = 1,N
        L = IPE(I)
        IPE(I) = IWFR
        DO 50 K = L,L+COUNT(I)-1
          J = IW(K)
          IF (FLAG(J).NE.I) THEN
            FLAG(J) = I
            IW(IWFR) = J
            IWFR = IWFR + 1
          ELSE
            IF (I.LT.J) INFO(4) = INFO(4) + 1
          END IF
   50   CONTINUE
        COUNT(I) = IWFR - IPE(I)
   60 CONTINUE
      IF (INFO(4).GT.0) THEN
        INFO(1) = INFO(1) + 2
        IF (LDIAG.GT.1 .AND. WP.GE.0) WRITE (WP,'(A/I10,A)')
     +      '*** Warning message from subroutine MA57AD ***',INFO(4),
     +      ' off-diagonal duplicate entries found'
      END IF
      END
      SUBROUTINE MA57JD(N,NE,IRN,JCN,PERM,IW,IPE,COUNT,FLAG,IWFR,
     +                  ICNTL,INFO)
      INTEGER N,NE,IRN(NE),JCN(NE),IW(NE+N),IPE(N),COUNT(N),
     +        PERM(N),FLAG(N),IWFR,ICNTL(20),INFO(40)
      INTRINSIC MAX,MIN
      INTEGER I,J,K,L,LDIAG,WP
      WP = ICNTL(2)
      LDIAG = ICNTL(5)
      IF (WP.LT.0) LDIAG = 0
      INFO(1) = 0
      DO 10 I = 1,N
        FLAG(I) = 0
        COUNT(I) = 0
   10 CONTINUE
      INFO(3) = 0
      DO 30 K = 1,NE
        I = IRN(K)
        J = JCN(K)
        IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) THEN
          IRN(K) = 0
          JCN(K) = 0
          INFO(3) = INFO(3) + 1
          INFO(1) = 1
          IF (INFO(3).EQ.1 .AND. LDIAG.GT.1) WRITE (WP,'(2A,I2)')
     +        '*** Warning message from subroutine MA57AD ***',
     +        ' INFO(1) =',INFO(1)
          IF (INFO(3).LE.10 .AND. LDIAG.GT.1) WRITE (WP,'(3(I10,A))')
     +        K,'th entry (in row',I,' and column',J,') ignored'
        ELSE IF (PERM(I).LE.PERM(J)) THEN
          COUNT(I) = COUNT(I) + 1
        ELSE
          COUNT(J) = COUNT(J) + 1
        END IF
   30 CONTINUE
      IPE(1) = COUNT(1) + 1
      DO 40 I = 2,N
        IPE(I) = IPE(I-1) + COUNT(I) + 1
   40 CONTINUE
      DO 50 K = 1,NE
        I = IRN(K)
        J = JCN(K)
        IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 50
        IF (PERM(I).LE.PERM(J)) THEN
          IW(IPE(I)) = J
          IPE(I) = IPE(I) - 1
        ELSE
          IW(IPE(J)) = I
          IPE(J) = IPE(J) - 1
        END IF
   50 CONTINUE
      IWFR = 1
      INFO(4) = 0
      DO 70 I = 1,N
        L = IPE(I)
        IPE(I) = IWFR
        DO 60 K = L + 1,L + COUNT(I)
          J = IW(K)
          IF (FLAG(J).NE.I) THEN
            FLAG(J) = I
            IWFR = IWFR + 1
            IW(IWFR) = J
          ELSE
            IF (I.LT.J) INFO(4) = INFO(4) + 1
          END IF
   60   CONTINUE
        IF (IWFR.GT.IPE(I)) THEN
          IW(IPE(I)) = IWFR - IPE(I)
          IWFR = IWFR + 1
        ELSE
          IPE(I) = 0
        END IF
   70 CONTINUE
      IF (INFO(4).GT.0) THEN
        INFO(1) = INFO(1) + 2
        IF (LDIAG.GT.1 .AND. WP.GE.0) WRITE (WP,'(A/I10,A)')
     +      '*** Warning message from subroutine MA57AD ***',
     +      INFO(4),' off-diagonal duplicate entries found'
      END IF
      END
      SUBROUTINE MA57KD(N, IPE, IW, LW, IWFR, PERM, IPS, NV, FLAG,
     *                  NCMPA)
      INTEGER N,LW,IWFR,NCMPA
      INTEGER IPE(N)
      INTEGER IW(LW), PERM(N), IPS(N), NV(N), FLAG(N)
      INTEGER I,J,ML,MS,ME,IP,MINJS,IE,KDUMMY,JP
      INTEGER LN,JP1,JS,LWFR,JP2,JE
      EXTERNAL MA57FD
      DO 10 I=1,N
        FLAG(I) = 0
        NV(I) = 0
        J = PERM(I)
        IPS(J) = I
   10 CONTINUE
      NCMPA = 0
      DO 100 ML=1,N
        MS = IPS(ML)
        ME = MS
        FLAG(MS) = ME
        IP = IWFR
        MINJS = N
        IE = ME
        DO 70 KDUMMY=1,N
          JP = IPE(IE)
          LN = 0
          IF (JP.LE.0) GO TO 60
          LN = IW(JP)
          DO 50 JP1=1,LN
            JP = JP + 1
            JS = IW(JP)
            IF (FLAG(JS).EQ.ME) GO TO 50
            FLAG(JS) = ME
            IF (IWFR.LT.LW) GO TO 40
            IPE(IE) = JP
            IW(JP) = LN - JP1
            CALL MA57FD(N, IPE, IW, IP-1, LWFR, NCMPA)
            JP2 = IWFR - 1
            IWFR = LWFR
            IF (IP.GT.JP2) GO TO 30
            DO 20 JP=IP,JP2
              IW(IWFR) = IW(JP)
              IWFR = IWFR + 1
   20       CONTINUE
   30       IP = LWFR
            JP = IPE(IE)
   40       IW(IWFR) = JS
            MINJS = MIN0(MINJS,PERM(JS)+0)
            IWFR = IWFR + 1
   50     CONTINUE
   60     IPE(IE) = -ME
          JE = NV(IE)
          NV(IE) = LN + 1
          IE = JE
          IF (IE.EQ.0) GO TO 80
   70   CONTINUE
   80   IF (IWFR.GT.IP) GO TO 90
        IPE(ME) = 0
        NV(ME) = 1
        GO TO 100
   90   MINJS = IPS(MINJS)
        NV(ME) = NV(MINJS)
        NV(MINJS) = ME
        IW(IWFR) = IW(IP)
        IW(IP) = IWFR - IP
        IPE(ME) = IP
        IWFR = IWFR + 1
  100 CONTINUE
      RETURN
      END
C** end of MA57KD**
      SUBROUTINE MA57FD(N, IPE, IW, LW, IWFR, NCMPA)
      INTEGER N,LW,IWFR,NCMPA
      INTEGER IPE(N)
      INTEGER   IW(LW)
      INTEGER I,K1,LWFR,IR,K,K2
      NCMPA = NCMPA + 1
      DO 10 I=1,N
        K1 = IPE(I)
        IF (K1.LE.0) GO TO 10
        IPE(I) = IW(K1)
        IW(K1) = -I
   10 CONTINUE
      IWFR = 1
      LWFR = IWFR
      DO 60 IR=1,N
        IF (LWFR.GT.LW) GO TO 70
        DO 20 K=LWFR,LW
          IF (IW(K).LT.0) GO TO 30
   20   CONTINUE
        GO TO 70
   30   I = -IW(K)
        IW(IWFR) = IPE(I)
        IPE(I) = IWFR
        K1 = K + 1
        K2 = K + IW(IWFR)
        IWFR = IWFR + 1
        IF (K1.GT.K2) GO TO 50
        DO 40 K=K1,K2
          IW(IWFR) = IW(K)
          IWFR = IWFR + 1
   40   CONTINUE
   50   LWFR = K2 + 1
   60 CONTINUE
   70 RETURN
      END
C--------------------------------------------------------------------
C-             Copyright CCLRC Rutherford Appleton Laboratory
C--------------------------------------------------------------------
      SUBROUTINE MA57LD(N, IPE, NV, IPS, NE, NA, NODE, PERM, NSTEPS,
     *                  FILS, FRERE, ND, NEMIN, SUBORD)
      INTEGER N, NSTEPS
      INTEGER ND(N)
      INTEGER IPE(N), FILS(N), FRERE(N), SUBORD(N)
      INTEGER NV(N), IPS(N), NE(N), NA(N), NODE(N), PERM(N)
      INTEGER NEMIN
      INTEGER I,IF,IS,NR,NR1,INS,INL,INB,INF,INFS,INSW
      INTEGER K,L,ISON,IN,IFSON,INO
      INTEGER INOS,IB,IL,INT
      INTEGER IPERM
      DO 10 I=1,N
        IPS(I) = 0
        NE(I) = 0
        NODE(I) = 0
        SUBORD(I) = 0
   10 CONTINUE
      NR = N + 1
      DO 50 I=1,N
        IF = -IPE(I)
        IF (NV(I).EQ.0) THEN
          IF (SUBORD(IF).NE.0) SUBORD(I) = SUBORD(IF)
          SUBORD(IF) = I
          NODE(IF) = NODE(IF)+1
        ELSE
          IF (IF.NE.0) THEN
            IS = -IPS(IF)
            IF (IS.GT.0) IPE(I) = IS
            IPS(IF) = -I
          ELSE
            NR = NR - 1
            NE(NR) = I
          ENDIF
        ENDIF
   50 CONTINUE
      DO 999 I=1,N
       FILS(I) = IPS(I)
 999  CONTINUE
      NR1 = NR
      INS = 0
 1000 IF (NR1.GT.N) GO TO 1151
      INS = NE(NR1)
      NR1 = NR1 + 1
 1070 INL = FILS(INS)
      IF (INL.LT.0) THEN
       INS = -INL
       GO TO 1070
      ENDIF
 1080 IF (IPE(INS).LT.0) THEN
       INS       = -IPE(INS)
       FILS(INS) = 0
       GO TO 1080
      ENDIF
      IF (IPE(INS).EQ.0) THEN
       INS = 0
       GO TO 1000
      ENDIF
      INB = IPE(INS)
C?? I think this test is the wrong way round
      IF (NV(INB).GE.NV(INS)) THEN
C?? So reversed
       INS = INB
       GO TO 1070
      ENDIF
      INF = INB
 1090 INF = IPE(INF)
      IF (INF.GT.0) GO TO 1090
      INF  = -INF
      INFS = -FILS(INF)
      IF (INFS.EQ.INS) THEN
        FILS(INF) = -INB
        IPS(INF)  = -INB
        IPE(INS)  = IPE(INB)
        IPE(INB)  = INS
      ELSE
        INSW = INFS
 1100   INFS = IPE(INSW)
        IF (INFS.NE.INS) THEN
          INSW = INFS
          GO TO 1100
        ENDIF
        IPE(INS) = IPE(INB)
        IPE(INB) = INS
        IPE(INSW)= INB
      ENDIF
        INS      = INB
        GO TO 1070
 1151 DO 51 I=1,N
       FRERE(I) = IPE(I)
       FILS(I) = IPS(I)
 51   CONTINUE
      IS = 1
      I = 0
      IPERM = 1
      DO 160 K=1,N
        IF (I.GT.0) GO TO 60
        IF (NR.GT.N) GO TO 161
        I = NE(NR)
        NE(NR) = 0
        NR = NR + 1
        IL = N
        NA(N) = 0
   60   CONTINUE
        DO 70 L=1,N
          IF (IPS(I).GE.0) GO TO 80
          ISON = -IPS(I)
          IPS(I) = 0
          I = ISON
          IL = IL - 1
          NA(IL) = 0
   70   CONTINUE
   80   CONTINUE
C?? Do we want to expand for subordinate variables
        IPS(I) = K
        NE(IS) = NE(IS) + NODE(I) + 1
        IF (IL.LT.N) NA(IL+1) = NA(IL+1) + 1
        NA(IS) = NA(IL)
        ND(IS) = NV(I)
        NODE(I) = IS
        PERM(I) = IPERM
        IPERM = IPERM + 1
        IN = I
  777   IF (SUBORD(IN).EQ.0) GO TO 778
          IN = SUBORD(IN)
          NODE(IN) = IS
          PERM(IN) = IPERM
          IPERM = IPERM + 1
          GO TO 777
  778   IF (NA(IS).NE.1) GO TO 90
        IF (ND(IS-1)-NE(IS-1).EQ.ND(IS)) GO TO 100
   90   IF (NE(IS).GE.NEMIN) GO TO 110
        IF (NA(IS).EQ.0) GO TO 110
        IF (NE(IS-1).GE.NEMIN) GO TO 110
  100   NA(IS-1) = NA(IS-1) + NA(IS) - 1
        ND(IS-1) = ND(IS) + NE(IS-1)
        NE(IS-1) = NE(IS) + NE(IS-1)
        NE(IS) = 0
        NODE(I) = IS-1
        IFSON = -FILS(I)
        IN = IFSON
 102    INO = IN
        IN =  FRERE(IN)
        IF (IN.GT.0) GO TO 102
        NV(INO) = 0
        IN = I
  888   IF (SUBORD(IN).EQ.0) GO TO 889
        IN = SUBORD(IN)
        NODE(IN) = IS-1
        GO TO 888
  889   SUBORD(IN) = INO
        IN = INO
        IF (SUBORD(IN).EQ.0) GO TO 887
        IN = SUBORD(IN)
        IPE(IN) = -I
  887   CONTINUE
      INOS = -FILS(INO)
      IF (IFSON.EQ.INO) GO TO 107
      IN = IFSON
 105  INS = IN
      IN =  FRERE(IN)
      IF (IN.NE.INO) GO TO 105
        IF (INOS.EQ.0) THEN
          FRERE(INS) = -I
          GO TO 120
        ELSE
          FRERE(INS) =  INOS
        ENDIF
 107    IN = INOS
        IF (IN.EQ.0) GO TO 120
 108    INT = IN
        IN =  FRERE(IN)
        IF (IN.GT.0) GO TO 108
        FRERE(INT) = -I
        GO TO 120
  110   IS = IS + 1
  120   IB = IPE(I)
        IF (IB.GE.0) THEN
          IF (IB.GT.0) NA(IL) = 0
          I = IB
          GO TO 160
        ELSE
          I = -IB
          IL = IL + 1
        ENDIF
  160 CONTINUE
  161 NSTEPS = IS - 1
      RETURN
      END
      SUBROUTINE MA57MD(N,NE,IRN,JCN,MAP,IRNPRM,
     +                  LROW,PERM,COUNT,IDIAG)
      INTEGER N,NE
      INTEGER IRN(NE),JCN(NE),MAP(NE),IRNPRM(N+NE),LROW(N),PERM(N),
     +        COUNT(N),
     +        IDIAG(N)
      INTEGER EXPNE,I,J,K
      DO 10 I = 1,N
        COUNT(I) = 1
        IDIAG(I) = 0
   10 CONTINUE
      EXPNE = NE + N
      DO 20 K = 1,NE
        I = IRN(K)
        J = JCN(K)
        IF (MAX(I,J).GT.N .OR. MIN(I,J).LT.1) THEN
          EXPNE = EXPNE - 1
          GO TO 20
        ENDIF
        IF (I.EQ.J) THEN
          I = PERM(I)
          IF (IDIAG(I).GE.1) THEN
            COUNT(I) = COUNT(I) + 1
            IDIAG(I) = IDIAG(I) + 1
          ELSE
            IDIAG(I) = 1
            EXPNE = EXPNE - 1
          ENDIF
          GO TO 20
        ENDIF
        IF (PERM(I).LT.PERM(J)) THEN
          I = PERM(I)
          COUNT(I) = COUNT(I) + 1
        ELSE
          J = PERM(J)
          COUNT(J) = COUNT(J) + 1
        END IF
   20 CONTINUE
      LROW(1) = COUNT(1)
      IDIAG(1) = MAX(IDIAG(1),1)
      DO 30 I = 2,N
        LROW(I) = COUNT(I)
        COUNT(I) = COUNT(I-1) + LROW(I)
        IDIAG(I) = COUNT(I-1) + MAX(IDIAG(I),1)
   30 CONTINUE
      DO 35 I = 1,N
        K = PERM(I)
        IRNPRM(IDIAG(K)) = I
   35 CONTINUE
      DO 40 K = 1,NE
        I = IRN(K)
        J = JCN(K)
        IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) THEN
          MAP(K) = 0
          GO TO 40
        ENDIF
        I = PERM(IRN(K))
        J = PERM(JCN(K))
        IF (I.EQ.J) THEN
          MAP(K) = IDIAG(I)
          IRNPRM(IDIAG(I)) = IRN(K)
          IDIAG(I) = IDIAG(I) - 1
        ELSE
          IF (I.GT.J) THEN
            MAP(K) = COUNT(J)
            IRNPRM(COUNT(J)) = IRN(K)
            COUNT(J) = COUNT(J) - 1
          ELSE
            MAP(K) = COUNT(I)
            IRNPRM(COUNT(I)) = JCN(K)
            COUNT(I) = COUNT(I) - 1
          ENDIF
        ENDIF
   40 CONTINUE
      IDIAG(1) = EXPNE
      RETURN
      END
      SUBROUTINE MA57ND(N,LENR,NA,NE,ND,NSTEPS,LSTKI,LSTKR,
     *                  INFO,RINFO)
      INTEGER N,NSTEPS
      INTEGER LENR(N),LSTKI(N),LSTKR(N),NA(NSTEPS),
     +        ND(NSTEPS),NE(NSTEPS),INFO(40)
      DOUBLE PRECISION RINFO(20)
      INTEGER I,IORG,ISTKI,ISTKR,ITOP,ITREE,JORG,K,
     +        LSTK,NASSR,NELIM,NFR,NSTK,NTOTPV,NZ1,NZ2
      DOUBLE PRECISION DELIM
      INTRINSIC MAX
      DOUBLE PRECISION OPS,OPSASS
      INTEGER NIRADU,NIRNEC,NIRTOT,NRLADU,NRLNEC,NRLTOT,MAXFRT
      NZ1 = 0
      DO 40 I = 1,N
        NZ1 = NZ1 + LENR(I)
   40 CONTINUE
      NZ2 = NZ1
      ISTKI = 0
      ISTKR = 0
      OPS = 0.0D0
      OPSASS = 0.0D0
      NRLADU = 0
      NIRADU = 3
      NIRTOT = NZ1+N+5
      NRLTOT = NZ1
      NIRNEC = NZ2+N+5
      NRLNEC = NZ2
      NTOTPV = 0
      ITOP = 0
      MAXFRT = 0
      DO 100 ITREE = 1,NSTEPS
        NELIM = NE(ITREE)
        DELIM = NELIM
        NFR = ND(ITREE)
        MAXFRT = MAX(MAXFRT,NFR)
        NSTK = NA(ITREE)
        NASSR = NELIM*(NELIM+1)/2 + NFR*NFR
        NRLTOT = MAX(NRLTOT,NRLADU+NASSR+ISTKR+NZ1)
        NRLNEC = MAX(NRLNEC,NRLADU+NASSR+ISTKR+NZ2)
        DO 70 IORG = 1,NELIM
          JORG = NTOTPV + IORG
          OPSASS = OPSASS + LENR(JORG)
          NZ2 = NZ2 - LENR(JORG)
   70   CONTINUE
        NTOTPV = NTOTPV + NELIM
        DO 80 K = 1,NSTK
          LSTK = LSTKR(ITOP)
          ISTKR = ISTKR - LSTK
          OPSASS = OPSASS + LSTK
          LSTK = LSTKI(ITOP)
          ISTKI = ISTKI - LSTK
          ITOP = ITOP - 1
   80   CONTINUE
        NRLADU = NRLADU + (NELIM*(NELIM+1))/2 + (NFR-NELIM)*NELIM
        NIRADU = NIRADU + 2 + NFR
        OPS = OPS + (DELIM* (12*NFR+6*NFR*NFR - (DELIM+1)*
     +        (6*NFR+6-(2*DELIM+1))))/6 + DELIM
        IF (NFR.GT.NELIM) THEN
          ITOP = ITOP + 1
          LSTKR(ITOP) = ((NFR-NELIM)*(NFR-NELIM+1))/2
          LSTKI(ITOP) = NFR - NELIM + 1
          ISTKI = ISTKI + LSTKI(ITOP)
          ISTKR = ISTKR + LSTKR(ITOP)
        ENDIF
        IF (ITREE.EQ.NSTEPS) THEN
          NIRTOT = MAX(NIRTOT,NIRADU+ISTKI+NZ1)
          NIRNEC = MAX(NIRNEC,NIRADU+ISTKI+NZ2)
        ELSE
          NIRTOT = MAX(NIRTOT,NIRADU+(N-NTOTPV+2)+ISTKI+NZ1)
          NIRNEC = MAX(NIRNEC,NIRADU+(N-NTOTPV+2)+ISTKI+NZ2)
        ENDIF
  100 CONTINUE
      INFO(5)   = NRLADU
      INFO(6)   = NIRADU
      INFO(7)   = MAXFRT
      INFO(8)   = NSTEPS
      INFO(9)   = NRLTOT
      INFO(10)  = NIRTOT
      INFO(11)  = NRLNEC
      INFO(12)  = NIRNEC
      RINFO(1)  = OPSASS
      RINFO(2)  = OPS
      RETURN
      END
      SUBROUTINE MA57OD(N,NE,A,LA,IW,LIW,LROW,PERM,NSTEPS,NSTK,NODE,
     +                  DIAG,SCHNAB,PPOS,CNTL,ICNTL,INFO,RINFO,HOLD,
     +                  BIGA)
      INTEGER N,NE,LA
      DOUBLE PRECISION A(LA),DIAG(N),SCHNAB(*),CNTL(5),RINFO(20),BIGA
      INTEGER LIW,IW(LIW),LROW(N),PERM(N),NSTEPS,NSTK(NSTEPS),
     +        NODE(N),PPOS(N),ICNTL(20),INFO(40),HOLD(40)
      INTEGER ZCOL,RPOS
      DOUBLE PRECISION ZERO,HALF,ONE
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0)
      INTEGER AINPUT
      DOUBLE PRECISION AMAX,AMULT1,AMULT2
      INTEGER APOS,APOSA,APOSB,APOSBB,APOSBK,APOSC,APOSI,APOSJ,APOSM,
     +        APOS1,APOS2,APOS3,APOS4,ASTK,ATRASH,BLK
      DOUBLE PRECISION DELTA,DETPIV
      INTEGER ELT
      DOUBLE PRECISION FLOPSA,FLOPSB,FLOPSX
      INTEGER I,I1,IASS,IBEG,IELL,IEND,IEXCH,IINPUT,INTSPA,
     +        IORG,IPIV,IPOS,IROW,ISNPIV,ISTK,ISWOP,IWNFS,IWPOS,
     +        J,JAY,JA1,JCOL,JJ,JJJ,JMAX,J1,J2,K,
     +        KB,KBLK,KCT,KR,KROW,K1,K2,L,LASPIV,LDIAG,LIELL,
     +        LP,LPIV, NBSTATIC
      LOGICAL LASTBK,LTWO
      INTEGER MAXFRT
      DOUBLE PRECISION MAXPIV
      INTEGER NASS,NBLK,NBLOC,NCMPBI,NCMPBR,NEIG,NELL,NFRONT,NIRBDU
      DOUBLE PRECISION NORMJ
      INTEGER NTWO
      LOGICAL SCHUR,LSTAT
      INTEGER MPIV,NPIV,NPOTPV,NRLBDU,NSC1,NST,
     +        NSTACK(2),NSTKAC(2),NTOTPV,
     +        NUMORG,OFFDAG,PHASE,PIVBLK
      DOUBLE PRECISION PIVOT
      INTEGER PIVSIZ,POSELT,POSPV1,POSPV2,PTRA,PTRIRN,RLSPA,
     +        SIZBLK,SIZC,SIZF,TRLSPA,TINSPA,TOTSTA(2),WP
      DOUBLE PRECISION RMAX,SWOP,TMAX,TOL,UU,ULOC,UTARG,STCTOL
      DOUBLE PRECISION FD05AD
C?? To identify bug
      INTRINSIC MIN,MAX,ABS
      EXTERNAL DGEMM,FD05AD,MA57PD,MA57WD
      NBLOC = ICNTL(11)
      TOL = CNTL(2)
      LP = ICNTL(1)
      WP = ICNTL(2)
      LDIAG = ICNTL(5)
      UU = MIN(CNTL(1),HALF)
      UU = MAX(UU,ZERO)
      LSTAT = .FALSE.
      IF (CNTL(4).GT.ZERO) THEN
        IF (CNTL(5).EQ.ZERO) LSTAT = .TRUE.
        UTARG = SQRT(UU/CNTL(4))*CNTL(4)
        STCTOL = BIGA*CNTL(4)
      ENDIF
      IF (HOLD(1).GT.0) THEN
        INFO(1) = 0
        NBLK = HOLD(2)
        NTWO = HOLD(3)
        INFO(23) = HOLD(4)
        NCMPBR = 0
        NCMPBI = 0
        NEIG   = HOLD(6)
        MAXFRT = HOLD(7)
        IWPOS  = HOLD(8)
        APOS   = HOLD(9)
        APOSBB = HOLD(10)
        NSTKAC(1) = HOLD(11)
        NSTKAC(2) = HOLD(12)
        AINPUT  = HOLD(13)
        IINPUT  = HOLD(14)
        ISTK    = HOLD(15)
        ASTK    = HOLD(16)
        INTSPA  = HOLD(17)
        RLSPA   = HOLD(18)
        PTRIRN  = HOLD(19)
        PTRA    = HOLD(20)
        NTOTPV  = HOLD(21)
        NPOTPV  = HOLD(22)
        NUMORG  = HOLD(23)
        NFRONT  = HOLD(24)
        NASS    = HOLD(25)
        IF (HOLD(1).EQ.1) NELL    = HOLD(27)
        IF (HOLD(1).EQ.2) NPIV    = HOLD(27)
        IASS    = HOLD(28)
        TINSPA  = HOLD(29)
        TRLSPA  = HOLD(30)
        TOTSTA(1) = HOLD(31)
        TOTSTA(2) = HOLD(32)
        NSTACK(1) = HOLD(33)
        NSTACK(2) = HOLD(34)
        INFO(32)  = HOLD(37)
        INFO(33)  = HOLD(38)
        INFO(34)  = HOLD(39)
        NBSTATIC  = HOLD(40)
        IF (ICNTL(7).EQ.2 .OR.ICNTL(7).EQ.3) ISNPIV = HOLD(35)
        IF (ICNTL(7).EQ.4) PHASE = HOLD(36)
        IF (HOLD(1).EQ.2) NSC1    = NFRONT-NPIV
        FLOPSA = RINFO(3)
        FLOPSB = RINFO(4)
        FLOPSX = RINFO(5)
        IF (HOLD(1).EQ.1) THEN
          HOLD(1) = 0
          GO TO 333
        ELSE
          HOLD(1) = 0
          GO TO 444
        ENDIF
      ENDIF
      NBSTATIC = 0
      NBLK = 0
      NTWO = 0
      NCMPBR = 0
      NCMPBI = 0
      FLOPSA = ZERO
      FLOPSB = ZERO
      FLOPSX = ZERO
      NEIG = 0
      MAXFRT  = 0
      INFO(1) = 0
      INFO(2) = 0
      INFO(17) = 0
      INFO(40) = 0
      INFO(26) = 0
      INFO(27) = 0
      INFO(32) = 0
      INFO(33) = 0
      INFO(34) = 0
      INFO(23) = 0
      RINFO(3) = ZERO
      RINFO(4) = ZERO
      RINFO(5) = ZERO
      RINFO(15) = ZERO
      DO 10 I = 1,N
        PPOS(I) = N + 1
   10 CONTINUE
      IWPOS = 6
      IW(1) = 0
      IW(2) = 0
      IW(3) = 0
      IW(4) = 0
      IW(5) = 0
      APOSBB = 1
      NSTACK(1) = 0
      NSTACK(2) = 0
      NSTKAC(1) = NE
      NSTKAC(2) = NE
      TOTSTA(1) = NE
      TOTSTA(2) = NE
      INTSPA = NE+5+N
      RLSPA = NE
      TINSPA = NE+5+N
      TRLSPA = NE
      PTRIRN = LIW - NE + 1
      PTRA = LA - NE + 1
      ISTK = PTRIRN - 1
      ASTK = PTRA - 1
      AINPUT = PTRA
      IINPUT = PTRIRN
      NTOTPV = 0
      NPOTPV = 0
      IF (ICNTL(7).EQ.2 .OR. ICNTL(7).EQ.3) ISNPIV = 0
      IF (ICNTL(7).EQ.4) THEN
        PHASE = 1
        DO 19 I = 1,N
          DIAG(I) = ZERO
   19   CONTINUE
        APOS1 = PTRA-1
        J1 = PTRIRN
        DO 20 I = 1,N
          J2 = J1 + LROW(I) - 1
          DO 25 JJ = J1,J2
            J = IW(JJ)
            APOS1 = APOS1 + 1
            IF (J.EQ.PERM(I)) DIAG(J) = DIAG(J) + A(APOS1)
   25     CONTINUE
          J1 = J2 + 1
   20   CONTINUE
        SCHNAB(1) = ZERO
        SCHNAB(5) = ZERO
        DO 21 I = 1,N
          SCHNAB(1) = MAX(SCHNAB(1),ABS(DIAG(I)))
          SCHNAB(5) = MIN(SCHNAB(5),DIAG(I))
   21   CONTINUE
        SCHNAB(4) = SCHNAB(1)
        SCHNAB(2) = FD05AD(1)**(1.0/3.0)
        SCHNAB(3) = 0.1
        RINFO(15) = FD05AD(5)
        DELTA     = ZERO
      ENDIF
      IASS = 1
 2160 CONTINUE
        NUMORG = 0
        DO 30 I = NPOTPV + 1,N
          J = PERM(I)
          IF (ABS(NODE(J)).GT.IASS) GO TO 40
          IW(IWPOS+NUMORG) = J
          NUMORG = NUMORG + 1
          PPOS(J) = NUMORG
   30   CONTINUE
   40   NASS = NUMORG
        NELL = NSTK(IASS)
        IELL = ISTK + 1
        DO 70 ELT = 1,NELL
          DO 50 JJ = IELL + 1,IELL + IW(IELL)
            J = IW(JJ)
            IF (NODE(J).GT.IASS) GO TO 50
            IF (PPOS(J).LE.N) GO TO 50
            IW(IWPOS+NASS) = J
            NASS = NASS + 1
            PPOS(J) = NASS
   50     CONTINUE
          IELL = IELL + IW(IELL) + 1
   70   CONTINUE
        IWNFS = IWPOS + NASS
        J1 = PTRIRN
        DO 90 IORG = 1,NUMORG
          J2 = J1 + LROW(NPOTPV+IORG) - 1
          DO 80 JJ = J1,J2
            J = IW(JJ)
            IF (PPOS(J).LE.N) GO TO 80
            IW(IWNFS) = J
            IWNFS = IWNFS + 1
            PPOS(J) = IWNFS - IWPOS
   80     CONTINUE
          J1 = J2 + 1
   90   CONTINUE
        IELL = ISTK + 1
        DO 170 ELT = 1,NELL
          J1 = IELL+1
          J2 = IELL+IW(IELL)
          DO 150 JJ = J1,J2
            J = IW(JJ)
            IF (PPOS(J).LE.N) GO TO 150
            IW(IWNFS) = J
            IWNFS = IWNFS + 1
            PPOS(J) = IWNFS - IWPOS
  150     CONTINUE
          IELL = J2 + 1
  170   CONTINUE
        NFRONT = IWNFS - IWPOS
        MAXFRT = MAX(MAXFRT,NFRONT)
        IF (INFO(1).NE.-3) THEN
          APOS = APOSBB + (NASS*(NASS+1))/2
        ELSE
          APOS = 1
        END IF
        RLSPA  = MAX(RLSPA,INFO(40)+APOS+NFRONT*NFRONT-1+NSTKAC(1))
        TRLSPA = MAX(TRLSPA,INFO(40)+APOS+NFRONT*NFRONT-1+TOTSTA(1))
  333   IF (APOS+NFRONT*NFRONT-1.GT.ASTK) THEN
          CALL MA57PD(A,IW,ASTK,AINPUT,PTRA,.TRUE.)
          NCMPBR = NCMPBR + 1
          IF (APOS+NFRONT*NFRONT-1.GT.ASTK) THEN
            IF (ICNTL(8).NE.0) THEN
              HOLD(1) = 1
              HOLD(2) = NBLK
              HOLD(3) = NTWO
              HOLD(4) = INFO(23)
              HOLD(5) = NCMPBI
              HOLD(6) = NEIG
              HOLD(7) = MAXFRT
              HOLD(8) = IWPOS
              HOLD(9) = APOS
              HOLD(10) = APOSBB
              HOLD(11) = NSTKAC(1)
              HOLD(12) = NSTKAC(2)
              HOLD(13) = AINPUT
              HOLD(14) = IINPUT
              HOLD(15) = ISTK
              HOLD(16) = ASTK
              HOLD(17) = INTSPA
              HOLD(18) = RLSPA
              HOLD(19) = PTRIRN
              HOLD(20) = PTRA
              HOLD(21) = NTOTPV
              HOLD(22) = NPOTPV
              HOLD(23) = NUMORG
              HOLD(24) = NFRONT
              HOLD(25) = NASS
              HOLD(27) = NELL
              HOLD(28) = IASS
              HOLD(29) = TINSPA
              HOLD(30) = TRLSPA
              HOLD(31) = TOTSTA(1)
              HOLD(32) = TOTSTA(2)
              HOLD(33) = NSTACK(1)
              HOLD(34) = NSTACK(2)
              IF (ICNTL(7).EQ.2 .OR.ICNTL(7).EQ.3) HOLD(35) = ISNPIV
              IF (ICNTL(7).EQ.4) HOLD(36) = PHASE
              HOLD(37) = INFO(32)
              HOLD(38) = INFO(33)
              HOLD(39) = INFO(34)
              RINFO(3) =FLOPSA
              RINFO(4) =FLOPSB
              RINFO(5) =FLOPSX
              HOLD(40) = NBSTATIC
              INFO(35) = HOLD(40)
              INFO(1) = 10
              RETURN
            ELSE
              INFO(40) = INFO(40) + APOS - 1
              APOS = 1
              APOSBB = 1
              INFO(1) = -3
              IF (NFRONT*NFRONT.GT.ASTK) THEN
                INFO(17) = MAX(INFO(17),RLSPA)
                IF (ICNTL(7).EQ.4) INFO(17) = MAX(INFO(17),RLSPA + N)
                INFO(2) = LA
                RETURN
              ENDIF
            ENDIF
          ENDIF
        END IF
        ATRASH = APOS + NFRONT*NFRONT - 1
        DO 210 JJ = APOS,ATRASH
          A(JJ) = ZERO
  210   CONTINUE
        J1 = PTRIRN
        DO 230 IORG = 1,NUMORG
          J = PERM(NPOTPV+IORG)
          APOSI = APOS + (PPOS(J)-1)*NFRONT - 1
          J2 = J1 + LROW(NPOTPV+IORG) - 1
          FLOPSA = FLOPSA + J2 - J1 + 1
          DO 220 JJ = J1,J2
            JAY = IW(JJ)
            IF (PPOS(JAY).GE.PPOS(J)) THEN
              APOS2 = APOSI + PPOS(JAY)
            ELSE
              APOS2 = APOS + (PPOS(JAY)-1)*NFRONT + PPOS(J) - 1
            ENDIF
            A(APOS2) = A(APOS2) + A(PTRA)
            PTRA = PTRA + 1
  220     CONTINUE
          NSTKAC(1) = NSTKAC(1) - J2 + J1 - 1
          J1 = J2 + 1
  230   CONTINUE
        NSTKAC(2) = NSTKAC(2) - J1 + PTRIRN
        PTRIRN = J1
        NPOTPV = NPOTPV + NUMORG
C???
C?? Depends if we need lower triangle and whether all entries are already
        DO 380 ELT = 1,NELL
          POSELT = ASTK + 1
          LIELL = IW(ISTK+1)
          J1 = ISTK + 2
          J2 = ISTK+1 + LIELL
          FLOPSA = FLOPSA + (LIELL*(LIELL+1))/2
          DO 250 JJ = J1,J2
            J = IW(JJ)
            APOS2 = APOS + (PPOS(J)-1)*NFRONT
            APOS1 = POSELT
            DO 240 JJJ=JJ,J2
              JAY = IW(JJJ)
C???          APOS3 = APOS2 + PPOS(JAY) - 1
C???          A(APOS3) = A(APOS3) + A(APOS1)
C???          APOS5 = APOS+(PPOS(JAY)-1)*NFRONT+PPOS(J)-1
C???          IF (APOS3.NE.APOS5) A(APOS5) = A(APOS5) + A(APOS1)
              IF (PPOS(JAY) .GE. PPOS(J)) THEN
                APOS3 = APOS2 + PPOS(JAY) - 1
              ELSE
                APOS3 = APOS+(PPOS(JAY)-1)*NFRONT+PPOS(J)-1
              ENDIF
              A(APOS3) = A(APOS3) + A(APOS1)
              APOS1 = APOS1 + 1
  240       CONTINUE
            POSELT = POSELT + LIELL - (JJ-J1)
  250     CONTINUE
          NSTKAC(2) = NSTKAC(2) - (J2-ISTK)
          NSTACK(2) = NSTACK(2) - (J2-ISTK)
          TOTSTA(2) = TOTSTA(2) - (J2-ISTK)
          ISTK = J2
          ASTK = ASTK + (LIELL*(LIELL+1))/2
          NSTKAC(1) = NSTKAC(1) - (LIELL*(LIELL+1))/2
          NSTACK(1) = NSTACK(1) - (LIELL*(LIELL+1))/2
          TOTSTA(1) = TOTSTA(1) - (LIELL*(LIELL+1))/2
  380   CONTINUE
C1122     CONTINUE
        PIVBLK = MIN(NBLOC,NASS)
        APOSBK = APOS
        NPIV = 0
        ULOC = UU
        DO 918 BLK = 1,NASS
        IF (NPIV+PIVBLK .GE. NASS) THEN
          LASTBK = .TRUE.
          SIZBLK = NASS - NPIV
        ELSE
          LASTBK = .FALSE.
          SIZBLK = PIVBLK
        ENDIF
        LASPIV = NPIV
        MPIV = 0
CCCCCCCC
        KR = 0
CCC Set to following to force 2 by 2 pivots in Nocedal examples
        KCT = SIZBLK + 1
  920   CONTINUE
          KR = KR + 1
          KCT = KCT - 1
          IF (KCT.EQ.0) GO TO 930
          IF (KR.GT.SIZBLK) KR = MPIV + 1
          IPIV = LASPIV + KR
            APOSI = APOS + (IPIV-1)*NFRONT
            POSPV1 = APOSI + IPIV - 1
            PIVOT = A(POSPV1)
   29       IF (ICNTL(7).EQ.4 .AND. PHASE.EQ.2) THEN
              IF (INFO(27).EQ.0) INFO(27) = NTOTPV + 1
              NORMJ = ZERO
              DO 28 I = POSPV1+1,POSPV1+NFRONT-NPIV-1
                NORMJ = NORMJ + ABS(A(I))
   28         CONTINUE
              DELTA = MAX(ZERO,
     *                    - A(POSPV1) + MAX(NORMJ,SCHNAB(2)*SCHNAB(1)))
              A(POSPV1) = A(POSPV1) + DELTA
              IF (A(POSPV1).EQ.ZERO) GO TO 970
              RINFO(15) = MIN(RINFO(15),A(POSPV1))
              DIAG(PERM(NTOTPV+1)) = DELTA
              PIVSIZ = 1
              GO TO 811
            ENDIF
            IF (ICNTL(7).GT.1) THEN
              IF (ABS(PIVOT).LE.CNTL(2)) THEN
                IF (ICNTL(7).LT.4) GO TO 970
                PHASE = 2
                GO TO 29
              ENDIF
              IF (NTOTPV.EQ.0) THEN
                IF (PIVOT.GT.ZERO) ISNPIV = 1
                IF (PIVOT.LT.ZERO) ISNPIV = -1
              ELSE
                IF (ICNTL(7).EQ.2 .AND. ISNPIV*PIVOT.LT.ZERO) GO TO 980
                IF (ICNTL(7).EQ.3 .AND. ISNPIV*PIVOT.LT.ZERO) THEN
                    INFO(26) = INFO(26) + 1
                    ISNPIV = -ISNPIV
                ENDIF
              ENDIF
              IF (ICNTL(7).EQ.4) THEN
                IF (PIVOT.GE.SCHNAB(1)*SCHNAB(2) .AND.
     *              SCHNAB(5).GE.-SCHNAB(3)*SCHNAB(4)) THEN
                  SCHNAB(5) = ZERO
                  SCHNAB(4) = ZERO
                  DO 22 I = POSPV1+1,POSPV1+NFRONT-NPIV-1
                    J = IW(IWPOS+NPIV+I-POSPV1)
                    DIAG(J) = DIAG(J) - A(I)*A(I)/PIVOT
                    SCHNAB(5) = MIN(DIAG(J),SCHNAB(5))
                    SCHNAB(4) = MAX(DIAG(J),SCHNAB(4))
                    IF (DIAG(J).LT.-SCHNAB(3)*SCHNAB(1)) THEN
                      PHASE = 2
                      GO TO 29
                    ENDIF
   22             CONTINUE
                  DIAG(PERM(NTOTPV+1)) = ZERO
                  RINFO(15) = MIN(RINFO(15),PIVOT)
                ELSE
                  PHASE = 2
                  GO TO 29
                ENDIF
              ENDIF
              PIVSIZ = 1
              GO TO 811
            ENDIF
            AMAX = ZERO
            JMAX = 0
            DO 110 K = 1, IPIV - NPIV - 1
              IF (ABS(A(POSPV1-K*NFRONT)).GT.AMAX) THEN
                AMAX = ABS(A(POSPV1-K*NFRONT))
                JMAX = IPIV - K
              ENDIF
  110       CONTINUE
            DO 111 K =  1, MIN(NASS,LASPIV+PIVBLK) - IPIV
              IF (ABS(A(POSPV1+K)).GT.AMAX) THEN
                AMAX = ABS(A(POSPV1+K))
                JMAX = IPIV + K
              ENDIF
  111       CONTINUE
            RMAX = ZERO
            DO 112 K = MIN(NASS,LASPIV+PIVBLK)-IPIV+1,NFRONT-IPIV
               RMAX = MAX(RMAX,ABS(A(POSPV1+K)))
 112        CONTINUE
            IF (MAX(AMAX,RMAX,ABS(PIVOT)).LE.TOL) THEN
              GO TO 920
            END IF
            IF (MAX(AMAX,ABS(PIVOT)).LE.TOL) GO TO 920
            PIVSIZ = 0
            IF (ABS(PIVOT).GT.ULOC*MAX(RMAX,AMAX)) THEN
              PIVSIZ = 1
              A(POSPV1) = PIVOT
              GO TO 810
            END IF
            IF (NPIV+1.EQ.NASS) THEN
              A(POSPV1) = PIVOT
              GO TO 920
            END IF
            IF (AMAX.LE.TOL) GO TO 920
            IF (RMAX.LT.AMAX) THEN
              RMAX = ZERO
              DO 113 K = 1, IPIV - NPIV - 1
                IF (IPIV-K.EQ.JMAX) GO TO 113
                RMAX=MAX(RMAX,ABS(A(POSPV1-K*NFRONT)))
  113         CONTINUE
              DO 114 K =  1, NFRONT - IPIV
                IF (IPIV+K.EQ.JMAX) GO TO 114
                RMAX = MAX(RMAX,ABS(A(POSPV1+K)))
  114         CONTINUE
            ENDIF
            APOSJ = APOS + (JMAX-1)*NFRONT
            POSPV2 = APOSJ + JMAX - 1
            IF (IPIV.GT.JMAX) THEN
              OFFDAG = APOSJ + IPIV - 1
            ELSE
              OFFDAG = APOSI + JMAX - 1
            END IF
            TMAX = ZERO
            DO 115 K = 1, JMAX - NPIV - 1
              IF (JMAX-K.EQ.IPIV) GO TO 115
              TMAX=MAX(TMAX,ABS(A(POSPV2-K*NFRONT)))
  115       CONTINUE
            DO 116 K =  1, NFRONT - JMAX
              IF (JMAX+K.EQ.IPIV) GO TO 116
              TMAX = MAX(TMAX,ABS(A(POSPV2+K)))
  116       CONTINUE
            DETPIV = A(POSPV1)*A(POSPV2) - AMAX*AMAX
            MAXPIV = MAX(ABS(A(POSPV1)),ABS(A(POSPV2)))
            IF (MAXPIV.EQ.ZERO) MAXPIV = ONE
            IF (ABS(DETPIV)/MAXPIV.LE.TOL) GO TO 920
            PIVSIZ = 2
            IF ((ABS(A(POSPV2))*RMAX+AMAX*TMAX)*ULOC.GT.
     +          ABS(DETPIV)) GO TO 920
            IF ((ABS(A(POSPV1))*TMAX+AMAX*RMAX)*ULOC.GT.
     +          ABS(DETPIV)) GO TO 920
  810       LPIV = IPIV
            IF (PIVSIZ.EQ.2) LPIV = MIN(IPIV,JMAX)
            KR = MAX(KR,NPIV+PIVSIZ)
            KCT = SIZBLK - MPIV - PIVSIZ + 1
            DO 860 KROW = NPIV,NPIV + PIVSIZ - 1
              IF (LPIV.EQ.KROW+1) GO TO 850
              JA1 = APOS + (LPIV-1)
              J1 = APOS + KROW
              DO 820 JJ = 1,KROW
                SWOP = A(JA1)
                A(JA1) = A(J1)
                A(J1) = SWOP
                JA1 = JA1 + NFRONT
                J1 = J1 + NFRONT
  820         CONTINUE
              JA1 = JA1 + NFRONT
              J1 = J1 + 1
              DO 830 JJ = 1,LPIV - KROW - 2
                SWOP = A(JA1)
                A(JA1) = A(J1)
                A(J1) = SWOP
                JA1 = JA1 + NFRONT
                J1 = J1 + 1
  830         CONTINUE
              SWOP = A(APOS+KROW* (NFRONT+1))
              A(APOS+KROW* (NFRONT+1)) = A(JA1)
              A(JA1) = SWOP
              DO 840 JJ = 1,NFRONT - LPIV
                JA1 = JA1 + 1
                J1 = J1 + 1
                SWOP = A(JA1)
                A(JA1) = A(J1)
                A(J1) = SWOP
  840         CONTINUE
              IPOS = IWPOS + KROW
              IEXCH = IWPOS + LPIV - 1
              ISWOP = IW(IPOS)
              IW(IPOS) = IW(IEXCH)
              IW(IEXCH) = ISWOP
  850         LPIV = MAX(IPIV,JMAX)
  860       CONTINUE
  811       POSPV1 = APOS + NPIV* (NFRONT+1)
            POSPV2 = POSPV1 + NFRONT + 1
            IF (PIVSIZ.EQ.1) THEN
              FLOPSB = FLOPSB + ONE
              A(POSPV1) = ONE/A(POSPV1)
              IF (A(POSPV1).LT.ZERO) NEIG = NEIG + 1
              J1 = POSPV1 + 1
              J2 = POSPV1 + NASS - (NPIV+1)
              IBEG = POSPV1 + NFRONT + 1
              IEND = APOS + (NPIV+1)*NFRONT + NFRONT - 1
              DO 880 JJ = J1,J2
                AMULT1 = -A(JJ)*A(POSPV1)
                IF (.NOT.LASTBK) A(POSPV1+(JJ-J1+1)*NFRONT) = A(JJ)
                JCOL = JJ
                FLOPSB = FLOPSB + (IEND-IBEG+1)*2 + 1
                IF (MPIV+JJ-J1+2.GT.PIVBLK) GO TO 871
CDIR$            IVDEP
                DO 870 IROW = IBEG,IEND
                  A(IROW) = A(IROW) + AMULT1*A(JCOL)
                  JCOL = JCOL + 1
  870           CONTINUE
  871           A(JJ) = AMULT1
                IBEG = IBEG + NFRONT + 1
                IEND = IEND + NFRONT
  880         CONTINUE
              NPIV = NPIV + 1
              MPIV = MPIV + 1
              NTOTPV = NTOTPV + 1
              IF (MPIV.EQ.SIZBLK) GO TO 930
            ELSE
              OFFDAG = POSPV1 + 1
              FLOPSB = FLOPSB + 6.0
              SWOP = A(POSPV2)
              IF (DETPIV.LT.ZERO) THEN
                NEIG = NEIG + 1
              ELSE
                IF (SWOP.LT.ZERO) NEIG = NEIG + 2
              END IF
              A(POSPV2) = A(POSPV1)/DETPIV
              A(POSPV1) = SWOP/DETPIV
              A(OFFDAG) = -A(OFFDAG)/DETPIV
              J1 = POSPV1 + 2
              J2 = POSPV1 + NASS - (NPIV+1)
              IBEG = POSPV2 + NFRONT + 1
              IEND = APOS + (NPIV+2)*NFRONT + NFRONT - 1
              DO 900 JJ = J1,J2
                K1 = JJ
                K2 = JJ + NFRONT
                AMULT1 = - (A(POSPV1)*A(K1)+A(POSPV1+1)*A(K2))
                AMULT2 = - (A(POSPV1+1)*A(K1)+A(POSPV2)*A(K2))
                IF (.NOT.LASTBK) THEN
                  A(POSPV1 + (JJ-J1+2)*NFRONT) = A(K1)
                  A(POSPV1 + (JJ-J1+2)*NFRONT + 1) = A(K2)
                ENDIF
                FLOPSB = FLOPSB + (IEND-IBEG+1)*4 + 6
                IF (MPIV+JJ-J1+3.GT.PIVBLK) GO TO 891
CDIR$            IVDEP
                DO 890 IROW = IBEG,IEND
                  A(IROW) = A(IROW) + AMULT1*A(K1) + AMULT2*A(K2)
                  K1 = K1 + 1
                  K2 = K2 + 1
  890           CONTINUE
  891           A(JJ) = AMULT1
                A(JJ+NFRONT) = AMULT2
                IBEG = IBEG + NFRONT + 1
                IEND = IEND + NFRONT
  900         CONTINUE
              IPOS = IWPOS + NPIV
              IW(IPOS) = -IW(IPOS)
              IW(IPOS+1) = -IW(IPOS+1)
              NPIV = NPIV + 2
              MPIV = MPIV + 2
              NTOTPV = NTOTPV + 2
              NTWO = NTWO + 1
              IF (MPIV.EQ.SIZBLK) GO TO 930
            END IF
        GO TO 920
 930    IF (LASTBK) THEN
          IF (NPIV.EQ.NASS) GO TO 935
          IF (.NOT. LSTAT)  GO TO 935
          ULOC = ULOC/10.0D0
          IF (ULOC.LT.UTARG) THEN
            ULOC = ULOC * 10.0D0
            GO TO 9919
          ENDIF
          KCT = SIZBLK + 1 - MPIV
          GO TO 920
        ENDIF
        IF (MPIV.EQ.0) THEN
          PIVBLK = 2*PIVBLK
          GO TO 918
        ENDIF
        KBLK = (NASS-(LASPIV+PIVBLK))/PIVBLK
        L = NASS - (LASPIV+PIVBLK)
        APOS4 = APOS+(LASPIV+PIVBLK)*(NFRONT+1)
        DO 931 KB = 1,KBLK
          FLOPSX = FLOPSX + PIVBLK*(PIVBLK-1)*MPIV
          CALL DGEMM('N','N',L-(KB-1)*PIVBLK,PIVBLK,MPIV,ONE,
     +               A(APOSBK+PIVBLK*KB),NFRONT,
     +               A(APOSBK+PIVBLK*KB*NFRONT),NFRONT,ONE,
     +               A(APOS4+PIVBLK*(KB-1)*(NFRONT+1)),NFRONT)
          IF (NFRONT.GT.NASS)
     +    CALL DGEMM('N','T',NFRONT-NASS,PIVBLK,MPIV,ONE,
     +               A(APOSBK+NASS-LASPIV),NFRONT,
     +               A(APOSBK+PIVBLK*KB),NFRONT,ONE,
     +               A(APOSBK+KB*NFRONT*PIVBLK+NASS-LASPIV),NFRONT)
  931   CONTINUE
       SIZC = NASS - (KBLK+1)*PIVBLK - LASPIV
       SIZF = NFRONT - (KBLK+1)*PIVBLK - LASPIV
       APOSA = APOSBK + (KBLK+1)*PIVBLK
       DO 934 K = 1,MPIV
         APOSB = APOSBK + NFRONT*PIVBLK*(KBLK+1) + K - 1
         APOSM = APOSBK + PIVBLK*(KBLK+1) + (K-1)*NFRONT
         APOSC = APOSBK + PIVBLK*(KBLK+1)*(NFRONT+1)
         DO 933 JJ = 1,SIZC
            DO 932 J = JJ,SIZC
              A(APOSC+J-1) = A(APOSC+J-1) + A(APOSA+J-1)*A(APOSB)
  932       CONTINUE
            DO 936 J = SIZC+1,SIZF
              A(APOSC+J-1) = A(APOSC+J-1) + A(APOSA+J-1)*A(APOSM)
  936       CONTINUE
            APOSC = APOSC + NFRONT
            APOSB = APOSB + NFRONT
            APOSM = APOSM + 1
  933     CONTINUE
          APOSA = APOSA + NFRONT
  934   CONTINUE
        APOSBK = APOSBK + MPIV*(NFRONT+1)
        LASPIV = NPIV
  918   CONTINUE
        IF (LP.GE.0) WRITE(LP,'(A)') '****** BE WORRIED LOOP 918'
 9919      IPIV = LASPIV+MPIV
 9920      IPIV = IPIV + 1
CADD Probably not needed .. use only IPIV
           APOSI = APOS + (IPIV-1)*NFRONT
           POSPV1 = APOSI + IPIV - 1
           PIVOT = A(POSPV1)
CADD
           PIVSIZ = 1
           LPIV = IPIV
           AMAX = ZERO
           DO 9876 K = 1, IPIV - NPIV - 1
             AMAX = MAX(AMAX,ABS(A(POSPV1-K*NFRONT)))
 9876      CONTINUE
           DO 9878 K =  1, NFRONT - IPIV
             AMAX = MAX(AMAX,ABS(A(POSPV1+K)))
 9878      CONTINUE
           IF (ABS(A(POSPV1)).LT.STCTOL) THEN
               PIVOT = STCTOL
              IF (A(POSPV1) .LT. ZERO) THEN
                 A(POSPV1) = -PIVOT
                 PIVOT     = -PIVOT
              ELSE
                 A(POSPV1) = PIVOT
              ENDIF
              NBSTATIC = NBSTATIC + 1
           ENDIF
           FLOPSB = FLOPSB + ONE
           A(POSPV1) = ONE/A(POSPV1)
           IF (A(POSPV1).LT.ZERO) NEIG = NEIG + 1
           J1 = POSPV1 + 1
           J2 = POSPV1 + NASS - (NPIV+1)
           IBEG = POSPV1 + NFRONT + 1
CBUG
           IEND = APOSI + 2*NFRONT - 1
           DO 9880 JJ = J1,J2
              AMULT1 = -A(JJ)*A(POSPV1)
CADD  Not needed since LASTBK always true
              JCOL = JJ
              FLOPSB = FLOPSB + (IEND-IBEG+1)*2 + 1
CADD Not necessary because control is on setting of J2
CDIR$            IVDEP
              DO 9870 IROW = IBEG,IEND
                 A(IROW) = A(IROW) + AMULT1*A(JCOL)
                 JCOL = JCOL + 1
 9870         CONTINUE
 9871         A(JJ) = AMULT1
              IBEG = IBEG + NFRONT + 1
              IEND = IEND + NFRONT
 9880      CONTINUE
           NPIV = NPIV + 1
           MPIV = MPIV + 1
           NTOTPV = NTOTPV + 1
           IF (MPIV.LT.SIZBLK) GO TO 9920
CADD
  935   SCHUR = (NBLOC.LT.(NFRONT-NASS) .AND. NPIV.GE.NBLOC)
        NSC1 = NFRONT - NPIV
        IF (IASS.NE.NSTEPS) INFO(23) = INFO(23) + NASS - NPIV
        IF (CNTL(4).GT.ZERO .AND. INFO(23).GT.CNTL(5)*N) LSTAT = .TRUE.
        IF (NSC1.EQ.0) GO TO 1830
        IF (.NOT.SCHUR) THEN
          RLSPA = MAX(RLSPA,INFO(40)+APOS+NFRONT*NFRONT-1+
     +                      NSTKAC(1))
          TRLSPA = MAX(TRLSPA,INFO(40)+APOS+NFRONT*NFRONT-1+
     +                      TOTSTA(1))
          NSTKAC(1) = NSTKAC(1) + ((NSC1+1)*NSC1)/2
          NSTACK(1) = NSTACK(1) + ((NSC1+1)*NSC1)/2
          TOTSTA(1) = TOTSTA(1) + ((NSC1+1)*NSC1)/2
          APOSI = APOS + NFRONT*NFRONT - 1
          DO 1370 JJ = 1,NFRONT-NPIV
            J = APOSI
            DO 1360 JJJ = 1,JJ
                A(ASTK) = A(J)
                ASTK = ASTK - 1
                J = J - 1
 1360       CONTINUE
            APOSI = APOSI - NFRONT
 1370     CONTINUE
          APOS4 = ASTK + 1
          J1 = IWPOS
          LTWO = .FALSE.
          POSPV1 = APOS
          DO 1450 I1 = 1,NPIV
            IF (LTWO) GO TO 1440
            APOSI = APOS + (I1-1)*NFRONT + NASS
            J2 = APOS + NFRONT* (I1-1) + NFRONT - 1
CCC What happens here ??
            APOSC = APOS4 + ((NASS-NPIV)*(2*NFRONT-NPIV-NASS+1))/2
            IF (IW(J1).GT.0) THEN
              FLOPSB = FLOPSB + (NFRONT-NASS) +
     *                          (NFRONT-NASS)* (NFRONT-NASS+1)
              DO 1410 JJ = APOSI,J2
                AMULT1 = -A(JJ)*A(POSPV1)
                DO 1400 JJJ = JJ,J2
                  A(APOSC) = A(APOSC) + AMULT1*A(JJJ)
                  APOSC = APOSC + 1
 1400           CONTINUE
                A(JJ) = AMULT1
 1410         CONTINUE
              J1 = J1 + 1
            ELSE
              POSPV2 = POSPV1 + NFRONT + 1
              OFFDAG = POSPV1 + 1
              FLOPSB = FLOPSB + 6* (NFRONT-NASS) +
     +                 2* (NFRONT-NASS)* (NFRONT-NASS+1)
              DO 1430 JJ = APOSI,J2
                AMULT1 = - (A(POSPV1)*A(JJ)+A(OFFDAG)*A(JJ+NFRONT))
                AMULT2 = -A(POSPV2)*A(JJ+NFRONT) - A(OFFDAG)*A(JJ)
                DO 1420 JJJ = JJ,J2
                  A(APOSC) = A(APOSC) + AMULT1*A(JJJ) +
     +                       AMULT2*A(JJJ+NFRONT)
                  APOSC = APOSC + 1
 1420           CONTINUE
                A(JJ) = AMULT1
                A(JJ+NFRONT) = AMULT2
 1430         CONTINUE
              J1 = J1 + 2
              POSPV1 = POSPV2
              LTWO = .TRUE.
              GO TO 1450
            END IF
 1440       LTWO = .FALSE.
            POSPV1 = POSPV1 + NFRONT + 1
 1450     CONTINUE
        ELSE
          APOS4 = APOS+NASS*(NFRONT+1)
        APOS3 = APOS+NASS*NFRONT
        J1 = IWPOS
        LTWO = .FALSE.
        POSPV1 = APOS
          DO 1490 I = 1,NPIV
            IF (LTWO) GO TO 1480
            APOSI = APOS + (I-1)*NFRONT + NASS
            POSELT = APOS3 + I - 1
            IF (IW(J1).GT.0) THEN
              FLOPSB = FLOPSB + (NFRONT-NASS)
              DO 1460 JJ = APOSI,APOS + NFRONT*I - 1
                A(POSELT) = A(JJ)
                A(JJ) = -A(JJ)*A(POSPV1)
                POSELT = POSELT + NFRONT
 1460         CONTINUE
              J1 = J1 + 1
            ELSE
              POSPV2 = POSPV1 + NFRONT + 1
              OFFDAG = POSPV1 + 1
              FLOPSB = FLOPSB + 6* (NFRONT-NASS)
              DO 1470 JJ = APOSI,APOS + NFRONT*I - 1
                A(POSELT) = A(JJ)
                A(POSELT+1) = A(JJ+NFRONT)
                A(JJ) = - (A(POSPV1)*A(JJ)+A(OFFDAG)*A(JJ+NFRONT))
                A(JJ+NFRONT) = -A(POSPV2)*A(JJ+NFRONT) -
     +                         A(OFFDAG)*A(POSELT)
                POSELT = POSELT + NFRONT
 1470         CONTINUE
              J1 = J1 + 2
              POSPV1 = POSPV2
              LTWO = .TRUE.
              GO TO 1490
            END IF
 1480       LTWO = .FALSE.
            POSPV1 = POSPV1 + NFRONT + 1
 1490     CONTINUE
          FLOPSB = FLOPSB + NPIV* (NFRONT-NASS)**2 +
     *                      NPIV* (NFRONT-NASS)
          KBLK = ( NFRONT-NASS)/NBLOC
          L =  NFRONT - NASS
          DO 1500 KB = 1,KBLK
            FLOPSX = FLOPSX + NBLOC* (NBLOC-1)* (NPIV)
            CALL DGEMM('N','N',L-(KB-1)*NBLOC,NBLOC,NPIV,ONE,
     +                 A(APOS+NASS+NBLOC*(KB-1)),NFRONT,
     +                 A(APOS3+NBLOC*(KB-1)*NFRONT),NFRONT,ONE,
     +                 A(APOS4+NBLOC*(NFRONT+1)*(KB-1)),NFRONT)
 1500     CONTINUE
          DO 1550 I = 1 + KBLK*NBLOC,L
            APOSA = APOS + NASS
            APOSB = APOS3 +(I-1)*NFRONT
            APOSC = APOS4 + (I-1)*NFRONT - 1
            DO 1540 K = 1,NPIV
              DO 1530 J = I,L
                A(APOSC+J) = A(APOSC+J) + A(APOSA+J-1)*A(APOSB)
 1530         CONTINUE
              APOSA = APOSA + NFRONT
              APOSB = APOSB + 1
 1540       CONTINUE
 1550     CONTINUE
          JA1 = APOS+NFRONT*NFRONT-1
          NSTKAC(1) = NSTKAC(1) + ((NSC1+1)* (NSC1))/2
          NSTACK(1) = NSTACK(1) + ((NSC1+1)* (NSC1))/2
          TOTSTA(1) = TOTSTA(1) + ((NSC1+1)* (NSC1))/2
          DO 1710 I = NSC1,1,-1
            DO 1700 JJ = JA1,JA1-(NSC1-I),-1
              A(ASTK) = A(JJ)
              ASTK = ASTK - 1
 1700       CONTINUE
            JA1 = JA1 - NFRONT
 1710     CONTINUE
        END IF
        NSTKAC(2) = NSTKAC(2) + NSC1 + 1
        NSTACK(2) = NSTACK(2) + NSC1 + 1
        TOTSTA(2) = TOTSTA(2) + NSC1 + 1
 1830   IF (IASS.EQ.NSTEPS) THEN
          INTSPA = MAX(INTSPA,IWPOS+NFRONT-1+NSTKAC(2))
          TINSPA = MAX(TINSPA,IWPOS+NFRONT-1+TOTSTA(2))
          GO TO 2158
        ELSE
          INTSPA = MAX(INTSPA,IWPOS+NFRONT-1+(N-NTOTPV+2)+NSTKAC(2))
          TINSPA = MAX(TINSPA,IWPOS+NFRONT-1+(N-NTOTPV+2)+TOTSTA(2))
        ENDIF
  444   NST = 0
        IF (NSC1.GT.0) NST = NSC1 + 1
        IF (IWPOS+NFRONT-1+(N-NTOTPV+2)+NST.GT.ISTK) THEN
          CALL MA57PD(A,IW,ISTK,IINPUT,PTRIRN,.FALSE.)
          NCMPBI = NCMPBI + 1
          IF (IWPOS+NFRONT-1+(N-NTOTPV+2)+NST.GT.ISTK) THEN
            IF (ICNTL(8).NE.0) THEN
              HOLD(1) = 2
              HOLD(2) = NBLK
              HOLD(3) = NTWO
              HOLD(4) = INFO(23)
              HOLD(5) = NCMPBI
              HOLD(6) = NEIG
              HOLD(7) = MAXFRT
              HOLD(8) = IWPOS
              HOLD(9) = APOS
              HOLD(10) = APOSBB
              HOLD(11) = NSTKAC(1)
              HOLD(12) = NSTKAC(2)
              HOLD(13) = AINPUT
              HOLD(14) = IINPUT
              HOLD(15) = ISTK
              HOLD(16) = ASTK
              HOLD(17) = INTSPA
              HOLD(18) = RLSPA
              HOLD(19) = PTRIRN
              HOLD(20) = PTRA
              HOLD(21) = NTOTPV
              HOLD(22) = NPOTPV
              HOLD(23) = NUMORG
              HOLD(24) = NFRONT
              HOLD(25) = NASS
              HOLD(27) = NPIV
              HOLD(28) = IASS
              HOLD(29) = TINSPA
              HOLD(30) = TRLSPA
              HOLD(31) = TOTSTA(1)
              HOLD(32) = TOTSTA(2)
              HOLD(33) = NSTACK(1)
              HOLD(34) = NSTACK(2)
              IF (ICNTL(7).EQ.2 .OR.ICNTL(7).EQ.3) HOLD(35) = ISNPIV
              IF (ICNTL(7).EQ.4) HOLD(36) = PHASE
              HOLD(37) = INFO(32)
              HOLD(38) = INFO(33)
              HOLD(39) = INFO(34)
              NSC1    = NFRONT-NPIV
              RINFO(3) =FLOPSA
              RINFO(4) =FLOPSB
              RINFO(5) =FLOPSX
              INFO(1) = 11
              HOLD(40) = NBSTATIC
              INFO(35) = HOLD(40)
            ELSE
              INFO(1)  = -4
              INFO(2)  = LIW
              INFO(18) = INTSPA
            ENDIF
            RETURN
          END IF
        END IF
        IF (NSC1.GT.0) THEN
          DO 1720 I = 1,NSC1
            IW(ISTK) = IW(IWPOS+NFRONT-I)
            ISTK = ISTK - 1
 1720     CONTINUE
          IW(ISTK) = NSC1
          ISTK = ISTK - 1
        ENDIF
        DO 1840 JJ = IWPOS + NPIV,IWPOS + NFRONT - 1
          J = ABS(IW(JJ))
          PPOS(J) = N + 1
 1840   CONTINUE
C********************************
C********************************
 2158   IF (NPIV.EQ.0) GO TO 2159
        NBLK = NBLK + 1
        IW(IWPOS-2) = NFRONT
        IW(IWPOS-1) = NPIV
        IWPOS = IWPOS + NFRONT + 2
        IF (INFO(1).EQ.-3) THEN
          INFO(40) = INFO(40) + (NPIV * (2*NFRONT-NPIV+1))/2
          GO TO 2159
        END IF
        APOS2 = APOSBB
        DO 2130 I = 1,NPIV
          JA1 = APOS + (I-1)* (NFRONT+1)
          DO 2120 J = I,NPIV
            A(APOS2) = A(JA1)
            IF (A(APOS2).EQ.ZERO) INFO(32) = INFO(32) + 1
            APOS2 = APOS2 + 1
            JA1 = JA1 + 1
 2120     CONTINUE
 2130   CONTINUE
        RPOS = APOS2
        DO 2150 I = 1,NPIV
          JA1 = APOS + (I-1)*NFRONT + NPIV
          DO 2140 J = 1,NFRONT - NPIV
            A(APOS2) = A(JA1)
            APOS2 = APOS2 + 1
            JA1 = JA1 + 1
 2140     CONTINUE
 2150   CONTINUE
        APOSBB = APOS2
        DO 2152 J = 1,NFRONT-NPIV
        APOS2 = RPOS+J-1
        ZCOL = 1
          DO 2151 I = 1,NPIV
            IF (A(APOS2).EQ.ZERO) INFO(33) = INFO(33)+1
            IF (A(APOS2).NE.ZERO) ZCOL = 0
            APOS2 = APOS2 + NFRONT - NPIV
 2151     CONTINUE
        IF (ZCOL.EQ.1) INFO(34) = INFO(34)+1
 2152   CONTINUE
 2159   IASS = IASS + 1
      IF (IASS.LE.NSTEPS) THEN
        IW(IWPOS-2) = 0
        IW(IWPOS-1) = 0
        GO TO 2160
      ENDIF
C2160 CONTINUE
      INFO(35) = NBSTATIC
      IF (INFO(1).EQ.-3) THEN
        INFO(2)  = LA
        INFO(17) = MAX(INFO(17),RLSPA)
        IF (ICNTL(7).EQ.4) INFO(17) = MAX(INFO(17),RLSPA + N)
        RETURN
      END IF
      GO TO 1000
 970  INFO(1) = -5
      INFO(2) = NTOTPV + 1
      IF (LDIAG.GT.0 .AND. LP.GE.0)
     *    WRITE(LP,99992) INFO(1),PIVOT,CNTL(2),INFO(2),ICNTL(7)
99992 FORMAT (/'*** Error message from routine MA57BD **',
     *       '   INFO(1) = ',I3/'Pivot has value ',D16.8,' when ',
     *       'CNTL(2) has value ',D16.8/
     *       'at stage',I11,2X,'when ICNTL(7) =',I3)
      RETURN
 980  INFO(1) = -6
      INFO(2) = NTOTPV + 1
      IF (LDIAG.GT.0 .AND. LP.GE.0)
     *    WRITE(LP,99993) INFO(1),INFO(2),ICNTL(7)
99993 FORMAT (/'*** Error message from routine MA57BD **',
     *       '   INFO(1) = ',I3/'Change in sign of pivot at stage',
     *       I10,2X,'when ICNTL(7) = ',I3)
      RETURN
 1000 NRLBDU = APOSBB - 1
      NIRBDU = IWPOS - 3
      IF (NTOTPV.NE.N) THEN
        INFO(1) = 4
        IF (LDIAG.GT.0 .AND. WP.GE.0)
     *      WRITE(WP,99994) INFO(1),NTOTPV
99994 FORMAT (/'*** Warning message from routine MA57BD **',
     *         '   INFO(1) =',I2/5X, 'Matrix is singular, rank =', I5)
        DO 3331 I = 1,N
          PPOS(I) = 0
 3331   CONTINUE
        IWPOS = 4
        DO 3332 I = 1,NBLK
          NFRONT = IW(IWPOS)
          NPIV = IW(IWPOS+1)
          DO 3330 J = IWPOS+2,IWPOS+NPIV+1
            PPOS(ABS(IW(J))) = 1
 3330     CONTINUE
          IWPOS = IWPOS + NFRONT + 2
 3332   CONTINUE
        K= 0
        DO 3333 I=1,N
          IF (PPOS(I).EQ.0) THEN
            K=K+1
            NBLK = NBLK + 1
            NRLBDU = NRLBDU+1
            A(NRLBDU) = ONE
            IW(NIRBDU+1) = 1
            IW(NIRBDU+2) = 1
            IW(NIRBDU+3) = I
            NIRBDU = NIRBDU+3
          ENDIF
 3333   CONTINUE
      ENDIF
      INFO(14) = NRLBDU
      IW(1) = NRLBDU + 1
      IW(2) = NRLBDU + NTWO
      INFO(15) = IW(2)
      IW(3) = NBLK
      INFO(31) = NBLK
      CALL MA57WD(A,LA,IW,LIW,NRLBDU)
      INFO(16) = NIRBDU
      INFO(18) = INTSPA
      INFO(20) = TINSPA
      INFO(17) = RLSPA
      INFO(19) = TRLSPA
      INFO(21) = MAXFRT
      INFO(22) = NTWO
      INFO(24) = NEIG
      INFO(25) = NTOTPV
      INFO(28) = NCMPBR
      INFO(29) = NCMPBI
      RINFO(3) = FLOPSA
      RINFO(4) = FLOPSB
      RINFO(5) = FLOPSX
      IF (INFO(27).GT.0) THEN
        RINFO(14) = ZERO
        DO 332 I = 1,N
          RINFO(14) = MAX(RINFO(14),DIAG(I))
 332    CONTINUE
      ENDIF
      RETURN
      END
      SUBROUTINE MA57PD(A,IW,J1,J2,ITOP,REAL)
      INTEGER ITOP,J1,J2
      LOGICAL REAL
      DOUBLE PRECISION A(*)
      INTEGER IW(*)
      INTEGER IPOS,JJ
      IF (J2.EQ.ITOP) GO TO 50
      IPOS = ITOP - 1
      IF (REAL) THEN
        DO 10 JJ = J2-1,J1+1,-1
          A(IPOS) = A(JJ)
          IPOS = IPOS - 1
   10   CONTINUE
      ELSE
        DO 20 JJ = J2-1,J1+1,-1
          IW(IPOS) = IW(JJ)
          IPOS = IPOS - 1
   20   CONTINUE
      ENDIF
      J2 = ITOP
      J1 = IPOS
   50 RETURN
      END
      SUBROUTINE MA57WD(A,LA,IW,LIW,NRLBDU)
      INTEGER LA,LIW
      DOUBLE PRECISION A(LA)
      INTEGER IW(LIW)
      INTEGER NRLBDU
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER APOS,IBLK,IROW,IWPOS,J,JPIV,NCOLS,NROWS
      APOS = 1
      IWPOS = 6
      DO 40 IBLK = 1,IW(3)
        NCOLS = IW(IWPOS-2)
        NROWS = IW(IWPOS-1)
        JPIV = 1
        DO 30 IROW = 1,NROWS
          JPIV = JPIV - 1
          IF (JPIV.EQ.1) GO TO 10
          IF (IW(IWPOS+IROW-1).LT.0) THEN
            JPIV = 2
            NRLBDU = NRLBDU + 1
            A(NRLBDU) = A(APOS+1)
            A(APOS+1) = ZERO
          END IF
   10     DO 20 J = APOS + 1,APOS + NROWS - IROW
            A(J) = -A(J)
   20     CONTINUE
          APOS = APOS + NROWS - IROW + 1
   30   CONTINUE
        APOS = APOS + NROWS* (NCOLS-NROWS)
        IWPOS = IWPOS + NCOLS + 2
   40 CONTINUE
      END
      SUBROUTINE MA57XD(N,FACT,LFACT,IFACT,LIFACT,RHS,LRHS,
     *                  W,LW,IW1,ICNTL)
      INTEGER N,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),LRHS,LW
      DOUBLE PRECISION W(LW),RHS(LRHS)
      INTEGER IW1(N),ICNTL(20)
      INTRINSIC ABS
      EXTERNAL DGEMV,DTPSV
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      INTEGER APOS,I,IBLK,II,IPIV,IRHS,IWPOS,J,J1,J2,K,K1,K2,
     +        NCOLS,NROWS
      DOUBLE PRECISION W1,W2
      APOS = 1
      IWPOS = 4
      DO 270 IBLK = 1,IFACT(3)
        IW1(IBLK) = IWPOS
        NCOLS = IFACT(IWPOS)
        NROWS = IFACT(IWPOS+1)
        IWPOS = IWPOS + 2
        IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN
          DO 10 I = 1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            W(I) = RHS(II)
   10     CONTINUE
          CALL DTPSV('L','N','U',NROWS,FACT(APOS),W,1)
          APOS = APOS + (NROWS* (NROWS+1))/2
          IF (NCOLS.GT.NROWS) CALL DGEMV('N',NCOLS-NROWS,NROWS,
     +                                  ONE,FACT(APOS),NCOLS-NROWS,
     +                                  W,1,ONE,W(NROWS+1),1)
          APOS = APOS + NROWS* (NCOLS-NROWS)
          DO 35 I = 1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            RHS(II) = W(I)
   35     CONTINUE
        ELSE
        J1 = IWPOS
        J2 = IWPOS + NROWS - 1
        DO 130 IPIV = 1,NROWS
          APOS = APOS + 1
          W1 = RHS(ABS(IFACT(J1)))
          K = APOS
          DO 100 J = J1+1,J2
            IRHS = ABS(IFACT(J))
            RHS(IRHS) = RHS(IRHS) - FACT(K)*W1
            K = K + 1
  100     CONTINUE
          APOS = K
          J1 = J1 + 1
  130   CONTINUE
        J2 = IWPOS + NCOLS - 1
        DO 136 IPIV = 1,NROWS-1,2
          K1 = APOS
          K2 = APOS+NCOLS-NROWS
          W1 = RHS(ABS(IFACT(IWPOS+IPIV-1)))
          W2 = RHS(ABS(IFACT(IWPOS+IPIV)))
          DO 133 J = J1,J2
            IRHS = ABS(IFACT(J))
            RHS(IRHS) = RHS(IRHS) + W1*FACT(K1) + W2*FACT(K2)
            K1 = K1 + 1
            K2 = K2 + 1
  133     CONTINUE
          APOS = K2
  136   CONTINUE
        IF (MOD(NROWS,2).EQ.1) THEN
          K = APOS
          W1 = RHS(ABS(IFACT(IWPOS+IPIV-1)))
          DO 137 J = J1,J2
            IRHS = ABS(IFACT(J))
            RHS(IRHS) = RHS(IRHS) + W1*FACT(K)
            K = K + 1
  137     CONTINUE
          APOS = K
        ENDIF
      END IF
      IWPOS = IWPOS + NCOLS
  270 CONTINUE
      END
      SUBROUTINE MA57YD(N,FACT,LFACT,IFACT,LIFACT,RHS,LRHS,
     *                  W,LW,IW1,ICNTL)
      INTEGER N,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),LRHS,LW
      DOUBLE PRECISION W(LW),RHS(LRHS)
      INTEGER IW1(N),ICNTL(20)
      INTRINSIC ABS
      EXTERNAL DGEMV,DTPSV
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      INTEGER APOS,APOS2,I,IBLK,II,IPIV,IRHS,IRHS1,
     +        IRHS2,IWPOS,J,JPIV,J1,J2,K,K2,LROW,NCOLS,NROWS
      DOUBLE PRECISION W1,W2
      APOS = IFACT(1)
      APOS2 = IFACT(2)
      DO 380 IBLK = IFACT(3),1,-1
        IWPOS = IW1(IBLK)
        NCOLS = ABS(IFACT(IWPOS))
        NROWS = ABS(IFACT(IWPOS+1))
        APOS = APOS - NROWS* (NCOLS-NROWS)
        IWPOS = IWPOS + 2
        IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN
          DO 5 I = NROWS + 1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            W(I) = RHS(II)
    5     CONTINUE
          DO 10 IPIV = NROWS,1,-1
            IRHS = ABS(IFACT(IWPOS+IPIV-1))
            APOS = APOS - (NROWS+1-IPIV)
            W(IPIV) = RHS(IRHS)*FACT(APOS)
   10     CONTINUE
          JPIV = -1
          DO 20 IPIV = NROWS,1,-1
            IRHS = IFACT(IWPOS+IPIV-1)
            IF (IRHS.LT.0) THEN
              IRHS1 = -IFACT(IWPOS+IPIV-1+JPIV)
              W(IPIV) = RHS(IRHS1)*FACT(APOS2) + W(IPIV)
              IF (JPIV.EQ.1) APOS2 = APOS2 - 1
              JPIV = -JPIV
            END IF
   20     CONTINUE
          K = NCOLS - NROWS
          IF (K.GT.0) CALL DGEMV('T',K,NROWS,ONE,
     +                           FACT(APOS+(NROWS*(NROWS+1))/2),K,
     +                           W(NROWS+1),1,ONE,W,1)
          CALL DTPSV('L','T','U',NROWS,FACT(APOS),W,1)
          DO 60 I = 1,NROWS
            II = ABS(IFACT(IWPOS+I-1))
            RHS(II) = W(I)
   60     CONTINUE
        ELSE
          J1 = IWPOS
          J2 = IWPOS + NCOLS - 1
          JPIV = -1
          DO 210 IPIV = NROWS,1,-1
            IRHS = IFACT(IWPOS+IPIV-1)
            LROW = NROWS + 1 - IPIV
            IF (IRHS.GT.0) THEN
              APOS = APOS - LROW
              RHS(IRHS) = RHS(IRHS)*FACT(APOS)
            ELSE
              IF (JPIV.EQ.-1) THEN
                IRHS1 = -IFACT(IWPOS+IPIV-2)
                IRHS2 = -IRHS
                APOS = APOS - LROW - LROW - 1
                W1 = RHS(IRHS1)*FACT(APOS) +
     +               RHS(IRHS2)*FACT(APOS2)
                RHS(IRHS2) = RHS(IRHS1)*FACT(APOS2) +
     +                         RHS(IRHS2)*FACT(APOS+LROW+1)
                RHS(IRHS1) = W1
                APOS2 = APOS2 - 1
              END IF
              JPIV = -JPIV
            END IF
  210     CONTINUE
          APOS = APOS + (NROWS* (NROWS+1))/2
          K = APOS
          J1 = IWPOS + NROWS
          DO 220 IPIV = 1,NROWS-1,2
            IRHS = ABS(IFACT(IWPOS+IPIV-1))
            W1 = RHS(IRHS)
            IRHS1 = ABS(IFACT(IWPOS+IPIV))
            W2 = RHS(IRHS1)
            K2 = K+(NCOLS-NROWS)
            DO 215 J = J1,J2
              II = ABS(IFACT(J))
              W1 = W1 + FACT(K)*RHS(II)
              W2 = W2 + FACT(K2)*RHS(II)
              K = K + 1
              K2 = K2 + 1
  215       CONTINUE
            RHS(IRHS) = W1
            RHS(IRHS1) = W2
            K = K2
  220     CONTINUE
          IF (MOD(NROWS,2).EQ.1) THEN
            IRHS = ABS(IFACT(IWPOS+IPIV-1))
            W1 = RHS(IRHS)
            DO 216 J = J1,J2
              W1 = W1 + FACT(K)*RHS(ABS(IFACT(J)))
              K = K + 1
  216       CONTINUE
            RHS(IRHS) = W1
          ENDIF
          J2 = IWPOS + NROWS - 1
          DO 260 IPIV = 1,NROWS
            IRHS = ABS(IFACT(J1-1))
            APOS = APOS - IPIV
            W1 = RHS(IRHS)
            K = APOS + 1
            DO 230 J = J1,J2
              W1 = W1 - FACT(K)*RHS(ABS(IFACT(J)))
              K = K + 1
  230       CONTINUE
            RHS(IRHS) = W1
            J1 = J1 - 1
  260     CONTINUE
        END IF
  380 CONTINUE
      END
      SUBROUTINE MA57VD(N,NZ,IRN,ICN,IW,LW,IPE,IQ,FLAG,IWFR,
     +                 ICNTL,INFO)
      INTEGER IWFR,LW,N,NZ
      INTEGER FLAG(N),ICN(*),IPE(N),IQ(N),IRN(*),IW(LW)
      INTEGER ICNTL(*),INFO(*)
      INTEGER I,ID,J,JN,K,K1,K2,L,LAST,LR,N1,NDUP
      INFO(2) = 0
      DO 10 I = 1,N
        IPE(I) = 0
   10 CONTINUE
      LR = NZ
      IF (NZ.EQ.0) GO TO 120
      DO 110 K = 1,NZ
        I = IRN(K)
        J = ICN(K)
        IF (I.LT.J) THEN
          IF (I.GE.1 .AND. J.LE.N) GO TO 90
        ELSE IF (I.GT.J) THEN
          IF (J.GE.1 .AND. I.LE.N) GO TO 90
        ELSE
          IF (I.GE.1 .AND. I.LE.N) GO TO 80
        END IF
        INFO(2) = INFO(2) + 1
        INFO(1) = 1
        IF (INFO(2).LE.1 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=60) INFO(1)
        END IF
   60   FORMAT (' *** WARNING MESSAGE FROM SUBROUTINE MA57AD',
     +          '  *** INFO(1) =',I2)
        IF (INFO(2).LE.10 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=70) K,I,J
        END IF
   70   FORMAT (I6,'TH NON-ZERO (IN ROW',I6,' AND COLUMN',I6,
     +         ') IGNORED')
   80   I = 0
        J = 0
        GO TO 100
   90   IPE(I) = IPE(I) + 1
        IPE(J) = IPE(J) + 1
  100   IW(K) = J
        LR = LR + 1
        IW(LR) = I
  110 CONTINUE
  120 IQ(1) = 1
      N1 = N - 1
      IF (N1.LE.0) GO TO 140
      DO 130 I = 1,N1
        FLAG(I) = 0
        IF (IPE(I).EQ.0) IPE(I) = -1
        IQ(I+1) = IPE(I) + IQ(I) + 1
        IPE(I) = IQ(I)
  130 CONTINUE
  140 LAST = IPE(N) + IQ(N)
      FLAG(N) = 0
      IF (LR.GE.LAST) GO TO 160
      K1 = LR + 1
      DO 150 K = K1,LAST
        IW(K) = 0
  150 CONTINUE
  160 IPE(N) = IQ(N)
      IWFR = LAST + 1
      IF (NZ.EQ.0) GO TO 230
      DO 220 K = 1,NZ
        J = IW(K)
        IF (J.LE.0) GO TO 220
        L = K
        IW(K) = 0
        DO 210 ID = 1,NZ
          IF (L.GT.NZ) GO TO 170
          L = L + NZ
          GO TO 180
  170     L = L - NZ
  180     I = IW(L)
          IW(L) = 0
          IF (I.LT.J) GO TO 190
          L = IQ(J) + 1
          IQ(J) = L
          JN = IW(L)
          IW(L) = -I
          GO TO 200
  190     L = IQ(I) + 1
          IQ(I) = L
          JN = IW(L)
          IW(L) = -J
  200     J = JN
          IF (J.LE.0) GO TO 220
  210   CONTINUE
  220 CONTINUE
  230 NDUP = 0
      DO 280 I = 1,N
        K1 = IPE(I) + 1
        K2 = IQ(I)
        IF (K1.LE.K2) GO TO 240
        IPE(I) = 0
        IQ(I) = 0
        GO TO 280
  240   DO 260 K = K1,K2
          J = -IW(K)
          IF (J.LE.0) GO TO 270
          L = IQ(J) + 1
          IQ(J) = L
          IW(L) = I
          IW(K) = J
          IF (FLAG(J).NE.I) GO TO 250
          NDUP = NDUP + 1
          IW(L) = 0
          IW(K) = 0
  250     FLAG(J) = I
  260   CONTINUE
  270   IQ(I) = IQ(I) - IPE(I)
        IF (NDUP.EQ.0) IW(K1-1) = IQ(I)
  280 CONTINUE
      IF (NDUP.EQ.0) GO TO 310
      IWFR = 1
      DO 300 I = 1,N
        K1 = IPE(I) + 1
        IF (K1.EQ.1) GO TO 300
        K2 = IQ(I) + IPE(I)
        L = IWFR
        IPE(I) = IWFR
        IWFR = IWFR + 1
        DO 290 K = K1,K2
          IF (IW(K).EQ.0) GO TO 290
          IW(IWFR) = IW(K)
          IWFR = IWFR + 1
  290   CONTINUE
        IW(L) = IWFR - L - 1
  300 CONTINUE
  310 RETURN
      END
      SUBROUTINE MA57HD(N,IPE,IW,LW,IWFR,NV,NXT,LST,IPD,FLAG,IOVFLO,
     +                 NCMPA,FRATIO)
      DOUBLE PRECISION FRATIO
      INTEGER IWFR,LW,N,IOVFLO,NCMPA
      INTEGER FLAG(N),IPD(N),IPE(N),IW(LW),LST(N),NV(N),NXT(N)
      INTEGER I,ID,IDL,IDN,IE,IP,IS,JP,JP1,JP2,JS,K,K1,K2,KE,KP,KP0,KP1,
     +        KP2,KS,L,LEN,LIMIT,LN,LS,LWFR,MD,ME,ML,MS,NEL,NFLG,NP,
     +        NP0,NS,NVPIV,NVROOT,ROOT
      EXTERNAL MA57ZD
      INTRINSIC ABS,MIN
      DO 10 I = 1,N
        IPD(I) = 0
        NV(I) = 1
        FLAG(I) = IOVFLO
   10 CONTINUE
      MD = 1
      NCMPA = 0
      NFLG = IOVFLO
      NEL = 0
      ROOT = N+1
      NVROOT = 0
      DO 30 IS = 1,N
        K = IPE(IS)
        IF (K.LE.0) GO TO 20
        ID = IW(K) + 1
        NS = IPD(ID)
        IF (NS.GT.0) LST(NS) = IS
        NXT(IS) = NS
        IPD(ID) = IS
        LST(IS) = -ID
        GO TO 30
   20   NEL = NEL + 1
        FLAG(IS) = -1
        NXT(IS) = 0
        LST(IS) = 0
   30 CONTINUE
      DO 340 ML = 1,N
        IF (NEL+NVROOT+1.GE.N) GO TO 350
        DO 40 ID = MD,N
          MS = IPD(ID)
          IF (MS.GT.0) GO TO 50
   40   CONTINUE
   50   MD = ID
        NVPIV = NV(MS)
        NS = NXT(MS)
        NXT(MS) = 0
        LST(MS) = 0
        IF (NS.GT.0) LST(NS) = -ID
        IPD(ID) = NS
        ME = MS
        NEL = NEL + NVPIV
        IDN = 0
        KP = IPE(ME)
        FLAG(MS) = -1
        IP = IWFR
        LEN = IW(KP)
        DO 140 KP1 = 1,LEN
          KP = KP + 1
          KE = IW(KP)
          IF (FLAG(KE).LE.-2) GO TO 60
          IF (FLAG(KE).LE.0) THEN
             IF (IPE(KE).NE.-ROOT) GO TO 140
             KE = ROOT
             IF (FLAG(KE).LE.0) GO TO 140
          END IF
          JP = KP - 1
          LN = LEN - KP1 + 1
          IE = MS
          GO TO 70
   60     IE = KE
          JP = IPE(IE)
          LN = IW(JP)
   70     DO 130 JP1 = 1,LN
            JP = JP + 1
            IS = IW(JP)
            IF (FLAG(IS).LE.0) THEN
               IF (IPE(IS).EQ.-ROOT) THEN
                  IS = ROOT
                  IW(JP) = ROOT
                  IF (FLAG(IS).LE.0) GO TO 130
               ELSE
                  GO TO 130
               END IF
            END IF
            FLAG(IS) = 0
            IF (IWFR.LT.LW) GO TO 100
            IPE(MS) = KP
            IW(KP) = LEN - KP1
            IPE(IE) = JP
            IW(JP) = LN - JP1
            CALL MA57ZD(N,IPE,IW,IP-1,LWFR,NCMPA)
            JP2 = IWFR - 1
            IWFR = LWFR
            IF (IP.GT.JP2) GO TO 90
            DO 80 JP = IP,JP2
              IW(IWFR) = IW(JP)
              IWFR = IWFR + 1
   80       CONTINUE
   90       IP = LWFR
            JP = IPE(IE)
            KP = IPE(ME)
  100       IW(IWFR) = IS
            IDN = IDN + NV(IS)
            IWFR = IWFR + 1
            LS = LST(IS)
            LST(IS) = 0
            NS = NXT(IS)
            NXT(IS) = 0
            IF (NS.GT.0) LST(NS) = LS
            IF (LS.LT.0) THEN
              LS = -LS
              IPD(LS) = NS
            ELSE IF (LS.GT.0) THEN
              NXT(LS) = NS
            END IF
  130     CONTINUE
          IF (IE.EQ.MS) GO TO 150
          IPE(IE) = -ME
          FLAG(IE) = -1
  140   CONTINUE
  150   NV(MS) = IDN + NVPIV
        IF (IWFR.EQ.IP) GO TO 330
        K1 = IP
        K2 = IWFR - 1
        LIMIT = NINT(FRATIO*(N-NEL))
        DO 310 K = K1,K2
          IS = IW(K)
          IF (IS.EQ.ROOT) GO TO 310
          IF (NFLG.GT.2) GO TO 170
          DO 160 I = 1,N
            IF (FLAG(I).GT.0) FLAG(I) = IOVFLO
            IF (FLAG(I).LE.-2) FLAG(I) = -IOVFLO
  160     CONTINUE
          NFLG = IOVFLO
  170     NFLG = NFLG - 1
          ID = IDN
          KP1 = IPE(IS) + 1
          NP = KP1
          KP2 = IW(KP1-1) + KP1 - 1
          DO 220 KP = KP1,KP2
            KE = IW(KP)
          IF (FLAG(KE).EQ.-1) THEN
             IF (IPE(KE).NE.-ROOT) GO TO 220
             KE = ROOT
             IW(KP) = ROOT
             IF (FLAG(KE).EQ.-1) GO TO 220
          END IF
          IF (FLAG(KE).GE.0) GO TO 230
            JP1 = IPE(KE) + 1
            JP2 = IW(JP1-1) + JP1 - 1
            IDL = ID
            DO 190 JP = JP1,JP2
              JS = IW(JP)
              IF (FLAG(JS).LE.NFLG) GO TO 190
              ID = ID + NV(JS)
              FLAG(JS) = NFLG
  190       CONTINUE
            IF (ID.GT.IDL) GO TO 210
            DO 200 JP = JP1,JP2
              JS = IW(JP)
              IF (FLAG(JS).NE.0) GO TO 210
  200       CONTINUE
            IPE(KE) = -ME
            FLAG(KE) = -1
            GO TO 220
  210       IW(NP) = KE
            FLAG(KE) = -NFLG
            NP = NP + 1
  220     CONTINUE
          NP0 = NP
          GO TO 250
  230     KP0 = KP
          NP0 = NP
          DO 240 KP = KP0,KP2
            KS = IW(KP)
            IF (FLAG(KS).LE.NFLG) THEN
               IF (IPE(KS).EQ.-ROOT) THEN
                  KS = ROOT
                  IW(KP) = ROOT
                  IF (FLAG(KS).LE.NFLG) GO TO 240
               ELSE
                  GO TO 240
               END IF
            END IF
            ID = ID + NV(KS)
            FLAG(KS) = NFLG
            IW(NP) = KS
            NP = NP + 1
  240     CONTINUE
  250     IF (ID.GE.LIMIT) GO TO 295
          IW(NP) = IW(NP0)
          IW(NP0) = IW(KP1)
          IW(KP1) = ME
          IW(KP1-1) = NP - KP1 + 1
          JS = IPD(ID)
          DO 280 L = 1,N
            IF (JS.LE.0) GO TO 300
            KP1 = IPE(JS) + 1
            IF (IW(KP1).NE.ME) GO TO 300
            KP2 = KP1 - 1 + IW(KP1-1)
            DO 260 KP = KP1,KP2
              IE = IW(KP)
              IF (ABS(FLAG(IE)+0).GT.NFLG) GO TO 270
  260       CONTINUE
            GO TO 290
  270       JS = NXT(JS)
  280     CONTINUE
  290     IPE(JS) = -IS
          NV(IS) = NV(IS) + NV(JS)
          NV(JS) = 0
          FLAG(JS) = -1
          NS = NXT(JS)
          LS = LST(JS)
          IF (NS.GT.0) LST(NS) = IS
          IF (LS.GT.0) NXT(LS) = IS
          LST(IS) = LS
          NXT(IS) = NS
          LST(JS) = 0
          NXT(JS) = 0
          IF (IPD(ID).EQ.JS) IPD(ID) = IS
          GO TO 310
  295     IF (NVROOT.EQ.0) THEN
            ROOT = IS
            IPE(IS) = 0
          ELSE
            IW(K) = ROOT
            IPE(IS) = -ROOT
            NV(ROOT) = NV(ROOT) + NV(IS)
            NV(IS) = 0
            FLAG(IS) = -1
          END IF
          NVROOT = NV(ROOT)
          GO TO 310
  300     NS = IPD(ID)
          IF (NS.GT.0) LST(NS) = IS
          NXT(IS) = NS
          IPD(ID) = IS
          LST(IS) = -ID
          MD = MIN(MD,ID)
  310   CONTINUE
        DO 320 K = K1,K2
          IS = IW(K)
          IF (NV(IS).EQ.0) GO TO 320
          FLAG(IS) = NFLG
          IW(IP) = IS
          IP = IP + 1
  320   CONTINUE
        IWFR = K1
        FLAG(ME) = -NFLG
        IW(IP) = IW(K1)
        IW(K1) = IP - K1
        IPE(ME) = K1
        IWFR = IP + 1
        GO TO 335
  330   IPE(ME) = 0
  335   CONTINUE
  340 CONTINUE
  350 DO 360 IS = 1,N
        IF(NXT(IS).NE.0 .OR. LST(IS).NE.0) THEN
          IF (NVROOT.EQ.0) THEN
            ROOT = IS
            IPE(IS) = 0
          ELSE
            IPE(IS) = -ROOT
          END IF
          NVROOT = NVROOT + NV(IS)
          NV(IS) = 0
         END IF
  360 CONTINUE
      DO 370 IE = 1,N
        IF (IPE(IE).GT.0) IPE(IE) = -ROOT
  370 CONTINUE
      IF(NVROOT.GT.0)NV(ROOT)=NVROOT
      END
      SUBROUTINE MA57ZD(N,IPE,IW,LW,IWFR,NCMPA)
      INTEGER IWFR,LW,N,NCMPA
      INTEGER IPE(N),IW(LW)
      INTEGER I,IR,K,K1,K2,LWFR
      NCMPA = NCMPA + 1
      DO 10 I = 1,N
        K1 = IPE(I)
        IF (K1.LE.0) GO TO 10
        IPE(I) = IW(K1)
        IW(K1) = -I
   10 CONTINUE
      IWFR = 1
      LWFR = IWFR
      DO 60 IR = 1,N
        IF (LWFR.GT.LW) GO TO 70
        DO 20 K = LWFR,LW
          IF (IW(K).LT.0) GO TO 30
   20   CONTINUE
        GO TO 70
   30   I = -IW(K)
        IW(IWFR) = IPE(I)
        IPE(I) = IWFR
        K1 = K + 1
        K2 = K + IW(IWFR)
        IWFR = IWFR + 1
        IF (K1.GT.K2) GO TO 50
        DO 40 K = K1,K2
          IW(IWFR) = IW(K)
          IWFR = IWFR + 1
   40   CONTINUE
   50   LWFR = K2 + 1
   60 CONTINUE
   70 RETURN
      END
Cdense    Version 1
CAMD      SUBROUTINE AMDD (N, IWLEN, PE, PFREE, LEN, IW, NV, ELEN,
CAMD     $                   LAST, NCMPA, DEGREE, HEAD, NEXT, W)
      SUBROUTINE MC50BD (THRESH, N, IWLEN, PE, PFREE, LEN, IW, NV,
     $                   ELEN, LAST, NCMPA, DEGREE, HEAD, NEXT, W)
Cdense
      INTEGER N, IWLEN, PE(N), PFREE, LEN(N), IW(IWLEN), NV(N),
     $        ELEN(N), LAST(N), NCMPA, DEGREE(N), HEAD(N), NEXT(N),
     $        W(N),NDENSE(N)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
Cdense
      INTEGER THRESH
      INTEGER THRESM, MINDEN, MAXDEN, NDME
      INTEGER NBD,NBED, NBDM, LASTD, NELME
      LOGICAL IDENSE
      DOUBLE PRECISION RELDEN
Cdense
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
Cdense
Cdense
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, HASH, HMOD, I,
     $        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     $        LENJ, LN, MAXMEM, ME, MEM, MINDEG, NEL, NEWMEM,
     $        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      INTEGER P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, PN, PSRC
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      INTRINSIC MAX, MIN, MOD
C=======================================================================
C=======================================================================
Cdense
      IF (THRESH.GT.0) THEN
         THRESM  = 0
C
         RELDEN = 0.0
         DO I=1,N
             RELDEN = RELDEN + FLOAT(LEN(I))/FLOAT(N)
             THRESM = MAX(THRESM, LEN(I))
          ENDDO
         THRESM =  INT(RELDEN)*10 + (THRESM-INT(RELDEN))/10 + 1
      ELSE
         THRESM = THRESH
      ENDIF
      IF (THRESM.GE.0) THEN
       IF ((THRESM.GE.N).OR.(THRESM.LT.2)) THEN
          THRESM = N
       ENDIF
      ENDIF
      LASTD = 0
      NBD   = 0
      NBED  = 0
      NBDM  = 0
Cdense
      WFLG = 2
      MINDEG = 1
      NCMPA = 0
      NEL = 0
      HMOD = MAX (1, N-1)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      DO 10 I = 1, N
Cdense
        NDENSE(I)= 0
Cdense
        LAST (I) = 0
        HEAD (I) = 0
        NV (I) = 1
        W (I) = 1
        ELEN (I) = 0
        DEGREE (I) = LEN (I)
   10 CONTINUE
      DO 20 I = 1, N
        DEG = DEGREE (I)
        IF (DEG .GT. 0) THEN
Cdense
          IF ( (THRESM.GE.0) .AND.
     &         (DEG+1.GE.THRESM) ) THEN
            NBD = NBD+1
            IF (DEG+1.NE.N-NEL) THEN
             DEGREE(I) = DEGREE(I)+N+1
             DEG = N
             INEXT = HEAD (DEG)
             IF (INEXT .NE. 0) LAST (INEXT) = I
             NEXT (I) = INEXT
             HEAD (DEG) = I
             LAST(I)  = 0
             IF (LASTD.EQ.0) LASTD=I
            ELSE
             NBED = NBED+1
             DEGREE(I) = N+1
             DEG = N
             IF (LASTD.EQ.0) THEN
               LASTD     = I
               HEAD(DEG) = I
               NEXT(I)   = 0
               LAST(I)   = 0
             ELSE
               NEXT(LASTD) = I
               LAST(I)     = LASTD
               LASTD       = I
               NEXT(I)     = 0
             ENDIF
            ENDIF
          ELSE
Cdense
            INEXT = HEAD (DEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            HEAD (DEG) = I
          ENDIF
        ELSE
          NEL = NEL + 1
          ELEN (I) = -NEL
          PE (I) = 0
          W (I) = 0
        ENDIF
   20 CONTINUE
          IF (NBD.EQ.0) THRESM = N
C=======================================================================
C=======================================================================
   30 IF (NEL .LT. N) THEN
C=======================================================================
C=======================================================================
        DO 40 DEG = MINDEG, N
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
Cdense
        IF (DEG.LT.N)  THEN
Cdense
          INEXT = NEXT (ME)
          IF (INEXT .NE. 0) LAST (INEXT) = 0
          HEAD (DEG) = INEXT
Cdense
        ELSE
          NBDM = MAX(NBDM,NBD)
          IF (DEGREE(ME).GT.N+1) THEN
            MINDEN = NBD
            MAXDEN = 0
            IF (WFLG+NBD+1 .LE. WFLG) THEN
             DO  52 X = 1, N
              IF (W (X) .NE. 0) W (X) = 1
  52         CONTINUE
             WFLG = 2
            ENDIF
            WFLG = WFLG + 1
  51        CONTINUE
            INEXT = NEXT (ME)
            IF (INEXT .NE. 0) THEN
               LAST (INEXT) = 0
            ELSE
               LASTD = 0
            ENDIF
            NDENSE(ME) = 0
            W(ME)      = WFLG
            P1 = PE(ME)
            P2 = P1 + LEN(ME) -1
            LN       = P1
            ELN      = P1
            DO 55 P=P1,P2
              E= IW(P)
              IF (W(E).EQ.WFLG) GOTO 55
              W(E) = WFLG
              IF (PE(E).LT.0) THEN
                X = E
  53            X = -PE(X)
                IF (W(X) .EQ.WFLG) GOTO 55
                W(X) = WFLG
                IF ( PE(X) .LT. 0 ) GOTO 53
                E = X
              ENDIF
              IF (ELEN(E).LT.0) THEN
               NDENSE(E) = NDENSE(E) - NV(ME)
               IW(LN) = IW(ELN)
               IW(ELN) = E
               LN  = LN+1
               ELN = ELN + 1
               PME1 = PE(E)
               DO 54 PME = PME1, PME1+LEN(E)-1
                X = IW(PME)
                IF ((ELEN(X).GE.0).AND.(W(X).NE.WFLG)) THEN
                 NDENSE(ME) = NDENSE(ME) + NV(X)
                 W(X) = WFLG
                ENDIF
 54            CONTINUE
              ELSE
               NDENSE(ME) = NDENSE(ME) + NV(E)
               IW(LN)=E
               LN = LN+1
              ENDIF
  55        CONTINUE
            WFLG     = WFLG + 1
            LEN(ME)  = LN-P1
            ELEN(ME) = ELN- P1
            NDME = NDENSE(ME)+NV(ME)
            MINDEN = MIN (MINDEN, NDME)
            MAXDEN = MAX (MAXDEN, NDME)
            IF (NDENSE(ME).EQ.0) NDENSE(ME) =1
            DEGREE(ME) = NDENSE(ME)
            DEG = DEGREE(ME)
            MINDEG = MIN(DEG,MINDEG)
            JNEXT = HEAD(DEG)
            IF (JNEXT.NE. 0) LAST (JNEXT) = ME
            NEXT(ME) = JNEXT
            HEAD(DEG) = ME
            ME    = INEXT
            IF (ME.NE.0) THEN
              IF (DEGREE(ME).GT.(N+1) ) GOTO 51
            ENDIF
            HEAD (N) = ME
            THRESM = MAX(THRESM*2, MINDEN+(MAXDEN-MINDEN)/2)
            THRESM = MIN(THRESM,NBD)
            IF (THRESM.GE.NBD) THRESM=N
            NBD    = NBED
C
            GOTO 30
          ENDIF
          IF (DEGREE(ME).EQ.N+1) THEN
           IF (NBD.NE.NBED) THEN
            write(6,*) ' Error -1 quasi dense rows remains'
            stop
           ENDIF
           NELME    = -(NEL+1)
           DO 59 X=1,N
            IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
             PE(X) = -ME
            ELSEIF (DEGREE(X).EQ.N+1) THEN
             NEL   = NEL + NV(X)
             PE(X) = -ME
             ELEN(X) = 0
             NV(X) = 0
            ENDIF
   59      CONTINUE
           ELEN(ME) = NELME
           NV(ME)   = NBD
           PE(ME)   = 0
           IF (NEL.NE.N) THEN
            write(6,*) 'ERROR 2 detected in AMDD'
            write(6,*) ' NEL not equal to N: N, NEL =',N,NEL
            stop
           ENDIF
           GOTO 265
          ENDIF
        ENDIF
Cdense
Cdense traces
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
Cdense
        NDENSE(ME) = 0
Cdense
C=======================================================================
C=======================================================================
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + LEN (ME) - 1
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
              DEGME = DEGME + NVI
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
cdense
              IF (DEGREE(I).LE.N) THEN
cdense
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                NEXT (ILAST) = INEXT
              ELSE
                HEAD (DEGREE (I)) = INEXT
              ENDIF
Cdense
              ELSE
               NDENSE(ME) = NDENSE(ME) + NVI
              ENDIF
Cdense
            ENDIF
   60     CONTINUE
          NEWMEM = 0
        ELSE
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          DO 120 KNT1 = 1, ELENME + 1
            IF (KNT1 .GT. ELENME) THEN
              E = ME
              PJ = P
              LN = SLENME
            ELSE
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
            DO 110 KNT2 = 1, LN
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
                IF (PFREE .GT. IWLEN) THEN
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0
                  PE (E) = PJ
                  LEN (E) = LN - KNT2
                  IF (LEN (E) .EQ. 0) PE (E) = 0
                  NCMPA = NCMPA + 1
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = IW (PN)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = PE (J)
                      PE (J) = PDST
                      PDST = PDST + 1
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
                DEGME = DEGME + NVI
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
cdense
                IF (DEGREE(I).LE.N) THEN
cdense
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
                  HEAD (DEGREE (I)) = INEXT
                ENDIF
Cdense
                ELSE
                 NDENSE(ME) = NDENSE(ME) + NVI
                ENDIF
Cdense
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
              PE (E) = -ME
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = MAX (MAXMEM, MEM)
        ENDIF
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = PME2 - PME1 + 1
        IF (WFLG+N .LE. WFLG) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
C=======================================================================
Cdense
Cdense
C=======================================================================
Cdense
Cdense
        DO 150 PME = PME1, PME2
          I = IW (PME)
Cdense
          IF (DEGREE(I).GT.N) GOTO 150
Cdense
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + ELN - 1
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
Cdense
                WE = DEGREE (E) + WNVI - NDENSE(E)
Cdense
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
C=======================================================================
C=======================================================================
Cdense Lme should be read Lme(G')
        DO 180 PME = PME1, PME2
          I = IW (PME)
Cdense
          IF (DEGREE(I).GT.N) GOTO 180
Cdense
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0
          DEG = 0
          DO 160 P = P1, P2
            E = IW (P)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + E
CAMD            ELSE IF (DEXT .EQ. 0) THEN
CAMDC             aggressive absorption: e is not adjacent to me, but
CAMDC             the |Le \ Lme| is 0, so absorb it into me
CAMD              PE (E) = -ME
CAMD              W (E) = 0
Cdense
              ELSE IF ((DEXT .EQ. 0) .AND.
     &                (NDENSE(ME).EQ.NBD)) THEN
                PE (E) = -ME
                W (E)  = 0
              ELSE IF (DEXT.EQ.0) THEN
                  IW(PN) = E
                  PN     = PN+1
                  HASH   = HASH + E
Cdense
            ENDIF
  160     CONTINUE
          ELEN (I) = PN - P1 + 1
          P3 = PN
          DO 170 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
CAMD              DEG = DEG + NVJ
Cdense
              IF (DEGREE(J).LE.N) DEG=DEG+NVJ
Cdense
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + J
            ENDIF
  170     CONTINUE
CAMD          IF (DEG .EQ. 0) THEN
Cdense
          IF ((DEG .EQ. 0).AND.(NDENSE(ME).EQ.NBD)) THEN
Cdense
            PE (I) = -ME
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
CAMD            DEGREE (I) = MIN (DEGREE (I), DEG)
CAMD        modified test moved to loop 260
Cdense
            DEGREE(I) = MIN (DEG+NBD-NDENSE(ME),
     &                       DEGREE(I))
Cdense
            IW (PN) = IW (P3)
            IW (P3) = IW (P1)
            IW (P1) = ME
            LEN (I) = PN - P1 + 1
            HASH = MOD (HASH, HMOD) + 1
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = HASH
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
        DMAX = MAX (DMAX, DEGME)
        WFLG = WFLG + DMAX
        IF (WFLG+N .LE. WFLG) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
C=======================================================================
C=======================================================================
        DO 250 PME = PME1, PME2
          I = IW (PME)
CAMD          IF (NV (I) .LT. 0) THEN
Cdense
          IF ( (NV(I).LT.0) .AND. (DEGREE(I).LE.N) ) THEN
Cdense
            HASH = LAST (I)
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
              I = -J
              HEAD (HASH) = 0
            ELSE
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
              LN = LEN (I)
              ELN = ELEN (I)
              DO 210 P = PE (I) + 1, PE (I) + LN - 1
                W (IW (P)) = WFLG
  210         CONTINUE
              JLAST = I
              J = NEXT (I)
  220         CONTINUE
              IF (J .NE. 0) THEN
                IF (LEN (J) .NE. LN) GO TO 240
                IF (ELEN (J) .NE. ELN) GO TO 240
                DO 230 P = PE (J) + 1, PE (J) + LN - 1
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
                PE (J) = -I
                NV (I) = NV (I) + NV (J)
                NV (J) = 0
                ELEN (J) = 0
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
  240           CONTINUE
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
C=======================================================================
C=======================================================================
        P = PME1
        NLEFT = N - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
            NV (I) = NVI
Cdense
            IF (DEGREE(I).LE.N) THEN
Cdense
            DEG = MIN (DEGREE (I)+ DEGME - NVI, NLEFT - NVI)
            DEGREE (I) = DEG
            IDENSE = .FALSE.
            IF (THRESM.GE.0) THEN
             IF (DEG+NVI .GE. THRESM) THEN
              IF (THRESM.EQ.N) THEN
               IF ((ELEN(I).LE.2) .AND. ((DEG+NVI).EQ.NLEFT) ) THEN
                DEGREE(I) = N+1
                IDENSE = .TRUE.
               ENDIF
              ELSE
               IDENSE = .TRUE.
               IF ((ELEN(I).LE.2).AND.((DEG+NVI).EQ.NLEFT) ) THEN
                 DEGREE(I) = N+1
               ELSE
                 DEGREE(I) = N+1+DEGREE(I)
               ENDIF
              ENDIF
             ENDIF
             IF (IDENSE) THEN
               P1 = PE(I)
               P2 = P1 + ELEN(I) - 1
               IF (P2.GE.P1) THEN
               DO 264 PJ=P1,P2
                 E= IW(PJ)
                 NDENSE (E) = NDENSE(E) + NVI
 264           CONTINUE
               ENDIF
               NBD = NBD+NVI
               DEG = N
               IF (DEGREE(I).EQ.N+1) THEN
                NBED = NBED +NVI
                IF (LASTD.EQ.0) THEN
                 LASTD     = I
                 HEAD(DEG) = I
                 NEXT(I)   = 0
                 LAST(I)   = 0
                ELSE
                 NEXT(LASTD) = I
                 LAST(I)     = LASTD
                 LASTD       = I
                 NEXT(I)     = 0
                ENDIF
               ELSE
                INEXT = HEAD(DEG)
                IF (INEXT .NE. 0) LAST (INEXT) = I
                NEXT (I) = INEXT
                HEAD (DEG) = I
                LAST(I)    = 0
                IF (LASTD.EQ.0) LASTD=I
               ENDIF
             ENDIF
            ENDIF
            IF (.NOT.IDENSE) THEN
Cdense
            INEXT = HEAD (DEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            LAST (I) = 0
            HEAD (DEG) = I
Cdense
            ENDIF
Cdense
            MINDEG = MIN (MINDEG, DEG)
Cdense
            ENDIF
Cdense
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
C=======================================================================
C=======================================================================
        NV (ME) = NVPIV + DEGME
        LEN (ME) = P - PME1
        IF (LEN (ME) .EQ. 0) THEN
          PE (ME) = 0
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
          PFREE = P
          MEM = MEM - NEWMEM + LEN (ME)
        ENDIF
Cdense
C=======================================================================
      GO TO 30
      ENDIF
C=======================================================================
Cdense
  265 CONTINUE
Cdense
C=======================================================================
C=======================================================================
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
          J = -PE (I)
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = -PE (J)
              GO TO 270
            ENDIF
            E = J
            K = -ELEN (E)
            J = I
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = -PE (J)
              PE (J) = -E
              IF (ELEN (J) .EQ. 0) THEN
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
      DO 300 I = 1, N
        K = ABS (ELEN (I))
        LAST (K) = I
        ELEN (I) = K
  300 CONTINUE
C=======================================================================
C=======================================================================
      PFREE = MAXMEM
      RETURN
      END
C      SUBROUTINE METIS_NODEND(N,IPTR,IRN,METFTN,METOPT,INVPRM,PERM)
CC Dummy routine that is called if MeTiS is not linked.
C      INTEGER N
C      INTEGER IPTR(N+1),IRN(*),METFTN,METOPT(8),INVPRM(N),PERM(N)
C      PERM(1) = -1
C      write(*,*) "DUMMY METIS FCN CALLED"
C      RETURN
C      END
C *******************************************************************
C COPYRIGHT (c) 1995 Timothy A. Davis, Patrick Amestoy and
C             Council for the Central Laboratory of the Research Councils
C All rights reserved.
C
C None of the comments in this Copyright notice between the lines
C of asterisks shall be removed or altered in any way.
C
C This Package is intended for compilation without modification,
C so most of the embedded comments have been removed.
C
C ALL USE IS SUBJECT TO LICENCE. For full details of a UK ACADEMIC
C Licence, see http://hsl.rl.ac.uk/acuk/cou.html
C
C Please note that for a UK ACADEMIC Licence:
C
C 1. The Packages may only be used for academic research or teaching
C    purposes by the Licensee, and must not be copied by the Licensee for
C    use by any other persons. Use of the Packages in any commercial
C    application shall be subject to prior written agreement between
C    Hyprotech UK Limited and the Licensee on suitable terms and
C    conditions, which will include financial conditions.
C 2. All information on the Package is provided to the Licensee on the
C    understanding that the details thereof are confidential.
C 3. All publications issued by the Licensee that include results obtained
C    with the help of one or more of the Packages shall acknowledge the
C    use of the Packages. The Licensee will notify the Numerical Analysis
C    Group at Rutherford Appleton Laboratory of any such publication.
C 4. The Packages may be modified by or on behalf of the Licensee
C    for such use in research applications but at no time shall such
C    Packages or modifications thereof become the property of the
C    Licensee. The Licensee shall make available free of charge to the
C    copyright holder for any purpose all information relating to
C    any modification.
C 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
C    direct or consequential loss or damage whatsoever arising out of
C    the use of Packages by the Licensee.
C *******************************************************************
C
C Original date 30 November 1995
C  April 2001: call to MC49 changed to MC59 to make routine threadsafe
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC47AD(N, NE, PE, IW, IWLEN, MP, INFO)
      INTEGER N, NE, PE(N+1), IWLEN, IW(IWLEN), MP, INFO(7)
      INTEGER DEGREE
      DOUBLE PRECISION DUMMY(1)
      INTEGER ELEN,HEAD,I,IFLAG,II,I1,I2,J,LAST,LEN,LENIW,NCMPA,
     *        NEXT,NV,PFREE,W
      INTEGER ICT59(10),INFO59(10),IOUT,JOUT,IDUP,NZOUT
      EXTERNAL MC59AD,MC34AD,MC47BD
      INFO(1) = 0
      IF (N.LT.1) THEN
        INFO(1) = -1
        GO TO 1000
      ENDIF
      IF (PE(1).LT.1) THEN
        IF (2*NE+N.GT.IWLEN) THEN
          INFO(1) = -2
          GO TO 1000
        ENDIF
      ELSE
        IF (NE+N.GT.IWLEN) THEN
          INFO(1) = -2
          GO TO 1000
        ENDIF
      ENDIF
      IF (MP.GT.0) THEN
        WRITE(MP,'(/A)') 'Entry to MC47A/AD'
        WRITE(MP,'(A,I10,A,I10,A)') 'Matrix of order',N,' with',NE,
     *                            ' entries'
        IF (PE(1).LT.0)  THEN
          WRITE(MP,'(A)') 'Matrix input in coordinate form'
          WRITE(MP,'(A/(4(I8,I8,4X)))') 'Row and column indices',
     *          (IW(I),IW(NE+I),I=1,NE)
        ELSE
          WRITE(MP,'(A)') 'Matrix input by columns'
          DO 10 J=1,N
            WRITE(MP,'(A,I4/(10I8))') 'Column',J,
     *                                (IW(I),I=PE(J),PE(J+1)-1)
   10     CONTINUE
        ENDIF
      ENDIF
      LAST   = IWLEN  - N + 1
      ELEN   = LAST   - N
      NV     = ELEN   - N
      W      = NV     - N
      DEGREE = W      - N
      HEAD   = DEGREE - N
      NEXT   = HEAD   - N
      LEN    = NEXT   - N
      LENIW = LEN-1
      INFO(6) = 0
      INFO(7) = 0
      IF (PE(1).LT.0) THEN
        DO 20 I=1,NE
          IF (IW(I).LE.IW(NE+I)) THEN
            IF (IW(I).EQ.IW(NE+I) .AND. IW(I).NE.0) THEN
              INFO(7) = INFO(7) + 1
            ELSE
              IF (IW(I).GT.0) INFO(6) = INFO(6) + 1
            ENDIF
            IW(I)=0
          ENDIF
   20   CONTINUE
        ICT59(1) = 0
        ICT59(2) = 1
        ICT59(3) = 1
        ICT59(4) = MP
        ICT59(5) = -1
        ICT59(6) = 0
        CALL MC59AD(ICT59,N,N,NE,IW,NE,IW(NE+1),1,DUMMY,
     *              N+1,PE,N+1,IW(2*NE+1),INFO59)
        IFLAG = INFO59(1)
        IDUP  = INFO59(3)
        IOUT  = INFO59(4)
        JOUT  = INFO59(5)
        NZOUT = INFO59(6)
      ELSE
        IDUP = 0
        IOUT = 0
        JOUT = 0
        DO 30 I = 1,N
          IW(NE+I) = 0
   30   CONTINUE
        DO 50 J=1,N
          I1 = PE(J)
          PE(J) = I1-(IOUT+IDUP)
          I2 = PE(J+1)-1
          IF (I2.LT.I1-1) THEN
            INFO(1) = -3
            GO TO 1000
          ENDIF
          DO 40 II = I1,I2
            I = IW(II)
            IF (I.LE.J .OR. I.GT.N) THEN
              IF (I.EQ.J) INFO(7) = INFO(7) + 1
              IF (I.GT.0 .AND. I.LT.J) INFO(6) = INFO(6) + 1
              IOUT = IOUT + 1
            ELSE
              IF (IW(NE+I).EQ.J) THEN
                IDUP = IDUP + 1
              ELSE
                IW(NE+I)=J
                IW(II-(IOUT+IDUP)) = I
              ENDIF
            ENDIF
   40     CONTINUE
   50   CONTINUE
        PE(N+1) = NE - (IOUT+IDUP) + 1
      ENDIF
      IF (IDUP.GT.0) THEN
        INFO(1) = 1
        INFO(4) = IDUP
      ELSE
        INFO(4) = 0
      ENDIF
      IF (IOUT.GT.0 .OR. JOUT.GT.0) THEN
        INFO(1) = 1
        INFO(5) = IOUT + JOUT - INFO(7)
      ELSE
        INFO(5) = 0
      ENDIF
      IF (INFO(6).GT.0 .OR. INFO(7).GT.0) INFO(1) = 1
      IF (NE-(IOUT+IDUP).EQ.0) THEN
        INFO(1) = -4
        GO TO 1000
      ENDIF
      IF (LENIW.LT.2*(PE(N+1)-1)) THEN
        INFO(1) = -2
        GO TO 1000
      ENDIF
      CALL MC34AD(N,IW,PE,.FALSE.,DUMMY,IW(W))
      PFREE = PE(N+1)
      DO 60 I=1,N
        IW(LEN+I-1) = PE(I+1) - PE(I)
   60 CONTINUE
      CALL MC47BD(N,LENIW,PE,PFREE,IW(LEN),IW,IW(NV),IW(ELEN),
     *            IW(LAST),NCMPA,IW(DEGREE),IW(HEAD),IW(NEXT),IW(W))
      INFO(2) = NCMPA
      INFO(3) = PFREE+8*N
      IF (MP.GT.0) THEN
        WRITE(MP,'(/A)') 'Exit from MC47A/AD'
        WRITE(MP,'(A/7I10)') 'INFO(1-7):',(INFO(I),I=1,7)
        WRITE(MP,'(A/(8I10))') 'Parent array',(PE(I),I=1,N)
        WRITE(MP,'(A/(8I10))') 'Permutation',(IW(ELEN+I-1),I=1,N)
        WRITE(MP,'(A/(8I10))') 'Inverse permutation',
     *                         (IW(LAST+I-1),I=1,N)
        WRITE(MP,'(A/(8I10))') 'Degree array',(IW(NV+I-1),I=1,N)
      ENDIF
 1000 RETURN
      END
      SUBROUTINE MC47BD (N, IWLEN, PE, PFREE, LEN, IW, NV, ELEN,
     $                   LAST, NCMPA, DEGREE, HEAD, NEXT, W)
      INTEGER N, IWLEN, PE(N), PFREE, LEN(N), IW(IWLEN), NV(N),
     $        ELEN(N), LAST(N), NCMPA, DEGREE(N), HEAD(N), NEXT(N),
     $        W(N)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, HASH, HMOD, I,
     $        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     $        LENJ, LN, MAXMEM, ME, MEM, MINDEG, NEL, NEWMEM,
     $        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      INTEGER P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, PN, PSRC
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      INTRINSIC MAX, MIN, MOD
C=======================================================================
C=======================================================================
      WFLG = 2
      MINDEG = 1
      NCMPA = 0
      NEL = 0
      HMOD = MAX (1, N-1)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      DO 10 I = 1, N
        LAST (I) = 0
        HEAD (I) = 0
        NV (I) = 1
        W (I) = 1
        ELEN (I) = 0
        DEGREE (I) = LEN (I)
   10 CONTINUE
      DO 20 I = 1, N
        DEG = DEGREE (I)
        IF (DEG .GT. 0) THEN
          INEXT = HEAD (DEG)
          IF (INEXT .NE. 0) LAST (INEXT) = I
          NEXT (I) = INEXT
          HEAD (DEG) = I
        ELSE
          NEL = NEL + 1
          ELEN (I) = -NEL
          PE (I) = 0
          W (I) = 0
        ENDIF
   20 CONTINUE
C=======================================================================
C=======================================================================
   30 IF (NEL .LT. N) THEN
C=======================================================================
C=======================================================================
        DO 40 DEG = MINDEG, N
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
        INEXT = NEXT (ME)
        IF (INEXT .NE. 0) LAST (INEXT) = 0
        HEAD (DEG) = INEXT
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
C=======================================================================
C=======================================================================
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + LEN (ME) - 1
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
              DEGME = DEGME + NVI
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                NEXT (ILAST) = INEXT
              ELSE
                HEAD (DEGREE (I)) = INEXT
              ENDIF
            ENDIF
   60     CONTINUE
          NEWMEM = 0
        ELSE
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          DO 120 KNT1 = 1, ELENME + 1
            IF (KNT1 .GT. ELENME) THEN
              E = ME
              PJ = P
              LN = SLENME
            ELSE
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
            DO 110 KNT2 = 1, LN
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
                IF (PFREE .GT. IWLEN) THEN
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0
                  PE (E) = PJ
                  LEN (E) = LN - KNT2
                  IF (LEN (E) .EQ. 0) PE (E) = 0
                  NCMPA = NCMPA + 1
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = IW (PN)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = PE (J)
                      PE (J) = PDST
                      PDST = PDST + 1
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
                DEGME = DEGME + NVI
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
                  HEAD (DEGREE (I)) = INEXT
                ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
              PE (E) = -ME
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = MAX (MAXMEM, MEM)
        ENDIF
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = PME2 - PME1 + 1
        IF (WFLG+N .LE. WFLG) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
C=======================================================================
C=======================================================================
        DO 150 PME = PME1, PME2
          I = IW (PME)
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + ELN - 1
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
                WE = DEGREE (E) + WNVI
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
C=======================================================================
C=======================================================================
        DO 180 PME = PME1, PME2
          I = IW (PME)
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0
          DEG = 0
          DO 160 P = P1, P2
            E = IW (P)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + E
            ELSE IF (DEXT .EQ. 0) THEN
              PE (E) = -ME
              W (E) = 0
            ENDIF
  160     CONTINUE
          ELEN (I) = PN - P1 + 1
          P3 = PN
          DO 170 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
              DEG = DEG + NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + J
            ENDIF
  170     CONTINUE
          IF (DEG .EQ. 0) THEN
            PE (I) = -ME
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
            DEGREE (I) = MIN (DEGREE (I), DEG)
            IW (PN) = IW (P3)
            IW (P3) = IW (P1)
            IW (P1) = ME
            LEN (I) = PN - P1 + 1
            HASH = MOD (HASH, HMOD) + 1
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = HASH
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
        DMAX = MAX (DMAX, DEGME)
        WFLG = WFLG + DMAX
        IF (WFLG+N .LE. WFLG) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
C=======================================================================
C=======================================================================
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF (NV (I) .LT. 0) THEN
            HASH = LAST (I)
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
              I = -J
              HEAD (HASH) = 0
            ELSE
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
              LN = LEN (I)
              ELN = ELEN (I)
              DO 210 P = PE (I) + 1, PE (I) + LN - 1
                W (IW (P)) = WFLG
  210         CONTINUE
              JLAST = I
              J = NEXT (I)
  220         CONTINUE
              IF (J .NE. 0) THEN
                IF (LEN (J) .NE. LN) GO TO 240
                IF (ELEN (J) .NE. ELN) GO TO 240
                DO 230 P = PE (J) + 1, PE (J) + LN - 1
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
                PE (J) = -I
                NV (I) = NV (I) + NV (J)
                NV (J) = 0
                ELEN (J) = 0
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
  240           CONTINUE
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
C=======================================================================
C=======================================================================
        P = PME1
        NLEFT = N - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
            NV (I) = NVI
            DEG = MIN (DEGREE (I) + DEGME - NVI, NLEFT - NVI)
            INEXT = HEAD (DEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            LAST (I) = 0
            HEAD (DEG) = I
            MINDEG = MIN (MINDEG, DEG)
            DEGREE (I) = DEG
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
C=======================================================================
C=======================================================================
        NV (ME) = NVPIV + DEGME
        LEN (ME) = P - PME1
        IF (LEN (ME) .EQ. 0) THEN
          PE (ME) = 0
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
          PFREE = P
          MEM = MEM - NEWMEM + LEN (ME)
        ENDIF
C=======================================================================
      GO TO 30
      ENDIF
C=======================================================================
C=======================================================================
C=======================================================================
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
          J = -PE (I)
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = -PE (J)
              GO TO 270
            ENDIF
            E = J
            K = -ELEN (E)
            J = I
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = -PE (J)
              PE (J) = -E
              IF (ELEN (J) .EQ. 0) THEN
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
      DO 300 I = 1, N
        K = ABS (ELEN (I))
        LAST (K) = I
        ELEN (I) = K
  300 CONTINUE
C=======================================================================
C=======================================================================
      PFREE = MAXMEM
      RETURN
      END
* *******************************************************************
* COPYRIGHT (c) 1988 Hyprotech UK and
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of a UK ACADEMIC
* Licence, see http://hsl.rl.ac.uk/acuk/cou.html
*
* Please note that for a UK ACADEMIC Licence:
*
* 1. The Packages may only be used for academic research or teaching
*    purposes by the Licensee, and must not be copied by the Licensee for
*    use by any other persons. Use of the Packages in any commercial
*    application shall be subject to prior written agreement between
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
* Council for the Central Laboratory of the Research Councils
C Original date 14 June 2001
C  June 2001: threadsafe version of MC41
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC71AD(N,KASE,X,EST,W,IW,KEEP)
      INTEGER ITMAX
      PARAMETER (ITMAX=5)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      DOUBLE PRECISION EST
      INTEGER KASE,N
      DOUBLE PRECISION W(*),X(*)
      INTEGER IW(*),KEEP(5)
      DOUBLE PRECISION ALTSGN,TEMP
      INTEGER I,ITER,J,JLAST,JUMP
      INTEGER IDAMAX
      EXTERNAL IDAMAX
      INTRINSIC ABS,SIGN,NINT,DBLE
      IF (N.LE.0) THEN
        KASE = -1
        RETURN
      END IF
      IF (KASE.EQ.0) THEN
        DO 10 I = 1,N
          X(I) = ONE/DBLE(N)
   10   CONTINUE
        KASE = 1
        JUMP = 1
        KEEP(1) = JUMP
        KEEP(2) = 0
        KEEP(3) = 0
        KEEP(4) = 0
        RETURN
      END IF
      JUMP  = KEEP(1)
      ITER  = KEEP(2)
      J     = KEEP(3)
      JLAST = KEEP(4)
      GO TO (100,200,300,400,500) JUMP
  100 CONTINUE
      IF (N.EQ.1) THEN
        W(1) = X(1)
        EST = ABS(W(1))
        GO TO 510
      END IF
      DO 110 I = 1,N
        X(I) = SIGN(ONE,X(I))
        IW(I) = NINT(X(I))
  110 CONTINUE
      KASE = 2
      JUMP = 2
      GO TO 1010
  200 CONTINUE
      J = IDAMAX(N,X,1)
      ITER = 2
  220 CONTINUE
      DO 230 I = 1,N
        X(I) = ZERO
  230 CONTINUE
      X(J) = ONE
      KASE = 1
      JUMP = 3
      GO TO 1010
  300 CONTINUE
      DO 310 I = 1,N
        W(I) = X(I)
  310 CONTINUE
      DO 320 I = 1,N
        IF (NINT(SIGN(ONE,X(I))).NE.IW(I)) GO TO 330
  320 CONTINUE
      GO TO 410
  330 CONTINUE
      DO 340 I = 1,N
        X(I) = SIGN(ONE,X(I))
        IW(I) = NINT(X(I))
  340 CONTINUE
      KASE = 2
      JUMP = 4
      GO TO 1010
  400 CONTINUE
      JLAST = J
      J = IDAMAX(N,X,1)
      IF ((ABS(X(JLAST)).NE.ABS(X(J))) .AND. (ITER.LT.ITMAX)) THEN
        ITER = ITER + 1
        GO TO 220
      END IF
  410 CONTINUE
      EST = ZERO
      DO 420 I = 1,N
        EST = EST + ABS(W(I))
  420 CONTINUE
      ALTSGN = ONE
      DO 430 I = 1,N
        X(I) = ALTSGN* (ONE+DBLE(I-1)/DBLE(N-1))
        ALTSGN = -ALTSGN
  430 CONTINUE
      KASE = 1
      JUMP = 5
      GO TO 1010
  500 CONTINUE
      TEMP = ZERO
      DO 520 I = 1,N
        TEMP = TEMP + ABS(X(I))
  520 CONTINUE
      TEMP = 2.0*TEMP/DBLE(3*N)
      IF (TEMP.GT.EST) THEN
        DO 530 I = 1,N
          W(I) = X(I)
  530   CONTINUE
        EST = TEMP
      END IF
  510 KASE = 0
 1010 CONTINUE
      KEEP(1) = JUMP
      KEEP(2) = ITER
      KEEP(3) = J
      KEEP(4) = JLAST
      RETURN
      END
* *******************************************************************
* COPYRIGHT (c) 1987 Hyprotech UK
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of a UK ACADEMIC
* Licence, see http://hsl.rl.ac.uk/acuk/cou.html
*
* Please note that for a UK ACADEMIC Licence:
*
* 1. The Packages may only be used for academic research or teaching
*    purposes by the Licensee, and must not be copied by the Licensee for
*    use by any other persons. Use of the Packages in any commercial
*    application shall be subject to prior written agreement between
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
* Original date 10 Feb 1993
C       Toolpack tool decs employed.
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC34AD(N,IRN,JCOLST,YESA,A,IW)
      INTEGER N
      LOGICAL YESA
      DOUBLE PRECISION A(*)
      INTEGER IRN(*),IW(*),JCOLST(*)
      INTEGER CKP1,I,I1,I2,II,IPKP1,IPOS,J,JSTART,LENK,NDIAG,NEWTAU,
     +        OLDTAU
      OLDTAU = JCOLST(N+1) - 1
      DO 5 I = 1,N
        IW(I) = 0
    5 CONTINUE
      NDIAG = 0
      DO 20 J = 1,N
        I1 = JCOLST(J)
        I2 = JCOLST(J+1) - 1
        IW(J) = IW(J) + I2 - I1 + 1
        DO 10 II = I1,I2
          I = IRN(II)
          IF (I.NE.J) THEN
            IW(I) = IW(I) + 1
          ELSE
            NDIAG = NDIAG + 1
          END IF
   10   CONTINUE
   20 CONTINUE
      NEWTAU = 2*OLDTAU - NDIAG
      IPKP1 = OLDTAU + 1
      CKP1 = NEWTAU + 1
      DO 40 J = N,1,-1
        I1 = JCOLST(J)
        I2 = IPKP1
        LENK = I2 - I1
        JSTART = CKP1
        IPKP1 = I1
        I2 = I2 - 1
        DO 30 II = I2,I1,-1
          JSTART = JSTART - 1
          IF (YESA) A(JSTART) = A(II)
          IRN(JSTART) = IRN(II)
   30   CONTINUE
        JCOLST(J) = JSTART
        CKP1 = CKP1 - IW(J)
        IW(J) = LENK
   40 CONTINUE
      DO 80 J = N,1,-1
        I1 = JCOLST(J)
        I2 = JCOLST(J) + IW(J) - 1
        DO 60 II = I1,I2
          I = IRN(II)
          IF (I.EQ.J) GO TO 60
          JCOLST(I) = JCOLST(I) - 1
          IPOS = JCOLST(I)
          IF (YESA) A(IPOS) = A(II)
          IRN(IPOS) = J
   60   CONTINUE
   80 CONTINUE
      JCOLST(N+1) = NEWTAU + 1
      RETURN
      END
* *******************************************************************
* COPYRIGHT (c) 1993 Council for the Central Laboratory
*                    of the Research Councils
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of a UK ACADEMIC
* Licence, see http://hsl.rl.ac.uk/acuk/cou.html
*
* Please note that for a UK ACADEMIC Licence:
*
* 1. The Packages may only be used for academic research or teaching
*    purposes by the Licensee, and must not be copied by the Licensee for
*    use by any other persons. Use of the Packages in any commercial
*    application shall be subject to prior written agreement between
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
C Original date 29 Jan 2001
C 29 January 2001. Modified from MC49 to be threadsafe.

C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC59AD(ICNTL,NC,NR,NE,IRN,LJCN,JCN,LA,A,LIP,IP,
     &                  LIW,IW,INFO)
      INTEGER LA,LIP,LIW,LJCN,NC,NE,NR
      DOUBLE PRECISION A(LA)
      INTEGER ICNTL(10),IP(LIP),INFO(10),IRN(NE),IW(LIW),JCN(LJCN)
      INTEGER I,ICNTL1,ICNTL2,ICNTL3,ICNTL6,LAA
      INTEGER IDUP,IOUT,IUP,JOUT,LP,MP,KNE,PART
      LOGICAL LCHECK
      EXTERNAL MC59BD,MC59CD,MC59DD,MC59ED,MC59FD
      INTRINSIC MAX
      DO 10 I = 1,10
         INFO(I) = 0
   10 CONTINUE
      ICNTL1 = ICNTL(1)
      ICNTL2 = ICNTL(2)
      ICNTL3 = ICNTL(3)
      ICNTL6 = ICNTL(6)
      LCHECK = (ICNTL1.EQ.0)
      LP = ICNTL(4)
      MP = ICNTL(5)
      IF (ICNTL2.GT.2 .OR. ICNTL2.LT.0) THEN
         INFO(1) = -1
         INFO(2) = ICNTL2
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9010) ICNTL2
         END IF
         GO TO 70
      END IF
      IF (ICNTL6.GT.2 .OR. ICNTL6.LT.-2) THEN
         INFO(1) = -11
         INFO(2) = ICNTL6
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9150) ICNTL6
         END IF
         GO TO 70
      END IF
      IF (NC.LT.1) THEN
        INFO(1) = -2
        INFO(2) = NC
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9020) NC
        END IF
        GO TO 70
      END IF
      IF (NR.LT.1) THEN
        INFO(1) = -3
        INFO(2) = NR
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9030) NR
        END IF
        GO TO 70
      END IF
      IF (ICNTL6.NE.0 .AND. NR.NE.NC) THEN
        INFO(1) = -3
        INFO(2) = NR
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9035) NC,NR
        END IF
        GO TO 70
      END IF
      IF (NE.LT.1) THEN
        INFO(1) = -4
        INFO(2) = NE
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9040) NE
        END IF
        GO TO 70
      END IF
      IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.1) THEN
        IF (LJCN.LT.NE) THEN
          INFO(1) = -5
          INFO(2) = NE
        END IF
      ELSE
        IF (LJCN.LT.1) THEN
          INFO(1) = -5
          INFO(2) = 1
        END IF
      END IF
      IF (INFO(1).EQ.-5) THEN
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9050) LJCN,INFO(2)
         END IF
         GO TO 70
      END IF
      IF (ICNTL3.EQ.0) THEN
        IF (LA.LT.NE) THEN
          INFO(1) = -6
          INFO(2) = NE
        END IF
      ELSE
        IF (LA.LT.1) THEN
          INFO(1) = -6
          INFO(2) = 1
        END IF
      END IF
      IF (INFO(1).EQ.-6) THEN
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9060) LA,INFO(2)
         END IF
         GO TO 70
      END IF
      IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.2) THEN
        IF (LIP.LT.NC+1) THEN
          INFO(1) = -7
          INFO(2) = NC+1
        END IF
      ELSE IF (LIP.LT.MAX(NR,NC)+1) THEN
        INFO(1) = -7
        INFO(2) = MAX(NR,NC)+1
      END IF
      IF (INFO(1).EQ.-7) THEN
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9065) LIP,INFO(2)
        END IF
        GO TO 70
      END IF
      IF (LIW.LT.MAX(NR,NC)+1) THEN
        INFO(1) = -8
        INFO(2) = MAX(NR,NC)+1
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9070) LIW,INFO(2)
        END IF
        GO TO 70
      END IF
      LAA = NE
      IF (ICNTL3.NE.0) LAA = 1
      IOUT = 0
      JOUT = 0
      IDUP = 0
      IUP = 0
      PART = 0
      IF (ICNTL6.NE.0) PART = 1
      IF (ICNTL2.EQ.0) THEN
        CALL MC59BD(LCHECK,PART,NC,NR,NE,IRN,JCN,LAA,A,IP,IW,
     +              IOUT,JOUT,KNE)
        IF (KNE.EQ.0) GO TO 50
        IF (LCHECK) CALL MC59ED(NC,NR,NE,IRN,LIP,IP,LAA,A,IW,IDUP,
     &                          KNE,ICNTL6)
      ELSE IF (ICNTL2.EQ.1) THEN
        IF (ICNTL6.NE.0) PART = -1
        CALL MC59BD(LCHECK,PART,NR,NC,NE,JCN,IRN,LAA,A,IW,IP,
     +              JOUT,IOUT,KNE)
        IF (KNE.EQ.0) GO TO 50
        IF (LCHECK) CALL MC59ED(NR,NC,NE,JCN,NR+1,IW,LAA,A,IP,
     &                          IDUP,KNE,ICNTL6)
        CALL MC59CD(NC,NR,KNE,IRN,JCN,LAA,A,IP,IW)
      ELSE IF (ICNTL2.EQ.2) THEN
        IF (LCHECK) THEN
          CALL MC59FD(NC,NR,NE,IRN,NC+1,IP,LAA,A,LIW,IW,IDUP,
     +                IOUT,IUP,KNE,ICNTL6,INFO)
          IF (INFO(1).EQ.-9) GO TO 40
          IF (KNE.EQ.0) GO TO 50
        ELSE
           KNE = NE
        END IF
        CALL MC59DD(NC,KNE,IRN,IP,LAA,A)
      END IF
      INFO(3) = IDUP
      INFO(4) = IOUT
      INFO(5) = JOUT
      INFO(6) = KNE
      INFO(7) = IUP
      IF (IDUP.GT.0) INFO(1) = INFO(1) + 1
      IF (IOUT.GT.0) INFO(1) = INFO(1) + 2
      IF (JOUT.GT.0) INFO(1) = INFO(1) + 4
      IF (INFO(1).GT.0 .AND. MP.GT.0) THEN
        WRITE (MP,FMT=9080) INFO(1)
        IF (IOUT.GT.0) WRITE (MP,FMT=9090) IOUT
        IF (JOUT.GT.0) WRITE (MP,FMT=9110) JOUT
        IF (IDUP.GT.0) WRITE (MP,FMT=9100) IDUP
        IF (IUP.GT.0)  WRITE (MP,FMT=9130) IUP
      END IF
      GO TO 70
   40 INFO(3) = IDUP
      INFO(4) = IOUT
      INFO(7) = IUP
      IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9140)
      END IF
      GO TO 70
   50 INFO(1) = -10
      INFO(4) = IOUT
      INFO(5) = JOUT
      INFO(2) = IOUT + JOUT
      IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9120)
      END IF
   70 RETURN
 9000 FORMAT (/,' *** Error return from MC59AD *** INFO(1) = ',I3)
 9010 FORMAT (1X,'ICNTL(2) = ',I2,' is out of range')
 9020 FORMAT (1X,'NC = ',I6,' is out of range')
 9030 FORMAT (1X,'NR = ',I6,' is out of range')
 9035 FORMAT (1X,'Symmetric case. NC = ',I6,' but NR = ',I6)
 9040 FORMAT (1X,'NE = ',I10,' is out of range')
 9050 FORMAT (1X,'Increase LJCN from ',I10,' to at least ',I10)
 9060 FORMAT (1X,'Increase LA from ',I10,' to at least ',I10)
 9065 FORMAT (1X,'Increase LIP from ',I8,' to at least ',I10)
 9070 FORMAT (1X,'Increase LIW from ',I8,' to at least ',I10)
 9080 FORMAT (/,' *** Warning message from MC59AD *** INFO(1) = ',I3)
 9090 FORMAT (1X,I8,' entries in IRN supplied by the user were ',
     +       /,'       out of range and were ignored by the routine')
 9100 FORMAT (1X,I8,' duplicate entries were supplied by the user')
 9110 FORMAT (1X,I8,' entries in JCN supplied by the user were ',
     +       /,'       out of range and were ignored by the routine')
 9120 FORMAT (1X,'All entries out of range')
 9130 FORMAT (1X,I8,' of these entries were in the upper triangular ',
     +       /,'       part of matrix')
 9140 FORMAT (1X,'Entries in IP are not monotonic increasing')
 9150 FORMAT (1X,'ICNTL(6) = ',I2,' is out of range')
      END
C***********************************************************************
      SUBROUTINE MC59BD(LCHECK,PART,NC,NR,NE,IRN,JCN,LA,A,IP,IW,IOUT,
     +                  JOUT,KNE)
      INTEGER LA,NC,NE,NR,IOUT,JOUT,KNE,PART
      LOGICAL LCHECK
      DOUBLE PRECISION A(LA)
      INTEGER IP(NC+1),IRN(NE),IW(NC+1),JCN(NE)
      DOUBLE PRECISION ACE,ACEP
      INTEGER I,ICE,ICEP,J,JCE,JCEP,K,L,LOC
      DO 10 J = 1,NC + 1
        IW(J) = 0
   10 CONTINUE
      KNE = 0
      IOUT = 0
      JOUT = 0
      IF (LCHECK) THEN
        IF (LA.GT.1) THEN
          IF (PART.EQ.0) THEN
            DO 20 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IRN(KNE) = I
                JCN(KNE) = J
                A(KNE) = A(K)
                IW(J) = IW(J) + 1
              END IF
   20       CONTINUE
          ELSE IF (PART.EQ.1) THEN
            DO 21 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.LT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
                A(KNE) = A(K)
              END IF
   21       CONTINUE
          ELSE IF (PART.EQ.-1) THEN
            DO 22 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.GT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
                A(KNE) = A(K)
              END IF
   22       CONTINUE
          END IF
        ELSE
          IF (PART.EQ.0) THEN
            DO 25 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IRN(KNE) = I
                JCN(KNE) = J
                IW(J) = IW(J) + 1
              END IF
   25       CONTINUE
          ELSE IF (PART.EQ.1) THEN
            DO 26 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.LT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
              END IF
   26       CONTINUE
          ELSE IF (PART.EQ.-1) THEN
            DO 27 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.GT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
              END IF
   27       CONTINUE
          END IF
        END IF
        IF (KNE.EQ.0) GO TO 130
      ELSE
        KNE = NE
        IF (PART.EQ.0) THEN
          DO 30 K = 1,NE
            J = JCN(K)
            IW(J) = IW(J) + 1
   30     CONTINUE
        ELSE IF (PART.EQ.1) THEN
          DO 35 K = 1,NE
            I = IRN(K)
            J = JCN(K)
            IF (I.LT.J) THEN
               IRN(K) = J
               JCN(K) = I
               IW(I) = IW(I) + 1
            ELSE
              IW(J) = IW(J) + 1
            END IF
   35     CONTINUE
        ELSE IF (PART.EQ.-1) THEN
          DO 36 K = 1,NE
            I = IRN(K)
            J = JCN(K)
            IF (I.GT.J) THEN
               IRN(K) = J
               JCN(K) = I
               IW(I) = IW(I) + 1
            ELSE
              IW(J) = IW(J) + 1
            END IF
   36     CONTINUE
        END IF
      END IF
      IP(1) = 1
      DO 37 J = 2,NC + 1
        IP(J) = IW(J-1) + IP(J-1)
        IW(J-1) = IP(J-1)
   37 CONTINUE
      IF (LA.EQ.1) THEN
        DO 70 L = 1,NC
          DO 60 K = IW(L),IP(L+1) - 1
            ICE = IRN(K)
            JCE = JCN(K)
            DO 40 J = 1,NE
              IF (JCE.EQ.L) GO TO 50
              LOC = IW(JCE)
              JCEP = JCN(LOC)
              ICEP = IRN(LOC)
              IW(JCE) = LOC + 1
              JCN(LOC) = JCE
              IRN(LOC) = ICE
              JCE = JCEP
              ICE = ICEP
   40       CONTINUE
   50       JCN(K) = JCE
            IRN(K) = ICE
   60     CONTINUE
   70   CONTINUE
      ELSE
        DO 120 L = 1,NC
          DO 110 K = IW(L),IP(L+1) - 1
            ICE = IRN(K)
            JCE = JCN(K)
            ACE = A(K)
            DO 90 J = 1,NE
              IF (JCE.EQ.L) GO TO 100
              LOC = IW(JCE)
              JCEP = JCN(LOC)
              ICEP = IRN(LOC)
              IW(JCE) = LOC + 1
              JCN(LOC) = JCE
              IRN(LOC) = ICE
              JCE = JCEP
              ICE = ICEP
              ACEP = A(LOC)
              A(LOC) = ACE
              ACE = ACEP
   90       CONTINUE
  100       JCN(K) = JCE
            IRN(K) = ICE
            A(K) = ACE
  110     CONTINUE
  120   CONTINUE
      END IF
  130 CONTINUE
      RETURN
      END
C**********************************************************
      SUBROUTINE MC59CD(NC,NR,NE,IRN,JCN,LA,A,IP,IW)
      INTEGER LA,NC,NE,NR
      DOUBLE PRECISION A(LA)
      INTEGER IP(NC+1),IRN(NE),IW(NR+1),JCN(NE)
      DOUBLE PRECISION ACE,ACEP
      INTEGER I,ICE,ICEP,J,J1,J2,K,L,LOC,LOCP
      DO 10 J = 1,NC
        IP(J) = 0
   10 CONTINUE
      IF (LA.GT.1) THEN
        DO 20 K = 1,NE
          I = JCN(K)
          IP(I) = IP(I) + 1
          IRN(K) = JCN(K)
   20   CONTINUE
        IP(NC+1) = NE + 1
        IP(1) = IP(1) + 1
        DO 30 J = 2,NC
          IP(J) = IP(J) + IP(J-1)
   30   CONTINUE
        DO 50 I = NR,1,-1
          J1 = IW(I)
          J2 = IW(I+1) - 1
          DO 40 J = J1,J2
            K = IRN(J)
            L = IP(K) - 1
            JCN(J) = L
            IRN(J) = I
            IP(K) = L
   40     CONTINUE
   50   CONTINUE
        IP(NC+1) = NE + 1
        DO 70 J = 1,NE
          LOC = JCN(J)
          IF (LOC.EQ.0) GO TO 70
          ICE = IRN(J)
          ACE = A(J)
          JCN(J) = 0
          DO 60 K = 1,NE
            LOCP = JCN(LOC)
            ICEP = IRN(LOC)
            ACEP = A(LOC)
            JCN(LOC) = 0
            IRN(LOC) = ICE
            A(LOC) = ACE
            IF (LOCP.EQ.0) GO TO 70
            ICE = ICEP
            ACE = ACEP
            LOC = LOCP
   60     CONTINUE
   70   CONTINUE
      ELSE
        DO 90 K = 1,NE
          I = JCN(K)
          IP(I) = IP(I) + 1
   90   CONTINUE
        IP(NC+1) = NE + 1
        IP(1) = IP(1) + 1
        DO 100 J = 2,NC
          IP(J) = IP(J) + IP(J-1)
  100   CONTINUE
        DO 120 I = NR,1,-1
          J1 = IW(I)
          J2 = IW(I+1) - 1
          DO 110 J = J1,J2
            K = JCN(J)
            L = IP(K) - 1
            IRN(L) = I
            IP(K) = L
  110     CONTINUE
  120   CONTINUE
      END IF
      RETURN
      END
C**********************************************************
      SUBROUTINE MC59DD(NC,NE,IRN,IP,LA,A)
      INTEGER LA,NC,NE
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(NC)
      DOUBLE PRECISION ACE
      INTEGER ICE,IK,J,JJ,K,KDUMMY,KLO,KMAX,KOR
      INTRINSIC ABS
      IF (LA.GT.1) THEN
        KMAX = NE
        DO 50 JJ = 1,NC
          J = NC + 1 - JJ
          KLO = IP(J) + 1
          IF (KLO.GT.KMAX) GO TO 40
          KOR = KMAX
          DO 30 KDUMMY = KLO,KMAX
            ACE = A(KOR-1)
            ICE = IRN(KOR-1)
            DO 10 K = KOR,KMAX
              IK = IRN(K)
              IF (ABS(ICE).LE.ABS(IK)) GO TO 20
              IRN(K-1) = IK
              A(K-1) = A(K)
   10       CONTINUE
            K = KMAX + 1
   20       IRN(K-1) = ICE
            A(K-1) = ACE
            KOR = KOR - 1
   30     CONTINUE
   40     KMAX = KLO - 2
   50   CONTINUE
      ELSE
        KMAX = NE
        DO 150 JJ = 1,NC
          J = NC + 1 - JJ
          KLO = IP(J) + 1
          IF (KLO.GT.KMAX) GO TO 140
          KOR = KMAX
          DO 130 KDUMMY = KLO,KMAX
            ICE = IRN(KOR-1)
            DO 110 K = KOR,KMAX
              IK = IRN(K)
              IF (ABS(ICE).LE.ABS(IK)) GO TO 120
              IRN(K-1) = IK
  110       CONTINUE
            K = KMAX + 1
  120       IRN(K-1) = ICE
            KOR = KOR - 1
  130     CONTINUE
  140     KMAX = KLO - 2
  150   CONTINUE
      END IF
      END
C***********************************************************************
      SUBROUTINE MC59ED(NC,NR,NE,IRN,LIP,IP,LA,A,IW,IDUP,KNE,ICNTL6)
      INTEGER ICNTL6,IDUP,KNE,LIP,LA,NC,NR,NE
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(LIP),IW(NR)
      INTEGER I,J,K,KSTART,KSTOP,NZJ
      IDUP = 0
      KNE = 0
      DO 10 I = 1,NR
        IW(I) = 0
   10 CONTINUE
      KSTART = IP(1)
      IF (LA.GT.1) THEN
        NZJ = 0
        DO 30 J = 1,NC
          KSTOP = IP(J+1)
          IP(J+1) = IP(J)
          DO 20 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (IW(I).LE.NZJ) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              A(KNE) = A(K)
              IP(J+1) = IP(J+1) + 1
              IW(I) = KNE
            ELSE
              IDUP = IDUP + 1
              IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
            END IF
   20     CONTINUE
          KSTART = KSTOP
          NZJ = KNE
   30   CONTINUE
      ELSE
        DO 50 J = 1,NC
          KSTOP = IP(J+1)
          IP(J+1) = IP(J)
          DO 40 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (IW(I).LT.J) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              IP(J+1) = IP(J+1) + 1
              IW(I) = J
            ELSE
              IDUP = IDUP + 1
            END IF
   40     CONTINUE
          KSTART = KSTOP
   50   CONTINUE
      END IF
      RETURN
      END
C***********************************************************************
      SUBROUTINE MC59FD(NC,NR,NE,IRN,LIP,IP,LA,A,LIW,IW,IDUP,IOUT,
     +                  IUP,KNE,ICNTL6,INFO)
      INTEGER ICNTL6,IDUP,IOUT,IUP,KNE,LA,LIP,LIW,NC,NR,NE
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(LIP),IW(LIW),INFO(2)
      INTEGER I,J,K,KSTART,KSTOP,NZJ,LOWER
      IDUP = 0
      IOUT = 0
      IUP = 0
      KNE = 0
      DO 10 I = 1,NR
        IW(I) = 0
   10 CONTINUE
      KSTART = IP(1)
      LOWER = 1
      IF (LA.GT.1) THEN
        NZJ = 0
        DO 30 J = 1,NC
          IF (ICNTL6.NE.0) LOWER = J
          KSTOP = IP(J+1)
          IF (KSTART.GT.KSTOP) THEN
            INFO(1) = -9
            INFO(2) = J
            RETURN
          END IF
          IP(J+1) = IP(J)
          DO 20 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (I.GT.NR .OR. I.LT.LOWER) THEN
              IOUT = IOUT + 1
              IF (ICNTL6.NE.0 .AND. I.LT.J) IUP = IUP + 1
            ELSE IF (IW(I).LE.NZJ) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              A(KNE) = A(K)
              IP(J+1) = IP(J+1) + 1
              IW(I) = KNE
            ELSE
              IDUP = IDUP + 1
              IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
            END IF
   20     CONTINUE
          KSTART = KSTOP
          NZJ = KNE
   30   CONTINUE
      ELSE
        DO 50 J = 1,NC
          IF (ICNTL6.NE.0) LOWER = J
          KSTOP = IP(J+1)
          IF (KSTART.GT.KSTOP) THEN
            INFO(1) = -9
            INFO(2) = J
            RETURN
          END IF
          IP(J+1) = IP(J)
          DO  40 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (I.GT.NR .OR. I.LT.LOWER) THEN
              IOUT = IOUT + 1
              IF (ICNTL6.NE.0 .AND. I.GT.1) IUP = IUP + 1
            ELSE IF (IW(I).LT.J) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              IP(J+1) = IP(J+1) + 1
              IW(I) = J
            ELSE
              IDUP = IDUP + 1
            END IF
   40     CONTINUE
          KSTART = KSTOP
   50   CONTINUE
      END IF
      RETURN
      END
* COPYRIGHT (c) 1982 AEA Technology
*######DATE 20 September 2001
C  September 2001: threadsafe version of MA27
C  19/3/03. Array ICNTL in MA27GD made assumed size. 

      SUBROUTINE MA27ID(ICNTL,CNTL)
      INTEGER ICNTL(30)
      DOUBLE PRECISION CNTL(5)

      INTEGER IFRLVL
      PARAMETER ( IFRLVL=5 )

C Stream number for error messages
      ICNTL(1) = 6
C Stream number for diagnostic messages
      ICNTL(2) = 6
C Control the level of diagnostic printing.
C   0 no printing
C   1 printing of scalar parameters and first parts of arrays.
C   2 printing of scalar parameters and whole of arrays.
      ICNTL(3) = 0
C The largest integer such that all integers I in the range
C -ICNTL(4).LE.I.LE.ICNTL(4) can be handled by the shortest integer
C type in use.
      ICNTL(4) = 2139062143
C Minimum number of eliminations in a step that is automatically
C accepted. if two adjacent steps can be combined and each has less
C eliminations then they are combined.
      ICNTL(5) = 1
C Control whether direct or indirect access is used by MA27C/CD.
C Indirect access is employed in forward and back substitution 
C respectively if the size of a block is less than
C ICNTL(IFRLVL+MIN(10,NPIV)) and ICNTL(IFRLVL+10+MIN(10,NPIV))
C respectively, where NPIV is the number of pivots in the block.
      ICNTL(IFRLVL+1)  = 32639
      ICNTL(IFRLVL+2)  = 32639
      ICNTL(IFRLVL+3)  = 32639
      ICNTL(IFRLVL+4)  = 32639
      ICNTL(IFRLVL+5)  = 14
      ICNTL(IFRLVL+6)  = 9
      ICNTL(IFRLVL+7)  = 8
      ICNTL(IFRLVL+8)  = 8
      ICNTL(IFRLVL+9)  = 9
      ICNTL(IFRLVL+10) = 10
      ICNTL(IFRLVL+11) = 32639
      ICNTL(IFRLVL+12) = 32639
      ICNTL(IFRLVL+13) = 32639
      ICNTL(IFRLVL+14) = 32689
      ICNTL(IFRLVL+15) = 24
      ICNTL(IFRLVL+16) = 11
      ICNTL(IFRLVL+17) = 9
      ICNTL(IFRLVL+18) = 8
      ICNTL(IFRLVL+19) = 9
      ICNTL(IFRLVL+20) = 10
C Not used
      ICNTL(26) = 0
      ICNTL(27) = 0
      ICNTL(28) = 0
      ICNTL(29) = 0
      ICNTL(30) = 0

C Control threshold pivoting.
      CNTL(1) = 0.1D0
C If a column of the reduced matrix has relative density greater than
C CNTL(2), it is forced into the root. All such columns are taken to
C have sparsity pattern equal to their merged patterns, so the fill
C and operation counts may be overestimated.
      CNTL(2) = 1.0D0
C An entry with absolute value less than CNTL(3) is never accepted as
C a 1x1 pivot or as the off-diagonal of a 2x2 pivot.
      CNTL(3) = 0.0D0
C Not used
      CNTL(4) = 0.0
      CNTL(5) = 0.0

      RETURN
      END

      SUBROUTINE MA27AD(N,NZ,IRN,ICN,IW,LIW,IKEEP,IW1,NSTEPS,IFLAG,
     +                 ICNTL,CNTL,INFO,OPS)
C THIS SUBROUTINE COMPUTES A MINIMUM DEGREE ORDERING OR ACCEPTS A GIVEN
C     ORDERING. IT COMPUTES ASSOCIATED ASSEMBLY AND ELIMINATION
C     INFORMATION FOR MA27B/BD.
C N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
C NZ MUST BE SET TO THE NUMBER OF NON-ZEROS INPUT. IT IS NOT
C     ALTERED.
C IRN(I),I=1,2,...,NZ MUST BE SET TO THE ROW NUMBERS OF THE
C     NON-ZEROS. IT IS NOT ALTERED UNLESS IT IS EQUIVALENCED
C     TO IW (SEE DESCRIPTION OF IW).
C ICN(I),I=1,2,...,NZ MUST BE SET TO THE COLUMN NUMBERS OF THE
C     NON-ZEROS. IT IS NOT ALTERED UNLESS IT IS EQUIVALENCED
C     TO IW (SEE DESCRIPTION OF IW).
C IW NEED NOT BE SET ON INPUT. IT IS USED FOR WORKSPACE.
C     IRN(1) MAY BE EQUIVALENCED TO IW(1) AND ICN(1) MAY BE
C     EQUIVALENCED TO IW(K), WHERE K.GT.NZ.
C LIW MUST BE SET TO THE LENGTH OF IW. IT MUST BE AT LEAST 2*NZ+3*N
C     FOR THE IFLAG=0 ENTRY AND AT LEAST NZ+3*N FOR THE IFLAG=1
C     ENTRY. IT IS NOT ALTERED.
C IKEEP NEED NOT BE SET UNLESS AN ORDERING IS GIVEN, IN WHICH CASE
C     IKEEP(I,1) MUST BE SET TO THE POSITION OF VARIABLE I IN THE
C     ORDER. ON OUTPUT IKEEP CONTAINS INFORMATION NEEDED BY MA27B/BD.
C     IKEEP(I,1) HOLDS THE POSITION OF VARIABLE I IN THE PIVOT ORDER.
C     IKEEP(I,2), IKEEP(I,3) HOLD THE NUMBER OF ELIMINATIONS, ASSEMBLIES
C     AT MAJOR STEP I, I=1,2,...,NSTEPS. NOTE THAT WHEN AN ORDER IS
C     GIVEN IT MAY BE REPLACED BY ANOTHER ORDER THAT GIVES IDENTICAL
C     NUMERICAL RESULTS.
C IW1 IS USED FOR WORKSPACE.
C NSTEPS NEED NOT BE SET. ON OUTPUT IT CONTAINS THE NUMBER OF MAJOR
C     STEPS NEEDED FOR A DEFINITE MATRIX AND MUST BE PASSED UNCHANGED
C     TO MA27B/BD.
C IFLAG MUST SET TO ZERO IF THE USER WANTS THE PIVOT ORDER CHOSEN
C     AUTOMATICALLY AND TO ONE IF HE WANTS TO SPECIFY IT IN IKEEP.
C ICNTL is an INTEGER array of length 30 containing user options
C     with integer values (defaults set in MA27I/ID)
C   ICNTL(1) (LP) MUST BE SET TO THE STREAM NUMBER FOR ERROR MESSAGES.
C     ERROR MESSAGES ARE SUPPRESSED IF ICNTL(1) IS NOT POSITIVE.
C     IT IS NOT ALTERED.
C   ICNTL(2) (MP) MUST BE SET TO THE STREAM NUMBER FOR DIAGNOSTIC
C     MESSAGES.  DIAGNOSTIC MESSAGES ARE SUPPRESSED IF ICNTL(2) IS NOT
C     POSITIVE.  IT IS NOT ALTERED.
C   ICNTL(3) (LDIAG) CONTROLS THE LEVEL OF DIAGNOSTIC PRINTING.
C     0 NO PRINTING
C     1 PRINTING OF SCALAR PARAMETERS AND FIRST PARTS OF ARRAYS.
C     2 PRINTING OF SCALAR PARAMETERS AND WHOLE OF ARRAYS.
C   ICNTL(4) (IOVFLO) IS THE LARGEST INTEGER SUCH THAT ALL INTEGERS
C     I IN THE RANGE -IOVFLO.LE.I.LE.IOVFLO CAN BE HANDLED BY THE
C     SHORTEST INTEGER TYPE IN USE.
C   ICNT(5) (NEMIN) MUST BE SET TO THE MINIMUM NUMBER OF ELIMINATIONS
C     IN A STEP THAT IS AUTOMATICALLY ACCEPTED. IF TWO ADJACENT STEPS
C     CAN BE COMBINED AND EACH HAS LESS ELIMINATIONS THEN THEY ARE
C     COMBINED.
C   ICNTL(IFRLVL+I) I=1,20, (IFRLVL) MUST BE SET TO CONTROL WHETHER
C     DIRECT OR INDIRECT ACCESS IS USED BY MA27C/CD.  INDIRECT ACCESS
C     IS EMPLOYED IN FORWARD AND BACK SUBSTITUTION RESPECTIVELY IF THE
C     SIZE OF A BLOCK IS LESS THAN ICNTL(IFRLVL+(MIN(10,NPIV)) AND
C     ICNTL(IFRLVL+10+MIN(10,NPIV)) RESPECTIVELY, WHERE NPIV IS THE
C     NUMBER OF PIVOTS IN THE BLOCK.
C   ICNTL(I) I=26,30 are not used.
C CNTL is an DOUBLE PRECISION array of length 5 containing user options
C     with real values (defaults set in MA27I/ID)
C   CNTL(1) (U) IS USED TO CONTROL THRESHOLD PIVOTING. IT IS NOT
C     ALTERED.
C   CNTL(2) (FRATIO) has default value 1.0.  If a column of the
C      reduced matrix has relative density greater than CNTL(2), it
C      is forced into the root. All such columns are taken to have
C      sparsity pattern equal to their merged patterns, so the fill
C      and operation counts may be overestimated.
C   CNTL(3) (PIVTOL) has default value 0.0. An entry with absolute
C      value less than CNTL(3) is never accepted as a 1x1 pivot or
C      as the off-diagonal of a 2x2 pivot.
C   CNTL(I) I=4,5 are not used.
C INFO is an INTEGER array of length 20 which is used to return
C     information to the user.
C   INFO(1) (IFLAG) is an error return code, zero for success, greater
C      than zero for a warning and less than zero for errors, see
C      INFO(2).
C   INFO(2) (IERROR) HOLDS ADDITIONAL INFORMATION IN THE EVENT OF ERRORS.
C     IF INFO(1)=-3 INFO(2) HOLDS A LENGTH THAT MAY SUFFICE FOR IW.
C     IF INFO(1)=-4 INFO(2) HOLDS A LENGTH THAT MAY SUFFICE FOR A.
C     IF INFO(1)=-5 INFO(2) IS SET TO THE PIVOT STEP AT WHICH SINGULARITY
C                 WAS DETECTED.
C     IF INFO(1)=-6 INFO(2) IS SET TO THE PIVOT STEP AT WHICH A CHANGE OF
C                 PIVOT SIGN WAS FOUND.
C     IF INFO(1)= 1 INFO(2) HOLDS THE NUMBER OF FAULTY ENTRIES.
C     IF INFO(1)= 2 INFO(2) IS SET TO THE NUMBER OF SIGNS CHANGES IN
C                 THE PIVOTS.
C     IF INFO(1)=3 INFO(2) IS SET TO THE RANK OF THE MATRIX.
C   INFO(3) and INFO(4) (NRLTOT and NIRTOT) REAL AND INTEGER STRORAGE
C     RESPECTIVELY REQUIRED FOR THE FACTORIZATION IF NO COMPRESSES ARE
C     ALLOWED.
C   INFO(5) and INFO(6) (NRLNEC and NIRNEC) REAL AND INTEGER STORAGE
C     RESPECTIVELY REQUIRED FOR THE FACTORIZATION IF COMPRESSES ARE
C     ALLOWED AND THE MATRIX IS DEFINITE.
C   INFO(7) and INFO(8) (NRLADU and NIRADU) REAL AND INTEGER STORAGE
C     RESPECTIVELY FOR THE MATRIX FACTORS AS CALCULATED BY MA27A/AD
C     FOR THE DEFINITE CASE.  
C   INFO(9) and INFO(10) (NRLBDU and NIRBDU) REAL AND INTEGER STORAGE
C     RESPECTIVELY FOR THE MATRIX FACTORS AS FOUND  BY MA27B/BD.
C   INFO(11) (NCMPA) ACCUMULATES THE NUMBER OF TIMES THE ARRAY IW IS
C     COMPRESSED BY MA27A/AD.
C   INFO(12) and INFO(13) (NCMPBR and NCMPBI) ACCUMULATE THE NUMBER
C     OF COMPRESSES OF THE REALS AND INTEGERS PERFORMED BY MA27B/BD.
C   INFO(14) (NTWO) IS USED BY MA27B/BD TO RECORD THE NUMBER OF 2*2
C     PIVOTS USED.
C   INFO(15) (NEIG) IS USED BY ME27B/BD TO RECORD THE NUMBER OF
C     NEGATIVE EIGENVALUES OF A.
C   INFO(16) to INFO(20) are not used. 
C OPS ACCUMULATES THE NO. OF MULTIPLY/ADD PAIRS NEEDED TO CREATE THE
C     TRIANGULAR FACTORIZATION, IN THE DEFINITE CASE.
C
C     .. Scalar Arguments ..
      INTEGER IFLAG,LIW,N,NSTEPS,NZ
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CNTL(5),OPS
      INTEGER ICNTL(30),INFO(20)
      INTEGER ICN(*),IKEEP(N,3),IRN(*),IW(LIW),IW1(N,2)
C     ..
C     .. Local Scalars ..
      INTEGER I,IWFR,K,L1,L2,LLIW
C     ..
C     .. External Subroutines ..
      EXTERNAL MA27GD,MA27HD,MA27JD,MA27KD,MA27LD,MA27MD,MA27UD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
C     .. Executable Statements ..
      DO 5 I = 1,15
        INFO(I) = 0
    5 CONTINUE

      IF (ICNTL(3).LE.0 .OR. ICNTL(2).LE.0) GO TO 40
C PRINT INPUT VARIABLES.
      WRITE (ICNTL(2),FMT=10) N,NZ,LIW,IFLAG

   10 FORMAT(/,/,' ENTERING MA27AD WITH      N     NZ      LIW  IFLAG',
     +       /,21X,I7,I7,I9,I7)

      NSTEPS = 0
      K = MIN(8,NZ)
      IF (ICNTL(3).GT.1) K = NZ
      IF (K.GT.0) WRITE (ICNTL(2),FMT=20) (IRN(I),ICN(I),I=1,K)

   20 FORMAT (' MATRIX NON-ZEROS',/,4 (I9,I6),/,
     +       (I9,I6,I9,I6,I9,I6,I9,I6))

      K = MIN(10,N)
      IF (ICNTL(3).GT.1) K = N
      IF (IFLAG.EQ.1 .AND. K.GT.0) THEN
        WRITE (ICNTL(2),FMT=30) (IKEEP(I,1),I=1,K)
      END IF

   30 FORMAT (' IKEEP(.,1)=',10I6,/, (12X,10I6))

   40 IF (N.LT.1 .OR. N.GT.ICNTL(4)) GO TO 70
      IF (NZ.LT.0) GO TO 100
      LLIW = LIW - 2*N
      L1 = LLIW + 1
      L2 = L1 + N
      IF (IFLAG.EQ.1) GO TO 50
      IF (LIW.LT.2*NZ+3*N+1) GO TO 130
C SORT
      CALL MA27GD(N,NZ,IRN,ICN,IW,LLIW,IW1,IW1(1,2),IW(L1),IWFR,
     +           ICNTL,INFO)
C ANALYZE USING MINIMUM DEGREE ORDERING
      CALL MA27HD(N,IW1,IW,LLIW,IWFR,IW(L1),IW(L2),IKEEP(1,2),
     +           IKEEP(1,3),IKEEP,ICNTL(4),INFO(11),CNTL(2))
      GO TO 60
C SORT USING GIVEN ORDER
   50 IF (LIW.LT.NZ+3*N+1) GO TO 120
      CALL MA27JD(N,NZ,IRN,ICN,IKEEP,IW,LLIW,IW1,IW1(1,2),IW(L1),IWFR,
     +           ICNTL,INFO)
C ANALYZE USING GIVEN ORDER
      CALL MA27KD(N,IW1,IW,LLIW,IWFR,IKEEP,IKEEP(1,2),IW(L1),IW(L2),
     +           INFO(11))
C PERFORM DEPTH-FIRST SEARCH OF ASSEMBLY TREE
   60 CALL MA27LD(N,IW1,IW(L1),IKEEP,IKEEP(1,2),IKEEP(1,3),IW(L2),
     +           NSTEPS,ICNTL(5))
C EVALUATE STORAGE AND OPERATION COUNT REQUIRED BY MA27B/BD IN THE
C     DEFINITE CASE.
C SET IW(1) SO THAT ARRAYS IW AND IRN CAN BE TESTED FOR EQUIVALENCE
C     IN MA27M/MD.
      IF(NZ.GE.1) IW(1) = IRN(1) + 1
      CALL MA27MD(N,NZ,IRN,ICN,IKEEP,IKEEP(1,3),IKEEP(1,2),IW(L2),
     +           NSTEPS,IW1,IW1(1,2),IW,INFO,OPS)
      GO TO 160

   70 INFO(1) = -1
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)

   80 FORMAT (' **** ERROR RETURN FROM MA27AD **** INFO(1)=',I3)

      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=90) N

   90 FORMAT (' VALUE OF N OUT OF RANGE ... =',I10)

      GO TO 160

  100 INFO(1) = -2
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=110) NZ

  110 FORMAT (' VALUE OF NZ OUT OF RANGE .. =',I10)

      GO TO 160

  120 INFO(2) = NZ + 3*N + 1
      GO TO 140

  130 INFO(2) = 2*NZ + 3*N + 1
  140 INFO(1) = -3
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=150) LIW,INFO(2)

  150 FORMAT (' LIW TOO SMALL, MUST BE INCREASED FROM',I10,
     +       ' TO AT LEAST',I10)

  160 IF (ICNTL(3).LE.0 .OR. ICNTL(2).LE.0 .OR. INFO(1).LT.0) GO TO 200
C PRINT PARAMETER VALUES ON EXIT.
      WRITE (ICNTL(2),FMT=170) NSTEPS,INFO(1),OPS,INFO(2),INFO(3),
     +  INFO(4),INFO(5),INFO(6),INFO(7),INFO(8),INFO(11)

  170 FORMAT (/,' LEAVING MA27AD WITH NSTEPS  INFO(1)    OPS IERROR',
     +          ' NRLTOT NIRTOT',
     +        /,20X,2I7,F7.0,3I7,
     +        /,20X,' NRLNEC NIRNEC NRLADU NIRADU  NCMPA',
     +        /,20X,6I7)

      K = MIN(9,N)
      IF (ICNTL(3).GT.1) K = N
      IF (K.GT.0) WRITE (ICNTL(2),FMT=30) (IKEEP(I,1),I=1,K)
      K = MIN(K,NSTEPS)
      IF (K.GT.0) WRITE (ICNTL(2),FMT=180) (IKEEP(I,2),I=1,K)

  180 FORMAT (' IKEEP(.,2)=',10I6,/, (12X,10I6))

      IF (K.GT.0) WRITE (ICNTL(2),FMT=190) (IKEEP(I,3),I=1,K)

  190 FORMAT (' IKEEP(.,3)=',10I6,/, (12X,10I6))

  200 CONTINUE

      RETURN
      END

      SUBROUTINE MA27BD(N,NZ,IRN,ICN,A,LA,IW,LIW,IKEEP,NSTEPS,MAXFRT,
     +                 IW1,ICNTL,CNTL,INFO)
C THIS SUBROUTINE COMPUTES THE FACTORISATION OF THE MATRIX INPUT IN
C     A,IRN,ICN USING INFORMATION (IN IKEEP) FROM MA27A/AD.
C N MUST BE SET TO THE ORDER OF THE MATRIX.  IT IS NOT ALTERED.
C NZ MUST BE SET TO THE NUMBER OF NON-ZEROS INPUT.  IT IS NOT
C     ALTERED.
C IRN,ICN,A.  ENTRY K (K=1,NZ) OF IRN,ICN MUST BE SET TO THE ROW
C     AND COLUMN INDEX RESPECTIVELY OF THE NON-ZERO IN A.
C     IRN AND ICN ARE UNALTERED BY MA27B/BD.
C     ON EXIT, ENTRIES 1 TO NRLBDU OF A HOLD REAL INFORMATION
C     ON THE FACTORS AND SHOULD BE PASSED UNCHANGED TO MA27C/CD.
C LA LENGTH OF ARRAY A.  AN INDICATION OF A SUITABLE VALUE,
C     SUFFICIENT FOR FACTORIZATION OF A DEFINITE MATRIX, WILL
C     HAVE BEEN PROVIDED IN NRLNEC AND NRLTOT BY MA27A/AD.
C     IT IS NOT ALTERED BY MA27B/BD.
C IW NEED NOT BE SET ON ENTRY.  USED AS A WORK ARRAY BY MA27B/BD.
C     ON EXIT, ENTRIES 1 TO NIRBDU HOLD INTEGER INFORMATION ON THE
C     FACTORS AND SHOULD BE PASSED UNCHANGED TO MA27C/CD.
C LIW LENGTH OF ARRAY IW.  AN INDICATION OF A SUITABLE VALUE WILL
C     HAVE BEEN PROVIDED IN NIRNEC AND NIRTOT BY MA27A/AD.
C     IT IS NOT ALTERED BY MA27B/BD.
C IKEEP MUST BE UNCHANGED SINCE THE CALL TO MA27A/AD.
C     IT IS NOT ALTERED BY MA27B/BD.
C NSTEPS MUST BE UNCHANGED SINCE THE CALL TO MA27A/AD.
C     IT IS NOT ALTERED BY MA27B/BD.
C MAXFRT NEED NOT BE SET AND MUST BE PASSED UNCHANGED TO MA27C/CD.
C     IT IS THE MAXIMUM SIZE OF THE FRONTAL MATRIX GENERATED DURING
C     THE DECOMPOSITION.
C IW1 USED AS WORKSPACE BY MA27B/BD.
C ICNTL is an INTEGER array of length 30, see MA27A/AD.
C CNTL is a DOUBLE PRECISION array of length 5, see MA27A/AD.
C INFO is an INTEGER array of length 20, see MA27A/AD.
C
C     .. Scalar Arguments ..
      INTEGER LA,LIW,MAXFRT,N,NSTEPS,NZ
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA),CNTL(5)
      INTEGER ICN(*),IKEEP(N,3),IRN(*),IW(LIW),IW1(N)
      INTEGER ICNTL(30),INFO(20)
C     ..
C     .. Local Scalars ..
      INTEGER I,IAPOS,IBLK,IPOS,IROWS,J1,J2,JJ,K,KBLK,KZ,LEN,NCOLS,
     +        NROWS,NZ1
C     ..
C     .. External Subroutines ..
      EXTERNAL MA27ND,MA27OD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN
C     ..
C     .. Executable Statements ..
      INFO(1) = 0

      IF (ICNTL(3).LE.0 .OR. ICNTL(2).LE.0) GO TO 60
C PRINT INPUT PARAMETERS.
      WRITE (ICNTL(2),FMT=10) N,NZ,LA,LIW,NSTEPS,CNTL(1)

   10 FORMAT (/,/,
     + ' ENTERING MA27BD WITH      N     NZ       LA      LIW',
     +       ' NSTEPS      U',/,21X,I7,I7,I9,I9,I7,1PD10.2)

      KZ = MIN(6,NZ)
      IF (ICNTL(3).GT.1) KZ = NZ
      IF (NZ.GT.0) WRITE (ICNTL(2),FMT=20) (A(K),IRN(K),ICN(K),K=1,KZ)

   20 FORMAT (' MATRIX NON-ZEROS',/,1X,2 (1P,D16.3,2I6),/,
     +       (1X,1P,D16.3,2I6,1P,D16.3,2I6))

      K = MIN(9,N)
      IF (ICNTL(3).GT.1) K = N
      IF (K.GT.0) WRITE (ICNTL(2),FMT=30) (IKEEP(I,1),I=1,K)

   30 FORMAT (' IKEEP(.,1)=',10I6,/, (12X,10I6))

      K = MIN(K,NSTEPS)
      IF (K.GT.0) WRITE (ICNTL(2),FMT=40) (IKEEP(I,2),I=1,K)

   40 FORMAT (' IKEEP(.,2)=',10I6,/, (12X,10I6))

      IF (K.GT.0) WRITE (ICNTL(2),FMT=50) (IKEEP(I,3),I=1,K)

   50 FORMAT (' IKEEP(.,3)=',10I6,/, (12X,10I6))

   60 IF (N.LT.1 .OR. N.GT.ICNTL(4)) GO TO 70
      IF (NZ.LT.0) GO TO 100
      IF (LIW.LT.NZ) GO TO 120
      IF (LA.LT.NZ+N) GO TO 150
      IF (NSTEPS.LT.1 .OR. NSTEPS.GT.N) GO TO 175
C SORT
      CALL MA27ND(N,NZ,NZ1,A,LA,IRN,ICN,IW,LIW,IKEEP,IW1,ICNTL,INFO)
      IF (INFO(1).EQ.-3) GO TO 130
      IF (INFO(1).EQ.-4) GO TO 160
C FACTORIZE
      CALL MA27OD(N,NZ1,A,LA,IW,LIW,IKEEP,IKEEP(1,3),NSTEPS,MAXFRT,
     +           IKEEP(1,2),IW1,ICNTL,CNTL,INFO)
      IF (INFO(1).EQ.-3) GO TO 130
      IF (INFO(1).EQ.-4) GO TO 160
      IF (INFO(1).EQ.-5) GO TO 180
      IF (INFO(1).EQ.-6) GO TO 200
C **** WARNING MESSAGE ****
      IF (INFO(1).EQ.3 .AND. ICNTL(2).GT.0) THEN
        WRITE (ICNTL(2),FMT=65) INFO(1),INFO(2)
      END IF

   65 FORMAT (' *** WARNING MESSAGE FROM SUBROUTINE MA27BD',
     +        '  *** INFO(1) =',I2,
     +        /,5X,'MATRIX IS SINGULAR. RANK=',I5)

      GO TO 220
C **** ERROR RETURNS ****
   70 INFO(1) = -1
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)

   80 FORMAT (' **** ERROR RETURN FROM MA27BD **** INFO(1)=',I3)

      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=90) N

   90 FORMAT (' VALUE OF N OUT OF RANGE ... =',I10)

      GO TO 220

  100 INFO(1) = -2
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=110) NZ

  110 FORMAT (' VALUE OF NZ OUT OF RANGE .. =',I10)

      GO TO 220

  120 INFO(1) = -3
      INFO(2) = NZ
  130 IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=140) LIW,INFO(2)

  140 FORMAT (' LIW TOO SMALL, MUST BE INCREASED FROM',I10,' TO',
     +       ' AT LEAST',I10)

      GO TO 220

  150 INFO(1) = -4
      INFO(2) = NZ + N
  160 IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=170) LA,INFO(2)

  170 FORMAT (' LA TOO SMALL, MUST BE INCREASED FROM ',I10,' TO',
     +       ' AT LEAST',I10)

      GO TO 220

  175 INFO(1) = -7
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) THEN
        WRITE (ICNTL(1),FMT='(A)') ' NSTEPS is out of range'
      END IF
      GO TO 220

  180 IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=190) INFO(2)

  190 FORMAT (' ZERO PIVOT AT STAGE',I10,
     +        ' WHEN INPUT MATRIX DECLARED DEFINITE')

      GO TO 220

  200 IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=80) INFO(1)
      IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=210)

  210 FORMAT (' CHANGE IN SIGN OF PIVOT ENCOUNTERED',
     +        ' WHEN FACTORING ALLEGEDLY DEFINITE MATRIX')

  220 IF (ICNTL(3).LE.0 .OR. ICNTL(2).LE.0 .OR. INFO(1).LT.0) GO TO 310
C PRINT OUTPUT PARAMETERS.
      WRITE (ICNTL(2),FMT=230) MAXFRT,INFO(1),INFO(9),INFO(10),INFO(12),
     +  INFO(13),INFO(14),INFO(2)

  230 FORMAT (/,' LEAVING MA27BD WITH',
     +        /,10X,'  MAXFRT  INFO(1) NRLBDU NIRBDU NCMPBR',
     +         ' NCMPBI   NTWO IERROR',
     +        /,11X,8I7)

      IF (INFO(1).LT.0) GO TO 310
C PRINT OUT MATRIX FACTORS FROM MA27B/BD.
      KBLK = ABS(IW(1)+0)
      IF (KBLK.EQ.0) GO TO 310
      IF (ICNTL(3).EQ.1) KBLK = 1
      IPOS = 2
      IAPOS = 1
      DO 300 IBLK = 1,KBLK
        NCOLS = IW(IPOS)
        NROWS = IW(IPOS+1)
        J1 = IPOS + 2
        IF (NCOLS.GT.0) GO TO 240
        NCOLS = -NCOLS
        NROWS = 1
        J1 = J1 - 1
  240   WRITE (ICNTL(2),FMT=250) IBLK,NROWS,NCOLS

  250   FORMAT (' BLOCK PIVOT =',I8,' NROWS =',I8,' NCOLS =',I8)

        J2 = J1 + NCOLS - 1
        IPOS = J2 + 1
        WRITE (ICNTL(2),FMT=260) (IW(JJ),JJ=J1,J2)

  260   FORMAT (' COLUMN INDICES =',10I6,/, (17X,10I6))

        WRITE (ICNTL(2),FMT=270)

  270   FORMAT (' REAL ENTRIES .. EACH ROW STARTS ON A NEW LINE')

        LEN = NCOLS
        DO 290 IROWS = 1,NROWS
          J1 = IAPOS
          J2 = IAPOS + LEN - 1
          WRITE (ICNTL(2),FMT=280) (A(JJ),JJ=J1,J2)

  280     FORMAT (1P,5D13.3)

          LEN = LEN - 1
          IAPOS = J2 + 1
  290   CONTINUE
  300 CONTINUE
  310 RETURN
      END

      SUBROUTINE MA27CD(N,A,LA,IW,LIW,W,MAXFRT,RHS,IW1,NSTEPS,
     + ICNTL,INFO)
C THIS SUBROUTINE USES THE FACTORISATION OF THE MATRIX IN A,IW TO
C     SOLVE A SYSTEM OF EQUATIONS.
C N MUST BE SET TO THE ORDER OF THE MATRIX.  IT IS NOT ALTERED.
C A,IW HOLD INFORMATION ON THE FACTORS AND MUST BE UNCHANGED SINCE
C     THE CALL TO MA27B/BD. THEY ARE NOT ALTERED BY MA27C/CDD.
C LA,LIW MUST BE SET TO THE LENGTHS OF A,IW RESPECTIVELY.  THEY
C     ARE NOT ALTERED.
C W USED AS A WORK ARRAY.
C MAXFRT IS THE LENGTH OF W AND MUST BE PASSED UNCHANGED FROM THE
C     CALL TO MA27B/BD.  IT IS NOT ALTERED.
C RHS MUST BE SET TO THE RIGHT HAND SIDE FOR THE EQUATIONS BEING
C     SOLVED.  ON EXIT, THIS ARRAY WILL HOLD THE SOLUTION.
C IW1 USED AS A WORK ARRAY.
C NSTEPS IS THE LENGTH OF IW1 WHICH MUST BE AT LEAST THE ABSOLUTE
C     VALUE OF IW(1).  IT IS NOT ALTERED.
C ICNTL is an INTEGER array of length 30, see MA27A/AD.
C INFO is an INTEGER array of length 20, see MA27A/AD.
C
C     .. Scalar Arguments ..
      INTEGER LA,LIW,MAXFRT,N,NSTEPS
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA),RHS(N),W(MAXFRT)
      INTEGER IW(LIW),IW1(NSTEPS),ICNTL(30),INFO(20)
C     ..
C     .. Local Scalars ..
      INTEGER I,IAPOS,IBLK,IPOS,IROWS,J1,J2,JJ,K,KBLK,LATOP,LEN,NBLK,
     +        NCOLS,NROWS
C     ..
C     .. External Subroutines ..
      EXTERNAL MA27QD,MA27RD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN
C     ..
C     .. Executable Statements ..
      INFO(1) = 0

      IF (ICNTL(3).LE.0 .OR. ICNTL(2).LE.0) GO TO 110
C PRINT INPUT PARAMETERS.
      WRITE (ICNTL(2),FMT=10) N,LA,LIW,MAXFRT,NSTEPS

   10 FORMAT (/,/,' ENTERING MA27CD WITH      N     LA    LIW MAXFRT',
     +       '  NSTEPS',/,21X,5I7)
C PRINT OUT MATRIX FACTORS FROM MA27B/BD.
      KBLK = ABS(IW(1)+0)
      IF (KBLK.EQ.0) GO TO 90
      IF (ICNTL(3).EQ.1) KBLK = 1
      IPOS = 2
      IAPOS = 1
      DO 80 IBLK = 1,KBLK
        NCOLS = IW(IPOS)
        NROWS = IW(IPOS+1)
        J1 = IPOS + 2
        IF (NCOLS.GT.0) GO TO 20
        NCOLS = -NCOLS
        NROWS = 1
        J1 = J1 - 1
   20   WRITE (ICNTL(2),FMT=30) IBLK,NROWS,NCOLS

   30   FORMAT (' BLOCK PIVOT =',I8,' NROWS =',I8,' NCOLS =',I8)

        J2 = J1 + NCOLS - 1
        IPOS = J2 + 1
        WRITE (ICNTL(2),FMT=40) (IW(JJ),JJ=J1,J2)

   40   FORMAT (' COLUMN INDICES =',10I6,/, (17X,10I6))

        WRITE (ICNTL(2),FMT=50)

   50   FORMAT (' REAL ENTRIES .. EACH ROW STARTS ON A NEW LINE')

        LEN = NCOLS
        DO 70 IROWS = 1,NROWS
          J1 = IAPOS
          J2 = IAPOS + LEN - 1
          WRITE (ICNTL(2),FMT=60) (A(JJ),JJ=J1,J2)

   60     FORMAT (1P,5D13.3)

          LEN = LEN - 1
          IAPOS = J2 + 1
   70   CONTINUE
   80 CONTINUE
   90 K = MIN(10,N)
      IF (ICNTL(3).GT.1) K = N
      IF (N.GT.0) WRITE (ICNTL(2),FMT=100) (RHS(I),I=1,K)

  100 FORMAT (' RHS',1P,5D13.3,/, (4X,1P,5D13.3))

  110 IF (IW(1).LT.0) GO TO 130
      NBLK = IW(1)
      IF (NBLK.GT.0) GO TO 140
C WE HAVE A ZERO MATRIX
      DO 120 I = 1,N
        RHS(I) = 0.0D0
  120 CONTINUE
      GO TO 150

  130 NBLK = -IW(1)
C FORWARD SUBSTITUTION
  140 CALL MA27QD(N,A,LA,IW(2),LIW-1,W,MAXFRT,RHS,IW1,NBLK,LATOP,ICNTL)
C BACK SUBSTITUTION.
      CALL MA27RD(N,A,LA,IW(2),LIW-1,W,MAXFRT,RHS,IW1,NBLK,LATOP,ICNTL)
  150 IF (ICNTL(3).LE.0 .OR. ICNTL(2).LE.0) GO TO 170
C PRINT OUTPUT PARAMETERS.
      WRITE (ICNTL(2),FMT=160)

  160 FORMAT (/,/,' LEAVING MA27CD WITH')

      IF (N.GT.0) WRITE (ICNTL(2),FMT=100) (RHS(I),I=1,K)
  170 CONTINUE

      RETURN
      END
      SUBROUTINE MA27GD(N,NZ,IRN,ICN,IW,LW,IPE,IQ,FLAG,IWFR,
     +                 ICNTL,INFO)
C
C SORT PRIOR TO CALLING ANALYSIS ROUTINE MA27H/HD.
C
C GIVEN THE POSITIONS OF THE OFF-DIAGONAL NON-ZEROS OF A SYMMETRIC
C     MATRIX, CONSTRUCT THE SPARSITY PATTERN OF THE OFF-DIAGONAL
C     PART OF THE WHOLE MATRIX (UPPER AND LOWER TRIANGULAR PARTS).
C EITHER ONE OF A PAIR (I,J),(J,I) MAY BE USED TO REPRESENT
C     THE PAIR. DIAGONAL ELEMENTS AND DUPLICATES ARE IGNORED.
C
C N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
C NZ MUST BE SET TO THE NUMBER OF NON-ZEROS INPUT. IT IS NOT
C     ALTERED.
C IRN(I),I=1,2,...,NZ MUST BE SET TO THE ROW NUMBERS OF THE
C     NON-ZEROS ON INPUT. IT IS NOT ALTERED UNLESS IT IS EQUIVALENCED
C     TO IW (SEE DESCRIPTION OF IW).
C ICN(I),I=1,2,...,NZ MUST BE SET TO THE COLUMN NUMBERS OF THE
C     NON-ZEROS ON INPUT. IT IS NOT ALTERED UNLESS IT IS EQUIVALENCED
C     TO IW (SEE DESCRIPTION OF IW).
C IW NEED NOT BE SET ON INPUT. ON OUTPUT IT CONTAINS LISTS OF
C     COLUMN INDICES, EACH LIST BEING HEADED BY ITS LENGTH.
C     IRN(1) MAY BE EQUIVALENCED TO IW(1) AND ICN(1) MAY BE
C     EQUIVALENCED TO IW(K), WHERE K.GT.NZ.
C LW MUST BE SET TO THE LENGTH OF IW. IT MUST BE AT LEAST 2*NZ+N.
C     IT IS NOT ALTERED.
C IPE NEED NOT BE SET ON INPUT. ON OUTPUT IPE(I) POINTS TO THE START OF
C     THE ENTRY IN IW FOR ROW I, OR IS ZERO IF THERE IS NO ENTRY.
C IQ NEED NOT BE SET.  ON OUTPUT IQ(I),I=1,N CONTAINS THE NUMBER OF
C     OFF-DIAGONAL N0N-ZEROS IN ROW I INCLUDING DUPLICATES.
C FLAG IS USED FOR WORKSPACE TO HOLD FLAGS TO PERMIT DUPLICATE ENTRIES
C     TO BE IDENTIFIED QUICKLY.
C IWFR NEED NOT BE SET ON INPUT. ON OUTPUT IT POINTS TO THE FIRST
C     UNUSED LOCATION IN IW.
C ICNTL is an INTEGER array of length 30, see MA27A/AD.
C INFO is an INTEGER array of length 20, see MA27A/AD.
C
C     .. Scalar Arguments ..
      INTEGER IWFR,LW,N,NZ
C     ..
C     .. Array Arguments ..
      INTEGER FLAG(N),ICN(*),IPE(N),IQ(N),IRN(*),IW(LW)
      INTEGER ICNTL(*),INFO(20)
C     ..
C     .. Local Scalars ..
      INTEGER I,ID,J,JN,K,K1,K2,L,LAST,LR,N1,NDUP
C     ..
C     .. Executable Statements ..
C
C INITIALIZE INFO(2) AND COUNT IN IPE THE
C     NUMBERS OF NON-ZEROS IN THE ROWS AND MOVE ROW AND COLUMN
C     NUMBERS INTO IW.
      INFO(2) = 0
      DO 10 I = 1,N
        IPE(I) = 0
   10 CONTINUE
      LR = NZ
      IF (NZ.EQ.0) GO TO 120
      DO 110 K = 1,NZ
        I = IRN(K)
        J = ICN(K)
        IF (I.LT.J) THEN
          IF (I.GE.1 .AND. J.LE.N) GO TO 90
        ELSE IF (I.GT.J) THEN
          IF (J.GE.1 .AND. I.LE.N) GO TO 90
        ELSE
          IF (I.GE.1 .AND. I.LE.N) GO TO 80
        END IF
        INFO(2) = INFO(2) + 1
        INFO(1) = 1
        IF (INFO(2).LE.1 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=60) INFO(1)
        END IF

   60   FORMAT (' *** WARNING MESSAGE FROM SUBROUTINE MA27AD',
     +          '  *** INFO(1) =',I2)

        IF (INFO(2).LE.10 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=70) K,I,J
        END IF

   70   FORMAT (I6,'TH NON-ZERO (IN ROW',I6,' AND COLUMN',I6,
     +         ') IGNORED')

   80   I = 0
        J = 0
        GO TO 100

   90   IPE(I) = IPE(I) + 1
        IPE(J) = IPE(J) + 1
  100   IW(K) = J
        LR = LR + 1
        IW(LR) = I
  110 CONTINUE
C
C ACCUMULATE ROW COUNTS TO GET POINTERS TO ROW STARTS IN BOTH IPE AND IQ
C     AND INITIALIZE FLAG
  120 IQ(1) = 1
      N1 = N - 1
      IF (N1.LE.0) GO TO 140
      DO 130 I = 1,N1
        FLAG(I) = 0
        IF (IPE(I).EQ.0) IPE(I) = -1
        IQ(I+1) = IPE(I) + IQ(I) + 1
        IPE(I) = IQ(I)
  130 CONTINUE
  140 LAST = IPE(N) + IQ(N)
      FLAG(N) = 0
      IF (LR.GE.LAST) GO TO 160
      K1 = LR + 1
      DO 150 K = K1,LAST
        IW(K) = 0
  150 CONTINUE
  160 IPE(N) = IQ(N)
      IWFR = LAST + 1
      IF (NZ.EQ.0) GO TO 230
C
C RUN THROUGH PUTTING THE MATRIX ELEMENTS IN THE RIGHT PLACE
C     BUT WITH SIGNS INVERTED. IQ IS USED FOR HOLDING RUNNING POINTERS
C     AND IS LEFT HOLDING POINTERS TO ROW ENDS.
      DO 220 K = 1,NZ
        J = IW(K)
        IF (J.LE.0) GO TO 220
        L = K
        IW(K) = 0
        DO 210 ID = 1,NZ
          IF (L.GT.NZ) GO TO 170
          L = L + NZ
          GO TO 180

  170     L = L - NZ
  180     I = IW(L)
          IW(L) = 0
          IF (I.LT.J) GO TO 190
          L = IQ(J) + 1
          IQ(J) = L
          JN = IW(L)
          IW(L) = -I
          GO TO 200

  190     L = IQ(I) + 1
          IQ(I) = L
          JN = IW(L)
          IW(L) = -J
  200     J = JN
          IF (J.LE.0) GO TO 220
  210   CONTINUE
  220 CONTINUE
C
C RUN THROUGH RESTORING SIGNS, REMOVING DUPLICATES AND SETTING THE
C     MATE OF EACH NON-ZERO.
C NDUP COUNTS THE NUMBER OF DUPLICATE ELEMENTS.
  230 NDUP = 0
      DO 280 I = 1,N
        K1 = IPE(I) + 1
        K2 = IQ(I)
        IF (K1.LE.K2) GO TO 240
C ROW IS EMPTY. SET POINTER TO ZERO.
        IPE(I) = 0
        IQ(I) = 0
        GO TO 280
C ON ENTRY TO THIS LOOP FLAG(J).LT.I FOR J=1,2,...,N. DURING THE LOOP
C     FLAG(J) IS SET TO I IF A NON-ZERO IN COLUMN J IS FOUND. THIS
C     PERMITS DUPLICATES TO BE RECOGNIZED QUICKLY.
  240   DO 260 K = K1,K2
          J = -IW(K)
          IF (J.LE.0) GO TO 270
          L = IQ(J) + 1
          IQ(J) = L
          IW(L) = I
          IW(K) = J
          IF (FLAG(J).NE.I) GO TO 250
          NDUP = NDUP + 1
          IW(L) = 0
          IW(K) = 0
  250     FLAG(J) = I
  260   CONTINUE
  270   IQ(I) = IQ(I) - IPE(I)
        IF (NDUP.EQ.0) IW(K1-1) = IQ(I)
  280 CONTINUE
      IF (NDUP.EQ.0) GO TO 310
C
C COMPRESS IW TO REMOVE DUMMY ENTRIES CAUSED BY DUPLICATES.
      IWFR = 1
      DO 300 I = 1,N
        K1 = IPE(I) + 1
        IF (K1.EQ.1) GO TO 300
        K2 = IQ(I) + IPE(I)
        L = IWFR
        IPE(I) = IWFR
        IWFR = IWFR + 1
        DO 290 K = K1,K2
          IF (IW(K).EQ.0) GO TO 290
          IW(IWFR) = IW(K)
          IWFR = IWFR + 1
  290   CONTINUE
        IW(L) = IWFR - L - 1
  300 CONTINUE
  310 RETURN

      END
      SUBROUTINE MA27HD(N,IPE,IW,LW,IWFR,NV,NXT,LST,IPD,FLAG,IOVFLO,
     +                 NCMPA,FRATIO)
C
C ANALYSIS SUBROUTINE
C
C GIVEN REPRESENTATION OF THE WHOLE MATRIX (EXCLUDING DIAGONAL)
C     PERFORM MINIMUM DEGREE ORDERING, CONSTRUCTING TREE POINTERS.
C     IT WORKS WITH SUPERVARIABLES WHICH ARE COLLECTIONS OF ONE OR MORE
C     VARIABLES, STARTING WITH SUPERVARIABLE I CONTAINING VARIABLE I FOR
C     I = 1,2,...,N. ALL VARIABLES IN A SUPERVARIABLE ARE ELIMINATED
C     TOGETHER. EACH SUPERVARIABLE HAS AS NUMERICAL NAME THAT OF ONE
C     OF ITS VARIABLES (ITS PRINCIPAL VARIABLE).
C
C N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
C IPE(I) MUST BE SET TO POINT TO THE POSITION IN IW OF THE
C     START OF ROW I OR HAVE THE VALUE ZERO IF ROW I HAS NO OFF-
C     DIAGONAL NON-ZEROS. DURING EXECUTION IT IS USED AS FOLLOWS. IF
C     SUPERVARIABLE I IS ABSORBED INTO SUPERVARIABLE J THEN IPE(I)=-J.
C     IF SUPERVARIABLE I IS ELIMINATED THEN IPE(I) EITHER POINTS TO THE
C     LIST OF SUPERVARIABLES FOR CREATED ELEMENT I OR IS ZERO IF
C     THE CREATED ELEMENT IS NULL. IF ELEMENT I
C     IS ABSORBED INTO ELEMENT J THEN IPE(I)=-J.
C IW MUST BE SET ON ENTRY TO HOLD LISTS OF VARIABLES BY
C     ROWS, EACH LIST BEING HEADED BY ITS LENGTH.
C     DURING EXECUTION THESE LISTS ARE REVISED AND HOLD
C     LISTS OF ELEMENTS AND SUPERVARIABLES. THE ELEMENTS
C     ALWAYS HEAD THE LISTS. WHEN A SUPERVARIABLE
C     IS ELIMINATED ITS LIST IS REPLACED BY A LIST OF SUPERVARIABLES
C     IN THE NEW ELEMENT.
C LW MUST BE SET TO THE LENGTH OF IW. IT IS NOT ALTERED.
C IWFR MUST BE SET TO THE POSITION IN IW OF THE FIRST FREE VARIABLE.
C     IT IS REVISED DURING EXECUTION AND CONTINUES TO HAVE THIS MEANING.
C NV(JS) NEED NOT BE SET. DURING EXECUTION IT IS ZERO IF
C     JS IS NOT A PRINCIPAL VARIABLE AND IF IT IS IT HOLDS
C     THE NUMBER OF VARIABLES IN SUPERVARIABLE JS. FOR ELIMINATED
C     VARIABLES IT IS SET TO THE DEGREE AT THE TIME OF ELIMINATION.
C NXT(JS) NEED NOT BE SET. DURING EXECUTION IT IS THE NEXT
C     SUPERVARIABLE HAVING THE SAME DEGREE AS JS, OR ZERO
C     IF IT IS LAST IN ITS LIST.
C LST(JS) NEED NOT BE SET. DURING EXECUTION IT IS THE
C     LAST SUPERVARIABLE HAVING THE SAME DEGREE AS JS OR
C     -(ITS DEGREE) IF IT IS FIRST IN ITS LIST.
C IPD(ID) NEED NOT BE SET. DURING EXECUTION IT
C     IS THE FIRST SUPERVARIABLE WITH DEGREE ID OR ZERO
C     IF THERE ARE NONE.
C FLAG IS USED AS WORKSPACE FOR ELEMENT AND SUPERVARIABLE FLAGS.
C     WHILE THE CODE IS FINDING THE DEGREE OF SUPERVARIABLE IS
C     FLAG HAS THE FOLLOWING VALUES.
C     A) FOR THE CURRENT PIVOT/NEW ELEMENT ME
C           FLAG(ME)=-1
C     B) FOR VARIABLES JS
C           FLAG(JS)=-1 IF JS IS NOT A PRINCIPAL VARIABLE
C           FLAG(JS)=0 IF JS IS A SUPERVARIABLE IN THE NEW ELEMENT
C           FLAG(JS)=NFLG IF JS IS A SUPERVARIABLE NOT IN THE NEW
C                 ELEMENT THAT HAS BEEN COUNTED IN THE DEGREE
C                 CALCULATION
C           FLAG(JS).GT.NFLG IF JS IS A SUPERVARIABLE NOT IN THE NEW
C                 ELEMENT THAT HAS NOT BEEN COUNTED IN THE DEGREE
C                 CALCULATION
C     C) FOR ELEMENTS IE
C           FLAG(IE)=-1 IF ELEMENT IE HAS BEEN MERGED INTO ANOTHER
C           FLAG(IE)=-NFLG IF ELEMENT IE HAS BEEN USED IN THE DEGREE
C                 CALCULATION FOR IS.
C           FLAG(IE).LT.-NFLG IF ELEMENT IE HAS NOT BEEN USED IN THE
C                 DEGREE CALCULATION FOR IS
C IOVFLO see ICNTL(4) in MA27A/AD.
C NCMPA see INFO(11) in MA27A/AD.
C FRATIO see CNTL(2) in MA27A/AD.
C 
C     .. Scalar Arguments ..
      DOUBLE PRECISION FRATIO
      INTEGER IWFR,LW,N,IOVFLO,NCMPA
C     ..
C     .. Array Arguments ..
      INTEGER FLAG(N),IPD(N),IPE(N),IW(LW),LST(N),NV(N),NXT(N)
C     ..
C     .. Local Scalars ..
      INTEGER I,ID,IDL,IDN,IE,IP,IS,JP,JP1,JP2,JS,K,K1,K2,KE,KP,KP0,KP1,
     +        KP2,KS,L,LEN,LIMIT,LN,LS,LWFR,MD,ME,ML,MS,NEL,NFLG,NP,
     +        NP0,NS,NVPIV,NVROOT,ROOT
C LIMIT  Limit on number of variables for putting node in root.
C NVROOT Number of variables in the root node
C ROOT   Index of the root node (N+1 if none chosen yet).
C     ..
C     .. External Subroutines ..
      EXTERNAL MA27UD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN
C     ..
C If a column of the reduced matrix has relative density greater than
C CNTL(2), it is forced into the root. All such columns are taken to
C have sparsity pattern equal to their merged patterns, so the fill
C and operation counts may be overestimated.
C
C IS,JS,KS,LS,MS,NS ARE USED TO REFER TO SUPERVARIABLES.
C IE,JE,KE ARE USED TO REFER TO ELEMENTS.
C IP,JP,KP,K,NP ARE USED TO POINT TO LISTS OF ELEMENTS.
C     OR SUPERVARIABLES.
C ID IS USED FOR THE DEGREE OF A SUPERVARIABLE.
C MD IS USED FOR THE CURRENT MINIMUM DEGREE.
C IDN IS USED FOR THE NO. OF VARIABLES IN A NEWLY CREATED ELEMENT
C NEL IS USED TO HOLD THE NO. OF VARIABLES THAT HAVE BEEN
C     ELIMINATED.
C ME=MS IS THE NAME OF THE SUPERVARIABLE ELIMINATED AND
C     OF THE ELEMENT CREATED IN THE MAIN LOOP.
C NFLG IS USED FOR THE CURRENT FLAG VALUE IN ARRAY FLAG. IT STARTS
C     WITH THE VALUE IOVFLO AND IS REDUCED BY 1 EACH TIME IT IS USED
C     UNTIL IT HAS THE VALUE 2 WHEN IT IS RESET TO THE VALUE IOVFLO.
C
C     .. Executable Statements ..
C INITIALIZATIONS
      DO 10 I = 1,N
        IPD(I) = 0
        NV(I) = 1
        FLAG(I) = IOVFLO
   10 CONTINUE
      MD = 1
      NCMPA = 0
      NFLG = IOVFLO
      NEL = 0
      ROOT = N+1
      NVROOT = 0
C
C LINK TOGETHER VARIABLES HAVING SAME DEGREE
      DO 30 IS = 1,N
        K = IPE(IS)
        IF (K.LE.0) GO TO 20
        ID = IW(K) + 1
        NS = IPD(ID)
        IF (NS.GT.0) LST(NS) = IS
        NXT(IS) = NS
        IPD(ID) = IS
        LST(IS) = -ID
        GO TO 30
C WE HAVE A VARIABLE THAT CAN BE ELIMINATED AT ONCE BECAUSE THERE IS
C     NO OFF-DIAGONAL NON-ZERO IN ITS ROW.
   20   NEL = NEL + 1
        FLAG(IS) = -1
        NXT(IS) = 0
        LST(IS) = 0
   30 CONTINUE
C
C START OF MAIN LOOP
C
      DO 340 ML = 1,N
C LEAVE LOOP IF ALL VARIABLES HAVE BEEN ELIMINATED.
        IF (NEL+NVROOT+1.GE.N) GO TO 350
C
C FIND NEXT SUPERVARIABLE FOR ELIMINATION.
        DO 40 ID = MD,N
          MS = IPD(ID)
          IF (MS.GT.0) GO TO 50
   40   CONTINUE
   50   MD = ID
C NVPIV HOLDS THE NUMBER OF VARIABLES IN THE PIVOT.
        NVPIV = NV(MS)
C
C REMOVE CHOSEN VARIABLE FROM LINKED LIST
        NS = NXT(MS)
        NXT(MS) = 0
        LST(MS) = 0
        IF (NS.GT.0) LST(NS) = -ID
        IPD(ID) = NS
        ME = MS
        NEL = NEL + NVPIV
C IDN HOLDS THE DEGREE OF THE NEW ELEMENT.
        IDN = 0
C
C RUN THROUGH THE LIST OF THE PIVOTAL SUPERVARIABLE, SETTING TREE
C     POINTERS AND CONSTRUCTING NEW LIST OF SUPERVARIABLES.
C KP IS A POINTER TO THE CURRENT POSITION IN THE OLD LIST.
        KP = IPE(ME)
        FLAG(MS) = -1
C IP POINTS TO THE START OF THE NEW LIST.
        IP = IWFR
C LEN HOLDS THE LENGTH OF THE LIST ASSOCIATED WITH THE PIVOT.
        LEN = IW(KP)
        DO 140 KP1 = 1,LEN
          KP = KP + 1
          KE = IW(KP)
C JUMP IF KE IS AN ELEMENT THAT HAS NOT BEEN MERGED INTO ANOTHER.
          IF (FLAG(KE).LE.-2) GO TO 60
C JUMP IF KE IS AN ELEMENT THAT HAS BEEN MERGED INTO ANOTHER OR IS
C     A SUPERVARIABLE THAT HAS BEEN ELIMINATED.
          IF (FLAG(KE).LE.0) THEN
             IF (IPE(KE).NE.-ROOT) GO TO 140
C KE has been merged into the root
             KE = ROOT
             IF (FLAG(KE).LE.0) GO TO 140
          END IF
C WE HAVE A SUPERVARIABLE. PREPARE TO SEARCH REST OF LIST.
          JP = KP - 1
          LN = LEN - KP1 + 1
          IE = MS
          GO TO 70
C SEARCH VARIABLE LIST OF ELEMENT KE, USING JP AS A POINTER TO IT.
   60     IE = KE
          JP = IPE(IE)
          LN = IW(JP)
C
C SEARCH FOR DIFFERENT SUPERVARIABLES AND ADD THEM TO THE NEW LIST,
C     COMPRESSING WHEN NECESSARY. THIS LOOP IS EXECUTED ONCE FOR
C     EACH ELEMENT IN THE LIST AND ONCE FOR ALL THE SUPERVARIABLES
C     IN THE LIST.
   70     DO 130 JP1 = 1,LN
            JP = JP + 1
            IS = IW(JP)
C JUMP IF IS IS NOT A PRINCIPAL VARIABLE OR HAS ALREADY BEEN COUNTED.
            IF (FLAG(IS).LE.0) THEN
               IF (IPE(IS).EQ.-ROOT) THEN
C IS has been merged into the root
                  IS = ROOT
                  IW(JP) = ROOT
                  IF (FLAG(IS).LE.0) GO TO 130
               ELSE
                  GO TO 130
               END IF
            END IF
            FLAG(IS) = 0
            IF (IWFR.LT.LW) GO TO 100
C PREPARE FOR COMPRESSING IW BY ADJUSTING POINTERS AND
C     LENGTHS SO THAT THE LISTS BEING SEARCHED IN THE INNER AND OUTER
C     LOOPS CONTAIN ONLY THE REMAINING ENTRIES.
            IPE(MS) = KP
            IW(KP) = LEN - KP1
            IPE(IE) = JP
            IW(JP) = LN - JP1
C COMPRESS IW
            CALL MA27UD(N,IPE,IW,IP-1,LWFR,NCMPA)
C COPY NEW LIST FORWARD
            JP2 = IWFR - 1
            IWFR = LWFR
            IF (IP.GT.JP2) GO TO 90
            DO 80 JP = IP,JP2
              IW(IWFR) = IW(JP)
              IWFR = IWFR + 1
   80       CONTINUE
C ADJUST POINTERS FOR THE NEW LIST AND THE LISTS BEING SEARCHED.
   90       IP = LWFR
            JP = IPE(IE)
            KP = IPE(ME)
C STORE IS IN NEW LIST.
  100       IW(IWFR) = IS
            IDN = IDN + NV(IS)
            IWFR = IWFR + 1
C REMOVE IS FROM DEGREE LINKED LIST
            LS = LST(IS)
            LST(IS) = 0
            NS = NXT(IS)
            NXT(IS) = 0
            IF (NS.GT.0) LST(NS) = LS
            IF (LS.LT.0) THEN
              LS = -LS
              IPD(LS) = NS
            ELSE IF (LS.GT.0) THEN
              NXT(LS) = NS
            END IF
  130     CONTINUE
C JUMP IF WE HAVE JUST BEEN SEARCHING THE VARIABLES AT THE END OF
C     THE LIST OF THE PIVOT.
          IF (IE.EQ.MS) GO TO 150
C SET TREE POINTER AND FLAG TO INDICATE ELEMENT IE IS ABSORBED INTO
C     NEW ELEMENT ME.
          IPE(IE) = -ME
          FLAG(IE) = -1
  140   CONTINUE

C STORE THE DEGREE OF THE PIVOT.
  150   NV(MS) = IDN + NVPIV
C JUMP IF NEW ELEMENT IS NULL.
        IF (IWFR.EQ.IP) GO TO 330
        K1 = IP
        K2 = IWFR - 1
C
C RUN THROUGH NEW LIST OF SUPERVARIABLES REVISING EACH ASSOCIATED LIST,
C     RECALCULATING DEGREES AND REMOVING DUPLICATES.
        LIMIT = NINT(FRATIO*(N-NEL))
        DO 310 K = K1,K2
          IS = IW(K)
          IF (IS.EQ.ROOT) GO TO 310
          IF (NFLG.GT.2) GO TO 170
C RESET FLAG VALUES TO +/-IOVFLO.
          DO 160 I = 1,N
            IF (FLAG(I).GT.0) FLAG(I) = IOVFLO
            IF (FLAG(I).LE.-2) FLAG(I) = -IOVFLO
  160     CONTINUE
          NFLG = IOVFLO
C REDUCE NFLG BY ONE TO CATER FOR THIS SUPERVARIABLE.
  170     NFLG = NFLG - 1
C BEGIN WITH THE DEGREE OF THE NEW ELEMENT. ITS VARIABLES MUST ALWAYS
C     BE COUNTED DURING THE DEGREE CALCULATION AND THEY ARE ALREADY
C     FLAGGED WITH THE VALUE 0.
          ID = IDN
C RUN THROUGH THE LIST ASSOCIATED WITH SUPERVARIABLE IS
          KP1 = IPE(IS) + 1
C NP POINTS TO THE NEXT ENTRY IN THE REVISED LIST.
          NP = KP1
          KP2 = IW(KP1-1) + KP1 - 1
          DO 220 KP = KP1,KP2
            KE = IW(KP)
C TEST WHETHER KE IS AN ELEMENT, A REDUNDANT ENTRY OR A SUPERVARIABLE.
          IF (FLAG(KE).EQ.-1) THEN
             IF (IPE(KE).NE.-ROOT) GO TO 220
C KE has been merged into the root
             KE = ROOT
             IW(KP) = ROOT
             IF (FLAG(KE).EQ.-1) GO TO 220
          END IF
          IF (FLAG(KE).GE.0) GO TO 230
C SEARCH LIST OF ELEMENT KE, REVISING THE DEGREE WHEN NEW VARIABLES
C     FOUND.
            JP1 = IPE(KE) + 1
            JP2 = IW(JP1-1) + JP1 - 1
            IDL = ID
            DO 190 JP = JP1,JP2
              JS = IW(JP)
C JUMP IF JS HAS ALREADY BEEN COUNTED.
              IF (FLAG(JS).LE.NFLG) GO TO 190
              ID = ID + NV(JS)
              FLAG(JS) = NFLG
  190       CONTINUE
C JUMP IF ONE OR MORE NEW SUPERVARIABLES WERE FOUND.
            IF (ID.GT.IDL) GO TO 210
C CHECK WHETHER EVERY VARIABLE OF ELEMENT KE IS IN NEW ELEMENT ME.
            DO 200 JP = JP1,JP2
              JS = IW(JP)
              IF (FLAG(JS).NE.0) GO TO 210
  200       CONTINUE
C SET TREE POINTER AND FLAG TO INDICATE THAT ELEMENT KE IS ABSORBED
C     INTO NEW ELEMENT ME.
            IPE(KE) = -ME
            FLAG(KE) = -1
            GO TO 220
C STORE ELEMENT KE IN THE REVISED LIST FOR SUPERVARIABLE IS AND FLAG IT.
  210       IW(NP) = KE
            FLAG(KE) = -NFLG
            NP = NP + 1
  220     CONTINUE
          NP0 = NP
          GO TO 250
C TREAT THE REST OF THE LIST ASSOCIATED WITH SUPERVARIABLE IS. IT
C     CONSISTS ENTIRELY OF SUPERVARIABLES.
  230     KP0 = KP
          NP0 = NP
          DO 240 KP = KP0,KP2
            KS = IW(KP)
            IF (FLAG(KS).LE.NFLG) THEN
               IF (IPE(KS).EQ.-ROOT) THEN
                  KS = ROOT
                  IW(KP) = ROOT
                  IF (FLAG(KS).LE.NFLG) GO TO 240
               ELSE
                  GO TO 240
               END IF
            END IF
C ADD TO DEGREE, FLAG SUPERVARIABLE KS AND ADD IT TO NEW LIST.
            ID = ID + NV(KS)
            FLAG(KS) = NFLG
            IW(NP) = KS
            NP = NP + 1
  240     CONTINUE
C MOVE FIRST SUPERVARIABLE TO END OF LIST, MOVE FIRST ELEMENT TO END
C     OF ELEMENT PART OF LIST AND ADD NEW ELEMENT TO FRONT OF LIST.
  250     IF (ID.GE.LIMIT) GO TO 295
          IW(NP) = IW(NP0)
          IW(NP0) = IW(KP1)
          IW(KP1) = ME
C STORE THE NEW LENGTH OF THE LIST.
          IW(KP1-1) = NP - KP1 + 1
C
C CHECK WHETHER ROW IS IS IDENTICAL TO ANOTHER BY LOOKING IN LINKED
C     LIST OF SUPERVARIABLES WITH DEGREE ID AT THOSE WHOSE LISTS HAVE
C     FIRST ENTRY ME. NOTE THAT THOSE CONTAINING ME COME FIRST SO THE
C     SEARCH CAN BE TERMINATED WHEN A LIST NOT STARTING WITH ME IS
C     FOUND.
          JS = IPD(ID)
          DO 280 L = 1,N
            IF (JS.LE.0) GO TO 300
            KP1 = IPE(JS) + 1
            IF (IW(KP1).NE.ME) GO TO 300
C JS HAS SAME DEGREE AND IS ACTIVE. CHECK IF IDENTICAL TO IS.
            KP2 = KP1 - 1 + IW(KP1-1)
            DO 260 KP = KP1,KP2
              IE = IW(KP)
C JUMP IF IE IS A SUPERVARIABLE OR AN ELEMENT NOT IN THE LIST OF IS.
              IF (ABS(FLAG(IE)+0).GT.NFLG) GO TO 270
  260       CONTINUE
            GO TO 290

  270       JS = NXT(JS)
  280     CONTINUE
C SUPERVARIABLE AMALGAMATION. ROW IS IS IDENTICAL TO ROW JS.
C REGARD ALL VARIABLES IN THE TWO SUPERVARIABLES AS BEING IN IS. SET
C     TREE POINTER, FLAG AND NV ENTRIES.
  290     IPE(JS) = -IS
          NV(IS) = NV(IS) + NV(JS)
          NV(JS) = 0
          FLAG(JS) = -1
C REPLACE JS BY IS IN LINKED LIST.
          NS = NXT(JS)
          LS = LST(JS)
          IF (NS.GT.0) LST(NS) = IS
          IF (LS.GT.0) NXT(LS) = IS
          LST(IS) = LS
          NXT(IS) = NS
          LST(JS) = 0
          NXT(JS) = 0
          IF (IPD(ID).EQ.JS) IPD(ID) = IS
          GO TO 310
C Treat IS as full. Merge it into the root node.
  295     IF (NVROOT.EQ.0) THEN
            ROOT = IS
            IPE(IS) = 0
          ELSE
            IW(K) = ROOT
            IPE(IS) = -ROOT
            NV(ROOT) = NV(ROOT) + NV(IS)
            NV(IS) = 0
            FLAG(IS) = -1
          END IF
          NVROOT = NV(ROOT)
          GO TO 310
C INSERT IS INTO LINKED LIST OF SUPERVARIABLES OF SAME DEGREE.
  300     NS = IPD(ID)
          IF (NS.GT.0) LST(NS) = IS
          NXT(IS) = NS
          IPD(ID) = IS
          LST(IS) = -ID
          MD = MIN(MD,ID)
  310   CONTINUE
C
C RESET FLAGS FOR SUPERVARIABLES IN NEWLY CREATED ELEMENT AND
C     REMOVE THOSE ABSORBED INTO OTHERS.
        DO 320 K = K1,K2
          IS = IW(K)
          IF (NV(IS).EQ.0) GO TO 320
          FLAG(IS) = NFLG
          IW(IP) = IS
          IP = IP + 1
  320   CONTINUE
        IWFR = K1
        FLAG(ME) = -NFLG
C MOVE FIRST ENTRY TO END TO MAKE ROOM FOR LENGTH.
        IW(IP) = IW(K1)
        IW(K1) = IP - K1
C SET POINTER FOR NEW ELEMENT AND RESET IWFR.
        IPE(ME) = K1
        IWFR = IP + 1
        GO TO 335

  330   IPE(ME) = 0
C
  335   CONTINUE
  340 CONTINUE
C

C Absorb any remaining variables into the root
  350 DO 360 IS = 1,N
        IF(NXT(IS).NE.0 .OR. LST(IS).NE.0) THEN
          IF (NVROOT.EQ.0) THEN
            ROOT = IS
            IPE(IS) = 0
          ELSE
            IPE(IS) = -ROOT
          END IF
          NVROOT = NVROOT + NV(IS)
          NV(IS) = 0
         END IF
  360 CONTINUE
C Link any remaining elements to the root
      DO 370 IE = 1,N
        IF (IPE(IE).GT.0) IPE(IE) = -ROOT
  370 CONTINUE
      IF(NVROOT.GT.0)NV(ROOT)=NVROOT
      END

      SUBROUTINE MA27UD(N,IPE,IW,LW,IWFR,NCMPA)
C COMPRESS LISTS HELD BY MA27H/HD AND MA27K/KD IN IW AND ADJUST POINTERS
C     IN IPE TO CORRESPOND.
C N IS THE MATRIX ORDER. IT IS NOT ALTERED.
C IPE(I) POINTS TO THE POSITION IN IW OF THE START OF LIST I OR IS
C     ZERO IF THERE IS NO LIST I. ON EXIT IT POINTS TO THE NEW POSITION.
C IW HOLDS THE LISTS, EACH HEADED BY ITS LENGTH. ON OUTPUT THE SAME
C     LISTS ARE HELD, BUT THEY ARE NOW COMPRESSED TOGETHER.
C LW HOLDS THE LENGTH OF IW. IT IS NOT ALTERED.
C IWFR NEED NOT BE SET ON ENTRY. ON EXIT IT POINTS TO THE FIRST FREE
C     LOCATION IN IW.
C     ON RETURN IT IS SET TO THE FIRST FREE LOCATION IN IW.
C NCMPA see INFO(11) in MA27A/AD.
C
C     .. Scalar Arguments ..
      INTEGER IWFR,LW,N,NCMPA
C     ..
C     .. Array Arguments ..
      INTEGER IPE(N),IW(LW)
C     ..
C     .. Local Scalars ..
      INTEGER I,IR,K,K1,K2,LWFR
C     ..
C     .. Executable Statements ..
      NCMPA = NCMPA + 1
C PREPARE FOR COMPRESSING BY STORING THE LENGTHS OF THE
C     LISTS IN IPE AND SETTING THE FIRST ENTRY OF EACH LIST TO
C     -(LIST NUMBER).
      DO 10 I = 1,N
        K1 = IPE(I)
        IF (K1.LE.0) GO TO 10
        IPE(I) = IW(K1)
        IW(K1) = -I
   10 CONTINUE
C
C COMPRESS
C IWFR POINTS JUST BEYOND THE END OF THE COMPRESSED FILE.
C LWFR POINTS JUST BEYOND THE END OF THE UNCOMPRESSED FILE.
      IWFR = 1
      LWFR = IWFR
      DO 60 IR = 1,N
        IF (LWFR.GT.LW) GO TO 70
C SEARCH FOR THE NEXT NEGATIVE ENTRY.
        DO 20 K = LWFR,LW
          IF (IW(K).LT.0) GO TO 30
   20   CONTINUE
        GO TO 70
C PICK UP ENTRY NUMBER, STORE LENGTH IN NEW POSITION, SET NEW POINTER
C     AND PREPARE TO COPY LIST.
   30   I = -IW(K)
        IW(IWFR) = IPE(I)
        IPE(I) = IWFR
        K1 = K + 1
        K2 = K + IW(IWFR)
        IWFR = IWFR + 1
        IF (K1.GT.K2) GO TO 50
C COPY LIST TO NEW POSITION.
        DO 40 K = K1,K2
          IW(IWFR) = IW(K)
          IWFR = IWFR + 1
   40   CONTINUE
   50   LWFR = K2 + 1
   60 CONTINUE
   70 RETURN

      END
      SUBROUTINE MA27JD(N,NZ,IRN,ICN,PERM,IW,LW,IPE,IQ,FLAG,IWFR,
     +                 ICNTL,INFO)
C
C SORT PRIOR TO CALLING ANALYSIS ROUTINE MA27K/KD.
C
C GIVEN THE POSITIONS OF THE OFF-DIAGONAL NON-ZEROS OF A SYMMETRIC
C     MATRIX AND A PERMUTATION, CONSTRUCT THE SPARSITY PATTERN
C     OF THE STRICTLY UPPER TRIANGULAR PART OF THE PERMUTED MATRIX.
C     EITHER ONE OF A PAIR (I,J),(J,I) MAY BE USED TO REPRESENT
C     THE PAIR. DIAGONAL ELEMENTS ARE IGNORED. NO CHECK IS MADE
C     FOR DUPLICATE ELEMENTS UNLESS ANY ROW HAS MORE THAN ICNTL(4)
C     NON-ZEROS, IN WHICH CASE DUPLICATES ARE REMOVED.
C
C N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
C NZ MUST BE SET TO THE NUMBER OF NON-ZEROS INPUT. IT IS NOT
C     ALTERED.
C IRN(I),I=1,2,...,NZ MUST BE SET TO THE ROW INDICES OF THE
C     NON-ZEROS ON INPUT. IT IS NOT ALTERED UNLESS EQUIVALENCED WITH IW.
C     IRN(1) MAY BE EQUIVALENCED WITH IW(1).
C ICN(I),I=1,2,...,NZ MUST BE SET TO THE COLUMN INDICES OF THE
C     NON-ZEROS ON INPUT. IT IS NOT ALTERED UNLESS EQUIVALENCED
C     WITH IW.ICN(1) MAY BE EQUIVELENCED WITH IW(K),K.GT.NZ.
C PERM(I) MUST BE SET TO HOLD THE POSITION OF VARIABLE I IN THE
C     PERMUTED ORDER. IT IS NOT ALTERED.
C IW NEED NOT BE SET ON INPUT. ON OUTPUT IT CONTAINS LISTS OF
C     COLUMN NUMBERS, EACH LIST BEING HEADED BY ITS LENGTH.
C LW MUST BE SET TO THE LENGTH OF IW. IT MUST BE AT LEAST
C     MAX(NZ,N+(NO. OF OFF-DIAGONAL NON-ZEROS)). IT IS NOT ALTERED.
C IPE NEED NOT BE SET ON INPUT. ON OUTPUT IPE(I) POINTS TO THE START OF
C     THE ENTRY IN IW FOR ROW I, OR IS ZERO IF THERE IS NO ENTRY.
C IQ NEED NOT BE SET. ON OUTPUT IQ(I) CONTAINS THE NUMBER OF
C     OFF-DIAGONAL NON-ZEROS IN ROW I, INCLUDING DUPLICATES.
C FLAG IS USED FOR WORKSPACE TO HOLD FLAGS TO PERMIT DUPLICATE
C     ENTRIES TO BE IDENTIFIED QUICKLY.
C IWFR NEED NOT BE SET ON INPUT. ON OUTPUT IT POINTS TO THE FIRST
C     UNUSED LOCATION IN IW.
C ICNTL is an INTEGER array of length 30, see MA27A/AD.
C INFO is an INTEGER array of length 20, see MA27A/AD.
C
C     .. Scalar Arguments ..
      INTEGER IWFR,LW,N,NZ
C     ..
C     .. Array Arguments ..
      INTEGER FLAG(N),ICN(*),IPE(N),IQ(N),IRN(*),IW(LW),PERM(N)
      INTEGER ICNTL(30),INFO(20)
C     ..
C     .. Local Scalars ..
      INTEGER I,ID,IN,J,JDUMMY,K,K1,K2,L,LBIG,LEN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX
C     ..
C     .. Executable Statements ..
C
C INITIALIZE INFO(1), INFO(2) AND IQ
      INFO(1) = 0
      INFO(2) = 0
      DO 10 I = 1,N
        IQ(I) = 0
   10 CONTINUE
C
C COUNT THE NUMBERS OF NON-ZEROS IN THE ROWS, PRINT WARNINGS ABOUT
C     OUT-OF-RANGE INDICES AND TRANSFER GENUINE ROW NUMBERS
C     (NEGATED) INTO IW.
      IF (NZ.EQ.0) GO TO 110
      DO 100 K = 1,NZ
        I = IRN(K)
        J = ICN(K)
        IW(K) = -I
        IF(I.LT.J) THEN
          IF (I.GE.1 .AND. J.LE.N) GO TO 80
        ELSE IF(I.GT.J) THEN
          IF (J.GE.1 .AND. I.LE.N) GO TO 80
        ELSE
          IW(K) = 0
          IF (I.GE.1 .AND. I.LE.N) GO TO 100
        END IF
        INFO(2) = INFO(2) + 1
        INFO(1) = 1
        IW(K) = 0
        IF (INFO(2).LE.1 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=60) INFO(1)
        END IF

   60   FORMAT (' *** WARNING MESSAGE FROM SUBROUTINE MA27AD',
     +          '  *** INFO(1) =',I2)

        IF (INFO(2).LE.10 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=70) K,I,J
        END IF

   70   FORMAT (I6,'TH NON-ZERO (IN ROW',I6,' AND COLUMN',I6,
     +         ') IGNORED')

        GO TO 100

   80   IF (PERM(J).GT.PERM(I)) GO TO 90
        IQ(J) = IQ(J) + 1
        GO TO 100

   90   IQ(I) = IQ(I) + 1
  100 CONTINUE
C
C ACCUMULATE ROW COUNTS TO GET POINTERS TO ROW ENDS
C     IN IPE.
  110 IWFR = 1
      LBIG = 0
      DO 120 I = 1,N
        L = IQ(I)
        LBIG = MAX(L,LBIG)
        IWFR = IWFR + L
        IPE(I) = IWFR - 1
  120 CONTINUE
C
C PERFORM IN-PLACE SORT
      IF (NZ.EQ.0) GO TO 250
      DO 160 K = 1,NZ
        I = -IW(K)
        IF (I.LE.0) GO TO 160
        L = K
        IW(K) = 0
        DO 150 ID = 1,NZ
          J = ICN(L)
          IF (PERM(I).LT.PERM(J)) GO TO 130
          L = IPE(J)
          IPE(J) = L - 1
          IN = IW(L)
          IW(L) = I
          GO TO 140

  130     L = IPE(I)
          IPE(I) = L - 1
          IN = IW(L)
          IW(L) = J
  140     I = -IN
          IF (I.LE.0) GO TO 160
  150   CONTINUE
  160 CONTINUE
C
C MAKE ROOM IN IW FOR ROW LENGTHS AND INITIALIZE FLAG.
      K = IWFR - 1
      L = K + N
      IWFR = L + 1
      DO 190 I = 1,N
        FLAG(I) = 0
        J = N + 1 - I
        LEN = IQ(J)
        IF (LEN.LE.0) GO TO 180
        DO 170 JDUMMY = 1,LEN
          IW(L) = IW(K)
          K = K - 1
          L = L - 1
  170   CONTINUE
  180   IPE(J) = L
        L = L - 1
  190 CONTINUE
      IF (LBIG.GE.ICNTL(4)) GO TO 210
C
C PLACE ROW LENGTHS IN IW
      DO 200 I = 1,N
        K = IPE(I)
        IW(K) = IQ(I)
        IF (IQ(I).EQ.0) IPE(I) = 0
  200 CONTINUE
      GO TO 250
C
C
C REMOVE DUPLICATE ENTRIES
  210 IWFR = 1
      DO 240 I = 1,N
        K1 = IPE(I) + 1
        K2 = IPE(I) + IQ(I)
        IF (K1.LE.K2) GO TO 220
        IPE(I) = 0
        GO TO 240

  220   IPE(I) = IWFR
        IWFR = IWFR + 1
        DO 230 K = K1,K2
          J = IW(K)
          IF (FLAG(J).EQ.I) GO TO 230
          IW(IWFR) = J
          IWFR = IWFR + 1
          FLAG(J) = I
  230   CONTINUE
        K = IPE(I)
        IW(K) = IWFR - K - 1
  240 CONTINUE
  250 RETURN

      END
      SUBROUTINE MA27KD(N,IPE,IW,LW,IWFR,IPS,IPV,NV,FLAG,NCMPA)
C
C USING A GIVEN PIVOTAL SEQUENCE AND A REPRESENTATION OF THE MATRIX THAT
C     INCLUDES ONLY NON-ZEROS OF THE STRICTLY UPPER-TRIANGULAR PART
C     OF THE PERMUTED MATRIX, CONSTRUCT TREE POINTERS.
C
C N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
C IPE(I) MUST BE SET TO POINT TO THE POSITION IN IW OF THE
C     START OF ROW I OR HAVE THE VALUE ZERO IF ROW I HAS NO OFF-
C     DIAGONAL NON-ZEROS. DURING EXECUTION IT IS USED AS FOLLOWS.
C     IF VARIABLE I IS ELIMINATED THEN IPE(I) POINTS TO THE LIST
C     OF VARIABLES FOR CREATED ELEMENT I. IF ELEMENT I IS
C     ABSORBED INTO NEWLY CREATED ELEMENT J THEN IPE(I)=-J.
C IW MUST BE SET ON ENTRY TO HOLD LISTS OF VARIABLES BY
C     ROWS, EACH LIST BEING HEADED BY ITS LENGTH. WHEN A VARIABLE
C     IS ELIMINATED ITS LIST IS REPLACED BY A LIST OF VARIABLES
C     IN THE NEW ELEMENT.
C LW MUST BE SET TO THE LENGTH OF IW. IT IS NOT ALTERED.
C IWFR MUST BE SET TO THE POSITION IN IW OF THE FIRST FREE VARIABLE.
C     IT IS REVISED DURING EXECUTION, CONTINUING TO HAVE THIS MEANING.
C IPS(I) MUST BE SET TO THE POSITION OF VARIABLE I IN THE REQUIRED
C     ORDERING. IT IS NOT ALTERED.
C IPV NEED NOT BE SET. IPV(K) IS SET TO HOLD THE K TH VARIABLE
C     IN PIVOT ORDER.
C NV NEED NOT BE SET. IF VARIABLE J HAS NOT BEEN ELIMINATED THEN
C     THE LAST ELEMENT WHOSE LEADING VARIABLE (VARIABLE EARLIEST
C     IN THE PIVOT SEQUENCE) IS J IS ELEMENT NV(J). IF ELEMENT J
C     EXISTS THEN THE LAST ELEMENT HAVING THE SAME LEADING
C     VARIABLE IS NV(J). IN BOTH CASES NV(J)=0 IF THERE IS NO SUCH
C     ELEMENT. IF ELEMENT J HAS BEEN MERGED INTO A LATER ELEMENT
C     THEN NV(J) IS THE DEGREE AT THE TIME OF ELIMINATION.
C FLAG IS USED AS WORKSPACE FOR VARIABLE FLAGS.
C     FLAG(JS)=ME IF JS HAS BEEN INCLUDED IN THE LIST FOR ME.
C NCMPA see INFO(11) in MA27A/AD.
C
C     .. Scalar Arguments ..
      INTEGER IWFR,LW,N,NCMPA
C     ..
C     .. Array Arguments ..
      INTEGER FLAG(N),IPE(N),IPS(N),IPV(N),IW(LW),NV(N)
C     ..
C     .. Local Scalars ..
      INTEGER I,IE,IP,J,JE,JP,JP1,JP2,JS,KDUMMY,LN,LWFR,ME,MINJS,ML,MS
C     ..
C     .. External Subroutines ..
      EXTERNAL MA27UD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
C     .. Executable Statements ..
C
C INITIALIZATIONS
      DO 10 I = 1,N
        FLAG(I) = 0
        NV(I) = 0
        J = IPS(I)
        IPV(J) = I
   10 CONTINUE
      NCMPA = 0
C
C START OF MAIN LOOP
C
      DO 100 ML = 1,N
C ME=MS IS THE NAME OF THE VARIABLE ELIMINATED AND
C     OF THE ELEMENT CREATED IN THE MAIN LOOP.
        MS = IPV(ML)
        ME = MS
        FLAG(MS) = ME
C
C MERGE ROW MS WITH ALL THE ELEMENTS HAVING MS AS LEADING VARIABLE.
C IP POINTS TO THE START OF THE NEW LIST.
        IP = IWFR
C MINJS IS SET TO THE POSITION IN THE ORDER OF THE LEADING VARIABLE
C     IN THE NEW LIST.
        MINJS = N
        IE = ME
        DO 70 KDUMMY = 1,N
C SEARCH VARIABLE LIST OF ELEMENT IE.
C JP POINTS TO THE CURRENT POSITION IN THE LIST BEING SEARCHED.
          JP = IPE(IE)
C LN IS THE LENGTH OF THE LIST BEING SEARCHED.
          LN = 0
          IF (JP.LE.0) GO TO 60
          LN = IW(JP)
C
C SEARCH FOR DIFFERENT VARIABLES AND ADD THEM TO LIST,
C     COMPRESSING WHEN NECESSARY
          DO 50 JP1 = 1,LN
            JP = JP + 1
C PLACE NEXT VARIABLE IN JS.
            JS = IW(JP)
C JUMP IF VARIABLE HAS ALREADY BEEN INCLUDED.
            IF (FLAG(JS).EQ.ME) GO TO 50
            FLAG(JS) = ME
            IF (IWFR.LT.LW) GO TO 40
C PREPARE FOR COMPRESSING IW BY ADJUSTING POINTER TO AND LENGTH OF
C     THE LIST FOR IE TO REFER TO THE REMAINING ENTRIES.
            IPE(IE) = JP
            IW(JP) = LN - JP1
C COMPRESS IW.
            CALL MA27UD(N,IPE,IW,IP-1,LWFR,NCMPA)
C COPY NEW LIST FORWARD
            JP2 = IWFR - 1
            IWFR = LWFR
            IF (IP.GT.JP2) GO TO 30
            DO 20 JP = IP,JP2
              IW(IWFR) = IW(JP)
              IWFR = IWFR + 1
   20       CONTINUE
   30       IP = LWFR
            JP = IPE(IE)
C ADD VARIABLE JS TO NEW LIST.
   40       IW(IWFR) = JS
            MINJS = MIN(MINJS,IPS(JS)+0)
            IWFR = IWFR + 1
   50     CONTINUE
C RECORD ABSORPTION OF ELEMENT IE INTO NEW ELEMENT.
   60     IPE(IE) = -ME
C PICK UP NEXT ELEMENT WITH LEADING VARIABLE MS.
          JE = NV(IE)
C STORE DEGREE OF IE.
          NV(IE) = LN + 1
          IE = JE
C LEAVE LOOP IF THERE ARE NO MORE ELEMENTS.
          IF (IE.EQ.0) GO TO 80
   70   CONTINUE
   80   IF (IWFR.GT.IP) GO TO 90
C DEAL WITH NULL NEW ELEMENT.
        IPE(ME) = 0
        NV(ME) = 1
        GO TO 100
C LINK NEW ELEMENT WITH OTHERS HAVING SAME LEADING VARIABLE.
   90   MINJS = IPV(MINJS)
        NV(ME) = NV(MINJS)
        NV(MINJS) = ME
C MOVE FIRST ENTRY IN NEW LIST TO END TO ALLOW ROOM FOR LENGTH AT
C     FRONT. SET POINTER TO FRONT.
        IW(IWFR) = IW(IP)
        IW(IP) = IWFR - IP
        IPE(ME) = IP
        IWFR = IWFR + 1
  100 CONTINUE
      RETURN

      END
      SUBROUTINE MA27LD(N,IPE,NV,IPS,NE,NA,ND,NSTEPS,NEMIN)
C
C TREE SEARCH
C
C GIVEN SON TO FATHER TREE POINTERS, PERFORM DEPTH-FIRST
C     SEARCH TO FIND PIVOT ORDER AND NUMBER OF ELIMINATIONS
C     AND ASSEMBLIES AT EACH STAGE.
C N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
C IPE(I) MUST BE SET EQUAL TO -(FATHER OF NODE I) OR ZERO IF
C      NODE IS A ROOT. IT IS ALTERED TO POINT TO ITS NEXT
C      YOUNGER BROTHER IF IT HAS ONE, BUT OTHERWISE IS NOT
C      CHANGED.
C NV(I) MUST BE SET TO ZERO IF NO VARIABLES ARE ELIMINATED AT NODE
C      I AND TO THE DEGREE OTHERWISE. ONLY LEAF NODES CAN HAVE
C      ZERO VALUES OF NV(I). NV IS NOT ALTERED.
C IPS(I) NEED NOT BE SET. IT IS USED TEMPORARILY TO HOLD
C      -(ELDEST SON OF NODE I) IF IT HAS ONE AND 0 OTHERWISE. IT IS
C      EVENTUALLY SET TO HOLD THE POSITION OF NODE I IN THE ORDER.
C NE(IS) NEED NOT BE SET. IT IS SET TO THE NUMBER OF VARIABLES
C      ELIMINATED AT STAGE IS OF THE ELIMINATION.
C NA(IS) NEED NOT BE SET. IT IS SET TO THE NUMBER OF ELEMENTS
C      ASSEMBLED AT STAGE IS OF THE ELIMINATION.
C ND(IS) NEED NOT BE SET. IT IS SET TO THE DEGREE AT STAGE IS OF
C     THE ELIMINATION.
C NSTEPS NEED NOT BE SET. IT IS SET TO  THE NUMBER OF ELIMINATION
C      STEPS.
C NEMIN see ICNTL(5) in MA27A/AD.
C
C     .. Scalar Arguments ..
      INTEGER N,NSTEPS,NEMIN
C     ..
C     .. Array Arguments ..
      INTEGER IPE(N),IPS(N),NA(N),ND(N),NE(N),NV(N)
C     ..
C     .. Local Scalars ..
      INTEGER I,IB,IF,IL,IS,ISON,K,L,NR
C     ..
C     .. Executable Statements ..
C INITIALIZE IPS AND NE.
      DO 10 I = 1,N
        IPS(I) = 0
        NE(I) = 0
   10 CONTINUE
C
C SET IPS(I) TO -(ELDEST SON OF NODE I) AND IPE(I) TO NEXT YOUNGER
C     BROTHER OF NODE I IF IT HAS ONE.
C FIRST PASS IS FOR NODES WITHOUT ELIMINATIONS.
      DO 20 I = 1,N
        IF (NV(I).GT.0) GO TO 20
        IF = -IPE(I)
        IS = -IPS(IF)
        IF (IS.GT.0) IPE(I) = IS
        IPS(IF) = -I
   20 CONTINUE
C NR IS DECREMENTED FOR EACH ROOT NODE. THESE ARE STORED IN
C     NE(I),I=NR,N.
      NR = N + 1
C SECOND PASS TO ADD NODES WITH ELIMINATIONS.
      DO 50 I = 1,N
        IF (NV(I).LE.0) GO TO 50
C NODE IF IS THE FATHER OF NODE I.
        IF = -IPE(I)
        IF (IF.EQ.0) GO TO 40
        IS = -IPS(IF)
C JUMP IF NODE IF HAS NO SONS YET.
        IF (IS.LE.0) GO TO 30
C SET POINTER TO NEXT BROTHER
        IPE(I) = IS
C NODE I IS ELDEST SON OF NODE IF.
   30   IPS(IF) = -I
        GO TO 50
C WE HAVE A ROOT
   40   NR = NR - 1
        NE(NR) = I
   50 CONTINUE
C
C DEPTH-FIRST SEARCH.
C IL HOLDS THE CURRENT TREE LEVEL. ROOTS ARE AT LEVEL N, THEIR SONS
C     ARE AT LEVEL N-1, ETC.
C IS HOLDS THE CURRENT ELIMINATION STAGE. WE ACCUMULATE THE NUMBER
C     OF ELIMINATIONS AT STAGE IS DIRECTLY IN NE(IS). THE NUMBER OF
C     ASSEMBLIES IS ACCUMULATED TEMPORARILY IN NA(IL), FOR TREE
C     LEVEL IL, AND IS TRANSFERED TO NA(IS) WHEN WE REACH THE
C     APPROPRIATE STAGE IS.
      IS = 1
C I IS THE CURRENT NODE.
      I = 0
      DO 160 K = 1,N
        IF (I.GT.0) GO TO 60
C PICK UP NEXT ROOT.
        I = NE(NR)
        NE(NR) = 0
        NR = NR + 1
        IL = N
        NA(N) = 0
C GO TO SON FOR AS LONG AS POSSIBLE, CLEARING FATHER-SON POINTERS
C     IN IPS AS EACH IS USED AND SETTING NA(IL)=0 FOR ALL LEVELS
C     REACHED.
   60   DO 70 L = 1,N
          IF (IPS(I).GE.0) GO TO 80
          ISON = -IPS(I)
          IPS(I) = 0
          I = ISON
          IL = IL - 1
          NA(IL) = 0
   70   CONTINUE
C RECORD POSITION OF NODE I IN THE ORDER.
   80   IPS(I) = K
        NE(IS) = NE(IS) + 1
C JUMP IF NODE HAS NO ELIMINATIONS.
        IF (NV(I).LE.0) GO TO 120
        IF (IL.LT.N) NA(IL+1) = NA(IL+1) + 1
        NA(IS) = NA(IL)
        ND(IS) = NV(I)
C CHECK FOR STATIC CONDENSATION
        IF (NA(IS).NE.1) GO TO 90
        IF (ND(IS-1)-NE(IS-1).EQ.ND(IS)) GO TO 100
C CHECK FOR SMALL NUMBERS OF ELIMINATIONS IN BOTH LAST TWO STEPS.
   90   IF (NE(IS).GE.NEMIN) GO TO 110
        IF (NA(IS).EQ.0) GO TO 110
        IF (NE(IS-1).GE.NEMIN) GO TO 110
C COMBINE THE LAST TWO STEPS
  100   NA(IS-1) = NA(IS-1) + NA(IS) - 1
        ND(IS-1) = ND(IS) + NE(IS-1)
        NE(IS-1) = NE(IS) + NE(IS-1)
        NE(IS) = 0
        GO TO 120

  110   IS = IS + 1
  120   IB = IPE(I)
        IF (IB.GE.0) THEN
C NODE I HAS A BROTHER OR IS A ROOT
          IF (IB.GT.0) NA(IL) = 0
          I = IB
        ELSE
C GO TO FATHER OF NODE I
          I = -IB
          IL = IL + 1
        END IF
  160 CONTINUE
      NSTEPS = IS - 1
      RETURN

      END
      SUBROUTINE MA27MD(N,NZ,IRN,ICN,PERM,NA,NE,ND,NSTEPS,LSTKI,LSTKR,
     +                 IW,INFO,OPS)
C
C STORAGE AND OPERATION COUNT EVALUATION.
C
C EVALUATE NUMBER OF OPERATIONS AND SPACE REQUIRED BY FACTORIZATION
C     USING MA27B/BD.  THE VALUES GIVEN ARE EXACT ONLY IF NO NUMERICAL
C     PIVOTING IS PERFORMED AND THEN ONLY IF IRN(1) WAS NOT
C     EQUIVALENCED TO IW(1) BY THE USER BEFORE CALLING MA27A/AD.  IF
C     THE EQUIVALENCE HAS BEEN MADE ONLY AN UPPER BOUND FOR NIRNEC
C     AND NRLNEC CAN BE CALCULATED ALTHOUGH THE OTHER COUNTS WILL
C     STILL BE EXACT.
C
C N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
C NZ MUST BE SET TO THE NUMBER OF NON-ZEROS INPUT.  IT IS NOT ALTERED.
C IRN,ICN.  UNLESS IRN(1) HAS BEEN EQUIVALENCED TO IW(1)
C     IRN,ICN MUST BE SET TO THE ROW AND COLUMN INDICES OF THE
C     NON-ZEROS ON INPUT.  THEY ARE NOT ALTERED BY MA27M/MD.
C PERM MUST BE SET TO THE POSITION IN THE PIVOT ORDER OF EACH ROW.
C     IT IS NOT ALTERED.
C NA,NE,ND MUST BE SET TO HOLD, FOR EACH TREE NODE, THE NUMBER OF STACK
C     ELEMENTS ASSEMBLED, THE NUMBER OF ELIMINATIONS AND THE SIZE OF
C     THE ASSEMBLED FRONT MATRIX RESPECTIVELY.  THEY ARE NOT ALTERED.
C NSTEPS MUST BE SET TO HOLD THE NUMBER OF TREE NODES. IT IS NOT
C     ALTERED.
C LSTKI IS USED AS A WORK ARRAY BY MA27M/MD.
C LSTKR.  IF IRN(1) IS EQUIVALENCED TO IW(1)  THEN LSTKR(I)
C     MUST HOLD THE TOTAL NUMBER OF OFF-DIAGONAL ENTRIES (INCLUDING
C     DUPLICATES) IN ROW I (I=1,..,N) OF THE ORIGINAL MATRIX.  IT
C     IS USED AS WORKSPACE BY MA27M/MD.
C IW IS A WORKSPACE ARRAY USED BY OTHER SUBROUTINES AND PASSED TO THIS
C     SUBROUTINE ONLY SO THAT A TEST FOR EQUIVALENCE WITH IRN CAN BE
C     MADE.
C
C COUNTS FOR OPERATIONS AND STORAGE ARE ACCUMULATED IN VARIABLES
C     OPS,NRLTOT,NIRTOT,NRLNEC,NIRNEC,NRLADU,NRLNEC,NIRNEC.
C OPS NUMBER OF MULTIPLICATIONS AND ADDITIONS DURING FACTORIZATION.
C NRLADU,NIRADU REAL AND INTEGER STORAGE RESPECTIVELY FOR THE
C     MATRIX FACTORS.
C NRLTOT,NIRTOT REAL AND INTEGER STRORAGE RESPECTIVELY REQUIRED
C     FOR THE FACTORIZATION IF NO COMPRESSES ARE ALLOWED.
C NRLNEC,NIRNEC REAL AND INTEGER STORAGE RESPECTIVELY REQUIRED FOR
C     THE FACTORIZATION IF COMPRESSES ARE ALLOWED.
C INFO is an INTEGER array of length 20, see MA27A/AD.
C OPS ACCUMULATES THE NO. OF MULTIPLY/ADD PAIRS NEEDED TO CREATE THE
C     TRIANGULAR FACTORIZATION, IN THE DEFINITE CASE.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION OPS
      INTEGER N,NSTEPS,NZ
C     ..
C     .. Array Arguments ..
      INTEGER ICN(*),IRN(*),IW(*),LSTKI(N),LSTKR(N),NA(NSTEPS),
     +        ND(NSTEPS),NE(NSTEPS),PERM(N),INFO(20)
C     ..
C     .. Local Scalars ..
      INTEGER I,INEW,IOLD,IORG,IROW,ISTKI,ISTKR,ITOP,ITREE,JOLD,JORG,K,
     +        LSTK,NASSR,NELIM,NFR,NSTK,NUMORG,NZ1,NZ2
      DOUBLE PRECISION DELIM
      INTEGER NRLADU,NIRADU,NIRTOT,NRLTOT,NIRNEC,NRLNEC
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX,MIN
C     ..
C     .. Executable Statements ..
C
      IF (NZ.EQ.0) GO TO 20
C JUMP IF IW AND IRN HAVE NOT BEEN EQUIVALENCED.
      IF (IRN(1).NE.IW(1)) GO TO 20
C RESET IRN(1).
      IRN(1) = IW(1) - 1
C THE TOTAL NUMBER OF OFF-DIAGONAL ENTRIES IS ACCUMULATED IN NZ2.
C LSTKI IS SET TO HOLD THE TOTAL NUMBER OF ENTRIES (INCUDING
C     THE DIAGONAL) IN EACH ROW IN PERMUTED ORDER.
      NZ2 = 0
      DO 10 IOLD = 1,N
        INEW = PERM(IOLD)
        LSTKI(INEW) = LSTKR(IOLD) + 1
        NZ2 = NZ2 + LSTKR(IOLD)
   10 CONTINUE
C NZ1 IS THE NUMBER OF ENTRIES IN ONE TRIANGLE INCLUDING THE DIAGONAL.
C NZ2 IS THE TOTAL NUMBER OF ENTRIES INCLUDING THE DIAGONAL.
      NZ1 = NZ2/2 + N
      NZ2 = NZ2 + N
      GO TO 60
C COUNT (IN LSTKI) NON-ZEROS IN ORIGINAL MATRIX BY PERMUTED ROW (UPPER
C     TRIANGLE ONLY). INITIALIZE COUNTS.
   20 DO 30 I = 1,N
        LSTKI(I) = 1
   30 CONTINUE
C ACCUMULATE NUMBER OF NON-ZEROS WITH INDICES IN RANGE IN NZ1
C     DUPLICATES ON THE DIAGONAL ARE IGNORED BUT NZ1 INCLUDES ANY
C     DIAGONALS NOT PRESENT ON INPUT.
C ACCUMULATE ROW COUNTS IN LSTKI.
      NZ1 = N
      IF (NZ.EQ.0) GO TO 50
      DO 40 I = 1,NZ
        IOLD = IRN(I)
        JOLD = ICN(I)
C JUMP IF INDEX IS OUT OF RANGE.
        IF (IOLD.LT.1 .OR. IOLD.GT.N) GO TO 40
        IF (JOLD.LT.1 .OR. JOLD.GT.N) GO TO 40
        IF (IOLD.EQ.JOLD) GO TO 40
        NZ1 = NZ1 + 1
        IROW = MIN(PERM(IOLD)+0,PERM(JOLD)+0)
        LSTKI(IROW) = LSTKI(IROW) + 1
   40 CONTINUE
   50 NZ2 = NZ1
C ISTKR,ISTKI CURRENT NUMBER OF STACK ENTRIES IN
C     REAL AND INTEGER STORAGE RESPECTIVELY.
C OPS,NRLADU,NIRADU,NIRTOT,NRLTOT,NIRNEC,NRLNEC,NZ2 ARE DEFINED ABOVE.
C NZ2 CURRENT NUMBER OF ORIGINAL MATRIX ENTRIES NOT YET PROCESSED.
C NUMORG CURRENT TOTAL NUMBER OF ROWS ELIMINATED.
C ITOP CURRENT NUMBER OF ELEMENTS ON THE STACK.
   60 ISTKI = 0
      ISTKR = 0
      OPS = 0.0D0
      NRLADU = 0
C     ONE LOCATION IS NEEDED TO RECORD THE NUMBER OF BLOCKS
C     ACTUALLY USED.
      NIRADU = 1
      NIRTOT = NZ1
      NRLTOT = NZ1
      NIRNEC = NZ2
      NRLNEC = NZ2
      NUMORG = 0
      ITOP = 0
C
C EACH PASS THROUGH THIS LOOP PROCESSES A NODE OF THE TREE.
      DO 100 ITREE = 1,NSTEPS
        NELIM = NE(ITREE)
        DELIM = NELIM
        NFR = ND(ITREE)
        NSTK = NA(ITREE)
C ADJUST STORAGE COUNTS ON ASSEMBLY OF CURRENT FRONTAL MATRIX.
        NASSR = NFR* (NFR+1)/2
        IF (NSTK.NE.0) NASSR = NASSR - LSTKR(ITOP) + 1
        NRLTOT = MAX(NRLTOT,NRLADU+NASSR+ISTKR+NZ1)
        NIRTOT = MAX(NIRTOT,NIRADU+NFR+2+ISTKI+NZ1)
        NRLNEC = MAX(NRLNEC,NRLADU+NASSR+ISTKR+NZ2)
        NIRNEC = MAX(NIRNEC,NIRADU+NFR+2+ISTKI+NZ2)
C DECREASE NZ2 BY THE NUMBER OF ENTRIES IN ROWS BEING ELIMINATED AT
C     THIS STAGE.
        DO 70 IORG = 1,NELIM
          JORG = NUMORG + IORG
          NZ2 = NZ2 - LSTKI(JORG)
   70   CONTINUE
        NUMORG = NUMORG + NELIM
C JUMP IF THERE ARE NO STACK ASSEMBLIES AT THIS NODE.
        IF (NSTK.LE.0) GO TO 90
C REMOVE ELEMENTS FROM THE STACK.  THERE ARE ITOP ELEMENTS ON THE
C     STACK WITH THE APPROPRIATE ENTRIES IN LSTKR,LSTKI GIVING
C     THE REAL AND INTEGER STORAGE RESPECTIVELY FOR EACH STACK
C     ELEMENT.
        DO 80 K = 1,NSTK
          LSTK = LSTKR(ITOP)
          ISTKR = ISTKR - LSTK
          LSTK = LSTKI(ITOP)
          ISTKI = ISTKI - LSTK
          ITOP = ITOP - 1
   80   CONTINUE
C ACCUMULATE NON-ZEROS IN FACTORS AND NUMBER OF OPERATIONS.
   90   NRLADU = NRLADU + (NELIM* (2*NFR-NELIM+1))/2
        NIRADU = NIRADU + 2 + NFR
        IF (NELIM.EQ.1) NIRADU = NIRADU - 1
        OPS = OPS + ((NFR*DELIM*(NFR+1)-(2*NFR+1)*DELIM*(DELIM+1)/2+
     +        DELIM* (DELIM+1)* (2*DELIM+1)/6)/2)
        IF (ITREE.EQ.NSTEPS) GO TO 100
C JUMP IF ALL OF FRONTAL MATRIX HAS BEEN ELIMINATED.
        IF (NFR.EQ.NELIM) GO TO 100
C STACK REMAINDER OF ELEMENT.
        ITOP = ITOP + 1
        LSTKR(ITOP) = (NFR-NELIM)* (NFR-NELIM+1)/2
        LSTKI(ITOP) = NFR - NELIM + 1
        ISTKI = ISTKI + LSTKI(ITOP)
        ISTKR = ISTKR + LSTKR(ITOP)
C WE DO NOT NEED TO ADJUST THE COUNTS FOR THE REAL STORAGE BECAUSE
C     THE REMAINDER OF THE FRONTAL MATRIX IS SIMPLY MOVED IN THE
C     STORAGE FROM FACTORS TO STACK AND NO EXTRA STORAGE IS REQUIRED.
        NIRTOT = MAX(NIRTOT,NIRADU+ISTKI+NZ1)
        NIRNEC = MAX(NIRNEC,NIRADU+ISTKI+NZ2)
  100 CONTINUE
C
C ADJUST THE STORAGE COUNTS TO ALLOW FOR THE USE OF THE REAL AND
C     INTEGER STORAGE FOR PURPOSES OTHER THAN PURELY THE
C     FACTORIZATION ITSELF.
C THE SECOND TWO TERMS ARE THE MINUMUM AMOUNT REQUIRED BY MA27N/ND.
      NRLNEC = MAX(NRLNEC,N+MAX(NZ,NZ1))
      NRLTOT = MAX(NRLTOT,N+MAX(NZ,NZ1))
      NRLNEC = MIN(NRLNEC,NRLTOT)
      NIRNEC = MAX(NZ,NIRNEC)
      NIRTOT = MAX(NZ,NIRTOT)
      NIRNEC = MIN(NIRNEC,NIRTOT)

      INFO(3) = NRLTOT
      INFO(4) = NIRTOT
      INFO(5) = NRLNEC
      INFO(6) = NIRNEC
      INFO(7) = NRLADU
      INFO(8) = NIRADU
      RETURN

      END 
      SUBROUTINE MA27ND(N,NZ,NZ1,A,LA,IRN,ICN,IW,LIW,PERM,IW2,ICNTL,
     +                 INFO)
C
C SORT PRIOR TO FACTORIZATION USING MA27O/OD.
C
C THIS SUBROUTINE REORDERS THE USER'S INPUT SO THAT THE UPPER TRIANGLE
C     OF THE PERMUTED MATRIX, INCLUDING THE DIAGONAL, IS
C     HELD ORDERED BY ROWS AT THE END OF THE STORAGE FOR A AND IW.
C     IT IGNORES ENTRIES WITH ONE OR BOTH INDICES OUT OF RANGE AND
C     ACCUMULATES DIAGONAL ENTRIES.
C     IT ADDS EXPLICIT ZEROS ON THE DIAGONAL WHERE NECESSARY.
C N      - MUST BE SET TO THE ORDER OF THE MATRIX.
C          IT IS NOT ALTERED BY MA27N/ND.
C NZ     - ON ENTRY NZ MUST BE SET TO THE NUMBER
C          OF NON-ZEROS INPUT.  NOT ALTERED BY MA27N/ND.
C NZ1    - ON EXIT NZ1 WILL BE EQUAL TO THE NUMBER OF ENTRIES IN THE
C          SORTED MATRIX.
C A      - ON ENTRY A(I) MUST
C          HOLD THE VALUE OF THE ORIGINAL MATRIX ELEMENT IN POSITION
C          (IRN(I),ICN(I)),I=1,NZ.  ON EXIT A(LA-NZ1+I),I=1,NZ1 HOLDS
C          THE UPPER TRIANGLE OF THE PERMUTED MATRIX BY ROWS WITH
C          THE DIAGONAL ENTRY FIRST ALTHOUGH THERE IS NO FURTHER
C          ORDERING WITHIN THE ROWS THEMSELVES.
C LA     - LENGTH OF ARRAY A. MUST BE AT LEAST N+MAX(NZ,NZ1).
C          IT IS NOT ALTERED BY MA27N/ND.
C IRN    - IRN(I) MUST BE SET TO
C          HOLD THE ROW INDEX OF ENTRY A(I),I=1,NZ.  IRN WILL BE
C          UNALTERED BY MA27N/ND, UNLESS IT IS EQUIVALENCED WITH IW.
C ICN    - ICN(I) MUST BE SET TO
C          HOLD THE COLUMN INDEX OF ENTRY A(I),I=1,NZ.  ICN WILL BE
C          UNALTERED BY MA27N/ND, UNLESS IT IS EQUIVALENCED WITH IW.
C IW     - USED AS WORKSPACE AND ON
C          EXIT, ENTRIES IW(LIW-NZ1+I),I=1,NZ1 HOLD THE COLUMN INDICES
C          (THE ORIGINAL UNPERMUTED INDICES) OF THE CORRESPONDING ENTRY
C          OF A WITH THE FIRST ENTRY FOR EACH ROW FLAGGED NEGATIVE.
C          IRN(1) MAY BE EQUIVALENCED WITH IW(1) AND ICN(1) MAY BE
C          EQUIVALENCED WITH IW(K) WHERE K.GT.NZ.
C LIW    - LENGTH OF ARRAY IW. MUST BE AT LEAST AS
C          GREAT AS THE MAXIMUM OF NZ AND NZ1.
C          NOT ALTERED BY MA27N/ND.
C PERM   - PERM(I) HOLDS THE
C          POSITION IN THE TENTATIVE PIVOT ORDER OF ROW I IN THE
C          ORIGINAL MATRIX,I=1,N. NOT ALTERED BY MA27N/ND.
C IW2    - USED AS WORKSPACE.
C          SEE COMMENTS IN CODE IMMEDIATELY PRIOR TO
C          EACH USE.
C ICNTL is an INTEGER array of length 30, see MA27A/AD.
C INFO is an INTEGER array of length 20, see MA27A/AD.
C   INFO(1)  - ON EXIT FROM MA27N/ND, A ZERO VALUE OF
C          INFO(1) INDICATES THAT NO ERROR HAS BEEN DETECTED.
C          POSSIBLE NON-ZERO VALUES ARE ..
C          +1  WARNING.  INDICES OUT OF RANGE.  THESE ARE IGNORED,
C              THEIR NUMBER IS RECORDED IN INFO(2) OF MA27E/ED AND
C              MESSAGES IDENTIFYING THE FIRST TEN ARE OUTPUT ON UNIT
C              ICNTL(2).
C          -3  INTEGER ARRAY IW IS TOO SMALL.
C          -4  DOUBLE PRECISION ARRAY A IS TOO SMALL.
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER LA,LIW,N,NZ,NZ1
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA)
      INTEGER ICN(*),IRN(*),IW(LIW),IW2(N),PERM(N),ICNTL(30),INFO(20)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ANEXT,ANOW
      INTEGER I,IA,ICH,II,IIW,INEW,IOLD,IPOS,J1,J2,JJ,JNEW,JOLD,JPOS,K
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
C     .. Executable Statements ..
      INFO(1) = 0
C INITIALIZE WORK ARRAY (IW2) IN PREPARATION FOR
C     COUNTING NUMBERS OF NON-ZEROS IN THE ROWS AND INITIALIZE
C     LAST N ENTRIES IN A WHICH WILL HOLD THE DIAGONAL ENTRIES
      IA = LA
      DO 10 IOLD = 1,N
        IW2(IOLD) = 1
        A(IA) = ZERO
        IA = IA - 1
   10 CONTINUE
C SCAN INPUT COPYING ROW INDICES FROM IRN TO THE FIRST NZ POSITIONS
C     IN IW.  THE NEGATIVE OF THE INDEX IS HELD TO FLAG ENTRIES FOR
C     THE IN-PLACE SORT.  ENTRIES IN IW CORRESPONDING TO DIAGONALS AND
C     ENTRIES WITH OUT-OF-RANGE INDICES ARE SET TO ZERO.
C     FOR DIAGONAL ENTRIES, REALS ARE ACCUMULATED IN THE LAST N
C     LOCATIONS OF A.
C     THE NUMBER OF ENTRIES IN EACH ROW OF THE PERMUTED MATRIX IS
C     ACCUMULATED IN IW2.
C INDICES OUT OF RANGE ARE IGNORED  AFTER BEING COUNTED AND
C     AFTER APPROPRIATE MESSAGES HAVE BEEN PRINTED.
      INFO(2) = 0
C NZ1 IS THE NUMBER OF NON-ZEROS HELD AFTER INDICES OUT OF RANGE HAVE
C     BEEN IGNORED AND DIAGONAL ENTRIES ACCUMULATED.
      NZ1 = N
      IF (NZ.EQ.0) GO TO 80
      DO 70 K = 1,NZ
        IOLD = IRN(K)
        IF (IOLD.GT.N .OR. IOLD.LE.0) GO TO 30
        JOLD = ICN(K)
        IF (JOLD.GT.N .OR. JOLD.LE.0) GO TO 30
        INEW = PERM(IOLD)
        JNEW = PERM(JOLD)
        IF (INEW.NE.JNEW) GO TO 20
        IA = LA - N + IOLD
        A(IA) = A(IA) + A(K)
        GO TO 60

   20   INEW = MIN(INEW,JNEW)
C INCREMENT NUMBER OF ENTRIES IN ROW INEW.
        IW2(INEW) = IW2(INEW) + 1
        IW(K) = -IOLD
        NZ1 = NZ1 + 1
        GO TO 70
C ENTRY OUT OF RANGE.  IT WILL BE IGNORED AND A FLAG SET.
   30   INFO(1) = 1
        INFO(2) = INFO(2) + 1
        IF (INFO(2).LE.1 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=40) INFO(1)
        ENDIF

   40   FORMAT (' *** WARNING MESSAGE FROM SUBROUTINE MA27BD',
     +          '  *** INFO(1) =',I2)

        IF (INFO(2).LE.10 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=50) K,IRN(K),ICN(K)
        END IF

   50   FORMAT (I6,'TH NON-ZERO (IN ROW',I6,' AND COLUMN',I6,
     +         ') IGNORED')

   60   IW(K) = 0
   70 CONTINUE
C CALCULATE POINTERS (IN IW2) TO THE POSITION IMMEDIATELY AFTER THE END
C     OF EACH ROW.
   80 IF (NZ.LT.NZ1 .AND. NZ1.NE.N) GO TO 100
C ROOM IS INCLUDED FOR THE DIAGONALS.
      K = 1
      DO 90 I = 1,N
        K = K + IW2(I)
        IW2(I) = K
   90 CONTINUE
      GO TO 120
C ROOM IS NOT INCLUDED FOR THE DIAGONALS.
  100 K = 1
      DO 110 I = 1,N
        K = K + IW2(I) - 1
        IW2(I) = K
  110 CONTINUE
C FAIL IF INSUFFICIENT SPACE IN ARRAYS A OR IW.
  120 IF (NZ1.GT.LIW) GO TO 210
      IF (NZ1+N.GT.LA) GO TO 220
C NOW RUN THROUGH NON-ZEROS IN ORDER PLACING THEM IN THEIR NEW
C POSITION AND DECREMENTING APPROPRIATE IW2 ENTRY.  IF WE ARE
C ABOUT TO OVERWRITE AN ENTRY NOT YET MOVED, WE MUST DEAL WITH
C THIS AT THIS TIME.
      IF (NZ1.EQ.N) GO TO 180
      DO 140 K = 1,NZ
        IOLD = -IW(K)
        IF (IOLD.LE.0) GO TO 140
        JOLD = ICN(K)
        ANOW = A(K)
        IW(K) = 0
        DO 130 ICH = 1,NZ
          INEW = PERM(IOLD)
          JNEW = PERM(JOLD)
          INEW = MIN(INEW,JNEW)
          IF (INEW.EQ.PERM(JOLD)) JOLD = IOLD
          JPOS = IW2(INEW) - 1
          IOLD = -IW(JPOS)
          ANEXT = A(JPOS)
          A(JPOS) = ANOW
          IW(JPOS) = JOLD
          IW2(INEW) = JPOS
          IF (IOLD.EQ.0) GO TO 140
          ANOW = ANEXT
          JOLD = ICN(JPOS)
  130   CONTINUE
  140 CONTINUE
      IF (NZ.GE.NZ1) GO TO 180
C MOVE UP ENTRIES TO ALLOW FOR DIAGONALS.
      IPOS = NZ1
      JPOS = NZ1 - N
      DO 170 II = 1,N
        I = N - II + 1
        J1 = IW2(I)
        J2 = JPOS
        IF (J1.GT.JPOS) GO TO 160
        DO 150 JJ = J1,J2
          IW(IPOS) = IW(JPOS)
          A(IPOS) = A(JPOS)
          IPOS = IPOS - 1
          JPOS = JPOS - 1
  150   CONTINUE
  160   IW2(I) = IPOS + 1
        IPOS = IPOS - 1
  170 CONTINUE
C RUN THROUGH ROWS INSERTING DIAGONAL ENTRIES AND FLAGGING BEGINNING
C     OF EACH ROW BY NEGATING FIRST COLUMN INDEX.
  180 DO 190 IOLD = 1,N
        INEW = PERM(IOLD)
        JPOS = IW2(INEW) - 1
        IA = LA - N + IOLD
        A(JPOS) = A(IA)
        IW(JPOS) = -IOLD
  190 CONTINUE
C MOVE SORTED MATRIX TO THE END OF THE ARRAYS.
      IPOS = NZ1
      IA = LA
      IIW = LIW
      DO 200 I = 1,NZ1
        A(IA) = A(IPOS)
        IW(IIW) = IW(IPOS)
        IPOS = IPOS - 1
        IA = IA - 1
        IIW = IIW - 1
  200 CONTINUE
      GO TO 230
C **** ERROR RETURN ****
  210 INFO(1) = -3
      INFO(2) = NZ1
      GO TO 230

  220 INFO(1) = -4
      INFO(2) = NZ1 + N
C
  230 RETURN

      END
      SUBROUTINE MA27OD(N,NZ,A,LA,IW,LIW,PERM,NSTK,NSTEPS,MAXFRT,NELIM,
     +                 IW2,ICNTL,CNTL,INFO)
C
C FACTORIZATION SUBROUTINE
C
C THIS SUBROUTINE OPERATES ON THE INPUT MATRIX ORDERED BY MA27N/ND AND
C     PRODUCES THE FACTORS OF U AND D ('A'=UTRANSPOSE*D*U) FOR USE IN
C     SUBSEQUENT SOLUTIONS.  GAUSSIAN ELIMINATION IS USED WITH PIVOTS
C     CHOSEN FROM THE DIAGONAL.  TO ENSURE STABILITY, BLOCK PIVOTS OF
C     ORDER 2 WILL BE USED IF THE DIAGONAL ENTRY IS NOT LARGE ENOUGH.
C
C N      - MUST BE SET TO THE ORDER OF THE MATRIX. IT IS NOT ALTERED.
C NZ     - MUST BE SET TO THE NUMBER OF NON-ZEROS IN UPPER TRIANGLE OF
C          PERMUTED MATRIX.  NOT ALTERED BY MA27O/OD.
C A      - MUST BE SET ON INPUT TO MATRIX HELD BY ROWS REORDERED BY
C          PERMUTATION FROM MA27A/AD IN A(LA-NZ+I),I=1,NZ.   ON
C          EXIT FROM MA27O/OD, THE FACTORS OF U AND D ARE HELD IN
C          POSITIONS 1 TO POSFAC-1.
C LA     - LENGTH OF ARRAY A.  A VALUE FOR LA
C          SUFFICIENT FOR DEFINITE SYSTEMS
C          WILL HAVE BEEN PROVIDED BY MA27A/AD. NOT ALTERED BY MA27O/OD.
C IW     - MUST BE SET SO THAT,ON INPUT, IW(LIW-NZ+I),I=1,NZ
C          HOLDS THE COLUMN INDEX OF THE ENTRY IN A(LA-NZ+I).  ON EXIT,
C          IW HOLDS INTEGER INDEXING INFORMATION ON THE FACTORS.
C          THE ABSOLUTE VALUE OF THE FIRST ENTRY IN IW WILL BE SET TO
C          THE NUMBER OF BLOCK PIVOTS ACTUALLY USED.  THIS MAY BE
C          DIFFERENT FROM NSTEPS SINCE NUMERICAL CONSIDERATIONS
C          MAY PREVENT US CHOOSING A PIVOT AT EACH STAGE.  IF THIS ENTRY
C          IN IW IS NEGATIVE, THEN AT LEAST ONE TWO BY TWO
C          PIVOT HAS BEEN USED DURING THE DECOMPOSITION.
C          INTEGER INFORMATION ON EACH BLOCK PIVOT ROW FOLLOWS.  FOR
C          EACH BLOCK PIVOT ROW THE COLUMN INDICES ARE PRECEDED BY A
C          COUNT OF THE NUMBER OF ROWS AND COLUMNS IN THE BLOCK PIVOT
C          WHERE, IF ONLY ONE ROW IS PRESENT, ONLY THE NUMBER OF
C          COLUMNS TOGETHER WITH A NEGATIVE FLAG IS HELD.  THE FIRST
C          COLUMN INDEX FOR A TWO BY TWO PIVOT IS FLAGGED NEGATIVE.
C LIW    - LENGTH OF ARRAY IW.  A VALUE FOR LIW SUFFICIENT FOR
C          DEFINITE SYSTEMS
C          WILL HAVE BEEN PROVIDED BY MA27A/AD.  NOT ALTERED BY MA27O/OD
C PERM   - PERM(I) MUST BE SET TO POSITION OF ROW/COLUMN I IN THE
C          TENTATIVE PIVOT ORDER GENERATED BY MA27A/AD.
C          IT IS NOT ALTERED BY MA27O/OD.
C NSTK   - MUST BE LEFT UNCHANGED SINCE OUTPUT FROM MA27A/AD. NSTK(I)
C          GIVES THE NUMBER OF GENERATED STACK ELEMENTS ASSEMBLED AT
C          STAGE I.  IT IS NOT ALTERED BY MA27O/OD.
C NSTEPS - LENGTH OF ARRAYS NSTK AND NELIM. VALUE IS GIVEN ON OUTPUT
C          FROM MA27A/AD (WILL NEVER EXCEED N). IT IS NOT ALTERED BY
C          MA27O/OD.
C MAXFRT - NEED NOT BE SET ON INPUT.  ON OUTPUT
C          MAXFRT WILL BE SET TO THE MAXIMUM FRONT SIZE ENCOUNTERED
C          DURING THE DECOMPOSITION.
C NELIM  - MUST BE UNCHANGED SINCE OUTPUT FROM MA27A/AD. NELIM(I)
C          GIVES THE NUMBER OF ORIGINAL ROWS ASSEMBLED AT STAGE I.
C          IT IS NOT ALTERED BY MA27O/OD.
C IW2    - INTEGER ARRAY OF LENGTH N. USED AS WORKSPACE BY MA27O/OD.
C          ALTHOUGH WE COULD HAVE USED A SHORT WORD INTEGER IN THE IBM
C          VERSION, WE HAVE NOT DONE SO BECAUSE THERE IS A SPARE
C          FULL INTEGER ARRAY (USED IN THE SORT BEFORE MA27O/OD)
C          AVAILABLE WHEN MA27O/OD IS CALLED FROM MA27B/BD.
C ICNTL is an INTEGER array of length 30, see MA27A/AD.
C CNTL is a DOUBLE PRECISION array of length 5, see MA27A/AD.
C INFO is an INTEGER array of length 20, see MA27A/AD.
C   INFO(1)  - INTEGER VARIABLE.  DIAGNOSTIC FLAG.  A ZERO VALUE ON EXIT
C          INDICATES SUCCESS.  POSSIBLE NEGATIVE VALUES ARE ...
C          -3  INSUFFICIENT STORAGE FOR IW.
C          -4  INSUFFICIENT STORAGE FOR A.
C          -5  ZERO PIVOT FOUND IN FACTORIZATION OF DEFINITE MATRIX.
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO,HALF,ONE
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER LA,LIW,MAXFRT,N,NSTEPS,NZ
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA),CNTL(5)
      INTEGER IW(LIW),IW2(N),NELIM(NSTEPS),NSTK(NSTEPS),PERM(N)
      INTEGER ICNTL(30),INFO(20)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AMAX,AMULT,AMULT1,AMULT2,DETPIV,RMAX,SWOP,
     +        THRESH,TMAX,UU
      INTEGER AINPUT,APOS,APOS1,APOS2,APOS3,ASTK,ASTK2,AZERO,I,IASS,
     +        IBEG,IDUMMY,IELL,IEND,IEXCH,IFR,IINPUT,IOLDPS,IORG,IPIV,
     +        IPMNP,IPOS,IROW,ISNPIV,ISTK,ISTK2,ISWOP,IWPOS,IX,IY,J,J1,
     +        J2,JCOL,JDUMMY,JFIRST,JJ,JJ1,JJJ,JLAST,JMAX,JMXMIP,JNEW,
     +        JNEXT,JPIV,JPOS,K,K1,K2,KDUMMY,KK,KMAX,KROW,LAELL,LAPOS2,
     +        LIELL,LNASS,LNPIV,LT,LTOPST,NASS,NBLK,NEWEL,NFRONT,NPIV,
     +        NPIVP1,NTOTPV,NUMASS,NUMORG,NUMSTK,PIVSIZ,POSFAC,POSPV1,
     +        POSPV2
      INTEGER NTWO,NEIG,NCMPBI,NCMPBR,NRLBDU,NIRBDU
C     ..
C     .. External Subroutines ..
      EXTERNAL MA27PD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN
C     ..
C     .. Statement Functions ..
      INTEGER IDIAG
C     ..
C     .. Statement Function definitions ..
C THE FOLLOWING ARITHMETIC FUNCTION GIVES THE DISPLACEMENT FROM
C     THE START OF THE ASSEMBLED MATRIX(OF ORDER IX) OF THE DIAGONAL
C     ENTRY IN ITS ROW IY.
      IDIAG(IX,IY) = ((IY-1)* (2*IX-IY+2))/2
C     ..
C     .. Executable Statements ..
C INITIALIZATION.
C NBLK IS THE NUMBER OF BLOCK PIVOTS USED.
      NBLK = 0
      NTWO = 0
      NEIG = 0
      NCMPBI = 0
      NCMPBR = 0
      MAXFRT = 0
      NRLBDU = 0
      NIRBDU = 0
C A PRIVATE VARIABLE UU IS SET TO CNTL(1), SO THAT CNTL(1) WILL REMAIN
C UNALTERED.
      UU = MIN(CNTL(1),HALF)
      UU = MAX(UU,-HALF)
      DO 10 I = 1,N
        IW2(I) = 0
   10 CONTINUE
C IWPOS IS POINTER TO FIRST FREE POSITION FOR FACTORS IN IW.
C POSFAC IS POINTER FOR FACTORS IN A. AT EACH PASS THROUGH THE
C     MAJOR LOOP POSFAC INITIALLY POINTS TO THE FIRST FREE LOCATION
C     IN A AND THEN IS SET TO THE POSITION OF THE CURRENT PIVOT IN A.
C ISTK IS POINTER TO TOP OF STACK IN IW.
C ISTK2 IS POINTER TO BOTTOM OF STACK IN IW (NEEDED BY COMPRESS).
C ASTK IS POINTER TO TOP OF STACK IN A.
C ASTK2 IS POINTER TO BOTTOM OF STACK IN A (NEEDED BY COMPRESS).
C IINPUT IS POINTER TO CURRENT POSITION IN ORIGINAL ROWS IN IW.
C AINPUT IS POINTER TO CURRENT POSITION IN ORIGINAL ROWS IN A.
C AZERO IS POINTER TO LAST POSITION ZEROED IN A.
C NTOTPV IS THE TOTAL NUMBER OF PIVOTS SELECTED. THIS IS USED
C     TO DETERMINE WHETHER THE MATRIX IS SINGULAR.
      IWPOS = 2
      POSFAC = 1
      ISTK = LIW - NZ + 1
      ISTK2 = ISTK - 1
      ASTK = LA - NZ + 1
      ASTK2 = ASTK - 1
      IINPUT = ISTK
      AINPUT = ASTK
      AZERO = 0
      NTOTPV = 0
C NUMASS IS THE ACCUMULATED NUMBER OF ROWS ASSEMBLED SO FAR.
      NUMASS = 0
C
C EACH PASS THROUGH THIS MAIN LOOP PERFORMS ALL THE OPERATIONS
C     ASSOCIATED WITH ONE SET OF ASSEMBLY/ELIMINATIONS.
      DO 760 IASS = 1,NSTEPS
C NASS WILL BE SET TO THE NUMBER OF FULLY ASSEMBLED VARIABLES IN
C     CURRENT NEWLY CREATED ELEMENT.
        NASS = NELIM(IASS)
C NEWEL IS A POINTER INTO IW TO CONTROL OUTPUT OF INTEGER INFORMATION
C     FOR NEWLY CREATED ELEMENT.
        NEWEL = IWPOS + 1
C SYMBOLICALLY ASSEMBLE INCOMING ROWS AND GENERATED STACK ELEMENTS
C     ORDERING THE RESULTANT ELEMENT ACCORDING TO PERMUTATION PERM.  WE
C     ASSEMBLE THE STACK ELEMENTS FIRST BECAUSE THESE WILL ALREADY BE
C     ORDERED.
C SET HEADER POINTER FOR MERGE OF INDEX LISTS.
        JFIRST = N + 1
C INITIALIZE NUMBER OF VARIABLES IN CURRENT FRONT.
        NFRONT = 0
        NUMSTK = NSTK(IASS)
        LTOPST = 1
        LNASS = 0
C JUMP IF NO STACK ELEMENTS ARE BEING ASSEMBLED AT THIS STAGE.
        IF (NUMSTK.EQ.0) GO TO 80
        J2 = ISTK - 1
        LNASS = NASS
        LTOPST = ((IW(ISTK)+1)*IW(ISTK))/2
        DO 70 IELL = 1,NUMSTK
C ASSEMBLE ELEMENT IELL PLACING
C     THE INDICES INTO A LINKED LIST IN IW2 ORDERED
C     ACCORDING TO PERM.
          JNEXT = JFIRST
          JLAST = N + 1
          J1 = J2 + 2
          J2 = J1 - 1 + IW(J1-1)
C RUN THROUGH INDEX LIST OF STACK ELEMENT IELL.
          DO 60 JJ = J1,J2
            J = IW(JJ)
C JUMP IF ALREADY ASSEMBLED
            IF (IW2(J).GT.0) GO TO 60
            JNEW = PERM(J)
C IF VARIABLE WAS PREVIOUSLY FULLY SUMMED BUT WAS NOT PIVOTED ON
C     EARLIER BECAUSE OF NUMERICAL TEST, INCREMENT NUMBER OF FULLY
C     SUMMED ROWS/COLUMNS IN FRONT.
            IF (JNEW.LE.NUMASS) NASS = NASS + 1
C FIND POSITION IN LINKED LIST FOR NEW VARIABLE.  NOTE THAT WE START
C     FROM WHERE WE LEFT OFF AFTER ASSEMBLY OF PREVIOUS VARIABLE.
            DO 20 IDUMMY = 1,N
              IF (JNEXT.EQ.N+1) GO TO 30
              IF (PERM(JNEXT).GT.JNEW) GO TO 30
              JLAST = JNEXT
              JNEXT = IW2(JLAST)
   20       CONTINUE
   30       IF (JLAST.NE.N+1) GO TO 40
            JFIRST = J
            GO TO 50

   40       IW2(JLAST) = J
   50       IW2(J) = JNEXT
            JLAST = J
C INCREMENT NUMBER OF VARIABLES IN THE FRONT.
            NFRONT = NFRONT + 1
   60     CONTINUE
   70   CONTINUE
        LNASS = NASS - LNASS
C NOW INCORPORATE ORIGINAL ROWS.  NOTE THAT THE COLUMNS IN THESE ROWS
C     NEED NOT BE IN ORDER. WE ALSO PERFORM
C     A SWOP SO THAT THE DIAGONAL ENTRY IS THE FIRST IN ITS
C     ROW.  THIS ALLOWS US TO AVOID STORING THE INVERSE OF ARRAY PERM.
   80   NUMORG = NELIM(IASS)
        J1 = IINPUT
        DO 150 IORG = 1,NUMORG
          J = -IW(J1)
          DO 140 IDUMMY = 1,LIW
            JNEW = PERM(J)
C JUMP IF VARIABLE ALREADY INCLUDED.
            IF (IW2(J).GT.0) GO TO 130
C HERE WE MUST ALWAYS START OUR SEARCH AT THE BEGINNING.
            JLAST = N + 1
            JNEXT = JFIRST
            DO 90 JDUMMY = 1,N
              IF (JNEXT.EQ.N+1) GO TO 100
              IF (PERM(JNEXT).GT.JNEW) GO TO 100
              JLAST = JNEXT
              JNEXT = IW2(JLAST)
   90       CONTINUE
  100       IF (JLAST.NE.N+1) GO TO 110
            JFIRST = J
            GO TO 120

  110       IW2(JLAST) = J
  120       IW2(J) = JNEXT
C INCREMENT NUMBER OF VARIABLES IN FRONT.
            NFRONT = NFRONT + 1
  130       J1 = J1 + 1
            IF (J1.GT.LIW) GO TO 150
            J = IW(J1)
            IF (J.LT.0) GO TO 150
  140     CONTINUE
  150   CONTINUE
C NOW RUN THROUGH LINKED LIST IW2 PUTTING INDICES OF VARIABLES IN NEW
C     ELEMENT INTO IW AND SETTING IW2 ENTRY TO POINT TO THE RELATIVE
C     POSITION OF THE VARIABLE IN THE NEW ELEMENT.
        IF (NEWEL+NFRONT.LT.ISTK) GO TO 160
C COMPRESS IW.
        CALL MA27PD(A,IW,ISTK,ISTK2,IINPUT,2,NCMPBR,NCMPBI)
        IF (NEWEL+NFRONT.LT.ISTK) GO TO 160
        INFO(2) = LIW + 1 + NEWEL + NFRONT - ISTK
        GO TO 770

  160   J = JFIRST
        DO 170 IFR = 1,NFRONT
          NEWEL = NEWEL + 1
          IW(NEWEL) = J
          JNEXT = IW2(J)
          IW2(J) = NEWEL - (IWPOS+1)
          J = JNEXT
  170   CONTINUE
C
C ASSEMBLE REALS INTO FRONTAL MATRIX.
        MAXFRT = MAX(MAXFRT,NFRONT)
        IW(IWPOS) = NFRONT
C FIRST ZERO OUT FRONTAL MATRIX AS APPROPRIATE FIRST CHECKING TO SEE
C     IF THERE IS SUFFICIENT SPACE.
        LAELL = ((NFRONT+1)*NFRONT)/2
        APOS2 = POSFAC + LAELL - 1
        IF (NUMSTK.NE.0) LNASS = LNASS* (2*NFRONT-LNASS+1)/2
        IF (POSFAC+LNASS-1.GE.ASTK) GO TO 180
        IF (APOS2.LT.ASTK+LTOPST-1) GO TO 190
C COMPRESS A.
  180   CALL MA27PD(A,IW,ASTK,ASTK2,AINPUT,1,NCMPBR,NCMPBI)
        IF (POSFAC+LNASS-1.GE.ASTK) GO TO 780
        IF (APOS2.GE.ASTK+LTOPST-1) GO TO 780
  190   IF (APOS2.LE.AZERO) GO TO 220
        APOS = AZERO + 1
        LAPOS2 = MIN(APOS2,ASTK-1)
        IF (LAPOS2.LT.APOS) GO TO 210
        DO 200 K = APOS,LAPOS2
          A(K) = ZERO
  200   CONTINUE
  210   AZERO = APOS2
C JUMP IF THERE ARE NO STACK ELEMENTS TO ASSEMBLE.
  220   IF (NUMSTK.EQ.0) GO TO 260
C PLACE REALS CORRESPONDING TO STACK ELEMENTS IN CORRECT POSITIONS IN A.
        DO 250 IELL = 1,NUMSTK
          J1 = ISTK + 1
          J2 = ISTK + IW(ISTK)
          DO 240 JJ = J1,J2
            IROW = IW(JJ)
            IROW = IW2(IROW)
            APOS = POSFAC + IDIAG(NFRONT,IROW)
            DO 230 JJJ = JJ,J2
              J = IW(JJJ)
              APOS2 = APOS + IW2(J) - IROW
              A(APOS2) = A(APOS2) + A(ASTK)
              A(ASTK) = ZERO
              ASTK = ASTK + 1
  230       CONTINUE
  240     CONTINUE
C INCREMENT STACK POINTER.
          ISTK = J2 + 1
  250   CONTINUE
C INCORPORATE REALS FROM ORIGINAL ROWS.
  260   DO 280 IORG = 1,NUMORG
          J = -IW(IINPUT)
C WE CAN DO THIS BECAUSE THE DIAGONAL IS NOW THE FIRST ENTRY.
          IROW = IW2(J)
          APOS = POSFAC + IDIAG(NFRONT,IROW)
C THE FOLLOWING LOOP GOES FROM 1 TO NZ BECAUSE THERE MAY BE DUPLICATES.
          DO 270 IDUMMY = 1,NZ
            APOS2 = APOS + IW2(J) - IROW
            A(APOS2) = A(APOS2) + A(AINPUT)
            AINPUT = AINPUT + 1
            IINPUT = IINPUT + 1
            IF (IINPUT.GT.LIW) GO TO 280
            J = IW(IINPUT)
            IF (J.LT.0) GO TO 280
  270     CONTINUE
  280   CONTINUE
C RESET IW2 AND NUMASS.
        NUMASS = NUMASS + NUMORG
        J1 = IWPOS + 2
        J2 = IWPOS + NFRONT + 1
        DO 290 K = J1,J2
          J = IW(K)
          IW2(J) = 0
  290   CONTINUE
C PERFORM PIVOTING ON ASSEMBLED ELEMENT.
C NPIV IS THE NUMBER OF PIVOTS SO FAR SELECTED.
C LNPIV IS THE NUMBER OF PIVOTS SELECTED AFTER THE LAST PASS THROUGH
C     THE THE FOLLOWING LOOP.
        LNPIV = -1
        NPIV = 0
        DO 650 KDUMMY = 1,NASS
          IF (NPIV.EQ.NASS) GO TO 660
          IF (NPIV.EQ.LNPIV) GO TO 660
          LNPIV = NPIV
          NPIVP1 = NPIV + 1
C JPIV IS USED AS A FLAG TO INDICATE WHEN 2 BY 2 PIVOTING HAS OCCURRED
C     SO THAT IPIV IS INCREMENTED CORRECTLY.
          JPIV = 1
C NASS IS MAXIMUM POSSIBLE NUMBER OF PIVOTS.
C WE EITHER TAKE THE DIAGONAL ENTRY OR THE 2 BY 2 PIVOT WITH THE
C     LARGEST OFF-DIAGONAL AT EACH STAGE.
C EACH PASS THROUGH THIS LOOP TRIES TO CHOOSE ONE PIVOT.
          DO 640 IPIV = NPIVP1,NASS
            JPIV = JPIV - 1
C JUMP IF WE HAVE JUST PROCESSED A 2 BY 2 PIVOT.
            IF (JPIV.EQ.1) GO TO 640
            APOS = POSFAC + IDIAG(NFRONT-NPIV,IPIV-NPIV)
C IF THE USER HAS INDICATED THAT THE MATRIX IS DEFINITE, WE
C     DO NOT NEED TO TEST FOR STABILITY BUT WE DO CHECK TO SEE IF THE
C     PIVOT IS NON-ZERO OR HAS CHANGED SIGN.
C     IF IT IS ZERO, WE EXIT WITH AN ERROR. IF IT HAS CHANGED SIGN
C     AND U WAS SET NEGATIVE, THEN WE AGAIN EXIT IMMEDIATELY.  IF THE
C     PIVOT CHANGES SIGN AND U WAS ZERO, WE CONTINUE WITH THE
C     FACTORIZATION BUT PRINT A WARNING MESSAGE ON UNIT ICNTL(2).
C ISNPIV HOLDS A FLAG FOR THE SIGN OF THE PIVOTS TO DATE SO THAT
C     A SIGN CHANGE WHEN DECOMPOSING AN ALLEGEDLY DEFINITE MATRIX CAN
C     BE DETECTED.
            IF (UU.GT.ZERO) GO TO 320
            IF (ABS(A(APOS)).LE.CNTL(3)) GO TO 790
C JUMP IF THIS IS NOT THE FIRST PIVOT TO BE SELECTED.
            IF (NTOTPV.GT.0) GO TO 300
C SET ISNPIV.
            IF (A(APOS).GT.ZERO) ISNPIV = 1
            IF (A(APOS).LT.ZERO) ISNPIV = -1
  300       IF (A(APOS).GT.ZERO .AND. ISNPIV.EQ.1) GO TO 560
            IF (A(APOS).LT.ZERO .AND. ISNPIV.EQ.-1) GO TO 560
            IF (INFO(1).NE.2) INFO(2) = 0
            INFO(2) = INFO(2) + 1
            INFO(1) = 2
            I = NTOTPV + 1
            IF (ICNTL(2).GT.0 .AND. INFO(2).LE.10) THEN
              WRITE (ICNTL(2),FMT=310) INFO(1),I
            END IF

  310       FORMAT (' *** WARNING MESSAGE FROM SUBROUTINE MA27BD',
     +              '  *** INFO(1) =',I2,/,' PIVOT',I6,
     +             ' HAS DIFFERENT SIGN FROM THE PREVIOUS ONE')

            ISNPIV = -ISNPIV
            IF (UU.EQ.ZERO) GO TO 560
            GO TO 800

  320       AMAX = ZERO
            TMAX = AMAX
C FIND LARGEST ENTRY TO RIGHT OF DIAGONAL IN ROW OF PROSPECTIVE PIVOT
C     IN THE FULLY-SUMMED PART.  ALSO RECORD COLUMN OF THIS LARGEST
C     ENTRY.
            J1 = APOS + 1
            J2 = APOS + NASS - IPIV
            IF (J2.LT.J1) GO TO 340
            DO 330 JJ = J1,J2
              IF (ABS(A(JJ)).LE.AMAX) GO TO 330
              JMAX = IPIV + JJ - J1 + 1
              AMAX = ABS(A(JJ))
  330       CONTINUE
C DO SAME AS ABOVE FOR NON-FULLY-SUMMED PART ONLY HERE WE DO NOT NEED
C     TO RECORD COLUMN SO LOOP IS SIMPLER.
  340       J1 = J2 + 1
            J2 = APOS + NFRONT - IPIV
            IF (J2.LT.J1) GO TO 360
            DO 350 JJ = J1,J2
              TMAX = MAX(ABS(A(JJ)),TMAX)
  350       CONTINUE
C NOW CALCULATE LARGEST ENTRY IN OTHER PART OF ROW.
  360       RMAX = MAX(TMAX,AMAX)
            APOS1 = APOS
            KK = NFRONT - IPIV
            LT = IPIV - (NPIV+1)
            IF (LT.EQ.0) GO TO 380
            DO 370 K = 1,LT
              KK = KK + 1
              APOS1 = APOS1 - KK
              RMAX = MAX(RMAX,ABS(A(APOS1)))
  370       CONTINUE
C JUMP IF STABILITY TEST SATISFIED.
  380       IF (ABS(A(APOS)).GT.MAX(CNTL(3),UU*RMAX)) GO TO 450
C CHECK BLOCK PIVOT OF ORDER 2 FOR STABILITY.
            IF (ABS(AMAX).LE.CNTL(3)) GO TO 640
            APOS2 = POSFAC + IDIAG(NFRONT-NPIV,JMAX-NPIV)
            DETPIV = A(APOS)*A(APOS2) - AMAX*AMAX
            THRESH = ABS(DETPIV)
C SET THRESH TO U TIMES THE RECIPROCAL OF THE MAX-NORM OF THE INVERSE
C     OF THE PROSPECTIVE BLOCK.
            THRESH = THRESH/ (UU*MAX(ABS(A(APOS))+AMAX,
     +               ABS(A(APOS2))+AMAX))
C CHECK 2 BY 2 PIVOT FOR STABILITY.
C FIRST CHECK AGAINST ROW IPIV.
            IF (THRESH.LE.RMAX) GO TO 640
C FIND LARGEST ENTRY IN ROW JMAX.
C FIND MAXIMUM TO THE RIGHT OF THE DIAGONAL.
            RMAX = ZERO
            J1 = APOS2 + 1
            J2 = APOS2 + NFRONT - JMAX
            IF (J2.LT.J1) GO TO 400
            DO 390 JJ = J1,J2
              RMAX = MAX(RMAX,ABS(A(JJ)))
  390       CONTINUE
C NOW CHECK TO THE LEFT OF THE DIAGONAL.
C WE USE TWO LOOPS TO AVOID TESTING FOR ROW IPIV INSIDE THE LOOP.
  400       KK = NFRONT - JMAX + 1
            APOS3 = APOS2
            JMXMIP = JMAX - IPIV - 1
            IF (JMXMIP.EQ.0) GO TO 420
            DO 410 K = 1,JMXMIP
              APOS2 = APOS2 - KK
              KK = KK + 1
              RMAX = MAX(RMAX,ABS(A(APOS2)))
  410       CONTINUE
  420       IPMNP = IPIV - NPIV - 1
            IF (IPMNP.EQ.0) GO TO 440
            APOS2 = APOS2 - KK
            KK = KK + 1
            DO 430 K = 1,IPMNP
              APOS2 = APOS2 - KK
              KK = KK + 1
              RMAX = MAX(RMAX,ABS(A(APOS2)))
  430       CONTINUE
  440       IF (THRESH.LE.RMAX) GO TO 640
            PIVSIZ = 2
            GO TO 460

  450       PIVSIZ = 1
  460       IROW = IPIV - NPIV
C
C PIVOT HAS BEEN CHOSEN.  IF BLOCK PIVOT OF ORDER 2, PIVSIZ IS EQUAL TO
C     TWO OTHERWISE PIVSIZ EQUALS ONE..
C THE FOLLOWING LOOP MOVES THE PIVOT BLOCK TO THE TOP LEFT HAND CORNER
C     OF THE FRONTAL MATRIX.
            DO 550 KROW = 1,PIVSIZ
C WE JUMP IF SWOP IS NOT NECESSARY.
              IF (IROW.EQ.1) GO TO 530
              J1 = POSFAC + IROW
              J2 = POSFAC + NFRONT - (NPIV+1)
              IF (J2.LT.J1) GO TO 480
              APOS2 = APOS + 1
C SWOP PORTION OF ROWS WHOSE COLUMN INDICES ARE GREATER THAN LATER ROW.
              DO 470 JJ = J1,J2
                SWOP = A(APOS2)
                A(APOS2) = A(JJ)
                A(JJ) = SWOP
                APOS2 = APOS2 + 1
  470         CONTINUE
  480         J1 = POSFAC + 1
              J2 = POSFAC + IROW - 2
              APOS2 = APOS
              KK = NFRONT - (IROW+NPIV)
              IF (J2.LT.J1) GO TO 500
C SWOP PORTION OF ROWS/COLUMNS WHOSE INDICES LIE BETWEEN THE TWO ROWS.
              DO 490 JJJ = J1,J2
                JJ = J2 - JJJ + J1
                KK = KK + 1
                APOS2 = APOS2 - KK
                SWOP = A(APOS2)
                A(APOS2) = A(JJ)
                A(JJ) = SWOP
  490         CONTINUE
  500         IF (NPIV.EQ.0) GO TO 520
              APOS1 = POSFAC
              KK = KK + 1
              APOS2 = APOS2 - KK
C SWOP PORTION OF COLUMNS WHOSE INDICES ARE LESS THAN EARLIER ROW.
              DO 510 JJ = 1,NPIV
                KK = KK + 1
                APOS1 = APOS1 - KK
                APOS2 = APOS2 - KK
                SWOP = A(APOS2)
                A(APOS2) = A(APOS1)
                A(APOS1) = SWOP
  510         CONTINUE
C SWOP DIAGONALS AND INTEGER INDEXING INFORMATION
  520         SWOP = A(APOS)
              A(APOS) = A(POSFAC)
              A(POSFAC) = SWOP
              IPOS = IWPOS + NPIV + 2
              IEXCH = IWPOS + IROW + NPIV + 1
              ISWOP = IW(IPOS)
              IW(IPOS) = IW(IEXCH)
              IW(IEXCH) = ISWOP
  530         IF (PIVSIZ.EQ.1) GO TO 550
C SET VARIABLES FOR THE SWOP OF SECOND ROW OF BLOCK PIVOT.
              IF (KROW.EQ.2) GO TO 540
              IROW = JMAX - (NPIV+1)
              JPOS = POSFAC
              POSFAC = POSFAC + NFRONT - NPIV
              NPIV = NPIV + 1
              APOS = APOS3
              GO TO 550
C RESET VARIABLES PREVIOUSLY SET FOR SECOND PASS.
  540         NPIV = NPIV - 1
              POSFAC = JPOS
  550       CONTINUE
C
            IF (PIVSIZ.EQ.2) GO TO 600
C PERFORM THE ELIMINATION USING ENTRY (IPIV,IPIV) AS PIVOT.
C WE STORE U AND DINVERSE.
  560       A(POSFAC) = ONE/A(POSFAC)
            IF (A(POSFAC).LT.ZERO) NEIG = NEIG + 1
            J1 = POSFAC + 1
            J2 = POSFAC + NFRONT - (NPIV+1)
            IF (J2.LT.J1) GO TO 590
            IBEG = J2 + 1
            DO 580 JJ = J1,J2
              AMULT = -A(JJ)*A(POSFAC)
              IEND = IBEG + NFRONT - (NPIV+JJ-J1+2)
C THE FOLLOWING SPECIAL COMMENT FORCES VECTORIZATION ON THE CRAY-1.
CDIR$ IVDEP
              DO 570 IROW = IBEG,IEND
                JCOL = JJ + IROW - IBEG
                A(IROW) = A(IROW) + AMULT*A(JCOL)
  570         CONTINUE
              IBEG = IEND + 1
              A(JJ) = AMULT
  580       CONTINUE
  590       NPIV = NPIV + 1
            NTOTPV = NTOTPV + 1
            JPIV = 1
            POSFAC = POSFAC + NFRONT - NPIV + 1
            GO TO 640
C PERFORM ELIMINATION USING BLOCK PIVOT OF ORDER TWO.
C REPLACE BLOCK PIVOT BY ITS INVERSE.
C SET FLAG TO INDICATE USE OF 2 BY 2 PIVOT IN IW.
  600       IPOS = IWPOS + NPIV + 2
            NTWO = NTWO + 1
            IW(IPOS) = -IW(IPOS)
            POSPV1 = POSFAC
            POSPV2 = POSFAC + NFRONT - NPIV
            SWOP = A(POSPV2)
            IF (DETPIV.LT.ZERO) NEIG = NEIG + 1
            IF (DETPIV.GT.ZERO .AND. SWOP.LT.ZERO) NEIG = NEIG + 2
            A(POSPV2) = A(POSPV1)/DETPIV
            A(POSPV1) = SWOP/DETPIV
            A(POSPV1+1) = -A(POSPV1+1)/DETPIV
            J1 = POSPV1 + 2
            J2 = POSPV1 + NFRONT - (NPIV+1)
            IF (J2.LT.J1) GO TO 630
            JJ1 = POSPV2
            IBEG = POSPV2 + NFRONT - (NPIV+1)
            DO 620 JJ = J1,J2
              JJ1 = JJ1 + 1
              AMULT1 = - (A(POSPV1)*A(JJ)+A(POSPV1+1)*A(JJ1))
              AMULT2 = - (A(POSPV1+1)*A(JJ)+A(POSPV2)*A(JJ1))
              IEND = IBEG + NFRONT - (NPIV+JJ-J1+3)
C THE FOLLOWING SPECIAL COMMENT FORCES VECTORIZATION ON THE CRAY-1.
CDIR$ IVDEP
              DO 610 IROW = IBEG,IEND
                K1 = JJ + IROW - IBEG
                K2 = JJ1 + IROW - IBEG
                A(IROW) = A(IROW) + AMULT1*A(K1) + AMULT2*A(K2)
  610         CONTINUE
              IBEG = IEND + 1
              A(JJ) = AMULT1
              A(JJ1) = AMULT2
  620       CONTINUE
  630       NPIV = NPIV + 2
            NTOTPV = NTOTPV + 2
            JPIV = 2
            POSFAC = POSPV2 + NFRONT - NPIV + 1
  640     CONTINUE
  650   CONTINUE
C END OF MAIN ELIMINATION LOOP.
C
  660   IF (NPIV.NE.0) NBLK = NBLK + 1
        IOLDPS = IWPOS
        IWPOS = IWPOS + NFRONT + 2
        IF (NPIV.EQ.0) GO TO 690
        IF (NPIV.GT.1) GO TO 680
        IW(IOLDPS) = -IW(IOLDPS)
        DO 670 K = 1,NFRONT
          J1 = IOLDPS + K
          IW(J1) = IW(J1+1)
  670   CONTINUE
        IWPOS = IWPOS - 1
        GO TO 690

  680   IW(IOLDPS+1) = NPIV
C COPY REMAINDER OF ELEMENT TO TOP OF STACK
  690   LIELL = NFRONT - NPIV
        IF (LIELL.EQ.0 .OR. IASS.EQ.NSTEPS) GO TO 750
        IF (IWPOS+LIELL.LT.ISTK) GO TO 700
        CALL MA27PD(A,IW,ISTK,ISTK2,IINPUT,2,NCMPBR,NCMPBI)
  700   ISTK = ISTK - LIELL - 1
        IW(ISTK) = LIELL
        J1 = ISTK
        KK = IWPOS - LIELL - 1
C THE FOLLOWING SPECIAL COMMENT FORCES VECTORIZATION ON THE CRAY-1.
CDIR$ IVDEP
        DO 710 K = 1,LIELL
          J1 = J1 + 1
          KK = KK + 1
          IW(J1) = IW(KK)
  710   CONTINUE
C WE COPY IN REVERSE DIRECTION TO AVOID OVERWRITE PROBLEMS.
        LAELL = ((LIELL+1)*LIELL)/2
        KK = POSFAC + LAELL
        IF (KK.NE.ASTK) GO TO 720
        ASTK = ASTK - LAELL
        GO TO 740
C THE MOVE AND ZEROING OF ARRAY A IS PERFORMED WITH TWO LOOPS SO
C THAT THE CRAY-1 WILL VECTORIZE THEM SAFELY.
  720   KMAX = KK - 1
C THE FOLLOWING SPECIAL COMMENT FORCES VECTORIZATION ON THE CRAY-1.
CDIR$ IVDEP
        DO 730 K = 1,LAELL
          KK = KK - 1
          ASTK = ASTK - 1
          A(ASTK) = A(KK)
  730   CONTINUE
        KMAX = MIN(KMAX,ASTK-1)
        DO 735 K = KK,KMAX
          A(K) = ZERO
  735   CONTINUE
  740   AZERO = MIN(AZERO,ASTK-1)
  750   IF (NPIV.EQ.0) IWPOS = IOLDPS
  760 CONTINUE
C
C END OF LOOP ON TREE NODES.
C
      IW(1) = NBLK
      IF (NTWO.GT.0) IW(1) = -NBLK
      NRLBDU = POSFAC - 1
      NIRBDU = IWPOS - 1
      IF (NTOTPV.EQ.N) GO TO 810
      INFO(1) = 3
      INFO(2) = NTOTPV
      GO TO 810
C **** ERROR RETURNS ****
  770 INFO(1) = -3
      GO TO 810

  780 INFO(1) = -4
      INFO(2) = LA + MAX(POSFAC+LNASS,APOS2-LTOPST+2) - ASTK
      GO TO 810

  790 INFO(1) = -5
      INFO(2) = NTOTPV + 1
      GO TO 810

  800 INFO(1) = -6
      INFO(2) = NTOTPV + 1
  810 CONTINUE
      INFO(9) = NRLBDU
      INFO(10) = NIRBDU
      INFO(12) = NCMPBR
      INFO(13) = NCMPBI
      INFO(14) = NTWO
      INFO(15) = NEIG

      RETURN
      END
      SUBROUTINE MA27PD(A,IW,J1,J2,ITOP,IREAL,NCMPBR,NCMPBI)
C THIS SUBROUTINE PERFORMS A VERY SIMPLE COMPRESS (BLOCK MOVE).
C     ENTRIES J1 TO J2 (INCL.) IN A OR IW AS APPROPRIATE ARE MOVED TO
C     OCCUPY THE POSITIONS IMMEDIATELY PRIOR TO POSITION ITOP.
C A/IW HOLD THE ARRAY BEING COMPRESSED.
C J1/J2 DEFINE THE ENTRIES BEING MOVED.
C ITOP DEFINES THE POSITION IMMEDIATELY AFTER THE POSITIONS TO WHICH
C     J1 TO J2 ARE MOVED.
C IREAL MUST BE SET BY THE USER TO 2 IF THE MOVE IS ON ARRAY IW,
C     ANY OTHER VALUE WILL PERFORM THE MOVE ON A.
C NCMPBR and NCMPBI, see INFO(12) and INFO(13) in MA27A/AD (ACCUMULATE
C     THE NUMBER OF COMPRESSES OF THE REALS AND INTEGERS PERFORMED BY
C     MA27B/BD.
C
C     .. Scalar Arguments ..
      INTEGER IREAL,ITOP,J1,J2,NCMPBR,NCMPBI
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*)
      INTEGER IW(*)
C     ..
C     .. Local Scalars ..
      INTEGER IPOS,JJ,JJJ
C     ..
C     .. Executable Statements ..
      IPOS = ITOP - 1
      IF (J2.EQ.IPOS) GO TO 50
      IF (IREAL.EQ.2) GO TO 20
      NCMPBR = NCMPBR + 1
      IF (J1.GT.J2) GO TO 40
      DO 10 JJJ = J1,J2
        JJ = J2 - JJJ + J1
        A(IPOS) = A(JJ)
        IPOS = IPOS - 1
   10 CONTINUE
      GO TO 40

   20 NCMPBI = NCMPBI + 1
      IF (J1.GT.J2) GO TO 40
      DO 30 JJJ = J1,J2
        JJ = J2 - JJJ + J1
        IW(IPOS) = IW(JJ)
        IPOS = IPOS - 1
   30 CONTINUE
   40 J2 = ITOP - 1
      J1 = IPOS + 1
   50 RETURN

      END
      SUBROUTINE MA27QD(N,A,LA,IW,LIW,W,MAXFNT,RHS,IW2,NBLK,LATOP,ICNTL)
C THIS SUBROUTINE PERFORMS FORWARD ELIMINATION
C     USING THE FACTOR U TRANSPOSE STORED IN A/IA AFTER MA27B/BD.
C
C N      - MUST BE SET TO THE ORDER OF THE MATRIX. NOT ALTERED
C          BY MA27Q/QD.
C A      - MUST BE SET TO HOLD THE REAL VALUES
C          CORRESPONDING TO THE FACTORS OF DINVERSE AND U.  THIS MUST BE
C          UNCHANGED SINCE THE PRECEDING CALL TO MA27B/BD.  NOT ALTERED
C          BY MA27Q/QD.
C LA     - LENGTH OF ARRAY A.  NOT ALTERED BY MA27Q/QD.
C IW     - HOLDS THE INTEGER INDEXING
C          INFORMATION FOR THE MATRIX FACTORS IN A.  THIS MUST BE
C          UNCHANGED SINCE THE PRECEDING CALL TO MA27B/BD.  NOT ALTERED
C          BY MA27Q/QD.
C LIW    - LENGTH OF ARRAY IW.  NOT ALTERED BY MA27Q/QD.
C W      - USED
C          AS WORKSPACE BY MA27Q/QD TO HOLD THE COMPONENTS OF THE RIGHT
C          HAND SIDES CORRESPONDING TO CURRENT BLOCK PIVOTAL ROWS.
C MAXFNT - MUST BE SET TO THE LARGEST NUMBER OF
C          VARIABLES IN ANY BLOCK PIVOT ROW.  THIS VALUE WILL HAVE
C          BEEN OUTPUT BY MA27B/BD.  NOT ALTERED BY MA27Q/QD.
C RHS    - ON INPUT,
C          MUST BE SET TO HOLD THE RIGHT HAND SIDES FOR THE EQUATIONS
C          WHICH THE USER DESIRES TO SOLVE.  ON OUTPUT, RHS WILL HOLD
C          THE MODIFIED VECTORS CORRESPONDING TO PERFORMING FORWARD
C          ELIMINATION ON THE RIGHT HAND SIDES.
C IW2    - NEED NOT BE SET ON ENTRY. ON EXIT IW2(I) (I = 1,NBLK)
C          WILL HOLD POINTERS TO THE
C          BEGINNING OF EACH BLOCK PIVOT IN ARRAY IW.
C NBLK   - NUMBER OF BLOCK PIVOT ROWS. NOT ALTERED BY MA27Q/QD.
C LATOP  - NEED NOT BE SET ON ENTRY. ON EXIT, IT IS THE POSITION IN
C          A OF THE LAST ENTRY IN THE  FACTORS. IT MUST BE PASSED
C          UNCHANGED TO MA27R/RD.
C ICNTL is an INTEGER array of length 30, see MA27A/AD.
C   ICNTL(IFRLVL+I) I=1,20 IS USED TO CONTROL WHETHER DIRECT OR INDIRECT
C     ACCESS IS USED BY MA27C/CD.  INDIRECT ACCESS IS EMPLOYED
C     IN FORWARD AND BACK SUBSTITUTION RESPECTIVELY IF THE SIZE OF
C     A BLOCK IS LESS THAN ICNTL(IFRLVL+MIN(10,NPIV)) AND
C     ICNTL(IFRLVL+10+MIN(10,NPIV)) RESPECTIVELY, WHERE NPIV IS THE
C     NUMBER OF PIVOTS IN THE BLOCK.
C
      INTEGER IFRLVL
      PARAMETER ( IFRLVL=5 )
C     .. Scalar Arguments ..
      INTEGER LA,LATOP,LIW,MAXFNT,N,NBLK
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA),RHS(N),W(MAXFNT)
      INTEGER IW(LIW),IW2(NBLK),ICNTL(30)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION W1,W2
      INTEGER APOS,IBLK,IFR,ILVL,IPIV,IPOS,IRHS,IROW,IST,J,J1,J2,J3,JJ,
     +        JPIV,K,K1,K2,K3,LIELL,NPIV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN
C     ..
C     .. Executable Statements ..
C APOS. RUNNING POINTER TO CURRENT PIVOT POSITION IN ARRAY A.
C IPOS. RUNNING POINTER TO BEGINNING OF BLOCK PIVOT ROW IN IW.
      APOS = 1
      IPOS = 1
      J2 = 0
      IBLK = 0
      NPIV = 0
      DO 140 IROW = 1,N
        IF (NPIV.GT.0) GO TO 90
        IBLK = IBLK + 1
        IF (IBLK.GT.NBLK) GO TO 150
        IPOS = J2 + 1
C SET UP POINTER IN PREPARATION FOR BACK SUBSTITUTION.
        IW2(IBLK) = IPOS
C ABS(LIELL) IS NUMBER OF VARIABLES (COLUMNS) IN BLOCK PIVOT ROW.
        LIELL = -IW(IPOS)
C NPIV IS NUMBER OF PIVOTS (ROWS) IN BLOCK PIVOT.
        NPIV = 1
        IF (LIELL.GT.0) GO TO 10
        LIELL = -LIELL
        IPOS = IPOS + 1
        NPIV = IW(IPOS)
   10   J1 = IPOS + 1
        J2 = IPOS + LIELL
        ILVL = MIN(NPIV,10)
        IF (LIELL.LT.ICNTL(IFRLVL+ILVL)) GO TO 90
C
C PERFORM OPERATIONS USING DIRECT ADDRESSING.
C
C LOAD APPROPRIATE COMPONENTS OF RIGHT HAND SIDES INTO ARRAY W.
        IFR = 0
        DO 20 JJ = J1,J2
          J = ABS(IW(JJ)+0)
          IFR = IFR + 1
          W(IFR) = RHS(J)
   20   CONTINUE
C JPIV IS USED AS A FLAG SO THAT IPIV IS INCREMENTED CORRECTLY AFTER
C THE USE OF A 2 BY 2 PIVOT.
        JPIV = 1
        J3 = J1
C PERFORM OPERATIONS.
        DO 70 IPIV = 1,NPIV
          JPIV = JPIV - 1
          IF (JPIV.EQ.1) GO TO 70
C JUMP IF WE HAVE A 2 BY 2 PIVOT.
          IF (IW(J3).LT.0) GO TO 40
C PERFORM FORWARD SUBSTITUTION USING 1 BY 1 PIVOT.
          JPIV = 1
          J3 = J3 + 1
          APOS = APOS + 1
          IST = IPIV + 1
          IF (LIELL.LT.IST) GO TO 70
          W1 = W(IPIV)
          K = APOS
          DO 30 J = IST,LIELL
            W(J) = W(J) + A(K)*W1
            K = K + 1
   30     CONTINUE
          APOS = APOS + LIELL - IST + 1
          GO TO 70
C PERFORM OPERATIONS WITH 2 BY 2 PIVOT.
   40     JPIV = 2
          J3 = J3 + 2
          APOS = APOS + 2
          IST = IPIV + 2
          IF (LIELL.LT.IST) GO TO 60
          W1 = W(IPIV)
          W2 = W(IPIV+1)
          K1 = APOS
          K2 = APOS + LIELL - IPIV
          DO 50 J = IST,LIELL
            W(J) = W(J) + W1*A(K1) + W2*A(K2)
            K1 = K1 + 1
            K2 = K2 + 1
   50     CONTINUE
   60     APOS = APOS + 2* (LIELL-IST+1) + 1
   70   CONTINUE
C RELOAD W BACK INTO RHS.
        IFR = 0
        DO 80 JJ = J1,J2
          J = ABS(IW(JJ)+0)
          IFR = IFR + 1
          RHS(J) = W(IFR)
   80   CONTINUE
        NPIV = 0
        GO TO 140
C
C PERFORM OPERATIONS USING INDIRECT ADDRESSING.
C
C JUMP IF WE HAVE A 2 BY 2 PIVOT.
   90   IF (IW(J1).LT.0) GO TO 110
C PERFORM FORWARD SUBSTITUTION USING 1 BY 1 PIVOT.
        NPIV = NPIV - 1
        APOS = APOS + 1
        J1 = J1 + 1
        IF (J1.GT.J2) GO TO 140
        IRHS = IW(J1-1)
        W1 = RHS(IRHS)
        K = APOS
        DO 100 J = J1,J2
          IRHS = ABS(IW(J)+0)
          RHS(IRHS) = RHS(IRHS) + A(K)*W1
          K = K + 1
  100   CONTINUE
        APOS = APOS + J2 - J1 + 1
        GO TO 140
C PERFORM OPERATIONS WITH 2 BY 2 PIVOT
  110   NPIV = NPIV - 2
        J1 = J1 + 2
        APOS = APOS + 2
        IF (J1.GT.J2) GO TO 130
        IRHS = -IW(J1-2)
        W1 = RHS(IRHS)
        IRHS = IW(J1-1)
        W2 = RHS(IRHS)
        K1 = APOS
        K3 = APOS + J2 - J1 + 2
        DO 120 J = J1,J2
          IRHS = ABS(IW(J)+0)
          RHS(IRHS) = RHS(IRHS) + W1*A(K1) + W2*A(K3)
          K1 = K1 + 1
          K3 = K3 + 1
  120   CONTINUE
  130   APOS = APOS + 2* (J2-J1+1) + 1
  140 CONTINUE
  150 LATOP = APOS - 1
      RETURN

      END
      SUBROUTINE MA27RD(N,A,LA,IW,LIW,W,MAXFNT,RHS,IW2,NBLK,LATOP,ICNTL)
C THIS SUBROUTINE PERFORMS BACKWARD ELIMINATION OPERATIONS
C     USING THE FACTORS DINVERSE AND U
C     STORED IN A/IW AFTER MA27B/BD.
C
C N      - MUST BE SET TO THE ORDER OF THE MATRIX. NOT ALTERED
C          BY MA27R/RD.
C A      - MUST BE SET TO HOLD THE REAL VALUES CORRESPONDING
C          TO THE FACTORS OF DINVERSE AND U.  THIS MUST BE
C          UNCHANGED SINCE THE PRECEDING CALL TO MA27B/BD.  NOT ALTERED
C          BY MA27R/RD.
C LA     - LENGTH OF ARRAY A. NOT ALTERED BY MA27R/RD.
C IW     - HOLDS THE INTEGER INDEXING
C          INFORMATION FOR THE MATRIX FACTORS IN A.  THIS MUST BE
C          UNCHANGED SINCE THE PRECEDING CALL TO MA27B/BD.  NOT ALTERED
C          BY MA27R/RD.
C LIW    - LENGTH OF ARRAY IW.  NOT ALTERED BY MA27R/RD.
C W      - USED
C          AS WORKSPACE BY MA27R/RD TO HOLD THE COMPONENTS OF THE RIGHT
C          HAND SIDES CORRESPONDING TO CURRENT BLOCK PIVOTAL ROWS.
C MAXFNT - INTEGER VARIABLE.  MUST BE SET TO THE LARGEST NUMBER OF
C          VARIABLES IN ANY BLOCK PIVOT ROW.  THIS VALUE WAS GIVEN
C          ON OUTPUT FROM MA27B/BD.  NOT ALTERED BY MA27R/RD.
C RHS    - ON INPUT,
C          MUST BE SET TO HOLD THE RIGHT HAND SIDE MODIFIED BY THE
C          FORWARD SUBSTITUTION OPERATIONS.  ON OUTPUT, RHS WILL HOLD
C          THE SOLUTION VECTOR.
C IW2    - ON ENTRY IW2(I) (I = 1,NBLK)
C          MUST HOLD POINTERS TO THE
C          BEGINNING OF EACH BLOCK PIVOT IN ARRAY IW, AS SET BY
C          MA27Q/QD.
C NBLK   - NUMBER OF BLOCK PIVOT ROWS. NOT ALTERED BY MA27R/RD.
C LATOP  - IT IS THE POSITION IN
C          A OF THE LAST ENTRY IN THE  FACTORS. IT MUST BE UNCHANGED
C          SINCE THE CALL TO MA27Q/QD.  IT IS NOT ALTERED BY MA27R/RD.
C ICNTL is an INTEGER array of length 30, see MA27A/AD.
C   ICNTL(IFRLVL+I) I=1,20 IS USED TO CONTROL WHETHER DIRECT OR INDIRECT
C     ACCESS IS USED BY MA27C/CD.  INDIRECT ACCESS IS EMPLOYED
C     IN FORWARD AND BACK SUBSTITUTION RESPECTIVELY IF THE SIZE OF
C     A BLOCK IS LESS THAN ICNTL(IFRLVL+MIN(10,NPIV)) AND
C     ICNTL(IFRLVL+10+MIN(10,NPIV)) RESPECTIVELY, WHERE NPIV IS THE
C     NUMBER OF PIVOTS IN THE BLOCK.
C
      INTEGER IFRLVL
      PARAMETER ( IFRLVL=5 )
C
C     .. Scalar Arguments ..
      INTEGER LA,LATOP,LIW,MAXFNT,N,NBLK
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(LA),RHS(N),W(MAXFNT)
      INTEGER IW(LIW),IW2(NBLK),ICNTL(30)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION W1,W2
      INTEGER APOS,APOS2,I1RHS,I2RHS,IBLK,IFR,IIPIV,IIRHS,ILVL,IPIV,
     +        IPOS,IRHS,IST,J,J1,J2,JJ,JJ1,JJ2,JPIV,JPOS,K,LIELL,LOOP,
     +        NPIV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN
C     ..
C     .. Executable Statements ..
C APOS. RUNNING POINTER TO CURRENT PIVOT POSITION IN ARRAY A.
C IPOS. RUNNING POINTER TO BEGINNING OF CURRENT BLOCK PIVOT ROW.
      APOS = LATOP + 1
      NPIV = 0
      IBLK = NBLK + 1
C RUN THROUGH BLOCK PIVOT ROWS IN THE REVERSE ORDER.
      DO 180 LOOP = 1,N
        IF (NPIV.GT.0) GO TO 110
        IBLK = IBLK - 1
        IF (IBLK.LT.1) GO TO 190
        IPOS = IW2(IBLK)
C ABS(LIELL) IS NUMBER OF VARIABLES (COLUMNS) IN BLOCK PIVOT ROW.
        LIELL = -IW(IPOS)
C NPIV IS NUMBER OF PIVOTS (ROWS) IN BLOCK PIVOT.
        NPIV = 1
        IF (LIELL.GT.0) GO TO 10
        LIELL = -LIELL
        IPOS = IPOS + 1
        NPIV = IW(IPOS)
   10   JPOS = IPOS + NPIV
        J2 = IPOS + LIELL
        ILVL = MIN(10,NPIV) + 10
        IF (LIELL.LT.ICNTL(IFRLVL+ILVL)) GO TO 110
C
C PERFORM OPERATIONS USING DIRECT ADDRESSING.
C
        J1 = IPOS + 1
C LOAD APPROPRIATE COMPONENTS OF RHS INTO W.
        IFR = 0
        DO 20 JJ = J1,J2
          J = ABS(IW(JJ)+0)
          IFR = IFR + 1
          W(IFR) = RHS(J)
   20   CONTINUE
C JPIV IS USED AS A FLAG SO THAT IPIV IS INCREMENTED CORRECTLY AFTER
C     THE USE OF A 2 BY 2 PIVOT.
        JPIV = 1
C PERFORM ELIMINATIONS.
        DO 90 IIPIV = 1,NPIV
          JPIV = JPIV - 1
          IF (JPIV.EQ.1) GO TO 90
          IPIV = NPIV - IIPIV + 1
          IF (IPIV.EQ.1) GO TO 30
C JUMP IF WE HAVE A 2 BY 2 PIVOT.
          IF (IW(JPOS-1).LT.0) GO TO 60
C PERFORM BACK-SUBSTITUTION USING 1 BY 1 PIVOT.
   30     JPIV = 1
          APOS = APOS - (LIELL+1-IPIV)
          IST = IPIV + 1
          W1 = W(IPIV)*A(APOS)
          IF (LIELL.LT.IST) GO TO 50
          JJ1 = APOS + 1
          DO 40 J = IST,LIELL
            W1 = W1 + A(JJ1)*W(J)
            JJ1 = JJ1 + 1
   40     CONTINUE
   50     W(IPIV) = W1
          JPOS = JPOS - 1
          GO TO 90
C PERFORM BACK-SUBSTITUTION OPERATIONS WITH 2 BY 2 PIVOT
   60     JPIV = 2
          APOS2 = APOS - (LIELL+1-IPIV)
          APOS = APOS2 - (LIELL+2-IPIV)
          IST = IPIV + 1
          W1 = W(IPIV-1)*A(APOS) + W(IPIV)*A(APOS+1)
          W2 = W(IPIV-1)*A(APOS+1) + W(IPIV)*A(APOS2)
          IF (LIELL.LT.IST) GO TO 80
          JJ1 = APOS + 2
          JJ2 = APOS2 + 1
          DO 70 J = IST,LIELL
            W1 = W1 + W(J)*A(JJ1)
            W2 = W2 + W(J)*A(JJ2)
            JJ1 = JJ1 + 1
            JJ2 = JJ2 + 1
   70     CONTINUE
   80     W(IPIV-1) = W1
          W(IPIV) = W2
          JPOS = JPOS - 2
   90   CONTINUE
C RELOAD WORKING VECTOR INTO SOLUTION VECTOR.
        IFR = 0
        DO 100 JJ = J1,J2
          J = ABS(IW(JJ)+0)
          IFR = IFR + 1
          RHS(J) = W(IFR)
  100   CONTINUE
        NPIV = 0
        GO TO 180
C
C PERFORM OPERATIONS USING INDIRECT ADDRESSING.
C
  110   IF (NPIV.EQ.1) GO TO 120
C JUMP IF WE HAVE A 2 BY 2 PIVOT.
        IF (IW(JPOS-1).LT.0) GO TO 150
C PERFORM BACK-SUBSTITUTION USING 1 BY 1 PIVOT.
  120   NPIV = NPIV - 1
        APOS = APOS - (J2-JPOS+1)
        IIRHS = IW(JPOS)
        W1 = RHS(IIRHS)*A(APOS)
        J1 = JPOS + 1
        IF (J1.GT.J2) GO TO 140
        K = APOS + 1
        DO 130 J = J1,J2
          IRHS = ABS(IW(J)+0)
          W1 = W1 + A(K)*RHS(IRHS)
          K = K + 1
  130   CONTINUE
  140   RHS(IIRHS) = W1
        JPOS = JPOS - 1
        GO TO 180
C PERFORM OPERATIONS WITH 2 BY 2 PIVOT
  150   NPIV = NPIV - 2
        APOS2 = APOS - (J2-JPOS+1)
        APOS = APOS2 - (J2-JPOS+2)
        I1RHS = -IW(JPOS-1)
        I2RHS = IW(JPOS)
        W1 = RHS(I1RHS)*A(APOS) + RHS(I2RHS)*A(APOS+1)
        W2 = RHS(I1RHS)*A(APOS+1) + RHS(I2RHS)*A(APOS2)
        J1 = JPOS + 1
        IF (J1.GT.J2) GO TO 170
        JJ1 = APOS + 2
        JJ2 = APOS2 + 1
        DO 160 J = J1,J2
          IRHS = ABS(IW(J)+0)
          W1 = W1 + RHS(IRHS)*A(JJ1)
          W2 = W2 + RHS(IRHS)*A(JJ2)
          JJ1 = JJ1 + 1
          JJ2 = JJ2 + 1
  160   CONTINUE
  170   RHS(I1RHS) = W1
        RHS(I2RHS) = W2
        JPOS = JPOS - 2
  180 CONTINUE
  190 RETURN

      END
C *******************************************************************
C COPYRIGHT (c) 1999 Council for the Central Laboratory
*                    of the Research Councils
C All rights reserved.
C
C None of the comments in this Copyright notice between the lines
C of asterisks shall be removed or altered in any way.
C
C This Package is intended for compilation without modification,
C so most of the embedded comments have been removed.
C
C ALL USE IS SUBJECT TO LICENCE. For full details of a UK ACADEMIC
C Licence, see http://hsl.rl.ac.uk/acuk/cou.html
C
C Please note that for a UK ACADEMIC Licence:
C
C 1. The Packages may only be used for academic research or teaching
C    purposes by the Licensee, and must not be copied by the Licensee for
C    use by any other persons. Use of the Packages in any commercial
C    application shall be subject to prior written agreement between
C    Hyprotech UK Limited and the Licensee on suitable terms and
C    conditions, which will include financial conditions.
C 2. All information on the Package is provided to the Licensee on the
C    understanding that the details thereof are confidential.
C 3. All publications issued by the Licensee that include results obtained
C    with the help of one or more of the Packages shall acknowledge the
C    use of the Packages. The Licensee will notify the Numerical Analysis
C    Group at Rutherford Appleton Laboratory of any such publication.
C 4. The Packages may be modified by or on behalf of the Licensee
C    for such use in research applications but at no time shall such
C    Packages or modifications thereof become the property of the
C    Licensee. The Licensee shall make available free of charge to the
C    copyright holder for any purpose all information relating to
C    any modification.
C 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
C    direct or consequential loss or damage whatsoever arising out of
C    the use of Packages by the Licensee.
C *******************************************************************
C
C Original date July 1999
CCCCC PACKAGE MC64A/AD
CCCCC AUTHORS Iain Duff (i.duff@rl.ac.uk) and Jacko Koster (jak@ii.uib.no)

C 12th July 2004 Version 1.0.0. Version numbering added.

C 30/07/04  Version 1.1.0. Permutation array flagged negative to indicate
C           dependent columns in singular case.  Calls to MC64F changed
C           to avoid unsafe reference to array L.


      SUBROUTINE MC64ID(ICNTL)
      IMPLICIT NONE
      INTEGER ICNTL(10)
      INTEGER I
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = -1
      DO 10 I = 4,10
        ICNTL(I) = 0
   10 CONTINUE
      RETURN
      END
C**********************************************************************
      SUBROUTINE MC64AD(JOB,N,NE,IP,IRN,A,NUM,CPERM,LIW,IW,LDW,DW,
     &           ICNTL,INFO)
      IMPLICIT NONE
      INTEGER JOB,N,NE,NUM,LIW,LDW
      INTEGER IP(N+1),IRN(NE),CPERM(N),IW(LIW),ICNTL(10),INFO(10)
      DOUBLE PRECISION A(NE),DW(LDW)
      INTEGER I,J,K
      DOUBLE PRECISION FACT,ZERO,RINF
      PARAMETER (ZERO=0.0D+00)
      EXTERNAL FD05AD,MC21AD,MC64BD,MC64RD,MC64SD,MC64WD
      DOUBLE PRECISION FD05AD
      INTRINSIC ABS,LOG
      RINF = FD05AD(5)
      IF (JOB.LT.1 .OR. JOB.GT.5) THEN
        INFO(1) = -1
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
        GO TO 99
      ENDIF
      IF (N.LT.1) THEN
        INFO(1) = -2
        INFO(2) = N
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
        GO TO 99
      ENDIF
      IF (NE.LT.1) THEN
        INFO(1) = -3
        INFO(2) = NE
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'NE',NE
        GO TO 99
      ENDIF
      IF (JOB.EQ.1) K = 5*N
      IF (JOB.EQ.2) K = 4*N
      IF (JOB.EQ.3) K = 10*N + NE
      IF (JOB.EQ.4) K = 5*N
      IF (JOB.EQ.5) K = 5*N
      IF (LIW.LT.K) THEN
        INFO(1) = -4
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
        GO TO 99
      ENDIF
      IF (JOB.GT.1) THEN
        IF (JOB.EQ.2) K = N
        IF (JOB.EQ.3) K = NE
        IF (JOB.EQ.4) K = 2*N + NE
        IF (JOB.EQ.5) K = 3*N + NE
        IF (LDW.LT.K) THEN
          INFO(1) = -5
          INFO(2) = K
          IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
          GO TO 99
        ENDIF
      ENDIF
      IF (ICNTL(4).EQ.0) THEN
        DO 3 I = 1,N
          IW(I) = 0
    3   CONTINUE
        DO 6 J = 1,N
          DO 4 K = IP(J),IP(J+1)-1
            I = IRN(K)
            IF (I.LT.1 .OR. I.GT.N) THEN
              INFO(1) = -6
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),J,I
              GO TO 99
            ENDIF
            IF (IW(I).EQ.J) THEN
              INFO(1) = -7
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9007) INFO(1),J,I
              GO TO 99
            ELSE
              IW(I) = J
            ENDIF
    4     CONTINUE
    6   CONTINUE
      ENDIF
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9020) JOB,N,NE
        WRITE(ICNTL(3),9021) (IP(J),J=1,N+1)
        WRITE(ICNTL(3),9022) (IRN(J),J=1,NE)
        IF (JOB.GT.1) WRITE(ICNTL(3),9023) (A(J),J=1,NE)
      ENDIF
      DO 8 I=1,10
        INFO(I) = 0
    8 CONTINUE
      IF (JOB.EQ.1) THEN
        DO 10 J = 1,N
          IW(J) = IP(J+1) - IP(J)
   10   CONTINUE
        CALL MC21AD(N,IRN,NE,IP,IW(1),CPERM,NUM,IW(N+1))
        GO TO 90
      ENDIF
      IF (JOB.EQ.2) THEN
        CALL MC64BD(N,NE,IP,IRN,A,CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),DW)
        GO TO 90
      ENDIF
      IF (JOB.EQ.3) THEN
        DO 20 K = 1,NE
          IW(K) = IRN(K)
          DW(K) = ABS(A(K))
   20   CONTINUE
        CALL MC64RD(N,NE,IP,IW,DW)
        CALL MC64SD(N,NE,IP,IW(1),DW,CPERM,NUM,IW(NE+1),
     &     IW(NE+N+1),IW(NE+2*N+1),IW(NE+3*N+1),IW(NE+4*N+1),
     &     IW(NE+5*N+1),IW(NE+6*N+1))
        GO TO 90
      ENDIF
      IF (JOB.EQ.4) THEN
        DO 50 J = 1,N
          FACT = ZERO
          DO 30 K = IP(J),IP(J+1)-1
            IF (ABS(A(K)).GT.FACT) FACT = ABS(A(K))
   30     CONTINUE
          DO 40 K = IP(J),IP(J+1)-1
            DW(2*N+K) = FACT - ABS(A(K))
   40     CONTINUE
   50   CONTINUE
        CALL MC64WD(N,NE,IP,IRN,DW(2*N+1),CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(4*N+1),
     &     DW(1),DW(N+1))
        GO TO 90
      ENDIF
      IF (JOB.EQ.5) THEN
        DO 75 J = 1,N
          FACT = ZERO
          DO 60 K = IP(J),IP(J+1)-1
            DW(3*N+K) = ABS(A(K))
            IF (DW(3*N+K).GT.FACT) FACT = DW(3*N+K)
   60     CONTINUE
          DW(2*N+J) = FACT
          IF (FACT.NE.ZERO) THEN
            FACT = LOG(FACT)
          ELSE
            FACT = RINF/N
          ENDIF
          DO 70 K = IP(J),IP(J+1)-1
            IF (DW(3*N+K).NE.ZERO) THEN
              DW(3*N+K) = FACT - LOG(DW(3*N+K))
            ELSE
              DW(3*N+K) = RINF/N
            ENDIF
   70     CONTINUE
   75   CONTINUE
        CALL MC64WD(N,NE,IP,IRN,DW(3*N+1),CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(4*N+1),
     &     DW(1),DW(N+1))
        IF (NUM.EQ.N) THEN
          DO 80 J = 1,N
            IF (DW(2*N+J).NE.ZERO) THEN
              DW(N+J) = DW(N+J) - LOG(DW(2*N+J))
            ELSE
              DW(N+J) = ZERO
            ENDIF
   80     CONTINUE
        ENDIF
        FACT = 0.5*LOG(RINF)
        DO 86 J = 1,N
          IF (DW(J).LT.FACT .AND. DW(N+J).LT.FACT) GO TO 86
          INFO(1) = 2
          GO TO 90
   86   CONTINUE
      ENDIF
   90 IF (INFO(1).EQ.0 .AND. NUM.LT.N) THEN
        INFO(1) = 1
        IF (ICNTL(2).GE.0) WRITE(ICNTL(2),9011) INFO(1)
      ENDIF
      IF (INFO(1).EQ.2) THEN
        IF (ICNTL(2).GE.0) WRITE(ICNTL(2),9012) INFO(1)
      ENDIF
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9030) (INFO(J),J=1,2)
        WRITE(ICNTL(3),9031) NUM
        WRITE(ICNTL(3),9032) (CPERM(J),J=1,N)
        IF (JOB.EQ.5) THEN
          WRITE(ICNTL(3),9033) (DW(J),J=1,N)
          WRITE(ICNTL(3),9034) (DW(N+J),J=1,N)
        ENDIF
      ENDIF
   99 RETURN
 9001 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2,
     &        ' because ',(A),' = ',I10)
 9004 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        Column ',I8,
     &        ' contains an entry with invalid row index ',I8)
 9007 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        Column ',I8,
     &        ' contains two or more entries with row index ',I8)
 9011 FORMAT (' ****** Warning from MC64A/AD. INFO(1) = ',I2/
     &        '        The matrix is structurally singular.')
 9012 FORMAT (' ****** Warning from MC64A/AD. INFO(1) = ',I2/
     &        '        Some scaling factors may be too large.')
 9020 FORMAT (' ****** Input parameters for MC64A/AD:'/
     &        ' JOB = ',I8/' N   = ',I8/' NE  = ',I8)
 9021 FORMAT (' IP(1:N+1)  = ',8I8/(14X,8I8))
 9022 FORMAT (' IRN(1:NE)  = ',8I8/(14X,8I8))
 9023 FORMAT (' A(1:NE)    = ',4(1PD14.4)/(14X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for MC64A/AD:'/
     &        ' INFO(1:2)  = ',2I8)
 9031 FORMAT (' NUM        = ',I8)
 9032 FORMAT (' CPERM(1:N) = ',8I8/(14X,8I8))
 9033 FORMAT (' DW(1:N)    = ',5(F11.3)/(14X,5(F11.3)))
 9034 FORMAT (' DW(N+1:2N) = ',5(F11.3)/(14X,5(F11.3)))
      END
C**********************************************************************
      SUBROUTINE MC64BD(N,NE,IP,IRN,A,IPERM,NUM,JPERM,PR,Q,L,D)
      IMPLICIT NONE
      INTEGER N,NE,NUM
      INTEGER IP(N+1),IRN(NE),IPERM(N),JPERM(N),PR(N),Q(N),L(N)
      DOUBLE PRECISION A(NE),D(N)
      INTEGER I,II,J,JJ,JORD,Q0,QLEN,IDUM,JDUM,ISP,JSP,
     &        K,KK,KK1,KK2,I0,UP,LOW,LPOS
      DOUBLE PRECISION CSP,DI,DNEW,DQ0,AI,A0,BV
      DOUBLE PRECISION RINF,ZERO,MINONE
      PARAMETER (ZERO=0.0D+0,MINONE=-1.0D+0)
      INTRINSIC ABS,MIN
      EXTERNAL FD05AD,MC64DD,MC64ED,MC64FD
      DOUBLE PRECISION FD05AD
      RINF = FD05AD(5)
      NUM = 0
      BV = RINF
      DO 10 K = 1,N
        IPERM(K) = 0
        JPERM(K) = 0
        PR(K) = IP(K)
        D(K) = ZERO
   10 CONTINUE
      DO 20 J = 1,N
        A0 = MINONE
        DO 30 K = IP(J),IP(J+1)-1
          I = IRN(K)
          AI = ABS(A(K))
          IF (AI.GT.D(I)) D(I) = AI
          IF (JPERM(J).NE.0) GO TO 30
          IF (AI.GE.BV) THEN
            A0 = BV
            IF (IPERM(I).NE.0) GO TO 30
            JPERM(J) = I
            IPERM(I) = J
            NUM = NUM + 1
          ELSE
            IF (AI.LE.A0) GO TO 30
            A0 = AI
            I0 = I
          ENDIF
   30   CONTINUE
        IF (A0.NE.MINONE .AND. A0.LT.BV) THEN
          BV = A0
          IF (IPERM(I0).NE.0) GO TO 20
          IPERM(I0) = J
          JPERM(J) = I0
          NUM = NUM + 1
        ENDIF
   20 CONTINUE
      DO 25 I = 1,N
        BV = MIN(BV,D(I))
   25 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 95 J = 1,N
        IF (JPERM(J).NE.0) GO TO 95
        DO 50 K = IP(J),IP(J+1)-1
          I = IRN(K)
          AI = ABS(A(K))
          IF (AI.LT.BV) GO TO 50
          IF (IPERM(I).EQ.0) GO TO 90
          JJ = IPERM(I)
          KK1 = PR(JJ)
          KK2 = IP(JJ+1) - 1
          IF (KK1.GT.KK2) GO TO 50
          DO 70 KK = KK1,KK2
            II = IRN(KK)
            IF (IPERM(II).NE.0) GO TO 70
            IF (ABS(A(KK)).GE.BV) GO TO 80
   70     CONTINUE
          PR(JJ) = KK2 + 1
   50   CONTINUE
        GO TO 95
   80   JPERM(JJ) = II
        IPERM(II) = JJ
        PR(JJ) = KK + 1
   90   NUM = NUM + 1
        JPERM(J) = I
        IPERM(I) = J
        PR(J) = K + 1
   95 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 99 I = 1,N
        D(I) = MINONE
        L(I) = 0
   99 CONTINUE
      DO 100 JORD = 1,N
        IF (JPERM(JORD).NE.0) GO TO 100
        QLEN = 0
        LOW = N + 1
        UP = N + 1
        CSP = MINONE
        J = JORD
        PR(J) = -1
        DO 115 K = IP(J),IP(J+1)-1
          I = IRN(K)
          DNEW = ABS(A(K))
          IF (CSP.GE.DNEW) GO TO 115
          IF (IPERM(I).EQ.0) THEN
            CSP = DNEW
            ISP = I
            JSP = J
            IF (CSP.GE.BV) GO TO 160
          ELSE
            D(I) = DNEW
            IF (DNEW.GE.BV) THEN
              LOW = LOW - 1
              Q(LOW) = I
            ELSE
              QLEN = QLEN + 1
              L(I) = QLEN
              CALL MC64DD(I,N,Q,D,L,1)
            ENDIF
            JJ = IPERM(I)
            PR(JJ) = J
          ENDIF
  115   CONTINUE
        DO 150 JDUM = 1,NUM
          IF (LOW.EQ.UP) THEN
            IF (QLEN.EQ.0) GO TO 160
            I = Q(1)
            IF (CSP.GE.D(I)) GO TO 160
            BV = D(I)
            DO 152 IDUM = 1,N
              CALL MC64ED(QLEN,N,Q,D,L,1)
              L(I) = 0
              LOW = LOW - 1
              Q(LOW) = I
              IF (QLEN.EQ.0) GO TO 153
              I = Q(1)
              IF (D(I).NE.BV) GO TO 153
  152       CONTINUE
          ENDIF
  153     UP = UP - 1
          Q0 = Q(UP)
          DQ0 = D(Q0)
          L(Q0) = UP
          J = IPERM(Q0)
          DO 155 K = IP(J),IP(J+1)-1
            I = IRN(K)
            IF (L(I).GE.UP) GO TO 155
            DNEW = MIN(DQ0,ABS(A(K)))
            IF (CSP.GE.DNEW) GO TO 155
            IF (IPERM(I).EQ.0) THEN
              CSP = DNEW
              ISP = I
              JSP = J
              IF (CSP.GE.BV) GO TO 160
            ELSE
              DI = D(I)
              IF (DI.GE.BV .OR. DI.GE.DNEW) GO TO 155
              D(I) = DNEW
              IF (DNEW.GE.BV) THEN
                IF (DI.NE.MINONE) THEN
                  LPOS = L(I)
                  CALL MC64FD(LPOS,QLEN,N,Q,D,L,1)
                ENDIF
                L(I) = 0
                LOW = LOW - 1
                Q(LOW) = I
              ELSE
                IF (DI.EQ.MINONE) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL MC64DD(I,N,Q,D,L,1)
              ENDIF
              JJ = IPERM(I)
              PR(JJ) = J
            ENDIF
  155     CONTINUE
  150   CONTINUE
  160   IF (CSP.EQ.MINONE) GO TO 190
        BV = MIN(BV,CSP)
        NUM = NUM + 1
        I = ISP
        J = JSP
        DO 170 JDUM = 1,NUM+1
          I0 = JPERM(J)
          JPERM(J) = I
          IPERM(I) = J
          J = PR(J)
          IF (J.EQ.-1) GO TO 190
          I = I0
  170   CONTINUE
  190   DO 191 KK = UP,N
          I = Q(KK)
          D(I) = MINONE
          L(I) = 0
  191   CONTINUE
        DO 192 KK = LOW,UP-1
          I = Q(KK)
          D(I) = MINONE
  192   CONTINUE
        DO 193 KK = 1,QLEN
          I = Q(KK)
          D(I) = MINONE
          L(I) = 0
  193   CONTINUE
  100 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 300 J = 1,N
        JPERM(J) = 0
  300 CONTINUE
      K = 0
      DO 310 I = 1,N
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          PR(K) = I
        ELSE
          J = IPERM(I)
          JPERM(J) = I
        ENDIF
  310 CONTINUE
      K = 0
      DO 320 I = 1,N
        IF (JPERM(I).NE.0) GO TO 320
        K = K + 1
        JDUM = PR(K)
        IPERM(JDUM) = I
  320 CONTINUE
 1000 RETURN
      END
C**********************************************************************
      SUBROUTINE MC64DD(I,N,Q,D,L,IWAY)
      IMPLICIT NONE
      INTEGER I,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)
      INTEGER IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      DOUBLE PRECISION DI
      DI = D(I)
      POS = L(I)
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          IF (POS.LE.1) GO TO 20
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.LE.D(QK)) GO TO 20
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   10   CONTINUE
      ELSE
        DO 15 IDUM = 1,N
          IF (POS.LE.1) GO TO 20
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.GE.D(QK)) GO TO 20
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   15   CONTINUE
      ENDIF
   20 Q(POS) = I
      L(I) = POS
      RETURN
      END
C**********************************************************************
      SUBROUTINE MC64ED(QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
      INTEGER QLEN,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)
      INTEGER I,IDUM,K,POS,POSK
      PARAMETER (K=2)
      DOUBLE PRECISION DK,DR,DI
      I = Q(QLEN)
      DI = D(I)
      QLEN = QLEN - 1
      POS = 1
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 20
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.LT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.GE.DK) GO TO 20
          Q(POS) = Q(POSK)
          L(Q(POS)) = POS
          POS = POSK
   10   CONTINUE
      ELSE
        DO 15 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 20
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.GT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.LE.DK) GO TO 20
          Q(POS) = Q(POSK)
          L(Q(POS)) = POS
          POS = POSK
   15   CONTINUE
      ENDIF
   20 Q(POS) = I
      L(I) = POS
      RETURN
      END
C**********************************************************************
      SUBROUTINE MC64FD(POS0,QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
      INTEGER POS0,QLEN,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)
      INTEGER I,IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      DOUBLE PRECISION DK,DR,DI
      IF (QLEN.EQ.POS0) THEN
        QLEN = QLEN - 1
        RETURN
      ENDIF
      I = Q(QLEN)
      DI = D(I)
      QLEN = QLEN - 1
      POS = POS0
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          IF (POS.LE.1) GO TO 20
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.LE.D(QK)) GO TO 20
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   10   CONTINUE
   20   Q(POS) = I
        L(I) = POS
        DO 30 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 40
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.LT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.GE.DK) GO TO 40
          QK = Q(POSK)
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   30   CONTINUE
      ELSE
        DO 32 IDUM = 1,N
          IF (POS.LE.1) GO TO 34
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.GE.D(QK)) GO TO 34
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   32   CONTINUE
   34   Q(POS) = I
        L(I) = POS
        DO 36 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 40
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.GT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.LE.DK) GO TO 40
          QK = Q(POSK)
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   36   CONTINUE
      ENDIF
   40 Q(POS) = I
      L(I) = POS
      RETURN
      END
C**********************************************************************
      SUBROUTINE MC64RD(N,NE,IP,IRN,A)
      IMPLICIT NONE
      INTEGER N,NE
      INTEGER IP(N+1),IRN(NE)
      DOUBLE PRECISION A(NE)
      INTEGER THRESH,TDLEN
      PARAMETER (THRESH=15,TDLEN=50)
      INTEGER J,IPJ,K,LEN,R,S,HI,FIRST,MID,LAST,TD
      DOUBLE PRECISION HA,KEY
      INTEGER TODO(TDLEN)
      DO 100 J = 1,N
        LEN = IP(J+1) - IP(J)
        IF (LEN.LE.1) GO TO 100
        IPJ = IP(J)
        IF (LEN.LT.THRESH) GO TO 400
        TODO(1) = IPJ
        TODO(2) = IPJ + LEN
        TD = 2
  500   CONTINUE
        FIRST = TODO(TD-1)
        LAST = TODO(TD)
        KEY = A((FIRST+LAST)/2)
        DO 475 K = FIRST,LAST-1
          HA = A(K)
          IF (HA.EQ.KEY) GO TO 475
          IF (HA.GT.KEY) GO TO 470
          KEY = HA
          GO TO 470
  475   CONTINUE
        TD = TD - 2
        GO TO 425
  470   MID = FIRST
        DO 450 K = FIRST,LAST-1
          IF (A(K).LE.KEY) GO TO 450
          HA = A(MID)
          A(MID) = A(K)
          A(K) = HA
          HI = IRN(MID)
          IRN(MID) = IRN(K)
          IRN(K) = HI
          MID = MID + 1
  450   CONTINUE
        IF (MID-FIRST.GE.LAST-MID) THEN
          TODO(TD+2) = LAST
          TODO(TD+1) = MID
          TODO(TD) = MID
        ELSE
          TODO(TD+2) = MID
          TODO(TD+1) = FIRST
          TODO(TD) = LAST
          TODO(TD-1) = MID
        ENDIF
        TD = TD + 2
  425   CONTINUE
        IF (TD.EQ.0) GO TO 400
        IF (TODO(TD)-TODO(TD-1).GE.THRESH) GO TO 500
        TD = TD - 2
        GO TO 425
  400   DO 200 R = IPJ+1,IPJ+LEN-1
          IF (A(R-1) .LT. A(R)) THEN
            HA = A(R)
            HI = IRN(R)
            A(R) = A(R-1)
            IRN(R) = IRN(R-1)
            DO 300 S = R-1,IPJ+1,-1
              IF (A(S-1) .LT. HA) THEN
                A(S) = A(S-1)
                IRN(S) = IRN(S-1)
              ELSE
                A(S) = HA
                IRN(S) = HI
                GO TO 200
              END IF
  300       CONTINUE
            A(IPJ) = HA
            IRN(IPJ) = HI
          END IF
  200   CONTINUE
  100 CONTINUE
      RETURN
      END
C**********************************************************************
      SUBROUTINE MC64SD(N,NE,IP,IRN,A,IPERM,NUMX,
     &           W,LEN,LENL,LENH,FC,IW,IW4)
      IMPLICIT NONE
      INTEGER N,NE,NUMX
      INTEGER IP(N+1),IRN(NE),IPERM(N),
     &        W(N),LEN(N),LENL(N),LENH(N),FC(N),IW(N),IW4(4*N)
      DOUBLE PRECISION A(NE)
      INTEGER NUM,NVAL,WLEN,II,I,J,K,L,CNT,MOD,IDUM1,IDUM2,IDUM3
      DOUBLE PRECISION BVAL,BMIN,BMAX,RINF
      EXTERNAL FD05AD,MC64QD,MC64UD
      DOUBLE PRECISION FD05AD
      RINF = FD05AD(5)
      DO 20 J = 1,N
        FC(J) = J
        IW(J) = 0
        LEN(J) = IP(J+1) - IP(J)
   20 CONTINUE
      CNT = 1
      MOD = 1
      NUMX = 0
      CALL MC64UD(CNT,MOD,N,IRN,NE,IP,LEN,FC,IW,NUMX,N,
     &            IW4(1),IW4(N+1),IW4(2*N+1),IW4(3*N+1))
      NUM = NUMX
      IF (NUM.NE.N) THEN
        BMAX = RINF
      ELSE
        BMAX = RINF
        DO 30 J = 1,N
          BVAL = 0.0
          DO 25 K = IP(J),IP(J+1)-1
            IF (A(K).GT.BVAL) BVAL = A(K)
   25     CONTINUE
          IF (BVAL.LT.BMAX) BMAX = BVAL
   30   CONTINUE
        BMAX = 1.001 * BMAX
      ENDIF
      BVAL = 0.0
      BMIN = 0.0
      WLEN = 0
      DO 48 J = 1,N
        L = IP(J+1) - IP(J)
        LENH(J) = L
        LEN(J) = L
        DO 45 K = IP(J),IP(J+1)-1
          IF (A(K).LT.BMAX) GO TO 46
   45   CONTINUE
        K = IP(J+1)
   46   LENL(J) = K - IP(J)
        IF (LENL(J).EQ.L) GO TO 48
        WLEN = WLEN + 1
        W(WLEN) = J
   48 CONTINUE
      DO 90 IDUM1 = 1,NE
        IF (NUM.EQ.NUMX) THEN
          DO 50 I = 1,N
            IPERM(I) = IW(I)
   50     CONTINUE
          DO 80 IDUM2 = 1,NE
            BMIN = BVAL
            IF (BMAX .EQ. BMIN) GO TO 99
            CALL MC64QD(IP,LENL,LEN,W,WLEN,A,NVAL,BVAL)
            IF (NVAL.LE.1) GO TO 99
            K = 1
            DO 70 IDUM3 = 1,N
              IF (K.GT.WLEN) GO TO 71
              J = W(K)
              DO 55 II = IP(J)+LEN(J)-1,IP(J)+LENL(J),-1
                IF (A(II).GE.BVAL) GO TO 60
                I = IRN(II)
                IF (IW(I).NE.J) GO TO 55
                IW(I) = 0
                NUM = NUM - 1
                FC(N-NUM) = J
   55         CONTINUE
   60         LENH(J) = LEN(J)
              LEN(J) = II - IP(J) + 1
              IF (LENL(J).EQ.LENH(J)) THEN
                W(K) = W(WLEN)
                WLEN = WLEN - 1
              ELSE
                K = K + 1
              ENDIF
   70       CONTINUE
   71       IF (NUM.LT.NUMX) GO TO 81
   80     CONTINUE
   81     MOD = 1
        ELSE
          BMAX = BVAL
          CALL MC64QD(IP,LEN,LENH,W,WLEN,A,NVAL,BVAL)
          IF (NVAL.EQ.0. OR. BVAL.EQ.BMIN) GO TO 99
          K = 1
          DO 87 IDUM3 = 1,N
            IF (K.GT.WLEN) GO TO 88
            J = W(K)
            DO 85 II = IP(J)+LEN(J),IP(J)+LENH(J)-1
              IF (A(II).LT.BVAL) GO TO 86
   85       CONTINUE
   86       LENL(J) = LEN(J)
            LEN(J) = II - IP(J)
            IF (LENL(J).EQ.LENH(J)) THEN
              W(K) = W(WLEN)
              WLEN = WLEN - 1
            ELSE
              K = K + 1
            ENDIF
   87     CONTINUE
   88     MOD = 0
        ENDIF
        CNT = CNT + 1
        CALL MC64UD(CNT,MOD,N,IRN,NE,IP,LEN,FC,IW,NUM,NUMX,
     &            IW4(1),IW4(N+1),IW4(2*N+1),IW4(3*N+1))
   90 CONTINUE
   99 IF (NUMX.EQ.N) GO TO 1000
      DO 300 J = 1,N
        W(J) = 0
  300 CONTINUE
      K = 0
      DO 310 I = 1,N
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          IW(K) = I
        ELSE
          J = IPERM(I)
          W(J) = I
        ENDIF
  310 CONTINUE
      K = 0
      DO 320 J = 1,N
        IF (W(J).NE.0) GO TO 320
        K = K + 1
        IDUM1 = IW(K)
        IPERM(IDUM1) = J
  320 CONTINUE
 1000 RETURN
      END
C**********************************************************************
      SUBROUTINE MC64QD(IP,LENL,LENH,W,WLEN,A,NVAL,VAL)
      IMPLICIT NONE
      INTEGER WLEN,NVAL
      INTEGER IP(*),LENL(*),LENH(*),W(*)
      DOUBLE PRECISION A(*),VAL
      INTEGER XX,J,K,II,S,POS
      PARAMETER (XX=10)
      DOUBLE PRECISION SPLIT(XX),HA
      NVAL = 0
      DO 10 K = 1,WLEN
        J = W(K)
        DO 15 II = IP(J)+LENL(J),IP(J)+LENH(J)-1
          HA = A(II)
          IF (NVAL.EQ.0) THEN
            SPLIT(1) = HA
            NVAL = 1
          ELSE
            DO 20 S = NVAL,1,-1
              IF (SPLIT(S).EQ.HA) GO TO 15
              IF (SPLIT(S).GT.HA) THEN
                POS = S + 1
                GO TO 21
              ENDIF
  20        CONTINUE
            POS = 1
  21        DO 22 S = NVAL,POS,-1
              SPLIT(S+1) = SPLIT(S)
  22        CONTINUE
            SPLIT(POS) = HA
            NVAL = NVAL + 1
          ENDIF
          IF (NVAL.EQ.XX) GO TO 11
  15    CONTINUE
  10  CONTINUE
  11  IF (NVAL.GT.0) VAL = SPLIT((NVAL+1)/2)
      RETURN
      END
C**********************************************************************
      SUBROUTINE MC64UD(ID,MOD,N,IRN,LIRN,IP,LENC,FC,IPERM,NUM,NUMX,
     &           PR,ARP,CV,OUT)
      IMPLICIT NONE
      INTEGER ID,MOD,N,LIRN,NUM,NUMX
      INTEGER ARP(N),CV(N),IRN(LIRN),IP(N),
     &        FC(N),IPERM(N),LENC(N),OUT(N),PR(N)
      INTEGER I,II,IN1,IN2,J,J1,JORD,K,KK,LAST,NFC,
     &        NUM0,NUM1,NUM2,ID0,ID1
      IF (ID.EQ.1) THEN
        DO 5 I = 1,N
          CV(I) = 0
          ARP(I) = 0
    5   CONTINUE
        NUM1 = N
        NUM2 = N
      ELSE
        IF (MOD.EQ.1) THEN
          DO 8 I = 1,N
            ARP(I) = 0
    8     CONTINUE
        ENDIF
        NUM1 = NUMX
        NUM2 = N - NUMX
      ENDIF
      NUM0 = NUM
      NFC = 0
      ID0 = (ID-1)*N
      DO 100 JORD = NUM0+1,N
        ID1 = ID0 + JORD
        J = FC(JORD-NUM0)
        PR(J) = -1
        DO 70 K = 1,JORD
          IF (ARP(J).GE.LENC(J)) GO TO 30
          IN1 = IP(J) + ARP(J)
          IN2 = IP(J) + LENC(J) - 1
          DO 20 II = IN1,IN2
            I = IRN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
          ARP(J) = LENC(J)
   30     OUT(J) = LENC(J) - 1
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENC(J) - 1
            IN1 = IN2 - IN1
            DO 40 II = IN1,IN2
              I = IRN(II)
              IF (CV(I).EQ.ID1) GO TO 40
              J1 = J
              J = IPERM(I)
              CV(I) = ID1
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
   40       CONTINUE
   50       J1 = PR(J)
            IF (J1.EQ.-1) THEN
              NFC = NFC + 1
              FC(NFC) = J
              IF (NFC.GT.NUM2) THEN
                LAST = JORD
                GO TO 101
              ENDIF
              GO TO 100
            ENDIF
            J = J1
   60     CONTINUE
   70   CONTINUE
   80   IPERM(I) = J
        ARP(J) = II - IP(J) + 1
        NUM = NUM + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 95
          II = IP(J) + LENC(J) - OUT(J) - 2
          I = IRN(II)
          IPERM(I) = J
   90   CONTINUE
   95   IF (NUM.EQ.NUM1) THEN
          LAST = JORD
          GO TO 101
        ENDIF
  100 CONTINUE
      LAST = N
  101 DO 110 JORD = LAST+1,N
        NFC = NFC + 1
        FC(NFC) = FC(JORD-NUM0)
  110 CONTINUE
      RETURN
      END
C**********************************************************************
      SUBROUTINE MC64WD(N,NE,IP,IRN,A,IPERM,NUM,
     &           JPERM,OUT,PR,Q,L,U,D)
      IMPLICIT NONE
      INTEGER N,NE,NUM
      INTEGER IP(N+1),IRN(NE),IPERM(N),
     &        JPERM(N),OUT(N),PR(N),Q(N),L(N)
      DOUBLE PRECISION A(NE),U(N),D(N)
      INTEGER I,I0,II,J,JJ,JORD,Q0,QLEN,JDUM,ISP,JSP,
     &        K,K0,K1,K2,KK,KK1,KK2,UP,LOW,LPOS
      DOUBLE PRECISION CSP,DI,DMIN,DNEW,DQ0,VJ
      DOUBLE PRECISION RINF,ZERO
      PARAMETER (ZERO=0.0D+0)
      EXTERNAL FD05AD,MC64DD,MC64ED,MC64FD
      DOUBLE PRECISION FD05AD
      RINF = FD05AD(5)
      NUM = 0
      DO 10 K = 1,N
        U(K) = RINF
        D(K) = ZERO
        IPERM(K) = 0
        JPERM(K) = 0
        PR(K) = IP(K)
        L(K) = 0
   10 CONTINUE
      DO 30 J = 1,N
        DO 20 K = IP(J),IP(J+1)-1
          I = IRN(K)
          IF (A(K).GT.U(I)) GO TO 20
          U(I) = A(K)
          IPERM(I) = J
          L(I) = K
   20   CONTINUE
   30 CONTINUE
      DO 40 I = 1,N
        J = IPERM(I)
        IF (J.EQ.0) GO TO 40
        IPERM(I) = 0
        IF (JPERM(J).NE.0) GO TO 40
        NUM = NUM + 1
        IPERM(I) = J
        JPERM(J) = L(I)
   40 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 95 J = 1,N
        IF (JPERM(J).NE.0) GO TO 95
        K1 = IP(J)
        K2 = IP(J+1) - 1
        IF (K1.GT.K2) GO TO 95
        VJ = RINF
        DO 50 K = K1,K2
          I = IRN(K)
          DI = A(K) - U(I)
          IF (DI.GT.VJ) GO TO 50
          IF (DI.LT.VJ .OR. DI.EQ.RINF) GO TO 55
          IF (IPERM(I).NE.0 .OR. IPERM(I0).EQ.0) GO TO 50
   55     VJ = DI
          I0 = I
          K0 = K
   50   CONTINUE
        D(J) = VJ
        K = K0
        I = I0
        IF (IPERM(I).EQ.0) GO TO 90
        DO 60 K = K0,K2
          I = IRN(K)
          IF (A(K)-U(I).GT.VJ) GO TO 60
          JJ = IPERM(I)
          KK1 = PR(JJ)
          KK2 = IP(JJ+1) - 1
          IF (KK1.GT.KK2) GO TO 60
          DO 70 KK = KK1,KK2
            II = IRN(KK)
            IF (IPERM(II).GT.0) GO TO 70
            IF (A(KK)-U(II).LE.D(JJ)) GO TO 80
   70     CONTINUE
          PR(JJ) = KK2 + 1
   60   CONTINUE
        GO TO 95
   80   JPERM(JJ) = KK
        IPERM(II) = JJ
        PR(JJ) = KK + 1
   90   NUM = NUM + 1
        JPERM(J) = K
        IPERM(I) = J
        PR(J) = K + 1
   95 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 99 I = 1,N
        D(I) = RINF
        L(I) = 0
   99 CONTINUE
      DO 100 JORD = 1,N
        IF (JPERM(JORD).NE.0) GO TO 100
        DMIN = RINF
        QLEN = 0
        LOW = N + 1
        UP = N + 1
        CSP = RINF
        J = JORD
        PR(J) = -1
        DO 115 K = IP(J),IP(J+1)-1
          I = IRN(K)
          DNEW = A(K) - U(I)
          IF (DNEW.GE.CSP) GO TO 115
          IF (IPERM(I).EQ.0) THEN
            CSP = DNEW
            ISP = K
            JSP = J
          ELSE
            IF (DNEW.LT.DMIN) DMIN = DNEW
            D(I) = DNEW
            QLEN = QLEN + 1
            Q(QLEN) = K
          ENDIF
  115   CONTINUE
        Q0 = QLEN
        QLEN = 0
        DO 120 KK = 1,Q0
          K = Q(KK)
          I = IRN(K)
          IF (CSP.LE.D(I)) THEN
            D(I) = RINF
            GO TO 120
          ENDIF
          IF (D(I).LE.DMIN) THEN
            LOW = LOW - 1
            Q(LOW) = I
            L(I) = LOW
          ELSE
            QLEN = QLEN + 1
            L(I) = QLEN
            CALL MC64DD(I,N,Q,D,L,2)
          ENDIF
          JJ = IPERM(I)
          OUT(JJ) = K
          PR(JJ) = J
  120   CONTINUE
        DO 150 JDUM = 1,NUM
          IF (LOW.EQ.UP) THEN
            IF (QLEN.EQ.0) GO TO 160
            I = Q(1)
            IF (D(I).GE.CSP) GO TO 160
            DMIN = D(I)
  152       CALL MC64ED(QLEN,N,Q,D,L,2)
            LOW = LOW - 1
            Q(LOW) = I
            L(I) = LOW
            IF (QLEN.EQ.0) GO TO 153
            I = Q(1)
            IF (D(I).GT.DMIN) GO TO 153
            GO TO 152
          ENDIF
  153     Q0 = Q(UP-1)
          DQ0 = D(Q0)
          IF (DQ0.GE.CSP) GO TO 160
          UP = UP - 1
          J = IPERM(Q0)
          VJ = DQ0 - A(JPERM(J)) + U(Q0)
          DO 155 K = IP(J),IP(J+1)-1
            I = IRN(K)
            IF (L(I).GE.UP) GO TO 155
            DNEW = VJ + A(K)-U(I)
            IF (DNEW.GE.CSP) GO TO 155
            IF (IPERM(I).EQ.0) THEN
              CSP = DNEW
              ISP = K
              JSP = J
            ELSE
              DI = D(I)
              IF (DI.LE.DNEW) GO TO 155
              IF (L(I).GE.LOW) GO TO 155
              D(I) = DNEW
              IF (DNEW.LE.DMIN) THEN
                LPOS = L(I)
                IF (LPOS.NE.0)
     *            CALL MC64FD(LPOS,QLEN,N,Q,D,L,2)
                LOW = LOW - 1
                Q(LOW) = I
                L(I) = LOW
              ELSE
                IF (L(I).EQ.0) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL MC64DD(I,N,Q,D,L,2)
              ENDIF
              JJ = IPERM(I)
              OUT(JJ) = K
              PR(JJ) = J
            ENDIF
  155     CONTINUE
  150   CONTINUE
  160   IF (CSP.EQ.RINF) GO TO 190
        NUM = NUM + 1
        I = IRN(ISP)
        IPERM(I) = JSP
        JPERM(JSP) = ISP
        J = JSP
        DO 170 JDUM = 1,NUM
          JJ = PR(J)
          IF (JJ.EQ.-1) GO TO 180
          K = OUT(J)
          I = IRN(K)
          IPERM(I) = JJ
          JPERM(JJ) = K
          J = JJ
  170   CONTINUE
  180   DO 185 KK = UP,N
          I = Q(KK)
          U(I) = U(I) + D(I) - CSP
  185   CONTINUE
  190   DO 191 KK = LOW,N
          I = Q(KK)
          D(I) = RINF
          L(I) = 0
  191   CONTINUE
        DO 193 KK = 1,QLEN
          I = Q(KK)
          D(I) = RINF
          L(I) = 0
  193   CONTINUE
  100 CONTINUE
 1000 DO 200 J = 1,N
        K = JPERM(J)
        IF (K.NE.0) THEN
          D(J) = A(K) - U(IRN(K))
        ELSE
          D(J) = ZERO
        ENDIF
        IF (IPERM(J).EQ.0) U(J) = ZERO
  200 CONTINUE
      IF (NUM.EQ.N) GO TO 1100
      DO 300 J = 1,N
        JPERM(J) = 0
  300 CONTINUE
      K = 0
      DO 310 I = 1,N
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          OUT(K) = I
        ELSE
          J = IPERM(I)
          JPERM(J) = I
        ENDIF
  310 CONTINUE
      K = 0
      DO 320 J = 1,N
        IF (JPERM(J).NE.0) GO TO 320
        K = K + 1
        JDUM = OUT(K)
        IPERM(JDUM) = - J
  320 CONTINUE
 1100 RETURN
      END
* *******************************************************************
* COPYRIGHT (c) 1977 Hyprotech UK
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of a UK ACADEMIC
* Licence, see http://hsl.rl.ac.uk/acuk/cou.html
*
* Please note that for a UK ACADEMIC Licence:
*
* 1. The Packages may only be used for academic research or teaching
*    purposes by the Licensee, and must not be copied by the Licensee for
*    use by any other persons. Use of the Packages in any commercial
*    application shall be subject to prior written agreement between
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
* Original date 8 Oct 1992
C######8/10/92 Toolpack tool decs employed.
C######8/10/92 D version created by name change only.
C 13/3/02 Cosmetic changes applied to reduce single/double differences
C
C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC21AD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW)
      INTEGER LICN,N,NUMNZ
      INTEGER ICN(LICN),IP(N),IPERM(N),IW(N,4),LENR(N)
      EXTERNAL MC21BD
      CALL MC21BD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW(1,1),IW(1,2),
     +            IW(1,3),IW(1,4))
      RETURN
      END
      SUBROUTINE MC21BD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,PR,ARP,CV,OUT)
      INTEGER LICN,N,NUMNZ
      INTEGER ARP(N),CV(N),ICN(LICN),IP(N),IPERM(N),LENR(N),OUT(N),PR(N)
      INTEGER I,II,IN1,IN2,IOUTK,J,J1,JORD,K,KK
      DO 10 I = 1,N
        ARP(I) = LENR(I) - 1
        CV(I) = 0
        IPERM(I) = 0
   10 CONTINUE
      NUMNZ = 0
      DO 100 JORD = 1,N
        J = JORD
        PR(J) = -1
        DO 70 K = 1,JORD
          IN1 = ARP(J)
          IF (IN1.LT.0) GO TO 30
          IN2 = IP(J) + LENR(J) - 1
          IN1 = IN2 - IN1
          DO 20 II = IN1,IN2
            I = ICN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
          ARP(J) = -1
   30     CONTINUE
          OUT(J) = LENR(J) - 1
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENR(J) - 1
            IN1 = IN2 - IN1
            DO 40 II = IN1,IN2
              I = ICN(II)
              IF (CV(I).EQ.JORD) GO TO 40
              J1 = J
              J = IPERM(I)
              CV(I) = JORD
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
   40       CONTINUE
   50       CONTINUE
            J = PR(J)
            IF (J.EQ.-1) GO TO 100
   60     CONTINUE
   70   CONTINUE
   80   CONTINUE
        IPERM(I) = J
        ARP(J) = IN2 - II - 1
        NUMNZ = NUMNZ + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 100
          II = IP(J) + LENR(J) - OUT(J) - 2
          I = ICN(II)
          IPERM(I) = J
   90   CONTINUE
  100 CONTINUE
      IF (NUMNZ.EQ.N) RETURN
      DO 110 I = 1,N
        ARP(I) = 0
  110 CONTINUE
      K = 0
      DO 130 I = 1,N
        IF (IPERM(I).NE.0) GO TO 120
        K = K + 1
        OUT(K) = I
        GO TO 130
  120   CONTINUE
        J = IPERM(I)
        ARP(J) = I
  130 CONTINUE
      K = 0
      DO 140 I = 1,N
        IF (ARP(I).NE.0) GO TO 140
        K = K + 1
        IOUTK = OUT(K)
        IPERM(IOUTK) = I
  140 CONTINUE
      RETURN
      END
C       Toolpack tool decs employed.
C       Arg dimension set to *.
C
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
C
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
C     .. Scalar Arguments ..
      INTEGER INCX,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DMAX
      INTEGER I,IX
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS
C     ..
C     .. Executable Statements ..
C
      IDAMAX = 0
      IF (N.LT.1) RETURN
      IDAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      DMAX = DABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
        IF (DABS(DX(IX)).LE.DMAX) GO TO 5
        IDAMAX = I
        DMAX = DABS(DX(IX))
    5   IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
        IF (DABS(DX(I)).LE.DMAX) GO TO 30
        IDAMAX = I
        DMAX = DABS(DX(I))
   30 CONTINUE
      RETURN

      END
C       A BLAS routine modified for use with HSL
C       Toolpack tool decs employed.
C       Contained comment lines which caused Decs to fail.
C
C      SUBROUTINE XERBLA(SRNAME,INFO)
C     .. Scalar Arguments ..
C      INTEGER INFO
C      CHARACTER SRNAME*6
C     ..
C
C  Purpose
C  =======
C
C  XERBLA  is an error handler for the Level 2 BLAS routines.
C
C  It is called by the Level 2 BLAS routines if an input parameter is
C  invalid.
C
C  Installers should consider modifying the STOP statement in order to
C  call system-specific exception-handling facilities.
C
C  Parameters
C  ==========
C
C  SRNAME - CHARACTER*6.
C           On entry, SRNAME specifies the name of the routine which
C           called XERBLA.
C
C  INFO   - INTEGER.
C           On entry, INFO specifies the position of the invalid
C           parameter in the parameter-list of the calling routine.
C
C
C  Auxiliary routine for Level 2 Blas.
C
C  Written on 20-July-1986.
C
C     .. Executable Statements ..
C
C      WRITE (*,FMT=99999) SRNAME,INFO
C
C      STOP
C
C99999 FORMAT (' ** On entry to ',A6,' parameter number ',I2,
C     +       ' had an illegal value')
C
C     End of XERBLA.
C
C      END
      LOGICAL FUNCTION LSAME ( CA, CB )
*     .. Scalar Arguments ..
      CHARACTER*1            CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME  tests if CA is the same letter as CB regardless of case.
*
*  N.B. This version of the routine is only correct for ASCII code.
*       Installers must modify the routine for other character-codes.
*
*       For EBCDIC systems the constant IOFF must be changed to -64.
*       For CDC systems using 6-12 bit representations, the system-
*       specific code in comments must be activated.
*
*  Parameters
*  ==========
*
*  CA     - CHARACTER*1
*  CB     - CHARACTER*1
*           On entry, CA and CB specify characters to be compared.
*           Unchanged on exit.
*
*
*  Auxiliary routine for Level 2 Blas.
*
*  -- Written on 11-October-1988.
*     Richard Hanson, Sandia National Labs.
*     Jeremy Du Croz, Nag Central Office.
*
*     .. Parameters ..
      INTEGER                IOFF
      PARAMETER            ( IOFF=32 )
*     .. Intrinsic Functions ..
      INTRINSIC              ICHAR
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA .EQ. CB
*
*     Now test for equivalence
*
      IF ( .NOT.LSAME ) THEN
         LSAME = ICHAR(CA) - IOFF .EQ. ICHAR(CB)
      END IF
      IF ( .NOT.LSAME ) THEN
         LSAME = ICHAR(CA) .EQ. ICHAR(CB) - IOFF
      END IF
*
      RETURN
*
*  The following comments contain code for CDC systems using 6-12 bit
*  representations.
*
*     .. Parameters ..
*     INTEGER                ICIRFX
*     PARAMETER            ( ICIRFX=62 )
*     .. Scalar Arguments ..
*     CHARACTER*1            CB
*     .. Array Arguments ..
*     CHARACTER*1            CA(*)
*     .. Local Scalars ..
*     INTEGER                IVAL
*     .. Intrinsic Functions ..
*     INTRINSIC              ICHAR, CHAR
*     .. Executable Statements ..
*
*     See if the first character in string CA equals string CB.
*
*     LSAME = CA(1) .EQ. CB .AND. CA(1) .NE. CHAR(ICIRFX)
*
*     IF (LSAME) RETURN
*
*     The characters are not identical. Now check them for equivalence.
*     Look for the 'escape' character, circumflex, followed by the
*     letter.
*
*     IVAL = ICHAR(CA(2))
*     IF (IVAL.GE.ICHAR('A') .AND. IVAL.LE.ICHAR('Z')) THEN
*        LSAME = CA(1) .EQ. CHAR(ICIRFX) .AND. CA(2) .EQ. CB
*     END IF
*
*     RETURN
*
*     End of LSAME.
*
      END
      SUBROUTINE MA57LF(N,FACT,LFACT,IFACT,LIFACT,IRN,JCN,FL,NNZ,
     *                  ID,JD,D,NNZD,IVP,IPERM,S64,INFO,ICNTL)
C This subroutine extract the factors
C     stored in FACT/IFACT by MA57BD.
C 
      INTEGER N,LFACT,LIFACT,NNZ,NNZD,ISCAL
      DOUBLE PRECISION FACT(LFACT),FL(LFACT),D(2*N),S64(N)
      INTEGER IFACT(LIFACT),IRN(LFACT),JCN(LFACT),ID(2*N),JD(2*N)
      INTEGER IPERM(N),IVP(N),INFO(40),ICNTL(20)
C N   must be set to the order of the matrix. It is not altered.
C FACT   must be set to hold the real values corresponding to the
C     factors. This must be unchanged since the preceding call to
C     MA57BD. It is not altered.
C LFACT  length of array FACT. It is not altered.
C IFACT  holds the integer indexing information for the matrix factors
C     in FACT. This must be unchanged since the preceding call to
C     MA57BD. It is not altered.
C LIFACT length of array IFACT. It is not altered.
C IRN on output, will hold the row indeces of the matrix L.
C ICN on output, will hold the column indeces of the matrix L.
C FL  on output, will hold the values of the matrix L.
C     sides corresponding to current block pivotal rows.
C NNZ is the number of nonzeros in L
C D   on output, will hold the value of the pivot (1x1 and 2x2)
C ID  on output, will hold the i rows coordinates of D
C JD  on output, will hold the j column coordinates of D
C IVP need not be set on entry. On exit IVP(IPERM(I)) = I
C IPERM on output holds the permutation of the matrix A
C S64 is the scaling vector produce by mc64
C INFO holds the informations (see spec of ma57)
C ICNTL holds the control parameters (see spec of ma57)
C Procedures
      INTRINSIC ABS
C Local variables
      INTEGER APOS,I,IBLK,II,IWPOS,J,K,K1,NCOLS,NROWS,IND,KD,APOS2
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
C APOS  Current position in array FACT.
C I     Temporary DO index
C IBLK  Index of block pivot.
C II    Temporary index.
C IWPOS Position in IFACT of start of current index list.
C J     Temporary DO index
C K     Temporary pointer to position in real array.
C NCOLS Number of columns in the block pivot.
C NROWS Number of rows in the block pivot.
      DO 10 I = 1,N
        IPERM(I) = 0
 10   CONTINUE
      APOS = 1
      APOS2 = 1
      IND = 0
      IWPOS = 4
      K = 0
      DO 270 IBLK = 1,IFACT(3)

C Find the number of rows and columns in the block.
        NCOLS = IFACT(IWPOS)
        NROWS = IFACT(IWPOS+1)
        IWPOS = IWPOS + 2
C Build the permutation
        DO 20 I = 1,NROWS
          II = IFACT(IWPOS+I-1)
          IPERM(I+IND) = II
 20     CONTINUE
        IWPOS = IWPOS + NCOLS
        IND = IND + NROWS
        APOS2 = APOS2 + (NROWS*(NROWS+1))/2 + NROWS*(NCOLS-NROWS)
  270 CONTINUE
C Build the inverse permutation
      DO 170 K1 = 1,N
        IVP(ABS(IPERM(K1))) = K1
  170 CONTINUE
      APOS = 1
      IND = 0
      IWPOS = 4
      K = 1
      KD = 1
      DO 370 IBLK = 1,IFACT(3)
C Find the number of rows and columns in the block.
        NCOLS = IFACT(IWPOS)
        NROWS = IFACT(IWPOS+1)
        IWPOS = IWPOS + 2
C Treat diagonal block and build L and D^(-1)
        DO 30 J = 1,NROWS
          DO 40 I = J,NROWS
             IF (I.EQ.J) THEN
               FL(K) = ONE
               IRN(K) = J + IND
               JCN(K) = J + IND
               IF ((I+IND).LE.INFO(25)) THEN
                 D(KD) = FACT(APOS)
                 ID(KD) = I + IND
                 JD(KD) = I + IND
                 KD = KD + 1
               ENDIF
C Build the matrix D^(-1) structure as a nonsymmetric matrix
               IF (IPERM(I+IND).LT.0) THEN
                 D(KD) = FACT(APOS2)
                 ID(KD) = I + IND
                 JD(KD) = I + 1 + IND
                 D(KD+1) = FACT(APOS2)
                 ID(KD+1) = JD(KD)
                 JD(KD+1) = ID(KD)
                 KD = KD + 2
                 APOS2 = APOS2 + 1
                 IPERM(I+IND) = ABS(IPERM(I+IND))
                 IPERM(I+1+IND) = ABS(IPERM(I+1+IND))
               ENDIF
             ELSE
               FL(K) = FACT(APOS)
               IRN(K) = J + IND
               JCN(K) = I + IND
             ENDIF
             APOS = APOS + 1
             K = K + 1
 40       CONTINUE
 30     CONTINUE
C Treat off-diagonal block
        IF (NCOLS.GT.NROWS) THEN
          DO 60 I = 1,NROWS
            DO 50 J = 1,NCOLS-NROWS
              II = IFACT(IWPOS+J-1+NROWS)
              IP = IVP(ABS(II))
              FL(K) = -FACT(APOS)
              IRN(K) = I + IND
              JCN(K) = ABS(IP)
              APOS = APOS + 1
              K = K + 1
 50         CONTINUE
 60       CONTINUE
        ENDIF
        IWPOS = IWPOS + NCOLS
        IND = IND + NROWS
  370 CONTINUE
      NNZ = K - 1
      NNZD = KD - 1
      IF (ICNTL(15) .EQ. 1) THEN
        ISCAL = LFACT - N
        DO 470 K1 = 1,N
          S64(K1) = FACT(ISCAL + K1 -1)
  470   CONTINUE
      ELSE
        DO 471 K1 = 1,N
          S64(K1) = 1.0d0
  471   CONTINUE
      ENDIF
c$$$      DO 1000 K1=1,NNZ
c$$$      WRITE(6,'(2I5XE20.14)') IRN(K1),JCN(K1),FL(K1)
c$$$ 1000 CONTINUE
c$$$      WRITE(6,'(a)') '----------------------'
c$$$      WRITE(6,*) (IPERM(I),I=1,N)
c$$$      WRITE(6,*) (IVP(I),I=1,N)
c$$$       WRITE(6,'(a)') '----------------------'
c$$$      DO 2000 K1 =1,NNZD
c$$$      WRITE(6,'(2i4xE20.14)') ID(K1),JD(K1),D(K1)
c$$$ 2000 CONTINUE
c$$$       WRITE(6,'(a)') '----------------------'
c$$$       WRITE(6,*) LFACT
c$$$       WRITE(6,'(a)') '----------------------'
c$$$       WRITE(6,'(E20.14)') (S64(K1),K1=1,N)
c$$$c       DO 1000 K1=1,lfact
c$$$c       WRITE(6,'(5E12.4)') fact(k1)
c$$$c  1000 CONTINUE

      END
      DOUBLE PRECISION FUNCTION FD05AD(INUM)
C----------------------------------------------------------------
C  Real constants for: IEEE double precision (8-byte arithmetic)
C
C  Obtained from H.S.L. subroutine ZE02AM.
C  Nick Gould and Sid Marlow, Harwell Laboratory, April 1988.
C----------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER INUM
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DC(5)
C     ..
C     .. Save statement ..
      SAVE DC
C     ..
C     .. Data statements ..
C
C  DC(1) THE SMALLEST POSITIVE NUMBER: 1.0 + DC(1) > 1.0.
C  DC(2) THE SMALLEST POSITIVE NUMBER: 1.0 - DC(2) < 1.0.
C  DC(3) THE SMALLEST NONZERO +VE REAL NUMBER.
C  DC(4) THE SMALLEST FULL PRECISION +VE REAL NUMBER.
C  DC(5) THE LARGEST FINITE +VE REAL NUMBER.
C
      DATA DC(1)/2.2204460492504D-16/
      DATA DC(2)/1.1102230246253D-16/
C     DATA DC(3)/4.9406564584126D-324/
      DATA DC(4)/2.2250738585073D-308/
      DATA DC(5)/1.7976931348622D+308/
C     ..
C     .. Executable Statements ..

      IF ( INUM .LE. 0 ) THEN
         FD05AD = DC( 1 )
      ELSE IF ( INUM .GE. 6 ) THEN
         FD05AD = DC( 5 )
      ELSE IF ( INUM .EQ. 3 ) THEN
         FD05AD = DC(4)/2.0D0**52
      ELSE
         FD05AD = DC( INUM )
      ENDIF
      RETURN
      END

C       Toolpack tool decs employed.
C       SAVE statements added.
C 16th October 2002: STOP and WRITE statements removed.
C
      INTEGER FUNCTION ID05AD(INUM)
C----------------------------------------------------------------
C  Integer constants for: IEEE double precision (8-byte arithmetic).
C
C  Obtained from H.S.L. subroutine ZE02AM.
C  Nick Gould and Sid Marlow, Harwell Laboratory, April 1988.
C----------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER INUM
C     ..
C     .. Local Arrays ..
      INTEGER IC(10)
C     ..
C     .. Save statement ..
      SAVE IC
C     ..
C     .. Data statements ..
C  IC(1) THE BASE (RADIX) OF THE FLOATING-POINT ARITHMETIC.
C  IC(2) THE NUMBER OF BASE IC(1) DIGITS IN THE SIGNIFICAND.
C  IC(3) THE NUMBER OF BITS USED FOR THE EXPONENT
C  IC(4) = 0 FLOATING-POINT ADDITION CHOPS, = 1 IT ROUNDS.
C  IC(5) = 0 A GUARD DIGIT IS NOT USED FOR *, = 1 IT IS.
C  IC(6) LARGEST -VE INTEGER:1.0 + DBLE(IC(1))**IC(6) > 1.0.
C  IC(7) LARGEST -VE INTEGER:1.0 - DBLE(IC(1))**IC(7) < 1.0.
C  IC(8) LARGEST -VE INTEGER: DBLE(IC(1))**IC(8) > 0.0.
C  IC(9) LARGEST -VE INTEGER: REAL(IC(1))**IC(9) IS NORMAL.
C  IC(10) LARGEST +VE INTEGER: REAL(IC(1))**IC(10) FINITE.
C
      DATA IC(1)/2/
      DATA IC(2)/53/
      DATA IC(3)/11/
      DATA IC(4)/1/
      DATA IC(5)/1/
      DATA IC(6)/-52/
      DATA IC(7)/-53/
      DATA IC(8)/-1074/
      DATA IC(9)/-1022/
      DATA IC(10)/1023/
C     ..
C     .. Executable Statements ..

      IF ( INUM .LE. 0 ) THEN
        ID05AD = IC( 1 )
      ELSE IF ( INUM .GE. 11 ) THEN
        ID05AD = IC( 10 )
      ELSE
        ID05AD = IC( INUM )
      ENDIF
      RETURN
      END
