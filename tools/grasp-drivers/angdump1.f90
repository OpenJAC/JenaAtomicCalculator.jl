!***********************************************************************
!   ANGDUMP1 -- dump the ONE-particle spin-angular coefficients that
!   GRASP2018's librang90 produces, for every CSF pair of a loaded list.
!
!   Reads  rcsf.c  (standard GRASP CSL file); writes a flat table to
!   stdout, one line per non-zero coefficient:
!
!       RANK  ICSF  JCSF  IA  IB  NP  NAK  COEFF
!
!   Built to compare against JAC's SpinAngular / SpinAngularNew.
!   See work/diag-grasp-angular.jl, which drives this.
!***********************************************************************
      PROGRAM ANGDUMP1
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNW
      USE orb_C,           ONLY: NW, NAK, NP, NCF
      USE onescalar_I
      USE oneparticlejj_I
      USE itjpo_I
      USE ispar_I
      IMPLICIT NONE

      INTEGER      :: NCORE, IC, IR, IA, IB, I, KT, IPT, MAXRANK
      REAL(DOUBLE) :: TSHELL(NNNW)
      REAL(DOUBLE), PARAMETER :: CUT = 1.0D-14
!  ... SETCSLA takes CHARACTER(LEN=24) and builds the file name via
!      INDEX(NAME,' '); a length-4 literal gives INDEX = 0, an EMPTY name,
!      and it then silently reads fort.21. The argument MUST be padded.
      CHARACTER(LEN=24) :: FNAME

      MAXRANK = 3

!  ... machine constants, physical constants and factorials. FACTT is
!      required by the angular routines and is easy to forget.
      CALL SETMC
      CALL SETCON
      CALL FACTT

      FNAME = 'rcsf'
      CALL SETCSLA(FNAME, NCORE)

      WRITE(*,'(A,I5,A,I5)') '# NCF = ', NCF, '   NW = ', NW
      DO I = 1, NW
         WRITE(*,'(A,I3,A,I3,A,I4)') '# subshell ', I, '  NP = ', NP(I), '  NAK = ', NAK(I)
      END DO
!  ... each CSF's 2J+1 and parity, so that a caller can match GRASP's CSF ordering to its own without
!      parsing the CSL file's coupling tree by hand.
      DO IC = 1, NCF
         WRITE(*,'(A,I5,A,I4,A,I3)') '# csf ', IC, '  ITJPO = ', ITJPO(IC), '  ISPAR = ', ISPAR(IC)
      END DO
      WRITE(*,'(A)') '# RANK ICSF JCSF IA IB NP NAK COEFF'

!  ... RANK 0 : ONESCALAR
      DO IC = 1, NCF
         DO IR = 1, NCF
            TSHELL = 0.0D0
            CALL ONESCALAR(IC, IR, IA, IB, TSHELL)
            IF (IA == 0) CYCLE
            IF (IA == IB) THEN
               DO I = 1, NW
                  IF (ABS(TSHELL(I)) > CUT)                                     &
                     WRITE(*,'(I4,2I5,2I4,I4,I5,ES26.17)')                      &
                           0, IC, IR, I, I, NP(I), NAK(I), TSHELL(I)
               END DO
            ELSE
               IF (ABS(TSHELL(1)) > CUT)                                        &
                  WRITE(*,'(I4,2I5,2I4,I4,I5,ES26.17)')                         &
                        0, IC, IR, IA, IB, NP(IA), NAK(IA), TSHELL(1)
            END IF
         END DO
      END DO

!  ... RANK > 0 : ONEPARTICLEJJ; IPT is the operator parity, as rhfs90 uses it
      DO KT = 1, MAXRANK
         IPT = 1
         DO IC = 1, NCF
            DO IR = 1, NCF
               TSHELL = 0.0D0
               CALL ONEPARTICLEJJ(KT, IPT, IC, IR, IA, IB, TSHELL)
               IF (IA == 0) CYCLE
               IF (IA == IB) THEN
                  DO I = 1, NW
                     IF (ABS(TSHELL(I)) > CUT)                                  &
                        WRITE(*,'(I4,2I5,2I4,I4,I5,ES26.17)')                   &
                              KT, IC, IR, I, I, NP(I), NAK(I), TSHELL(I)
                  END DO
               ELSE
                  IF (ABS(TSHELL(1)) > CUT)                                     &
                     WRITE(*,'(I4,2I5,2I4,I4,I5,ES26.17)')                      &
                           KT, IC, IR, IA, IB, NP(IA), NAK(IA), TSHELL(1)
               END IF
            END DO
         END DO
      END DO

      END PROGRAM ANGDUMP1
