!***********************************************************************
!   ANGDUMP2 -- dump the TWO-particle (electron-electron) spin-angular
!   coefficients that GRASP2018's librang90 produces, for every CSF pair
!   of a loaded list.  Reads rcsf.c; writes one line per coefficient:
!
!       ICSF  JCSF  IA1 IB1 IA2 IB2  K   COEFF
!
!   Takes INCOR as its first command-line argument (default 1 = include the core).
!
!   RKCO_GG deposits its results in BUFFER_C via SPEAK, so the buffer is
!   simply zeroed before each pair and read out afterwards.
!***********************************************************************
      PROGRAM ANGDUMP2
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNW
      USE orb_C,           ONLY: NW, NAK, NP, NCF
      USE buffer_C,        ONLY: LABEL, COEFF, NVCOEF
      USE cord_I
      USE alcbuf_I
      IMPLICIT NONE

      INTEGER      :: NCORE, IC, IR, I, INCOR, IARGLEN, IST
      CHARACTER(LEN=24) :: FNAME
      CHARACTER(LEN=8)  :: ARG
      REAL(DOUBLE), PARAMETER :: CUT = 1.0D-14

!  ... INCOR gates the CLOSED-SHELL (core) contributions: 0 = valence-valence only, as GRASP's own mcp uses
!      because it treats the core separately; 1 = the complete set. Default 1, overridden by the first argument.
      INCOR = 1
      CALL GET_COMMAND_ARGUMENT(1, ARG, IARGLEN, IST)
      IF (IST == 0 .AND. IARGLEN > 0) READ(ARG,*) INCOR

      CALL SETMC
      CALL SETCON
      CALL FACTT
      FNAME = 'rcsf'
      CALL SETCSLA(FNAME, NCORE)
      CALL ALCBUF(1)

      WRITE(*,'(A,I5,A,I5)') '# NCF = ', NCF, '   NW = ', NW
      DO I = 1, NW
         WRITE(*,'(A,I3,A,I3,A,I4)') '# subshell ', I, '  NP = ', NP(I), '  NAK = ', NAK(I)
      END DO
      WRITE(*,'(A)') '# ICSF JCSF IA1 IB1 IA2 IB2  K  COEFF'

      DO IC = 1, NCF
         DO IR = 1, NCF
            NVCOEF = 0
            CALL RKCO_GG(IC, IR, CORD, INCOR, 1)
            DO I = 1, NVCOEF
               IF (ABS(COEFF(I)) > CUT)                                        &
                  WRITE(*,'(2I5,4I4,I4,ES26.17)') IC, IR, LABEL(1,I),          &
                        LABEL(2,I), LABEL(3,I), LABEL(4,I), LABEL(5,I), COEFF(I)
            END DO
         END DO
      END DO

      CALL ALCBUF(3)
      END PROGRAM ANGDUMP2
