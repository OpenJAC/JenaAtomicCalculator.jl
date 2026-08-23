!   Minimal DLAMCH, sufficient for GRASP's SETMC, which asks only for
!   'L' (max base-10 exponent), 'O' (overflow), 'U' (underflow), 'E' (eps).
      DOUBLE PRECISION FUNCTION DLAMCH(CMACH)
      IMPLICIT NONE
      CHARACTER :: CMACH, C
      C = CMACH
      IF (C == 'l') C = 'L'
      IF (C == 'o') C = 'O'
      IF (C == 'u') C = 'U'
      IF (C == 'e') C = 'E'
      SELECT CASE (C)
      CASE ('L');  DLAMCH = DBLE(MAXEXPONENT(1.0D0)) * LOG10(2.0D0)
      CASE ('O');  DLAMCH = HUGE(1.0D0)
      CASE ('U');  DLAMCH = TINY(1.0D0)
      CASE ('E');  DLAMCH = EPSILON(1.0D0)
      CASE DEFAULT; DLAMCH = 0.0D0
      END SELECT
      RETURN
      END FUNCTION DLAMCH
