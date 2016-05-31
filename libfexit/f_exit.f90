SUBROUTINE f_exit(exit_code)
  IMPLICIT NONE
  INTEGER :: exit_code
  WRITE(*, *) 'CALL F_EXIT(', exit_code, ') - exiting with exit code', exit_code
  call c_exit(exit_code)
END SUBROUTINE f_exit

SUBROUTINE f_exit_message(exit_code, exit_message)
  IMPLICIT NONE
  INTEGER :: exit_code
  CHARACTER*(*) :: exit_message
  WRITE(*, *) TRIM(exit_message), ' - exiting with exit code', exit_code
  call c_exit(exit_code)
END SUBROUTINE f_exit_message
