PROGRAM test_f_exit_message
  INTEGER exit_code
  CHARACTER*(200) exit_message

  WRITE(*, *) 'Enter desired exit code:'
  READ(*, *) exit_code
  WRITE(*, *) 'Enter desired exit message:'
  READ(*, *) exit_message
  WRITE(*, *) 'Exiting with exit code', exit_code
  call f_exit_message(exit_code, exit_message)
END PROGRAM test_f_exit_message
