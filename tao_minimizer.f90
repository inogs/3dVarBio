subroutine tao_minimizer

  use drv_str

  implicit none
#include "finclude/petsctao.h"
! #include "finclude/petscsys.h"

  PetscErrorCode       ierr
  Tao                  tao
  
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  
  write(drv%dia,*) "Within tao_minimizer subroutine!"
  
  call PetscFinalize(ierr)

end subroutine tao_minimizer
