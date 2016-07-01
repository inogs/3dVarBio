subroutine tao_minimizer

  use drv_str

  implicit none
#include "finclude/petsctao.h"
#include "finclude/petscsys.h"

  PetscErrorCode       ierr
  Tao                  tao
  
  write(drv%dia,*) ''
  write(drv%dia,*) 'call PetscInitialize'
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  write(drv%dia,*) 'PETSC init done :)'

  write(drv%dia,*) "Within tao_minimizer subroutine!"
   
  write(drv%dia,*) 'call PetscFinalize'
  call PetscFinalize(ierr)
  write(drv%dia,*) 'PetscFinalize done'
  write(drv%dia,*) ''

end subroutine tao_minimizer
