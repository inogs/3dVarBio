subroutine tao_minimizer

  use drv_str
  implicit none

! include for petsc and tao stuffs
#include "tao_minimizer.h"

  PetscErrorCode  ::   ierr
  Tao             ::   tao
  Vec             ::   x
  Mat             ::   A
  integer         ::   size, rank
  call MPI_Comm_size(PETSC_COMM_WORLD, size, ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)

  write(drv%dia,*) ''
  write(drv%dia,*) 'call PetscInitialize'
  print*,'call PetscInitialize'
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  print*, 'PetscInitialize() done'
  write(drv%dia,*) 'PETSC init done :)'

  write(drv%dia,*) "Within tao_minimizer subroutine!"
   
  write(drv%dia,*) 'call PetscFinalize'
  call PetscFinalize(ierr)
  write(drv%dia,*) 'PetscFinalize done'
  write(drv%dia,*) ''

end subroutine tao_minimizer
