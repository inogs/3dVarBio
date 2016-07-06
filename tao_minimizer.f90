subroutine tao_minimizer
  
  use drv_str
  use ctl_str
  
  implicit none
  
  ! include for petsc and tao stuffs
#include "tao_minimizer.h"
  
  PetscErrorCode  ::   ierr
  Tao             ::   tao
  Vec             ::   x
  PetscInt        ::   n, M
  integer         ::   size, rank, j
  
  ! Working arrays
  PetscInt, allocatable, dimension(:)     :: loc
  PetscScalar, allocatable, dimension(:)  :: myvalues
  
  external MyFuncAndGradient
  
  write(drv%dia,*) ''
  write(drv%dia,*) 'call PetscInitialize'
  print*,'call PetscInitialize'
  
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  
  call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  
  ! Fill working arrays
  n = ctl%n
  M = ctl%n
  ALLOCATE(loc(n), myvalues(n))
  
  do j = 1, ctl%n
     loc(j) = j-1
     myvalues(j) = ctl%x_c(j)
  end do
  
  call VecCreateMPI(MPI_COMM_WORLD, n, M, x, ierr)
  
  call VecSetValues(x, ctl%n, loc, myvalues, INSERT_VALUES, ierr)
  call VecAssemblyBegin(x, ierr)
  call VecAssemblyEnd(x, ierr)
  
  print*, 'PetscInitialize() done by rank ', rank
  write(drv%dia,*) 'PETSC init done :)'
  
  write(drv%dia,*) "Within tao_minimizer subroutine!"
  ! call VecView(x, PETSC_VIEWER_STDOUT_WORLD, ierr)
  
  call TaoCreate(MPI_COMM_WORLD, tao, ierr)
  CHKERRQ(ierr)
  call TaoSetType(tao,"blmvm",ierr)
  CHKERRQ(ierr)
  
  call TaoSetInitialVector(tao, x, ierr)
  CHKERRQ(ierr)

  call TaoSetObjectiveAndGradientRoutine(tao, MyFuncAndGradient, PETSC_NULL_OBJECT, ierr)
  CHKERRQ(ierr)
  
  call TaoSolve(tao, ierr)
  CHKERRQ(ierr)

  write(drv%dia,*) 'call PetscFinalize'
  
  DEALLOCATE(loc, myvalues)
  
  call TaoDestroy(tao, ierr)
  CHKERRQ(ierr)

  call VecDestroy(x, ierr)
  CHKERRQ(ierr)

  call PetscFinalize(ierr)
  write(drv%dia,*) 'PetscFinalize done'
  write(drv%dia,*) ''

end subroutine tao_minimizer

subroutine MyFuncAndGradient(tao, x, f, g, dummy, ierr)
  
  use set_knd
  use obs_str
  use grd_str
  use eof_str
  use ctl_str
  
  implicit none
  
#include "tao_minimizer.h"
  Tao             ::   tao
  Vec             ::   x, g
  PetscReal       ::   f
  integer         ::   dummy, ierr

  print*, ""
  print*, "Im here, within MyFuncAndGradient :)"
  print*, ""

end subroutine MyFuncAndGradient
