subroutine tao_minimizer
  
  use drv_str
  use ctl_str
  
  implicit none
  
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

  drv%MyCounter = 0
  
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

  write(drv%dia,*) 'Finalizing...'

  !
  ! TaoGetSolutionVector() ???
  !
  
  DEALLOCATE(loc, myvalues)
  
  call TaoDestroy(tao, ierr)
  CHKERRQ(ierr)

  call VecDestroy(x, ierr)
  CHKERRQ(ierr)

  call PetscFinalize(ierr)
  write(drv%dia,*) 'PetscFinalize done'
  write(drv%dia,*) ''
  
  print*, "Minimization done with ", drv%MyCounter
  print*, "iterations"

end subroutine tao_minimizer

subroutine MyFuncAndGradient(tao, x, f, g, dummy, ierr)
  
  use set_knd
  use drv_str
  use obs_str
  use grd_str
  use eof_str
  use ctl_str
  
  implicit none
  
#include "tao_minimizer.h"
  Tao             ::   tao
  Vec             ::   x, g
  PetscReal       ::   f
  integer         ::   dummy, ierr, j

  ! Working arrays
  PetscInt, allocatable, dimension(:)     :: loc
  PetscScalar, allocatable, dimension(:)  :: my_grad
  PetscScalar, pointer                    :: xtmp(:)

  ALLOCATE(loc(ctl%n), my_grad(ctl%n))

  ! print*, ""
  ! print*, "Im here, within MyFuncAndGradient :)"
  ! print*, ""

  call VecGetArrayReadF90(x, xtmp, ierr)
  CHKERRQ(ierr)

  do j=1,ctl%n
     ctl%x_c(j) = xtmp(j)
  end do

  ! compute function and gradient
  call costf

  ! assign to f the value computed by costf
  f = ctl%f_c

  do j = 1, ctl%n
     loc(j) = j-1
     my_grad(j) = ctl%g_c(j)
  end do

  call VecSetValues(g, ctl%n, loc, my_grad, INSERT_VALUES, ierr)
  CHKERRQ(ierr)
  call VecAssemblyBegin(g, ierr)
  CHKERRQ(ierr)
  call VecAssemblyEnd(g, ierr)
  CHKERRQ(ierr)

  DEALLOCATE(loc, my_grad)

  ! Update counter
  drv%MyCounter = drv%MyCounter + 1
  ! Exit without errors
  ierr = 0

end subroutine MyFuncAndGradient
