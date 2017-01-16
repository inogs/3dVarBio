subroutine tao_minimizer

#include "petsc/finclude/petscvecdef.h"

  use drv_str
  use ctl_str
  use mpi_str
  use petscvec
  use tao_str

  implicit none

#include "petsc/finclude/petsctao.h"

  PetscErrorCode  ::   ierr
  Tao             ::   tao
  Vec             ::   MyState    ! array that stores the (temporary) state
  PetscInt        ::   n, M, GlobalStart, MyEnd
  PetscScalar     ::   MyTolerance
  integer(i4)     ::   j
  real(8)         ::   MaxGrad

  ! Working arrays
  PetscInt, allocatable, dimension(:)     :: loc
  PetscScalar, allocatable, dimension(:)  :: MyValues
  PetscScalar, pointer                    :: xtmp(:)

  external MyFuncAndGradient, MyBounds, MyConvTest

  if(MyRank .eq. 0) print*,'Initialize Petsc and Tao stuffs'

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  CHKERRQ(ierr)

  print*, 'PetscInitialize(...) done by MyRank ', MyRank, ctl%n, ctl%n_global

  if(MyRank .eq. 0) then
     write(drv%dia,*) ''
     write(drv%dia,*) "Within tao_minimizer subroutine!"
  endif

  ! Allocate working arrays
  n = ctl%n
  M = ctl%n_global
  NewCtl%n_global = ctl%n_global

  ALLOCATE(loc(n), MyValues(n))

  ! Create MyState array and fill it
  call VecCreateMPI(MPI_COMM_WORLD, n, M, MyState, ierr)
  call VecGetOwnershipRange(MyState, GlobalStart, MyEnd, ierr)

  print*, "MyState initialization by MyRank ", MyRank, "with indices: ", GlobalStart, MyEnd

  if( ctl%n .ne. MyEnd - GlobalStart ) then
     print*, ""
     print*, "WARNING!!"
     print*, "ctl%n .ne. GlobalStart - MyEnd"
     print*, "ctl%n = ", ctl%n
     print*, "GlobalStart = ", GlobalStart
     print*, "MyEnd = ", MyEnd
     print*, ""
  endif

  ! Take values from ctl%x_c in order to initialize
  ! the solution array for Tao solver
  do j = 1, ctl%n
     loc(j) = GlobalStart + j - 1
     MyValues(j) = 0.
  end do

  ! Setting only local values (since each process can access at all entries of MyState)
  call VecSetValues(MyState, ctl%n, loc, MyValues, INSERT_VALUES, ierr)
  call VecAssemblyBegin(MyState, ierr)
  call VecAssemblyEnd(MyState, ierr)

  ! Counter init
  drv%MyCounter = 0

  ! Create Tao object and set type BLMVM (ones that use BFGS minimization algorithm)
  call TaoCreate(MPI_COMM_WORLD, tao, ierr)
  CHKERRQ(ierr)
  ! call TaoSetType(tao,"blmvm",ierr)
  call TaoSetType(tao,"lmvm",ierr)
  CHKERRQ(ierr)

  ! Set initial solution array, MyBounds and MyFuncAndGradient routines
  call TaoSetInitialVector(tao, MyState, ierr)
  CHKERRQ(ierr)
  ! call TaoSetVariableBoundsRoutine(tao, MyBounds, PETSC_NULL_OBJECT, ierr)
  ! CHKERRQ(ierr)
  call TaoSetObjectiveAndGradientRoutine(tao, MyFuncAndGradient, PETSC_NULL_OBJECT, ierr)
  CHKERRQ(ierr)

  ! Set MyTolerance and ConvergenceTest
  call parallel_costf
  MaxGrad = 0
  do j=1,ctl%n
     MaxGrad = max(MaxGrad, abs(ctl%g_c(j)))
  end do

  call MPI_Allreduce(MPI_IN_PLACE, MaxGrad, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
  MyTolerance = ctl%pgper * MaxGrad
  if(MyRank .eq. 0) then
     print*, "Setting MyTolerance", MyTolerance
     write(drv%dia,*) "Setting MyTolerance", MyTolerance
  endif

  call TaoSetTolerances(tao, MyTolerance, 2.3d-2, ctl%pgper, ierr) !PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, ierr) !
  CHKERRQ(ierr)

  ! Perform minimization
  call TaoSolve(tao, ierr)
  CHKERRQ(ierr)

  if(MyRank .eq. 0) then
     print*, ''
     print*, 'Tao Solver Info:'
     print*, ''
  endif

  call TaoView(tao, PETSC_VIEWER_STDOUT_WORLD, ierr)

  ! Take computed solution and copy into ctl%x_c array
  call TaoGetSolutionVector(tao, MyState, ierr)
  CHKERRQ(ierr)
  call VecGetArrayReadF90(MyState, xtmp, ierr)
  CHKERRQ(ierr)

  do j = 1, ctl%n
     ctl%x_c(j) = xtmp(j)
  end do

  call VecRestoreArrayReadF90(MyState, xtmp, ierr)
  CHKERRQ(ierr)

  ! Deallocating variables
  DEALLOCATE(loc, MyValues)

  call TaoDestroy(tao, ierr)
  CHKERRQ(ierr)

  call VecDestroy(MyState, ierr)
  CHKERRQ(ierr)

  call PetscFinalize(ierr)
  if(MyRank .eq. 0) then
     write(drv%dia,*) 'Minimization done with ', drv%MyCounter
     write(drv%dia,*) 'iterations'
     write(drv%dia,*) ''

     print*, ""
     print*, "Minimization done with ", drv%MyCounter, "iterations"
     print*, ""
  endif

end subroutine tao_minimizer

!-------------------------------------------------!
! subroutine that, given a state                  !
! MyState (provided by Tao solver),               !
! computes the value of cost function             !
! in MyState and its gradient                     !
!-------------------------------------------------!

subroutine MyFuncAndGradient(tao, MyState, CostFunc, Grad, dummy, ierr)

#include "petsc/finclude/petscvecdef.h"

  use set_knd
  use drv_str
  use ctl_str
  use tao_str
  use petscvec
  use mpi_str
  use tao_str

  implicit none

#include "petsc/finclude/petsctao.h"

  Tao             ::   tao
  Vec             ::   MyState, Grad
  PetscScalar     ::   CostFunc
  PetscErrorCode  ::   ierr
  integer(i4)     ::   dummy, j

  ! Working arrays
  PetscScalar, pointer, dimension(:)  :: my_grad
  PetscScalar, pointer, dimension(:)  :: xtmp

  ! read temporary state provided by Tao Solver
  ! and set it in ctl%x_c array in order to compute
  ! the actual value of Cost Function and the gradient
  call VecGetArrayReadF90(MyState, xtmp, ierr)
  CHKERRQ(ierr)

  do j=1,ctl%n
     ctl%x_c(j) = xtmp(j)
  end do

  call VecRestoreArrayReadF90(MyState, xtmp, ierr)
  CHKERRQ(ierr)

  ! compute function and gradient
  call parallel_costf

  ! assign the Cost Function value computed by costf to CostFunc
  CostFunc = ctl%f_c

  call VecGetArrayF90(Grad, my_grad, ierr)

  do j = 1, ctl%n
     my_grad(j) = ctl%g_c(j)
  end do

  call VecRestoreArrayF90(Grad, my_grad, ierr)

  ! Update counter
  drv%MyCounter = drv%MyCounter + 1

  ! Exit without errors
  ierr = 0

end subroutine MyFuncAndGradient
