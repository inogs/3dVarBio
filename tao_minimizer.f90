subroutine tao_minimizer

  use drv_str
  use ctl_str
  
  implicit none
  
#include "tao_minimizer.h"
  
  PetscErrorCode  ::   ierr
  Tao             ::   tao
  Vec             ::   MyState ! array that stores the (temporary) state
  PetscInt        ::   n, M, MyStart, MyEnd
  PetscReal       ::   MyTolerance
  integer         ::   size, rank, j
  
  ! Working arrays
  PetscInt, allocatable, dimension(:)     :: loc
  PetscScalar, allocatable, dimension(:)  :: MyValues
  PetscScalar, pointer                    :: xtmp(:)
  
  external MyFuncAndGradient, MyBounds, MyConvTest
  
  print*,'Initialize Petsc and Tao stuffs'  
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  CHKERRQ(ierr)

  call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  
  print*, 'PetscInitialize() done by rank ', rank, ctl%n, ctl%n_global

  write(drv%dia,*) ''
  write(drv%dia,*) "Within tao_minimizer subroutine!"

  ! Allocate working arrays
  n = ctl%n
  M = ctl%n_global
  ALLOCATE(loc(n), MyValues(n))
  
  ! Take values from ctl%x_c in order to initialize 
  ! the solution array for Tao solver
  do j = 1, ctl%n
     loc(j) = j-1
     MyValues(j) = 0. !ctl%x_c(j)
  end do
  
  ! Create MyState array and fill it
  call VecCreateMPI(MPI_COMM_WORLD, n, M, MyState, ierr)
  call VecGetOwnershipRange(MyState, MyStart, MyEnd, ierr)
  print*, "MyState initialization by rank ", rank, "with indices: ", MyStart, MyEnd

  ! Setting only local values (since each process can access at all entries of MyState)
  do j=MyStart, MyEnd-1
     call VecSetValues(MyState, 1, j, MyValues(j-MyStart+1), INSERT_VALUES, ierr)
  end do

  ! call VecSetValues(MyState, ctl%n, loc, MyValues, INSERT_VALUES, ierr)
  call VecAssemblyBegin(MyState, ierr)
  call VecAssemblyEnd(MyState, ierr)

  ! Counter init
  drv%MyCounter = 0

  ! Create Tao object and set type BLMVM (ones that use BFGS minimization algorithm)
  call TaoCreate(MPI_COMM_WORLD, tao, ierr)
  CHKERRQ(ierr)
  call TaoSetType(tao,"blmvm",ierr)
  CHKERRQ(ierr)
  
  ! Set initial solution array, MyBounds and MyFuncAndGradient routines
  call TaoSetInitialVector(tao, MyState, ierr)
  CHKERRQ(ierr)
  call TaoSetVariableBoundsRoutine(tao, MyBounds, PETSC_NULL_OBJECT, ierr)
  CHKERRQ(ierr)
  call TaoSetObjectiveAndGradientRoutine(tao, MyFuncAndGradient, PETSC_NULL_OBJECT, ierr)
  CHKERRQ(ierr)

  ! Set MyTolerance and ConvergenceTest
  MyTolerance = 1.0d1
  ! MyTolerance = 2.0d-2
  call TaoSetTolerances(tao, MyTolerance, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, ierr)
  CHKERRQ(ierr)
  call TaoSetConvergenceTest(tao, MyConvTest, PETSC_NULL_OBJECT, ierr)
  CHKERRQ(ierr)

  ! Perform minimization
  call TaoSolve(tao, ierr)
  CHKERRQ(ierr)
  
  print*, ''
  print*, 'Tao Solver Info:'
  print*, ''
  call TaoView(tao, PETSC_VIEWER_STDOUT_WORLD, ierr)
  print*, ''

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
  write(drv%dia,*) 'Minimization done with ', drv%MyCounter
  write(drv%dia,*) 'iterations'
  write(drv%dia,*) ''
  
  print*, "Minimization done with ", drv%MyCounter
  print*, "iterations"
  print*, ""

end subroutine tao_minimizer

!-------------------------------------------------!
! subroutine that, given a state                  !
! MyState (provided by Tao solver),               !
! computes the value of cost function             !
! in MyState and its gradient                     !
!-------------------------------------------------!

subroutine MyFuncAndGradient(tao, MyState, CostFunc, Grad, dummy, ierr)
  
  use set_knd
  use drv_str
  use obs_str
  use grd_str
  use eof_str
  use ctl_str
  
  implicit none
  
#include "tao_minimizer.h"
  Tao             ::   tao
  Vec             ::   MyState, Grad
  PetscReal       ::   CostFunc
  integer         ::   dummy, ierr, j

  ! Working arrays
  PetscInt, allocatable, dimension(:)     :: loc
  PetscScalar, allocatable, dimension(:)  :: my_grad
  PetscScalar, pointer                    :: xtmp(:)

  ALLOCATE(loc(ctl%n), my_grad(ctl%n))

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

  ! assign the gradient value computed by costf to Grad
  do j = 1, ctl%n
     loc(j) = j-1
     my_grad(j) = ctl%g_c(j)
  end do

  call VecSetValues(Grad, ctl%n, loc, my_grad, INSERT_VALUES, ierr)
  CHKERRQ(ierr)
  call VecAssemblyBegin(Grad, ierr)
  CHKERRQ(ierr)
  call VecAssemblyEnd(Grad, ierr)
  CHKERRQ(ierr)

  DEALLOCATE(loc, my_grad)

  ! Update counter
  drv%MyCounter = drv%MyCounter + 1

  ! Exit without errors
  ierr = 0

end subroutine MyFuncAndGradient

!-------------------------------------------------!
! Subroutine that sets upper and lower            !
! bounds of the solution state array              !
! This routine is called only one time at the     !
! beginning of the iteration                      !
!-------------------------------------------------!

subroutine MyBounds(tao, lb, ub, dummy, ierr)
  use ctl_str

  implicit none
#include "tao_minimizer.h"
  Tao        :: tao
  Vec        :: lb, ub
  integer    :: dummy, ierr, j

  PetscInt, allocatable, dimension(:)       :: loc
  PetscScalar, allocatable, dimension(:)    :: lbound, ubound

  ALLOCATE(loc(ctl%n), lbound(ctl%n), ubound(ctl%n))
  
  do j = 1, ctl%n
     loc(j) = j-1
     lbound(j) = ctl%l_c(j)
     ubound(j) = ctl%u_c(j)
  end do

  call VecSetValues(lb, ctl%n, loc, lbound, INSERT_VALUES, ierr)
  CHKERRQ(ierr)
  call VecAssemblyBegin(lb, ierr)
  CHKERRQ(ierr)
  call VecAssemblyEnd(lb, ierr)
  CHKERRQ(ierr)

  call VecSetValues(ub, ctl%n, loc, ubound, INSERT_VALUES, ierr)
  CHKERRQ(ierr)
  call VecAssemblyBegin(ub, ierr)
  CHKERRQ(ierr)
  call VecAssemblyEnd(ub, ierr)
  CHKERRQ(ierr)

  DEALLOCATE(loc, lbound, ubound)
  
end subroutine MyBounds


!-------------------------------------------------!
! subroutine that performs the computation of     !
! the infinity norm of the gradient. If that norm !
! is less than the provided tolerance             !
! the solution is convergent                      !
!-------------------------------------------------!

subroutine MyConvTest(tao, dummy, ierr)

  use ctl_str

  implicit none
#include "tao_minimizer.h"
  Tao                  :: tao
  integer              :: dummy, ierr, j, n, M, CheckVal
  Vec                  :: TmpGrad
  PetscScalar, pointer :: ReadGrad(:)
  PetscReal            :: MyTol, grtol, gttol

  ! set useful variables
  n = ctl%n
  M = ctl%n
  CheckVal = 0

  ! taking tolerance value (skipping useless values)
  call TaoGetTolerances(tao, MyTol, PETSC_NULL_REAL, PETSC_NULL_REAL, ierr)
  CHKERRQ(ierr)

  call TaoGetGradientVector(tao, TmpGrad, ierr)
  CHKERRQ(ierr)  

  call VecGetArrayReadF90(TmpGrad, ReadGrad, ierr)
  CHKERRQ(ierr)

  ! check with infinity norm of gradient...
  do j=1, ctl%n
     if( ReadGrad(j) .gt. MyTol ) then
        CheckVal = 1
        EXIT
     end if
  end do

  call VecRestoreArrayReadF90(TmpGrad, ReadGrad, ierr)
  CHKERRQ(ierr)

  if( CheckVal .eq. 1) then
     call TaoSetConvergedReason(tao, TAO_CONTINUE_ITERATING, ierr)
     CHKERRQ(ierr)
  else
     call TaoSetConvergedReason(tao, TAO_CONVERGED_USER, ierr)
     CHKERRQ(ierr)
  end if

end subroutine MyConvTest
