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
     MyValues(j) = 0. !ctl%x_c(j)
  end do
  
  ! Setting only local values (since each process can access at all entries of MyState)
  do j=1, ctl%n
     call VecSetValues(MyState, 1, loc(j), MyValues(j), INSERT_VALUES, ierr)
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
     print*, "Minimization done with ", drv%MyCounter
     print*, "iterations"
     print*, ""
  endif
  print*,"MyRank ", MyRank, " exiting from tao_minimizer"
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
  PetscInt, allocatable, dimension(:)     :: loc
  ! PetscScalar, allocatable, dimension(:)  :: my_grad
  PetscScalar, pointer, dimension(:)  :: my_grad
  PetscScalar, pointer                    :: xtmp(:)
  PetscInt                                :: GlobalStart(1), MyEnd(1)
  ALLOCATE(loc(ctl%n)) !, my_grad(ctl%n))
  ! ALLOCATE(my_grad(ctl%n))

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
  call VecGetOwnershipRange(Grad, GlobalStart, MyEnd, ierr)

  ! print*,""
  print*,"MyRank ", MyRank, "GlobStart ", GlobalStart, "MyEnd ", MyEnd
  ! print*,""

  call VecGetArrayF90(Grad, my_grad, ierr)

  do j = 1, ctl%n
     my_grad(j) = ctl%g_c(j)
  end do

  call VecRestoreArrayF90(Grad, my_grad, ierr)

  ! call VecView(Grad, PETSC_VIEWER_STDOUT_SELF, ierr)

  ! do j = 1, ctl%n
  !    loc(j) = GlobalStart(1) + j - 1
  !    my_grad(j) = ctl%g_c(j)
  ! end do
  ! do j=1, ctl%n
  !    call VecSetValues(Grad, 1, loc(j), my_grad(j), INSERT_VALUES, ierr)
  ! end do
  ! call VecSetValues(Grad, ctl%n, loc, my_grad, INSERT_VALUES, ierr)
  ! CHKERRQ(ierr)
  ! call VecAssemblyBegin(Grad, ierr)
  ! CHKERRQ(ierr)
  ! call VecAssemblyEnd(Grad, ierr)
  ! CHKERRQ(ierr)

  DEALLOCATE(loc) !, my_grad)
  ! DEALLOCATE(my_grad)

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
#include "petsc/finclude/petscvecdef.h"

  use ctl_str
  use tao_str
  use petscvec
  use tao_str

  implicit none
  
#include "petsc/finclude/petsctao.h"

  Tao        :: tao
  Vec        :: lb, ub
  integer    :: dummy, ierr, j

  PetscInt, allocatable, dimension(:)       :: loc
  PetscScalar, allocatable, dimension(:)    :: lbound, ubound
  PetscInt                                  :: GlobalStart, MyEnd

  ALLOCATE(loc(ctl%n), lbound(ctl%n), ubound(ctl%n))
  call VecGetOwnershipRange(lb, GlobalStart, MyEnd, ierr)
  do j = 1, ctl%n
     loc(j) = GlobalStart + j - 1
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

#include "petsc/finclude/petscvecdef.h"

  use ctl_str
  use tao_str
  use petscvec

  implicit none
  
#include "petsc/finclude/petsctao.h"

  Tao                  :: tao
  integer              :: dummy, ierr, j, n, M, CheckVal
  Vec                  :: TmpGrad
  PetscScalar, pointer :: ReadGrad(:)
  PetscScalar          :: MyTol, grtol, gttol

  ! set useful variables
  n = ctl%n
  M = ctl%n_global
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
     print*,"AAAAAAAAAAAAAAAAAAAAA"
     call TaoSetConvergedReason(tao, TAO_CONVERGED_USER, ierr)
     CHKERRQ(ierr)
  end if

end subroutine MyConvTest
